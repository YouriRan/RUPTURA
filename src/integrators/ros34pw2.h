#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

#include "block_banded.h"
#include "column.h"
#include "jacobian.h"
#include "timing.h"

/**
 * \brief Adaptive column integrator with ROS34PW2, the four-stage Rosenbrock-W method of Rang and Angermann
 *        (BIT Numer. Math. 45 (2005) 761-787): third order, L-stable, and stiffly accurate.
 *
 * The stages advance the whole column state with the right-hand side of RK3 (precompute, then computeDerivatives).
 * Their linear systems use W = I - gamma h J~ with an approximate Jacobian J~ (jacobian.h): the sorption Jacobian on
 * the concentrations and physisorbed loadings, whose transport rows are exact for both convection schemes; the heat
 * exchange, convection, and conduction on the temperatures; zero on every other state and between the two blocks. A
 * W-method keeps its order for any matrix in the stage systems, so what J~ leaves out (velocity closures,
 * chemisorption, pore transport, reactions, heats of sorption) costs no accuracy: it is treated explicitly and can only
 * limit the stable step. Each step assembles J~ once. The loadings are eliminated node by node, which leaves W as a
 * block-banded matrix on the concentrations, and one on the temperatures. The step size follows from the embedded
 * second-order solution in the weighted RMS norm of CVODE, with a relative tolerance and per-block absolute
 * tolerances; a step with an error estimate above one is repeated with a smaller step.
 *
 * The continuity closures of InletPressureInletVelocity and InletVelocityOutletPressure march the molar balance of
 * each cell, so that the right-hand side keeps the total concentration sum_j c_ij of a cell at P/(R T_g): its
 * increment is -(sum_j c_ij / T_g) times that of the gas temperature, plus that of a pressure ramp. J~ leaves out
 * the closure, so that W would mix the total concentration into the stage increments, and over many steps it drifts
 * from P/(R T_g). W therefore acts as the identity on sum_j c_ij + (sum_j c_ij / T_g) T_g of each cell, at the
 * ratio of the step start: the stage increment of the total concentration is that of the right-hand side, less the
 * change that W makes to the temperature increment times sum_j c_ij / T_g, and the correction goes to the components
 * in proportion to their mole fractions at the step start. This is a W-method with J~ projected on the closure, which
 * keeps its order and the invariants of the right-hand side. InletPressureOutletPressure balances the molar flux of
 * the pressure profile instead, and there W keeps acting on the total concentration: the projection held the steps
 * of a dilute BEA breakthrough at a fiftieth of their length.
 */
struct ROS34PW2
{
  /// Coefficients in the stage form (I - gamma h J~) K_i = gamma h F(t + alpha_i h, U_i) + gamma sum_j c_ij K_j with
  /// U_i = u + sum_j a_ij K_j, the step u + sum_i m_i K_i, and the error estimate sum_i e_i K_i. They follow from the
  /// tables A, Gamma, b, and the embedded b^ of PETSc's TSROSWRA34PW2 (src/ts/impls/rosw/rosw.c) by the substitution
  /// K = Gamma k: a = A Gamma^-1, c = -Gamma^-1 below the diagonal, m = b Gamma^-1, e = (b - b^) Gamma^-1, and
  /// alpha_i = sum_j A_ij. Lower-triangular entries are stored row by row: (2,1), (3,1), (3,2), (4,1), (4,2), (4,3).
  /// The time derivative dF/dt of a Rosenbrock method is part of the Jacobian that a W-method may approximate; it is
  /// left out.
  static constexpr size_t stages = 4;
  static constexpr double gamma = 4.3586652150845900e-01;
  static constexpr std::array<double, 6> a{2.0,
                                           1.4192173174557646500,
                                           -0.25923221167296971378,
                                           4.1847604823191607312,
                                           -0.28519201735549591370,
                                           2.2942803602790417167};
  static constexpr std::array<double, 6> c{-4.5885607205580834861, -4.1847604823191607312, 0.28519201735549591370,
                                           -6.3681792001283577635, -6.7956209444668361844, 2.8700986043310560892};
  static constexpr std::array<double, 4> m{4.1847604823191601453, -0.28519201735549587377, 2.2942803602790413955, 1.0};
  static constexpr std::array<double, 4> e{0.27774994764796811038, -1.4032398951759990242, 1.7726301276675507452, 0.5};
  static constexpr std::array<double, 4> alpha{0.0, 0.87173304301691801, 0.73157995778885238, 1.0};

  /**
   * \brief Computes the stage increments K_i of one step from caller-supplied stage operations.
   *
   * On entry, function holds F(time, state); it is overwritten. evaluate(t, u, F) writes the right-hand side at a
   * stage state into F, and solve(x) overwrites x with W^-1 x. stageState is workspace of the state size.
   */
  template <typename Evaluate, typename Solve>
  static void stageIncrements(double time, double timeStep, std::span<const double> state, std::span<double> function,
                              std::span<double> stageState, const std::array<std::span<double>, stages>& increments,
                              const Evaluate& evaluate, const Solve& solve)
  {
    for (size_t stage = 0; stage < stages; ++stage)
    {
      const size_t first = stage * (stage - 1) / 2;  // start of this stage's row in a and c
      if (stage > 0)
      {
        for (size_t row = 0; row < state.size(); ++row)
        {
          double sum = state[row];
          for (size_t j = 0; j < stage; ++j) sum += a[first + j] * increments[j][row];
          stageState[row] = sum;
        }
        evaluate(time + alpha[stage] * timeStep, std::span<const double>(stageState), function);
      }
      for (size_t row = 0; row < state.size(); ++row)
      {
        double coupling = 0.0;
        for (size_t j = 0; j < stage; ++j) coupling += c[first + j] * increments[j][row];
        increments[stage][row] = gamma * timeStep * function[row] + gamma * coupling;
      }
      solve(increments[stage]);
    }
  }

  explicit ROS34PW2(const InputReader& inputReader);

  /**
   * \brief Constructs the integrator with the given tolerances; maximumTimeStep 0 leaves the step uncapped.
   */
  ROS34PW2(double relativeTolerance, double absoluteToleranceConcentration, double absoluteToleranceLoading,
           double absoluteToleranceTemperature, double maximumTimeStep = 0.0);

  /**
   * \brief Advances the column from startTime to endTime in s. On entry and on return the column holds the fields of
   *        a precompute at its state.
   */
  void advance(Column& column, double startTime, double endTime, Timing& timings);

  /**
   * \brief Forgets the step size, after a swing phase changes the boundary conditions: the next step starts small.
   */
  void reinitialize() noexcept { timeStep = 0.0; }

  /**
   * \brief Prints the step statistics.
   */
  void printStatistics() const;

  double relativeTolerance;               ///< Relative tolerance, applied to every state component.
  double absoluteToleranceConcentration;  ///< Absolute tolerance for bulk, surface, and pore concentrations in mol/m^3.
  double absoluteToleranceLoading;        ///< Absolute tolerance for physisorbed and chemisorbed loadings in mol/kg.
  double absoluteToleranceTemperature;    ///< Absolute tolerance for gas, solid, and wall temperatures in K.
  double maximumTimeStep;                 ///< Largest step in s; 0 leaves the step uncapped.

  double timeStep{0.0};                     ///< Step for the next attempt in s; 0 before the first step.
  size_t numberOfSteps{0};                  ///< Accepted steps.
  size_t numberOfRejectedSteps{0};          ///< Steps repeated with a smaller step.
  size_t numberOfEvaluations{0};            ///< Right-hand-side evaluations.
  size_t numberOfJacobians{0};              ///< Assemblies of J~.
  size_t numberOfMixturePredictions{0};     ///< Column-wide equilibrium evaluations (precompute calls).

 private:
  /**
   * \brief Sizes the workspace for the column, and sets the absolute tolerance of every state component.
   */
  void allocate(const Column& column);

  /**
   * \brief Factorizes W at gamma timeStep from the assembled J~; false if a block of it is singular.
   */
  bool factorize(double timeStep);

  /**
   * \brief Overwrites the state-sized vector x with W^-1 x; W is the identity outside the concentrations, physisorbed
   *        loadings, and (with an energy balance) temperatures.
   */
  void solve(std::span<double> x) const;

  ColumnStateLayout layout;                    ///< Blocks of the column state.
  bool energyBalance{false};                   ///< The temperatures are integrated and have a block in W.
  double stageGammaStep{0.0};                  ///< gamma h of the factorized W.
  std::optional<SorptionJacobian> sorptionJacobian;  ///< J~ on the concentrations and physisorbed loadings.
  std::vector<double> loadingDiagonal;         ///< 1 - gamma h dr_i/dq_i per node and component.
  BlockBandedMatrix concentrationMatrix;       ///< W on the concentrations with the loadings eliminated, factorized.
  /// Blocks above the diagonal of W that the Jacobian fills: two for WENO3 with some face velocity below zero, else one.
  size_t concentrationUpperBlocks{1};
  BlockBandedMatrix thermalJacobian;           ///< J~ on the temperatures, node blocks (T_g, T_s, T_w).
  BlockBandedMatrix thermalMatrix;             ///< W on the temperatures, factorized.
  mutable std::vector<double> thermalVector;   ///< Node-major temperatures of a solve.
  bool keepTotalConcentration{false};          ///< A gas with a marching continuity closure: W keeps sum_j c_ij.
  std::vector<double> stepMoleFractions;       ///< Mole fractions at the step start, grid-major.
  std::vector<double> stepExpansion;           ///< sum_j c_ij / T_g of each cell at the step start (energy balance).
  mutable std::vector<double> totalIncrement;  ///< Right-hand-side total concentration of each cell in a solve.
  std::vector<double> state;                   ///< State at the step start.
  std::vector<double> rates;                   ///< Right-hand side at the step start.
  std::vector<double> function;                ///< Right-hand side of the current stage.
  std::vector<double> stageState;              ///< State of the current stage, then the new state.
  std::vector<double> absoluteTolerance;       ///< Absolute tolerance of every state component.
  std::array<std::vector<double>, stages> increments;  ///< Stage increments K_i.
};

#pragma once

#include <array>
#include <functional>
#include <span>
#include <vector>

/**
 * \brief Coefficients of the three-stage, third-order, L-stable ROS3 method.
 *
 * Numerical coefficients are from KPP's Ros3 subroutine, at commit
 * c0a9b02466355461ba096facd5b0d9d83ef00b68 in int/rosenbrock.f90.
 * See docs/SIRK3_METHOD_VALIDATION.md for source links and the stage convention.
 * The lower-triangular entries are stored in the order (2,1), (3,1), (3,2).
 */
struct Rosenbrock3Coefficients
{
  static constexpr double gamma = 0.43586652150845899941601945119356;
  static constexpr std::array<double, 3> a{1.0, 1.0, 0.0};
  static constexpr std::array<double, 3> c{-1.0156171083877702091975600115545, 4.0759956452537699824805835358067,
                                           9.2076794298330791242156818474003};
  static constexpr std::array<double, 3> m{1.0, 6.1697947043828245592553615689730, -0.42772256543218573326238373806514};
  static constexpr std::array<double, 3> e{0.5, -2.9079558716805469821718236208017, 0.22354069897811569627360909276199};
  static constexpr std::array<double, 3> alpha{0.0, gamma, gamma};
  static constexpr std::array<double, 3> timeDerivativeWeights{gamma, 0.24291996454816804366592249683314,
                                                               2.1851380027664058511513169485832};
};

using Rosenbrock3Rhs = std::function<void(double, std::span<const double>, std::span<double>)>;
using Rosenbrock3Jacobian = std::function<void(double, std::span<const double>, std::span<double>)>;

/**
 * \brief Trial result of a ROS3 step, without modifying the supplied initial state.
 */
struct Rosenbrock3Step
{
  std::vector<double> state;  ///< Third-order solution at t + h.
  std::vector<double> error;  ///< Third-order solution minus the embedded second-order solution.
};

/**
 * \brief Advances a small ODE system by one ROS3 step using analytical callbacks.
 *
 * The Jacobian callback fills J[row + column * n] = dF[row]/du[column], in column-major order.
 * Callback outputs are initially zero; callbacks must set all nonzero entries.
 * partialTimeDerivative supplies dF/dt at fixed state, evaluated at the step start.
 * Omit it only for an autonomous problem. No numerical differentiation is performed.
 *
 * A dense pivoted LAPACK LU factorization of I - gamma*h*J is reused for all three stage solves.
 * This reference kernel does not perform error control, clipping, or column boundary handling.
 * The initial state is unchanged on success or failure. Invalid inputs, nonfinite evaluations,
 * and failed linear solves throw an exception, so a caller can reject the trial step.
 */
Rosenbrock3Step rosenbrock3Step(double time, double timeStep, std::span<const double> state, const Rosenbrock3Rhs& rhs,
                                const Rosenbrock3Jacobian& jacobian, const Rosenbrock3Rhs& partialTimeDerivative = {});

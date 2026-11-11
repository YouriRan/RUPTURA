#pragma once

#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

#include "column.h"
#include "jacobian.h"
#include "timing.h"

#if BUILD_SUNDIALS
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

/**
 * \brief CVODE-based variable-step integrator wrapper.
 *
 * Owns the SUNDIALS objects needed to integrate the Column ODE system when
 * the project is built with SUNDIALS support.
 */
struct CVODE
{
  /**
   * \brief Linear solver backing CVODE's modified Newton iteration.
   *
   * Dense assembles the full Jacobian by difference quotients and factors it directly: O(N^3) per setup and
   * O(N^2) storage, which becomes the dominant cost well before a few hundred grid points. SPGMR is
   * matrix-free, forming only Jacobian-vector products, so it costs one right-hand-side evaluation per Krylov
   * iteration and stores nothing. Neither choice affects the computed solution, only how fast and how
   * reliably the Newton iteration converges.
   */
  enum class LinearSolverType
  {
    Dense = 0,  ///< Direct dense solve with a difference-quotient Jacobian.
    SPGMR = 1,  ///< Matrix-free scaled preconditioned GMRES.
    KLU = 2     ///< Sparse direct solve with the analytical sorption Jacobian of jacobian.h.
  };

  /**
   * \brief Left preconditioner of SPGMR.
   *
   * SorptionJacobian factorizes I - gamma J with the banded sorption Jacobian of jacobian.h (the transport of the
   * convection scheme, velocity held fixed) on the concentration and physisorbed-loading blocks, and is the identity on
   * the other states
   * (chemisorption, pore transport, temperatures). It is an approximation of CVODE's iteration matrix, so it changes
   * how fast the Krylov iteration converges, not the solution.
   */
  enum class PreconditionerType
  {
    None = 0,             ///< Unpreconditioned.
    SorptionJacobian = 1  ///< Banded LU of I - gamma J on the concentrations and physisorbed loadings.
  };

  /**
   * \brief CVODE user data: the column and this integrator, for the right-hand side and the preconditioner, and the
   *        timing of the running advance, which every right-hand side evaluation is part of.
   */
  struct UserData
  {
    Column* column{nullptr};
    CVODE* integrator{nullptr};
    Timing* timings{nullptr};
    std::string lastError{};  ///< Message of the last exception that a callback turned into a CVODE error flag.
  };

  /**
   * \brief Constructs the CVODE wrapper from parsed input settings.
   */
  CVODE([[maybe_unused]] const InputReader& inputReader)
#if BUILD_SUNDIALS
      : relativeTolerance(inputReader.cvodeRelativeTolerance),
        absoluteToleranceConcentration(inputReader.cvodeAbsoluteToleranceConcentration),
        absoluteToleranceLoading(inputReader.cvodeAbsoluteToleranceLoading),
        absoluteToleranceTemperature(inputReader.cvodeAbsoluteToleranceTemperature),
        maximumTimeStep(inputReader.cvodeMaximumTimeStep),
        maximumOrder(static_cast<int>(inputReader.cvodeMaximumOrder)),
        linearSolverType(inputReader.cvodeLinearSolver == 2   ? LinearSolverType::KLU
                         : inputReader.cvodeLinearSolver == 1 ? LinearSolverType::SPGMR
                                                              : LinearSolverType::Dense),
        krylovDimension(static_cast<int>(inputReader.cvodeKrylovDimension)),
        preconditionerType(inputReader.cvodePreconditioner == 1 ? PreconditionerType::SorptionJacobian
                                                                 : PreconditionerType::None)
#endif
  {
  }

  /**
   * \brief Constructs the CVODE wrapper with the default tolerances and solver settings.
   */
  CVODE() = default;

  /**
   * \brief Releases any SUNDIALS resources owned by this wrapper.
   */
  ~CVODE();

  /**
   * \brief Copy construction is disabled because SUNDIALS handles are uniquely owned.
   */
  CVODE(const CVODE&) = delete;

  /**
   * \brief Copy assignment is disabled because SUNDIALS handles are uniquely owned.
   */
  CVODE& operator=(const CVODE&) = delete;

#if BUILD_SUNDIALS
  SUNContext sunContext = nullptr;           ///< SUNDIALS context handle.
  SUNLogger sunLogger = nullptr;             ///< SUNDIALS logger handle.
  N_Vector stateVector = nullptr;            ///< CVODE state vector.
  N_Vector stateDerivativeVector = nullptr;  ///< CVODE derivative vector.
  N_Vector absoluteToleranceVector = nullptr;  ///< Per-component absolute tolerances.
  SUNMatrix linearMatrix = nullptr;          ///< Linear-system matrix.
  void* cvodeMem = nullptr;                  ///< CVODE solver memory block.
  SUNNonlinearSolver solver = nullptr;       ///< Nonlinear solver handle.
  SUNLinearSolver linSolver = nullptr;       ///< Linear solver handle.
  sunrealtype currentTime = 0.0;             ///< Current absolute CVODE time.

  sunrealtype relativeTolerance = 1.0e-5;  ///< Relative integration tolerance, applied to every component.

  // The state vector mixes units (mol/m^3, mol/kg, K), so a single scalar absolute tolerance cannot serve
  // all blocks. These are applied per ColumnStateLayout block via CVodeSVtolerances.
  sunrealtype absoluteToleranceConcentration = 1.0e-6;  ///< Absolute tolerance for bulk/surface/pore concentrations.
  sunrealtype absoluteToleranceLoading = 1.0e-7;       ///< Absolute tolerance for physisorption/chemisorption loadings.
  sunrealtype absoluteToleranceTemperature = 1.0e-4;    ///< Absolute tolerance for gas/solid/wall temperatures.

  /// Maximum internal step in s: negative selects the automatic value (four grid-cell advection times, 4 dz/v),
  /// zero leaves CVODE uncapped, positive is used verbatim.
  double maximumTimeStep = -1.0;
  double appliedMaximumTimeStep = 0.0;  ///< Step cap actually handed to CVODE in s; 0 when uncapped.
  int maximumOrder = 5;                 ///< Maximum BDF order (1-5); CVODE's default is 5.

  LinearSolverType linearSolverType = LinearSolverType::SPGMR;  ///< Selected linear solver.
  int krylovDimension = 30;                                     ///< SPGMR Krylov subspace dimension.
  PreconditionerType preconditionerType = PreconditionerType::None;  ///< SPGMR preconditioner.

  UserData userData;                                ///< Handed to CVODE as user data.
  std::optional<SorptionJacobian> sorptionJacobian;  ///< Jacobian of the preconditioner, kept while CVODE reuses it.
  bool sorptionJacobianAssembled = false;            ///< sorptionJacobian holds a Jacobian that CVODE may reuse.
  std::vector<double> preconditionerBand;            ///< LU factors of I - gamma J, band storage of SorptionJacobian.
  std::vector<int> preconditionerPivots;             ///< Row interchanges of that LU.
  std::vector<double> preconditionerVector;          ///< Node-major work vector of the preconditioner solve.
  N_Vector preconditionerScratch = nullptr;          ///< Receives the right-hand side that a preconditioner setup
                                                     ///< evaluates for the column fields; never read.
  std::vector<sunindextype> sparseColumnStarts;      ///< KLU: CSC column pointers of the Jacobian pattern.
  std::vector<sunindextype> sparseRowIndices;        ///< KLU: CSC row indices of the Jacobian pattern.
  std::vector<size_t> sparseBandIndex;               ///< KLU: band-storage index of each nonzero, or SIZE_MAX.

  /**
   * \brief Evaluates the Column ODE right-hand side for CVODE; user_data is a CVODE::UserData.
   */
  static int evaluateDerivatives(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector, void* user_data);

  /**
   * \brief Assembles (unless CVODE allows reuse) and factorizes the SorptionJacobian preconditioner I - gamma J.
   */
  static int setupPreconditioner(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector,
                                 sunbooleantype jacobianOk, sunbooleantype* jacobianCurrent, sunrealtype gamma,
                                 void* user_data);

  /**
   * \brief KLU: fills the sparse Jacobian with the sorption Jacobian at the trial state; zero rows elsewhere.
   */
  static int evaluateJacobian(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector, SUNMatrix matrix,
                              void* user_data, N_Vector scratch1, N_Vector scratch2, N_Vector scratch3);

  /**
   * \brief Solves (I - gamma J) z = r on the concentration and loading blocks; z = r on the others.
   */
  static int solvePreconditioner(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector, N_Vector r,
                                 N_Vector z, sunrealtype gamma, sunrealtype delta, int side, void* user_data);
#endif

  /**
   * \brief Advances the column from the current CVODE time to endTime in s with one CVode call.
   *
   * CVODE chooses its own internal steps (capped by maximumTimeStep) and interpolates the state at endTime.
   */
  void advance(Column& column, double endTime, Timing& timings);

  /**
   * \brief Initializes CVODE state and solver data for the given column.
   */
  void initialize(Column& column);

  /**
   * \brief Resets CVODE integration history after an external state or boundary change.
   */
  void reinitialize();

  /**
   * \brief Resets CVODE integration history and re-derives the automatic step cap (one grid-cell advection time)
   * from the column's current velocity, e.g. after a swing phase changes the inlet velocity.
   */
  void reinitialize(const Column& column);

  /**
   * \brief Prints CVODE solver statistics, or nothing when the solver was never initialized.
   */
  void printStatistics() const;
};

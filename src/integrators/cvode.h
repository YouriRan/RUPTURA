#pragma once

#include <span>
#include <stdexcept>

#include "column.h"
#include "component.h"
#include "compute.h"
#include "mixture_prediction.h"
#include "sorption.h"
#include "timing.h"
#include "utils.h"

struct MultibedColumn;

#if BUILD_SUNDIALS
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_dense.h>
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_klu.h>
#include <sunlinsol/sunlinsol_spbcgs.h>
#include <sunlinsol/sunlinsol_spfgmr.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunlinsol/sunlinsol_sptfqmr.h>
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
    SPGMR = 1   ///< Matrix-free scaled preconditioned GMRES, currently unpreconditioned.
  };

  /**
   * \brief Constructs the CVODE wrapper from parsed input settings.
   */
  CVODE(const InputReader& inputReader)
      : timeStep(inputReader.timeStep),
        autoNumberOfSteps(inputReader.autoNumberOfTimeSteps),
        numberOfSteps(inputReader.numberOfTimeSteps)
#if BUILD_SUNDIALS
        ,
        relativeTolerance(inputReader.cvodeRelativeTolerance),
        absoluteToleranceConcentration(inputReader.cvodeAbsoluteToleranceConcentration),
        absoluteToleranceLoading(inputReader.cvodeAbsoluteToleranceLoading),
        absoluteToleranceTemperature(inputReader.cvodeAbsoluteToleranceTemperature),
        maximumTimeStep(inputReader.cvodeMaximumTimeStep),
        linearSolverType(inputReader.cvodeLinearSolver == 1 ? LinearSolverType::SPGMR : LinearSolverType::Dense),
        krylovDimension(static_cast<int>(inputReader.cvodeKrylovDimension))
#endif
  {
  }

  /**
   * \brief Constructs the CVODE wrapper from explicit time-step settings.
   */
  CVODE(double timeStep, bool autoNumberOfSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoNumberOfSteps(autoNumberOfSteps), numberOfSteps(numberOfSteps)
  {
  }

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

  double timeStep;         ///< Integration time step in s.
  bool autoNumberOfSteps;  ///< Continue until breakthrough criterion when true.
  size_t numberOfSteps;    ///< Requested number of integration steps.

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
  sunrealtype absoluteToleranceLoading = 1.0e-8;       ///< Absolute tolerance for physisorption/chemisorption loadings.
  sunrealtype absoluteToleranceTemperature = 1.0e-4;    ///< Absolute tolerance for gas/solid/wall temperatures.

  /// Maximum internal step in s: negative selects the automatic value (one grid-cell advection time, dz/v),
  /// zero leaves CVODE uncapped, positive is used verbatim.
  double maximumTimeStep = -1.0;
  double appliedMaximumTimeStep = 0.0;  ///< Step cap actually handed to CVODE in s; 0 when uncapped.

  LinearSolverType linearSolverType = LinearSolverType::Dense;  ///< Selected linear solver.
  int krylovDimension = 30;                                     ///< SPGMR Krylov subspace dimension.

  /**
   * \brief Evaluates the Column ODE right-hand side for CVODE.
   */
  static int evaluateDerivatives(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector, void* user_data);

  /**
   * \brief Evaluates the MultibedColumn ODE right-hand side for CVODE.
   */
  static int evaluateMultibedDerivatives(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector,
                                         void* user_data);
#endif

  /**
   * \brief Advances the column by one CVODE output step.
   */
  bool propagate(Column& column, size_t step, Timing& timings);

  /**
   * \brief Advances a multibed column by one CVODE output step.
   */
  bool propagate(MultibedColumn& column, size_t step, Timing& timings);

  /**
   * \brief Initializes CVODE state and solver data for the given column.
   */
  void initialize(Column& column);

  /**
   * \brief Initializes CVODE state and solver data for the given multibed column.
   */
  void initialize(MultibedColumn& column);

  /**
   * \brief Resets CVODE integration history after an external state or boundary change.
   */
  void reinitialize();

  /**
   * \brief Prints CVODE solver statistics, or nothing when the solver was never initialized.
   */
  void printStatistics() const;
};

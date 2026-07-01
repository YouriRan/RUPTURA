#pragma once

#include <span>
#include <stdexcept>

#include "column.h"
#include "component.h"
#include "compute.h"
#include "mixture_prediction.h"
#include "timing.h"
#include "utils.h"

#if BUILD_SUNDIALS
#include <cvode/cvode.h>
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
   * \brief Constructs the CVODE wrapper from parsed input settings.
   */
  CVODE(const InputReader& inputReader)
      : timeStep(inputReader.timeStep),
        autoNumberOfSteps(inputReader.autoNumberOfTimeSteps),
        numberOfSteps(inputReader.numberOfTimeSteps)
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
  SUNMatrix linearMatrix = nullptr;          ///< Linear-system matrix.
  void* cvodeMem = nullptr;                  ///< CVODE solver memory block.
  SUNNonlinearSolver solver = nullptr;       ///< Nonlinear solver handle.
  SUNLinearSolver linSolver = nullptr;       ///< Linear solver handle.

  const sunrealtype relativeTolerance = 1.0e-3;  ///< Relative integration tolerance.
  const sunrealtype absoluteTolerance = 1.0e-6;  ///< Absolute integration tolerance.

  /**
   * \brief Evaluates the Column ODE right-hand side for CVODE.
   */
  static int evaluateDerivatives(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector, void* user_data);
#endif

  /**
   * \brief Advances the column by one CVODE output step.
   */
  bool propagate(Column& column, size_t step, Timing& timings);

  /**
   * \brief Initializes CVODE state and solver data for the given column.
   */
  void initialize(Column& column);
};

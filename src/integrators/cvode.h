#pragma once

#include <span>
#include <stdexcept>

#include "column.h"
#include "component.h"
#include "compute.h"
#include "mixture_prediction.h"
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

struct CVODE
{
  CVODE(const InputReader& inputReader)
      : timeStep(inputReader.timeStep),
        autoSteps(inputReader.autoNumberOfTimeSteps),
        numberOfSteps(inputReader.numberOfTimeSteps)
  {
  }

  CVODE(double timeStep, bool autoSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoSteps(autoSteps), numberOfSteps(numberOfSteps)
  {
  }

  ~CVODE();

  CVODE(const CVODE&) = delete;
  CVODE& operator=(const CVODE&) = delete;

  double timeStep;
  bool autoSteps;
  size_t numberOfSteps;

#if BUILD_SUNDIALS
  SUNContext sunContext = nullptr;
  SUNLogger sunLogger = nullptr;
  N_Vector u = nullptr;
  N_Vector uDot = nullptr;
  SUNMatrix A = nullptr;
  void* cvodeMem = nullptr;
  SUNNonlinearSolver solver = nullptr;
  SUNLinearSolver linSolver = nullptr;

  const sunrealtype relativeTolerance = 1.0e-3;
  const sunrealtype absoluteTolerance = 1.0e-6;

  static int f(sunrealtype t, N_Vector u, N_Vector uDot, void* user_data);
#endif

  bool propagate(Column& column, size_t step);
  void initialize(Column& column);
};
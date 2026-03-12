#pragma once

#include <span>
#include <stdexcept>

#include "column.h"
#include "component.h"
#include "compute.h"
#include "mixture_prediction.h"

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
  CVODE(double timeStep, bool autoSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoSteps(autoSteps), numberOfSteps(numberOfSteps)
  {
  }

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
#endif

  bool propagate(Column& column, size_t step);
  void initialize(Column& column);
};

#if BUILD_SUNDIALS
inline void copyFromcolumn(const Column& column, N_Vector u);
inline void copyFromcolumnDot(const Column& column, N_Vector uDot);
inline void copyIntocolumn(Column& column, N_Vector u);
inline void copyIntocolumnDot(Column& column, N_Vector uDot);

inline std::span<double> getPartialPressureSpan(N_Vector u, size_t Ngrid, size_t Ncomp);
inline std::span<double> getAdsorptionSpan(N_Vector u, size_t Ngrid, size_t Ncomp);
inline std::span<double> getGasTemperatureSpan(N_Vector u, size_t Ngrid, size_t Ncomp);
inline std::span<double> getSolidTemperatureSpan(N_Vector u, size_t Ngrid, size_t Ncomp);
inline std::span<double> getWallTemperatureSpan(N_Vector u, size_t Ngrid, size_t Ncomp);

static int f(sunrealtype t, N_Vector u, N_Vector uDot, void* user_data);
#endif
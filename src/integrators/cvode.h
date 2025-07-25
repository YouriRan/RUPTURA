#pragma once

#include <cvode/cvode.h>             /* main integrator header file                 */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fct. and macros      */
#include <sundials/sundials_dense.h> /* use generic DENSE solver in preconditioning */
#include <sundials/sundials_logger.h>
#include <sundials/sundials_types.h> /* definition of sunrealtype                      */
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver            */
#include <sunlinsol/sunlinsol_spfgmr.h>  /* access to SPFGMR SUNLinearSolver            */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver             */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver           */
#include <sunmatrix/sunmatrix_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver         */

#include <span>

#include "breakthrough_state.h"
#include "component.h"
#include "compute.h"
#include "mixture_prediction.h"

struct CVODE
{
  CVODE(double timeStep, bool autoSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoSteps(autoSteps), numberOfSteps(numberOfSteps) {};

  double timeStep;
  bool autoSteps;
  size_t numberOfSteps;

  SUNContext sunContext;
  SUNLogger sunLogger;
  N_Vector u;
  N_Vector uDot;
  SUNMatrix A;
  void* cvodeMem;
  SUNNonlinearSolver solver;
  SUNLinearSolver linSolver;

  const sunrealtype relativeTolerance = 1.0e-3;
  const sunrealtype absoluteTolerance = 1.0e-6;

  bool propagate(BreakthroughState& state, size_t step);
  void initialize(BreakthroughState& state);
};

inline void copyFromState(const BreakthroughState& state, N_Vector u);
inline void copyFromStateDot(const BreakthroughState& state, N_Vector uDot);
inline void copyIntoState(BreakthroughState& state, N_Vector u);
inline void copyIntoStateDot(BreakthroughState& state, N_Vector uDot);

inline std::span<double> getTotalPressureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/);
inline std::span<double> getTemperatureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/);
inline std::span<double> getWallTemperatureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/);
inline std::span<double> getPartialPressureSpan(N_Vector u, size_t Ngrid, size_t Ncomp);
inline std::span<double> getAdsorptionSpan(N_Vector u, size_t Ngrid, size_t Ncomp);
inline std::span<double> getMoleFractionSpan(N_Vector u, size_t Ngrid, size_t Ncomp);

static int f(sunrealtype t, N_Vector u, N_Vector uDot, void* user_data);
#include "cvode.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <span>
#include <vector>

bool CVODE::propagate(BreakthroughState &state, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  double tNext = static_cast<double>(step + 1) * timeStep;
  size_t Ngrid = state.Ngrid;
  size_t Ncomp = state.Ncomp;

  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance = std::max(
          tolerance,
          std::abs((state.partialPressure[Ngrid * Ncomp + j] / (state.exitPressure * state.components[j].Yi0)) - 1.0));
    }

    // consider 1% as being visibily indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    if (tolerance < 0.01)
    {
      std::cout << "\nConvergence criteria reached, running 10% longer\n\n" << std::endl;
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoSteps = false;
    }
  }

  sunrealtype tReturn = t;
  CVode(cvodeMem, tNext, u, &tReturn, CV_NORMAL);
  unpackState(state, u);
  unpackStateDot(state, uDot);

  return (!autoSteps && step >= numberOfSteps - 1);
}

void CVODE::initialize(BreakthroughState &state)
{
  // initialize logger and context
  SUNContext_Create(SUN_COMM_NULL, &sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &sunLogger);
  SUNContext_SetLogger(sunContext, sunLogger);

  // create vector that cvode will operate on
  const sunindextype totalSize = static_cast<sunindextype>(3 * (state.Ncomp + 1) * (state.Ngrid + 1));
  u = N_VNew_Serial(totalSize, sunContext);
  packState(state, u);
  uDot = N_VNew_Serial(totalSize, sunContext);
  packStateDot(state, uDot);

  cvodeMem = CVodeCreate(CV_BDF, sunContext);
  CVodeSetMaxNumSteps(cvodeMem, 1e8);

  CVodeSetUserData(cvodeMem, &state);

  const sunrealtype t0 = 0.0;
  CVodeInit(cvodeMem, f, t0, u);

  CVodeSStolerances(cvodeMem, relativeTolerance, absoluteTolerance);
  solver = SUNNonlinSol_Newton(u, sunContext);
  CVodeSetNonlinearSolver(cvodeMem, solver);

  A = SUNDenseMatrix(totalSize, totalSize, sunContext);
  linSolver = SUNLinSol_Dense(u, A, sunContext);
  CVodeSetLinearSolver(cvodeMem, linSolver, A);

  CVodeSetJacFn(cvodeMem, nullptr);
}

inline std::span<double> getTotalPressureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  return {base + 0 * small, small};
}

inline std::span<double> getTemperatureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  return {base + 1 * small, small};
}

inline std::span<double> getWallTemperatureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  return {base + 2 * small, small};
}

inline std::span<double> getPartialPressureSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 3 * small + 0 * big, big};
}

inline std::span<double> getAdsorptionSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 3 * small + 1 * big, big};
}

inline std::span<double> getMoleFractionSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 3 * small + 2 * big, big};
}

inline void packState(const BreakthroughState &state, N_Vector u)
{
  std::copy(state.totalPressure.begin(), state.totalPressure.end(),
            getTotalPressureSpan(u, state.Ngrid, state.Ncomp).begin());
  std::copy(state.temperature.begin(), state.temperature.end(),
            getTemperatureSpan(u, state.Ngrid, state.Ncomp).begin());
  std::copy(state.wallTemperature.begin(), state.wallTemperature.end(),
            getWallTemperatureSpan(u, state.Ngrid, state.Ncomp).begin());

  std::copy(state.partialPressure.begin(), state.partialPressure.end(),
            getPartialPressureSpan(u, state.Ngrid, state.Ncomp).begin());
  std::copy(state.adsorption.begin(), state.adsorption.end(), getAdsorptionSpan(u, state.Ngrid, state.Ncomp).begin());
  std::copy(state.moleFraction.begin(), state.moleFraction.end(),
            getMoleFractionSpan(u, state.Ngrid, state.Ncomp).begin());
}

inline void packStateDot(const BreakthroughState &state, N_Vector uDot)
{
  std::copy(state.totalPressureDot.begin(), state.totalPressureDot.end(),
            getTotalPressureSpan(uDot, state.Ngrid, state.Ncomp).begin());
  std::copy(state.temperatureDot.begin(), state.temperatureDot.end(),
            getTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin());
  std::copy(state.wallTemperatureDot.begin(), state.wallTemperatureDot.end(),
            getWallTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin());

  std::copy(state.pressureDot.begin(), state.pressureDot.end(),
            getPartialPressureSpan(uDot, state.Ngrid, state.Ncomp).begin());
  std::copy(state.adsorptionDot.begin(), state.adsorptionDot.end(),
            getAdsorptionSpan(uDot, state.Ngrid, state.Ncomp).begin());
  std::copy(state.moleFractionDot.begin(), state.moleFractionDot.end(),
            getMoleFractionSpan(uDot, state.Ngrid, state.Ncomp).begin());
}

inline void unpackState(BreakthroughState &state, N_Vector u)
{
  std::copy(getTotalPressureSpan(u, state.Ngrid, state.Ncomp).begin(),
            getTotalPressureSpan(u, state.Ngrid, state.Ncomp).end(), state.totalPressure.begin());
  std::copy(getTemperatureSpan(u, state.Ngrid, state.Ncomp).begin(),
            getTemperatureSpan(u, state.Ngrid, state.Ncomp).end(), state.temperature.begin());
  std::copy(getWallTemperatureSpan(u, state.Ngrid, state.Ncomp).begin(),
            getWallTemperatureSpan(u, state.Ngrid, state.Ncomp).end(), state.wallTemperature.begin());

  std::copy(getPartialPressureSpan(u, state.Ngrid, state.Ncomp).begin(),
            getPartialPressureSpan(u, state.Ngrid, state.Ncomp).end(), state.partialPressure.begin());
  std::copy(getAdsorptionSpan(u, state.Ngrid, state.Ncomp).begin(),
            getAdsorptionSpan(u, state.Ngrid, state.Ncomp).end(), state.adsorption.begin());
  std::copy(getMoleFractionSpan(u, state.Ngrid, state.Ncomp).begin(),
            getMoleFractionSpan(u, state.Ngrid, state.Ncomp).end(), state.moleFraction.begin());
}

inline void unpackStateDot(BreakthroughState &state, N_Vector uDot)
{
  std::copy(getTotalPressureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getTotalPressureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.totalPressureDot.begin());
  std::copy(getTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getTemperatureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.temperatureDot.begin());
  std::copy(getWallTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getWallTemperatureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.wallTemperatureDot.begin());

  std::copy(getPartialPressureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getPartialPressureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.pressureDot.begin());
  std::copy(getAdsorptionSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getAdsorptionSpan(uDot, state.Ngrid, state.Ncomp).end(), state.adsorptionDot.begin());
  std::copy(getMoleFractionSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getMoleFractionSpan(uDot, state.Ngrid, state.Ncomp).end(), state.moleFractionDot.begin());
}

static int f(sunrealtype t, N_Vector u, N_Vector uDot, void *user_data)
{
  auto *state = reinterpret_cast<BreakthroughState *>(user_data);
  unpackState(*state, u);
  computeEquilibriumLoadings(*state);
  computeVelocity(*state);
  computeFirstDerivatives(*state);
  packStateDot(*state, uDot);

  return 0;
}

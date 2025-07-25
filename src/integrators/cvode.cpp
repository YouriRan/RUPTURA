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

  // maybe not necessary? Are there external things influencing these values?
  // copyFromState(state, u);
  // copyFromStateDot(state, u);

  CVode(cvodeMem, tNext, u, &tReturn, CV_NORMAL);

  copyIntoState(state, u);
  copyIntoStateDot(state, uDot);

  return (!autoSteps && step >= numberOfSteps - 1);
}

void CVODE::initialize(BreakthroughState &state)
{
  // initialize logger and context
  SUNContext_Create(SUN_COMM_NULL, &sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &sunLogger);
  SUNContext_SetLogger(sunContext, sunLogger);

  // create vector that cvode will operate on
  // const sunindextype totalSize = static_cast<sunindextype>(3 * (state.Ncomp + 1) * (state.Ngrid + 1));
  const sunindextype totalSize = static_cast<sunindextype>((2 * state.Ncomp) * (state.Ngrid + 1));
  u = N_VNew_Serial(totalSize, sunContext);
  copyFromState(state, u);
  uDot = N_VNew_Serial(totalSize, sunContext);
  copyFromStateDot(state, uDot);

  cvodeMem = CVodeCreate(CV_BDF, sunContext);
  CVodeSetMaxNumSteps(cvodeMem, 1e6);

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

// inline std::span<double> getTotalPressureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/)
// {
//   double *base = static_cast<double *>(N_VGetArrayPointer(u));
//   size_t small = Ngrid + 1;
//   return {base + 0 * small, small};
// }

// inline std::span<double> getTemperatureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/)
// {
//   double *base = static_cast<double *>(N_VGetArrayPointer(u));
//   size_t small = Ngrid + 1;
//   return {base + 1 * small, small};
// }

// inline std::span<double> getWallTemperatureSpan(N_Vector u, size_t Ngrid, size_t /*Ncomp*/)
// {
//   double *base = static_cast<double *>(N_VGetArrayPointer(u));
//   size_t small = Ngrid + 1;
//   return {base + 2 * small, small};
// }

inline std::span<double> getPartialPressureSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 0 * small + 0 * big, big};
}

inline std::span<double> getAdsorptionSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double *base = static_cast<double *>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 0 * small + 1 * big, big};
}

// inline std::span<double> getMoleFractionSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
// {
//   double *base = static_cast<double *>(N_VGetArrayPointer(u));
//   size_t small = Ngrid + 1;
//   size_t big = small * Ncomp;
//   return {base + 3 * small + 2 * big, big};
// }

inline void copyFromState(const BreakthroughState &state, N_Vector u)
{
  //   std::copy(state.totalPressure.begin(), state.totalPressure.end(),
  //             getTotalPressureSpan(u, state.Ngrid, state.Ncomp).begin());
  //   std::copy(state.temperature.begin(), state.temperature.end(),
  //             getTemperatureSpan(u, state.Ngrid, state.Ncomp).begin());
  //   std::copy(state.wallTemperature.begin(), state.wallTemperature.end(),
  //             getWallTemperatureSpan(u, state.Ngrid, state.Ncomp).begin());

  std::copy(state.partialPressure.begin(), state.partialPressure.end(),
            getPartialPressureSpan(u, state.Ngrid, state.Ncomp).begin());
  std::copy(state.adsorption.begin(), state.adsorption.end(), getAdsorptionSpan(u, state.Ngrid, state.Ncomp).begin());
  //   std::copy(state.moleFraction.begin(), state.moleFraction.end(),
  //             getMoleFractionSpan(u, state.Ngrid, state.Ncomp).begin());
}

inline void copyFromStateDot(const BreakthroughState &state, N_Vector uDot)
{
  //   std::copy(state.totalpartialPressureDot.begin(), state.totalpartialPressureDot.end(),
  //             getTotalPressureSpan(uDot, state.Ngrid, state.Ncomp).begin());
  //   std::copy(state.temperatureDot.begin(), state.temperatureDot.end(),
  //             getTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin());
  //   std::copy(state.wallTemperatureDot.begin(), state.wallTemperatureDot.end(),
  //             getWallTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin());

  std::copy(state.partialPressureDot.begin(), state.partialPressureDot.end(),
            getPartialPressureSpan(uDot, state.Ngrid, state.Ncomp).begin());
  std::copy(state.adsorptionDot.begin(), state.adsorptionDot.end(),
            getAdsorptionSpan(uDot, state.Ngrid, state.Ncomp).begin());
  //   std::copy(state.moleFractionDot.begin(), state.moleFractionDot.end(),
  //             getMoleFractionSpan(uDot, state.Ngrid, state.Ncomp).begin());
}

inline void copyIntoState(BreakthroughState &state, N_Vector u)
{
  //   std::copy(getTotalPressureSpan(u, state.Ngrid, state.Ncomp).begin(),
  //             getTotalPressureSpan(u, state.Ngrid, state.Ncomp).end(), state.totalPressure.begin());
  //   std::copy(getTemperatureSpan(u, state.Ngrid, state.Ncomp).begin(),
  //             getTemperatureSpan(u, state.Ngrid, state.Ncomp).end(), state.temperature.begin());
  //   std::copy(getWallTemperatureSpan(u, state.Ngrid, state.Ncomp).begin(),
  //             getWallTemperatureSpan(u, state.Ngrid, state.Ncomp).end(), state.wallTemperature.begin());

  std::copy(getPartialPressureSpan(u, state.Ngrid, state.Ncomp).begin(),
            getPartialPressureSpan(u, state.Ngrid, state.Ncomp).end(), state.partialPressure.begin());
  std::copy(getAdsorptionSpan(u, state.Ngrid, state.Ncomp).begin(),
            getAdsorptionSpan(u, state.Ngrid, state.Ncomp).end(), state.adsorption.begin());
  //   std::copy(getMoleFractionSpan(u, state.Ngrid, state.Ncomp).begin(),
  //             getMoleFractionSpan(u, state.Ngrid, state.Ncomp).end(), state.moleFraction.begin());
}

inline void copyIntoStateDot(BreakthroughState &state, N_Vector uDot)
{
  //   std::copy(getTotalPressureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
  //             getTotalPressureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.totalpartialPressureDot.begin());
  //   std::copy(getTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
  //             getTemperatureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.temperatureDot.begin());
  //   std::copy(getWallTemperatureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
  //             getWallTemperatureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.wallTemperatureDot.begin());

  std::copy(getPartialPressureSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getPartialPressureSpan(uDot, state.Ngrid, state.Ncomp).end(), state.partialPressureDot.begin());
  std::copy(getAdsorptionSpan(uDot, state.Ngrid, state.Ncomp).begin(),
            getAdsorptionSpan(uDot, state.Ngrid, state.Ncomp).end(), state.adsorptionDot.begin());
  //   std::copy(getMoleFractionSpan(uDot, state.Ngrid, state.Ncomp).begin(),
  //             getMoleFractionSpan(uDot, state.Ngrid, state.Ncomp).end(), state.moleFractionDot.begin());
}

static int f(sunrealtype t, N_Vector u, N_Vector uDot, void *user_data)
{
  auto *state = reinterpret_cast<BreakthroughState *>(user_data);

  // auto spanTotalPressure = getTotalPressureSpan(u, state->Ngrid, state->Ncomp);
  //   auto spanTemperature = getTemperatureSpan(u, state->Ngrid, state->Ncomp);
  //   auto spanWallTemperature = getWallTemperatureSpan(u, state->Ngrid, state->Ncomp);
  auto spanPartialPressure = getPartialPressureSpan(u, state->Ngrid, state->Ncomp);
  auto spanAdsorption = getAdsorptionSpan(u, state->Ngrid, state->Ncomp);
  //   auto spanMoleFraction = getMoleFractionSpan(u, state->Ngrid, state->Ncomp);

  auto spanAdsorptionDot = getAdsorptionSpan(uDot, state->Ngrid, state->Ncomp);
  auto spanpartialPressureDot = getPartialPressureSpan(uDot, state->Ngrid, state->Ncomp);

  computeEquilibriumLoadings(state->Ncomp, state->Ngrid, state->totalPressure, spanPartialPressure,
                             state->idealGasMolFractions, state->adsorbedMolFractions, state->numberOfMolecules,
                             state->cachedPressure, state->cachedGrandPotential, state->equilibriumAdsorption,
                             state->mixture, state->iastPerformance, state->pressureGradient, state->columnLength,
                             state->maxIsothermTerms);

  computeVelocity(state->Ncomp, state->Ngrid, state->resolution, state->interstitialGasVelocity,
                  state->columnEntranceVelocity, state->pressureGradient, state->totalPressure,
                  state->prefactorMassTransfer, state->equilibriumAdsorption, spanAdsorption, state->components,
                  spanPartialPressure);

  computeFirstDerivatives(state->Ncomp, state->Ngrid, state->resolution, spanPartialPressure,
                          state->equilibriumAdsorption, spanAdsorption, spanAdsorptionDot, spanpartialPressureDot,
                          state->interstitialGasVelocity, state->prefactorMassTransfer, state->components);

  return 0;
}

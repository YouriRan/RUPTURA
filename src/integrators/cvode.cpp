#include "cvode.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <span>
#include <vector>

bool CVODE::propagate(Column& column, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  double tNext = static_cast<double>(step + 1) * timeStep;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance = std::max(tolerance, std::abs((column.partialPressure[Ngrid * Ncomp + j] /
                                                (column.exitPressure * column.components[j].Yi0)) -
                                               1.0));
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
  // copyFromcolumn(column, u);
  // copyFromcolumnDot(column, u);

  CVode(cvodeMem, tNext, u, &tReturn, CV_NORMAL);

  copyIntocolumn(column, u);
  copyIntocolumnDot(column, uDot);

  return (!autoSteps && step >= numberOfSteps - 1);
}

void CVODE::initialize(Column& column)
{
  // initialize logger and context
  SUNContext_Create(SUN_COMM_NULL, &sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &sunLogger);
  SUNContext_SetLogger(sunContext, sunLogger);

  // create vector that cvode will operate on
  // const sunindextype totalSize = static_cast<sunindextype>(3 * (column.Ncomp + 1) * (column.Ngrid + 1));
  const sunindextype totalSize = static_cast<sunindextype>((2 * column.Ncomp) * (column.Ngrid + 1));
  u = N_VNew_Serial(totalSize, sunContext);
  copyFromcolumn(column, u);
  uDot = N_VNew_Serial(totalSize, sunContext);
  copyFromcolumnDot(column, uDot);

  cvodeMem = CVodeCreate(CV_BDF, sunContext);
  CVodeSetMaxNumSteps(cvodeMem, 1e6);

  CVodeSetUserData(cvodeMem, &column);

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
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 0 * small + 0 * big, big};
}

inline std::span<double> getAdsorptionSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
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

inline void copyFromcolumn(const Column& column, N_Vector u)
{
  //   std::copy(column.totalPressure.begin(), column.totalPressure.end(),
  //             getTotalPressureSpan(u, column.Ngrid, column.Ncomp).begin());
  //   std::copy(column.temperature.begin(), column.temperature.end(),
  //             getTemperatureSpan(u, column.Ngrid, column.Ncomp).begin());
  //   std::copy(column.wallTemperature.begin(), column.wallTemperature.end(),
  //             getWallTemperatureSpan(u, column.Ngrid, column.Ncomp).begin());

  std::copy(column.partialPressure.begin(), column.partialPressure.end(),
            getPartialPressureSpan(u, column.Ngrid, column.Ncomp).begin());
  std::copy(column.adsorption.begin(), column.adsorption.end(),
            getAdsorptionSpan(u, column.Ngrid, column.Ncomp).begin());
  //   std::copy(column.moleFraction.begin(), column.moleFraction.end(),
  //             getMoleFractionSpan(u, column.Ngrid, column.Ncomp).begin());
}

inline void copyFromcolumnDot(const Column& column, N_Vector uDot)
{
  //   std::copy(column.totalpartialPressureDot.begin(), column.totalpartialPressureDot.end(),
  //             getTotalPressureSpan(uDot, column.Ngrid, column.Ncomp).begin());
  //   std::copy(column.temperatureDot.begin(), column.temperatureDot.end(),
  //             getTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin());
  //   std::copy(column.wallTemperatureDot.begin(), column.wallTemperatureDot.end(),
  //             getWallTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin());

  std::copy(column.partialPressureDot.begin(), column.partialPressureDot.end(),
            getPartialPressureSpan(uDot, column.Ngrid, column.Ncomp).begin());
  std::copy(column.adsorptionDot.begin(), column.adsorptionDot.end(),
            getAdsorptionSpan(uDot, column.Ngrid, column.Ncomp).begin());
  //   std::copy(column.moleFractionDot.begin(), column.moleFractionDot.end(),
  //             getMoleFractionSpan(uDot, column.Ngrid, column.Ncomp).begin());
}

inline void copyIntocolumn(Column& column, N_Vector u)
{
  //   std::copy(getTotalPressureSpan(u, column.Ngrid, column.Ncomp).begin(),
  //             getTotalPressureSpan(u, column.Ngrid, column.Ncomp).end(), column.totalPressure.begin());
  //   std::copy(getTemperatureSpan(u, column.Ngrid, column.Ncomp).begin(),
  //             getTemperatureSpan(u, column.Ngrid, column.Ncomp).end(), column.temperature.begin());
  //   std::copy(getWallTemperatureSpan(u, column.Ngrid, column.Ncomp).begin(),
  //             getWallTemperatureSpan(u, column.Ngrid, column.Ncomp).end(), column.wallTemperature.begin());

  std::copy(getPartialPressureSpan(u, column.Ngrid, column.Ncomp).begin(),
            getPartialPressureSpan(u, column.Ngrid, column.Ncomp).end(), column.partialPressure.begin());
  std::copy(getAdsorptionSpan(u, column.Ngrid, column.Ncomp).begin(),
            getAdsorptionSpan(u, column.Ngrid, column.Ncomp).end(), column.adsorption.begin());
  //   std::copy(getMoleFractionSpan(u, column.Ngrid, column.Ncomp).begin(),
  //             getMoleFractionSpan(u, column.Ngrid, column.Ncomp).end(), column.moleFraction.begin());
}

inline void copyIntocolumnDot(Column& column, N_Vector uDot)
{
  //   std::copy(getTotalPressureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
  //             getTotalPressureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.totalpartialPressureDot.begin());
  //   std::copy(getTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
  //             getTemperatureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.temperatureDot.begin());
  //   std::copy(getWallTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
  //             getWallTemperatureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.wallTemperatureDot.begin());

  std::copy(getPartialPressureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
            getPartialPressureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.partialPressureDot.begin());
  std::copy(getAdsorptionSpan(uDot, column.Ngrid, column.Ncomp).begin(),
            getAdsorptionSpan(uDot, column.Ngrid, column.Ncomp).end(), column.adsorptionDot.begin());
  //   std::copy(getMoleFractionSpan(uDot, column.Ngrid, column.Ncomp).begin(),
  //             getMoleFractionSpan(uDot, column.Ngrid, column.Ncomp).end(), column.moleFractionDot.begin());
}

static int f(sunrealtype t, N_Vector u, N_Vector uDot, void* user_data)
{
  auto* column = reinterpret_cast<Column*>(user_data);

  // auto spanTotalPressure = getTotalPressureSpan(u, column->Ngrid, column->Ncomp);
  //   auto spanTemperature = getTemperatureSpan(u, column->Ngrid, column->Ncomp);
  //   auto spanWallTemperature = getWallTemperatureSpan(u, column->Ngrid, column->Ncomp);
  auto spanPartialPressure = getPartialPressureSpan(u, column->Ngrid, column->Ncomp);
  auto spanAdsorption = getAdsorptionSpan(u, column->Ngrid, column->Ncomp);
  //   auto spanMoleFraction = getMoleFractionSpan(u, column->Ngrid, column->Ncomp);

  auto spanAdsorptionDot = getAdsorptionSpan(uDot, column->Ngrid, column->Ncomp);
  auto spanpartialPressureDot = getPartialPressureSpan(uDot, column->Ngrid, column->Ncomp);

  computeEquilibriumLoadings(column->Ncomp, column->Ngrid, column->totalPressure, spanPartialPressure,
                             column->idealGasMolFractions, column->adsorbedMolFractions, column->numberOfMolecules,
                             column->cachedPressure, column->cachedGrandPotential, column->equilibriumAdsorption,
                             column->mixture, column->iastPerformance, column->pressureGradient, column->columnLength,
                             column->maxIsothermTerms);

  computeVelocity(column->Ncomp, column->Ngrid, column->resolution, column->interstitialGasVelocity,
                  column->columnEntranceVelocity, column->pressureGradient, column->totalPressure,
                  column->prefactorMassTransfer, column->equilibriumAdsorption, spanAdsorption, column->components,
                  spanPartialPressure);

  computeFirstDerivatives(column->Ncomp, column->Ngrid, column->resolution, spanPartialPressure,
                          column->equilibriumAdsorption, spanAdsorption, spanAdsorptionDot, spanpartialPressureDot,
                          column->interstitialGasVelocity, column->prefactorMassTransfer, column->components);

  return 0;
}

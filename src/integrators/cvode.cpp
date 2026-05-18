// cvode.cpp
#include "cvode.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#if BUILD_SUNDIALS

CVODE::~CVODE()
{
  if (linSolver) SUNLinSolFree(linSolver);
  if (A) SUNMatDestroy(A);
  if (solver) SUNNonlinSolFree(solver);
  if (uDot) N_VDestroy(uDot);
  if (u) N_VDestroy(u);
  if (cvodeMem) CVodeFree(&cvodeMem);
  if (sunLogger) SUNLogger_Destroy(&sunLogger);
  if (sunContext) SUNContext_Free(&sunContext);
}

inline std::span<double> getConcentrationSpan(N_Vector v, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = Ngrid + 1;
  const size_t big = small * Ncomp;
  return {base, big};
}

inline std::span<double> getAdsorptionSpan(N_Vector v, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = Ngrid + 1;
  const size_t big = small * Ncomp;
  return {base + big, big};
}

inline std::span<double> getGasTemperatureSpan(N_Vector v, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = Ngrid + 1;
  const size_t big = small * Ncomp;
  return {base + 2 * big, small};
}

inline std::span<double> getSolidTemperatureSpan(N_Vector v, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = Ngrid + 1;
  const size_t big = small * Ncomp;
  return {base + 2 * big + small, small};
}

inline std::span<double> getWallTemperatureSpan(N_Vector v, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = Ngrid + 1;
  const size_t big = small * Ncomp;
  return {base + 2 * big + 2 * small, small};
}

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
      tolerance = std::max(tolerance, std::abs((column.moleFraction[Ngrid * Ncomp + j] / column.components[j].Yi0) -
                                               1.0));
    }

    if (tolerance < 0.01)
    {
      std::cout << "\nConvergence criteria reached, running 10% longer\n\n" << std::endl;
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoSteps = false;
    }
  }

  sunrealtype tReturn = t;

  int flag = CVode(cvodeMem, tNext, u, &tReturn, CV_NORMAL);
  if (flag < 0)
  {
    throw std::runtime_error("CVode failed during propagation");
  }

  // Refresh derived quantities and stateDot at the accepted state.
  flag = CVODE::f(tReturn, u, uDot, &column);
  if (flag != 0)
  {
    throw std::runtime_error("CVODE::f failed after propagation");
  }

  return (!autoSteps && step >= numberOfSteps - 1);
}

void CVODE::initialize(Column& column)
{
  SUNContext_Create(SUN_COMM_NULL, &sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &sunLogger);
  SUNContext_SetLogger(sunContext, sunLogger);

  const sunindextype totalSize = static_cast<sunindextype>(column.state.size());

  u = N_VMake_Serial(totalSize, column.state.data(), sunContext);
  uDot = N_VMake_Serial(totalSize, column.stateDot.data(), sunContext);

  cvodeMem = CVodeCreate(CV_BDF, sunContext);
  CVodeSetMaxNumSteps(cvodeMem, 1e6);
  CVodeSetUserData(cvodeMem, &column);

  const sunrealtype t0 = 0.0;
  int flag = CVodeInit(cvodeMem, CVODE::f, t0, u);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeInit failed");

  flag = CVodeSStolerances(cvodeMem, relativeTolerance, absoluteTolerance);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSStolerances failed");

  solver = SUNNonlinSol_Newton(u, sunContext);
  flag = CVodeSetNonlinearSolver(cvodeMem, solver);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetNonlinearSolver failed");

  A = SUNDenseMatrix(totalSize, totalSize, sunContext);
  linSolver = SUNLinSol_Dense(u, A, sunContext);
  flag = CVodeSetLinearSolver(cvodeMem, linSolver, A);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");

  CVodeSetJacFn(cvodeMem, nullptr);
}

int CVODE::f(sunrealtype /*t*/, N_Vector u, N_Vector uDot, void* user_data)
{
  auto* column = static_cast<Column*>(user_data);

  auto spanConcentration = getConcentrationSpan(u, column->Ngrid, column->Ncomp);
  auto spanAdsorption = getAdsorptionSpan(u, column->Ngrid, column->Ncomp);

  auto spanConcentrationDot = getConcentrationSpan(uDot, column->Ngrid, column->Ncomp);
  auto spanAdsorptionDot = getAdsorptionSpan(uDot, column->Ngrid, column->Ncomp);

  auto spanGasTemperature = (column->energyBalance) ? getGasTemperatureSpan(u, column->Ngrid, column->Ncomp) : column->gasTemperature;
  auto spanSolidTemperature = (column->energyBalance) ? getSolidTemperatureSpan(u, column->Ngrid, column->Ncomp) : column->solidTemperature;
  auto spanWallTemperature = (column->energyBalance) ? getWallTemperatureSpan(u, column->Ngrid, column->Ncomp) : column->wallTemperature;
  auto spanGasTemperatureDot = (column->energyBalance) ? getGasTemperatureSpan(uDot, column->Ngrid, column->Ncomp) : column->gasTemperatureDot;
  auto spanSolidTemperatureDot = (column->energyBalance) ? getSolidTemperatureSpan(uDot, column->Ngrid, column->Ncomp) : column->solidTemperatureDot;
  auto spanWallTemperatureDot = (column->energyBalance) ? getWallTemperatureSpan(uDot, column->Ngrid, column->Ncomp) : column->wallTemperatureDot;
  
  computePressure(column->velocityProfile, column->boundaryCondition, column->components, column->Ngrid, column->Ncomp,
                  column->externalPressure, column->voidFraction, column->dynamicViscosity, column->particleDiameter,
                  column->resolution, column->interstitialGasVelocity, column->gasDensity, column->totalConcentration,
                  column->totalPressure, spanGasTemperature, spanConcentration, column->partialPressure,
                  column->moleFraction);

  computeEquilibriumLoadings(column->mixture, column->Ngrid, column->Ncomp, column->maxIsothermTerms,
                             column->iastPerformance, column->totalPressure, column->idealGasMolFractions,
                             column->adsorbedMolFractions, column->numberOfMolecules, column->equilibriumAdsorption,
                             column->moleFraction, column->cachedPressure, column->cachedGrandPotential);

  switch (column->velocityProfile)
  {
    case Column::VelocityProfile::FixedPressureGradient:
    {
      computeVelocityFixedGradient(column->boundaryCondition, column->components, column->Ngrid, column->Ncomp,
                                   column->pressureGradient, column->columnEntranceVelocity, column->resolution,
                                   column->prefactorMassTransfer, column->interstitialGasVelocity,
                                   column->totalPressure, column->totalConcentration, spanAdsorption,
                                   column->equilibriumAdsorption, spanConcentration);
      break;
    }
    case Column::VelocityProfile::Ergun:
    {
      computeVelocityErgun(column->boundaryCondition, column->Ngrid, column->voidFraction,
                           column->columnEntranceVelocity, column->columnLength, column->dynamicViscosity,
                           column->particleDiameter, column->resolution, column->interstitialGasVelocity,
                           column->gasDensity, column->totalPressure);
      break;
    }
    case Column::VelocityProfile::FixedVelocity:
      break;
    default:
      break;
  }

  if (column->energyBalance)
  {
    computeFirstDerivativesEnergyBalance(
        column->components, column->Ngrid, column->Ncomp, column->externalTemperature, column->voidFraction,
        column->particleDensity, column->particleDiameter, column->internalDiameter, column->outerDiameter,
        column->wallDensity, column->gasThermalConductivity, column->wallThermalConductivity,
        column->heatTransferGasSolid, column->heatTransferGasWall, column->heatTransferWallExternal,
        column->heatCapacityGas, column->heatCapacitySolid, column->heatCapacityWall, column->resolution,
        column->prefactorMassTransfer, column->interstitialGasVelocity, column->totalPressure, spanGasTemperature,
        spanGasTemperatureDot, spanSolidTemperature, spanSolidTemperatureDot, spanWallTemperature,
        spanWallTemperatureDot, spanConcentration, spanConcentrationDot, spanAdsorption, spanAdsorptionDot,
        column->equilibriumAdsorption, column->moleFraction, column->gasDensity, column->coeffGasGas,
        column->coeffGasSolid, column->coeffGasWall, column->coeffDiffusion, column->facePressures, column->massFlux);
  }
  else
  {
    computeFirstDerivatives(column->components, column->Ngrid, column->Ncomp, column->resolution,
                            column->prefactorMassTransfer, column->interstitialGasVelocity, spanConcentration,
                            spanConcentrationDot, spanAdsorption, spanAdsorptionDot, column->equilibriumAdsorption,
                            column->massFlux);
  }

  return 0;
}

#else

CVODE::~CVODE() = default;

bool CVODE::propagate(Column&, size_t)
{
  throw std::runtime_error("CVODE::propagate() called, but this build was compiled without SUNDIALS support");
}

void CVODE::initialize(Column&) {}

#endif
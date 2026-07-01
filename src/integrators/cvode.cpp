// cvode.cpp
#include "cvode.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <print>

#if BUILD_SUNDIALS

CVODE::~CVODE()
{
  if (linSolver) SUNLinSolFree(linSolver);
  if (linearMatrix) SUNMatDestroy(linearMatrix);
  if (solver) SUNNonlinSolFree(solver);
  if (stateDerivativeVector) N_VDestroy(stateDerivativeVector);
  if (stateVector) N_VDestroy(stateVector);
  if (cvodeMem) CVodeFree(&cvodeMem);
  if (sunLogger) SUNLogger_Destroy(&sunLogger);
  if (sunContext) SUNContext_Free(&sunContext);
}

inline std::span<double> getMoleFractionSpan(N_Vector v, size_t numberOfGridPoints, size_t numberOfComponents)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;
  return {base, big};
}

inline std::span<double> getAdsorptionSpan(N_Vector v, size_t numberOfGridPoints, size_t numberOfComponents)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;
  return {base + big, big};
}

inline std::span<double> getGasTemperatureSpan(N_Vector v, size_t numberOfGridPoints, size_t numberOfComponents)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;
  return {base + 2 * big, small};
}

inline std::span<double> getSolidTemperatureSpan(N_Vector v, size_t numberOfGridPoints, size_t numberOfComponents)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;
  return {base + 2 * big + small, small};
}

inline std::span<double> getWallTemperatureSpan(N_Vector v, size_t numberOfGridPoints, size_t numberOfComponents)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(v));
  const size_t small = numberOfGridPoints + 1;
  const size_t big = small * numberOfComponents;
  return {base + 2 * big + 2 * small, small};
}

bool CVODE::propagate(Column& column, size_t step, Timing& timings)
{
  double t = static_cast<double>(step) * timeStep;
  double tNext = static_cast<double>(step + 1) * timeStep;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  if (autoNumberOfSteps)
  {
    double tolerance = 0.0;

    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      tolerance = std::max(tolerance, std::abs((column.moleFraction[numberOfGridPoints * numberOfComponents + j] /
                                                column.components[j].initialGasMoleFraction) -
                                               1.0));
    }

    if (tolerance < 0.01)
    {
      std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoNumberOfSteps = false;
    }
  }

  sunrealtype tReturn = t;

  auto timer = timings.scoped(timings.total);

  int flag = CVode(cvodeMem, tNext, stateVector, &tReturn, CV_NORMAL);
  if (flag < 0)
  {
    throw std::runtime_error("CVode failed during propagation");
  }

  // Refresh derived quantities and stateDot at the accepted state.
  flag = CVODE::evaluateDerivatives(tReturn, stateVector, stateDerivativeVector, &column);
  if (flag != 0)
  {
    throw std::runtime_error("CVODE::evaluateDerivatives failed after propagation");
  }

  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}

void CVODE::initialize(Column& column)
{
  SUNContext_Create(SUN_COMM_NULL, &sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &sunLogger);
  SUNContext_SetLogger(sunContext, sunLogger);

  const sunindextype totalSize = static_cast<sunindextype>(column.state.size());

  stateVector = N_VMake_Serial(totalSize, column.state.data(), sunContext);
  stateDerivativeVector = N_VMake_Serial(totalSize, column.stateDot.data(), sunContext);

  cvodeMem = CVodeCreate(CV_BDF, sunContext);
  CVodeSetMaxNumSteps(cvodeMem, 1e6);
  CVodeSetUserData(cvodeMem, &column);

  const sunrealtype t0 = 0.0;
  int flag = CVodeInit(cvodeMem, CVODE::evaluateDerivatives, t0, stateVector);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeInit failed");

  flag = CVodeSStolerances(cvodeMem, relativeTolerance, absoluteTolerance);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSStolerances failed");

  solver = SUNNonlinSol_Newton(stateVector, sunContext);
  flag = CVodeSetNonlinearSolver(cvodeMem, solver);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetNonlinearSolver failed");

  linearMatrix = SUNDenseMatrix(totalSize, totalSize, sunContext);
  linSolver = SUNLinSol_Dense(stateVector, linearMatrix, sunContext);
  flag = CVodeSetLinearSolver(cvodeMem, linSolver, linearMatrix);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");

  CVodeSetJacFn(cvodeMem, nullptr);
}

int CVODE::evaluateDerivatives(sunrealtype /*t*/, N_Vector stateVector, N_Vector stateDerivativeVector, void* user_data)
{
  auto* column = static_cast<Column*>(user_data);

  auto spanMoleFraction = getMoleFractionSpan(stateVector, column->numberOfGridPoints, column->numberOfComponents);
  auto spanAdsorption = getAdsorptionSpan(stateVector, column->numberOfGridPoints, column->numberOfComponents);

  auto spanMoleFractionDot =
      getMoleFractionSpan(stateDerivativeVector, column->numberOfGridPoints, column->numberOfComponents);
  auto spanAdsorptionDot =
      getAdsorptionSpan(stateDerivativeVector, column->numberOfGridPoints, column->numberOfComponents);

  auto spanGasTemperature = (column->energyBalance) ? getGasTemperatureSpan(stateVector, column->numberOfGridPoints,
                                                                            column->numberOfComponents)
                                                    : column->gasTemperature;
  auto spanSolidTemperature = (column->energyBalance) ? getSolidTemperatureSpan(stateVector, column->numberOfGridPoints,
                                                                                column->numberOfComponents)
                                                      : column->solidTemperature;
  auto spanWallTemperature = (column->energyBalance) ? getWallTemperatureSpan(stateVector, column->numberOfGridPoints,
                                                                              column->numberOfComponents)
                                                     : column->wallTemperature;
  auto spanGasTemperatureDot =
      (column->energyBalance)
          ? getGasTemperatureSpan(stateDerivativeVector, column->numberOfGridPoints, column->numberOfComponents)
          : column->gasTemperatureDot;
  auto spanSolidTemperatureDot =
      (column->energyBalance)
          ? getSolidTemperatureSpan(stateDerivativeVector, column->numberOfGridPoints, column->numberOfComponents)
          : column->solidTemperatureDot;
  auto spanWallTemperatureDot =
      (column->energyBalance)
          ? getWallTemperatureSpan(stateDerivativeVector, column->numberOfGridPoints, column->numberOfComponents)
          : column->wallTemperatureDot;

  computePressure(column->components, column->velocityProfile, column->boundaryCondition, column->numberOfGridPoints,
                  column->numberOfComponents, column->inletPressure, column->outletPressure, column->voidFraction,
                  column->dynamicViscosity, column->particleDiameter, column->resolution,
                  column->interstitialGasVelocity, column->gasDensity, column->totalConcentration,
                  column->totalPressure, column->concentration, column->partialPressure, spanMoleFraction,
                  spanGasTemperature);

  computeEquilibriumLoadings(column->mixture, column->numberOfGridPoints, column->numberOfComponents,
                             column->maxIsothermTerms, column->iastPerformance, column->idealGasMolFractions,
                             column->adsorbedMolFractions, column->numberOfMolecules, column->totalPressure,
                             column->equilibriumAdsorption, column->cachedPressure, column->cachedGrandPotential,
                             spanMoleFraction, spanGasTemperature);

  switch (column->velocityProfile)
  {
    case Column::VelocityProfile::FixedPressureGradient:
    {
      computeVelocityFixedGradient(column->boundaryCondition, column->components, column->numberOfGridPoints,
                                   column->numberOfComponents, column->pressureGradient, column->columnEntranceVelocity,
                                   column->resolution, column->prefactorMassTransfer, column->interstitialGasVelocity,
                                   column->totalConcentration, column->totalPressure, column->equilibriumAdsorption,
                                   spanMoleFraction, spanAdsorption);
      break;
    }
    case Column::VelocityProfile::Ergun:
    {
      computeVelocityErgun(column->boundaryCondition, column->numberOfGridPoints, column->voidFraction,
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
    computeDerivativesEnergyBalance(
        column->components, column->numberOfGridPoints, column->numberOfComponents, column->externalTemperature,
        column->voidFraction, column->particleDensity, column->particleDiameter, column->internalDiameter,
        column->outerDiameter, column->wallDensity, column->gasThermalConductivity, column->wallThermalConductivity,
        column->heatTransferGasSolid, column->heatTransferGasWall, column->heatTransferWallExternal,
        column->heatCapacityGas, column->heatCapacitySolid, column->heatCapacityWall, column->resolution,
        column->prefactorMassTransfer, column->interstitialGasVelocity, column->gasDensity, column->totalConcentration,
        column->equilibriumAdsorption, column->coeffGasGas, column->coeffGasSolid, column->coeffGasWall,
        column->coeffDiffusion, spanMoleFraction, spanMoleFractionDot, spanAdsorption, spanAdsorptionDot,
        spanGasTemperature, spanGasTemperatureDot, spanSolidTemperature, spanSolidTemperatureDot, spanWallTemperature,
        spanWallTemperatureDot);
  }
  else
  {
    computeDerivatives(column->components, column->numberOfGridPoints, column->numberOfComponents, column->resolution,
                       column->prefactorMassTransfer, column->interstitialGasVelocity, column->totalConcentration,
                       column->equilibriumAdsorption, spanMoleFraction, spanMoleFractionDot, spanAdsorption,
                       spanAdsorptionDot);
  }

  return 0;
}

#else

CVODE::~CVODE() = default;

bool CVODE::propagate(Column&, size_t, Timing&)
{
  throw std::runtime_error("CVODE::propagate() called, but this build was compiled without SUNDIALS support");
}

void CVODE::initialize(Column&) {}

#endif

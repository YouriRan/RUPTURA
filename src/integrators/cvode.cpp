// cvode.cpp
#include "cvode.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <print>
#include <type_traits>

#include "transport.h"

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

bool projectState(Column& column)
{
  bool changed = false;

  for (double& concentration : column.concentration)
  {
    const double projected = std::max(0.0, concentration);
    changed = changed || projected != concentration;
    concentration = projected;
  }

  for (double& loading : column.physisorption)
  {
    const double projected = std::max(0.0, loading);
    changed = changed || projected != loading;
    loading = projected;
  }

  for (double& loading : column.chemisorption)
  {
    const double projected = std::max(0.0, loading);
    changed = changed || projected != loading;
    loading = projected;
  }

  for (double& concentration : column.surfaceConcentration)
  {
    const double projected = std::max(0.0, concentration);
    changed = changed || projected != concentration;
    concentration = projected;
  }

  for (double& concentration : column.poreConcentration)
  {
    const double projected = std::max(0.0, concentration);
    changed = changed || projected != concentration;
    concentration = projected;
  }

  if (column.energyBalance)
  {
    for (double& temperature : column.gasTemperature)
    {
      const double projected = std::max(1e-10, temperature);
      changed = changed || projected != temperature;
      temperature = projected;
    }
    for (double& temperature : column.solidTemperature)
    {
      const double projected = std::max(1e-10, temperature);
      changed = changed || projected != temperature;
      temperature = projected;
    }
    for (double& temperature : column.wallTemperature)
    {
      const double projected = std::max(1e-10, temperature);
      changed = changed || projected != temperature;
      temperature = projected;
    }
  }

  return changed;
}

bool CVODE::propagate(Column& column, size_t step, Timing& timings)
{
  double tNext = currentTime + timeStep;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  if (autoNumberOfSteps)
  {
    double tolerance = 0.0;

    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      if (column.components[j].initialGasMoleFraction <= 0.0) continue;

      const size_t outlet = numberOfGridPoints * numberOfComponents + j;
      const double feed = column.components[j].initialGasMoleFraction;
      tolerance = std::max(tolerance, std::abs((column.moleFraction[outlet] / feed) - 1.0));
    }

    if (tolerance < 0.01)
    {
      std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoNumberOfSteps = false;
    }
  }

  sunrealtype tReturn = currentTime;

  auto timer = timings.scoped(timings.total);

  int flag = CVode(cvodeMem, tNext, stateVector, &tReturn, CV_NORMAL);
  if (flag < 0)
  {
    throw std::runtime_error("CVode failed during propagation");
  }
  currentTime = tReturn;

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
  static_assert(std::is_same_v<sunrealtype, double>, "CVODE currently requires SUNDIALS sunrealtype to be double");

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
  N_VConst(0.0, stateDerivativeVector);

  const ColumnStateLayout layout = column->stateLayout();
  double* stateBase = static_cast<double*>(N_VGetArrayPointer(stateVector));
  double* derivativeBase = static_cast<double*>(N_VGetArrayPointer(stateDerivativeVector));

  auto spanConcentration = layout.concentration(stateBase);
  auto spanPhysisorption = layout.physisorption(stateBase);
  auto spanChemisorption = layout.chemisorption(stateBase);
  const bool usePoreSurfaceTransport = column->surfacePoreTransportEnabled;
  auto spanSurfaceConcentration =
      usePoreSurfaceTransport ? layout.surfaceConcentration(stateBase) : column->surfaceConcentration;
  auto spanPoreConcentration =
      usePoreSurfaceTransport ? layout.poreConcentration(stateBase) : column->poreConcentration;

  auto spanConcentrationDot = layout.concentration(derivativeBase);
  auto spanPhysisorptionDot = layout.physisorption(derivativeBase);
  auto spanChemisorptionDot = layout.chemisorption(derivativeBase);
  auto spanSurfaceConcentrationDot =
      usePoreSurfaceTransport ? layout.surfaceConcentration(derivativeBase) : column->surfaceConcentrationDot;
  auto spanPoreConcentrationDot =
      usePoreSurfaceTransport ? layout.poreConcentration(derivativeBase) : column->poreConcentrationDot;

  auto spanGasTemperature = (column->energyBalance) ? layout.gasTemperature(stateBase) : column->gasTemperature;
  auto spanSolidTemperature = (column->energyBalance) ? layout.solidTemperature(stateBase) : column->solidTemperature;
  auto spanWallTemperature = (column->energyBalance) ? layout.wallTemperature(stateBase) : column->wallTemperature;
  auto spanGasTemperatureDot =
      (column->energyBalance) ? layout.gasTemperature(derivativeBase) : column->gasTemperatureDot;
  auto spanSolidTemperatureDot =
      (column->energyBalance) ? layout.solidTemperature(derivativeBase) : column->solidTemperatureDot;
  auto spanWallTemperatureDot =
      (column->energyBalance) ? layout.wallTemperature(derivativeBase) : column->wallTemperatureDot;

  for (size_t iteration = 0; iteration < 2; ++iteration)
  {
    computeBulkSpeciesSink(
        column->components, column->numberOfGridPoints, column->numberOfComponents,
        column->maxChemisorptionSites, column->voidFraction, column->particleDensity,
        column->particleDiameter, spanConcentration, spanPhysisorptionDot, spanChemisorptionDot,
        spanSurfaceConcentration, column->bulkSpeciesSink);

    updateVelocityAndPressure(
        column->components, column->boundaryCondition, column->numberOfGridPoints, column->numberOfComponents,
        column->inletPressure, column->outletPressure, column->pressureGradient, column->columnLength,
        column->voidFraction, column->particleDensity, column->columnEntranceVelocity, column->dynamicViscosity,
        column->particleDiameter, column->resolution,
        column->interstitialGasVelocity, column->gasDensity, column->totalConcentration, column->totalPressure,
        spanConcentration, column->partialPressure, column->moleFraction, column->bulkSpeciesSink,
        spanGasTemperature);

    computeEquilibriumLoadings(column->mixture, column->numberOfGridPoints, column->numberOfComponents,
                               column->maxIsothermTerms, column->iastPerformance, column->idealGasMolFractions,
                               column->adsorbedMolFractions, column->numberOfMolecules, column->totalPressure,
                               column->equilibriumAdsorption, column->cachedPressure, column->cachedGrandPotential,
                               column->moleFraction, spanGasTemperature);

    computeSorptionDerivatives(
        column->components, column->numberOfGridPoints, column->numberOfComponents,
        column->maxChemisorptionSites, column->externalTemperature, column->voidFraction,
        column->particleDensity, column->particleDiameter, column->equilibriumAdsorption, spanConcentration,
        spanPhysisorption, spanPhysisorptionDot, spanChemisorption, spanChemisorptionDot,
        spanSurfaceConcentration, spanSurfaceConcentrationDot, spanPoreConcentration,
        spanPoreConcentrationDot, spanSolidTemperature, column->bulkSpeciesSink);
  }

  computeMassDerivatives(column->components, column->numberOfGridPoints, column->numberOfComponents,
                         column->resolution, column->interstitialGasVelocity, spanConcentration, spanConcentrationDot,
                         column->bulkSpeciesSink);

  if (column->energyBalance)
  {
    computeEnergyDerivatives(
        column->components, column->numberOfGridPoints, column->numberOfComponents,
        column->maxChemisorptionSites, column->externalTemperature, column->voidFraction,
        column->particleDensity, column->particleDiameter, column->internalDiameter, column->outerDiameter,
        column->wallDensity, column->gasThermalConductivity, column->wallThermalConductivity,
        column->heatTransferGasSolid, column->heatTransferGasWall, column->heatTransferWallExternal,
        column->heatCapacityGas, column->heatCapacitySolid, column->heatCapacityWall, column->resolution,
        column->interstitialGasVelocity, column->gasDensity, column->coeffDiffusion, spanPhysisorptionDot,
        spanChemisorptionDot, spanGasTemperature, spanGasTemperatureDot, spanSolidTemperature,
        spanSolidTemperatureDot, spanWallTemperature, spanWallTemperatureDot);
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

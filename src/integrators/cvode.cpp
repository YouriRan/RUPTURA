// cvode.cpp
#include "cvode.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <print>
#include <type_traits>

#include "column_multibed.h"
#include "compute_multibed.h"
#include "rk3.h"
#include "sorption.h"
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

namespace
{
template <typename ColumnType>
bool breakthroughConverged(const ColumnType& column)
{
  double tolerance = 0.0;
  for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
  {
    const bool gasPhase = static_cast<size_t>(column.fluidPhase) == 0;
    const double feed = gasPhase ? column.components[comp].initialGasMoleFraction
                                 : column.components[comp].inletLiquidConcentration;
    if (feed <= 0.0) continue;

    const size_t outlet = column.numberOfGridPoints * column.numberOfComponents + comp;
    const double outletValue = gasPhase ? column.moleFraction[outlet] : column.concentration[outlet];
    tolerance = std::max(tolerance, std::abs((outletValue / feed) - 1.0));
  }
  return tolerance < 0.01;
}

template <typename ColumnType, typename EvaluationFunction>
bool propagateCVODE(CVODE& integrator, ColumnType& column, size_t step, Timing& timings,
                    EvaluationFunction evaluateDerivatives)
{
  const double tNext = integrator.currentTime + integrator.timeStep;

  if (integrator.autoNumberOfSteps && column.reactions.empty() && breakthroughConverged(column))
  {
    std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
    integrator.numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
    integrator.autoNumberOfSteps = false;
  }

  sunrealtype tReturn = integrator.currentTime;
  auto timer = timings.scoped(timings.total);

  int flag = CVode(integrator.cvodeMem, tNext, integrator.stateVector, &tReturn, CV_NORMAL);
  if (flag < 0)
  {
    throw std::runtime_error("CVode failed during propagation");
  }
  integrator.currentTime = tReturn;
  clampNonnegative(column.state);

  flag = evaluateDerivatives(tReturn, integrator.stateVector, integrator.stateDerivativeVector, &column);
  if (flag != 0)
  {
    throw std::runtime_error("CVODE derivative evaluation failed after propagation");
  }

  if (integrator.autoNumberOfSteps && !column.reactions.empty() && reactionAutoStopReached(column, integrator.timeStep))
  {
    std::print("\nReaction convergence criteria reached, running 10% longer\n\n\n");

    const size_t minimumSteps = std::max<size_t>(step + 1, 1);
    integrator.numberOfSteps =
        std::max<size_t>(static_cast<size_t>(std::ceil(1.1 * static_cast<double>(minimumSteps))), 2);
    integrator.autoNumberOfSteps = false;
  }

  return (!integrator.autoNumberOfSteps && step >= integrator.numberOfSteps - 1);
}

template <typename ColumnType, typename EvaluationFunction>
void initializeCVODE(CVODE& integrator, ColumnType& column, EvaluationFunction evaluateDerivatives)
{
  static_assert(std::is_same_v<sunrealtype, double>, "CVODE currently requires SUNDIALS sunrealtype to be double");

  SUNContext_Create(SUN_COMM_NULL, &integrator.sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &integrator.sunLogger);
  SUNContext_SetLogger(integrator.sunContext, integrator.sunLogger);

  const sunindextype totalSize = static_cast<sunindextype>(column.state.size());
  integrator.stateVector = N_VMake_Serial(totalSize, column.state.data(), integrator.sunContext);
  integrator.stateDerivativeVector = N_VMake_Serial(totalSize, column.stateDot.data(), integrator.sunContext);
  if (!integrator.stateVector || !integrator.stateDerivativeVector)
  {
    throw std::runtime_error("Failed to allocate CVODE state vectors");
  }

  integrator.cvodeMem = CVodeCreate(CV_BDF, integrator.sunContext);
  CVodeSetMaxNumSteps(integrator.cvodeMem, 1000000);
  CVodeSetUserData(integrator.cvodeMem, &column);
  integrator.currentTime = 0.0;

  int flag = CVodeInit(integrator.cvodeMem, evaluateDerivatives, integrator.currentTime, integrator.stateVector);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeInit failed");

  flag = CVodeSStolerances(integrator.cvodeMem, integrator.relativeTolerance, integrator.absoluteTolerance);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSStolerances failed");

  integrator.solver = SUNNonlinSol_Newton(integrator.stateVector, integrator.sunContext);
  flag = CVodeSetNonlinearSolver(integrator.cvodeMem, integrator.solver);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetNonlinearSolver failed");

  integrator.linearMatrix = SUNDenseMatrix(totalSize, totalSize, integrator.sunContext);
  integrator.linSolver = SUNLinSol_Dense(integrator.stateVector, integrator.linearMatrix, integrator.sunContext);
  flag = CVodeSetLinearSolver(integrator.cvodeMem, integrator.linSolver, integrator.linearMatrix);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");

  CVodeSetJacFn(integrator.cvodeMem, nullptr);
}
}  // namespace

bool CVODE::propagate(Column& column, size_t step, Timing& timings)
{
  return propagateCVODE(*this, column, step, timings, CVODE::evaluateDerivatives);
}

bool CVODE::propagate(MultibedColumn& column, size_t step, Timing& timings)
{
  return propagateCVODE(*this, column, step, timings, CVODE::evaluateMultibedDerivatives);
}

void CVODE::initialize(Column& column) { initializeCVODE(*this, column, CVODE::evaluateDerivatives); }

void CVODE::initialize(MultibedColumn& column) { initializeCVODE(*this, column, CVODE::evaluateMultibedDerivatives); }

void CVODE::reinitialize()
{
  const int flag = CVodeReInit(cvodeMem, currentTime, stateVector);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeReInit failed");
}

/*
 * CVODE evaluates intermediate trial states in storage owned by SUNDIALS. The
 * callbacks below therefore bind layout views to the supplied N_Vector rather
 * than using the column's accepted-state spans directly.
 */
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

  std::fill(column->reactionPhysisorptionSource.begin(), column->reactionPhysisorptionSource.end(), 0.0);
  std::fill(column->reactionChemisorptionSource.begin(), column->reactionChemisorptionSource.end(), 0.0);
  std::fill(column->reactionPoreConcentrationSource.begin(), column->reactionPoreConcentrationSource.end(), 0.0);
  std::fill(column->reactionHeat.begin(), column->reactionHeat.end(), 0.0);
  for (size_t iteration = 0; iteration < 2; ++iteration)
  {
    computeBulkSpeciesSink(column->components, column->numberOfGridPoints, column->numberOfComponents,
                           column->maxChemisorptionSites, column->geometry, column->particleDensity, spanConcentration,
                           spanPhysisorptionDot, spanChemisorptionDot, spanSurfaceConcentration,
                           column->bulkSpeciesSink, column->reactionPhysisorptionSource,
                           column->reactionChemisorptionSource);

    updateVelocityAndPressure(
        column->components, column->boundaryCondition, column->numberOfGridPoints, column->numberOfComponents,
        column->inletPressure, column->outletPressure, column->pressureGradient, column->columnLength, column->geometry,
        column->columnEntranceVelocity, column->dynamicViscosity, column->resolution, column->interstitialGasVelocity,
        column->gasDensity, column->totalConcentration, column->totalPressure, spanConcentration,
        column->partialPressure, column->moleFraction, column->bulkSpeciesSink, spanGasTemperature,
        column->fluidPhase, column->liquidDensity, column->pHMode, column->pHValue, column->pKw,
        column->pHComponent, column->pH);

    computePhysisorptionEquilibriumLoadings(
        column->physisorptionMixture, column->numberOfGridPoints, column->numberOfComponents, column->maxIsothermTerms,
        column->iastPerformance, column->idealGasMolFractions, column->adsorbedMolFractions, column->numberOfMolecules,
        column->totalPressure, column->equilibriumPhysisorption, column->cachedPressure, column->cachedGrandPotential,
        column->moleFraction, spanGasTemperature,
        column->fluidPhase == Column::FluidPhase::Gas
            ? MixturePrediction::DrivingForceInput::MoleFraction
            : MixturePrediction::DrivingForceInput::Concentration,
        spanConcentration, column->pH);

    computeChemisorptionEquilibriumLoadings(
        column->chemisorptionMixture, column->numberOfGridPoints, column->numberOfComponents,
        column->maxChemisorptionSites, column->iastPerformance, column->idealGasMolFractions,
        column->adsorbedMolFractions, column->numberOfMolecules, column->totalPressure,
        column->equilibriumChemisorption, column->cachedChemisorptionPressure,
        column->cachedChemisorptionGrandPotential, column->moleFraction, spanGasTemperature,
        column->fluidPhase == Column::FluidPhase::Gas
            ? MixturePrediction::DrivingForceInput::MoleFraction
            : MixturePrediction::DrivingForceInput::Concentration,
        spanConcentration, column->pH);

    std::fill(spanPhysisorptionDot.begin(), spanPhysisorptionDot.end(), 0.0);
    std::fill(spanChemisorptionDot.begin(), spanChemisorptionDot.end(), 0.0);
    std::fill(spanSurfaceConcentrationDot.begin(), spanSurfaceConcentrationDot.end(), 0.0);
    std::fill(spanPoreConcentrationDot.begin(), spanPoreConcentrationDot.end(), 0.0);

    computePhysisorption(column->components, column->numberOfGridPoints, column->numberOfComponents,
                         column->equilibriumPhysisorption, spanPhysisorption, spanPhysisorptionDot);
    computeChemisorption(column->components, column->numberOfGridPoints, column->numberOfComponents,
                         column->maxChemisorptionSites, column->externalTemperature, column->geometry,
                         column->particleDensity, column->equilibriumChemisorption, spanConcentration,
                         spanChemisorption, spanChemisorptionDot, spanPoreConcentration, spanSolidTemperature);
    computeChemisorptionTransportDerivatives(
        column->components, column->numberOfGridPoints, column->numberOfComponents, column->maxChemisorptionSites,
        column->geometry, column->particleDensity, spanConcentration, spanChemisorptionDot, spanSurfaceConcentration,
        spanSurfaceConcentrationDot, spanPoreConcentration, spanPoreConcentrationDot);
    computeBulkSpeciesSink(column->components, column->numberOfGridPoints, column->numberOfComponents,
                           column->maxChemisorptionSites, column->geometry, column->particleDensity, spanConcentration,
                           spanPhysisorptionDot, spanChemisorptionDot, spanSurfaceConcentration,
                           column->bulkSpeciesSink);
    computeReactionDerivatives(column->components, column->reactions, column->numberOfGridPoints,
                               column->numberOfComponents, column->maxChemisorptionSites, column->externalTemperature,
                               spanPhysisorption, spanPhysisorptionDot, spanChemisorption, spanChemisorptionDot,
                               spanPoreConcentration, spanPoreConcentrationDot, spanSolidTemperature,
                               column->reactionPhysisorptionSource, column->reactionChemisorptionSource,
                               column->reactionPoreConcentrationSource, column->reactionHeat);
  }

  computeMassDerivatives(column->components, column->numberOfGridPoints, column->numberOfComponents, column->resolution,
                         column->interstitialGasVelocity, spanConcentration, spanConcentrationDot,
                         column->bulkSpeciesSink);

  if (column->energyBalance)
  {
    computeEnergyDerivatives(column->components, column->numberOfGridPoints, column->numberOfComponents,
                             column->maxChemisorptionSites, column->externalTemperature, column->geometry,
                             column->particleDensity, column->wallDensity, column->gasThermalConductivity,
                             column->wallThermalConductivity, column->heatTransferGasSolid, column->heatTransferGasWall,
                             column->heatTransferWallExternal, column->heatCapacityGas, column->heatCapacitySolid,
                             column->heatCapacityWall, column->resolution, column->interstitialGasVelocity,
                             column->gasDensity, column->coeffDiffusion, spanPhysisorptionDot, spanChemisorptionDot,
                             spanGasTemperature, spanGasTemperatureDot, spanSolidTemperature, spanSolidTemperatureDot,
                             spanWallTemperature, spanWallTemperatureDot, column->reactionPhysisorptionSource,
                             column->reactionChemisorptionSource, column->reactionHeat);
  }

  return 0;
}

int CVODE::evaluateMultibedDerivatives(sunrealtype /*t*/, N_Vector stateVector, N_Vector stateDerivativeVector,
                                       void* user_data)
{
  auto* column = static_cast<MultibedColumn*>(user_data);
  N_VConst(0.0, stateDerivativeVector);

  const MultibedColumnStateLayout layout = column->stateLayout();
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

  auto spanGasTemperature = column->energyBalance ? layout.gasTemperature(stateBase) : column->gasTemperature;
  auto spanSolidTemperature = column->energyBalance ? layout.solidTemperature(stateBase) : column->solidTemperature;
  auto spanWallTemperature = column->energyBalance ? layout.wallTemperature(stateBase) : column->wallTemperature;
  auto spanGasTemperatureDot =
      column->energyBalance ? layout.gasTemperature(derivativeBase) : column->gasTemperatureDot;
  auto spanSolidTemperatureDot =
      column->energyBalance ? layout.solidTemperature(derivativeBase) : column->solidTemperatureDot;
  auto spanWallTemperatureDot =
      column->energyBalance ? layout.wallTemperature(derivativeBase) : column->wallTemperatureDot;

  std::fill(column->reactionPhysisorptionSource.begin(), column->reactionPhysisorptionSource.end(), 0.0);
  std::fill(column->reactionChemisorptionSource.begin(), column->reactionChemisorptionSource.end(), 0.0);
  std::fill(column->reactionPoreConcentrationSource.begin(), column->reactionPoreConcentrationSource.end(), 0.0);
  std::fill(column->reactionHeat.begin(), column->reactionHeat.end(), 0.0);

  for (size_t iteration = 0; iteration < 2; ++iteration)
  {
    computeBulkSpeciesSink(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                           column->numberOfAdsorbents, column->maxChemisorptionSites, column->geometries,
                           column->adsorbentVoidFractions,
                           column->particleDensities, column->particleDiameters, column->fractionOfAdsorbent,
                           column->totalVoidFraction, spanConcentration, spanPhysisorptionDot, spanChemisorptionDot,
                           spanSurfaceConcentration, column->bulkSpeciesSink, column->reactionPhysisorptionSource,
                           column->reactionChemisorptionSource);

    updateVelocityAndPressure(
        column->components, column->boundaryCondition, column->numberOfGridPoints, column->numberOfComponents,
        column->inletPressure, column->outletPressure, column->pressureGradient, column->columnLength,
        column->numberOfAdsorbents, column->columnEntranceVelocity, column->dynamicViscosity, column->columnDistances,
        column->fractionOfAdsorbent, column->geometries, column->adsorbentScaledVoidFraction,
        column->totalVoidFraction,
        column->interstitialGasVelocity,
        column->gasDensity, column->totalConcentration, column->totalPressure, spanConcentration,
        column->partialPressure, column->moleFraction, column->bulkSpeciesSink, spanGasTemperature,
        column->fluidPhase, column->liquidDensity, column->pHMode, column->pHValue, column->pKw,
        column->pHComponent, column->pH);

    computePhysisorptionEquilibriumLoadings(
        column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->fractionOfAdsorbent, column->hasAdsorbentOfType, column->maxIsothermTerms,
        column->iastPerformance, column->idealGasMolFractions, column->adsorbedMolFractions, column->numberOfMolecules,
        column->totalPressure, column->equilibriumPhysisorption, column->cachedPressure, column->cachedGrandPotential,
        column->moleFraction, spanGasTemperature,
        column->fluidPhase == MultibedColumn::FluidPhase::Gas
            ? MixturePrediction::DrivingForceInput::MoleFraction
            : MixturePrediction::DrivingForceInput::Concentration,
        spanConcentration, column->pH);

    computeChemisorptionEquilibriumLoadings(
        column->chemisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->fractionOfAdsorbent, column->hasAdsorbentOfType,
        column->maxChemisorptionSites, column->iastPerformance, column->idealGasMolFractions,
        column->adsorbedMolFractions, column->numberOfMolecules, column->totalPressure,
        column->equilibriumChemisorption, column->cachedChemisorptionPressure,
        column->cachedChemisorptionGrandPotential, column->moleFraction, spanGasTemperature,
        column->fluidPhase == MultibedColumn::FluidPhase::Gas
            ? MixturePrediction::DrivingForceInput::MoleFraction
            : MixturePrediction::DrivingForceInput::Concentration,
        spanConcentration, column->pH);

    computePhysisorption(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                         column->numberOfAdsorbents, column->fractionOfAdsorbent, column->equilibriumPhysisorption,
                         spanPhysisorption, spanPhysisorptionDot);
    computeChemisorption(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                         column->numberOfAdsorbents, column->maxChemisorptionSites, column->externalTemperature,
                         column->fractionOfAdsorbent, column->adsorbentVoidFractions, column->particleDensities,
                         column->equilibriumChemisorption, spanConcentration, spanChemisorption, spanChemisorptionDot,
                         spanPoreConcentration, spanSolidTemperature);
    computeChemisorptionTransportDerivatives(
        column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->maxChemisorptionSites, column->fractionOfAdsorbent,
        column->geometries, column->adsorbentVoidFractions, column->particleDensities,
        column->particleDiameters, column->totalVoidFraction,
        spanConcentration, spanChemisorptionDot, spanSurfaceConcentration, spanSurfaceConcentrationDot,
        spanPoreConcentration, spanPoreConcentrationDot);

    computeReactionDerivatives(
        column->physisorptionMixtures.front().components, column->reactions, column->numberOfGridPoints,
        column->numberOfComponents, column->maxChemisorptionSites, column->externalTemperature, spanPhysisorption,
        spanPhysisorptionDot, spanChemisorption, spanChemisorptionDot, spanPoreConcentration, spanPoreConcentrationDot,
        spanSolidTemperature, column->reactionPhysisorptionSource, column->reactionChemisorptionSource,
        column->reactionPoreConcentrationSource, column->reactionHeat);

    computeBulkSpeciesSink(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                           column->numberOfAdsorbents, column->maxChemisorptionSites, column->geometries,
                           column->adsorbentVoidFractions,
                           column->particleDensities, column->particleDiameters, column->fractionOfAdsorbent,
                           column->totalVoidFraction, spanConcentration, spanPhysisorptionDot, spanChemisorptionDot,
                           spanSurfaceConcentration, column->bulkSpeciesSink, column->reactionPhysisorptionSource,
                           column->reactionChemisorptionSource);
  }

  computeMassDerivatives(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                         column->numberOfAdsorbents, column->columnDistances, column->fractionOfAdsorbent,
                         column->interstitialGasVelocity, spanConcentration, spanConcentrationDot,
                         column->bulkSpeciesSink, column->fluidPhase, column->totalVoidFraction);

  if (column->energyBalance)
  {
    computeEnergyDerivatives(
        column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->externalTemperature, column->totalVoidFraction, column->geometries,
        column->particleDensities,
        column->particleDiameters, column->fractionOfAdsorbent, column->internalDiameter, column->outerDiameter,
        column->wallDensity, column->gasThermalConductivity, column->wallThermalConductivity,
        column->heatTransferGasSolid, column->heatTransferGasWall, column->heatTransferWallExternal,
        column->heatCapacityGas, column->heatCapacitySolid, column->heatCapacityWall, column->columnDistances,
        column->interstitialGasVelocity, column->gasDensity, column->coeffDiffusion, column->maxChemisorptionSites,
        spanPhysisorptionDot, spanChemisorptionDot, spanGasTemperature, spanGasTemperatureDot, spanSolidTemperature,
        spanSolidTemperatureDot, spanWallTemperature, spanWallTemperatureDot, column->reactionPhysisorptionSource,
        column->reactionChemisorptionSource, column->reactionHeat);
  }

  return 0;
}

#else

CVODE::~CVODE() = default;

bool CVODE::propagate(Column&, size_t, Timing&)
{
  throw std::runtime_error("CVODE::propagate() called, but this build was compiled without SUNDIALS support");
}

bool CVODE::propagate(MultibedColumn&, size_t, Timing&)
{
  throw std::runtime_error("CVODE::propagate() called, but this build was compiled without SUNDIALS support");
}

void CVODE::initialize(Column&) {}

void CVODE::initialize(MultibedColumn&) {}

void CVODE::reinitialize() {}

#endif

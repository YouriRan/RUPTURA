#include "rk3.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mdspan>
#include <print>
#include <type_traits>
#include <vector>

#include "column_multibed.h"
#include "compute.h"
#include "compute_multibed.h"
#include "sorption.h"
#include "transport.h"

namespace
{
bool reactionStepChangeSmall(std::span<const double> state, std::span<const double> derivative,
                             double timeStep) noexcept
{
  if (state.empty() || derivative.empty()) return true;
  if (state.size() != derivative.size()) return false;

  constexpr double absoluteTolerance = 1.0e-8;
  constexpr double relativeTolerance = 1.0e-5;
  const double dt = std::max(0.0, timeStep);

  for (size_t i = 0; i < state.size(); ++i)
  {
    const double scale = std::max(1.0, std::abs(state[i]));
    const double allowedChange = absoluteTolerance + relativeTolerance * scale;
    if (std::abs(derivative[i] * dt) > allowedChange) return false;
  }
  return true;
}
}  // namespace

bool reactionAutoStopReached(const Column& column, double timeStep) noexcept
{
  if (column.reactions.empty()) return false;

  return reactionStepChangeSmall(column.physisorption, column.physisorptionDot, timeStep) &&
         reactionStepChangeSmall(column.chemisorption, column.chemisorptionDot, timeStep) &&
         reactionStepChangeSmall(column.surfaceConcentration, column.surfaceConcentrationDot, timeStep) &&
         reactionStepChangeSmall(column.poreConcentration, column.poreConcentrationDot, timeStep) &&
         reactionStepChangeSmall(column.physisorption, column.reactionPhysisorptionSource, timeStep) &&
         reactionStepChangeSmall(column.chemisorption, column.reactionChemisorptionSource, timeStep) &&
         reactionStepChangeSmall(column.poreConcentration, column.reactionPoreConcentrationSource, timeStep);
}

void precompute(Column& column, Timing& timings)
{
  timings.measure(timings.updateVelocityAndPressure,
                  [&]
                  {
                    ::computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                                             column.maxChemisorptionSites, column.geometry, column.particleDensity,
                                             column.concentration, column.physisorptionDot, column.chemisorptionDot,
                                             column.surfaceConcentration, column.bulkSpeciesSink,
                                             column.reactionPhysisorptionSource, column.reactionChemisorptionSource);
                    ::updateVelocityAndPressure(
                        column.components, column.boundaryCondition, column.numberOfGridPoints,
                        column.numberOfComponents, column.inletPressure, column.outletPressure, column.pressureGradient,
                        column.columnLength, column.geometry, column.columnEntranceVelocity, column.dynamicViscosity,
                        column.resolution, column.interstitialGasVelocity, column.gasDensity, column.totalConcentration,
                        column.totalPressure, column.concentration, column.partialPressure, column.moleFraction,
                        column.bulkSpeciesSink, column.gasTemperature, column.fluidPhase, column.liquidDensity,
                        column.pHMode, column.pHValue, column.pKw, column.pHComponent, column.pH);
                  });
  timings.measure(
      timings.computeEquilibriumLoadings,
      [&]
      {
        ::computePhysisorptionEquilibriumLoadings(
            column.physisorptionMixture, column.numberOfGridPoints, column.numberOfComponents, column.maxIsothermTerms,
            column.iastPerformance, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
            column.totalPressure, column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential,
            column.moleFraction, column.gasTemperature,
            column.fluidPhase == Column::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                         : MixturePrediction::DrivingForceInput::Concentration,
            column.concentration, column.pH);
        ::computeChemisorptionEquilibriumLoadings(
            column.chemisorptionMixture, column.numberOfGridPoints, column.numberOfComponents,
            column.maxChemisorptionSites, column.iastPerformance, column.idealGasMolFractions,
            column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
            column.equilibriumChemisorption, column.cachedChemisorptionPressure,
            column.cachedChemisorptionGrandPotential, column.moleFraction, column.gasTemperature,
            column.fluidPhase == Column::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                         : MixturePrediction::DrivingForceInput::Concentration,
            column.concentration, column.pH);
      });
}

void precompute(MultibedColumn& column, Timing& timings)
{
  timings.measure(
      timings.updateVelocityAndPressure,
      [&]
      {
        ::computeBulkSpeciesSink(
            column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
            column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries, column.adsorbentVoidFractions,
            column.particleDensities, column.particleDiameters, column.fractionOfAdsorbent, column.totalVoidFraction,
            column.concentration, column.physisorptionDot, column.chemisorptionDot, column.surfaceConcentration,
            column.bulkSpeciesSink, column.reactionPhysisorptionSource, column.reactionChemisorptionSource);
        ::updateVelocityAndPressure(
            column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
            column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength,
            column.numberOfAdsorbents, column.columnEntranceVelocity, column.dynamicViscosity, column.columnDistances,
            column.fractionOfAdsorbent, column.geometries, column.adsorbentScaledVoidFraction, column.totalVoidFraction,
            column.interstitialGasVelocity, column.gasDensity, column.totalConcentration, column.totalPressure,
            column.concentration, column.partialPressure, column.moleFraction, column.bulkSpeciesSink,
            column.gasTemperature, column.fluidPhase, column.liquidDensity, column.pHMode, column.pHValue, column.pKw,
            column.pHComponent, column.pH);
      });
  timings.measure(
      timings.computeEquilibriumLoadings,
      [&]
      {
        ::computePhysisorptionEquilibriumLoadings(
            column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
            column.numberOfAdsorbents, column.fractionOfAdsorbent, column.hasAdsorbentOfType, column.maxIsothermTerms,
            column.iastPerformance, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
            column.totalPressure, column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential,
            column.moleFraction, column.gasTemperature,
            column.fluidPhase == MultibedColumn::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                                 : MixturePrediction::DrivingForceInput::Concentration,
            column.concentration, column.pH);
        ::computeChemisorptionEquilibriumLoadings(
            column.chemisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
            column.numberOfAdsorbents, column.fractionOfAdsorbent, column.hasAdsorbentOfType,
            column.maxChemisorptionSites, column.iastPerformance, column.idealGasMolFractions,
            column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
            column.equilibriumChemisorption, column.cachedChemisorptionPressure,
            column.cachedChemisorptionGrandPotential, column.moleFraction, column.gasTemperature,
            column.fluidPhase == MultibedColumn::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                                 : MixturePrediction::DrivingForceInput::Concentration,
            column.concentration, column.pH);
      });
}

void computeDerivatives(Column& column, double elapsedTime)
{
  std::fill(column.physisorptionDot.begin(), column.physisorptionDot.end(), 0.0);
  std::fill(column.chemisorptionDot.begin(), column.chemisorptionDot.end(), 0.0);
  std::fill(column.surfaceConcentrationDot.begin(), column.surfaceConcentrationDot.end(), 0.0);
  std::fill(column.poreConcentrationDot.begin(), column.poreConcentrationDot.end(), 0.0);

  ::computePhysisorption(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.equilibriumPhysisorption, column.physisorption, column.physisorptionDot);

  ::computeChemisorption(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.maxChemisorptionSites, column.externalTemperature, column.geometry,
                         column.particleDensity, elapsedTime, column.equilibriumChemisorption, column.concentration,
                         column.chemisorption, column.chemisorptionDot, column.poreConcentration,
                         column.solidTemperature);

  ::computeChemisorptionTransportDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents,
                                             column.maxChemisorptionSites, column.geometry, column.particleDensity,
                                             column.concentration, column.chemisorptionDot, column.surfaceConcentration,
                                             column.surfaceConcentrationDot, column.poreConcentration,
                                             column.poreConcentrationDot);

  ::computeReactionDerivatives(
      column.components, column.reactions, column.numberOfGridPoints, column.numberOfComponents,
      column.maxChemisorptionSites, column.externalTemperature, column.physisorption, column.physisorptionDot,
      column.chemisorption, column.chemisorptionDot, column.poreConcentration, column.poreConcentrationDot,
      column.solidTemperature, column.reactionPhysisorptionSource, column.reactionChemisorptionSource,
      column.reactionPoreConcentrationSource, column.reactionHeat);

  ::computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                           column.maxChemisorptionSites, column.geometry, column.particleDensity, column.concentration,
                           column.physisorptionDot, column.chemisorptionDot, column.surfaceConcentration,
                           column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                           column.reactionChemisorptionSource);

  ::computeMassDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents, column.resolution,
                           column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                           column.bulkSpeciesSink);

  if (column.energyBalance)
  {
    ::computeEnergyDerivatives(
        column.components, column.numberOfGridPoints, column.numberOfComponents, column.maxChemisorptionSites,
        column.externalTemperature, column.geometry, column.particleDensity, column.wallDensity,
        column.gasThermalConductivity, column.wallThermalConductivity, column.heatTransferGasSolid,
        column.heatTransferGasWall, column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
        column.heatCapacityWall, column.resolution, column.interstitialGasVelocity, column.gasDensity,
        column.coeffDiffusion, column.physisorptionDot, column.chemisorptionDot, column.gasTemperature,
        column.gasTemperatureDot, column.solidTemperature, column.solidTemperatureDot, column.wallTemperature,
        column.wallTemperatureDot, column.reactionPhysisorptionSource, column.reactionChemisorptionSource,
        column.reactionHeat);
  }
}

void computeDerivatives(MultibedColumn& column, double elapsedTime)
{
  ::computePhysisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.fractionOfAdsorbent, column.equilibriumPhysisorption,
                         column.physisorption, column.physisorptionDot);

  ::computeChemisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.maxChemisorptionSites, column.externalTemperature,
                         elapsedTime, column.fractionOfAdsorbent, column.adsorbentVoidFractions, column.particleDensities,
                         column.equilibriumChemisorption, column.concentration, column.chemisorption,
                         column.chemisorptionDot, column.poreConcentration, column.solidTemperature);

  ::computeChemisorptionTransportDerivatives(
      column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.maxChemisorptionSites, column.fractionOfAdsorbent, column.geometries,
      column.adsorbentVoidFractions, column.particleDensities,
      column.particleDiameters, column.totalVoidFraction, column.concentration, column.chemisorptionDot,
      column.surfaceConcentration, column.surfaceConcentrationDot, column.poreConcentration,
      column.poreConcentrationDot);

  ::computeReactionDerivatives(
      column.physisorptionMixtures.front().components, column.reactions, column.numberOfGridPoints,
      column.numberOfComponents, column.maxChemisorptionSites, column.externalTemperature, column.physisorption,
      column.physisorptionDot, column.chemisorption, column.chemisorptionDot, column.poreConcentration,
      column.poreConcentrationDot, column.solidTemperature, column.reactionPhysisorptionSource,
      column.reactionChemisorptionSource, column.reactionPoreConcentrationSource, column.reactionHeat);

  ::computeBulkSpeciesSink(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                           column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries,
                           column.adsorbentVoidFractions,
                           column.particleDensities, column.particleDiameters, column.fractionOfAdsorbent,
                           column.totalVoidFraction, column.concentration, column.physisorptionDot,
                           column.chemisorptionDot, column.surfaceConcentration, column.bulkSpeciesSink,
                           column.reactionPhysisorptionSource, column.reactionChemisorptionSource);

  ::computeMassDerivatives(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                           column.numberOfAdsorbents, column.columnDistances, column.fractionOfAdsorbent,
                           column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                           column.bulkSpeciesSink, column.fluidPhase, column.totalVoidFraction);

  if (column.energyBalance)
  {
    ::computeEnergyDerivatives(
        column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
        column.externalTemperature, column.totalVoidFraction, column.geometries, column.particleDensities,
        column.particleDiameters,
        column.fractionOfAdsorbent, column.internalDiameter, column.outerDiameter, column.wallDensity,
        column.gasThermalConductivity, column.wallThermalConductivity, column.heatTransferGasSolid,
        column.heatTransferGasWall, column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
        column.heatCapacityWall, column.columnDistances, column.interstitialGasVelocity, column.gasDensity,
        column.coeffDiffusion, column.maxChemisorptionSites, column.physisorptionDot, column.chemisorptionDot,
        column.gasTemperature, column.gasTemperatureDot, column.solidTemperature, column.solidTemperatureDot,
        column.wallTemperature, column.wallTemperatureDot, column.reactionPhysisorptionSource,
        column.reactionChemisorptionSource, column.reactionHeat);
  }
}

bool reactionAutoStopReached(const MultibedColumn& column, double timeStep) noexcept
{
  if (column.reactions.empty()) return false;

  return reactionStepChangeSmall(column.physisorption, column.physisorptionDot, timeStep) &&
         reactionStepChangeSmall(column.chemisorption, column.chemisorptionDot, timeStep) &&
         reactionStepChangeSmall(column.surfaceConcentration, column.surfaceConcentrationDot, timeStep) &&
         reactionStepChangeSmall(column.poreConcentration, column.poreConcentrationDot, timeStep) &&
         reactionStepChangeSmall(column.physisorption, column.reactionPhysisorptionSource, timeStep) &&
         reactionStepChangeSmall(column.chemisorption, column.reactionChemisorptionSource, timeStep) &&
         reactionStepChangeSmall(column.poreConcentration, column.reactionPoreConcentrationSource, timeStep);
}

void updateStateRK(Column& column, Column& newColumn, double alpha, double beta, double timeStep)
{
  for (size_t i = 0; i < column.physisorption.size(); i++)
  {
    newColumn.physisorption[i] = alpha * column.physisorption[i] +
                                 beta * (newColumn.physisorption[i] + timeStep * newColumn.physisorptionDot[i]);
    newColumn.concentration[i] =
        std::max(0.0, alpha * column.concentration[i] +
                          beta * (newColumn.concentration[i] + timeStep * newColumn.concentrationDot[i]));
  }

  for (size_t i = 0; i < column.chemisorption.size(); ++i)
  {
    newColumn.chemisorption[i] = alpha * column.chemisorption[i] +
                                 beta * (newColumn.chemisorption[i] + timeStep * newColumn.chemisorptionDot[i]);
  }

  for (size_t i = 0; i < column.surfaceConcentration.size(); ++i)
  {
    newColumn.surfaceConcentration[i] =
        std::max(0.0, alpha * column.surfaceConcentration[i] +
                          beta * (newColumn.surfaceConcentration[i] + timeStep * newColumn.surfaceConcentrationDot[i]));
    newColumn.poreConcentration[i] =
        std::max(0.0, alpha * column.poreConcentration[i] +
                          beta * (newColumn.poreConcentration[i] + timeStep * newColumn.poreConcentrationDot[i]));
  }
  if (column.energyBalance)
  {
    for (size_t grid = 0; grid < column.numberOfGridPoints + 1; grid++)
    {
      newColumn.gasTemperature[grid] =
          alpha * column.gasTemperature[grid] +
          beta * (newColumn.gasTemperature[grid] + timeStep * newColumn.gasTemperatureDot[grid]);
      newColumn.solidTemperature[grid] =
          alpha * column.solidTemperature[grid] +
          beta * (newColumn.solidTemperature[grid] + timeStep * newColumn.solidTemperatureDot[grid]);
      newColumn.wallTemperature[grid] =
          alpha * column.wallTemperature[grid] +
          beta * (newColumn.wallTemperature[grid] + timeStep * newColumn.wallTemperatureDot[grid]);
    }
  }
}

void updateStateRK(MultibedColumn& column, MultibedColumn& newColumn, double alpha, double beta, double timeStep)
{
  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    newColumn.physisorption[i] = alpha * column.physisorption[i] +
                                 beta * (newColumn.physisorption[i] + timeStep * newColumn.physisorptionDot[i]);
    newColumn.concentration[i] =
        std::max(0.0, alpha * column.concentration[i] +
                          beta * (newColumn.concentration[i] + timeStep * newColumn.concentrationDot[i]));
  }

  for (size_t i = 0; i < column.chemisorption.size(); ++i)
  {
    newColumn.chemisorption[i] = alpha * column.chemisorption[i] +
                                 beta * (newColumn.chemisorption[i] + timeStep * newColumn.chemisorptionDot[i]);
  }

  for (size_t i = 0; i < column.surfaceConcentration.size(); ++i)
  {
    newColumn.surfaceConcentration[i] =
        std::max(0.0, alpha * column.surfaceConcentration[i] +
                          beta * (newColumn.surfaceConcentration[i] + timeStep * newColumn.surfaceConcentrationDot[i]));
    newColumn.poreConcentration[i] =
        std::max(0.0, alpha * column.poreConcentration[i] +
                          beta * (newColumn.poreConcentration[i] + timeStep * newColumn.poreConcentrationDot[i]));
  }

  if (column.energyBalance)
  {
    for (size_t grid = 0; grid < column.numberOfGridPoints + 1; ++grid)
    {
      newColumn.gasTemperature[grid] =
          alpha * column.gasTemperature[grid] +
          beta * (newColumn.gasTemperature[grid] + timeStep * newColumn.gasTemperatureDot[grid]);
      newColumn.solidTemperature[grid] =
          alpha * column.solidTemperature[grid] +
          beta * (newColumn.solidTemperature[grid] + timeStep * newColumn.solidTemperatureDot[grid]);
      newColumn.wallTemperature[grid] =
          alpha * column.wallTemperature[grid] +
          beta * (newColumn.wallTemperature[grid] + timeStep * newColumn.wallTemperatureDot[grid]);
    }
  }
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

}  // namespace

template <typename ColumnType>
bool RungeKutta3::propagate(ColumnType& column, size_t step, Timing& timings)
{
  static_assert(std::is_same_v<ColumnType, Column> || std::is_same_v<ColumnType, MultibedColumn>);

  auto totalTimer = timings.scoped(timings.total);

  if (autoNumberOfSteps && column.reactions.empty() && breakthroughConverged(column))
  {
    // consider 1% as being visibly indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
    numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
    autoNumberOfSteps = false;
  }

  const double startTime = static_cast<double>(step) * timeStep;
  auto evaluateDerivatives = [&](ColumnType& stage, double elapsedTime)
  { timings.measure(timings.computeDerivatives, [&] { computeDerivatives(stage, elapsedTime); }); };

  auto finalizeStage = [&](ColumnType& stage) { precompute(stage, timings); };

  evaluateDerivatives(column, startTime);
  ColumnType newColumn(column);

  updateStateRK(column, newColumn, 0.0, 1.0, timeStep);
  finalizeStage(newColumn);

  evaluateDerivatives(newColumn, startTime + timeStep);
  updateStateRK(column, newColumn, 0.75, 0.25, timeStep);
  finalizeStage(newColumn);

  evaluateDerivatives(newColumn, startTime + 0.5 * timeStep);
  updateStateRK(column, newColumn, 1.0 / 3.0, 2.0 / 3.0, timeStep);
  clampNonnegative(newColumn.state);
  finalizeStage(newColumn);

  column = newColumn;

  if (autoNumberOfSteps && !column.reactions.empty())
  {
    computeDerivatives(column, startTime + timeStep);
    if (reactionAutoStopReached(column, timeStep))
    {
      std::print("\nReaction convergence criteria reached, running 10% longer\n\n\n");

      const size_t minimumSteps = std::max<size_t>(step + 1, 1);
      numberOfSteps = std::max<size_t>(static_cast<size_t>(std::ceil(1.1 * static_cast<double>(minimumSteps))), 2);
      autoNumberOfSteps = false;
    }
  }
  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}

template bool RungeKutta3::propagate(Column& column, size_t step, Timing& timings);
template bool RungeKutta3::propagate(MultibedColumn& column, size_t step, Timing& timings);

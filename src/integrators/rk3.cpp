#include "rk3.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mdspan>
#include <print>
#include <vector>

#include "column_multibed.h"
#include "compute.h"
#include "compute_multibed.h"
#include "sorption.h"
#include "transport.h"

namespace RK3Helpers
{
void computeSorptionDerivatives(Column& column)
{
  ::computeSorptionDerivatives(
      column.components, column.numberOfGridPoints, column.numberOfComponents, column.maxChemisorptionSites,
      column.externalTemperature, column.geometry, column.particleDensity, column.equilibriumPhysisorption,
      column.equilibriumChemisorption, column.concentration, column.physisorption, column.physisorptionDot,
      column.chemisorption, column.chemisorptionDot, column.surfaceConcentration, column.surfaceConcentrationDot,
      column.poreConcentration, column.poreConcentrationDot, column.solidTemperature, column.bulkSpeciesSink);
  computeReactionDerivatives(column);
}

void computePhysisorption(Column& column)
{
  ::computePhysisorption(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.equilibriumPhysisorption, column.physisorption, column.physisorptionDot);
}

void computeChemisorption(Column& column)
{
  ::computeChemisorption(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.maxChemisorptionSites, column.externalTemperature, column.geometry,
                         column.particleDensity, column.equilibriumChemisorption, column.concentration,
                         column.chemisorption, column.chemisorptionDot, column.poreConcentration,
                         column.solidTemperature);
}

void computeChemisorptionTransportDerivatives(Column& column)
{
  ::computeChemisorptionTransportDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents,
                                             column.maxChemisorptionSites, column.geometry, column.particleDensity,
                                             column.concentration, column.chemisorptionDot, column.surfaceConcentration,
                                             column.surfaceConcentrationDot, column.poreConcentration,
                                             column.poreConcentrationDot);
}

void computeBulkSpeciesSink(Column& column)
{
  ::computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                           column.maxChemisorptionSites, column.geometry, column.particleDensity, column.concentration,
                           column.physisorptionDot, column.chemisorptionDot, column.surfaceConcentration,
                           column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                           column.reactionChemisorptionSource);
}

void updateVelocityAndPressure(Column& column)
{
  computeBulkSpeciesSink(column);
  ::updateVelocityAndPressure(
      column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
      column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength, column.geometry,
      column.columnEntranceVelocity, column.dynamicViscosity, column.resolution, column.interstitialGasVelocity,
      column.gasDensity, column.totalConcentration, column.totalPressure, column.concentration, column.partialPressure,
      column.moleFraction, column.bulkSpeciesSink, column.gasTemperature);
}

void computeEquilibriumLoadings(Column& column)
{
  computePhysisorptionEquilibriumLoadings(
      column.physisorptionMixture, column.numberOfGridPoints, column.numberOfComponents, column.maxIsothermTerms,
      column.iastPerformance, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
      column.totalPressure, column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential,
      column.moleFraction, column.gasTemperature);

  computeChemisorptionEquilibriumLoadings(
      column.chemisorptionMixture, column.numberOfGridPoints, column.numberOfComponents, column.maxChemisorptionSites,
      column.iastPerformance, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
      column.totalPressure, column.equilibriumChemisorption, column.cachedChemisorptionPressure,
      column.cachedChemisorptionGrandPotential, column.moleFraction, column.gasTemperature);
}

void computeReactionDerivatives(Column& column)
{
  ::computeReactionDerivatives(
      column.components, column.reactions, column.numberOfGridPoints, column.numberOfComponents,
      column.maxChemisorptionSites, column.externalTemperature, column.physisorption, column.physisorptionDot,
      column.chemisorption, column.chemisorptionDot, column.poreConcentration, column.poreConcentrationDot,
      column.solidTemperature, column.reactionPhysisorptionSource, column.reactionChemisorptionSource,
      column.reactionPoreConcentrationSource, column.reactionHeat);
}

void computeMassDerivatives(Column& column)
{
  computeBulkSpeciesSink(column);
  ::computeMassDerivatives(column.components, column.numberOfGridPoints, column.numberOfComponents, column.resolution,
                           column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                           column.bulkSpeciesSink);
}

void computeEnergyDerivatives(Column& column)
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

void computeDerivatives(Column& column)
{
  computeMassDerivatives(column);
  if (column.energyBalance)
  {
    computeEnergyDerivatives(column);
  }
}

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
}  // namespace RK3Helpers

namespace RK3MultibedHelpers
{
void computeSorptionDerivatives(ColumnMultibed& column) { computePhysisorption(column); }

void computePhysisorption(ColumnMultibed& column)
{
  ::computePhysisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.fractionOfAdsorbent, column.equilibriumPhysisorption,
                         column.physisorption, column.physisorptionDot);
}

void computeBulkSpeciesSink(ColumnMultibed& column)
{
  ::computeBulkSpeciesSink(column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
                           column.adsorbentVoidFractions, column.particleDensities, column.fractionOfAdsorbent,
                           column.totalVoidFraction, column.physisorptionDot, column.bulkSpeciesSink);
}

void updateVelocityAndPressure(ColumnMultibed& column)
{
  computeBulkSpeciesSink(column);
  ::updateVelocityAndPressure(
      column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
      column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength,
      column.numberOfAdsorbents, column.columnEntranceVelocity, column.dynamicViscosity, column.columnDistances,
      column.fractionOfAdsorbent, column.adsorbentScaledVoidFraction, column.interstitialGasVelocity, column.gasDensity,
      column.totalConcentration, column.totalPressure, column.concentration, column.partialPressure,
      column.moleFraction, column.bulkSpeciesSink, column.gasTemperature);
}

void computeEquilibriumLoadings(ColumnMultibed& column)
{
  computePhysisorptionEquilibriumLoadings(
      column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.fractionOfAdsorbent, column.hasAdsorbentOfType, column.maxIsothermTerms, column.iastPerformance,
      column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
      column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential, column.moleFraction,
      column.gasTemperature);
}

void computeMassDerivatives(ColumnMultibed& column)
{
  computeBulkSpeciesSink(column);
  ::computeMassDerivatives(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                           column.numberOfAdsorbents, column.columnDistances, column.fractionOfAdsorbent,
                           column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                           column.bulkSpeciesSink);
}

void computeEnergyDerivatives(ColumnMultibed& column)
{
  ::computeEnergyDerivatives(
      column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.externalTemperature, column.totalVoidFraction, column.particleDensities, column.particleDiameters,
      column.fractionOfAdsorbent, column.internalDiameter, column.outerDiameter, column.wallDensity,
      column.gasThermalConductivity, column.wallThermalConductivity, column.heatTransferGasSolid,
      column.heatTransferGasWall, column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
      column.heatCapacityWall, column.columnDistances, column.interstitialGasVelocity, column.gasDensity,
      column.coeffDiffusion, column.physisorptionDot, column.gasTemperature, column.gasTemperatureDot,
      column.solidTemperature, column.solidTemperatureDot, column.wallTemperature, column.wallTemperatureDot);
}

void computeDerivatives(ColumnMultibed& column)
{
  computeMassDerivatives(column);
  if (column.energyBalance)
  {
    computeEnergyDerivatives(column);
  }
}
}  // namespace RK3MultibedHelpers

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

void updateStateRK(ColumnMultibed& column, ColumnMultibed& newColumn, double alpha, double beta, double timeStep)
{
  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    newColumn.physisorption[i] = alpha * column.physisorption[i] +
                                 beta * (newColumn.physisorption[i] + timeStep * newColumn.physisorptionDot[i]);
    newColumn.concentration[i] =
        std::max(0.0, alpha * column.concentration[i] +
                          beta * (newColumn.concentration[i] + timeStep * newColumn.concentrationDot[i]));
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
    const double feed = column.components[comp].initialGasMoleFraction;
    if (feed <= 0.0) continue;

    const size_t outlet = column.numberOfGridPoints * column.numberOfComponents + comp;
    tolerance = std::max(tolerance, std::abs((column.moleFraction[outlet] / feed) - 1.0));
  }
  return tolerance < 0.01;
}

template <typename ColumnType, typename SorptionFunction, typename DerivativeFunction, typename VelocityFunction,
          typename EquilibriumFunction>
void advanceSSPRK3(ColumnType& column, double timeStep, Timing& timings, SorptionFunction computeSorptionDerivatives,
                   DerivativeFunction computeDerivatives, VelocityFunction updateVelocityAndPressure,
                   EquilibriumFunction computeEquilibriumLoadings)
{
  auto evaluateDerivatives = [&](ColumnType& stage)
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeSorptionDerivatives(stage);
                      computeDerivatives(stage);
                    });
  };

  auto finalizeStage = [&](ColumnType& stage)
  {
    timings.measure(timings.updateVelocityAndPressure, [&] { updateVelocityAndPressure(stage); });
    timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(stage); });
  };

  evaluateDerivatives(column);
  ColumnType newColumn(column);

  updateStateRK(column, newColumn, 0.0, 1.0, timeStep);
  finalizeStage(newColumn);

  evaluateDerivatives(newColumn);
  updateStateRK(column, newColumn, 0.75, 0.25, timeStep);
  finalizeStage(newColumn);

  evaluateDerivatives(newColumn);
  updateStateRK(column, newColumn, 1.0 / 3.0, 2.0 / 3.0, timeStep);
  finalizeStage(newColumn);

  column = newColumn;
}
}  // namespace

bool RungeKutta3::propagate(Column& column, size_t step, Timing& timings)
{
  auto totalTimer = timings.scoped(timings.total);

  if (autoNumberOfSteps && column.reactions.empty() && breakthroughConverged(column))
  {
    // consider 1% as being visibily indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
    numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
    autoNumberOfSteps = false;
  }

  advanceSSPRK3(column, timeStep, timings, RK3Helpers::computeSorptionDerivatives, RK3Helpers::computeDerivatives,
                RK3Helpers::updateVelocityAndPressure, RK3Helpers::computeEquilibriumLoadings);

  if (autoNumberOfSteps && !column.reactions.empty())
  {
    RK3Helpers::computeSorptionDerivatives(column);
    RK3Helpers::computeDerivatives(column);
    if (RK3Helpers::reactionAutoStopReached(column, timeStep))
    {
      std::print("\nReaction convergence criteria reached, running 10% longer\n\n\n");

      const size_t minimumSteps = std::max<size_t>(step + 1, 1);
      numberOfSteps = std::max<size_t>(static_cast<size_t>(std::ceil(1.1 * static_cast<double>(minimumSteps))), 2);
      autoNumberOfSteps = false;
    }
  }
  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}

bool RungeKutta3::propagate(ColumnMultibed& column, size_t step, Timing& timings)
{
  auto totalTimer = timings.scoped(timings.total);

  if (autoNumberOfSteps && breakthroughConverged(column))
  {
    // consider 1% as being visibly indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
    numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
    autoNumberOfSteps = false;
  }

  advanceSSPRK3(column, timeStep, timings, RK3MultibedHelpers::computeSorptionDerivatives,
                RK3MultibedHelpers::computeDerivatives, RK3MultibedHelpers::updateVelocityAndPressure,
                RK3MultibedHelpers::computeEquilibriumLoadings);

  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}

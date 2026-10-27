#include "rk3.h"

#include <algorithm>
#include <cmath>
#include <span>
#include <vector>

#include "compute.h"
#include "sorption.h"
#include "transport.h"

namespace
{
/**
 * \brief True when every state component changes by less than 1e-7 + 1e-4 max(1, |x|) per second.
 *
 * The bound is a rate, so the criterion does not depend on the time step or the output interval.
 */
bool reactionRateSmall(std::span<const double> state, std::span<const double> derivative) noexcept
{
  if (state.empty() || derivative.empty()) return true;
  if (state.size() != derivative.size()) return false;

  constexpr double absoluteRateTolerance = 1.0e-7;  // state units per s
  constexpr double relativeRateTolerance = 1.0e-4;  // 1/s

  for (size_t i = 0; i < state.size(); ++i)
  {
    const double scale = std::max(1.0, std::abs(state[i]));
    const double allowedRate = absoluteRateTolerance + relativeRateTolerance * scale;
    if (std::abs(derivative[i]) > allowedRate) return false;
  }
  return true;
}

void updateStateRK(std::vector<double>& initialState, Column& column, double alpha, double beta,
                   double timeStep)
{
  // The stage is advanced in place; initialState holds the state at the start of the step, with the same layout.
  const auto layout = column.stateLayout();
  const std::span<double> concentration = layout.concentration(initialState.data());
  const std::span<double> physisorption = layout.physisorption(initialState.data());
  const std::span<double> chemisorption = layout.chemisorption(initialState.data());
  const std::span<double> surfaceConcentration = layout.surfaceConcentration(initialState.data());
  const std::span<double> poreConcentration = layout.poreConcentration(initialState.data());
  const std::span<double> gasTemperature = layout.gasTemperature(initialState.data());
  const std::span<double> solidTemperature = layout.solidTemperature(initialState.data());
  const std::span<double> wallTemperature = layout.wallTemperature(initialState.data());

  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    column.physisorption[i] =
        alpha * physisorption[i] + beta * (column.physisorption[i] + timeStep * column.physisorptionDot[i]);
    column.concentration[i] = std::max(
        0.0, alpha * concentration[i] + beta * (column.concentration[i] + timeStep * column.concentrationDot[i]));
  }

  for (size_t i = 0; i < column.chemisorption.size(); ++i)
  {
    column.chemisorption[i] =
        alpha * chemisorption[i] + beta * (column.chemisorption[i] + timeStep * column.chemisorptionDot[i]);
  }

  for (size_t i = 0; i < column.surfaceConcentration.size(); ++i)
  {
    column.surfaceConcentration[i] =
        std::max(0.0, alpha * surfaceConcentration[i] +
                          beta * (column.surfaceConcentration[i] + timeStep * column.surfaceConcentrationDot[i]));
    column.poreConcentration[i] =
        std::max(0.0, alpha * poreConcentration[i] +
                          beta * (column.poreConcentration[i] + timeStep * column.poreConcentrationDot[i]));
  }

  if (column.energyBalance)
  {
    for (size_t grid = 0; grid < column.numberOfGridPoints + 1; ++grid)
    {
      column.gasTemperature[grid] = alpha * gasTemperature[grid] +
                                    beta * (column.gasTemperature[grid] + timeStep * column.gasTemperatureDot[grid]);
      column.solidTemperature[grid] =
          alpha * solidTemperature[grid] +
          beta * (column.solidTemperature[grid] + timeStep * column.solidTemperatureDot[grid]);
      column.wallTemperature[grid] = alpha * wallTemperature[grid] +
                                     beta * (column.wallTemperature[grid] + timeStep * column.wallTemperatureDot[grid]);
    }
  }
}
}  // namespace

bool reactionAutoStopReached(const Column& column) noexcept
{
  if (column.reactions.empty()) return false;

  // A bulk-concentration reaction changes no adsorbed or pore state, so it is converged only once the bulk is steady.
  // The reaction sources are not tested: a catalytic column at steady state keeps reacting.
  const bool bulkReaction =
      std::any_of(column.reactions.begin(), column.reactions.end(), [](const Reaction& reaction)
                  { return reaction.phase == Reaction::Phase::BulkConcentration; });

  return reactionRateSmall(column.physisorption, column.physisorptionDot) &&
         reactionRateSmall(column.chemisorption, column.chemisorptionDot) &&
         reactionRateSmall(column.surfaceConcentration, column.surfaceConcentrationDot) &&
         reactionRateSmall(column.poreConcentration, column.poreConcentrationDot) &&
         (!bulkReaction || reactionRateSmall(column.concentration, column.concentrationDot));
}

void precompute(Column& column, Timing& timings)
{
  timings.measure(timings.updateVelocityAndPressure,
                  [&]
                  {
                    ::computeBulkSpeciesSink(
                        column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                        column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries,
                        column.particleDensities, column.fractionOfAdsorbent, column.totalVoidFraction,
                        column.concentration, column.physisorptionDot, column.chemisorptionDot,
                        column.surfaceConcentration, column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                        column.reactionChemisorptionSource, column.reactionBulkSource);
                    ::updateVelocityAndPressure(
                        column.components, column.boundaryCondition, column.numberOfGridPoints,
                        column.numberOfComponents, column.inletPressure, column.outletPressure, column.pressureGradient,
                        column.columnLength, column.numberOfAdsorbents, column.columnEntranceVelocity,
                        column.dynamicViscosity, column.columnDistances, column.fractionOfAdsorbent, column.geometries,
                        column.totalVoidFraction, column.interstitialGasVelocity, column.gasDensity,
                        column.totalConcentration, column.totalPressure, column.concentration, column.partialPressure,
                        column.moleFraction, column.bulkSpeciesSink, column.gasTemperature, column.fluidPhase,
                        column.liquidDensity, column.pHMode, column.pHValue, column.pKw, column.pHComponent, column.pH);
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
            column.fluidPhase == Column::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                                 : MixturePrediction::DrivingForceInput::Concentration,
            column.concentration, column.pH);
        ::computeChemisorptionEquilibriumLoadings(
            column.chemisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
            column.numberOfAdsorbents, column.fractionOfAdsorbent, column.hasAdsorbentOfType,
            column.maxChemisorptionSites, column.iastPerformance, column.idealGasMolFractions,
            column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
            column.equilibriumChemisorption, column.cachedChemisorptionPressure,
            column.cachedChemisorptionGrandPotential, column.moleFraction, column.gasTemperature,
            column.fluidPhase == Column::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                                 : MixturePrediction::DrivingForceInput::Concentration,
            column.concentration, column.pH);
      });
}

void computeDerivatives(Column& column, double elapsedTime)
{
  ::computePhysisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.fractionOfAdsorbent, column.equilibriumPhysisorption,
                         column.physisorption, column.physisorptionDot);

  ::computeChemisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.maxChemisorptionSites, column.externalTemperature,
                         elapsedTime, column.fractionOfAdsorbent, column.particleDensities,
                         column.equilibriumChemisorption, column.concentration, column.chemisorption,
                         column.chemisorptionDot, column.poreConcentration, column.solidTemperature);

  ::computeChemisorptionTransportDerivatives(
      column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.maxChemisorptionSites, column.fractionOfAdsorbent, column.geometries, column.particleDensities,
      column.concentration, column.chemisorptionDot, column.surfaceConcentration, column.surfaceConcentrationDot,
      column.poreConcentration, column.poreConcentrationDot);

  ::computeReactionDerivatives(
      column.physisorptionMixtures.front().components, column.reactions, column.numberOfGridPoints,
      column.numberOfComponents, column.maxChemisorptionSites, column.externalTemperature, column.physisorption,
      column.physisorptionDot, column.chemisorption, column.chemisorptionDot, column.poreConcentration,
      column.poreConcentrationDot, column.solidTemperature, column.reactionPhysisorptionSource,
      column.reactionChemisorptionSource, column.reactionPoreConcentrationSource, column.reactionHeat,
      column.concentration, column.reactionBulkSource, column.fractionOfAdsorbent, column.numberOfAdsorbents);

  ::computeBulkSpeciesSink(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                           column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries,
                           column.particleDensities, column.fractionOfAdsorbent, column.totalVoidFraction,
                           column.concentration, column.physisorptionDot, column.chemisorptionDot,
                           column.surfaceConcentration, column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                           column.reactionChemisorptionSource, column.reactionBulkSource);

  ::computeMassDerivatives(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                           column.numberOfAdsorbents, column.columnDistances, column.fractionOfAdsorbent,
                           column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                           column.bulkSpeciesSink, column.fluidPhase, column.totalVoidFraction,
                           column.convectionScheme);

  if (column.energyBalance)
  {
    ::computeEnergyDerivatives(
        column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
        column.externalTemperature, column.geometries, column.particleDensities, column.fractionOfAdsorbent,
        column.wallDensity, column.gasThermalConductivity, column.wallThermalConductivity, column.heatTransferGasSolid,
        column.heatTransferGasWall, column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
        column.heatCapacityWall, column.columnDistances, column.interstitialGasVelocity, column.gasDensity,
        column.coeffDiffusion, column.maxChemisorptionSites, column.physisorptionDot, column.chemisorptionDot,
        column.gasTemperature, column.gasTemperatureDot, column.solidTemperature, column.solidTemperatureDot,
        column.wallTemperature, column.wallTemperatureDot, column.reactionPhysisorptionSource,
        column.reactionChemisorptionSource, column.reactionHeat);
  }
}

void RungeKutta3::step(Column& column, double startTime, double dt, Timing& timings)
{
  auto totalTimer = timings.scoped(timings.total);

  // The stages advance in place from a copy of the initial state, not of the whole column (mixture models,
  // components, caches).
  std::vector<double> initialState(column.state);

  auto stage = [&](double alpha, double beta, double elapsedTime, bool finalStage)
  {
    timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column, elapsedTime); });
    updateStateRK(initialState, column, alpha, beta, dt);
    if (finalStage) clampNonnegative(column.state);
    precompute(column, timings);
  };

  stage(0.0, 1.0, startTime, false);
  stage(0.75, 0.25, startTime + dt, false);
  stage(1.0 / 3.0, 2.0 / 3.0, startTime + 0.5 * dt, true);
}

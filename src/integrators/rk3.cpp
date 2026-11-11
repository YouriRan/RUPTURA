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

void precompute(ColumnData& column, const ColumnStateViews& views, Timing& timings)
{
  timings.measure(
      timings.updateVelocityAndPressure,
      [&]
      {
        ::computeBulkSpeciesSink(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                                 column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries,
                                 column.fractionOfAdsorbent, column.loadingPrefactor, views.concentration,
                                 views.physisorptionDot, views.chemisorptionDot, views.surfaceConcentration,
                                 column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                                 column.reactionChemisorptionSource, column.reactionBulkSource);
        ::updateVelocityAndPressure(
            column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
            column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength,
            column.numberOfAdsorbents, column.columnEntranceVelocity, column.dynamicViscosity, column.columnDistances,
            column.fractionOfAdsorbent, column.geometries, column.pressureDropLaw, column.totalVoidFraction,
            column.interstitialGasVelocity, column.gasDensity, column.totalConcentration, column.totalPressure,
            views.concentration, column.partialPressure, column.moleFraction, column.bulkSpeciesSink,
            views.gasTemperature, column.fluidPhase, column.liquidDensity, column.pHMode, column.pHValue, column.pKw,
            column.pHComponent, column.pH, column.molecularDispersionCoefficient, column.dispersionLaw,
            column.dispersionReferenceTemperature, column.dispersionReferencePressure,
            column.axialDispersionCoefficient, column.convectionScheme, column.inletCondition, nullptr,
            column.flowEnds());
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
            column.moleFraction, views.gasTemperature,
            column.fluidPhase == Column::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                                 : MixturePrediction::DrivingForceInput::Concentration,
            views.concentration, column.pH);
        ::computeChemisorptionEquilibriumLoadings(
            column.chemisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
            column.numberOfAdsorbents, column.fractionOfAdsorbent, column.hasAdsorbentOfType,
            column.maxChemisorptionSites, column.iastPerformance, column.idealGasMolFractions,
            column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
            column.equilibriumChemisorption, column.cachedChemisorptionPressure,
            column.cachedChemisorptionGrandPotential, column.moleFraction, views.gasTemperature,
            column.fluidPhase == Column::FluidPhase::Gas ? MixturePrediction::DrivingForceInput::MoleFraction
                                                                 : MixturePrediction::DrivingForceInput::Concentration,
            views.concentration, column.pH);
      });
}

void computeDerivatives(ColumnData& column, const ColumnStateViews& views, bool sorptionOnly)
{
  ::computePhysisorption(column.numberOfGridPoints, column.numberOfComponents, column.massTransferCoefficient,
                         column.equilibriumPhysisorption, views.physisorption, views.physisorptionDot);

  ::computeChemisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.maxChemisorptionSites, column.externalTemperature,
                         column.fractionOfAdsorbent, column.particleDensities,
                         column.equilibriumChemisorption, views.concentration, views.chemisorption,
                         views.chemisorptionDot, views.poreConcentration, views.solidTemperature);

  ::computeChemisorptionTransportDerivatives(
      column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.maxChemisorptionSites, column.fractionOfAdsorbent, column.geometries, column.poreLoadingConversion,
      views.concentration, views.chemisorptionDot, views.surfaceConcentration, views.surfaceConcentrationDot,
      views.poreConcentration, views.poreConcentrationDot);

  ::computeReactionDerivatives(
      column.physisorptionMixtures.front().components, column.reactions, column.numberOfGridPoints,
      column.numberOfComponents, column.maxChemisorptionSites, column.externalTemperature, views.physisorption,
      views.physisorptionDot, views.chemisorption, views.chemisorptionDot, views.poreConcentration,
      views.poreConcentrationDot, views.solidTemperature, column.reactionPhysisorptionSource,
      column.reactionChemisorptionSource, column.reactionPoreConcentrationSource, column.reactionHeat,
      views.concentration, column.reactionBulkSource, column.fractionOfAdsorbent, column.numberOfAdsorbents,
      column.loadingPrefactor, column.reactionCoverage);

  ::computeBulkSpeciesSink(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                           column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries,
                           column.fractionOfAdsorbent, column.loadingPrefactor, views.concentration,
                           views.physisorptionDot, views.chemisorptionDot, views.surfaceConcentration,
                           column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                           column.reactionChemisorptionSource, column.reactionBulkSource);
  if (sorptionOnly) return;

  // precompute built the velocity from the sink of the previous stage; rebuild it from this sink so that the
  // closure and the species balances remove the same uptake and sum(c_i) stays at P/RT (B-BC-2). With the energy
  // balance the closure also carries the thermal expansion of the gas.
  if (column.fluidPhase == Column::FluidPhase::Gas &&
      column.boundaryCondition != Column::BoundaryCondition::FixedVelocity &&
      column.boundaryCondition != Column::BoundaryCondition::FixedPressureInletVelocity)
  {
    const GasTemperatureBalance gasTemperatureBalance{
        column.columnDistances,
        views.gasTemperature,
        views.solidTemperature,
        views.wallTemperature,
        column.gasDensity,
        column.fluidSolidContactArea,
        column.totalVoidFraction,
        column.gasThermalConductivity,
        column.heatTransferGasSolid,
        column.heatCapacityGas,
        column.geometries.front().contactAreas.fluidWallPerFluidVolume * column.heatTransferGasWall /
            column.heatCapacityGas,
        column.inletCondition,
        column.flowEnds()};
    ::updateVelocityAndPressure(
        column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
        column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength,
        column.numberOfAdsorbents, column.columnEntranceVelocity, column.dynamicViscosity, column.columnDistances,
        column.fractionOfAdsorbent, column.geometries, column.pressureDropLaw, column.totalVoidFraction,
        column.interstitialGasVelocity, column.gasDensity, column.totalConcentration, column.totalPressure,
        views.concentration, column.partialPressure, column.moleFraction, column.bulkSpeciesSink,
        views.gasTemperature, column.fluidPhase, column.liquidDensity, column.pHMode, column.pHValue, column.pKw,
        column.pHComponent, column.pH, column.molecularDispersionCoefficient, column.dispersionLaw,
        column.dispersionReferenceTemperature, column.dispersionReferencePressure,
        column.axialDispersionCoefficient, column.convectionScheme, column.inletCondition,
        column.energyBalance ? &gasTemperatureBalance : nullptr, column.flowEnds());
  }

  ::computeMassDerivatives(column.numberOfGridPoints, column.numberOfComponents, column.columnDistances,
                           column.axialDispersionCoefficient, column.interstitialGasVelocity, views.concentration,
                           views.concentrationDot, column.bulkSpeciesSink, column.fluidPhase, column.totalVoidFraction,
                           column.convectionScheme, column.inletCondition, column.wenoScratch, column.flowEnds());

  if (column.energyBalance)
  {
    ::computeEnergyDerivatives(
        column.numberOfGridPoints, column.numberOfComponents, column.externalTemperature, column.geometries,
        column.particleDensity, column.solidFluidContactArea, column.fluidSolidContactArea, column.heatOfAdsorption,
        column.heatOfChemisorption, column.wallDensity, column.gasThermalConductivity, column.wallThermalConductivity,
        column.heatTransferGasSolid, column.heatTransferGasWall, column.heatTransferWallExternal,
        column.heatCapacityGas, column.heatCapacitySolid, column.heatCapacityWall, column.columnDistances,
        column.interstitialGasVelocity, column.totalVoidFraction, column.gasDensity, column.coeffDiffusion,
        column.maxChemisorptionSites,
        views.physisorptionDot, views.chemisorptionDot, views.gasTemperature, views.gasTemperatureDot,
        views.solidTemperature, views.solidTemperatureDot, views.wallTemperature, views.wallTemperatureDot,
        column.reactionPhysisorptionSource, column.reactionChemisorptionSource, column.reactionHeat,
        column.inletCondition, column.flowEnds());
  }
}

void precompute(Column& column, Timing& timings) { precompute(column, column, timings); }

void computeDerivatives(Column& column) { computeDerivatives(column, column); }

void RungeKutta3::step(Column& column, double dt, Timing& timings)
{
  // The stages advance in place from a copy of the initial state, not of the whole column (mixture models,
  // components, caches).
  std::vector<double> initialState(column.state);
  const double startTime = column.pressureRamp.time;

  // nextStageTime: the time of the next stage, or of the step end, in units of dt.
  auto stage = [&](double alpha, double beta, bool finalStage, double nextStageTime)
  {
    timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column); });
    updateStateRK(initialState, column, alpha, beta, dt);
    if (finalStage) clampNonnegative(column.state);
    column.advancePressureRamp(startTime + nextStageTime * dt);
    precompute(column, timings);
  };

  stage(0.0, 1.0, false, 1.0);
  stage(0.75, 0.25, false, 0.5);
  stage(1.0 / 3.0, 2.0 / 3.0, true, 1.0);
}

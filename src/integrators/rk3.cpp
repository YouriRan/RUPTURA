#include "rk3.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>

#include "compute.h"
#include "sorption.h"
#include "transport.h"
#include "utils.h"

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

  // The optional blocks after the temperatures: coolant, radial cells of the particle, stream amounts.
  for (size_t i = layout.temperatureOffset() + 3 * layout.scalarBlockSize(); i < column.state.size(); ++i)
  {
    column.state[i] = alpha * initialState[i] + beta * (column.state[i] + timeStep * column.stateDot[i]);
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
                                 column.reactionChemisorptionSource, column.reactionBulkSource,
                                 column.chemisorptionSites);
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
            column.flowEnds(), column.gasHeatCapacity);
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

  // Homogeneous surface diffusion in a spherical particle of N radial cells of one width, x = r / R: the loading of
  // cell k changes by 3 (D_s / R^2) (F_{k+1/2} - F_{k-1/2}) / (x_{k+1/2}^3 - x_{k-1/2}^3), with the flux
  // F = x^2 dq/dx through the faces between cells, none through the centre, and (q* - q_N) / (1 / (2 N)) through the
  // surface, which is at equilibrium with the gas. The loading of the particle, the physisorption state, takes the
  // sum of these rates over its volume, 3 (D_s / R^2) F_{N+1/2}. D_s / R^2 is k / 15 of the linear driving force,
  // which is this model at long times (Glueckauf).
  const size_t block = views.physisorption.size();
  const double cells = static_cast<double>(column.numberOfShells);
  for (size_t index = column.numberOfComponents; index < block && !views.shellLoading.empty(); ++index)
  {
    const double rate = column.massTransferCoefficient[index] / 5.0;
    double inner = 0.0;
    for (size_t shell = 0; shell < column.numberOfShells; ++shell)
    {
      const double lower = static_cast<double>(shell) / cells;
      const double upper = static_cast<double>(shell + 1) / cells;
      const double loading = views.shellLoading[shell * block + index];
      const double outer = shell + 1 < column.numberOfShells
                               ? upper * upper * cells * (views.shellLoading[(shell + 1) * block + index] - loading)
                               : 2.0 * cells * (column.equilibriumPhysisorption[index] - loading);
      views.shellLoadingDot[shell * block + index] =
          rate * (outer - inner) / (upper * upper * upper - lower * lower * lower);
      inner = outer;
    }
    views.physisorptionDot[index] = rate * inner;
  }

  ::computeChemisorption(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                         column.numberOfAdsorbents, column.maxChemisorptionSites, column.externalTemperature,
                         column.fractionOfAdsorbent, column.particleDensities,
                         column.equilibriumChemisorption, views.concentration, views.chemisorption,
                         views.chemisorptionDot, views.poreConcentration, views.solidTemperature,
                         column.chemisorptionSites);

  ::computeChemisorptionTransportDerivatives(
      column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
      column.maxChemisorptionSites, column.fractionOfAdsorbent, column.geometries, column.poreLoadingConversion,
      views.concentration, views.chemisorptionDot, views.surfaceConcentration, views.surfaceConcentrationDot,
      views.poreConcentration, views.poreConcentrationDot, column.chemisorptionSites);

  ::computeReactionDerivatives(
      column.physisorptionMixtures, column.reactions, column.numberOfGridPoints,
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
                           column.reactionChemisorptionSource, column.reactionBulkSource,
                           column.chemisorptionSites);
  if (sorptionOnly) return;

  // precompute built the velocity from the sink of the previous stage; rebuild it from this sink so that the
  // closure and the species balances remove the same uptake and sum(c_i) stays at P/RT (B-BC-2). With the energy
  // balance the closure also carries the thermal expansion of the gas. The velocity of MomentumBalance does not depend
  // on the sink.
  if (column.fluidPhase == Column::FluidPhase::Gas &&
      column.boundaryCondition != Column::BoundaryCondition::FixedVelocity &&
      column.boundaryCondition != Column::BoundaryCondition::FixedPressureInletVelocity &&
      column.boundaryCondition != Column::BoundaryCondition::MomentumBalance)
  {
    const GasTemperatureBalance gasTemperatureBalance{column.columnDistances,
                                                      views.gasTemperature,
                                                      views.solidTemperature,
                                                      views.wallTemperature,
                                                      column.gasDensity,
                                                      column.fluidSolidContactArea,
                                                      column.totalVoidFraction,
                                                      column.gasThermalConductivity,
                                                      column.heatTransferGasSolid,
                                                      column.heatCapacityGas,
                                                      column.fluidWallContactArea,
                                                      column.heatTransferGasWall,
                                                      column.inletCondition,
                                                      column.flowEnds(),
                                                      column.temperatureRateScratch,
                                                      column.gasHeatCapacity};
    ::updateVelocityAndPressure(
        column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
        column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength,
        column.numberOfAdsorbents, column.columnEntranceVelocity, column.dynamicViscosity, column.columnDistances,
        column.fractionOfAdsorbent, column.geometries, column.pressureDropLaw, column.totalVoidFraction,
        column.interstitialGasVelocity, column.gasDensity, column.totalConcentration, column.totalPressure,
        views.concentration, column.partialPressure, column.moleFraction, column.bulkSpeciesSink, views.gasTemperature,
        column.fluidPhase, column.liquidDensity, column.pHMode, column.pHValue, column.pKw, column.pHComponent,
        column.pH, column.molecularDispersionCoefficient, column.dispersionLaw, column.dispersionReferenceTemperature,
        column.dispersionReferencePressure, column.axialDispersionCoefficient, column.convectionScheme,
        column.inletCondition, column.energyBalance ? &gasTemperatureBalance : nullptr, column.flowEnds(),
        column.gasHeatCapacity);
  }

  ::computeMassDerivatives(column.numberOfGridPoints, column.numberOfComponents, column.columnDistances,
                           column.axialDispersionCoefficient, column.interstitialGasVelocity, views.concentration,
                           views.concentrationDot, column.bulkSpeciesSink, column.fluidPhase, column.totalVoidFraction,
                           column.convectionScheme, column.inletCondition, column.wenoScratch, column.flowEnds(),
                           views.streamAmountDot);

  // After the stream amounts: the heat through the outer wall, summed over the cells below with the energy balance,
  // and the work of an isentropic compression to the exhaust pressure of the gas that leaves through a pressure port
  // below it, n c_p T ((P_exhaust / P)^(R / c_p) - 1) / efficiency with the molar heat capacity of the end cell.
  const size_t numberOfStreams = 2 * column.numberOfComponents;
  for (size_t end = 0; end < 2 && !views.streamAmountDot.empty(); ++end)
  {
    if (end == 0) views.streamAmountDot[numberOfStreams] = views.streamAmountDot[numberOfStreams + 1] = 0.0;
    const Column::Port& port = column.ports[end];
    const size_t cell = end == 0 ? 1 : column.numberOfGridPoints;
    const auto first = views.streamAmountDot.begin() + static_cast<std::ptrdiff_t>(end * column.numberOfComponents);
    const double leaving = -std::reduce(first, first + static_cast<std::ptrdiff_t>(column.numberOfComponents));
    if (column.boundaryCondition != Column::BoundaryCondition::MomentumBalance ||
        port.mode != Column::Port::Mode::Pressure || port.pressure >= column.exhaustPressure || leaving <= 0.0)
    {
      continue;
    }
    const double heatCapacity =
        (column.gasHeatCapacity.empty() ? column.heatCapacityGas : column.gasHeatCapacity[cell]) *
        column.gasDensity[cell] / std::max(1e-30, column.totalConcentration[cell]);
    views.streamAmountDot[numberOfStreams + 1] +=
        leaving * heatCapacity * views.gasTemperature[cell] *
        (std::pow(column.exhaustPressure / port.pressure, R / heatCapacity) - 1.0) / column.compressionEfficiency;
  }

  if (column.energyBalance)
  {
    ::computeEnergyDerivatives(
        column.numberOfGridPoints, column.numberOfComponents, column.externalTemperature, column.geometries,
        column.particleDensity, column.solidFluidContactArea, column.fluidSolidContactArea, column.fluidWallContactArea,
        column.heatOfAdsorption, column.heatOfChemisorption, column.wallDensity, column.gasThermalConductivity,
        column.wallThermalConductivity, column.heatTransferGasSolid, column.heatTransferGasWall,
        column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid, column.heatCapacityWall,
        column.columnDistances, column.interstitialGasVelocity, column.totalVoidFraction, column.gasDensity,
        column.coeffDiffusion, column.maxChemisorptionSites, views.physisorptionDot, views.chemisorptionDot,
        views.gasTemperature, views.gasTemperatureDot, views.solidTemperature, views.solidTemperatureDot,
        views.wallTemperature, views.wallTemperatureDot, column.reactionPhysisorptionSource,
        column.reactionChemisorptionSource, column.reactionHeat, column.inletCondition, column.flowEnds(),
        column.gasHeatCapacity);

    // The outer wall area per unit length is that per wall volume times the wall cross-section; per unit cross-section
    // of the column, as the stream amounts (Column::crossSection).
    const Geometry& geometry = column.geometries.front();
    const double wallPerimeter = geometry.contactAreas.wallOuterPerWallVolume * geometry.dimensions.solidArea;
    const double crossSection =
        geometry.kind == GeometryKind::Monolith ? geometry.dimensions.outerArea : geometry.dimensions.openArea;

    // A heat-transfer fluid along the outer wall, donor cell from its inlet at the temperature beyond the wall: each
    // cell of it exchanges heat with the wall of that cell, which then sees its temperature instead of the one beyond.
    const size_t last = column.numberOfGridPoints;
    const double wallRate = column.heatTransferWallExternal * geometry.contactAreas.wallOuterPerWallVolume /
                            (column.heatCapacityWall * column.wallDensity);
    for (size_t grid = 1; grid <= last && !views.coolantTemperature.empty(); ++grid)
    {
      const bool inlet = column.coolantCounterCurrent ? grid == last : grid == 1;
      const double entering = inlet ? column.externalTemperature
                                    : views.coolantTemperature[column.coolantCounterCurrent ? grid + 1 : grid - 1];
      const double fluid = views.coolantTemperature[grid];
      views.coolantTemperatureDot[grid] =
          (column.coolantMassFlow * column.coolantHeatCapacity * (entering - fluid) /
               gridSpacing(column.columnDistances, grid) -
           column.heatTransferWallExternal * wallPerimeter * (fluid - views.wallTemperature[grid])) /
          column.coolantHeatCapacityPerLength;
      views.wallTemperatureDot[grid] += wallRate * (fluid - column.externalTemperature);
    }
    for (size_t grid = 1; grid <= last && !views.streamAmountDot.empty() && crossSection > 0.0; ++grid)
    {
      const double beyond =
          views.coolantTemperature.empty() ? column.externalTemperature : views.coolantTemperature[grid];
      views.streamAmountDot[numberOfStreams] += gridSpacing(column.columnDistances, grid) *
                                                column.heatTransferWallExternal * wallPerimeter / crossSection *
                                                (beyond - views.wallTemperature[grid]);
    }

    // MomentumBalance: the pressure of a cell is that of its gas, P = R T C with C = sum_j c_j, so that the pressure
    // work (dP/dt) / (rho c_p) of the gas-temperature balance holds the rate itself: with kappa = R C / (rho c_p),
    // dT/dt = (rate without the pressure work + kappa T (dC/dt) / C) / (1 - kappa), the balance at constant volume.
    for (size_t grid = 1;
         column.boundaryCondition == Column::BoundaryCondition::MomentumBalance && grid <= column.numberOfGridPoints;
         ++grid)
    {
      const auto first = views.concentrationDot.begin() + static_cast<std::ptrdiff_t>(grid * column.numberOfComponents);
      const double totalRate = std::reduce(first, first + static_cast<std::ptrdiff_t>(column.numberOfComponents));
      const double kappa = R * column.totalConcentration[grid] /
                           (std::max(1e-10, column.gasDensity[grid]) *
                            (column.gasHeatCapacity.empty() ? column.heatCapacityGas : column.gasHeatCapacity[grid]));
      if (kappa >= 1.0)
      {
        throw std::runtime_error(
            "Error: MomentumBalance with the energy balance needs a gas heat capacity above R "
            "per mole of gas; heatCapacityGas (per kg) is too small for the lightest gas. Give "
            "every component its GasHeatCapacity in J/(mol K)\n");
      }
      views.gasTemperatureDot[grid] =
          (views.gasTemperatureDot[grid] +
           kappa * views.gasTemperature[grid] * totalRate / std::max(1e-30, column.totalConcentration[grid])) /
          (1.0 - kappa);
    }
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
    if (finalStage)
      clampNonnegative(std::span<double>(column.state).first(column.state.size() - column.streamAmount.size()));
    column.advancePressureRamp(startTime + nextStageTime * dt);
    precompute(column, timings);
  };

  stage(0.0, 1.0, false, 1.0);
  stage(0.75, 0.25, false, 0.5);
  stage(1.0 / 3.0, 2.0 / 3.0, true, 1.0);
}

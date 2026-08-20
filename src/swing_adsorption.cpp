#include "swing_adsorption.h"

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <format>
#include <fstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "compute.h"
#include "compute_multibed.h"
#include "rk3.h"
#include "sorption.h"
#include "transport.h"
#include "utils.h"

namespace
{
template <typename ColumnType>
void refreshColumnState(ColumnType& column)
{
  auto updateDerivedState = [&]
  {
    if constexpr (std::is_same_v<ColumnType, Column>)
    {
      computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                             column.maxChemisorptionSites, column.geometry, column.particleDensity,
                             column.concentration, column.physisorptionDot, column.chemisorptionDot,
                             column.surfaceConcentration, column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                             column.reactionChemisorptionSource);
      updateVelocityAndPressure(
          column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
          column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength, column.geometry,
          column.columnEntranceVelocity, column.dynamicViscosity, column.resolution, column.interstitialGasVelocity,
          column.gasDensity, column.totalConcentration, column.totalPressure, column.concentration,
          column.partialPressure, column.moleFraction, column.bulkSpeciesSink, column.gasTemperature,
          column.fluidPhase, column.liquidDensity, column.pHMode, column.pHValue, column.pKw,
          column.pHComponent, column.pH);
      computePhysisorptionEquilibriumLoadings(
          column.physisorptionMixture, column.numberOfGridPoints, column.numberOfComponents, column.maxIsothermTerms,
          column.iastPerformance, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
          column.totalPressure, column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential,
          column.moleFraction, column.gasTemperature,
          column.fluidPhase == Column::FluidPhase::Gas
              ? MixturePrediction::DrivingForceInput::MoleFraction
              : MixturePrediction::DrivingForceInput::Concentration,
          column.concentration, column.pH);
      computeChemisorptionEquilibriumLoadings(
          column.chemisorptionMixture, column.numberOfGridPoints, column.numberOfComponents,
          column.maxChemisorptionSites, column.iastPerformance, column.idealGasMolFractions,
          column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure, column.equilibriumChemisorption,
          column.cachedChemisorptionPressure, column.cachedChemisorptionGrandPotential, column.moleFraction,
          column.gasTemperature,
          column.fluidPhase == Column::FluidPhase::Gas
              ? MixturePrediction::DrivingForceInput::MoleFraction
              : MixturePrediction::DrivingForceInput::Concentration,
          column.concentration, column.pH);
    }
    else
    {
      computeBulkSpeciesSink(column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents,
                             column.numberOfAdsorbents, column.maxChemisorptionSites, column.geometries,
                             column.adsorbentVoidFractions,
                             column.particleDensities, column.particleDiameters, column.fractionOfAdsorbent,
                             column.totalVoidFraction, column.concentration, column.physisorptionDot,
                             column.chemisorptionDot, column.surfaceConcentration, column.bulkSpeciesSink,
                             column.reactionPhysisorptionSource, column.reactionChemisorptionSource);
      updateVelocityAndPressure(
          column.components, column.boundaryCondition, column.numberOfGridPoints, column.numberOfComponents,
          column.inletPressure, column.outletPressure, column.pressureGradient, column.columnLength,
          column.numberOfAdsorbents, column.columnEntranceVelocity, column.dynamicViscosity, column.columnDistances,
          column.fractionOfAdsorbent, column.geometries, column.adsorbentScaledVoidFraction,
          column.totalVoidFraction,
          column.interstitialGasVelocity,
          column.gasDensity, column.totalConcentration, column.totalPressure, column.concentration,
          column.partialPressure, column.moleFraction, column.bulkSpeciesSink, column.gasTemperature,
          column.fluidPhase, column.liquidDensity, column.pHMode, column.pHValue, column.pKw,
          column.pHComponent, column.pH);
      computePhysisorptionEquilibriumLoadings(
          column.physisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
          column.fractionOfAdsorbent, column.hasAdsorbentOfType, column.maxIsothermTerms, column.iastPerformance,
          column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
          column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential, column.moleFraction,
          column.gasTemperature,
          column.fluidPhase == MultibedColumn::FluidPhase::Gas
              ? MixturePrediction::DrivingForceInput::MoleFraction
              : MixturePrediction::DrivingForceInput::Concentration,
          column.concentration, column.pH);
      computeChemisorptionEquilibriumLoadings(
          column.chemisorptionMixtures, column.numberOfGridPoints, column.numberOfComponents, column.numberOfAdsorbents,
          column.fractionOfAdsorbent, column.hasAdsorbentOfType, column.maxChemisorptionSites, column.iastPerformance,
          column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure,
          column.equilibriumChemisorption, column.cachedChemisorptionPressure, column.cachedChemisorptionGrandPotential,
          column.moleFraction, column.gasTemperature,
          column.fluidPhase == MultibedColumn::FluidPhase::Gas
              ? MixturePrediction::DrivingForceInput::MoleFraction
              : MixturePrediction::DrivingForceInput::Concentration,
          column.concentration, column.pH);
    }
  };

  updateDerivedState();
  computeDerivatives(column);
  updateDerivedState();
}
}  // namespace

template <typename ColumnType>
SwingAdsorption<ColumnType>::SwingAdsorption(const InputReader& inputReader) : breakthrough(inputReader)
{
  if (inputReader.swingAdsorptionPhases.empty())
  {
    throw std::invalid_argument("SwingAdsorption: at least one swing adsorption phase must be provided.");
  }

  subStages.reserve(inputReader.swingAdsorptionPhases.size());
  for (const InputReader::SwingAdsorptionPhase& phase : inputReader.swingAdsorptionPhases)
  {
    const double temperature = phase.temperature.value_or(inputReader.temperature);
    const double pressure = phase.inletPressure.value_or(inputReader.inletPressure);
    if (phase.numberOfSteps == 0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has zero time steps.");
    }
    if (temperature < 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid temperature.");
    }
    if (pressure <= 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid inlet pressure.");
    }
    subStages.emplace_back(SwingAdsorptionSubStage{phase.name, temperature, pressure, phase.numberOfSteps});
  }
}

template <typename ColumnType>
void SwingAdsorption<ColumnType>::run()
{
  ColumnType& column = breakthrough.column;

  std::vector<std::ofstream> streams;
  for (size_t comp = 0; comp < breakthrough.numberOfComponents; ++comp)
  {
    const std::string fileName = "component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data";
    streams.emplace_back(std::ofstream{fileName, std::ios_base::app});
  }
  std::ofstream columnStream("column.data", std::ios_base::app);

  size_t step = 0;
  double realTime = 0.0;
  constexpr double initialTimeStepFraction = 1.0e-4;

  const size_t totalNumberOfSteps =
      std::accumulate(subStages.begin(), subStages.end(), size_t{0},
                      [](size_t total, const SwingAdsorptionSubStage& stage) { return total + stage.numberOfSteps; });

  breakthrough.rk3.autoNumberOfSteps = false;
  breakthrough.sirk3.autoNumberOfSteps = false;
  breakthrough.cvode.autoNumberOfSteps = false;
  breakthrough.rk3.numberOfSteps = totalNumberOfSteps;
  breakthrough.sirk3.numberOfSteps = totalNumberOfSteps;
  breakthrough.cvode.numberOfSteps = totalNumberOfSteps;

  if (breakthrough.numberOfInitTimeSteps > 0)
  {
    breakthrough.rk3.timeStep = breakthrough.timeStep * initialTimeStepFraction;
    breakthrough.sirk3.timeStep = breakthrough.timeStep * initialTimeStepFraction;
    breakthrough.cvode.timeStep = breakthrough.timeStep * initialTimeStepFraction;
  }

  auto applyStageConditions = [&](const SwingAdsorptionSubStage& stage)
  {
    column.setTemperature(stage.temperature);
    column.inletPressure = stage.pressure;
    if (!column.gasTemperature.empty())
    {
      column.gasTemperature[0] = stage.temperature;
    }
    if (column.fluidPhase == decltype(column.fluidPhase)::Liquid)
    {
      for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
      {
        column.concentration[comp] = column.components[comp].inletLiquidConcentration;
      }
    }
    else if (column.inletPressure > 0.0)
    {
      const double inletTemperature =
          std::max(1e-10, column.gasTemperature.empty() ? stage.temperature : column.gasTemperature[0]);
      const double inletConcentration = column.inletPressure / (R * inletTemperature);
      for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
      {
        column.concentration[comp] = column.components[comp].initialGasMoleFraction * inletConcentration;
      }
    }
    refreshColumnState(column);
  };

  try
  {
    for (size_t stageIndex = 0; stageIndex < subStages.size(); ++stageIndex)
    {
      const SwingAdsorptionSubStage& stage = subStages[stageIndex];
      applyStageConditions(stage);
      if (breakthrough.integrationScheme == BreakthroughIntegrationScheme::CVODE)
      {
        breakthrough.cvode.reinitialize();
      }
      std::print("\nSwing phase {}: {} (T={} K, inlet pressure={} Pa, steps={})\n", stageIndex + 1, stage.name,
                 stage.temperature, stage.pressure, stage.numberOfSteps);

      for (size_t subStep = 0; subStep < stage.numberOfSteps; ++subStep)
      {
        const double stepTime = realTime;

        if (step % breakthrough.writeEvery == 0)
        {
          column.writeJSON("column.json");
        }

        switch (breakthrough.integrationScheme)
        {
          case BreakthroughIntegrationScheme::SSP_RK:
          {
            static_cast<void>(breakthrough.rk3.propagate(column, step, timings));
            realTime += breakthrough.rk3.timeStep;
            break;
          }
          case BreakthroughIntegrationScheme::CVODE:
          {
            static_cast<void>(breakthrough.cvode.propagate(column, step, timings));
            realTime += breakthrough.cvode.timeStep;
            break;
          }
          case BreakthroughIntegrationScheme::SIRK3:
          {
            if constexpr (std::is_same_v<ColumnType, Column>)
            {
              static_cast<void>(breakthrough.sirk3.propagate(column, step, timings));
              realTime += breakthrough.sirk3.timeStep;
            }
            else
            {
              throw std::runtime_error("Error: SIRK3 is not supported by MultibedColumn");
            }
            break;
          }
        }

        if (step < breakthrough.numberOfInitTimeSteps)
        {
          const double progress = static_cast<double>(step) / static_cast<double>(breakthrough.numberOfInitTimeSteps);
          const double nextTime =
              breakthrough.timeStep *
              (initialTimeStepFraction +
               (1.0 - initialTimeStepFraction) * (3.0 * progress * progress - 2.0 * progress * progress * progress));
          breakthrough.rk3.timeStep = nextTime;
          breakthrough.sirk3.timeStep = nextTime;
          breakthrough.cvode.timeStep = nextTime;
        }
        else if (step == breakthrough.numberOfInitTimeSteps)
        {
          breakthrough.rk3.timeStep = breakthrough.timeStep;
          breakthrough.sirk3.timeStep = breakthrough.timeStep;
          breakthrough.cvode.timeStep = breakthrough.timeStep;
        }

        if (step % breakthrough.writeEvery == 0)
        {
          column.writeOutput(streams, columnStream, stepTime);
        }
        if (step % breakthrough.printEvery == 0)
        {
          const double averageIastSteps = column.iastPerformance.second == 0
                                              ? 0.0
                                              : static_cast<double>(column.iastPerformance.first) /
                                                    static_cast<double>(column.iastPerformance.second);
          std::print("Timestep {}, time: {:6.5f} [s]\n", step, stepTime);
          std::print("    Average number of mixture-prediction steps: {:6.5f}\n", averageIastSteps);
          std::fflush(stdout);
        }
        ++step;
      }
    }
  }
  catch (const std::runtime_error&)
  {
    column.writeJSON("failed_state.json");
    throw;
  }

  std::print("Final timestep {}, time: {:6.5f} [s]\n\n", step, realTime);
  timings.print();
}

template <typename ColumnType>
void SwingAdsorption<ColumnType>::print() const
{
  std::print("{}", repr());
}

template <typename ColumnType>
std::string SwingAdsorption<ColumnType>::repr() const
{
  std::string result = breakthrough.repr();
  result += std::format(
      "Swing adsorption phases\n"
      "=======================================================\n");
  for (size_t i = 0; i < subStages.size(); ++i)
  {
    const SwingAdsorptionSubStage& stage = subStages[i];
    result += std::format("{}: {}  T={} [K], inlet pressure={} [Pa], steps={}\n", i + 1, stage.name, stage.temperature,
                          stage.pressure, stage.numberOfSteps);
  }
  result += "\n\n";
  return result;
}

template struct SwingAdsorption<Column>;
template struct SwingAdsorption<MultibedColumn>;

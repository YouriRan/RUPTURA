#pragma once

#include <array>
#include <cstddef>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "chemisorption.h"
#include "component.h"
#include "inputreader.h"
#include "isotherm.h"
#include "mixture_prediction.h"
#include "multi_site_chemisorption.h"
#include "utils.h"

namespace multibed_benchmark
{
constexpr std::size_t simulationSteps = 100;
constexpr std::size_t gridPoints = 50;
constexpr double absoluteTolerance = 1.0e-6;
constexpr double relativeTolerance = 1.0e-6;

enum class Workload
{
  Physisorption,
  ChemisorptionDirect,
  ChemisorptionTransport,
};

struct Method
{
  MixturePrediction::PredictionMethod value;
  std::string_view name;
};

struct Setting
{
  std::size_t componentCount;
  bool energyBalance;
  Workload workload;
  Method method;
};

constexpr std::array componentCounts{std::size_t{2}, std::size_t{3}, std::size_t{4}, std::size_t{5}, std::size_t{6}};
constexpr std::array energyBalances{false, true};
constexpr std::array workloads{Workload::Physisorption, Workload::ChemisorptionDirect,
                               Workload::ChemisorptionTransport};
constexpr std::array methods{
    Method{MixturePrediction::PredictionMethod::IAST, "IAST"},
    Method{MixturePrediction::PredictionMethod::SIAST, "SIAST"},
    Method{MixturePrediction::PredictionMethod::SEI, "SEI"},
    Method{MixturePrediction::PredictionMethod::SPI, "SPI"},
};

constexpr std::size_t settingCount = componentCounts.size() * energyBalances.size() * workloads.size() * methods.size();
static_assert(settingCount == 120);

constexpr std::array<Setting, settingCount> makeSettings()
{
  std::array<Setting, settingCount> result{};
  std::size_t index = 0;
  for (const std::size_t componentCount : componentCounts)
  {
    for (const bool energyBalance : energyBalances)
    {
      for (const Workload workload : workloads)
      {
        for (const Method method : methods)
        {
          result[index++] = Setting{componentCount, energyBalance, workload, method};
        }
      }
    }
  }
  return result;
}

constexpr auto settings = makeSettings();

constexpr std::string_view workloadName(Workload workload)
{
  switch (workload)
  {
    case Workload::Physisorption:
      return "physisorption";
    case Workload::ChemisorptionDirect:
      return "chemisorption-direct";
    case Workload::ChemisorptionTransport:
      return "chemisorption-transport";
  }
  return "unknown";
}

constexpr bool hasChemisorption(Workload workload) { return workload != Workload::Physisorption; }

constexpr bool hasSurfacePoreTransport(Workload workload) { return workload == Workload::ChemisorptionTransport; }

inline std::string settingName(const Setting& setting)
{
  return "components=" + std::to_string(setting.componentCount) + "/energy=" + (setting.energyBalance ? "on" : "off") +
         "/workload=" + std::string{workloadName(setting.workload)} + "/method=" + std::string{setting.method.name};
}

inline std::string settingDirectoryName(const Setting& setting)
{
  return "components-" + std::to_string(setting.componentCount) + "_energy-" + (setting.energyBalance ? "on" : "off") +
         "_workload-" + std::string{workloadName(setting.workload)} + "_method-" + std::string{setting.method.name};
}

inline Chemisorption makeChemisorptionSite(std::size_t component, bool surfacePoreTransport)
{
  const double componentScale = 1.0 + 0.1 * static_cast<double>(component);
  Chemisorption site;
  site.type = Chemisorption::Type::General;
  site.maximumLoading = 0.4 * componentScale;
  site.heatOfChemisorption = 40000.0 + 1000.0 * static_cast<double>(component);
  site.adsorptionRateCoefficient = 2.5e-5 * componentScale;
  site.desorptionRateCoefficient = 1.0e-6;
  site.filmMassTransferCoefficient = 1.5e-3;
  site.poreDiffusivity = 8.0e-11;
  site.usePoreSurfaceTransport = surfacePoreTransport;
  site.isotherm = Isotherm(Isotherm::Type::Langmuir, {site.maximumLoading, 1.0e-5 * componentScale}, true);
  return site;
}

inline std::vector<Component> makeComponents(const InputReader& base, const Setting& setting)
{
  std::vector<Component> components;
  components.reserve(setting.componentCount);

  Component carrier = base.components.at(0);
  carrier.id = 0;
  carrier.name = "Carrier";
  carrier.initialGasMoleFraction = 0.7;
  carrier.isCarrierGas = true;
  carrier.chemisorption = MultiSiteChemisorption{};
  components.push_back(std::move(carrier));

  const double adsorbateMoleFraction = 0.3 / static_cast<double>(setting.componentCount - 1);
  for (std::size_t component = 1; component < setting.componentCount; ++component)
  {
    Component adsorbate = base.components.at(1 + (component - 1) % (base.components.size() - 1));
    adsorbate.id = component;
    adsorbate.name = "Adsorbate" + std::to_string(component);
    adsorbate.initialGasMoleFraction = adsorbateMoleFraction;
    adsorbate.isCarrierGas = false;
    adsorbate.nonIsothermal = true;
    adsorbate.referenceTemperature = base.temperature;
    adsorbate.heatOfAdsorption = 18000.0 + 2000.0 * static_cast<double>(component);
    adsorbate.chemisorption = MultiSiteChemisorption{};

    const double componentScale = 1.0 + 0.05 * static_cast<double>(component);
    for (Isotherm& site : adsorbate.isotherm.sites)
    {
      site.parameters[0] *= componentScale;
      site.parameters[1] /= componentScale;
      site.nonIsothermal = true;
    }

    if (hasChemisorption(setting.workload))
    {
      adsorbate.chemisorption.add(makeChemisorptionSite(component, hasSurfacePoreTransport(setting.workload)));
    }
    components.push_back(std::move(adsorbate));
  }
  return components;
}

inline InputReader makeInput(const Setting& setting)
{
  const auto inputPath =
      std::filesystem::path{RUPTURA_SOURCE_DIR} / "examples/MOR-CO2-C3H8/breakthrough/simulation.json";
  InputReader input(inputPath.string());
  input.breakthroughIntegrator = 0;
  input.autoNumberOfTimeSteps = false;
  input.numberOfTimeSteps = simulationSteps;
  input.energyBalance = setting.energyBalance;
  input.numberOfGridPoints = gridPoints;
  input.columnDistances = makeUniformColumnDistances(input.numberOfGridPoints, input.columnLength);
  input.mixturePredictionMethod = static_cast<std::size_t>(setting.method.value);
  input.components = makeComponents(input, setting);
  input.adsorbentComponents = {input.components};
  input.numberOfCarrierGases = 1;
  input.carrierGasComponent = 0;
  input.maxIsothermTerms = 2;

  input.internalDiameter = 1.0e-2;
  input.outerDiameter = 1.2e-2;
  input.wallDensity = 7800.0;
  input.gasThermalConductivity = 0.025;
  input.wallThermalConductivity = 15.0;
  input.heatTransferGasSolid = 20.0;
  input.heatTransferGasWall = 10.0;
  input.heatTransferWallExternal = 5.0;
  input.heatCapacityGas = 1000.0;
  input.heatCapacitySolid = 900.0;
  input.heatCapacityWall = 500.0;
  input.geometry = makeGeometry(PackedBedTubeSpec{.voidFraction = input.columnVoidFraction,
                                                  .particleDiameter = input.particleDiameter,
                                                  .internalDiameter = input.internalDiameter,
                                                  .outerDiameter = input.outerDiameter});
  input.adsorbentGeometries = {input.geometry};
  return input;
}
}  // namespace multibed_benchmark

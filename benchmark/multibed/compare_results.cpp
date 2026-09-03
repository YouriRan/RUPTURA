#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "benchmark_cases.h"
#include "column.h"
#include "column_multibed.h"
#include "inputreader.h"
#include "integrators/rk3.h"
#include "json.h"
#include "mixture_prediction.h"
#include "timing.h"
#include "utils.h"

namespace
{
using multibed_benchmark::Setting;
using multibed_benchmark::settings;
using multibed_benchmark::simulationSteps;

struct Errors
{
  double maximumAbsolute{0.0};
  double maximumRelative{0.0};
  bool equal{true};
};

void writeSimulationInput(const Setting& setting, const InputReader& configuredInput,
                          const std::filesystem::path& outputDirectory, bool multibed)
{
  const auto source = std::filesystem::path{RUPTURA_SOURCE_DIR} / "examples/MOR-CO2-C3H8/breakthrough/simulation.json";
  std::ifstream inputStream(source);
  nlohmann::json input;
  inputStream >> input;
  input["energyBalance"] = setting.energyBalance;
  input["NumberOfGridPoints"] = configuredInput.numberOfGridPoints;
  input["MixturePredictionMethod"] = setting.method.name;
  input["BreakthroughIntegrator"] = "RungeKutta3";
  input["NumberOfTimeSteps"] = simulationSteps;
  input["PrintEvery"] = simulationSteps;
  input["WriteEvery"] = 1;
  input["wallDensity"] = configuredInput.wallDensity;
  input["gasThermalConductivity"] = configuredInput.gasThermalConductivity;
  input["wallThermalConductivity"] = configuredInput.wallThermalConductivity;
  input["heatTransferGasSolid"] = configuredInput.heatTransferGasSolid;
  input["heatTransferGasWall"] = configuredInput.heatTransferGasWall;
  input["heatTransferWallExternal"] = configuredInput.heatTransferWallExternal;
  input["heatCapacityGas"] = configuredInput.heatCapacityGas;
  input["heatCapacitySolid"] = configuredInput.heatCapacitySolid;
  input["heatCapacityWall"] = configuredInput.heatCapacityWall;
  input["Geometry"]["InternalDiameter"] = configuredInput.internalDiameter;
  input["Geometry"]["OuterDiameter"] = configuredInput.outerDiameter;

  input["Components"] = nlohmann::json::array();
  for (const Component& component : configuredInput.components)
  {
    nlohmann::json componentJson{{"Name", component.name},
                                 {"GasPhaseMolFraction", component.initialGasMoleFraction},
                                 {"MolecularWeight", component.molecularWeight}};
    if (component.isCarrierGas)
    {
      componentJson["CarrierGas"] = true;
    }
    else
    {
      componentJson["MassTransferCoefficient"] = component.massTransferCoefficient;
      componentJson["AxialDispersionCoefficient"] = component.axialDispersionCoefficient;
      componentJson["HeatOfAdsorption"] = component.heatOfAdsorption;
      componentJson["referenceTemperature"] = component.referenceTemperature.value_or(configuredInput.temperature);
      componentJson["PhysisorptionSites"] = nlohmann::json::array();
      for (const Isotherm& site : component.isotherm.sites)
      {
        componentJson["PhysisorptionSites"].push_back({{"Type", "Langmuir"}, {"Parameters", site.parameters}});
      }

      if (component.chemisorption.enabled())
      {
        componentJson["ChemisorptionSites"] = nlohmann::json::array();
        for (const Chemisorption& site : component.chemisorption.sites)
        {
          componentJson["ChemisorptionSites"].push_back(
              {{"Type", "General"},
               {"Parameters",
                {{"maximumLoading", site.maximumLoading},
                 {"heatOfChemisorption", site.heatOfChemisorption},
                 {"adsorptionRateCoefficient", site.adsorptionRateCoefficient},
                 {"adsorptionActivationEnergy", site.adsorptionActivationEnergy},
                 {"desorptionRateCoefficient", site.desorptionRateCoefficient},
                 {"desorptionActivationEnergy", site.desorptionActivationEnergy},
                 {"poreConcentrationOrder", site.poreConcentrationOrder},
                 {"capacityOrder", site.capacityOrder},
                 {"desorptionOrder", site.desorptionOrder},
                 {"filmMassTransferCoefficient", site.filmMassTransferCoefficient},
                 {"poreDiffusivity", site.poreDiffusivity},
                 {"usePoreSurfaceTransport", site.usePoreSurfaceTransport},
                 {"Isotherm", {{"Type", "Langmuir"}, {"Parameters", site.isotherm.value().parameters}}}}}});
        }
      }
    }
    input["Components"].push_back(std::move(componentJson));
  }
  if (multibed)
  {
    input["Adsorbents"] = nlohmann::json::array({{{"Name", "MOR"},
                                                  {"AdsorbentLength", input["ColumnLength"]},
                                                  {"ColumnVoidFraction", input["ColumnVoidFraction"]},
                                                  {"ParticleDensity", input["ParticleDensity"]}}});
  }
  else
  {
    input.erase("Adsorbents");
    input.erase("ColumnSections");
  }

  std::ofstream output(outputDirectory / "simulation.json", std::ios::trunc);
  output << std::setw(2) << input << '\n';
}

template <typename ColumnType>
ColumnType run(const InputReader& input, const std::filesystem::path& outputDirectory)
{
  std::filesystem::create_directories(outputDirectory);
  ColumnType column(input);
  column.initialize();
  std::vector<std::ofstream> componentStreams(column.numberOfComponents);
  for (std::size_t component = 0; component < column.numberOfComponents; ++component)
  {
    const std::string fileName =
        "component_" + std::to_string(component) + "_" + column.components[component].name + ".data";
    componentStreams[component].open(outputDirectory / fileName, std::ios::trunc);
  }
  std::ofstream columnStream(outputDirectory / "column.data", std::ios::trunc);
  column.writeOutputHeader(componentStreams, columnStream);
  column.writeOutput(componentStreams, columnStream, 0.0);

  RungeKutta3 integrator(input.timeStep, false, simulationSteps);
  Timing timing;
  for (std::size_t step = 0; step < simulationSteps; ++step)
  {
    static_cast<void>(integrator.propagate(column, step, timing));
  }
  column.writeOutput(componentStreams, columnStream, static_cast<double>(simulationSteps) * input.timeStep);
  column.writeJSON((outputDirectory / "column.json").string());
  return column;
}

Errors writeValues(std::ofstream& output, std::string_view field, std::span<const double> singleBed,
                   std::span<const double> multibed)
{
  Errors errors;
  for (std::size_t index = 0; index < singleBed.size(); ++index)
  {
    const double absolute = std::abs(singleBed[index] - multibed[index]);
    const double scale = std::max({1.0, std::abs(singleBed[index]), std::abs(multibed[index])});
    const double relative = absolute / scale;
    errors.maximumAbsolute = std::max(errors.maximumAbsolute, absolute);
    errors.maximumRelative = std::max(errors.maximumRelative, relative);
    errors.equal = errors.equal && (absolute <= multibed_benchmark::absoluteTolerance ||
                                    relative <= multibed_benchmark::relativeTolerance);
    output << field << ',' << index << ',' << singleBed[index] << ',' << multibed[index] << ',' << absolute << ','
           << relative << '\n';
  }
  return errors;
}
}  // namespace

int main(int argc, char** argv)
{
  const std::filesystem::path outputRoot =
      argc > 1 ? std::filesystem::path{argv[1]} : std::filesystem::path{RUPTURA_SOURCE_DIR} / "multibed";
  std::filesystem::create_directories(outputRoot);

  std::ofstream summary(outputRoot / "summary.csv");
  summary << "setting,component_count,energy_balance,workload,chemisorption,surface_pore_transport,method,"
             "grid_points,max_state_abs,max_state_rel,max_pressure_abs,"
             "max_pressure_rel,max_velocity_abs,max_velocity_rel,results_equal\n";
  summary << std::setprecision(17);

  std::cout << std::left << std::setw(43) << "setting" << std::right << std::setw(17) << "state abs" << std::setw(17)
            << "state rel" << std::setw(17) << "pressure abs" << std::setw(17) << "velocity abs" << '\n';

  std::size_t failures = 0;
  for (const Setting& setting : settings)
  {
    const InputReader input = multibed_benchmark::makeInput(setting);
    const std::string name = multibed_benchmark::settingDirectoryName(setting);
    const std::filesystem::path directory = outputRoot / name;
    std::filesystem::create_directories(directory);
    const Column singleBed = run<Column>(input, directory / "singlebed");
    const MultibedColumn multibed = run<MultibedColumn>(input, directory / "multibed");
    writeSimulationInput(setting, input, directory / "singlebed", false);
    writeSimulationInput(setting, input, directory / "multibed", true);

    std::ofstream values(directory / "comparison.csv");
    values << "field,index,singlebed,multibed,absolute_difference,relative_difference\n";
    values << std::setprecision(17);
    const Errors stateErrors = writeValues(values, "state", singleBed.state, multibed.state);
    const Errors pressureErrors =
        writeValues(values, "total_pressure", singleBed.totalPressure, multibed.totalPressure);
    const Errors velocityErrors = writeValues(values, "interstitial_gas_velocity", singleBed.interstitialGasVelocity,
                                              multibed.interstitialGasVelocity);
    const bool resultsEqual = stateErrors.equal && pressureErrors.equal && velocityErrors.equal;
    failures += resultsEqual ? 0 : 1;

    summary << name << ',' << setting.componentCount << ',' << (setting.energyBalance ? "on" : "off") << ','
            << multibed_benchmark::workloadName(setting.workload) << ','
            << (multibed_benchmark::hasChemisorption(setting.workload) ? "on" : "off") << ','
            << (multibed_benchmark::hasSurfacePoreTransport(setting.workload) ? "on" : "off") << ','
            << setting.method.name << ',' << input.numberOfGridPoints << ',' << stateErrors.maximumAbsolute << ','
            << stateErrors.maximumRelative << ',' << pressureErrors.maximumAbsolute << ','
            << pressureErrors.maximumRelative << ',' << velocityErrors.maximumAbsolute << ','
            << velocityErrors.maximumRelative << ',' << (resultsEqual ? "true" : "false") << '\n';

    std::cout << std::left << std::setw(43) << name << std::right << std::scientific << std::setprecision(6)
              << std::setw(17) << stateErrors.maximumAbsolute << std::setw(17) << stateErrors.maximumRelative
              << std::setw(17) << pressureErrors.maximumAbsolute << std::setw(17) << velocityErrors.maximumAbsolute
              << '\n';
  }

  std::cout << "\nDetailed results written to " << outputRoot << '\n';
  std::cout << (settings.size() - failures) << '/' << settings.size() << " settings agree within tolerance\n";
  return failures == 0 ? 0 : 1;
}

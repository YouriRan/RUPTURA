#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

#include "column.h"
#include "component.h"
#include "fitting.h"
#include "inputreader.h"
#include "integrators/compute.h"
#include "integrators/cvode.h"
#include "integrators/rk3.h"
#include "integrators/rk3_si.h"
#include "integrators/sorption.h"
#include "integrators/transport.h"
#include "isotherm.h"
#include "json.h"
#include "swing_adsorption.h"
#include "utils.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

double langmuir(double saturationLoading, double affinity, double pressure)
{
  const double bp = affinity * pressure;
  return saturationLoading * bp / (1.0 + bp);
}

double langmuirFreundlich(double saturationLoading, double affinity, double exponent, double pressure)
{
  const double bp = affinity * std::pow(pressure, exponent);
  return saturationLoading * bp / (1.0 + bp);
}

std::vector<double> toVector(std::span<const double> values) { return {values.begin(), values.end()}; }

double maxAbsDifference(std::span<const double> left, const std::vector<double>& right)
{
  double maxDifference = 0.0;
  for (size_t i = 0; i < left.size(); ++i)
  {
    maxDifference = std::max(maxDifference, std::abs(left[i] - right[i]));
  }
  return maxDifference;
}

double expectedMassDerivative(const Column& column, size_t grid, size_t comp)
{
  const size_t n = column.numberOfComponents;
  const size_t componentBlockSize = column.concentration.size();
  const double invDz = 1.0 / column.resolution;
  const double prefactor = column.geometry.loadingPrefactor(column.particleDensity);

  if (grid == 0)
  {
    return 0.0;
  }

  const double advectiveFluxDifference =
      (column.interstitialGasVelocity[grid] * column.concentration[grid * n + comp] -
       column.interstitialGasVelocity[grid - 1] * column.concentration[(grid - 1) * n + comp]) *
      invDz;

  const size_t index = grid * n + comp;
  double chemisorptionDot = 0.0;
  for (size_t site = 0; site < column.components[comp].chemisorption.numberOfSites; ++site)
  {
    chemisorptionDot += column.chemisorptionDot[site * componentBlockSize + index];
  }
  return -advectiveFluxDifference - prefactor * (column.physisorptionDot[index] + chemisorptionDot);
}

Column loadBeaBreakthroughColumn()
{
  InputReader reader(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
  Column column(reader);
  column.initialize();
  column.readJSON(examplePath("examples/BEA-alkanes-C7/breakthrough/column.json"));
  return column;
}

void writeGeneralChemisorptionInput(const std::filesystem::path& path, bool usePoreSurfaceTransport,
                                    const std::string& concentrationOrder = "1")
{
  std::ofstream out(path);
  out << R"json(
{
  "SimulationType": "Breakthrough",
  "BreakthroughIntegrator": "RungeKutta3",
  "BoundaryCondition": "InletPressureInletVelocity",
  "Temperature": 300.0,
  "ColumnVoidFraction": 0.4,
  "DynamicViscosity": 1.0e-5,
  "ParticleDiameter": 1.0e-3,
  "Geometry": {
    "Type": "HollowTube",
    "ColumnVoidFraction": 0.4,
    "ParticleDiameter": 1.0e-3
  },
  "ParticleDensity": 1000.0,
  "InletPressure": 100000.0,
  "ColumnEntranceVelocity": 0.1,
  "ColumnLength": 0.1,
  "NumberOfGridPoints": 2,
  "NumberOfTimeSteps": 1,
  "TimeStep": 0.01,
  "Components": [
    {
      "Name": "He",
      "CarrierGas": true,
      "GasPhaseMolFraction": 0.9,
      "MolecularWeight": 0.004
    },
    {
      "Name": "CO2",
      "GasPhaseMolFraction": 0.1,
      "MassTransferCoefficient": 0.0,
      "AxialDispersionCoefficient": 0.0,
      "MolecularWeight": 0.044,
      "ChemisorptionSites": [{
        "Type": "General",
        "Parameters": {
          "maximumLoading": 2.0,
          "heatOfChemisorption": 42000.0,
          "adsorptionRateCoefficient": 1.0e-4,
          "adsorptionActivationEnergy": 0.0,
          "desorptionRateCoefficient": 0.0,
          "desorptionActivationEnergy": 0.0,
          "poreConcentrationOrder": )json"
      << concentrationOrder << R"json(,
          "capacityOrder": 1,
          "desorptionOrder": 1,
          "filmMassTransferCoefficient": 2.0e-3,
          "poreDiffusivity": 1.0e-10,
          "usePoreSurfaceTransport": )json"
      << (usePoreSurfaceTransport ? "true" : "false") << R"json(,
          "Isotherm": {
            "Type": "Langmuir",
            "Parameters": [2.0, 1.0e-5]
          }
        }
      }],
      "PhysisorptionSites": [
        {"Type": "Langmuir", "Parameters": [1.0, 1.0e-5]}
      ]
    }
  ]
}
)json";
}

void writeElovichChemisorptionInput(const std::filesystem::path& path, bool usePoreSurfaceTransport)
{
  writeGeneralChemisorptionInput(path, usePoreSurfaceTransport);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  auto& chemisorption = data["Components"][1]["ChemisorptionSites"][0];
  chemisorption["Type"] = "Elovich";
  chemisorption["Parameters"] = {
      {"maximumLoading", 2.0},
      {"heatOfChemisorption", 42000.0},
      {"alpha", 0.25},
      {"beta", 0.5},
      {"filmMassTransferCoefficient", 2.0e-3},
      {"poreDiffusivity", 1.0e-10},
      {"usePoreSurfaceTransport", usePoreSurfaceTransport},
      {"Isotherm", {{"Type", "Langmuir"}, {"Parameters", {2.0, 1.0e-5}}}},
  };
  input.close();

  std::ofstream output(path);
  output << data;
}

void writeMultiSiteChemisorptionInput(const std::filesystem::path& path, bool secondSiteUsesPoreSurfaceTransport)
{
  writeGeneralChemisorptionInput(path, false);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  auto& component = data["Components"][1];
  nlohmann::json general = component["ChemisorptionSites"][0];
  general["Parameters"]["maximumLoading"] = 1.0;
  general["Parameters"]["Isotherm"]["Parameters"] = {1.0, 1.0e-5};
  nlohmann::json elovich = {
      {"Type", "Elovich"},
      {"Parameters",
       {
           {"maximumLoading", 2.0},
           {"heatOfChemisorption", 42000.0},
           {"alpha", 0.25},
           {"beta", 0.5},
           {"filmMassTransferCoefficient", 2.0e-3},
           {"poreDiffusivity", 1.0e-10},
           {"usePoreSurfaceTransport", secondSiteUsesPoreSurfaceTransport},
           {"Isotherm", {{"Type", "Langmuir"}, {"Parameters", {2.0, 2.0e-5}}}},
       }},
  };
  component["ChemisorptionSites"] = nlohmann::json::array({general, elovich});
  input.close();

  std::ofstream output(path);
  output << data;
}

void writeReactionInput(const std::filesystem::path& path, const std::string& phase = "Physisorbed",
                        const std::string& style = "GeneralPowerLaw", bool energyBalance = true,
                        const std::string& reactants = "[1]", const std::string& products = "[2]")
{
  std::ofstream out(path);
  out << R"json(
{
  "SimulationType": "Breakthrough",
  "BreakthroughIntegrator": "RungeKutta3",
  "BoundaryCondition": "InletPressureInletVelocity",
  "Temperature": 300.0,
  "ColumnVoidFraction": 0.4,
  "DynamicViscosity": 1.0e-5,
  "ParticleDiameter": 1.0e-3,
  "Geometry": {
    "Type": "HollowTube",
    "ColumnVoidFraction": 0.4,
    "ParticleDiameter": 1.0e-3
  },
  "ParticleDensity": 1000.0,
  "InletPressure": 100000.0,
  "ColumnEntranceVelocity": 0.1,
  "ColumnLength": 0.1,
  "NumberOfGridPoints": 2,
  "NumberOfTimeSteps": 1,
  "TimeStep": 0.01,
  "energyBalance": )json"
      << (energyBalance ? "true" : "false") << R"json(,
  "heatCapacityGas": 1.0,
  "heatCapacitySolid": 1000.0,
  "heatCapacityWall": 1.0,
  "wallDensity": 1.0,
  "Components": [
    {
      "Name": "He",
      "CarrierGas": true,
      "GasPhaseMolFraction": 0.8,
      "MolecularWeight": 0.004
    },
    {
      "Name": "A",
      "GasPhaseMolFraction": 0.2,
      "MassTransferCoefficient": 0.0,
      "AxialDispersionCoefficient": 0.0,
      "MolecularWeight": 0.044,
      "HeatOfAdsorption": 0.0,
      "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [1.0, 1.0e-5]}]
    },
    {
      "Name": "B",
      "GasPhaseMolFraction": 0.0,
      "MassTransferCoefficient": 0.0,
      "AxialDispersionCoefficient": 0.0,
      "MolecularWeight": 0.028,
      "HeatOfAdsorption": 0.0,
      "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [1.0, 1.0e-5]}]
    }
  ],
  "Reactions": [
    {
      "Phase": ")json"
      << phase << R"json(",
      "Style": ")json"
      << style << R"json(",
      "Reactants": )json"
      << reactants << R"json(,
      "Products": )json"
      << products << R"json(,
      "Stoichiometry": {
        "Reactants": [1.0],
        "Products": [1.0]
      },
      "Kinetics": {
        "forwardRateCoefficient": 1.0,
        "forwardActivationEnergy": 0.0,
        "equilibriumConstant": 1.0e9,
        "gibbsFreeEnergy": -1000.0
      }
    }
  ]
}
)json";
}
}  // namespace

TEST(ExampleInputReader, LoadsBeaAlkanesMixturePrediction)
{
  InputReader reader(examplePath("examples/BEA-alkanes-C7/iast/simulation.json"));

  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::MixturePrediction);
  EXPECT_EQ(reader.displayName, "BEA");
  EXPECT_DOUBLE_EQ(reader.temperature, 552.0);
  EXPECT_DOUBLE_EQ(reader.pressureStart, 1000.0);
  EXPECT_DOUBLE_EQ(reader.pressureEnd, 10000000.0);
  EXPECT_EQ(reader.numberOfPressurePoints, 100U);
  EXPECT_EQ(reader.pressureScale, 0U);

  ASSERT_EQ(reader.components.size(), 8U);
  EXPECT_EQ(reader.numberOfCarrierGases, 1U);
  EXPECT_EQ(reader.carrierGasComponent, 0U);
  EXPECT_EQ(reader.maxIsothermTerms, 2U);

  const Component& helium = reader.components[0];
  EXPECT_EQ(helium.name, "Helium");
  EXPECT_TRUE(helium.isCarrierGas);
  EXPECT_DOUBLE_EQ(helium.initialGasMoleFraction, 0.93);
  ASSERT_EQ(helium.isotherm.sites.size(), 1U);
  EXPECT_EQ(helium.isotherm.sites[0].type, Isotherm::Type::Langmuir);
  EXPECT_DOUBLE_EQ(helium.isotherm.sites[0].parameters[0], 1.0);
  EXPECT_DOUBLE_EQ(helium.isotherm.sites[0].parameters[1], 0.0);

  const Component& nC7 = reader.components[1];
  EXPECT_EQ(nC7.name, "nC7");
  EXPECT_DOUBLE_EQ(nC7.initialGasMoleFraction, 0.01);
  ASSERT_EQ(nC7.isotherm.sites.size(), 2U);
  EXPECT_EQ(nC7.isotherm.numberOfParameters, 4U);
  EXPECT_EQ(nC7.isotherm.sites[0].type, Isotherm::Type::Langmuir);
  EXPECT_DOUBLE_EQ(nC7.isotherm.sites[0].parameters[0], 1.09984);
  EXPECT_DOUBLE_EQ(nC7.isotherm.sites[0].parameters[1], 6.55857e-05);
  EXPECT_DOUBLE_EQ(nC7.isotherm.sites[1].parameters[0], 0.19466);
  EXPECT_DOUBLE_EQ(nC7.isotherm.sites[1].parameters[1], 8.90731e-07);

  constexpr double pressure = 100000.0;
  const double expectedLoading = langmuir(1.09984, 6.55857e-05, pressure) + langmuir(0.19466, 8.90731e-07, pressure);
  EXPECT_NEAR(nC7.isotherm.value(pressure, 1.0), expectedLoading, 1e-12);

  ASSERT_EQ(reader.adsorbentComponents.size(), 1U);
  EXPECT_EQ(reader.adsorbentComponents[0].size(), reader.components.size());
  EXPECT_EQ(reader.columnDistances.size(), reader.numberOfGridPoints + 1);
  EXPECT_DOUBLE_EQ(reader.columnDistances.front(), 0.0);
  EXPECT_DOUBLE_EQ(reader.columnDistances.back(), reader.columnLength);
}

TEST(ExampleInputReader, LoadsCoBdpAlkanesBreakthrough)
{
  InputReader reader(examplePath("examples/CoBDP-alkanes-C6/breakthrough/simulation.json"));

  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::Breakthrough);
  EXPECT_EQ(reader.displayName, "CoBDP");
  EXPECT_DOUBLE_EQ(reader.temperature, 443.0);
  EXPECT_DOUBLE_EQ(reader.columnVoidFraction, 0.4);
  EXPECT_DOUBLE_EQ(reader.particleDensity, 703.4);
  EXPECT_DOUBLE_EQ(reader.inletPressure, 2000000.0);
  EXPECT_DOUBLE_EQ(reader.pressureGradient, 0.0);
  EXPECT_DOUBLE_EQ(reader.columnEntranceVelocity, 0.1);
  EXPECT_DOUBLE_EQ(reader.columnLength, 0.3);
  EXPECT_EQ(reader.breakthroughIntegrator, 1U);
  EXPECT_TRUE(reader.autoNumberOfTimeSteps);
  EXPECT_DOUBLE_EQ(reader.timeStep, 0.005);
  EXPECT_EQ(reader.numberOfGridPoints, 200U);

  ASSERT_EQ(reader.components.size(), 6U);
  EXPECT_EQ(reader.numberOfCarrierGases, 1U);
  EXPECT_EQ(reader.carrierGasComponent, 0U);
  EXPECT_EQ(reader.maxIsothermTerms, 2U);

  const Component& nC6 = reader.components[1];
  EXPECT_EQ(nC6.name, "nC6");
  EXPECT_DOUBLE_EQ(nC6.initialGasMoleFraction, 0.01);
  EXPECT_DOUBLE_EQ(nC6.massTransferCoefficient, 0.06);
  EXPECT_DOUBLE_EQ(nC6.axialDispersionCoefficient, 0.0);
  ASSERT_EQ(nC6.isotherm.sites.size(), 2U);
  EXPECT_EQ(nC6.isotherm.sites[0].type, Isotherm::Type::Langmuir_Freundlich);
  EXPECT_DOUBLE_EQ(nC6.isotherm.sites[0].parameters[0], 1.1);
  EXPECT_DOUBLE_EQ(nC6.isotherm.sites[0].parameters[1], 0.000214);
  EXPECT_DOUBLE_EQ(nC6.isotherm.sites[0].parameters[2], 0.6);
  EXPECT_DOUBLE_EQ(nC6.isotherm.sites[1].parameters[0], 4.8);
  EXPECT_DOUBLE_EQ(nC6.isotherm.sites[1].parameters[1], 2.39e-05);
  EXPECT_DOUBLE_EQ(nC6.isotherm.sites[1].parameters[2], 1.35);

  constexpr double pressure = 250000.0;
  const double expectedLoading =
      langmuirFreundlich(1.1, 0.000214, 0.6, pressure) + langmuirFreundlich(4.8, 2.39e-05, 1.35, pressure);
  EXPECT_NEAR(nC6.isotherm.value(pressure, 1.0), expectedLoading, 1e-12);

  ASSERT_EQ(reader.adsorbentComponents.size(), 1U);
  EXPECT_EQ(reader.adsorbentComponents[0].size(), reader.components.size());
  EXPECT_EQ(reader.columnDistances.size(), 201U);
  EXPECT_DOUBLE_EQ(reader.columnDistances.front(), 0.0);
  EXPECT_DOUBLE_EQ(reader.columnDistances.back(), 0.3);
}

TEST(ExampleInputReader, ParsesSemiImplicitRungeKuttaIntegrator)
{
  std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();
  data["BreakthroughIntegrator"] = "SIRK3";
  data["NumberOfTimeSteps"] = 1;

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_sirk3_integrator.json";
  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  EXPECT_EQ(reader.breakthroughIntegrator, 2U);
  EXPECT_EQ(static_cast<BreakthroughIntegrationScheme>(reader.breakthroughIntegrator),
            BreakthroughIntegrationScheme::SIRK3);

  SemiImplicitRungeKutta3 integrator(reader);
  EXPECT_DOUBLE_EQ(integrator.timeStep, reader.timeStep);
  EXPECT_EQ(integrator.numberOfSteps, reader.numberOfTimeSteps);

  Column column(reader);
  column.initialize();
  Timing timings;
  EXPECT_TRUE(integrator.propagate(column, 0, timings));
}

TEST(SwingAdsorptionInput, ParsesStructuredPhasesAndBuildsSubStages)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_swing_structured_test.json";
  writeGeneralChemisorptionInput(path, false);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();

  data["SimulationType"] = "SwingAdsorption";
  data.erase("InletPressure");
  data["NumberOfTimeSteps"] = "auto";
  data["SwingAdsorptionPhases"] = nlohmann::json::array({
      {{"Name", "Adsorption"}, {"Temperature", 310.0}, {"InletPressure", 200000.0}, {"NumberOfTimeSteps", 3}},
      {{"Name", "Depressurization"}, {"InletPressure", 50000.0}, {"NumberOfTimeSteps", 2}},
  });

  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::SwingAdsorption);
  EXPECT_FALSE(reader.autoNumberOfTimeSteps);
  EXPECT_EQ(reader.numberOfTimeSteps, 5U);
  EXPECT_DOUBLE_EQ(reader.inletPressure, 200000.0);

  ASSERT_EQ(reader.swingAdsorptionPhases.size(), 2U);
  EXPECT_EQ(reader.swingAdsorptionPhases[0].name, "Adsorption");
  ASSERT_TRUE(reader.swingAdsorptionPhases[0].temperature.has_value());
  EXPECT_DOUBLE_EQ(*reader.swingAdsorptionPhases[0].temperature, 310.0);
  ASSERT_TRUE(reader.swingAdsorptionPhases[0].inletPressure.has_value());
  EXPECT_DOUBLE_EQ(*reader.swingAdsorptionPhases[0].inletPressure, 200000.0);
  EXPECT_EQ(reader.swingAdsorptionPhases[0].numberOfSteps, 3U);
  EXPECT_FALSE(reader.swingAdsorptionPhases[1].temperature.has_value());
  ASSERT_TRUE(reader.swingAdsorptionPhases[1].inletPressure.has_value());
  EXPECT_DOUBLE_EQ(*reader.swingAdsorptionPhases[1].inletPressure, 50000.0);

  SwingAdsorption<Column> swing(reader);
  ASSERT_EQ(swing.subStages.size(), 2U);
  EXPECT_EQ(swing.subStages[0].name, "Adsorption");
  EXPECT_DOUBLE_EQ(swing.subStages[0].temperature, 310.0);
  EXPECT_DOUBLE_EQ(swing.subStages[0].pressure, 200000.0);
  EXPECT_EQ(swing.subStages[0].numberOfSteps, 3U);
  EXPECT_EQ(swing.subStages[1].name, "Depressurization");
  EXPECT_DOUBLE_EQ(swing.subStages[1].temperature, reader.temperature);
  EXPECT_DOUBLE_EQ(swing.subStages[1].pressure, 50000.0);
  EXPECT_EQ(swing.subStages[1].numberOfSteps, 2U);
}

TEST(SwingAdsorptionInput, LoadsPsaAndTsaExamples)
{
  InputReader psaReader(examplePath("examples/BEA-alkanes-C7/psa/simulation.json"));
  EXPECT_EQ(psaReader.simulationType, InputReader::SimulationType::SwingAdsorption);
  ASSERT_EQ(psaReader.swingAdsorptionPhases.size(), 2U);
  EXPECT_EQ(psaReader.swingAdsorptionPhases[0].name, "Adsorption");
  EXPECT_EQ(psaReader.swingAdsorptionPhases[1].name, "Blowdown");
  ASSERT_TRUE(psaReader.swingAdsorptionPhases[1].inletPressure.has_value());
  EXPECT_DOUBLE_EQ(*psaReader.swingAdsorptionPhases[1].inletPressure, 1000.0);

  InputReader tsaReader(examplePath("examples/BEA-alkanes-C7/tsa/simulation.json"));
  EXPECT_EQ(tsaReader.simulationType, InputReader::SimulationType::SwingAdsorption);
  ASSERT_EQ(tsaReader.swingAdsorptionPhases.size(), 2U);
  EXPECT_EQ(tsaReader.swingAdsorptionPhases[1].name, "Thermal regeneration");
  ASSERT_TRUE(tsaReader.swingAdsorptionPhases[1].temperature.has_value());
  EXPECT_DOUBLE_EQ(*tsaReader.swingAdsorptionPhases[1].temperature, 673.0);
}

TEST(GeometryInput, LoadsExplicitHollowTubeGeometry)
{
  InputReader reader(examplePath("examples/CoBDP-alkanes-C6/breakthrough/simulation.json"));
  const Geometry& shape = reader.geometry;

  EXPECT_EQ(shape.kind, GeometryKind::HollowTube);
  EXPECT_DOUBLE_EQ(shape.voidFraction, reader.columnVoidFraction);
  EXPECT_DOUBLE_EQ(shape.solidToFluidVolumeRatio, 1.5);
  EXPECT_DOUBLE_EQ(shape.contactAreas.solidFluidPerSolidVolume, 6000.0);
  EXPECT_DOUBLE_EQ(shape.contactAreas.fluidSolidPerFluidVolume, 9000.0);

  Column column(reader);
  EXPECT_EQ(column.geometry.kind, GeometryKind::HollowTube);
}

TEST(GeometryInput, RequiresGeometryForBreakthrough)
{
  std::ifstream input(examplePath("examples/CoBDP-alkanes-C6/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  data.erase("Geometry");
  input.close();

  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_missing_geometry_test.json";
  std::ofstream output(path);
  output << data;
  output.close();

  EXPECT_THROW(
      {
        try
        {
          static_cast<void>(InputReader(path.string()));
        }
        catch (const std::runtime_error& error)
        {
          EXPECT_NE(std::string(error.what()).find("Geometry is required"), std::string::npos);
          throw;
        }
      },
      std::runtime_error);
}

TEST(GeometryInput, ParsesMonolithSubJsonAndPrecomputesChannelTerms)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_monolith_geometry_test.json";
  writeGeneralChemisorptionInput(path, false);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  data["Geometry"] = {
      {"Type", "Monolith"},      {"ChannelShape", "square"}, {"InternalChannelDimension", 1.0e-3},
      {"OuterDiameter", 1.0e-2}, {"NumberOfChannels", 10},   {"WashcoatThickness", 1.0e-4},
  };
  input.close();

  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  const Geometry& shape = reader.geometry;
  const double outerArea = std::acos(-1.0) * 1.0e-2 * 1.0e-2 / 4.0;
  const double openArea = 10.0e-6;

  EXPECT_EQ(shape.kind, GeometryKind::Monolith);
  EXPECT_EQ(channelShapeName(shape.dimensions.channelShape), "square");
  EXPECT_EQ(shape.dimensions.numberOfChannels, 10U);
  EXPECT_DOUBLE_EQ(shape.dimensions.channelArea, 1.0e-6);
  EXPECT_DOUBLE_EQ(shape.dimensions.channelPerimeter, 4.0e-3);
  EXPECT_DOUBLE_EQ(shape.dimensions.hydraulicDiameter, 1.0e-3);
  EXPECT_NEAR(shape.voidFraction, openArea / outerArea, 1.0e-14);
  EXPECT_DOUBLE_EQ(shape.contactAreas.fluidSolidPerFluidVolume, 4000.0);
  EXPECT_DOUBLE_EQ(shape.solidToFluidVolumeRatio, 0.4);
  EXPECT_DOUBLE_EQ(shape.contactAreas.solidFluidPerSolidVolume, 10000.0);
  EXPECT_DOUBLE_EQ(shape.dimensions.poreDiffusionLength, 1.0e-4);
  EXPECT_DOUBLE_EQ(reader.columnVoidFraction, shape.voidFraction);

  Column column(reader);
  column.initialize();
  EXPECT_EQ(column.geometry.kind, GeometryKind::Monolith);
}

TEST(GeometryInput, LoadsChemisorptionCo2N2MonolithExample)
{
  InputReader reader(examplePath("examples/Chemisorption-CO2-N2/monolith/simulation.json"));
  const Geometry& shape = reader.geometry;

  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::Breakthrough);
  EXPECT_EQ(reader.displayName, "Chemisorption CO2/N2 monolith");
  EXPECT_EQ(shape.kind, GeometryKind::Monolith);
  EXPECT_EQ(channelShapeName(shape.dimensions.channelShape), "square");
  EXPECT_EQ(shape.dimensions.numberOfChannels, 10U);
  EXPECT_DOUBLE_EQ(shape.dimensions.channelArea, 1.0e-6);
  EXPECT_DOUBLE_EQ(shape.dimensions.channelPerimeter, 4.0e-3);
  EXPECT_DOUBLE_EQ(shape.dimensions.hydraulicDiameter, 1.0e-3);
  EXPECT_DOUBLE_EQ(shape.contactAreas.fluidSolidPerFluidVolume, 4000.0);
  EXPECT_DOUBLE_EQ(shape.solidToFluidVolumeRatio, 0.4);
  EXPECT_DOUBLE_EQ(shape.dimensions.poreDiffusionLength, 1.0e-4);
  EXPECT_DOUBLE_EQ(reader.columnVoidFraction, shape.voidFraction);

  ASSERT_EQ(reader.components.size(), 3U);
  EXPECT_EQ(reader.numberOfCarrierGases, 1U);
  EXPECT_EQ(reader.carrierGasComponent, 0U);
  EXPECT_EQ(reader.components[1].name, "CO2");
  EXPECT_EQ(reader.components[2].name, "N2");

  Column column(reader);
  column.initialize();
  EXPECT_EQ(column.geometry.kind, GeometryKind::Monolith);
}

TEST(ExampleMassBalance, RecomputesBeaColumnStateDerivatives)
{
  Column column = loadBeaBreakthroughColumn();
  const std::vector<double> savedConcentrationDot = toVector(column.concentrationDot);

  RK3Helpers::computeMassDerivatives(column);

  EXPECT_LE(maxAbsDifference(column.concentrationDot, savedConcentrationDot), 1e-10);

  for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
  {
    EXPECT_DOUBLE_EQ(column.concentrationDot[comp], 0.0);
  }

  const size_t n = column.numberOfComponents;
  EXPECT_NEAR(column.concentrationDot[1 * n + 1], expectedMassDerivative(column, 1, 1), 1e-14);
  EXPECT_NEAR(column.concentrationDot[10 * n + 1], expectedMassDerivative(column, 10, 1), 1e-12);
  EXPECT_NEAR(column.concentrationDot[25 * n + 0], expectedMassDerivative(column, 25, 0), 1e-12);
  EXPECT_NEAR(column.concentrationDot[25 * n + 2], expectedMassDerivative(column, 25, 2), 1e-12);
  EXPECT_NEAR(column.concentrationDot[50 * n + 5], expectedMassDerivative(column, 50, 5), 1e-12);

  EXPECT_NEAR(column.concentrationDot[25 * n + 0], -2.4276051365836119e-06, 1e-14);
  EXPECT_NEAR(column.concentrationDot[25 * n + 2], -6.1706769254940182e-08, 1e-14);
  EXPECT_NEAR(column.concentrationDot[50 * n + 5], 1.0247150372002965e-07, 1e-14);
}

TEST(ExampleMixturePrediction, RecomputesBeaBreakthroughEquilibriumLoadings)
{
  Column column = loadBeaBreakthroughColumn();
  const std::vector<double> savedEquilibriumAdsorption = column.equilibriumPhysisorption;

  RK3Helpers::computeEquilibriumLoadings(column);

  EXPECT_LE(maxAbsDifference(column.equilibriumPhysisorption, savedEquilibriumAdsorption), 1e-10);

  const size_t n = column.numberOfComponents;
  EXPECT_DOUBLE_EQ(column.equilibriumPhysisorption[0], 0.0);
  EXPECT_NEAR(column.equilibriumPhysisorption[0 * n + 1], 0.070539163250492337, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[10 * n + 2], 0.046853375510025136, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[25 * n + 3], 0.031885705858180735, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[50 * n + 5], 0.021244052749882759, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[50 * n + 7], 0.0066257061765186923, 1e-14);
}

TEST(ChemisorptionGeneral, DefaultsToBulkConcentrationWithoutPoreSurfaceTransport)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_general_chem_direct_test.json";
  writeGeneralChemisorptionInput(path, false);

  InputReader reader(path.string());
  ASSERT_EQ(reader.components.size(), 2U);
  ASSERT_EQ(reader.components[1].chemisorption.numberOfSites, 1U);
  const Chemisorption& kinetics = reader.components[1].chemisorption.sites[0];
  EXPECT_EQ(kinetics.type, Chemisorption::Type::General);
  EXPECT_DOUBLE_EQ(kinetics.maximumLoading, 2.0);
  EXPECT_DOUBLE_EQ(kinetics.adsorptionRateCoefficient, 1.0e-4);
  EXPECT_FALSE(kinetics.usePoreSurfaceTransport);
  EXPECT_FALSE(reader.components[1].chemisorption.usesSurfacePoreTransport());

  Column column(reader);
  column.initialize();

  EXPECT_FALSE(column.surfacePoreTransportEnabled);
  EXPECT_EQ(column.stateSize(), (3 * column.numberOfComponents + 3) * (column.numberOfGridPoints + 1));
  EXPECT_TRUE(column.surfaceConcentration.empty());
  EXPECT_TRUE(column.poreConcentration.empty());

  RK3Helpers::computeSorptionDerivatives(column);
  const size_t adsorbateInlet = 1U;
  EXPECT_GT(column.chemisorptionDot[adsorbateInlet], 0.0);
  EXPECT_TRUE(column.poreConcentrationDot.empty());
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet], 0.0);
}

TEST(ChemisorptionGeneral, InitializesSurfacePoreStateAndFilmSinkWhenEnabled)
{
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ruptura_general_chem_transport_test.json";
  writeGeneralChemisorptionInput(path, true);

  InputReader reader(path.string());
  ASSERT_EQ(reader.components.size(), 2U);
  ASSERT_EQ(reader.components[1].chemisorption.numberOfSites, 1U);
  EXPECT_TRUE(reader.components[1].chemisorption.sites[0].usePoreSurfaceTransport);
  EXPECT_TRUE(reader.components[1].chemisorption.usesSurfacePoreTransport());

  Column column(reader);
  column.initialize();

  EXPECT_TRUE(column.surfacePoreTransportEnabled);
  EXPECT_EQ(column.stateSize(), (5 * column.numberOfComponents + 3) * (column.numberOfGridPoints + 1));
  EXPECT_EQ(column.surfaceConcentration.size(), column.concentration.size());
  EXPECT_EQ(column.poreConcentration.size(), column.concentration.size());
  EXPECT_LE(maxAbsDifference(column.surfaceConcentration, toVector(column.concentration)), 0.0);
  EXPECT_LE(maxAbsDifference(column.poreConcentration, toVector(column.concentration)), 0.0);

  RK3Helpers::computeSorptionDerivatives(column);
  const size_t adsorbateInlet = 1U;
  EXPECT_GT(column.chemisorptionDot[adsorbateInlet], 0.0);
  const double poreLoadingConversion = column.geometry.loadingPrefactor(column.particleDensity);
  EXPECT_NEAR(column.poreConcentrationDot[adsorbateInlet],
              -poreLoadingConversion * column.chemisorptionDot[adsorbateInlet], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[adsorbateInlet], 0.0, 1.0e-12);

  column.surfaceConcentration[adsorbateInlet] = 0.5 * column.concentration[adsorbateInlet];
  RK3Helpers::computeBulkSpeciesSink(column);
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet], 0.0);
}

TEST(Arrhenius, SupportsAbsoluteAndReferenceTemperatureForms)
{
  constexpr double preExponentialFactor = 2.5;
  constexpr double activationEnergy = 12000.0;
  constexpr double temperature = 350.0;
  constexpr double referenceTemperature = 300.0;

  EXPECT_NEAR(arrhenius(preExponentialFactor, activationEnergy, temperature),
              preExponentialFactor * std::exp(-activationEnergy / (R * temperature)), 1.0e-14);
  EXPECT_NEAR(arrhenius(preExponentialFactor, activationEnergy, temperature, referenceTemperature),
              preExponentialFactor * std::exp(-activationEnergy / R * (1.0 / temperature - 1.0 / referenceTemperature)),
              1.0e-14);

  Component component(0, "test");
  component.nonIsothermal = true;
  component.heatOfAdsorption = 12000.0;
  component.referenceTemperature = referenceTemperature;
  EXPECT_NEAR(component.scale(temperature),
              arrhenius(1.0, -component.heatOfAdsorption, temperature, component.referenceTemperature), 1.0e-14);
}

TEST(ChemisorptionRate, GeneralAndElovichUseUnifiedRate)
{
  Chemisorption general;
  general.type = Chemisorption::Type::General;
  general.maximumLoading = 3.0;
  general.adsorptionRateCoefficient = 2.0;
  general.desorptionRateCoefficient = 0.5;
  general.poreConcentrationOrder = 2;
  general.capacityOrder = 1;
  general.desorptionOrder = 1;
  EXPECT_DOUBLE_EQ(general.rate(3.0, 1.0, 4.0, 300.0), 63.5);

  Chemisorption elovich;
  elovich.type = Chemisorption::Type::Elovich;
  elovich.elovichAlpha = 2.0;
  elovich.elovichBeta = 0.5;
  EXPECT_NEAR(elovich.rate(0.0, 4.0, 3.0), 6.0 * std::exp(-2.0), 1.0e-14);
}

TEST(ReactionsInput, ParsesGeneralPowerLawReaction)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_parse.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  ASSERT_EQ(reader.reactions.size(), 1U);
  const Reaction& reaction = reader.reactions[0];
  EXPECT_EQ(reaction.phase, Reaction::Phase::Physisorbed);
  EXPECT_EQ(reaction.style, Reaction::Style::GeneralPowerLaw);
  EXPECT_EQ(reaction.reactants, std::vector<size_t>{1});
  EXPECT_EQ(reaction.products, std::vector<size_t>{2});
  EXPECT_DOUBLE_EQ(reaction.reactantStoichiometry[0], 1.0);
  EXPECT_DOUBLE_EQ(reaction.productStoichiometry[0], 1.0);
  EXPECT_DOUBLE_EQ(reaction.forwardRateCoefficient, 1.0);
  EXPECT_DOUBLE_EQ(reaction.equilibriumConstant, 1.0e9);
  EXPECT_DOUBLE_EQ(reaction.gibbsFreeEnergy, -1000.0);
}

TEST(ReactionsInput, ParsesParticipantNamesAsComponentIndices)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_names.json";
  writeReactionInput(path, "Physisorbed", "GeneralPowerLaw", true, R"json(["A"])json", R"json(["B"])json");

  InputReader reader(path.string());
  ASSERT_EQ(reader.reactions.size(), 1U);
  const Reaction& reaction = reader.reactions[0];
  EXPECT_EQ(reaction.reactants, std::vector<size_t>{1});
  EXPECT_EQ(reaction.products, std::vector<size_t>{2});
}

TEST(ReactionsAdsorbed, AddsAdsorbedSourcesWithoutBulkTransferAndLimitsQmax)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_adsorbed.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  const size_t reactant = 1U;
  const size_t product = 2U;
  column.physisorption[reactant] = 0.4;
  column.physisorption[product] = 0.0;

  RK3Helpers::computeSorptionDerivatives(column);
  EXPECT_LT(column.physisorptionDot[reactant], 0.0);
  EXPECT_GT(column.physisorptionDot[product], 0.0);
  EXPECT_NEAR(column.physisorptionDot[reactant], -column.physisorptionDot[product], 1.0e-12);
  EXPECT_NEAR(column.reactionPhysisorptionSource[reactant], column.physisorptionDot[reactant], 1.0e-12);
  EXPECT_NEAR(column.reactionPhysisorptionSource[product], column.physisorptionDot[product], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);

  RK3Helpers::computeEnergyDerivatives(column);
  EXPECT_NEAR(column.solidTemperatureDot[0], column.reactionHeat[0] / column.heatCapacitySolid, 1.0e-12);
  EXPECT_GT(column.solidTemperatureDot[0], 0.0);

  column.physisorption[product] = 1.0;
  RK3Helpers::computeSorptionDerivatives(column);
  EXPECT_NEAR(column.reactionPhysisorptionSource[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.reactionPhysisorptionSource[product], 0.0, 1.0e-12);
}

TEST(ReactionsPore, LangmuirHinshelwoodRequiresAndUsesPoreConcentration)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_lh_pore.json";
  writeReactionInput(path, "PoreConcentration", "LangmuirHinshelwood", false);

  InputReader reader(path.string());
  ASSERT_EQ(reader.reactions.size(), 1U);
  EXPECT_TRUE(reactionsRequirePoreConcentration(reader.reactions));

  Column column(reader);
  column.initialize();
  EXPECT_TRUE(column.surfacePoreTransportEnabled);
  EXPECT_EQ(column.poreConcentration.size(), column.concentration.size());

  RK3Helpers::computeSorptionDerivatives(column);
  const size_t reactant = 1U;
  const size_t product = 2U;
  EXPECT_LT(column.poreConcentrationDot[reactant], 0.0);
  EXPECT_GT(column.poreConcentrationDot[product], 0.0);
  EXPECT_LT(column.reactionPoreConcentrationSource[reactant], 0.0);
  EXPECT_GT(column.reactionPoreConcentrationSource[product], 0.0);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);
}

TEST(ReactionsAutoStop, RequiresLowReactionSourcesAndStablePhaseStates)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_auto_stop.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  const size_t reactant = 1U;
  const size_t product = 2U;
  column.physisorption[reactant] = 0.4;
  RK3Helpers::computeSorptionDerivatives(column);
  EXPECT_FALSE(RK3Helpers::reactionAutoStopReached(column, 0.01));

  column.physisorption[product] = 1.0;
  RK3Helpers::computeSorptionDerivatives(column);
  EXPECT_TRUE(RK3Helpers::reactionAutoStopReached(column, 0.01));

  column.physisorptionDot[reactant] = 1.0;
  EXPECT_FALSE(RK3Helpers::reactionAutoStopReached(column, 0.01));
}

TEST(ReactionsAutoStop, RungeKuttaSchedulesAutoStopForStableReactionState)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_auto_stop_rk3.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();
  column.physisorption[2U] = 1.0;

  RungeKutta3 integrator(0.01, true, 0);
  Timing timings;
  EXPECT_FALSE(integrator.propagate(column, 0, timings));
  EXPECT_FALSE(integrator.autoNumberOfSteps);
  EXPECT_EQ(integrator.numberOfSteps, 2U);
}

TEST(ReactionsInput, RejectsLHOutsidePoreConcentration)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_reaction_lh_bad_phase.json";
  writeReactionInput(path, "Physisorbed", "LangmuirHinshelwood", false);

  EXPECT_THROW(
      {
        try
        {
          static_cast<void>(InputReader(path.string()));
        }
        catch (const std::runtime_error& error)
        {
          EXPECT_NE(std::string(error.what()).find("LH/LHHW"), std::string::npos);
          throw;
        }
      },
      std::runtime_error);
}

TEST(ChemisorptionInput, RejectsFractionalOrders)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_fractional_chem_order.json";
  writeGeneralChemisorptionInput(path, false, "1.5");

  EXPECT_THROW(
      {
        try
        {
          static_cast<void>(InputReader(path.string()));
        }
        catch (const std::runtime_error& error)
        {
          EXPECT_NE(std::string(error.what()).find("non-negative integer"), std::string::npos);
          throw;
        }
      },
      std::runtime_error);
}

TEST(ChemisorptionInput, RequiresExactlyTheParametersForTheSelectedModel)
{
  const std::filesystem::path validPath = std::filesystem::temp_directory_path() / "ruptura_strict_chem_valid.json";
  writeGeneralChemisorptionInput(validPath, false);

  std::ifstream input(validPath);
  const nlohmann::json valid = nlohmann::json::parse(input);
  input.close();

  auto expectRejected = [&](nlohmann::json data, const std::string& suffix, const std::string& expectedMessage)
  {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / ("ruptura_strict_chem_" + suffix + ".json");
    std::ofstream output(path);
    output << data;
    output.close();

    EXPECT_THROW(
        {
          try
          {
            static_cast<void>(InputReader(path.string()));
          }
          catch (const std::runtime_error& error)
          {
            EXPECT_NE(std::string(error.what()).find(expectedMessage), std::string::npos);
            throw;
          }
        },
        std::runtime_error);
  };

  nlohmann::json missing = valid;
  missing["Components"][1]["ChemisorptionSites"][0]["Parameters"].erase("capacityOrder");
  expectRejected(std::move(missing), "missing", "required key 'capacityOrder'");

  nlohmann::json extra = valid;
  extra["Components"][1]["ChemisorptionSites"][0]["Parameters"]["alpha"] = 1.0;
  expectRejected(std::move(extra), "extra", "unknown key 'alpha'");

  nlohmann::json oldAlias = valid;
  auto& oldAliasParameters = oldAlias["Components"][1]["ChemisorptionSites"][0]["Parameters"];
  oldAliasParameters["Ea_ka"] = oldAliasParameters["adsorptionActivationEnergy"];
  oldAliasParameters.erase("adsorptionActivationEnergy");
  expectRejected(std::move(oldAlias), "old_alias", "unknown key 'Ea_ka'");

  nlohmann::json explicitCompleteMixture = valid;
  explicitCompleteMixture["MixturePredictionMethod"] = "EI";
  expectRejected(std::move(explicitCompleteMixture), "ei", "EI is not supported");
}

TEST(ChemisorptionInput, LoadsUpdatedSurfacePoreExample)
{
  InputReader reader(examplePath("examples/Chemisorption-CO2/breakthrough/simulation.json"));
  ASSERT_EQ(reader.components[1].chemisorption.numberOfSites, 2U);
  EXPECT_EQ(reader.components[1].chemisorption.sites[0].type, Chemisorption::Type::General);
  EXPECT_EQ(reader.components[1].chemisorption.sites[1].type, Chemisorption::Type::Elovich);
  EXPECT_TRUE(reader.components[1].chemisorption.sites[0].isotherm.has_value());
  EXPECT_TRUE(reader.components[1].chemisorption.sites[1].isotherm.has_value());
  EXPECT_TRUE(reader.components[1].chemisorption.usesSurfacePoreTransport());
}

TEST(ChemisorptionEquilibrium, ComputesCompetitiveSiteLoadingsSeparatelyFromPhysisorption)
{
  std::ifstream input(examplePath("examples/Chemisorption-CO2-N2/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();

  auto& nitrogen = data["Components"][2];
  nitrogen["ChemisorptionSites"] = nlohmann::json::array({
      {
          {"Type", "FirstOrder"},
          {"Parameters",
           {
               {"rateCoefficient", 0.02},
               {"maximumLoading", 0.40},
               {"heatOfChemisorption", 30000.0},
               {"Isotherm", {{"Type", "Langmuir"}, {"Parameters", {0.40, 1.0e-5}}}},
           }},
      },
  });

  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ruptura_competitive_chem_equilibrium.json";
  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  const size_t carbonDioxide = 1;
  const size_t nitrogenIndex = 2;
  const double pureCarbonDioxide = reader.components[carbonDioxide].chemisorption.sites[0].isotherm->value(
      column.partialPressure[carbonDioxide], 1.0);
  const double pureNitrogen = reader.components[nitrogenIndex].chemisorption.sites[0].isotherm->value(
      column.partialPressure[nitrogenIndex], 1.0);

  EXPECT_GT(column.equilibriumChemisorption[carbonDioxide], 0.0);
  EXPECT_GT(column.equilibriumChemisorption[nitrogenIndex], 0.0);
  EXPECT_LT(column.equilibriumChemisorption[carbonDioxide], pureCarbonDioxide);
  EXPECT_LT(column.equilibriumChemisorption[nitrogenIndex], pureNitrogen);
  EXPECT_NE(column.equilibriumChemisorption[carbonDioxide], column.equilibriumPhysisorption[carbonDioxide]);

  auto secondCarbonDioxideSite = data["Components"][1]["ChemisorptionSites"][0];
  secondCarbonDioxideSite["Parameters"]["maximumLoading"] = 0.20;
  secondCarbonDioxideSite["Parameters"]["Isotherm"]["Parameters"] = {0.20, 1.5e-5};
  data["Components"][1]["ChemisorptionSites"].push_back(secondCarbonDioxideSite);

  const std::filesystem::path missingSitePath =
      std::filesystem::temp_directory_path() / "ruptura_component_missing_chem_site.json";
  std::ofstream missingSiteOutput(missingSitePath);
  missingSiteOutput << data;
  missingSiteOutput.close();

  InputReader missingSiteReader(missingSitePath.string());
  Column missingSiteColumn(missingSiteReader);
  missingSiteColumn.initialize();
  const size_t secondSiteOffset = missingSiteColumn.concentration.size();
  EXPECT_EQ(missingSiteColumn.chemisorptionMixture.predictionMethod, MixturePrediction::PredictionMethod::IAST);
  EXPECT_EQ(missingSiteColumn.chemisorptionMixture.maxIsothermTerms, 2U);
  EXPECT_GT(missingSiteColumn.equilibriumChemisorption[secondSiteOffset + carbonDioxide], 0.0);
  EXPECT_DOUBLE_EQ(missingSiteColumn.equilibriumChemisorption[secondSiteOffset + nitrogenIndex], 0.0);

  const double completeMixtureSecondSite = missingSiteColumn.equilibriumChemisorption[secondSiteOffset + carbonDioxide];

  data["MixturePredictionMethod"] = "IAST";
  data["IASTMethod"] = "NestedLoopBisection";
  const std::filesystem::path completeBisectionPath =
      std::filesystem::temp_directory_path() / "ruptura_complete_chem_iast_bisection.json";
  std::ofstream completeBisectionOutput(completeBisectionPath);
  completeBisectionOutput << data;
  completeBisectionOutput.close();

  InputReader completeBisectionReader(completeBisectionPath.string());
  Column completeBisectionColumn(completeBisectionReader);
  completeBisectionColumn.initialize();
  EXPECT_EQ(completeBisectionColumn.chemisorptionMixture.predictionMethod, MixturePrediction::PredictionMethod::IAST);
  EXPECT_EQ(completeBisectionColumn.chemisorptionMixture.iastMethod,
            MixturePrediction::IASTMethod::NestedLoopBisection);
  EXPECT_GT(completeBisectionColumn.equilibriumChemisorption[carbonDioxide], 0.0);
  EXPECT_GT(completeBisectionColumn.equilibriumChemisorption[nitrogenIndex], 0.0);
  EXPECT_GT(completeBisectionColumn.equilibriumChemisorption[secondSiteOffset + carbonDioxide], 0.0);
  EXPECT_DOUBLE_EQ(completeBisectionColumn.equilibriumChemisorption[secondSiteOffset + nitrogenIndex], 0.0);

  const std::vector<std::pair<std::string, std::string>> segregatedMethods{{"SIAST", "FastIAST"},
                                                                           {"SIAST", "NestedLoopBisection"},
                                                                           {"SEI", "FastIAST"},
                                                                           {"SCI", "FastIAST"},
                                                                           {"SPI", "FastIAST"}};
  for (const auto& [method, iastMethod] : segregatedMethods)
  {
    data["MixturePredictionMethod"] = method;
    data["IASTMethod"] = iastMethod;
    const std::filesystem::path segregatedPath =
        std::filesystem::temp_directory_path() / ("ruptura_segregated_chem_" + method + "_" + iastMethod + ".json");
    std::ofstream segregatedOutput(segregatedPath);
    segregatedOutput << data;
    segregatedOutput.close();

    InputReader segregatedReader(segregatedPath.string());
    Column segregatedColumn(segregatedReader);
    segregatedColumn.initialize();

    MixturePrediction::PredictionMethod expectedMethod = MixturePrediction::PredictionMethod::SEI;
    if (method == "SIAST") expectedMethod = MixturePrediction::PredictionMethod::SIAST;
    if (method == "SCI") expectedMethod = MixturePrediction::PredictionMethod::SCI;
    if (method == "SPI") expectedMethod = MixturePrediction::PredictionMethod::SPI;
    EXPECT_EQ(segregatedColumn.chemisorptionMixture.predictionMethod, expectedMethod);
    EXPECT_GT(segregatedColumn.equilibriumChemisorption[carbonDioxide], 0.0);
    EXPECT_GT(segregatedColumn.equilibriumChemisorption[nitrogenIndex], 0.0);
    EXPECT_GT(segregatedColumn.equilibriumChemisorption[secondSiteOffset + carbonDioxide], 0.0);
    EXPECT_DOUBLE_EQ(segregatedColumn.equilibriumChemisorption[secondSiteOffset + nitrogenIndex], 0.0);
    EXPECT_NE(segregatedColumn.equilibriumChemisorption[secondSiteOffset + carbonDioxide], completeMixtureSecondSite);
  }
}

TEST(ChemisorptionElovich, ParsesAndSubtractsMassFromBulk)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_elovich_chem_test.json";
  writeElovichChemisorptionInput(path, false);

  InputReader reader(path.string());
  ASSERT_EQ(reader.components[1].chemisorption.numberOfSites, 1U);
  const Chemisorption& kinetics = reader.components[1].chemisorption.sites[0];
  EXPECT_EQ(kinetics.type, Chemisorption::Type::Elovich);
  EXPECT_DOUBLE_EQ(kinetics.elovichAlpha, 0.25);
  EXPECT_DOUBLE_EQ(kinetics.elovichBeta, 0.5);
  EXPECT_FALSE(kinetics.usesSurfacePoreTransport());

  Column column(reader);
  column.initialize();
  RK3Helpers::computeSorptionDerivatives(column);

  const size_t adsorbateInlet = 1U;
  EXPECT_NEAR(column.chemisorptionDot[adsorbateInlet], 0.25 * column.concentration[adsorbateInlet], 1.0e-12);
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet], 0.0);
}

TEST(MultiSiteChemisorption, ParsesIndependentSitesAndSumsTheirBulkSink)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_multisite_chem_direct_test.json";
  writeMultiSiteChemisorptionInput(path, false);

  InputReader reader(path.string());
  const MultiSiteChemisorption& kinetics = reader.components[1].chemisorption;
  ASSERT_EQ(kinetics.numberOfSites, 2U);
  ASSERT_EQ(kinetics.sites.size(), 2U);
  EXPECT_EQ(kinetics.sites[0].type, Chemisorption::Type::General);
  EXPECT_EQ(kinetics.sites[1].type, Chemisorption::Type::Elovich);
  EXPECT_DOUBLE_EQ(kinetics.maximumLoading(), 3.0);
  ASSERT_TRUE(kinetics.sites[0].isotherm.has_value());
  ASSERT_TRUE(kinetics.sites[1].isotherm.has_value());

  Column column(reader);
  column.initialize();
  ASSERT_EQ(column.maxChemisorptionSites, 2U);
  EXPECT_FALSE(column.surfacePoreTransportEnabled);
  EXPECT_EQ(column.chemisorption.size(), 2 * column.concentration.size());
  EXPECT_EQ(column.stateSize(), (4 * column.numberOfComponents + 3) * (column.numberOfGridPoints + 1));

  RK3Helpers::computeSorptionDerivatives(column);

  const size_t componentBlockSize = column.concentration.size();
  const size_t adsorbateInlet = 1U;
  const size_t secondSiteInlet = componentBlockSize + adsorbateInlet;
  ASSERT_GT(column.equilibriumChemisorption[adsorbateInlet], 0.0);
  ASSERT_GT(column.equilibriumChemisorption[secondSiteInlet], 0.0);
  EXPECT_NEAR(column.chemisorptionDot[adsorbateInlet],
              1.0e-4 * column.concentration[adsorbateInlet] * column.equilibriumChemisorption[adsorbateInlet], 1.0e-12);
  EXPECT_NEAR(column.chemisorptionDot[secondSiteInlet], 0.25 * column.concentration[adsorbateInlet], 1.0e-12);

  const double loadingPrefactor = column.geometry.loadingPrefactor(column.particleDensity);
  EXPECT_NEAR(column.bulkSpeciesSink[adsorbateInlet],
              loadingPrefactor * (column.physisorptionDot[adsorbateInlet] + column.chemisorptionDot[adsorbateInlet] +
                                  column.chemisorptionDot[secondSiteInlet]),
              1.0e-12);
}

TEST(MultiSiteChemisorption, KeepsDirectAndPoreTransportSiteBalancesIndependent)
{
  const std::filesystem::path path =
      std::filesystem::temp_directory_path() / "ruptura_multisite_chem_mixed_transport_test.json";
  writeMultiSiteChemisorptionInput(path, true);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  const size_t componentBlockSize = column.concentration.size();
  const size_t adsorbateInlet = 1U;
  const size_t secondSiteInlet = componentBlockSize + adsorbateInlet;
  ASSERT_TRUE(column.surfacePoreTransportEnabled);
  ASSERT_EQ(column.surfaceConcentration.size(), 2 * componentBlockSize);
  ASSERT_EQ(column.poreConcentration.size(), 2 * componentBlockSize);
  EXPECT_EQ(column.stateSize(), (8 * column.numberOfComponents + 3) * (column.numberOfGridPoints + 1));

  RK3Helpers::computeSorptionDerivatives(column);

  EXPECT_GT(column.chemisorptionDot[adsorbateInlet], 0.0);
  EXPECT_GT(column.chemisorptionDot[secondSiteInlet], 0.0);
  EXPECT_DOUBLE_EQ(column.surfaceConcentrationDot[adsorbateInlet], 0.0);
  EXPECT_DOUBLE_EQ(column.poreConcentrationDot[adsorbateInlet], 0.0);

  const double loadingPrefactor = column.geometry.loadingPrefactor(column.particleDensity);
  EXPECT_NEAR(column.poreConcentrationDot[secondSiteInlet],
              -loadingPrefactor * column.chemisorptionDot[secondSiteInlet], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[adsorbateInlet],
              loadingPrefactor * (column.physisorptionDot[adsorbateInlet] + column.chemisorptionDot[adsorbateInlet]),
              1.0e-12);

  column.surfaceConcentration[secondSiteInlet] = 0.5 * column.concentration[adsorbateInlet];
  RK3Helpers::computeBulkSpeciesSink(column);
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet],
            loadingPrefactor * (column.physisorptionDot[adsorbateInlet] + column.chemisorptionDot[adsorbateInlet]));

  const std::filesystem::path restartPath =
      std::filesystem::temp_directory_path() / "ruptura_multisite_chem_restart_test.json";
  column.writeJSON(restartPath.string());
  Column restored(reader);
  restored.initialize();
  restored.readJSON(restartPath.string());
  EXPECT_EQ(toVector(restored.chemisorption), toVector(column.chemisorption));
  EXPECT_EQ(toVector(restored.equilibriumChemisorption), toVector(column.equilibriumChemisorption));
  EXPECT_EQ(toVector(restored.surfaceConcentration), toVector(column.surfaceConcentration));
  EXPECT_EQ(toVector(restored.poreConcentration), toVector(column.poreConcentration));
}

#if BUILD_SUNDIALS
TEST(MultiSiteChemisorption, PropagatesMixedSitesWithCvode)
{
  const std::filesystem::path path = std::filesystem::temp_directory_path() / "ruptura_multisite_chem_cvode_test.json";
  writeMultiSiteChemisorptionInput(path, true);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  CVODE integrator(1.0e-5, false, 1);
  Timing timings;
  integrator.initialize(column);
  EXPECT_NO_THROW(integrator.propagate(column, 0, timings));

  const size_t adsorbateInlet = 1U;
  const size_t secondSiteInlet = column.concentration.size() + adsorbateInlet;
  EXPECT_GT(column.chemisorption[adsorbateInlet], 0.0);
  EXPECT_GT(column.chemisorption[secondSiteInlet], 0.0);
  EXPECT_GE(column.poreConcentration[secondSiteInlet], 0.0);
}
#endif

TEST(ExampleFitting, ScoresKnownBeaLangmuirFitWithoutRunningOptimizer)
{
  InputReader reader(examplePath("examples/BEA-alkanes-C7/fitting/simulation.json"));
  Fitting fitting(reader);
  fitting.filename[0] = examplePath("examples/BEA-alkanes-C7/fitting/Results.dat-BEA-Repeat-552K-nC7");
  fitting.readData(0);

  ASSERT_EQ(fitting.rawData.size(), 12U);
  EXPECT_DOUBLE_EQ(fitting.pressureRange.first, 10.0);
  EXPECT_DOUBLE_EQ(fitting.pressureRange.second, 3000000.0);
  EXPECT_NEAR(fitting.maximumLoading, 1.235003511, 1e-12);

  MultiSiteIsotherm fittedIsotherm({
      Isotherm(Isotherm::Type::Langmuir, {0.19466024781253777, 8.90730767575818e-07}, false),
      Isotherm(Isotherm::Type::Langmuir, {1.0998368974010369, 6.558568046086507e-05}, false),
  });

  EXPECT_NEAR(fitting.fitness(fittedIsotherm), 0.005825848539239016, 1e-14);
  EXPECT_GT(fitting.RCorrelation(fittedIsotherm), 0.9999);
}

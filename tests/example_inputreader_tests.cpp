#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numbers>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

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
#include "temporary_directory.h"

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
  // Column state after 100 s of RungeKutta3 on this example (Δt = 1e-3 s); running the example overwrites its own
  // column.json, so the snapshot is kept with the tests.
  column.readJSON(examplePath("tests/data/bea_breakthrough_column.json"));
  column.convectionScheme = Column::ConvectionScheme::Upwind;  // the snapshot and expectedMassDerivative are upwind
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
    "Type": "PackedBed",
    "ColumnVoidFraction": 0.4,
    "ParticleDiameter": 1.0e-3
  },
  "ParticleDensity": 1000.0,
  "InletPressure": 100000.0,
  "ColumnEntranceVelocity": 0.1,
  "ColumnLength": 0.1,
  "NumberOfGridPoints": 2,
  "SimulationTime": 0.01,
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
        "Type": "Langmuir",
        "Parameters": [2.0, 1.0e-5],
        "RateEquation": "General",
        "RateEquationParameters": {
          "AdsorptionRateCoefficient": 1.0e-4,
          "AdsorptionActivationEnergy": 0.0,
          "DesorptionRateCoefficient": 0.0,
          "DesorptionActivationEnergy": 0.0,
          "PoreConcentrationOrder": )json"
      << concentrationOrder << R"json(,
          "CapacityOrder": 1,
          "DesorptionOrder": 1,
          "FilmMassTransferCoefficient": 2.0e-3,
          "PoreDiffusivity": 1.0e-10,
          "UsePoreSurfaceTransport": )json"
      << (usePoreSurfaceTransport ? "true" : "false") << R"json(
        },
        "HeatOfAdsorption": 42000.0,
        "MaximumLoading": 2.0
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
  chemisorption["RateEquation"] = "Elovich";
  chemisorption["RateEquationParameters"] = {
      {"Alpha", 0.25},
      {"Beta", 0.5},
      {"FilmMassTransferCoefficient", 2.0e-3},
      {"PoreDiffusivity", 1.0e-10},
      {"UsePoreSurfaceTransport", usePoreSurfaceTransport},
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
  general["MaximumLoading"] = 1.0;
  general["Parameters"] = {1.0, 1.0e-5};
  nlohmann::json elovich = {
      {"Type", "Langmuir"},
      {"Parameters", {2.0, 2.0e-5}},
      {"RateEquation", "Elovich"},
      {"RateEquationParameters",
       {
           {"Alpha", 0.25},
           {"Beta", 0.5},
           {"FilmMassTransferCoefficient", 2.0e-3},
           {"PoreDiffusivity", 1.0e-10},
           {"UsePoreSurfaceTransport", secondSiteUsesPoreSurfaceTransport},
       }},
      {"HeatOfAdsorption", 42000.0},
      {"MaximumLoading", 2.0},
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
    "Type": "PackedBed",
    "ColumnVoidFraction": 0.4,
    "ParticleDiameter": 1.0e-3
  },
  "ParticleDensity": 1000.0,
  "InletPressure": 100000.0,
  "ColumnEntranceVelocity": 0.1,
  "ColumnLength": 0.1,
  "NumberOfGridPoints": 2,
  "SimulationTime": 0.01,
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

TEST(ExampleInputReader, LoadsCompleteBeaMpdBreakthroughBenchmarkMatrix)
{
  const std::filesystem::path benchmarkDirectory =
      std::filesystem::path(RUPTURA_SOURCE_DIR) / "examples" / "MPD" / "BEA-alkanes" / "breakthrough";
  std::set<std::tuple<size_t, size_t, size_t, int>> combinations;

  for (const std::filesystem::directory_entry& entry : std::filesystem::directory_iterator(benchmarkDirectory))
  {
    const std::filesystem::path inputPath = entry.path() / "simulation.json";
    if (!entry.is_directory() || !std::filesystem::is_regular_file(inputPath)) continue;

    InputReader reader(inputPath.string());
    EXPECT_EQ(reader.simulationType, InputReader::SimulationType::Breakthrough);
    EXPECT_EQ(reader.components.size(), 8U);
    EXPECT_EQ(reader.numberOfCarrierGases, 1U);
    EXPECT_FALSE(reader.autoSimulationTime);
    EXPECT_GT(reader.simulationTime, 0.0);
    EXPECT_TRUE(reader.breakthroughIntegrator == 0U || reader.breakthroughIntegrator == 1U);
    EXPECT_TRUE(reader.numberOfGridPoints == 20U || reader.numberOfGridPoints == 50U ||
                reader.numberOfGridPoints == 100U);

    const int velocityMilliMetresPerSecond = static_cast<int>(std::lround(1000.0 * reader.columnEntranceVelocity));
    EXPECT_TRUE(velocityMilliMetresPerSecond == 19 || velocityMilliMetresPerSecond == 190);
    const auto [_, inserted] = combinations.emplace(reader.mixturePredictionMethod, reader.breakthroughIntegrator,
                                                     reader.numberOfGridPoints, velocityMilliMetresPerSecond);
    EXPECT_TRUE(inserted) << "duplicate benchmark combination in " << inputPath;
  }

  EXPECT_EQ(combinations.size(), 24U);
  for (size_t method : {static_cast<size_t>(MixturePrediction::PredictionMethod::SCI),
                        static_cast<size_t>(MixturePrediction::PredictionMethod::MPD)})
  {
    for (size_t integrator : {0U, 1U})
    {
      for (size_t gridPoints : {20U, 50U, 100U})
      {
        for (int velocity : {19, 190})
        {
          EXPECT_TRUE(combinations.contains({method, integrator, gridPoints, velocity}));
        }
      }
    }
  }
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
  EXPECT_TRUE(reader.autoSimulationTime);
  EXPECT_DOUBLE_EQ(reader.writeInterval, 50.0);
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
  data["SimulationTime"] = 5.0e-4;

  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_sirk3_integrator.json";
  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  EXPECT_EQ(reader.breakthroughIntegrator, 2U);
  EXPECT_EQ(static_cast<BreakthroughIntegrationScheme>(reader.breakthroughIntegrator),
            BreakthroughIntegrationScheme::SIRK3);

  SemiImplicitRungeKutta3 integrator(reader);
  EXPECT_DOUBLE_EQ(integrator.timeStep, reader.timeStep);
  EXPECT_EQ(integrator.numberOfInitTimeSteps, reader.numberOfInitTimeSteps);

  Column column(reader);
  column.initialize();
  Timing timings;
  integrator.advance(column, 0.0, reader.simulationTime, timings);
  EXPECT_EQ(integrator.numberOfStepsTaken, 1U);
}

TEST(ExampleInputReader, Sirk3RejectsReactionsAndSurfacePoreTransport)
{
  auto expectSirk3Rejected = [](const std::filesystem::path& path, const std::string& expectedMessage)
  {
    std::ifstream input(path);
    nlohmann::json data = nlohmann::json::parse(input);
    input.close();
    data["BreakthroughIntegrator"] = "SIRK3";
    std::ofstream(path) << data;

    EXPECT_THROW(
        {
          try
          {
            static_cast<void>(InputReader(path.string()));
          }
          catch (const std::runtime_error& error)
          {
            EXPECT_NE(std::string(error.what()).find(expectedMessage), std::string::npos) << error.what();
            throw;
          }
        },
        std::runtime_error);
  };

  const std::filesystem::path reactionPath = testTemporaryDirectory() / "ruptura_sirk3_physisorbed_reaction.json";
  writeReactionInput(reactionPath, "Physisorbed", "GeneralPowerLaw", false);
  expectSirk3Rejected(reactionPath, "does not evaluate reactions");

  const std::filesystem::path porePath = testTemporaryDirectory() / "ruptura_sirk3_pore_transport.json";
  writeGeneralChemisorptionInput(porePath, true);
  expectSirk3Rejected(porePath, "does not evaluate surface-pore transport");
}

TEST(ExampleInputReader, ParsesAndValidatesCvodeMaximumOrder)
{
  std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_cvode_maximum_order.json";

  std::ofstream(path) << data;
  EXPECT_EQ(InputReader(path.string()).cvodeMaximumOrder, 5U);

  data["CVODEMaximumOrder"] = 3;
  std::ofstream(path) << data;
  const InputReader reader(path.string());
  EXPECT_EQ(reader.cvodeMaximumOrder, 3U);
#if BUILD_SUNDIALS
  EXPECT_EQ(CVODE(reader).maximumOrder, 3);
#endif

  for (int order : {0, 6})
  {
    data["CVODEMaximumOrder"] = order;
    std::ofstream(path) << data;
    EXPECT_THROW(InputReader{path.string()}, std::runtime_error);
  }
}

TEST(SwingAdsorptionInput, ParsesStructuredPhasesAndBuildsSubStages)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_swing_structured_test.json";
  writeGeneralChemisorptionInput(path, false);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();

  data["SimulationType"] = "SwingAdsorption";
  data.erase("InletPressure");
  data.erase("SimulationTime");
  data["SwingAdsorptionPhases"] = nlohmann::json::array({
      {{"Name", "Adsorption"}, {"Temperature", 310.0}, {"InletPressure", 200000.0}, {"Duration", 0.03}},
      {{"Name", "Depressurization"}, {"InletPressure", 50000.0}, {"Duration", 0.02}},
  });

  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  EXPECT_EQ(reader.simulationType, InputReader::SimulationType::SwingAdsorption);
  EXPECT_FALSE(reader.autoSimulationTime);
  EXPECT_DOUBLE_EQ(reader.simulationTime, 0.05);
  EXPECT_DOUBLE_EQ(reader.inletPressure, 200000.0);

  ASSERT_EQ(reader.swingAdsorptionPhases.size(), 2U);
  EXPECT_EQ(reader.swingAdsorptionPhases[0].name, "Adsorption");
  ASSERT_TRUE(reader.swingAdsorptionPhases[0].temperature.has_value());
  EXPECT_DOUBLE_EQ(*reader.swingAdsorptionPhases[0].temperature, 310.0);
  ASSERT_TRUE(reader.swingAdsorptionPhases[0].inletPressure.has_value());
  EXPECT_DOUBLE_EQ(*reader.swingAdsorptionPhases[0].inletPressure, 200000.0);
  EXPECT_DOUBLE_EQ(reader.swingAdsorptionPhases[0].duration, 0.03);
  EXPECT_FALSE(reader.swingAdsorptionPhases[1].temperature.has_value());
  ASSERT_TRUE(reader.swingAdsorptionPhases[1].inletPressure.has_value());
  EXPECT_DOUBLE_EQ(*reader.swingAdsorptionPhases[1].inletPressure, 50000.0);

  SwingAdsorption<Column> swing(reader);
  ASSERT_EQ(swing.subStages.size(), 2U);
  EXPECT_EQ(swing.subStages[0].name, "Adsorption");
  EXPECT_DOUBLE_EQ(swing.subStages[0].temperature, 310.0);
  EXPECT_DOUBLE_EQ(swing.subStages[0].pressure, 200000.0);
  EXPECT_DOUBLE_EQ(swing.subStages[0].duration, 0.03);
  EXPECT_EQ(swing.subStages[1].name, "Depressurization");
  EXPECT_DOUBLE_EQ(swing.subStages[1].temperature, reader.temperature);
  EXPECT_DOUBLE_EQ(swing.subStages[1].pressure, 50000.0);
  EXPECT_DOUBLE_EQ(swing.subStages[1].duration, 0.02);
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

TEST(GeometryInput, LoadsExplicitPackedBedGeometry)
{
  InputReader reader(examplePath("examples/CoBDP-alkanes-C6/breakthrough/simulation.json"));
  const Geometry& shape = reader.geometry;

  EXPECT_EQ(shape.kind, GeometryKind::PackedBed);
  EXPECT_DOUBLE_EQ(shape.voidFraction, reader.columnVoidFraction);
  EXPECT_DOUBLE_EQ(shape.solidToFluidVolumeRatio, 1.5);
  EXPECT_DOUBLE_EQ(shape.contactAreas.solidFluidPerSolidVolume, 6000.0);
  EXPECT_DOUBLE_EQ(shape.contactAreas.fluidSolidPerFluidVolume, 9000.0);
  EXPECT_DOUBLE_EQ(shape.pressureDrop.viscousCoefficient,
                   150.0 * shape.solidToFluidVolumeRatio * shape.solidToFluidVolumeRatio /
                       (reader.particleDiameter * reader.particleDiameter));
  EXPECT_DOUBLE_EQ(shape.pressureDrop.inertialCoefficient,
                   1.75 * shape.solidToFluidVolumeRatio / reader.particleDiameter);

  Column column(reader);
  EXPECT_EQ(column.geometry.kind, GeometryKind::PackedBed);
}

TEST(PressureDropLaw, SolvesSignedDarcyForchheimerVelocity)
{
  const PressureDropLaw law{2.0, 3.0};
  const double viscosity = 0.5;
  const double density = 2.0;
  const double pressureGradient = -7.0;

  const double velocity = law.velocity(viscosity, density, pressureGradient);
  EXPECT_GT(velocity, 0.0);
  EXPECT_NEAR(law.gradient(viscosity, density, velocity), -pressureGradient, 1.0e-12);

  const double reverseVelocity = law.velocity(viscosity, density, -pressureGradient);
  EXPECT_NEAR(reverseVelocity, -velocity, 1.0e-12);
  EXPECT_NEAR(law.gradient(viscosity, density, reverseVelocity), pressureGradient, 1.0e-12);
}

TEST(PressureDropLaw, MonolithDarcyVelocityUsesHydraulicPermeability)
{
  MonolithSpec specification;
  specification.channelShape = ChannelShape::Square;
  specification.internalChannelDimension = 1.0e-3;
  specification.outerDiameter = 1.0e-2;
  specification.numberOfChannels = 10;
  const Geometry geometry = makeGeometry(specification);

  const double viscosity = 1.8e-5;
  const double pressureGradient = -576.0;
  EXPECT_NEAR(geometry.pressureDrop.velocity(viscosity, 0.0, pressureGradient), 1.0, 1.0e-12);
}

TEST(GeometryInput, RequiresGeometryForBreakthrough)
{
  std::ifstream input(examplePath("examples/CoBDP-alkanes-C6/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  data.erase("Geometry");
  input.close();

  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_missing_geometry_test.json";
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
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_monolith_geometry_test.json";
  writeGeneralChemisorptionInput(path, false);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  data["Geometry"] = {
      {"Type", "Monolith"},      {"ChannelShape", "square"}, {"InternalChannelDimension", 1.0e-3},
      {"OuterDiameter", 1.0e-2}, {"NumberOfChannels", 10},   {"WashcoatThickness", 1.0e-4},
      {"ForchheimerCoefficient", 12.5},
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
  EXPECT_DOUBLE_EQ(shape.pressureDrop.viscousCoefficient, 32.0 / 1.0e-6);
  EXPECT_DOUBLE_EQ(shape.pressureDrop.inertialCoefficient, 12.5);
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

  computeDerivatives(column);

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

  EXPECT_NEAR(column.concentrationDot[25 * n + 0], 4.4401129708813942e-04, 1e-14);
  EXPECT_NEAR(column.concentrationDot[25 * n + 2], 1.9394951851919610e-03, 1e-14);
  EXPECT_NEAR(column.concentrationDot[50 * n + 5], 2.0346075380352105e-03, 1e-14);
}

TEST(ExampleMixturePrediction, RecomputesBeaBreakthroughEquilibriumLoadings)
{
  Column column = loadBeaBreakthroughColumn();
  const std::vector<double> savedEquilibriumAdsorption = column.equilibriumPhysisorption;

  computePhysisorptionEquilibriumLoadings(
      column.physisorptionMixture, column.numberOfGridPoints, column.numberOfComponents, column.maxIsothermTerms,
      column.iastPerformance, column.idealGasMolFractions, column.adsorbedMolFractions, column.numberOfMolecules,
      column.totalPressure, column.equilibriumPhysisorption, column.cachedPressure, column.cachedGrandPotential,
      column.moleFraction, column.gasTemperature, MixturePrediction::DrivingForceInput::MoleFraction,
      column.concentration, column.pH);

  EXPECT_LE(maxAbsDifference(column.equilibriumPhysisorption, savedEquilibriumAdsorption), 1e-10);

  const size_t n = column.numberOfComponents;
  EXPECT_DOUBLE_EQ(column.equilibriumPhysisorption[0], 0.0);
  EXPECT_NEAR(column.equilibriumPhysisorption[0 * n + 1], 0.070539163250446624, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[10 * n + 2], 0.044465087898088139, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[25 * n + 3], 0.027130378763624827, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[50 * n + 5], 0.01358645524192777, 1e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[50 * n + 7], 0.0080228059270143568, 1e-14);
}

TEST(ChemisorptionGeneral, DefaultsToBulkConcentrationWithoutPoreSurfaceTransport)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_general_chem_direct_test.json";
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

  computeDerivatives(column);
  const size_t adsorbateInlet = 1U;
  EXPECT_GT(column.chemisorptionDot[adsorbateInlet], 0.0);
  EXPECT_TRUE(column.poreConcentrationDot.empty());
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet], 0.0);
}

TEST(ChemisorptionGeneral, InitializesSurfacePoreStateAndFilmSinkWhenEnabled)
{
  const std::filesystem::path path =
      testTemporaryDirectory() / "ruptura_general_chem_transport_test.json";
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

  computeDerivatives(column);
  const size_t adsorbateInlet = 1U;
  EXPECT_GT(column.chemisorptionDot[adsorbateInlet], 0.0);
  // The pore reservoir holds the solid volume, so it loses rho_p per unit chemisorbed loading.
  const double poreLoadingConversion = column.particleDensity;
  EXPECT_NEAR(column.poreConcentrationDot[adsorbateInlet],
              -poreLoadingConversion * column.chemisorptionDot[adsorbateInlet], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[adsorbateInlet], 0.0, 1.0e-12);

  column.surfaceConcentration[adsorbateInlet] = 0.5 * column.concentration[adsorbateInlet];
  computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.maxChemisorptionSites, column.geometry, column.particleDensity, column.concentration,
                         column.physisorptionDot, column.chemisorptionDot, column.surfaceConcentration,
                         column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                         column.reactionChemisorptionSource);
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

TEST(ChemisorptionRate, AvramiMatchesDerivativeOfIntegratedModel)
{
  Chemisorption avrami;
  avrami.type = Chemisorption::Type::Avrami;
  avrami.rateCoefficient = 3.92e-2;
  avrami.order = 1.46;

  constexpr double equilibriumLoading = 2.06;
  constexpr double elapsedTime = 25.0;
  const double exponent = std::pow(avrami.rateCoefficient * elapsedTime, avrami.order);
  const double loading = equilibriumLoading * (1.0 - std::exp(-exponent));
  const double expected = avrami.order * std::pow(avrami.rateCoefficient, avrami.order) *
                          std::pow(elapsedTime, avrami.order - 1.0) * (equilibriumLoading - loading);

  EXPECT_NEAR(avrami.rate(equilibriumLoading, loading, 0.0, 298.15, elapsedTime), expected, 1.0e-14);
  EXPECT_DOUBLE_EQ(avrami.rate(equilibriumLoading, 0.0, 0.0, 298.15, 0.0), 0.0);
}

TEST(ReactionsInput, ParsesGeneralPowerLawReaction)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_parse.json";
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
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_names.json";
  writeReactionInput(path, "Physisorbed", "GeneralPowerLaw", true, R"json(["A"])json", R"json(["B"])json");

  InputReader reader(path.string());
  ASSERT_EQ(reader.reactions.size(), 1U);
  const Reaction& reaction = reader.reactions[0];
  EXPECT_EQ(reaction.reactants, std::vector<size_t>{1});
  EXPECT_EQ(reaction.products, std::vector<size_t>{2});
}

TEST(ReactionsAdsorbed, AddsAdsorbedSourcesWithoutBulkTransferAndLimitsQmax)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_adsorbed.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  const size_t reactant = 1U;
  const size_t product = 2U;
  column.physisorption[reactant] = 0.4;
  column.physisorption[product] = 0.0;

  computeDerivatives(column);
  EXPECT_LT(column.physisorptionDot[reactant], 0.0);
  EXPECT_GT(column.physisorptionDot[product], 0.0);
  EXPECT_NEAR(column.physisorptionDot[reactant], -column.physisorptionDot[product], 1.0e-12);
  EXPECT_NEAR(column.reactionPhysisorptionSource[reactant], column.physisorptionDot[reactant], 1.0e-12);
  EXPECT_NEAR(column.reactionPhysisorptionSource[product], column.physisorptionDot[product], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);

  EXPECT_NEAR(column.solidTemperatureDot[0], column.reactionHeat[0] / column.heatCapacitySolid, 1.0e-12);
  EXPECT_GT(column.solidTemperatureDot[0], 0.0);

  column.physisorption[product] = 1.0;
  computeDerivatives(column);
  EXPECT_NEAR(column.reactionPhysisorptionSource[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.reactionPhysisorptionSource[product], 0.0, 1.0e-12);
}

TEST(ReactionsPore, LangmuirHinshelwoodRequiresAndUsesPoreConcentration)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_lh_pore.json";
  writeReactionInput(path, "PoreConcentration", "LangmuirHinshelwood", false);

  InputReader reader(path.string());
  ASSERT_EQ(reader.reactions.size(), 1U);
  EXPECT_TRUE(reactionsRequirePoreConcentration(reader.reactions));

  Column column(reader);
  column.initialize();
  EXPECT_TRUE(column.surfacePoreTransportEnabled);
  EXPECT_EQ(column.poreConcentration.size(), column.concentration.size());

  computeDerivatives(column);
  const size_t reactant = 1U;
  const size_t product = 2U;
  EXPECT_LT(column.poreConcentrationDot[reactant], 0.0);
  EXPECT_GT(column.poreConcentrationDot[product], 0.0);
  EXPECT_LT(column.reactionPoreConcentrationSource[reactant], 0.0);
  EXPECT_GT(column.reactionPoreConcentrationSource[product], 0.0);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);
}

TEST(ReactionsAutoStop, RequiresSlowlyChangingPhaseStates)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_auto_stop.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  const size_t reactant = 1U;
  const size_t product = 2U;
  column.physisorption[reactant] = 0.4;
  computeDerivatives(column);
  EXPECT_FALSE(reactionAutoStopReached(column));

  column.physisorption[product] = 1.0;
  computeDerivatives(column);
  EXPECT_TRUE(reactionAutoStopReached(column));

  // 1e-4 max(1, |x|) per second is the bound: just above it the state is still changing.
  column.physisorptionDot[reactant] = 1.1e-4;
  EXPECT_FALSE(reactionAutoStopReached(column));

  // A reaction source alone does not stop the auto-stop: a catalytic column at steady state keeps reacting.
  computeDerivatives(column);
  column.reactionPhysisorptionSource[reactant] = 1.0;
  EXPECT_TRUE(reactionAutoStopReached(column));
}

TEST(ReactionsAutoStop, BreakthroughSchedulesAutoStopForStableReactionState)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_auto_stop_rk3.json";
  writeReactionInput(path);

  InputReader reader(path.string());
  reader.autoSimulationTime = true;
  reader.writeInterval = reader.timeStep;  // one RK3 step per output

  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_reaction_auto_stop_rk3");
  Breakthrough<Column> breakthrough(reader);
  breakthrough.column.physisorption[2U] = 1.0;

  EXPECT_FALSE(breakthrough.advanceToNextOutput());
  EXPECT_FALSE(breakthrough.autoSimulationTime);
  // 10% longer than the first output, rounded up to whole outputs: two outputs.
  EXPECT_DOUBLE_EQ(breakthrough.simulationTime, 2.0 * reader.writeInterval);
  EXPECT_TRUE(breakthrough.advanceToNextOutput());
  EXPECT_DOUBLE_EQ(breakthrough.time, 2.0 * reader.writeInterval);
}

TEST(ReactionsInput, RejectsLHOutsidePoreConcentration)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reaction_lh_bad_phase.json";
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

TEST(ChemisorptionInput, RejectsFractionalGeneralOrders)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_fractional_chem_order.json";
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

TEST(ChemisorptionInput, AcceptsFractionalAvramiOrder)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_fractional_avrami_order.json";
  writeGeneralChemisorptionInput(path, false);

  std::ifstream input(path);
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();

  auto& site = data["Components"][1]["ChemisorptionSites"][0];
  site["RateEquation"] = "Avrami";
  site["RateEquationParameters"] = {{"RateCoefficient", 3.92e-2}, {"Order", 1.46}};

  std::ofstream output(path);
  output << data;
  output.close();

  InputReader reader(path.string());
  ASSERT_EQ(reader.components[1].chemisorption.numberOfSites, 1U);
  EXPECT_EQ(reader.components[1].chemisorption.sites[0].type, Chemisorption::Type::Avrami);
  EXPECT_DOUBLE_EQ(reader.components[1].chemisorption.sites[0].order, 1.46);
}

TEST(ChemisorptionInput, RequiresExactlyTheParametersForTheSelectedModel)
{
  const std::filesystem::path validPath = testTemporaryDirectory() / "ruptura_strict_chem_valid.json";
  writeGeneralChemisorptionInput(validPath, false);

  std::ifstream input(validPath);
  const nlohmann::json valid = nlohmann::json::parse(input);
  input.close();

  auto expectRejected = [&](nlohmann::json data, const std::string& suffix, const std::string& expectedMessage)
  {
    const std::filesystem::path path =
        testTemporaryDirectory() / ("ruptura_strict_chem_" + suffix + ".json");
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
  missing["Components"][1]["ChemisorptionSites"][0]["RateEquationParameters"].erase("CapacityOrder");
  expectRejected(std::move(missing), "missing", "required key 'CapacityOrder'");

  nlohmann::json extra = valid;
  extra["Components"][1]["ChemisorptionSites"][0]["RateEquationParameters"]["Alpha"] = 1.0;
  expectRejected(std::move(extra), "extra", "unknown key 'Alpha'");

  nlohmann::json oldAlias = valid;
  auto& oldAliasParameters = oldAlias["Components"][1]["ChemisorptionSites"][0]["RateEquationParameters"];
  oldAliasParameters["Ea_ka"] = oldAliasParameters["AdsorptionActivationEnergy"];
  oldAliasParameters.erase("AdsorptionActivationEnergy");
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

TEST(ChemisorptionTransport, FilmPoreAndChemisorptionExchangeConservesMass)
{
  // Per unit bed volume the bulk gas (volume eps), the surface and pore reservoirs (each the solid volume 1 - eps)
  // and the solid (mass rho_p (1 - eps)) only exchange material, so their interphase rates must sum to zero.
  InputReader reader(examplePath("examples/Chemisorption-CO2/breakthrough/simulation.json"));
  Column column(reader);
  column.initialize();
  ASSERT_TRUE(column.surfacePoreTransportEnabled);

  const size_t carbonDioxide = 1;
  const size_t grid = 2;
  const size_t nodes = column.numberOfGridPoints + 1;
  const size_t components = column.numberOfComponents;
  const size_t bulkIndex = grid * components + carbonDioxide;
  column.concentration[bulkIndex] = column.concentration[carbonDioxide];  // the column starts with carrier only
  ASSERT_GT(column.concentration[bulkIndex], 0.0);
  for (size_t site = 0; site < column.maxChemisorptionSites; ++site)
  {
    const size_t siteIndex = site * nodes * components + bulkIndex;
    column.surfaceConcentration[siteIndex] = 0.6 * column.concentration[bulkIndex];
    column.poreConcentration[siteIndex] = 0.3 * column.concentration[bulkIndex];
    column.chemisorption[siteIndex] = 0.0;
  }

  computeDerivatives(column);

  const double voidFraction = column.voidFraction;
  const double solidMass = column.particleDensity * (1.0 - voidFraction);
  double reservoirRate = 0.0;
  double solidRate = column.physisorptionDot[bulkIndex];
  double scale = std::abs(voidFraction * column.bulkSpeciesSink[bulkIndex]);
  for (size_t site = 0; site < column.maxChemisorptionSites; ++site)
  {
    const size_t siteIndex = site * nodes * components + bulkIndex;
    reservoirRate += column.surfaceConcentrationDot[siteIndex] + column.poreConcentrationDot[siteIndex];
    solidRate += column.chemisorptionDot[siteIndex];
    scale += std::abs(column.surfaceConcentrationDot[siteIndex]) + std::abs(column.poreConcentrationDot[siteIndex]) +
             solidMass * std::abs(column.chemisorptionDot[siteIndex]);
  }
  ASSERT_GT(scale, 0.0);

  const double residual =
      -voidFraction * column.bulkSpeciesSink[bulkIndex] + (1.0 - voidFraction) * reservoirRate + solidMass * solidRate;
  EXPECT_LE(std::abs(residual), 1.0e-12 * scale);
}

TEST(EnergyBalance, GasWallExchangeConservesEnergyInPackedBed)
{
  // Per unit column length the gas (volume eps pi d_in^2 / 4) and the wall (volume pi (d_out^2 - d_in^2) / 4) only
  // exchange heat when the temperatures are uniform, the gas is at rest, the solid is at the gas temperature, and
  // the outer wall is insulated.
  const PackedBedTubeSpec specification{0.4, 3.0e-3, 0.025, 0.04};
  const Geometry geometry = makeGeometry(specification);
  const std::vector<Component> components{Component(0, "A", std::vector<Isotherm>{}, 1.0, 0.0, 0.0)};

  constexpr size_t numberOfGridPoints = 2;
  constexpr size_t nodes = numberOfGridPoints + 1;
  constexpr double gasDensityValue = 1.2;
  constexpr double heatCapacityGas = 1000.0;
  constexpr double wallDensity = 4000.0;
  constexpr double heatCapacityWall = 500.0;
  const std::vector<double> velocity(nodes, 0.0);
  const std::vector<double> gasDensity(nodes, gasDensityValue);
  std::vector<double> coeffDiffusion(nodes, 0.0);
  const std::vector<double> physisorptionDot(nodes, 0.0);
  const std::vector<double> gasTemperature(nodes, 320.0);
  const std::vector<double> wallTemperature(nodes, 300.0);
  std::vector<double> gasTemperatureDot(nodes, 0.0);
  std::vector<double> solidTemperatureDot(nodes, 0.0);
  std::vector<double> wallTemperatureDot(nodes, 0.0);

  computeEnergyDerivatives(components, numberOfGridPoints, 1, 0, 300.0, geometry, 1000.0, wallDensity, 0.0, 0.0, 10.0,
                           50.0, 0.0, heatCapacityGas, 900.0, heatCapacityWall, 0.05, velocity, gasDensity,
                           coeffDiffusion, physisorptionDot, {}, gasTemperature, gasTemperatureDot, gasTemperature,
                           solidTemperatureDot, wallTemperature, wallTemperatureDot);

  const double innerArea = std::numbers::pi * specification.internalDiameter * specification.internalDiameter / 4.0;
  const double outerArea = std::numbers::pi * specification.outerDiameter * specification.outerDiameter / 4.0;
  const double wallHeatPerLength = 50.0 * std::numbers::pi * specification.internalDiameter * (320.0 - 300.0);
  for (size_t grid = 1; grid < nodes; ++grid)
  {
    const double gasRate =
        specification.voidFraction * innerArea * gasDensityValue * heatCapacityGas * gasTemperatureDot[grid];
    const double wallRate = (outerArea - innerArea) * wallDensity * heatCapacityWall * wallTemperatureDot[grid];
    EXPECT_NEAR(wallRate, wallHeatPerLength, 1.0e-12 * wallHeatPerLength);
    EXPECT_NEAR(gasRate + wallRate, 0.0, 1.0e-12 * wallHeatPerLength);
  }
}

TEST(ChemisorptionEquilibrium, ComputesCompetitiveSiteLoadingsSeparatelyFromPhysisorption)
{
  std::ifstream input(examplePath("examples/Chemisorption-CO2-N2/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();

  auto& nitrogen = data["Components"][2];
  nitrogen["ChemisorptionSites"] = nlohmann::json::array({
      {
          {"Type", "Langmuir"},
          {"Parameters", {0.40, 1.0e-5}},
          {"RateEquation", "FirstOrder"},
          {"RateEquationParameters", {{"MassTransferCoefficient", 0.02}}},
          {"HeatOfAdsorption", 30000.0},
          {"MaximumLoading", 0.40},
      },
  });

  const std::filesystem::path path =
      testTemporaryDirectory() / "ruptura_competitive_chem_equilibrium.json";
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
  secondCarbonDioxideSite["MaximumLoading"] = 0.20;
  secondCarbonDioxideSite["Parameters"] = {0.20, 1.5e-5};
  data["Components"][1]["ChemisorptionSites"].push_back(secondCarbonDioxideSite);

  const std::filesystem::path missingSitePath =
      testTemporaryDirectory() / "ruptura_component_missing_chem_site.json";
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
      testTemporaryDirectory() / "ruptura_complete_chem_iast_bisection.json";
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
        testTemporaryDirectory() / ("ruptura_segregated_chem_" + method + "_" + iastMethod + ".json");
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
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_elovich_chem_test.json";
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
  computeDerivatives(column);

  const size_t adsorbateInlet = 1U;
  EXPECT_NEAR(column.chemisorptionDot[adsorbateInlet], 0.25 * column.concentration[adsorbateInlet], 1.0e-12);
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet], 0.0);
}

TEST(MultiSiteChemisorption, ParsesIndependentSitesAndSumsTheirBulkSink)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_multisite_chem_direct_test.json";
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

  computeDerivatives(column);

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
      testTemporaryDirectory() / "ruptura_multisite_chem_mixed_transport_test.json";
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

  computeDerivatives(column);

  EXPECT_GT(column.chemisorptionDot[adsorbateInlet], 0.0);
  EXPECT_GT(column.chemisorptionDot[secondSiteInlet], 0.0);
  EXPECT_DOUBLE_EQ(column.surfaceConcentrationDot[adsorbateInlet], 0.0);
  EXPECT_DOUBLE_EQ(column.poreConcentrationDot[adsorbateInlet], 0.0);

  const double loadingPrefactor = column.geometry.loadingPrefactor(column.particleDensity);
  // The pore reservoir holds the solid volume, so it loses rho_p per unit chemisorbed loading.
  EXPECT_NEAR(column.poreConcentrationDot[secondSiteInlet],
              -column.particleDensity * column.chemisorptionDot[secondSiteInlet], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[adsorbateInlet],
              loadingPrefactor * (column.physisorptionDot[adsorbateInlet] + column.chemisorptionDot[adsorbateInlet]),
              1.0e-12);

  column.surfaceConcentration[secondSiteInlet] = 0.5 * column.concentration[adsorbateInlet];
  computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.maxChemisorptionSites, column.geometry, column.particleDensity, column.concentration,
                         column.physisorptionDot, column.chemisorptionDot, column.surfaceConcentration,
                         column.bulkSpeciesSink, column.reactionPhysisorptionSource,
                         column.reactionChemisorptionSource);
  EXPECT_GT(column.bulkSpeciesSink[adsorbateInlet],
            loadingPrefactor * (column.physisorptionDot[adsorbateInlet] + column.chemisorptionDot[adsorbateInlet]));

  const std::filesystem::path restartPath =
      testTemporaryDirectory() / "ruptura_multisite_chem_restart_test.json";
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
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_multisite_chem_cvode_test.json";
  writeMultiSiteChemisorptionInput(path, true);

  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  CVODE integrator;
  Timing timings;
  integrator.initialize(column);
  EXPECT_NO_THROW(integrator.advance(column, 1.0e-5, timings));

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

TEST(ChemisorptionRate, FirstOrderUsesCoverageDependentDesorptionAboveEquilibrium)
{
  // Stadler et al., Front. Chem. Eng. 4 (2022) 1000064, Eqs. 6-10 and Table 3.
  Chemisorption site;
  site.type = Chemisorption::Type::FirstOrder;
  site.rateCoefficient = 0.1125;
  site.desorptionMassTransferCoefficient = 5.38e-4;
  site.desorptionCoverageEnergy = 74000.0;
  site.coverageReferenceLoading = 0.3421;
  const double temperature = 523.15;

  EXPECT_DOUBLE_EQ(site.rate(0.20, 0.10, 0.0, temperature), 0.1125 * (0.20 - 0.10));

  const double desorptionCoefficient = 5.38e-4 * std::exp(74000.0 * (0.15 / 0.3421) / (R * temperature));
  const double expected = desorptionCoefficient * (0.05 - 0.15);
  EXPECT_NEAR(site.rate(0.05, 0.15, 0.0, temperature), expected, 1.0e-12 * std::abs(expected));

  site.desorptionMassTransferCoefficient = 0.0;
  EXPECT_DOUBLE_EQ(site.rate(0.05, 0.15, 0.0, temperature), 0.1125 * (0.05 - 0.15));
}

TEST(ChemisorptionInput, ParsesFirstOrderDesorptionBranchAndRequiresReferenceLoading)
{
  std::ifstream input(examplePath("examples/Chemisorption-CO2-N2/breakthrough/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  input.close();

  nlohmann::json site = {
      {"Type", "Freundlich"},
      {"Parameters", {0.0140596, 4.54545}},
      {"MaximumLoading", 0.3421},
      {"RateEquation", "FirstOrder"},
      {"RateEquationParameters",
       {{"MassTransferCoefficient", 0.1125},
        {"DesorptionMassTransferCoefficient", 5.38e-4},
        {"DesorptionCoverageEnergy", 74000.0},
        {"CoverageReferenceLoading", 0.3421}}},
  };
  data["Components"][1]["ChemisorptionSites"] = nlohmann::json::array({site});

  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_first_order_desorption.json";
  std::ofstream(path) << data;
  InputReader reader(path.string());
  const Chemisorption& kinetics = reader.components[1].chemisorption.sites[0];
  EXPECT_EQ(kinetics.type, Chemisorption::Type::FirstOrder);
  EXPECT_DOUBLE_EQ(kinetics.rateCoefficient, 0.1125);
  EXPECT_DOUBLE_EQ(kinetics.desorptionMassTransferCoefficient, 5.38e-4);
  EXPECT_DOUBLE_EQ(kinetics.desorptionCoverageEnergy, 74000.0);
  EXPECT_DOUBLE_EQ(kinetics.coverageReferenceLoading, 0.3421);

  site["RateEquationParameters"].erase("CoverageReferenceLoading");
  data["Components"][1]["ChemisorptionSites"] = nlohmann::json::array({site});
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);
}

static nlohmann::json makeBulkWaterGasShiftInput()
{
  const auto inert = [](const std::string& name, double fraction, double molarMass)
  {
    return nlohmann::json{{"Name", name},
                          {"GasPhaseMolFraction", fraction},
                          {"MolecularWeight", molarMass},
                          {"AxialDispersionCoefficient", 0.0},
                          {"PhysisorptionSites", {{{"Type", "Henry"}, {"Parameters", {0.0}}}}}};
  };
  return {
      {"SimulationType", "Breakthrough"},
      {"BreakthroughIntegrator", "RungeKutta3"},
      {"BoundaryCondition", "FixedVelocity"},
      {"Temperature", 523.15},
      {"ColumnVoidFraction", 0.4},
      {"DynamicViscosity", 2.65e-5},
      {"ParticleDiameter", 2.0e-4},
      {"ParticleDensity", 800.0},
      {"InletPressure", 8.0e5},
      {"PressureGradient", 0.0},
      {"ColumnEntranceVelocity", 0.1},
      {"ColumnLength", 0.025},
      {"NumberOfGridPoints", 10},
      {"SimulationTime", 1.0e-4},
      {"TimeStep", 1.0e-4},
      {"Geometry", {{"Type", "PackedBed"}, {"ColumnVoidFraction", 0.4}, {"ParticleDiameter", 2.0e-4}}},
      {"Components",
       {{{"Name", "N2"}, {"GasPhaseMolFraction", 0.8}, {"CarrierGas", true}, {"MolecularWeight", 0.028}},
        inert("CO", 0.1, 0.028), inert("H2O", 0.1, 0.018), inert("CO2", 0.0, 0.044), inert("H2", 0.0, 0.002)}},
      {"Reactions",
       {{{"Phase", "BulkConcentration"},
         {"Style", "GeneralPowerLaw"},
         {"Reactants", {"CO", "H2O"}},
         {"Products", {"CO2", "H2"}},
         {"ForwardOrders", {1, 1}},
         {"BackwardOrders", {1, 1}},
         {"Kinetics",
          {{"forwardRateCoefficient", 2.0e-3},
           {"forwardActivationEnergy", 0.0},
           {"equilibriumConstant", 80.0},
           {"gibbsFreeEnergy", 0.0}}}}}},
  };
}

TEST(BulkReaction, WaterGasShiftActsOnBulkConcentrationsThroughTheSpeciesSink)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_bulk_wgs.json";
  std::ofstream(path) << makeBulkWaterGasShiftInput();
  InputReader reader(path.string());
  ASSERT_EQ(reader.reactions.size(), 1U);
  EXPECT_EQ(reader.reactions[0].phase, Reaction::Phase::BulkConcentration);

  Column column(reader);
  column.initialize();

  // Put feed gas plus some product at an interior node so both directions of the reaction contribute.
  const size_t components = column.numberOfComponents;
  const size_t grid = 3;
  const double total = 8.0e5 / (R * 523.15);
  const std::array<double, 5> fractions{0.78, 0.08, 0.09, 0.03, 0.02};  // N2, CO, H2O, CO2, H2
  for (size_t comp = 0; comp < components; ++comp)
  {
    column.concentration[grid * components + comp] = fractions[comp] * total;
  }

  computeDerivatives(column);

  const auto c = [&](size_t comp) { return column.concentration[grid * components + comp]; };
  const double rate = 2.0e-3 * (c(1) * c(2) - c(3) * c(4) / 80.0);
  ASSERT_GT(rate, 0.0);
  const std::array<double, 5> stoichiometry{0.0, -1.0, -1.0, 1.0, 1.0};
  for (size_t comp = 0; comp < components; ++comp)
  {
    const size_t index = grid * components + comp;
    EXPECT_NEAR(column.reactionBulkSource[index], stoichiometry[comp] * rate, 1.0e-12 * rate);
    // No sorption: the bulk sink is the negative reaction source, and loadings are untouched.
    EXPECT_NEAR(column.bulkSpeciesSink[index], -stoichiometry[comp] * rate, 1.0e-12 * rate);
    EXPECT_DOUBLE_EQ(column.physisorptionDot[index], 0.0);
  }
}

TEST(BulkReaction, RejectsEnergyBalanceAndSirk3)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_bulk_wgs_rejected.json";

  nlohmann::json withEnergy = makeBulkWaterGasShiftInput();
  withEnergy["energyBalance"] = true;
  std::ofstream(path) << withEnergy;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  nlohmann::json withSirk3 = makeBulkWaterGasShiftInput();
  withSirk3["BreakthroughIntegrator"] = "SIRK3";
  std::ofstream(path) << withSirk3;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  nlohmann::json restrictedOnSingleBed = makeBulkWaterGasShiftInput();
  restrictedOnSingleBed["Reactions"][0]["Adsorbents"] = {"catalyst"};
  std::ofstream(path) << restrictedOnSingleBed;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);
}

TEST(SwingAdsorptionInput, RepeatsCyclesWithPhaseFeedCompositionAndVelocity)
{
  nlohmann::json data = makeBulkWaterGasShiftInput();
  data.erase("Reactions");
  data.erase("SimulationTime");
  data["SimulationType"] = "SwingAdsorption";
  data["NumberOfCycles"] = 3;
  data["SwingAdsorptionPhases"] = {
      {{"Name", "Adsorption"}, {"InletPressure", 8.0e5}, {"Duration", 4.0e-4}},
      {{"Name", "Regeneration"},
       {"InletPressure", 1.0e5},
       {"ColumnEntranceVelocity", 0.4},
       {"GasPhaseMolFractions", {{"N2", 3.0}, {"H2O", 2.0}}},
       {"Duration", 6.0e-4}},
  };

  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_swing_cycles.json";
  std::ofstream(path) << data;
  InputReader reader(path.string());

  ASSERT_EQ(reader.swingAdsorptionPhases.size(), 6U);
  EXPECT_NEAR(reader.simulationTime, 3.0e-3, 1.0e-15);
  EXPECT_EQ(reader.swingAdsorptionPhases[0].name, "Adsorption (cycle 1)");
  EXPECT_EQ(reader.swingAdsorptionPhases[5].name, "Regeneration (cycle 3)");
  EXPECT_TRUE(reader.swingAdsorptionPhases[4].gasPhaseMolFractions.empty());
  EXPECT_FALSE(reader.swingAdsorptionPhases[4].columnEntranceVelocity.has_value());

  const InputReader::SwingAdsorptionPhase& regeneration = reader.swingAdsorptionPhases[3];
  ASSERT_EQ(regeneration.gasPhaseMolFractions.size(), reader.components.size());
  EXPECT_DOUBLE_EQ(regeneration.gasPhaseMolFractions[0], 0.6);  // N2, normalized from 3 : 2
  EXPECT_DOUBLE_EQ(regeneration.gasPhaseMolFractions[1], 0.0);  // CO not listed
  EXPECT_DOUBLE_EQ(regeneration.gasPhaseMolFractions[2], 0.4);  // H2O
  EXPECT_DOUBLE_EQ(*regeneration.columnEntranceVelocity, 0.4);

  SwingAdsorption<Column> swing(reader);
  ASSERT_EQ(swing.subStages.size(), 6U);
  EXPECT_DOUBLE_EQ(swing.subStages[0].velocity, 0.1);
  EXPECT_DOUBLE_EQ(swing.subStages[1].velocity, 0.4);
  EXPECT_EQ(swing.subStages[1].moleFractions, regeneration.gasPhaseMolFractions);
}

/**
 * \brief Writes `data` to a temporary input file and returns the message the InputReader rejects it with.
 */
static std::string inputRejectionMessage(const nlohmann::json& data, const std::string& fileName)
{
  const std::filesystem::path path = testTemporaryDirectory() / fileName;
  std::ofstream(path) << data;
  try
  {
    InputReader reader(path.string());
  }
  catch (const std::runtime_error& error)
  {
    return error.what();
  }
  return {};
}

static nlohmann::json makeTimeTestInput()
{
  nlohmann::json data = makeBulkWaterGasShiftInput();
  data.erase("Reactions");
  return data;
}

TEST(TimeInput, RejectsStepBasedKeysWithTheConvertedTimes)
{
  // The same cases as tests/test_convert.py, which checks that convert.py writes these times.
  nlohmann::json rk3 = makeTimeTestInput();
  rk3.erase("SimulationTime");
  rk3["TimeStep"] = 0.002;
  rk3["NumberOfTimeSteps"] = 750000;
  rk3["WriteEvery"] = 1000;
  rk3["PrintEvery"] = 50000;
  const std::string rk3Message = inputRejectionMessage(rk3, "ruptura_step_keys_rk3.json");
  EXPECT_NE(rk3Message.find(R"(with:    "TimeStep": 0.002, "SimulationTime": 1500, "WriteInterval": 2)"),
            std::string::npos)
      << rk3Message;
  EXPECT_NE(rk3Message.find("PrintEvery has no replacement"), std::string::npos) << rk3Message;
  EXPECT_NE(rk3Message.find("python convert.py"), std::string::npos) << rk3Message;

  nlohmann::json cvode = rk3;
  cvode["BreakthroughIntegrator"] = "CVODE";
  const std::string cvodeMessage = inputRejectionMessage(cvode, "ruptura_step_keys_cvode.json");
  EXPECT_NE(cvodeMessage.find(R"(with:    "SimulationTime": 1500, "WriteInterval": 2)"), std::string::npos)
      << cvodeMessage;
  EXPECT_NE(cvodeMessage.find("TimeStep is not used with CVODE"), std::string::npos) << cvodeMessage;

  nlohmann::json swing = makeTimeTestInput();
  swing.erase("SimulationTime");
  swing.erase("TimeStep");
  swing["SimulationType"] = "SwingAdsorption";
  swing["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"NumberOfTimeSteps", 1000000}},
                                    {{"Name", "Regeneration"}, {"InletPressure", 1.0e5}, {"NumberOfTimeSteps", 2000000}}};
  const std::string swingMessage = inputRejectionMessage(swing, "ruptura_step_keys_swing.json");
  EXPECT_NE(swingMessage.find(R"(Add:     "WriteInterval": 5)"), std::string::npos) << swingMessage;
  EXPECT_NE(swingMessage.find(R"(SwingAdsorptionPhases[0]: replace "NumberOfTimeSteps": 1000000 with "Duration": 500)"),
            std::string::npos)
      << swingMessage;
  EXPECT_NE(swingMessage.find(R"(SwingAdsorptionPhases[1]: replace "NumberOfTimeSteps": 2000000 with "Duration": 1000)"),
            std::string::npos)
      << swingMessage;
  EXPECT_NE(swingMessage.find("former default TimeStep of 5e-4 s"), std::string::npos) << swingMessage;
}

TEST(TimeInput, CvodeRejectsTimeStepAndRampSteps)
{
  nlohmann::json data = makeTimeTestInput();
  data["BreakthroughIntegrator"] = "CVODE";
  EXPECT_NE(inputRejectionMessage(data, "ruptura_cvode_time_step.json").find("TimeStep is not used by CVODE"),
            std::string::npos);

  data.erase("TimeStep");
  data["NumberOfInitTimeSteps"] = 10;
  EXPECT_NE(inputRejectionMessage(data, "ruptura_cvode_ramp.json").find("NumberOfInitTimeSteps is not used by CVODE"),
            std::string::npos);
}

TEST(TimeInput, SwingAdsorptionTakesItsTimeFromThePhases)
{
  nlohmann::json data = makeTimeTestInput();
  data["SimulationType"] = "SwingAdsorption";
  data["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"Duration", 3.0e-4}},
                                   {{"Name", "Regeneration"}, {"InletPressure", 1.0e5}, {"Duration", 2.0e-4}}};
  EXPECT_NE(inputRejectionMessage(data, "ruptura_swing_simulation_time.json").find("remove SimulationTime"),
            std::string::npos);

  data.erase("SimulationTime");
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_swing_durations.json";
  std::ofstream(path) << data;
  InputReader reader(path.string());
  EXPECT_FALSE(reader.autoSimulationTime);
  EXPECT_DOUBLE_EQ(reader.simulationTime, 5.0e-4);
}

TEST(TimeInput, WarnsWhenAnIntervalIsNotAWholeNumberOfSteps)
{
  nlohmann::json data = makeTimeTestInput();
  data["WriteInterval"] = 2.5e-4;
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_write_interval_warning.json";
  std::ofstream(path) << data;

  testing::internal::CaptureStderr();
  InputReader reader(path.string());
  const std::string warnings = testing::internal::GetCapturedStderr();
  EXPECT_NE(warnings.find("WriteInterval (0.00025 s) is not a multiple of TimeStep (0.0001 s)"), std::string::npos)
      << warnings;
}

TEST(TimeInput, RungeKuttaShortensOnlyTheLastStepOfAnInterval)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_rk3_intervals.json";
  std::ofstream(path) << makeTimeTestInput();
  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  RungeKutta3 integrator(3.0e-5);
  Timing timings;
  integrator.advance(column, 0.0, 1.0e-4, timings);  // three whole steps and a shortened one of 1e-5 s
  EXPECT_EQ(integrator.numberOfStepsTaken, 4U);
  integrator.advance(column, 1.0e-4, 1.9e-4, timings);  // three whole steps, up to rounding
  EXPECT_EQ(integrator.numberOfStepsTaken, 7U);
  EXPECT_TRUE(std::ranges::all_of(column.state, [](double value) { return std::isfinite(value); }));
}

TEST(TimeInput, BreakthroughOutputsAtWholeWriteIntervalsAndAtTheEnd)
{
  nlohmann::json data = makeTimeTestInput();
  data["WriteInterval"] = 2.0e-4;
  data["SimulationTime"] = 5.0e-4;
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_output_times.json";
  std::ofstream(path) << data;
  InputReader reader(path.string());

  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_output_times");
  Breakthrough<Column> breakthrough(reader);
  EXPECT_FALSE(breakthrough.advanceToNextOutput());
  EXPECT_DOUBLE_EQ(breakthrough.time, 2.0e-4);
  EXPECT_FALSE(breakthrough.advanceToNextOutput());
  EXPECT_DOUBLE_EQ(breakthrough.time, 4.0e-4);
  EXPECT_TRUE(breakthrough.advanceToNextOutput());
  EXPECT_DOUBLE_EQ(breakthrough.time, 5.0e-4);
  EXPECT_EQ(breakthrough.rk3.numberOfStepsTaken, 5U);
}

TEST(TimeInput, SwingWritesAnOutputAtAPhaseStartUnderThatPhasesConditions)
{
  nlohmann::json data = makeTimeTestInput();
  data.erase("SimulationTime");
  data["SimulationType"] = "SwingAdsorption";
  data["WriteInterval"] = 2.0e-4;
  data["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"Duration", 2.0e-4}},
                                   {{"Name", "Regeneration"}, {"Temperature", 600.0}, {"Duration", 2.0e-4}}};
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_swing_output_order.json";
  std::ofstream(path) << data;
  InputReader reader(path.string());

  const std::filesystem::path outputDirectory = testTemporaryDirectory() / "ruptura_swing_output_order";
  {
    WorkingDirectoryGuard workingDirectory(outputDirectory);
    SwingAdsorption<Column> swing(reader);
    swing.run();
  }

  // column.data: time [min] in column 2, position in column 3, gas temperature in column 6.
  std::ifstream columnData(outputDirectory / "column.data");
  std::vector<std::array<double, 3>> inletRows;
  for (std::string line; std::getline(columnData, line);)
  {
    if (line.empty() || line.front() == '#') continue;
    std::istringstream fields(line);
    std::array<double, 6> values{};
    for (double& value : values) fields >> value;
    if (values[2] == 0.0) inletRows.push_back({values[1] * 60.0, values[2], values[5]});
  }
  ASSERT_EQ(inletRows.size(), 3U);  // the start, the phase switch, and the end
  EXPECT_NEAR(inletRows[1][0], 2.0e-4, 1.0e-12);
  EXPECT_DOUBLE_EQ(inletRows[0][2], 523.15);
  EXPECT_DOUBLE_EQ(inletRows[1][2], 600.0);  // written after the second phase's conditions are applied
  EXPECT_NEAR(inletRows[2][0], 4.0e-4, 1.0e-12);
}

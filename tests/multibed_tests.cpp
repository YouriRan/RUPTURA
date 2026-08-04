#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "column_multibed.h"
#include "inputreader.h"
#include "integrators/compute_multibed.h"
#include "integrators/rk3.h"
#include "json.h"
#include "swing_adsorption.h"
#include "timing.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

bool allFinite(std::span<const double> values)
{
  return std::ranges::all_of(values, [](double value) { return std::isfinite(value); });
}

Reaction makeReaction(Reaction::Phase phase)
{
  Reaction reaction;
  reaction.phase = phase;
  reaction.style = Reaction::Style::GeneralPowerLaw;
  reaction.reactants = {1U};
  reaction.products = {2U};
  reaction.reactantStoichiometry = {1.0};
  reaction.productStoichiometry = {1.0};
  reaction.forwardOrders = {1.0};
  reaction.backwardOrders = {1.0};
  reaction.forwardRateCoefficient = 1.0;
  reaction.equilibriumConstant = 1.0e9;
  reaction.gibbsFreeEnergy = -1000.0;
  return reaction;
}

void disablePhysisorptionKinetics(MultibedColumn& column)
{
  for (MixturePrediction& mixture : column.physisorptionMixtures)
  {
    for (Component& component : mixture.components)
    {
      component.massTransferCoefficient = 0.0;
    }
  }
}

class WorkingDirectoryGuard
{
 public:
  explicit WorkingDirectoryGuard(const std::filesystem::path& path) : originalPath(std::filesystem::current_path())
  {
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }

  ~WorkingDirectoryGuard()
  {
    std::error_code error;
    std::filesystem::current_path(originalPath, error);
  }

 private:
  std::filesystem::path originalPath;
};
}  // namespace

TEST(MultibedExamples, ParseAndAdvanceWithAdsorbentSpecificParameters)
{
  const std::vector<std::pair<std::string, size_t>> examples = {
      {"examples/Multibed/C6-alkanes-layered/simulation.json", 2},
      {"examples/Multibed/CO2-C3H8-gradient/simulation.json", 3},
      {"examples/Multibed/CO2-N2-polishing/simulation.json", 2},
  };

  for (const auto& [path, expectedAdsorbents] : examples)
  {
    SCOPED_TRACE(path);
    InputReader reader(examplePath(path));
    ASSERT_EQ(reader.adsorbentComponents.size(), expectedAdsorbents);
    ASSERT_EQ(reader.adsorbentLengths.size(), expectedAdsorbents);
    ASSERT_EQ(reader.columnDistances.size(), reader.numberOfGridPoints + 1);

    MultibedColumn column(reader);
    column.initialize();
    EXPECT_EQ(column.numberOfAdsorbents, expectedAdsorbents);
    EXPECT_DOUBLE_EQ(column.columnDistances.front(), 0.0);
    EXPECT_DOUBLE_EQ(column.columnDistances.back(), column.columnLength);

    for (size_t ads = 0; ads < expectedAdsorbents; ++ads)
    {
      ASSERT_EQ(column.physisorptionMixtures[ads].components.size(), column.numberOfComponents);
      for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
      {
        if (!column.components[comp].isCarrierGas)
        {
          EXPECT_GT(column.physisorptionMixtures[ads].components[comp].massTransferCoefficient, 0.0);
        }
      }
    }

    RungeKutta3 integrator(reader.timeStep, false, 1);
    Timing timings;
    EXPECT_TRUE(integrator.propagate(column, 0, timings));
    EXPECT_TRUE(allFinite(column.state));
    EXPECT_TRUE(allFinite(column.totalPressure));
    EXPECT_TRUE(allFinite(column.interstitialGasVelocity));
    EXPECT_GT(*std::ranges::max_element(column.physisorption), 0.0);
  }
}

TEST(MultibedExamples, LayerFractionsFollowConfiguredSectionsAndInterfaces)
{
  InputReader reader(examplePath("examples/Multibed/C6-alkanes-layered/simulation.json"));
  MultibedColumn column(reader);
  column.initialize();

  ASSERT_EQ(column.numberOfAdsorbents, 2U);
  EXPECT_DOUBLE_EQ(column.fractionOfAdsorbent[0], 1.0);
  EXPECT_DOUBLE_EQ(column.fractionOfAdsorbent[1], 0.0);

  const size_t outlet = column.numberOfGridPoints * column.numberOfAdsorbents;
  EXPECT_DOUBLE_EQ(column.fractionOfAdsorbent[outlet], 0.0);
  EXPECT_DOUBLE_EQ(column.fractionOfAdsorbent[outlet + 1], 1.0);

  bool foundBlendedInterface = false;
  for (size_t grid = 0; grid < column.numberOfGridPoints + 1; ++grid)
  {
    const double left = column.fractionOfAdsorbent[grid * column.numberOfAdsorbents];
    const double right = column.fractionOfAdsorbent[grid * column.numberOfAdsorbents + 1];
    EXPECT_NEAR(left + right, 1.0, 1e-12);
    foundBlendedInterface = foundBlendedInterface || (left > 0.0 && left < 1.0 && right > 0.0 && right < 1.0);
  }
  EXPECT_TRUE(foundBlendedInterface);
}

#if BUILD_SUNDIALS
TEST(MultibedCvode, BreakthroughParsesInitializesAndPropagates)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  reader.breakthroughIntegrator = static_cast<size_t>(BreakthroughIntegrationScheme::CVODE);
  reader.autoNumberOfTimeSteps = false;
  reader.numberOfTimeSteps = 1;

  const std::filesystem::path outputDirectory =
      std::filesystem::path{::testing::TempDir()} / "ruptura_multibed_cvode_breakthrough";
  WorkingDirectoryGuard workingDirectory(outputDirectory);

  Breakthrough<MultibedColumn> breakthrough(reader);
  EXPECT_EQ(breakthrough.integrationScheme, BreakthroughIntegrationScheme::CVODE);
  EXPECT_NO_THROW(breakthrough.computeStep(0));
  EXPECT_TRUE(allFinite(breakthrough.column.state));
  EXPECT_TRUE(allFinite(breakthrough.column.stateDot));
  EXPECT_TRUE(allFinite(breakthrough.column.totalPressure));
  EXPECT_GT(*std::ranges::max_element(breakthrough.column.physisorption), 0.0);
}

TEST(MultibedCvode, PropagatesPoreReactionState)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  reader.reactions = {makeReaction(Reaction::Phase::PoreConcentration)};

  MultibedColumn column(reader);
  column.initialize();
  disablePhysisorptionKinetics(column);
  std::fill(column.poreConcentration.begin(), column.poreConcentration.end(), 0.0);
  column.poreConcentration[1U] = 0.4;

  CVODE integrator(1.0e-5, false, 1);
  Timing timings;
  integrator.initialize(column);
  EXPECT_NO_THROW(integrator.propagate(column, 0, timings));
  EXPECT_TRUE(allFinite(column.state));
  EXPECT_GT(column.poreConcentration[2U], 0.0);
  EXPECT_LT(column.reactionPoreConcentrationSource[1U], 0.0);
  EXPECT_GT(column.reactionPoreConcentrationSource[2U], 0.0);
}
#endif

TEST(MultibedColumn, StateLayoutMatchesColumnConventionAndRebindsCopies)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  MultibedColumn column(reader);
  column.initialize();

  const MultibedColumnStateLayout layout = column.stateLayout();
  const size_t componentBlockSize = (column.numberOfGridPoints + 1) * column.numberOfComponents;
  const size_t temperatureBlockSize = column.numberOfGridPoints + 1;

  EXPECT_EQ(layout.componentBlockSize(), componentBlockSize);
  EXPECT_EQ(column.stateSize(), 3 * componentBlockSize + 3 * temperatureBlockSize);
  EXPECT_EQ(column.concentration.data(), column.state.data());
  EXPECT_EQ(column.physisorption.data(), column.state.data() + componentBlockSize);
  EXPECT_EQ(column.chemisorption.data(), column.state.data() + 2 * componentBlockSize);
  EXPECT_EQ(column.gasTemperature.data(), column.state.data() + 3 * componentBlockSize);
  EXPECT_EQ(column.solidTemperature.data(), column.gasTemperature.data() + temperatureBlockSize);
  EXPECT_EQ(column.wallTemperature.data(), column.solidTemperature.data() + temperatureBlockSize);

  MultibedColumn copy(column);
  EXPECT_NE(copy.state.data(), column.state.data());
  EXPECT_EQ(copy.concentration.data(), copy.state.data());
  EXPECT_EQ(copy.physisorption.data(), copy.state.data() + componentBlockSize);
  EXPECT_EQ(copy.chemisorption.data(), copy.state.data() + 2 * componentBlockSize);
  ASSERT_FALSE(copy.concentration.empty());
  copy.concentration.front() += 1.0;
  EXPECT_NE(copy.concentration.front(), column.concentration.front());

  MultibedColumn assigned(reader);
  assigned = column;
  EXPECT_NE(assigned.state.data(), column.state.data());
  EXPECT_EQ(assigned.concentration.data(), assigned.state.data());
  EXPECT_EQ(assigned.physisorption.data(), assigned.state.data() + componentBlockSize);
  EXPECT_EQ(assigned.chemisorption.data(), assigned.state.data() + 2 * componentBlockSize);
}

TEST(MultibedReactions, PhysisorbedSourcesDoNotTransferMassToTheBulk)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  reader.energyBalance = true;
  reader.reactions = {makeReaction(Reaction::Phase::Physisorbed)};

  MultibedColumn column(reader);
  column.initialize();
  disablePhysisorptionKinetics(column);

  const size_t reactant = 1U;
  const size_t product = 2U;
  std::fill(column.physisorption.begin(), column.physisorption.end(), 0.0);
  column.physisorption[reactant] = 0.4;

  RK3MultibedHelpers::computeSorptionDerivatives(column);

  EXPECT_LT(column.reactionPhysisorptionSource[reactant], 0.0);
  EXPECT_GT(column.reactionPhysisorptionSource[product], 0.0);
  EXPECT_NEAR(column.physisorptionDot[reactant], -column.physisorptionDot[product], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);

  RK3MultibedHelpers::computeEnergyDerivatives(column);
  EXPECT_NEAR(column.solidTemperatureDot[0], column.reactionHeat[0] / column.heatCapacitySolid, 1.0e-12);
  EXPECT_GT(column.solidTemperatureDot[0], 0.0);
}

TEST(MultibedReactions, ParsesReactionTogetherWithAdsorbents)
{
  std::ifstream source(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  ASSERT_TRUE(source);
  nlohmann::json input;
  source >> input;
  input["NumberOfTimeSteps"] = 1;
  input["Reactions"] = nlohmann::json::array({{{"Phase", "Physisorbed"},
                                               {"Style", "GeneralPowerLaw"},
                                               {"Reactants", nlohmann::json::array({"CO2"})},
                                               {"Products", nlohmann::json::array({"N2"})},
                                               {"Stoichiometry", {{"Reactants", {1.0}}, {"Products", {1.0}}}},
                                               {"Kinetics",
                                                {{"forwardRateCoefficient", 1.0},
                                                 {"forwardActivationEnergy", 0.0},
                                                 {"equilibriumConstant", 1.0e9},
                                                 {"gibbsFreeEnergy", -1000.0}}}}});

  const std::filesystem::path path =
      std::filesystem::path{::testing::TempDir()} / "ruptura_multibed_reaction_input.json";
  std::ofstream output(path, std::ios::trunc);
  ASSERT_TRUE(output);
  output << input;
  output.close();

  InputReader reader(path.string());
  ASSERT_EQ(reader.adsorbentComponents.size(), 2U);
  ASSERT_EQ(reader.reactions.size(), 1U);
  EXPECT_EQ(reader.reactions.front().reactants, std::vector<size_t>{1U});
  EXPECT_EQ(reader.reactions.front().products, std::vector<size_t>{2U});

  const std::filesystem::path runDirectory =
      std::filesystem::path{::testing::TempDir()} / "ruptura_multibed_reaction_breakthrough";
  WorkingDirectoryGuard workingDirectory(runDirectory);
  Breakthrough<MultibedColumn> breakthrough(reader);
  breakthrough.computeStep(0);
  EXPECT_TRUE(allFinite(breakthrough.column.state));
}

TEST(MultibedReactions, ChemisorbedSourcesUseTheSharedSiteState)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  Chemisorption site;
  site.type = Chemisorption::Type::General;
  site.maximumLoading = 1.0;
  for (std::vector<Component>& components : reader.adsorbentComponents)
  {
    components[1].chemisorption.add(site);
    components[2].chemisorption.add(site);
  }
  reader.reactions = {makeReaction(Reaction::Phase::Chemisorbed)};

  MultibedColumn column(reader);
  column.initialize();
  disablePhysisorptionKinetics(column);

  const size_t reactant = 1U;
  const size_t product = 2U;
  std::fill(column.chemisorption.begin(), column.chemisorption.end(), 0.0);
  column.chemisorption[reactant] = 0.4;

  RK3MultibedHelpers::computeSorptionDerivatives(column);

  EXPECT_LT(column.reactionChemisorptionSource[reactant], 0.0);
  EXPECT_GT(column.reactionChemisorptionSource[product], 0.0);
  EXPECT_NEAR(column.chemisorptionDot[reactant], -column.chemisorptionDot[product], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);
}

TEST(MultibedReactions, PoreSourcesAreIntegratedWithoutDirectBulkTransfer)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  reader.reactions = {makeReaction(Reaction::Phase::PoreConcentration)};

  MultibedColumn column(reader);
  column.initialize();
  disablePhysisorptionKinetics(column);
  ASSERT_TRUE(column.surfacePoreTransportEnabled);
  ASSERT_EQ(column.poreConcentration.size(), column.concentration.size());

  const size_t reactant = 1U;
  const size_t product = 2U;
  std::fill(column.poreConcentration.begin(), column.poreConcentration.end(), 0.0);
  column.poreConcentration[reactant] = 0.4;

  RK3MultibedHelpers::computeSorptionDerivatives(column);

  EXPECT_LT(column.reactionPoreConcentrationSource[reactant], 0.0);
  EXPECT_GT(column.reactionPoreConcentrationSource[product], 0.0);
  EXPECT_NEAR(column.poreConcentrationDot[reactant], -column.poreConcentrationDot[product], 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[reactant], 0.0, 1.0e-12);
  EXPECT_NEAR(column.bulkSpeciesSink[product], 0.0, 1.0e-12);

  RungeKutta3 integrator(reader.timeStep, false, 1);
  Timing timings;
  EXPECT_TRUE(integrator.propagate(column, 0, timings));
  EXPECT_TRUE(allFinite(column.state));
  EXPECT_GT(column.poreConcentration[product], 0.0);
}

TEST(MultibedCompute, BulkSpeciesSinkUsesBedWeightedSolidLoading)
{
  InputReader reader(examplePath("examples/Multibed/C6-alkanes-layered/simulation.json"));
  MultibedColumn column(reader);
  column.initialize();

  size_t blendedGrid = column.numberOfGridPoints + 1;
  for (size_t grid = 0; grid < column.numberOfGridPoints + 1; ++grid)
  {
    const double fraction = column.fractionOfAdsorbent[grid * column.numberOfAdsorbents];
    if (fraction > 0.0 && fraction < 1.0)
    {
      blendedGrid = grid;
      break;
    }
  }
  ASSERT_LT(blendedGrid, column.numberOfGridPoints + 1);

  std::fill(column.physisorptionDot.begin(), column.physisorptionDot.end(), 0.0);
  for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
  {
    column.physisorptionDot[blendedGrid * column.numberOfComponents + comp] = 0.25 * static_cast<double>(comp + 1);
  }

  RK3MultibedHelpers::computeBulkSpeciesSink(column);

  double solidLoadingDensity = 0.0;
  for (size_t ads = 0; ads < column.numberOfAdsorbents; ++ads)
  {
    const double fraction = column.fractionOfAdsorbent[blendedGrid * column.numberOfAdsorbents + ads];
    solidLoadingDensity += fraction * (1.0 - column.adsorbentVoidFractions[ads]) * column.particleDensities[ads];
  }
  const double loadingPrefactor = solidLoadingDensity / column.totalVoidFraction[blendedGrid];

  for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
  {
    const size_t index = blendedGrid * column.numberOfComponents + comp;
    EXPECT_NEAR(column.bulkSpeciesSink[index], loadingPrefactor * column.physisorptionDot[index], 1e-12);
  }
}

TEST(MultibedSwingAdsorption, AdvancesConcreteMultibedColumnAcrossStages)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  reader.simulationType = InputReader::SimulationType::SwingAdsorption;
  reader.numberOfInitTimeSteps = 0;
  reader.printEvery = 1;
  reader.writeEvery = 1;

  InputReader::SwingAdsorptionPhase adsorption;
  adsorption.name = "adsorption";
  adsorption.temperature = reader.temperature;
  adsorption.inletPressure = reader.inletPressure;
  adsorption.numberOfSteps = 2;

  InputReader::SwingAdsorptionPhase regeneration;
  regeneration.name = "regeneration";
  regeneration.temperature = reader.temperature + 10.0;
  regeneration.inletPressure = 0.8 * reader.inletPressure;
  regeneration.numberOfSteps = 2;
  reader.swingAdsorptionPhases = {adsorption, regeneration};

  const std::filesystem::path outputDirectory = std::filesystem::path{::testing::TempDir()} / "ruptura_multibed_swing";
  WorkingDirectoryGuard workingDirectory(outputDirectory);

  SwingAdsorption<MultibedColumn> swing(reader);
  EXPECT_EQ(swing.breakthrough.column.numberOfAdsorbents, 2U);
  ASSERT_EQ(swing.subStages.size(), 2U);

  swing.run();

  EXPECT_TRUE(allFinite(swing.breakthrough.column.state));
  EXPECT_TRUE(allFinite(swing.breakthrough.column.totalPressure));
  EXPECT_TRUE(allFinite(swing.breakthrough.column.interstitialGasVelocity));
  EXPECT_DOUBLE_EQ(swing.breakthrough.column.inletPressure, *regeneration.inletPressure);
}

#if BUILD_SUNDIALS
TEST(MultibedSwingAdsorption, AdvancesWithCvode)
{
  InputReader reader(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  reader.simulationType = InputReader::SimulationType::SwingAdsorption;
  reader.breakthroughIntegrator = static_cast<size_t>(BreakthroughIntegrationScheme::CVODE);
  reader.numberOfInitTimeSteps = 0;
  reader.printEvery = 1;
  reader.writeEvery = 1;

  InputReader::SwingAdsorptionPhase adsorption;
  adsorption.name = "adsorption";
  adsorption.temperature = reader.temperature;
  adsorption.inletPressure = reader.inletPressure;
  adsorption.numberOfSteps = 2;

  InputReader::SwingAdsorptionPhase regeneration;
  regeneration.name = "regeneration";
  regeneration.temperature = reader.temperature + 10.0;
  regeneration.inletPressure = 0.9 * reader.inletPressure;
  regeneration.numberOfSteps = 1;
  reader.swingAdsorptionPhases = {adsorption, regeneration};

  const std::filesystem::path outputDirectory =
      std::filesystem::path{::testing::TempDir()} / "ruptura_multibed_cvode_swing";
  WorkingDirectoryGuard workingDirectory(outputDirectory);

  SwingAdsorption<MultibedColumn> swing(reader);
  EXPECT_EQ(swing.breakthrough.integrationScheme, BreakthroughIntegrationScheme::CVODE);
  EXPECT_NO_THROW(swing.run());
  EXPECT_TRUE(allFinite(swing.breakthrough.column.state));
  EXPECT_TRUE(allFinite(swing.breakthrough.column.totalPressure));
  EXPECT_DOUBLE_EQ(swing.breakthrough.column.inletPressure, *regeneration.inletPressure);
}
#endif

TEST(MultibedOutput, UsesSingleBedComponentColumnSchema)
{
  InputReader reader(examplePath("examples/Multibed/C6-alkanes-layered/simulation.json"));
  MultibedColumn column(reader);
  column.initialize();

  const std::filesystem::path outputDirectory = ::testing::TempDir();
  std::vector<std::filesystem::path> componentPaths(column.numberOfComponents);
  std::vector<std::ofstream> componentStreams(column.numberOfComponents);
  for (size_t comp = 0; comp < column.numberOfComponents; ++comp)
  {
    componentPaths[comp] = outputDirectory / ("ruptura_multibed_component_" + std::to_string(comp) + ".data");
    componentStreams[comp].open(componentPaths[comp], std::ios::trunc);
    ASSERT_TRUE(componentStreams[comp]);
  }

  const std::filesystem::path columnPath = outputDirectory / "ruptura_multibed_column.data";
  std::ofstream columnStream(columnPath, std::ios::trunc);
  ASSERT_TRUE(columnStream);

  column.writeOutputHeader(componentStreams, columnStream);
  column.writeOutput(componentStreams, columnStream, 0.0);
  for (std::ofstream& stream : componentStreams) stream.close();
  columnStream.close();

  std::ifstream componentInput(componentPaths[1]);
  ASSERT_TRUE(componentInput);
  std::string line;
  std::vector<double> values;
  while (std::getline(componentInput, line))
  {
    if (line.empty() || line.front() == '#') continue;
    std::istringstream row(line);
    double value = 0.0;
    while (row >> value) values.push_back(value);
    break;
  }

  ASSERT_EQ(values.size(), 14U);
  EXPECT_DOUBLE_EQ(values[5], column.moleFraction[1]);
  EXPECT_DOUBLE_EQ(values[6], column.physisorption[1]);
  EXPECT_DOUBLE_EQ(values[7], column.physisorptionDot[1]);
  EXPECT_DOUBLE_EQ(values[8], 0.0);
  EXPECT_DOUBLE_EQ(values[9], 0.0);
  EXPECT_DOUBLE_EQ(values[10], column.partialPressure[1]);
  EXPECT_DOUBLE_EQ(values[11], column.equilibriumPhysisorption[1]);
  EXPECT_DOUBLE_EQ(values[13], 0.0);
}

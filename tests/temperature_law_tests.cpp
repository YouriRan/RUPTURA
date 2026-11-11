#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "inputreader.h"
#include "json.h"
#include "mixture_prediction.h"
#include "temperature_law.h"
#include "temporary_directory.h"
#include "utils.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

/**
 * \brief The MOR breakthrough example: helium carrier, CO2 and C3H8 with dual-site Langmuir isotherms, 300 K.
 */
nlohmann::json morInput()
{
  std::ifstream source(examplePath("examples/MOR-CO2-C3H8/breakthrough/simulation.json"));
  nlohmann::json input;
  source >> input;
  input["NumberOfGridPoints"] = 4;
  return input;
}

/**
 * \brief Adds what the energy balance needs to an input.
 */
void enableEnergyBalance(nlohmann::json& input)
{
  input["energyBalance"] = true;
  input["wallDensity"] = 7800.0;
  input["heatCapacityWall"] = 500.0;
}

InputReader readInput(const nlohmann::json& input)
{
  const std::filesystem::path path =
      testTemporaryDirectory() /
      (std::string{"ruptura_temperature_"} + ::testing::UnitTest::GetInstance()->current_test_info()->name() + ".json");
  std::ofstream output(path);
  output << input;
  output.close();
  return InputReader(path.string());
}

void expectRejected(const nlohmann::json& input, const std::string& expectedMessage)
{
  try
  {
    static_cast<void>(readInput(input));
    ADD_FAILURE() << "input accepted; expected an error containing: " << expectedMessage;
  }
  catch (const std::runtime_error& error)
  {
    EXPECT_NE(std::string(error.what()).find(expectedMessage), std::string::npos) << error.what();
  }
}

nlohmann::json vantHoff(double value) { return {{"Value", value}, {"Law", "VantHoff"}}; }
}  // namespace

TEST(TemperatureInput, AcceptsLawsAtThePreselectedPositions)
{
  nlohmann::json input = morInput();
  nlohmann::json& co2 = input["Components"][1];
  co2["HeatOfAdsorption"] = 25000.0;
  co2["ReferenceTemperature"] = 300.0;
  co2["PhysisorptionSites"] = {{{"Type", "Sips"},
                                {"Parameters",
                                 {{{"Value", 4.4}, {"Law", "Exponential"}, {"Chi", 0.4}},
                                  vantHoff(2.91e-4),
                                  {{"Value", 1.2}, {"Law", "Inverse"}, {"Alpha", 0.3}, {"ReferenceTemperature", 310.0}}}}}};

  const InputReader reader = readInput(input);
  const Isotherm& sips = reader.components[1].isotherm.sites[0];
  ASSERT_EQ(sips.temperatureLaws.size(), 3U);
  EXPECT_EQ(sips.temperatureLaw(0)->type, TemperatureLaw::Type::Exponential);
  EXPECT_DOUBLE_EQ(sips.temperatureLaw(0)->coefficient, 0.4);
  // The affinity law inherits the heat of adsorption and the component's reference temperature.
  EXPECT_EQ(sips.temperatureLaw(1)->type, TemperatureLaw::Type::VantHoff);
  EXPECT_TRUE(sips.temperatureLaw(1)->energyInherited);
  EXPECT_DOUBLE_EQ(sips.temperatureLaw(1)->coefficient, 25000.0);
  EXPECT_EQ(sips.temperatureLaw(1)->referenceTemperature, 300.0);
  EXPECT_EQ(sips.temperatureLaw(2)->referenceTemperature, 310.0);
  EXPECT_DOUBLE_EQ(sips.parameters[1], 2.91e-4);
  EXPECT_DOUBLE_EQ(reader.components[1].heatOfAdsorption, 25000.0);
}

TEST(TemperatureInput, RejectsLawsOutsideThePreselectedPositionsAndIncompleteLaws)
{
  const auto withParameters = [](nlohmann::json parameters, const std::string& type = "Langmuir")
  {
    nlohmann::json input = morInput();
    input["Components"][1]["HeatOfAdsorption"] = 25000.0;
    input["Components"][1]["PhysisorptionSites"] = {{{"Type", type}, {"Parameters", std::move(parameters)}}};
    return input;
  };

  expectRejected(withParameters({{{"Value", 4.4}, {"Law", "VantHoff"}}, 2.91e-4}),
                 "the saturation loading cannot use the temperature law 'VantHoff'; use Exponential or Linear");
  expectRejected(withParameters({4.4, {{"Value", 2.91e-4}, {"Law", "Polynomial"}}}), "temperature law 'Polynomial'");
  expectRejected(withParameters({4.4, 2.91e-4, {{"Value", 7.0}, {"Law", "Linear"}, {"Slope", 0.1}}}, "pH-Langmuir"),
                 "the pH_0 cannot depend on temperature");
  expectRejected(withParameters({4.4, 1.0e-5, {{"Value", 1.0}, {"Law", "Linear"}, {"Slope", 0.1}}}, "BET"),
                 "cannot depend on temperature");
  expectRejected(withParameters({{{"Value", 4.4}, {"Law", "Exponential"}, {"Chi", 0.4}}, 2.91e-4}),
                 "needs a ReferenceTemperature");
  expectRejected(withParameters({{{"Value", 4.4}, {"Law", "Linear"}, {"ReferenceTemperature", 300.0}}, 2.91e-4}),
                 "needs \"Slope\"");
  expectRejected(withParameters({4.4, {{"Value", 2.91e-4}, {"Law", "VantHoff"}, {"Chi", 1.0}}}),
                 "unknown key 'Chi' in the VantHoff law");
  // Redlich-Peterson has two affinity-like parameters, so its VantHoff laws need their own Energy.
  expectRejected(withParameters({1.0e-5, vantHoff(1.0e-4), 0.8}, "Redlich-Peterson"), "needs \"Energy\"");
}

TEST(TemperatureInput, HeatOfAdsorptionAndAffinityEnergyLinkBothWays)
{
  // No heat given: the energy of the affinity law is the energy-balance heat.
  nlohmann::json input = morInput();
  input["Components"][1]["PhysisorptionSites"][0]["Parameters"][1] =
      {{"Value", 2.91e-4}, {"Law", "VantHoff"}, {"Energy", 30000.0}};
  EXPECT_DOUBLE_EQ(readInput(input).components[1].heatOfAdsorption, 30000.0);

  // Explicit heat wins: a law with its own Energy keeps it, and the heat stays as given.
  input["Components"][1]["HeatOfAdsorption"] = 25000.0;
  {
    const InputReader reader = readInput(input);
    EXPECT_DOUBLE_EQ(reader.components[1].heatOfAdsorption, 25000.0);
    EXPECT_DOUBLE_EQ(reader.components[1].isotherm.sites[0].temperatureLaw(1)->coefficient, 30000.0);
  }

  // An affinity law without Energy needs a heat to take it from.
  nlohmann::json withoutHeat = morInput();
  withoutHeat["Components"][1]["PhysisorptionSites"][0]["Parameters"][1] = vantHoff(2.91e-4);
  expectRejected(withoutHeat, "there is no HeatOfAdsorption to take it from");
}

TEST(TemperatureInput, SitesWithDifferentEnergiesNeedAHeatOnlyForTheEnergyBalance)
{
  nlohmann::json input = morInput();
  input["Components"][1]["PhysisorptionSites"][0]["Parameters"][1] =
      {{"Value", 2.91e-4}, {"Law", "VantHoff"}, {"Energy", 30000.0}};
  input["Components"][1]["PhysisorptionSites"][1]["Parameters"][1] =
      {{"Value", 6.96e-7}, {"Law", "VantHoff"}, {"Energy", 18000.0}};

  // Isothermal: the heat is not used.
  EXPECT_DOUBLE_EQ(readInput(input).components[1].heatOfAdsorption, 0.0);

  enableEnergyBalance(input);
  expectRejected(input, "CO2: 30000 J/mol, 18000 J/mol");

  input["Components"][1]["HeatOfAdsorption"] = 26000.0;
  EXPECT_DOUBLE_EQ(readInput(input).components[1].heatOfAdsorption, 26000.0);
}

TEST(TemperatureInput, ChemisorptionSiteHeatFollowsTheMostSpecificSetting)
{
  nlohmann::json input = morInput();
  nlohmann::json& co2 = input["Components"][1];
  co2["HeatOfAdsorption"] = 25000.0;
  const nlohmann::json firstOrder = {{"MassTransferCoefficient", 0.1}};
  co2["ChemisorptionSites"] = {
      // Inherits the component heat.
      {{"Type", "Langmuir"}, {"Parameters", {1.0, vantHoff(1.0e-5)}}, {"RateEquationParameters", firstOrder},
       {"MaximumLoading", 1.0}},
      // The site's own law energy is more specific than the component heat.
      {{"Type", "Langmuir"},
       {"Parameters", {1.0, {{"Value", 1.0e-5}, {"Law", "VantHoff"}, {"Energy", 80000.0}}}},
       {"RateEquationParameters", firstOrder},
       {"MaximumLoading", 1.0}},
      // A site heat is inherited by its law.
      {{"Type", "Langmuir"}, {"Parameters", {1.0, vantHoff(1.0e-5)}}, {"RateEquationParameters", firstOrder},
       {"HeatOfAdsorption", 60000.0}, {"MaximumLoading", 1.0}}};

  const InputReader reader = readInput(input);
  const std::vector<Chemisorption>& sites = reader.components[1].chemisorption.sites;
  ASSERT_EQ(sites.size(), 3U);
  EXPECT_DOUBLE_EQ(sites[0].heatOfChemisorption, 25000.0);
  EXPECT_DOUBLE_EQ(sites[0].isotherm->temperatureLaw(1)->coefficient, 25000.0);
  EXPECT_DOUBLE_EQ(sites[1].heatOfChemisorption, 80000.0);
  EXPECT_DOUBLE_EQ(sites[2].heatOfChemisorption, 60000.0);
  EXPECT_DOUBLE_EQ(sites[2].isotherm->temperatureLaw(1)->coefficient, 60000.0);
}

TEST(TemperatureInput, KineticCoefficientsAcceptArrheniusLaws)
{
  nlohmann::json input = morInput();
  input["Components"][1]["ChemisorptionSites"] = {
      {{"Type", "Langmuir"},
       {"Parameters", {1.0, 1.0e-5}},
       {"RateEquation", "PseudoNth"},
       {"RateEquationParameters",
        {{"RateCoefficient", {{"Value", 0.05}, {"Law", "Arrhenius"}, {"ActivationEnergy", 20000.0},
                              {"ReferenceTemperature", 300.0}}},
         {"Order", 1.5}}},
       {"MaximumLoading", 3.0}}};
  const InputReader reader = readInput(input);
  const Chemisorption& site = reader.components[1].chemisorption.sites[0];
  EXPECT_EQ(site.rateCoefficient.type, TemperatureLaw::Type::Arrhenius);
  const double k = 0.05 * std::exp(-20000.0 / R * (1.0 / 330.0 - 1.0 / 300.0));
  EXPECT_NEAR(site.rate(2.0, 1.75, 0.0, 330.0), k * std::pow(0.25, 1.5), 1.0e-14);

  input["Components"][1]["ChemisorptionSites"][0]["RateEquationParameters"]["RateCoefficient"] = vantHoff(0.05);
  expectRejected(input, "the RateCoefficient cannot use the temperature law 'VantHoff'; use Arrhenius");
}

TEST(TemperatureInput, ReactionHeatAndEquilibriumConstantLinkBothWays)
{
  nlohmann::json input = morInput();
  input["Reactions"] = {{{"Phase", "Physisorbed"},
                         {"Style", "GeneralPowerLaw"},
                         {"Reactants", {"CO2"}},
                         {"Products", {"C3H8"}},
                         {"Kinetics",
                          {{"ForwardRateCoefficient",
                            {{"Value", 2.0}, {"Law", "Arrhenius"}, {"ActivationEnergy", 50000.0}}},
                           {"EquilibriumConstant", {{"Value", 3.0}, {"Law", "VantHoff"}, {"Energy", 40000.0}}}}}}};
  {
    const InputReader reader = readInput(input);
    EXPECT_DOUBLE_EQ(reader.reactions[0].heatOfReaction, 40000.0);
    EXPECT_NEAR(reader.reactions[0].forwardRateConstant(350.0), 2.0 * std::exp(-50000.0 / (R * 350.0)), 1.0e-20);
  }

  input["Reactions"][0]["Kinetics"]["HeatOfReaction"] = 35000.0;
  EXPECT_DOUBLE_EQ(readInput(input).reactions[0].heatOfReaction, 35000.0);

  input["Reactions"][0]["Kinetics"].erase("HeatOfReaction");
  input["Reactions"][0]["Kinetics"]["EquilibriumConstant"] = vantHoff(3.0);
  expectRejected(input, "there is no HeatOfReaction to take it from");
}

TEST(TemperatureInput, AdsorbentOverrideOfTheHeatMovesTheInheritedEnergy)
{
  nlohmann::json input = morInput();
  input["Components"][1]["HeatOfAdsorption"] = 25000.0;
  input["Components"][1]["PhysisorptionSites"][0]["Parameters"][1] = vantHoff(2.91e-4);
  input["Adsorbents"] = {{{"Name", "front"}, {"ComponentParameters", {{"CO2", {{"HeatOfAdsorption", 40000.0}}}}}},
                         {{"Name", "tail"}}};
  input["ColumnSections"] = {{{"Adsorbent", "front"}, {"Length", 0.14}, {"NumberOfGridPoints", 2}},
                             {{"InterfaceLength", 0.02}, {"NumberOfGridPoints", 2}},
                             {{"Adsorbent", "tail"}, {"Length", 0.14}, {"NumberOfGridPoints", 2}}};

  const InputReader reader = readInput(input);
  ASSERT_EQ(reader.adsorbentComponents.size(), 2U);
  const Component& front = reader.adsorbentComponents[0][1];
  const Component& tail = reader.adsorbentComponents[1][1];
  EXPECT_DOUBLE_EQ(front.heatOfAdsorption, 40000.0);
  EXPECT_DOUBLE_EQ(front.isotherm.sites[0].temperatureLaw(1)->coefficient, 40000.0);
  EXPECT_DOUBLE_EQ(tail.heatOfAdsorption, 25000.0);
  EXPECT_DOUBLE_EQ(tail.isotherm.sites[0].temperatureLaw(1)->coefficient, 25000.0);
}

TEST(TemperatureInput, FittingRejectsTemperatureLaws)
{
  nlohmann::json input = morInput();
  input["SimulationType"] = "Fitting";
  input["Components"][1]["HeatOfAdsorption"] = 25000.0;
  input["Components"][1]["PhysisorptionSites"][0]["Parameters"][1] = vantHoff(2.91e-4);
  expectRejected(input, "Fitting fits isotherm parameters at one temperature");
}

TEST(TemperatureInput, RejectsRemovedKeysWithTheirReplacement)
{
  nlohmann::json input = morInput();
  input["Components"][1]["nonIsothermal"] = true;
  input["Components"][2]["PhysisorptionSites"] = {{{"Type", "GAB"}, {"Parameters", {35.0, 0.15, 0.51, 6600.0, 0.0}}}};
  input["Components"][1]["ChemisorptionSites"] = {
      {{"Type", "Langmuir"},
       {"Parameters", {1.0, 1.0e-5}},
       {"RateEquation", "General"},
       {"RateEquationParameters", {{"AdsorptionRateCoefficient", 1.0e-4}, {"AdsorptionActivationEnergy", 30000.0}}},
       {"MaximumLoading", 1.0}}};
  input["Reactions"] = {{{"Kinetics",
                          {{"forwardRateCoefficient", 1.0},
                           {"forwardActivationEnergy", 0.0},
                           {"equilibriumConstant", 2.0},
                           {"gibbsFreeEnergy", -1000.0}}}}};

  try
  {
    static_cast<void>(readInput(input));
    FAIL() << "the removed keys were accepted";
  }
  catch (const std::runtime_error& error)
  {
    const std::string message = error.what();
    for (const std::string expected :
         {"Components[1]: replace \"nonIsothermal\": true by writing the affinity of each site as {\"Value\": b, "
          "\"Law\": \"VantHoff\"}",
          "Components[2].PhysisorptionSites[0]: replace the GAB Parameters [35.0,0.15,0.51,6600.0,0.0] with [35.0, "
          "{\"Value\": 0.15, \"Law\": \"VantHoff\", \"Energy\": 6600.0}, {\"Value\": 0.51, \"Law\": \"VantHoff\", "
          "\"Energy\": 0.0}]",
          "replace \"AdsorptionActivationEnergy\": 30000.0 by \"AdsorptionRateCoefficient\": {\"Value\": k, \"Law\": "
          "\"Arrhenius\", \"ActivationEnergy\": 30000.0}",
          "Reactions[0].Kinetics: remove \"forwardActivationEnergy\": 0.0",
          "Reactions[0].Kinetics: replace \"gibbsFreeEnergy\": -1000.0 by \"HeatOfReaction\": 1000",
          "Run 'python convert.py"})
    {
      EXPECT_NE(message.find(expected), std::string::npos) << "missing: " << expected << "\nin: " << message;
    }
  }
}

TEST(TemperatureInput, ChemisorbedGABUsesTheTemperatureOfTheColumn)
{
  // The chemisorption mixture once took the first site's heat as its scaling energy, so GAB saw T/50700.
  nlohmann::json input = morInput();
  input["Components"][1]["ChemisorptionSites"] = {
      {{"Type", "GAB"},
       {"Parameters", {2.58, {{"Value", 0.155}, {"Law", "VantHoff"}, {"Energy", 6600.0}}, 0.871}},
       {"RateEquationParameters", {{"MassTransferCoefficient", 0.25}}},
       {"HeatOfAdsorption", 50700.0},
       {"MaximumLoading", 2.58}}};
  input["MixturePredictionMethod"] = "SPI";
  const InputReader reader = readInput(input);

  MixturePrediction physisorption(reader);
  MixturePrediction chemisorption = MixturePrediction::makeChemisorptionPrediction(physisorption, reader.components);
  std::vector<double> loadings(reader.components.size(), 0.0);
  constexpr double pressure = 1500.0;
  chemisorption.predictPureComponentLoadings(pressure, loadings, 298.15);

  const double saturationPressure = std::pow(10.0, 8.07131 - 1730.63 / (298.15 - 39.724)) * 133.32236842105263;
  const double x = pressure / saturationPressure;
  const double c = 0.155 * std::exp(6600.0 / (R * 298.15));
  const double k = 0.871;
  EXPECT_NEAR(loadings[1], 2.58 * c * k * x / ((1.0 - k * x) * (1.0 + k * x * (c - 1.0))), 1.0e-12);
  EXPECT_GT(loadings[1], 0.1);
}

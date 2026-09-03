#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "isotherm.h"
#include "json.h"
#include "macrostate_particle_distribution.h"
#include "mixture_prediction.h"
#include "multi_site_isotherm.h"

namespace
{
std::filesystem::path mpdTemporaryPath(const std::string& name)
{
  return std::filesystem::temp_directory_path() / ("ruptura_mpd_" + name);
}

void writeMPDFile(const std::filesystem::path& path, const std::string& contents)
{
  std::ofstream output(path);
  output << contents;
}
}  // namespace

TEST(Isotherm, GABImplementsElfvingEquationsFiveThroughSeven)
{
  const Isotherm gab(Isotherm::Type::GAB, {2.58, 0.155, 0.871, 6600.0, 0.0}, true);
  constexpr double temperature = 298.15;
  constexpr double pressure = 0.019 * 105000.0;
  const double scale = std::exp(1.0 / (R * temperature));
  const double temperatureFromScale = 1.0 / (R * std::log(scale));
  const double saturationPressure =
      std::pow(10.0, 8.07131 - 1730.63 / (temperatureFromScale - 39.724)) * 133.32236842105263;
  const double relativePressure = pressure / saturationPressure;
  const double c = 0.155 * std::pow(scale, 6600.0);
  const double k = 0.871 * std::pow(scale, 0.0);
  const double expected =
      2.58 * c * k * relativePressure /
      ((1.0 - k * relativePressure) * (1.0 + k * relativePressure * (c - 1.0)));

  EXPECT_DOUBLE_EQ(gab.value(0.0, scale), 0.0);
  EXPECT_NEAR(gab.value(pressure, scale), expected, 1.0e-13);
  EXPECT_THROW(gab.value(pressure, 1.0), std::domain_error);

  const Isotherm gabWithTemperatureDependentK(Isotherm::Type::GAB, {2.58, 0.155, 0.8, 6600.0, -1200.0}, true);
  constexpr double secondTemperature = 308.15;
  constexpr double secondPressure = 1200.0;
  const double secondScale = std::exp(1.0 / (R * secondTemperature));
  const double secondTemperatureFromScale = 1.0 / (R * std::log(secondScale));
  const double secondSaturationPressure =
      std::pow(10.0, 8.07131 - 1730.63 / (secondTemperatureFromScale - 39.724)) * 133.32236842105263;
  const double secondRelativePressure = secondPressure / secondSaturationPressure;
  const double secondC = 0.155 * std::pow(secondScale, 6600.0);
  const double secondK = 0.8 * std::pow(secondScale, -1200.0);
  const double secondExpected =
      2.58 * secondC * secondK * secondRelativePressure /
      ((1.0 - secondK * secondRelativePressure) *
       (1.0 + secondK * secondRelativePressure * (secondC - 1.0)));
  EXPECT_NEAR(gabWithTemperatureDependentK.value(secondPressure, secondScale), secondExpected, 1.0e-13);
}

TEST(MixturePrediction, SPIPassesTemperatureScaleToGAB)
{
  const Isotherm gab(Isotherm::Type::GAB, {2.58, 0.155, 0.871, 6600.0, 0.0}, true);
  std::vector<Component> components;
  components.emplace_back(0, "H2O", std::vector<Isotherm>{gab}, 1.0, 1.0, 0.0);
  components.front().nonIsothermal = true;
  components.front().heatOfAdsorption = 1.0;
  components.emplace_back(1, "N2", std::vector<Isotherm>{}, 0.0, 0.0, 0.0, true);
  MixturePrediction prediction("GAB temperature", components, 1, 1, 298.15, 1.0e3, 1.0e5, 2, 0,
                               static_cast<size_t>(MixturePrediction::PredictionMethod::SPI), 0);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(2, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  const std::vector<double> gasFractions{1.0, 0.0};
  constexpr double pressure = 1500.0;
  double temperature = 308.15;

  prediction.predictMixture(gasFractions, pressure, adsorbedFractions, loadings, cachedPressure,
                            cachedGrandPotential, temperature);

  EXPECT_NEAR(loadings[0], gab.value(pressure, components.front().scale(temperature)), 1.0e-13);
  EXPECT_NE(loadings[0], gab.value(pressure, components.front().scale(298.15)));
}

TEST(MixturePrediction, BuildsEquilibriumPredictionFromChemisorptionSites)
{
  const Isotherm physical(Isotherm::Type::Langmuir, {1.0, 1.0e-5}, false);
  const Isotherm firstChemical(Isotherm::Type::Langmuir, {2.0, 2.0e-5}, false);
  const Isotherm secondChemical(Isotherm::Type::Henry, {1.0e-6}, false);

  std::vector<Component> components;
  components.emplace_back(0, "CO2", std::vector<Isotherm>{physical}, 1.0, 1.0, 0.0);
  Chemisorption firstSite;
  firstSite.type = Chemisorption::Type::FirstOrder;
  firstSite.isotherm = firstChemical;
  Chemisorption secondSite;
  secondSite.type = Chemisorption::Type::Avrami;
  secondSite.isotherm = secondChemical;
  components.front().chemisorption.add(firstSite);
  components.front().chemisorption.add(secondSite);
  components.emplace_back(1, "N2", std::vector<Isotherm>{}, 0.0, 0.0, 0.0, true);

  MixturePrediction physicalPrediction(
      "combined sorption", components, 1, 1, 298.15, 1.0e3, 1.0e5, 2, 0,
      static_cast<size_t>(MixturePrediction::PredictionMethod::SPI), 0);
  MixturePrediction chemicalPrediction =
      MixturePrediction::makeChemisorptionPrediction(physicalPrediction, components);

  ASSERT_EQ(chemicalPrediction.components[0].isotherm.sites.size(), 2U);
  EXPECT_FALSE(chemicalPrediction.components[0].isCarrierGas);
  EXPECT_TRUE(chemicalPrediction.components[1].isCarrierGas);

  std::vector<double> pureLoadings(2, 0.0);
  constexpr double pressure = 1.0e5;
  chemicalPrediction.predictPureComponentLoadings(pressure, pureLoadings, 298.15);
  EXPECT_NEAR(pureLoadings[0], firstChemical.value(pressure, 1.0) + secondChemical.value(pressure, 1.0), 1.0e-14);
  EXPECT_DOUBLE_EQ(pureLoadings[1], 0.0);
}

TEST(MacrostateParticleDistribution, ReweightsCOrderedTwoComponentDistribution)
{
  const std::filesystem::path path = mpdTemporaryPath("c_order.data");
  writeMPDFile(path, "0.1\n0.2\n0.3\n0.4\n");

  MPDSettings settings;
  settings.fileName = path.string();
  settings.referenceTemperature = 300.0;
  settings.referenceFugacity = 1.0e5;
  settings.referenceFrameworkMass = 1.0;
  settings.componentBounds = {{"A", 0, 1, 1}, {"B", 0, 1, 1}};
  MacrostateParticleDistribution distribution(settings);

  ASSERT_EQ(distribution.rank(), 2U);
  EXPECT_EQ(distribution.numberOfMacrostates(), 4U);
  EXPECT_FALSE(distribution.hasMeanEnergies());

  const std::vector<double> gasFractions{0.25, 0.75};
  const std::vector<double> means = distribution.meanParticleNumbers(gasFractions, 2.0e5, 300.0);
  ASSERT_EQ(means.size(), 2U);
  EXPECT_NEAR(means[0], 0.45 / 0.85, 1.0e-14);
  EXPECT_NEAR(means[1], 0.60 / 0.85, 1.0e-14);
}

TEST(MacrostateParticleDistribution, AppliesFirstOrderTemperatureReweightingWhenEnergyColumnExists)
{
  const std::filesystem::path path = mpdTemporaryPath("energy.data");
  constexpr double boltzmannConstant = 1.380649e-23;
  {
    std::ofstream output(path);
    output << "0.5 0\n";
    output << "0.5 " << std::scientific << boltzmannConstant * 600.0 << "\n";
  }

  MPDSettings settings;
  settings.fileName = path.string();
  settings.referenceTemperature = 300.0;
  settings.referenceFugacity = 1.0e5;
  settings.referenceFrameworkMass = 1.0;
  settings.componentBounds = {{"A", 0, 1, 1}};
  MacrostateParticleDistribution distribution(settings);

  ASSERT_TRUE(distribution.hasMeanEnergies());
  const std::vector<double> gasFractions{1.0};
  const std::vector<double> means = distribution.meanParticleNumbers(gasFractions, 1.0e5, 600.0);
  const double odds = 0.5 * std::exp(1.0);
  EXPECT_NEAR(means[0], odds / (1.0 + odds), 1.0e-7);
}

TEST(MacrostateParticleDistribution, UsesInclusiveBoundsAndParticleNumberDeltas)
{
  const std::filesystem::path path = mpdTemporaryPath("delta.data");
  writeMPDFile(path, "0.2\n0.3\n0.5\n");

  MPDSettings settings;
  settings.fileName = path.string();
  settings.referenceTemperature = 300.0;
  settings.referenceFugacity = 1.0e5;
  settings.referenceFrameworkMass = 1.0;
  settings.componentBounds = {{"A", 2, 6, 2}};
  MacrostateParticleDistribution distribution(settings);

  const std::vector<double> gasFractions{1.0};
  const std::vector<double> means = distribution.meanParticleNumbers(gasFractions, 1.0e5, 300.0);
  EXPECT_NEAR(means[0], 4.6, 1.0e-14);
}

TEST(MacrostateParticleDistribution, RejectsWrongEntryCount)
{
  const std::filesystem::path path = mpdTemporaryPath("wrong_count.data");
  writeMPDFile(path, "0.25\n0.25\n0.5\n");

  MPDSettings settings;
  settings.fileName = path.string();
  settings.referenceTemperature = 300.0;
  settings.referenceFugacity = 1.0e5;
  settings.referenceFrameworkMass = 1.0;
  settings.componentBounds = {{"A", 0, 1, 1}, {"B", 0, 1, 1}};
  EXPECT_THROW(static_cast<void>(MacrostateParticleDistribution(settings)), std::runtime_error);
}

TEST(MixturePrediction, LoadsAndUsesMPDSettings)
{
  const std::filesystem::path distributionPath = mpdTemporaryPath("integration.data");
  const std::filesystem::path inputPath = mpdTemporaryPath("integration.json");
  writeMPDFile(distributionPath, "0.1\n0.2\n0.3\n0.4\n");

  constexpr double avogadroConstant = 6.02214076e23;
  nlohmann::json input = {
      {"SimulationType", "MixturePrediction"},
      {"MixturePredictionMethod", "MPD"},
      {"Temperature", 300.0},
      {"PressureStart", 2.0e5},
      {"PressureEnd", 2.0e5},
      {"NumberOfPressurePoints", 1},
      {"MPDSettings",
       {{"FileName", distributionPath.filename().string()},
        {"ReferenceTemperature", 300.0},
        {"ReferenceFugacity", 1.0e5},
        {"ReferenceFrameworkMass", 1.0 / avogadroConstant},
        {"ComponentBounds",
         nlohmann::json::array({{{"Component", "A"}, {"NMin", 0}, {"NMax", 1}, {"DeltaN", 1}},
                                {{"Component", "B"}, {"NMin", 0}, {"NMax", 1}, {"DeltaN", 1}}})}}},
      {"Components",
       nlohmann::json::array({{{"Name", "A"}, {"GasPhaseMolFraction", 0.25}},
                              {{"Name", "B"}, {"GasPhaseMolFraction", 0.75}}})},
  };
  {
    std::ofstream output(inputPath);
    output << input;
  }

  InputReader reader(inputPath.string());
  ASSERT_TRUE(reader.mpdSettings.has_value());
  EXPECT_EQ(reader.mixturePredictionMethod, static_cast<size_t>(MixturePrediction::PredictionMethod::MPD));
  EXPECT_EQ(reader.maxIsothermTerms, 1U);
  EXPECT_EQ(reader.mpdSettings->fileName, distributionPath.string());

  MixturePrediction prediction(reader);
  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(2, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  const std::vector<double> gasFractions{0.25, 0.75};
  double temperature = 300.0;
  prediction.predictMixture(gasFractions, 2.0e5, adsorbedFractions, loadings, cachedPressure,
                            cachedGrandPotential, temperature);

  EXPECT_NEAR(loadings[0], 0.45 / 0.85, 1.0e-14);
  EXPECT_NEAR(loadings[1], 0.60 / 0.85, 1.0e-14);
  EXPECT_NEAR(adsorbedFractions[0] + adsorbedFractions[1], 1.0, 1.0e-14);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[0], loadings[0], 1.0e-14);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[1], loadings[1], 1.0e-14);
}

TEST(MixturePrediction, BeaMultinomialMPDMatchesCompetitiveMultisiteLangmuir)
{
  const std::filesystem::path exampleDirectory =
      std::filesystem::path(RUPTURA_SOURCE_DIR) / "examples" / "MPD" / "BEA-alkanes";
  InputReader mpdReader((exampleDirectory / "mpd" / "simulation.json").string());
  InputReader langmuirReader((exampleDirectory / "langmuir" / "simulation.json").string());
  MixturePrediction mpd(mpdReader);
  MixturePrediction langmuir(langmuirReader);

  std::vector<double> gasFractions;
  gasFractions.reserve(mpdReader.components.size());
  for (const Component& component : mpdReader.components)
  {
    gasFractions.push_back(component.initialGasMoleFraction);
  }

  const auto predict = [&gasFractions](MixturePrediction& prediction, double fugacity)
  {
    std::vector<double> adsorbedFractions(prediction.numberOfComponents, 0.0);
    std::vector<double> loadings(prediction.numberOfComponents, 0.0);
    std::vector<double> cachedPressure(prediction.numberOfComponents * prediction.maxIsothermTerms, 0.0);
    std::vector<double> cachedGrandPotential(prediction.maxIsothermTerms, 0.0);
    double temperature = 552.0;
    prediction.predictMixture(gasFractions, fugacity, adsorbedFractions, loadings, cachedPressure,
                              cachedGrandPotential, temperature);
    return loadings;
  };

  for (double fugacity : {1.0e3, 1.0e5, 1.0e7})
  {
    const std::vector<double> mpdLoadings = predict(mpd, fugacity);
    const std::vector<double> langmuirLoadings = predict(langmuir, fugacity);
    ASSERT_EQ(mpdLoadings.size(), langmuirLoadings.size());
    for (std::size_t component = 0; component < mpdLoadings.size(); ++component)
    {
      EXPECT_NEAR(mpdLoadings[component], langmuirLoadings[component], 2.0e-13)
          << "component " << component << " at fugacity " << fugacity;
    }
  }
}

TEST(MixturePrediction, EvaluatesMPDPureComponentsThroughOneHotMixtures)
{
  const std::filesystem::path distributionPath = mpdTemporaryPath("pure_components.data");
  const std::filesystem::path inputPath = mpdTemporaryPath("pure_components.json");
  writeMPDFile(distributionPath, "0.1\n0.2\n0.3\n0.4\n");

  constexpr double avogadroConstant = 6.02214076e23;
  nlohmann::json input = {
      {"SimulationType", "MixturePrediction"},
      {"MixturePredictionMethod", "MPD"},
      {"Temperature", 300.0},
      {"PressureStart", 1.0e5},
      {"PressureEnd", 1.0e5},
      {"NumberOfPressurePoints", 1},
      {"MPDSettings",
       {{"FileName", distributionPath.filename().string()},
        {"ReferenceTemperature", 300.0},
        {"ReferenceFugacity", 1.0e5},
        {"ReferenceFrameworkMass", 1.0 / avogadroConstant},
        {"ComponentBounds",
         nlohmann::json::array({{{"Component", "A"}, {"NMin", 0}, {"NMax", 1}, {"DeltaN", 1}},
                                {{"Component", "B"}, {"NMin", 0}, {"NMax", 1}, {"DeltaN", 1}}})}}},
      {"Components",
       nlohmann::json::array({{{"Name", "A"}, {"GasPhaseMolFraction", 0.25}},
                              {{"Name", "B"}, {"GasPhaseMolFraction", 0.75}}})},
  };
  {
    std::ofstream output(inputPath);
    output << input;
  }

  InputReader reader(inputPath.string());
  MixturePrediction prediction(reader);
  std::vector<double> pureLoadings(2, 0.0);
  prediction.predictPureComponentLoadings(1.0e5, pureLoadings, 300.0);

  // For pure A only states (0,0) and (1,0) remain; for pure B only
  // states (0,0) and (0,1) remain.
  EXPECT_NEAR(pureLoadings[0], 0.3 / (0.1 + 0.3), 1.0e-14);
  EXPECT_NEAR(pureLoadings[1], 0.2 / (0.1 + 0.2), 1.0e-14);
}

TEST(MPDInput, WarnsWhenReferenceTemperatureHasNoTargetTemperature)
{
  const std::filesystem::path distributionPath = mpdTemporaryPath("warning.data");
  const std::filesystem::path inputPath = mpdTemporaryPath("warning.json");
  writeMPDFile(distributionPath, "1.0\n");
  nlohmann::json input = {
      {"SimulationType", "MixturePrediction"},
      {"MixturePredictionMethod", "MPD"},
      {"MPDSettings",
       {{"FileName", distributionPath.string()},
        {"ReferenceTemperature", 300.0},
        {"ReferenceFugacity", 1.0e5},
        {"FrameworkMass", 1.0},
        {"ComponentBounds",
         nlohmann::json::array({{{"Component", "A"}, {"NMin", 0}, {"NMax", 0}, {"DeltaN", 1}}})}}},
      {"Components", nlohmann::json::array({{{"Name", "A"}, {"GasPhaseMolFraction", 1.0}}})},
  };
  {
    std::ofstream output(inputPath);
    output << input;
  }

  testing::internal::CaptureStderr();
  InputReader reader(inputPath.string());
  const std::string warning = testing::internal::GetCapturedStderr();
  EXPECT_NE(warning.find("ReferenceTemperature is set, but target Temperature is not"), std::string::npos);
}

TEST(MultiSiteIsotherm, RetainsButIgnoresZeroLoadingSites)
{
  const Isotherm inactive(Isotherm::Type::Langmuir, {0.0, 0.0}, false);
  const Isotherm active(Isotherm::Type::Langmuir, {1.5, 2.0e-5}, false);
  MultiSiteIsotherm isotherm({inactive, active});

  ASSERT_EQ(isotherm.sites.size(), 2U);
  EXPECT_FALSE(isotherm.sites[0].enabled());
  EXPECT_TRUE(isotherm.sites[1].enabled());
  EXPECT_DOUBLE_EQ(isotherm.value(size_t{0}, 1.0e5, 1.0), 0.0);
  EXPECT_DOUBLE_EQ(isotherm.psiForPressure(0, 1.0e5, 1.0), 0.0);

  double cachedPressure = 0.0;
  EXPECT_DOUBLE_EQ(inactive.inversePressureForPsi(1.0, cachedPressure, 1.0), 0.0);
  const double psi = active.psiForPressure(1.0e5, 1.0);
  EXPECT_DOUBLE_EQ(isotherm.inversePressureForPsi(psi, cachedPressure, 1.0),
                   active.inversePressureForPsi(psi, cachedPressure, 1.0));
}

TEST(MixturePrediction, SIASTKeepsSiteOrderingWithZeroLoadingSites)
{
  const std::vector<Isotherm> firstSites{
      Isotherm(Isotherm::Type::Langmuir, {0.0, 0.0}, false),
      Isotherm(Isotherm::Type::Langmuir, {1.0, 2.0e-5}, false),
  };
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}, false),
      Isotherm(Isotherm::Type::Langmuir, {0.0, 0.0}, false),
  };

  for (const size_t iastMethod : {size_t{0}, size_t{1}})
  {
    std::vector<Component> components;
    components.emplace_back(0, "first", firstSites, 0.5, 1.0, 0.0);
    components.emplace_back(1, "second", secondSites, 0.5, 1.0, 0.0);

    MixturePrediction prediction("zero-loading sites", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                                 static_cast<size_t>(MixturePrediction::PredictionMethod::SIAST), iastMethod);

    ASSERT_EQ(prediction.maxIsothermTerms, 2U);
    ASSERT_EQ(prediction.components[0].isotherm.sites.size(), 2U);
    ASSERT_EQ(prediction.components[1].isotherm.sites.size(), 2U);
    EXPECT_EQ(prediction.segregatedNumberOfSortedComponents, (std::vector<size_t>{1U, 1U}));
    ASSERT_EQ(prediction.segregatedSortedComponents[0][1].isotherm.sites.size(), 1U);
    ASSERT_EQ(prediction.segregatedSortedComponents[1][1].isotherm.sites.size(), 1U);
    EXPECT_FALSE(prediction.segregatedSortedComponents[0][1].isotherm.sites[0].enabled());
    EXPECT_FALSE(prediction.segregatedSortedComponents[1][1].isotherm.sites[0].enabled());

    std::vector<double> adsorbedFractions(2, 0.0);
    std::vector<double> loadings(2, 0.0);
    std::vector<double> cachedPressure(4, 0.0);
    std::vector<double> cachedGrandPotential(2, 0.0);
    const std::vector<double> gasFractions{0.5, 0.5};
    double temperature = 300.0;

    EXPECT_NO_THROW(prediction.predictMixture(gasFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure,
                                              cachedGrandPotential, temperature));
    EXPECT_TRUE(std::isfinite(loadings[0]));
    EXPECT_TRUE(std::isfinite(loadings[1]));
    EXPECT_GT(loadings[0], 0.0);
    EXPECT_GT(loadings[1], 0.0);

    ASSERT_EQ(prediction.equilibriumSiteLoadings.size(), 4U);
    EXPECT_DOUBLE_EQ(prediction.equilibriumSiteLoadings[0], 0.0);
    EXPECT_GT(prediction.equilibriumSiteLoadings[1], 0.0);
    EXPECT_GT(prediction.equilibriumSiteLoadings[2], 0.0);
    EXPECT_DOUBLE_EQ(prediction.equilibriumSiteLoadings[3], 0.0);
  }
}

TEST(MixturePrediction, SCICompetesWithinEachAlignedSite)
{
  const std::vector<Isotherm> firstSites{
      Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}, false),
      Isotherm(Isotherm::Type::Langmuir, {0.0, 0.0}, false),
  };
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Langmuir, {4.0, 2.0e-5}, false),
      Isotherm(Isotherm::Type::Langmuir, {1.0, 5.0e-6}, false),
  };

  std::vector<Component> components;
  components.emplace_back(0, "first", firstSites, 0.25, 1.0, 0.0);
  components.emplace_back(1, "second", secondSites, 0.75, 1.0, 0.0);

  MixturePrediction prediction("SCI", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                               static_cast<size_t>(MixturePrediction::PredictionMethod::SCI), 0);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(4, 0.0);
  std::vector<double> cachedGrandPotential(2, 0.0);
  const std::vector<double> gasFractions{0.25, 0.75};
  double temperature = 300.0;

  prediction.predictMixture(gasFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                            temperature);

  const double firstSiteDenominator = 1.0 + 1.0e-5 * 0.25 * 1.0e5 + 2.0e-5 * 0.75 * 1.0e5;
  const double expectedFirstSiteFirst = 2.0 * (1.0e-5 * 0.25 * 1.0e5) / firstSiteDenominator;
  const double expectedFirstSiteSecond = 4.0 * (2.0e-5 * 0.75 * 1.0e5) / firstSiteDenominator;
  const double secondSiteActivity = 5.0e-6 * 0.75 * 1.0e5;
  const double expectedSecondSiteSecond = secondSiteActivity / (1.0 + secondSiteActivity);

  EXPECT_NEAR(loadings[0], expectedFirstSiteFirst, 1.0e-14);
  EXPECT_NEAR(loadings[1], expectedFirstSiteSecond + expectedSecondSiteSecond, 1.0e-14);
  ASSERT_EQ(prediction.equilibriumSiteLoadings.size(), 4U);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[0], expectedFirstSiteFirst, 1.0e-14);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[1], expectedFirstSiteSecond, 1.0e-14);
  EXPECT_DOUBLE_EQ(prediction.equilibriumSiteLoadings[2], 0.0);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[3], expectedSecondSiteSecond, 1.0e-14);
  EXPECT_NEAR(adsorbedFractions[0] + adsorbedFractions[1], 1.0, 1.0e-14);
}

TEST(MixturePrediction, SPISumsIndependentPureSiteLoadings)
{
  const std::vector<Isotherm> firstSites{
      Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0}, false),
      Isotherm(Isotherm::Type::Henry, {2.0e-6}, false),
  };
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Langmuir_Freundlich, {3.0, 2.0e-5, 0.5}, false),
      Isotherm(Isotherm::Type::Toth, {1.5, 1.0e-5, 0.8}, false),
  };

  std::vector<Component> components;
  components.emplace_back(0, "first", firstSites, 0.25, 1.0, 0.0);
  components.emplace_back(1, "second", secondSites, 0.75, 1.0, 0.0);

  MixturePrediction prediction("SPI", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                               static_cast<size_t>(MixturePrediction::PredictionMethod::SPI), 0);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(4, 0.0);
  std::vector<double> cachedGrandPotential(2, 0.0);
  const std::vector<double> gasFractions{0.25, 0.75};
  double temperature = 300.0;

  prediction.predictMixture(gasFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                            temperature);

  const double firstPressure = gasFractions[0] * 1.0e5;
  const double secondPressure = gasFractions[1] * 1.0e5;
  const double expectedFirstSiteFirst = firstSites[0].value(firstPressure, 1.0);
  const double expectedSecondSiteFirst = firstSites[1].value(firstPressure, 1.0);
  const double expectedFirstSiteSecond = secondSites[0].value(secondPressure, 1.0);
  const double expectedSecondSiteSecond = secondSites[1].value(secondPressure, 1.0);

  EXPECT_NEAR(loadings[0], expectedFirstSiteFirst + expectedSecondSiteFirst, 1.0e-13);
  EXPECT_NEAR(loadings[1], expectedFirstSiteSecond + expectedSecondSiteSecond, 1.0e-13);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[0], expectedFirstSiteFirst, 1.0e-13);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[1], expectedFirstSiteSecond, 1.0e-13);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[2], expectedSecondSiteFirst, 1.0e-13);
  EXPECT_NEAR(prediction.equilibriumSiteLoadings[3], expectedSecondSiteSecond, 1.0e-13);
  EXPECT_NEAR(adsorbedFractions[0] + adsorbedFractions[1], 1.0, 1.0e-14);
}

TEST(MixturePrediction, SCISafelyHandlesNegativeFractionalPowerDrivingForce)
{
  const std::vector<Isotherm> firstSites{Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0}, false)};
  const std::vector<Isotherm> secondSites{Isotherm(Isotherm::Type::Sips, {3.0, 2.0e-5, 2.0}, false)};

  std::vector<Component> components;
  components.emplace_back(0, "first", firstSites, 0.5, 1.0, 0.0);
  components.emplace_back(1, "second", secondSites, 0.5, 1.0, 0.0);
  MixturePrediction prediction("safe SCI", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                               static_cast<size_t>(MixturePrediction::PredictionMethod::SCI), 0);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(2, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  const std::vector<double> gasFractions{-1.0e-12, 1.0 + 1.0e-12};
  double temperature = 300.0;

  EXPECT_NO_THROW(prediction.predictMixture(gasFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure,
                                            cachedGrandPotential, temperature));
  EXPECT_DOUBLE_EQ(loadings[0], 0.0);
  EXPECT_TRUE(std::isfinite(loadings[1]));
  EXPECT_GT(loadings[1], 0.0);
  EXPECT_TRUE(std::isfinite(adsorbedFractions[0]));
  EXPECT_TRUE(std::isfinite(adsorbedFractions[1]));
}

TEST(MixturePrediction, SCISupportedModelsReduceToPureSiteLoadingWithoutCompetitors)
{
  const std::vector<Isotherm> supportedModels{
      Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}, false),
      Isotherm(Isotherm::Type::Anti_Langmuir, {2.0e-5, 1.0e-7}, false),
      Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0}, false),
      Isotherm(Isotherm::Type::Langmuir_Freundlich, {2.0, 1.0e-5, 0.75}, false),
      Isotherm(Isotherm::Type::Redlich_Peterson, {2.0e-5, 1.0e-5, 0.75}, false),
      Isotherm(Isotherm::Type::Toth, {2.0, 1.0e-5, 0.75}, false),
  };

  for (const Isotherm& isotherm : supportedModels)
  {
    std::vector<Component> components;
    components.emplace_back(0, "present", std::vector<Isotherm>{isotherm}, 1.0, 1.0, 0.0);
    components.emplace_back(1, "absent", std::vector<Isotherm>{isotherm}, 0.0, 1.0, 0.0);
    MixturePrediction prediction("SCI pure limit", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                                 static_cast<size_t>(MixturePrediction::PredictionMethod::SCI), 0);

    std::vector<double> adsorbedFractions(2, 0.0);
    std::vector<double> loadings(2, 0.0);
    std::vector<double> cachedPressure(2, 0.0);
    std::vector<double> cachedGrandPotential(1, 0.0);
    const std::vector<double> gasFractions{1.0, 0.0};
    double temperature = 300.0;

    prediction.predictMixture(gasFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                              temperature);

    const double expected = isotherm.value(1.0e5, 1.0);
    EXPECT_NEAR(loadings[0], expected, 1.0e-12 * std::max(1.0, std::abs(expected)));
    EXPECT_DOUBLE_EQ(loadings[1], 0.0);
  }
}

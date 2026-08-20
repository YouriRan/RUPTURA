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

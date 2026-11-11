#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "integrators/sorption.h"
#include "isotherm.h"
#include "json.h"
#include "macrostate_particle_distribution.h"
#include "mixture_prediction.h"
#include "multi_site_isotherm.h"
#include "temporary_directory.h"

namespace
{
std::filesystem::path mpdTemporaryPath(const std::string& name)
{
  return testTemporaryDirectory() / ("ruptura_mpd_" + name);
}

void writeMPDFile(const std::filesystem::path& path, const std::string& contents)
{
  std::ofstream output(path);
  output << contents;
}
}  // namespace

TEST(Isotherm, GABImplementsElfvingEquationsFiveThroughSeven)
{
  // C = C_0 exp(H_C/(R T)) as an uncentred VantHoff law; p_sat follows the Antoine equation at the temperature.
  Isotherm gab(Isotherm::Type::GAB, {2.58, 0.155, 0.871},
               {{1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 0.155, 6600.0)}});
  EXPECT_THROW(gab.value(1000.0), std::domain_error);

  constexpr double temperature = 298.15;
  constexpr double pressure = 0.019 * 105000.0;
  gab.setTemperature(temperature);
  const double saturationPressure = std::pow(10.0, 8.07131 - 1730.63 / (temperature - 39.724)) * 133.32236842105263;
  const double relativePressure = pressure / saturationPressure;
  const double c = 0.155 * std::exp(6600.0 / (R * temperature));
  const double k = 0.871;
  const double expected =
      2.58 * c * k * relativePressure /
      ((1.0 - k * relativePressure) * (1.0 + k * relativePressure * (c - 1.0)));

  EXPECT_DOUBLE_EQ(gab.value(0.0), 0.0);
  EXPECT_NEAR(gab.value(pressure), expected, 1.0e-13);

  Isotherm gabWithTemperatureDependentK(Isotherm::Type::GAB, {2.58, 0.155, 0.8},
                                        {{1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 0.155, 6600.0)},
                                         {2, TemperatureLaw(TemperatureLaw::Type::VantHoff, 0.8, -1200.0)}});
  constexpr double secondTemperature = 308.15;
  constexpr double secondPressure = 1200.0;
  gabWithTemperatureDependentK.setTemperature(secondTemperature);
  const double secondSaturationPressure =
      std::pow(10.0, 8.07131 - 1730.63 / (secondTemperature - 39.724)) * 133.32236842105263;
  const double secondRelativePressure = secondPressure / secondSaturationPressure;
  const double secondC = 0.155 * std::exp(6600.0 / (R * secondTemperature));
  const double secondK = 0.8 * std::exp(-1200.0 / (R * secondTemperature));
  const double secondExpected =
      2.58 * secondC * secondK * secondRelativePressure /
      ((1.0 - secondK * secondRelativePressure) *
       (1.0 + secondK * secondRelativePressure * (secondC - 1.0)));
  EXPECT_NEAR(gabWithTemperatureDependentK.value(secondPressure), secondExpected, 1.0e-13);
}

TEST(Isotherm, UnilanLoadingAndSpreadingPressureFollowTheUnilanForm)
{
  constexpr double saturationLoading = 3.0;
  constexpr double affinity = 2.0e-5;
  constexpr double heterogeneity = 1.5;
  const Isotherm unilan(Isotherm::Type::Unilan, {saturationLoading, affinity, heterogeneity});

  constexpr double pressure = 4.0e4;
  const double expected = saturationLoading / (2.0 * heterogeneity) *
                          std::log((1.0 + affinity * std::exp(heterogeneity) * pressure) /
                                   (1.0 + affinity * std::exp(-heterogeneity) * pressure));
  EXPECT_NEAR(unilan.value(pressure), expected, 1.0e-13);
  EXPECT_GT(unilan.value(pressure), 0.0);

  // psi(p) = integral of q(p')/p' dp' from 0 to p, here by the midpoint rule in ln p'.
  constexpr size_t intervals = 200000;
  const double lowerLogPressure = std::log(pressure) - 40.0;
  const double width = (std::log(pressure) - lowerLogPressure) / static_cast<double>(intervals);
  double integral = 0.0;
  for (size_t i = 0; i < intervals; ++i)
  {
    integral += unilan.value(std::exp(lowerLogPressure + (static_cast<double>(i) + 0.5) * width)) * width;
  }
  EXPECT_NEAR(unilan.psiForPressure(pressure), integral, 1.0e-8 * integral);
}

TEST(Isotherm, AntiLangmuirInversePressureInvertsTheSpreadingPressure)
{
  const Isotherm antiLangmuir(Isotherm::Type::Anti_Langmuir, {2.0, 4.0e-6});
  double cachedPressure = 0.0;
  for (double pressure : {1.0e3, 5.0e4, 2.0e5})
  {
    const double psi = antiLangmuir.psiForPressure(pressure);
    EXPECT_NEAR(antiLangmuir.inversePressureForPsi(psi, cachedPressure), 1.0 / pressure, 1.0e-12 / pressure);
  }
}

TEST(Isotherm, OBrienMyersLoadingAndSpreadingPressureFollowTheOBrienMyersForm)
{
  constexpr double saturationLoading = 2.0;
  constexpr double affinity = 3.0e-5;
  constexpr double sigma = 0.8;
  const Isotherm obrienMyers(Isotherm::Type::OBrien_Myers, {saturationLoading, affinity, sigma});

  constexpr double pressure = 4.0e4;
  const double kp = affinity * pressure;
  const double expected = saturationLoading * (kp / (1.0 + kp) + sigma * sigma * kp * (1.0 - kp) /
                                                                      (2.0 * std::pow(1.0 + kp, 3)));
  EXPECT_NEAR(obrienMyers.value(pressure), expected, 1.0e-13);

  // psi(p) = integral of q(p')/p' dp' from 0 to p, here by the midpoint rule in ln p'.
  constexpr size_t intervals = 200000;
  const double lowerLogPressure = std::log(pressure) - 40.0;
  const double width = (std::log(pressure) - lowerLogPressure) / static_cast<double>(intervals);
  double integral = 0.0;
  for (size_t i = 0; i < intervals; ++i)
  {
    integral += obrienMyers.value(std::exp(lowerLogPressure + (static_cast<double>(i) + 0.5) * width)) * width;
  }
  EXPECT_NEAR(obrienMyers.psiForPressure(pressure), integral, 1.0e-8 * integral);
}

TEST(Isotherm, FreundlichPhysicalityChecksOnlyItsTwoParameters)
{
  EXPECT_FALSE(Isotherm(Isotherm::Type::Freundlich, {0.5, 2.0}).isUnphysical());
  EXPECT_TRUE(Isotherm(Isotherm::Type::Freundlich, {-0.5, 2.0}).isUnphysical());
  EXPECT_TRUE(Isotherm(Isotherm::Type::Freundlich, {0.5, -2.0}).isUnphysical());
}

TEST(MixturePrediction, SPIAndSCIEvaluateTothAtTheTemperature)
{
  // VantHoff affinities b(T) = b_ref exp[Q/R (1/T - 1/T_ref)] with T_ref = 300 K.
  const std::vector<Isotherm> firstSites{
      Isotherm(Isotherm::Type::Toth, {2.0, 1.0e-5, 0.6},
               {{1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 1.0e-5, 30000.0, 300.0)}})};
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Toth, {1.5, 4.0e-6, 0.6},
               {{1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 4.0e-6, 20000.0, 300.0)}})};
  std::vector<Component> components;
  components.emplace_back(0, "first", firstSites, 0.25, 1.0, 0.0);
  components.emplace_back(1, "second", secondSites, 0.75, 1.0, 0.0);

  constexpr double pressure = 1.0e5;
  const std::vector<double> gasFractions{0.25, 0.75};
  double temperature = 340.0;
  const auto scale = [&](double heat) { return std::exp(heat / R * (1.0 / temperature - 1.0 / 300.0)); };
  const std::array<double, 2> scales{scale(30000.0), scale(20000.0)};
  ASSERT_LT(scales[0], 0.5);
  ASSERT_LT(scales[1], 0.5);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(2, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);

  MixturePrediction spi("SPI Toth", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                        static_cast<size_t>(MixturePrediction::PredictionMethod::SPI), 0);
  spi.predictMixture(gasFractions, pressure, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                     temperature);
  const std::array<double, 2> activities{scales[0] * 1.0e-5 * gasFractions[0] * pressure,
                                         scales[1] * 4.0e-6 * gasFractions[1] * pressure};
  EXPECT_NEAR(loadings[0], 2.0 * activities[0] / std::pow(1.0 + std::pow(activities[0], 0.6), 1.0 / 0.6), 1.0e-13);
  EXPECT_NEAR(loadings[1], 1.5 * activities[1] / std::pow(1.0 + std::pow(activities[1], 0.6), 1.0 / 0.6), 1.0e-13);

  MixturePrediction sci("SCI Toth", components, 0, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                        static_cast<size_t>(MixturePrediction::PredictionMethod::SCI), 0);
  sci.predictMixture(gasFractions, pressure, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                     temperature);
  const double denominator =
      std::pow(1.0 + std::pow(activities[0], 0.6) + std::pow(activities[1], 0.6), 1.0 / 0.6);
  EXPECT_NEAR(loadings[0], 2.0 * activities[0] / denominator, 1.0e-13);
  EXPECT_NEAR(loadings[1], 1.5 * activities[1] / denominator, 1.0e-13);
}

TEST(MixturePrediction, IASTEvaluatesSaturationLoadingAndExponentLaws)
{
  // q_sat(T) = 3 exp[0.4 (1 - T/300)], b(T) VantHoff, and t(T) = 0.7 + 0.2 (1 - 300/T): the IAST loading of a
  // single component equals its pure Toth isotherm at the temperature.
  const Isotherm toth(Isotherm::Type::Toth, {3.0, 2.0e-5, 0.7},
                      {{0, TemperatureLaw(TemperatureLaw::Type::Exponential, 3.0, 0.4, 300.0)},
                       {1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 2.0e-5, 25000.0, 300.0)},
                       {2, TemperatureLaw(TemperatureLaw::Type::Inverse, 0.7, 0.2, 300.0)}});
  std::vector<Component> components;
  components.emplace_back(0, "CO2", std::vector<Isotherm>{toth}, 0.5, 1.0, 0.0);
  components.emplace_back(1, "N2", std::vector<Isotherm>{}, 0.5, 0.0, 0.0, true);
  MixturePrediction iast("IAST Toth", components, 1, 1, 300.0, 1.0e3, 1.0e5, 2, 0,
                         static_cast<size_t>(MixturePrediction::PredictionMethod::IAST), 0);

  double temperature = 330.0;
  const double qsat = 3.0 * std::exp(0.4 * (1.0 - temperature / 300.0));
  const double b = 2.0e-5 * std::exp(25000.0 / R * (1.0 / temperature - 1.0 / 300.0));
  const double t = 0.7 + 0.2 * (1.0 - 300.0 / temperature);
  const double partialPressure = 0.5 * 2.0e5;
  const double expected = qsat * b * partialPressure / std::pow(1.0 + std::pow(b * partialPressure, t), 1.0 / t);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(2, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  const std::vector<double> gasFractions{0.5, 0.5};
  iast.predictMixture(gasFractions, 2.0e5, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                      temperature);
  EXPECT_NEAR(loadings[0], expected, 1.0e-12 * expected);

  // At the reference temperature the laws return their values.
  double referenceTemperature = 300.0;
  iast.predictMixture(gasFractions, 2.0e5, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                      referenceTemperature);
  EXPECT_NEAR(loadings[0], Isotherm(Isotherm::Type::Toth, {3.0, 2.0e-5, 0.7}).value(partialPressure), 1.0e-13);
}

TEST(MixturePrediction, AcceptsRoundingInTheMoleFractionSumAndRejectsAWrongSum)
{
  // Normalized mole fractions sum to 1 only up to rounding, which grows with the number of components; a sum off by
  // 5e-13 is accepted, one of 0.999 is not (B-PHY-13).
  std::vector<Component> components;
  components.emplace_back(0, "CO2", std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {3.0, 2.0e-5})}, 0.5,
                          1.0, 0.0);
  components.emplace_back(1, "N2", std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-6})}, 0.5,
                          1.0, 0.0);
  MixturePrediction iast("IAST Langmuir", components, 1, 1, 300.0, 1.0e3, 1.0e5, 2, 0,
                         static_cast<size_t>(MixturePrediction::PredictionMethod::IAST), 0);

  std::vector<double> adsorbedFractions(2, 0.0);
  std::vector<double> loadings(2, 0.0);
  std::vector<double> cachedPressure(2, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  double temperature = 300.0;
  const std::vector<double> roundedFractions{0.5 + 5.0e-13, 0.5};
  EXPECT_NO_THROW(iast.predictMixture(roundedFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure,
                                      cachedGrandPotential, temperature));

  const std::vector<double> wrongFractions{0.5, 0.499};
  EXPECT_THROW(
      {
        try
        {
          iast.predictMixture(wrongFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure,
                              cachedGrandPotential, temperature);
        }
        catch (const std::runtime_error& error)
        {
          EXPECT_NE(std::string(error.what()).find("sum to 0.999 instead of 1"), std::string::npos) << error.what();
          throw;
        }
      },
      std::runtime_error);
}

TEST(MixturePrediction, SPIEvaluatesGABAtTheTemperature)
{
  const Isotherm gab(Isotherm::Type::GAB, {2.58, 0.155, 0.871},
                     {{1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 0.155, 6600.0)}});
  std::vector<Component> components;
  components.emplace_back(0, "H2O", std::vector<Isotherm>{gab}, 1.0, 1.0, 0.0);
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

  Isotherm reference = gab;
  reference.setTemperature(temperature);
  EXPECT_NEAR(loadings[0], reference.value(pressure), 1.0e-13);
  reference.setTemperature(298.15);
  EXPECT_NE(loadings[0], reference.value(pressure));
}

TEST(MixturePrediction, BuildsEquilibriumPredictionFromChemisorptionSites)
{
  const Isotherm physical(Isotherm::Type::Langmuir, {1.0, 1.0e-5});
  const Isotherm firstChemical(Isotherm::Type::Langmuir, {2.0, 2.0e-5});
  const Isotherm secondChemical(Isotherm::Type::Henry, {1.0e-6});

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
  EXPECT_NEAR(pureLoadings[0], firstChemical.value(pressure) + secondChemical.value(pressure), 1.0e-14);
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
  const Isotherm inactive(Isotherm::Type::Langmuir, {0.0, 0.0});
  const Isotherm active(Isotherm::Type::Langmuir, {1.5, 2.0e-5});
  MultiSiteIsotherm isotherm({inactive, active});

  ASSERT_EQ(isotherm.sites.size(), 2U);
  EXPECT_FALSE(isotherm.sites[0].enabled());
  EXPECT_TRUE(isotherm.sites[1].enabled());
  EXPECT_DOUBLE_EQ(isotherm.value(size_t{0}, 1.0e5), 0.0);
  EXPECT_DOUBLE_EQ(isotherm.psiForPressure(0, 1.0e5), 0.0);

  double cachedPressure = 0.0;
  EXPECT_DOUBLE_EQ(inactive.inversePressureForPsi(1.0, cachedPressure), 0.0);
  const double psi = active.psiForPressure(1.0e5);
  EXPECT_DOUBLE_EQ(isotherm.inversePressureForPsi(psi, cachedPressure),
                   active.inversePressureForPsi(psi, cachedPressure));
}

TEST(MixturePrediction, SIASTKeepsSiteOrderingWithZeroLoadingSites)
{
  const std::vector<Isotherm> firstSites{
      Isotherm(Isotherm::Type::Langmuir, {0.0, 0.0}),
      Isotherm(Isotherm::Type::Langmuir, {1.0, 2.0e-5}),
  };
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}),
      Isotherm(Isotherm::Type::Langmuir, {0.0, 0.0}),
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
      Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}),
      Isotherm(Isotherm::Type::Langmuir, {0.0, 0.0}),
  };
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Langmuir, {4.0, 2.0e-5}),
      Isotherm(Isotherm::Type::Langmuir, {1.0, 5.0e-6}),
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
      Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0}),
      Isotherm(Isotherm::Type::Henry, {2.0e-6}),
  };
  const std::vector<Isotherm> secondSites{
      Isotherm(Isotherm::Type::Langmuir_Freundlich, {3.0, 2.0e-5, 0.5}),
      Isotherm(Isotherm::Type::Toth, {1.5, 1.0e-5, 0.8}),
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
  const double expectedFirstSiteFirst = firstSites[0].value(firstPressure);
  const double expectedSecondSiteFirst = firstSites[1].value(firstPressure);
  const double expectedFirstSiteSecond = secondSites[0].value(secondPressure);
  const double expectedSecondSiteSecond = secondSites[1].value(secondPressure);

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
  const std::vector<Isotherm> firstSites{Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0})};
  const std::vector<Isotherm> secondSites{Isotherm(Isotherm::Type::Sips, {3.0, 2.0e-5, 2.0})};

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
      Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}),
      Isotherm(Isotherm::Type::Anti_Langmuir, {2.0e-5, 1.0e-7}),
      Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0}),
      Isotherm(Isotherm::Type::Langmuir_Freundlich, {2.0, 1.0e-5, 0.75}),
      Isotherm(Isotherm::Type::Redlich_Peterson, {2.0e-5, 1.0e-5, 0.75}),
      Isotherm(Isotherm::Type::Toth, {2.0, 1.0e-5, 0.75}),
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

    const double expected = isotherm.value(1.0e5);
    EXPECT_NEAR(loadings[0], expected, 1.0e-12 * std::max(1.0, std::abs(expected)));
    EXPECT_DOUBLE_EQ(loadings[1], 0.0);
  }
}

TEST(MixturePrediction, AdsorbingCarrierParticipatesInEveryAnalyticalMethod)
{
  const Isotherm physical(Isotherm::Type::Langmuir, {2.0, 1.0e-5});
  for (size_t method = 0; method <= 5; ++method)
  {
    SCOPED_TRACE(method);
    std::vector<Component> components;
    components.emplace_back(0, "carrier", std::vector<Isotherm>{physical}, 0.8, 0.1, 0.0, true);
    components.emplace_back(1, "adsorbate", std::vector<Isotherm>{physical}, 0.2, 0.1, 0.0);
    MixturePrediction prediction("adsorbing carrier", components, 1, 0, 300, 1, 1e5, 2, 0, method, 0);
    std::vector<double> fractions(2), loadings(2), cachedPressure(2), cachedPotential(1);
    double temperature = 300;
    prediction.predictMixture(std::vector<double>{0.8, 0.2}, 1e5, fractions, loadings, cachedPressure,
                              cachedPotential, temperature);
    EXPECT_NEAR(loadings[0], method == 5 ? 2.0 * 0.8 / 1.8 : 0.8, 1e-10);
    EXPECT_NEAR(loadings[1], method == 5 ? 2.0 * 0.2 / 1.2 : 0.2, 1e-10);
  }
}

TEST(MixturePrediction, CarrierWithOnlyChemisorptionHasChemicalEquilibrium)
{
  std::vector<Component> components;
  components.emplace_back(0, "carrier", std::vector<Isotherm>{}, 0.8, 0.1, 0.0, true);
  components.emplace_back(1, "adsorbate", std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {2.0, 1e-5})},
                          0.2, 0.1, 0.0);
  Chemisorption site;
  site.type = Chemisorption::Type::FirstOrder;
  site.isotherm = Isotherm(Isotherm::Type::Langmuir, {1.0, 2e-5});
  components[0].chemisorption.add(site);
  MixturePrediction physical("chemical carrier", components, 1, 0, 300, 1, 1e5, 2, 0, 5, 0);
  EXPECT_TRUE(physical.hasChemisorptionEquilibrium());
  auto chemical = MixturePrediction::makeChemisorptionPrediction(physical, components);
  std::vector<double> fractions(2), loadings(2), cachedPressure(2), cachedPotential(1);
  double temperature = 300;
  chemical.predictMixture(std::vector<double>{0.8, 0.2}, 1e5, fractions, loadings, cachedPressure,
                           cachedPotential, temperature);
  EXPECT_NEAR(loadings[0], 1.6 / 2.6, 1e-12);
  EXPECT_DOUBLE_EQ(loadings[1], 0.0);
}

namespace
{
// Central differences of the predictMixture loadings with respect to the partial pressure of every adsorbing
// component at a fixed total pressure. The carrier (component 0) takes up the change in mole fraction; the steps are
// y_j / 4096, so the mole fractions keep summing to exactly one. Column-major, like predictMixtureSensitivities.
std::vector<double> centralDifferenceSensitivities(MixturePrediction& prediction, const std::vector<double>& fractions,
                                                   double pressure)
{
  const size_t n = fractions.size();
  std::vector<double> adsorbed(n), loadings(n), cachedPressure(n * prediction.maxIsothermTerms);
  std::vector<double> cachedPotential(prediction.maxIsothermTerms);
  std::vector<double> result(n * n, 0.0);
  for (size_t j = 1; j < n; ++j)
  {
    if (fractions[j] == 0.0) continue;
    const double step = std::ldexp(fractions[j], -12);
    std::array<std::vector<double>, 2> perturbedLoadings;
    for (size_t side = 0; side < 2; ++side)
    {
      std::vector<double> perturbed = fractions;
      const double signedStep = side == 0 ? step : -step;
      perturbed[j] += signedStep;
      perturbed[0] -= signedStep;
      double temperature = 300.0;
      std::fill(cachedPotential.begin(), cachedPotential.end(), 0.0);
      prediction.predictMixture(perturbed, pressure, adsorbed, loadings, cachedPressure, cachedPotential, temperature);
      perturbedLoadings[side] = loadings;
    }
    for (size_t i = 0; i < n; ++i)
    {
      result[i + j * n] = (perturbedLoadings[0][i] - perturbedLoadings[1][i]) / (2.0 * step * pressure);
    }
  }
  return result;
}

// Largest entrywise difference of the analytical sensitivities from central differences, relative to the largest
// central-difference entry.
double relativeSensitivityError(MixturePrediction& prediction, const std::vector<double>& fractions, double pressure)
{
  const size_t n = fractions.size();
  std::vector<double> adsorbed(n), loadings(n), cachedPressure(n * prediction.maxIsothermTerms);
  std::vector<double> cachedPotential(prediction.maxIsothermTerms), sensitivities(n * n);
  double temperature = 300.0;
  prediction.predictMixture(fractions, pressure, adsorbed, loadings, cachedPressure, cachedPotential, temperature);
  prediction.predictMixtureSensitivities(fractions, pressure, cachedPressure, temperature, 7.0, sensitivities);

  const std::vector<double> reference = centralDifferenceSensitivities(prediction, fractions, pressure);
  double scale = 0.0;
  double error = 0.0;
  for (size_t entry = 0; entry < n * n; ++entry)
  {
    scale = std::max(scale, std::abs(reference[entry]));
    error = std::max(error, std::abs(sensitivities[entry] - reference[entry]));
  }
  return error / scale;
}

// A non-adsorbing carrier and three adsorbing components with the given sites.
std::vector<Component> carrierAndThreeAdsorbates(const std::array<std::vector<Isotherm>, 3>& sites)
{
  std::vector<Component> components;
  components.emplace_back(0, "carrier", std::vector<Isotherm>{}, 0.5, 0.0, 0.0, true);
  for (size_t i = 0; i < 3; ++i)
  {
    components.emplace_back(i + 1, "adsorbate " + std::to_string(i + 1), sites[i], 0.0, 0.06, 0.0);
  }
  return components;
}

// A feed-like and a front-toe composition (carrier first), both summing to exactly one.
const std::array<std::vector<double>, 2> sensitivityCompositions{
    std::vector<double>{0.5, 0.25, 0.125, 0.125},
    std::vector<double>{1.0 - 0x1p-10 - 0x1p-20 - 0x1p-30, 0x1p-10, 0x1p-20, 0x1p-30}};
}  // namespace

TEST(Isotherm, SlopeMatchesCentralDifferencesOfTheLoading)
{
  Isotherm gab(Isotherm::Type::GAB, {2.58, 0.155, 0.871},
               {{1, TemperatureLaw(TemperatureLaw::Type::VantHoff, 0.155, 6600.0)}});
  gab.setTemperature(298.15);
  const std::vector<std::pair<Isotherm, std::vector<double>>> cases{
      {Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}), {1.0e2, 1.0e4, 1.0e6}},
      {Isotherm(Isotherm::Type::Anti_Langmuir, {2.0e-5, 1.0e-7}), {1.0e2, 1.0e4, 1.0e6}},
      {Isotherm(Isotherm::Type::BET, {1.5, 1.0e-4, 1.0e-6}), {1.0e2, 1.0e4, 2.0e5}},
      {Isotherm(Isotherm::Type::Henry, {2.0e-6}), {1.0e2, 1.0e5}},
      {Isotherm(Isotherm::Type::Freundlich, {1.0e-3, 2.0}), {1.0e-2, 1.0e2, 1.0e5}},
      {Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 2.0}), {1.0e-2, 1.0e2, 1.0e5}},
      {Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 0.7}), {1.0e-2, 1.0e2, 1.0e5}},
      {Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.4, 0.0024, 0.8}), {1.0e-4, 1.0e2, 2.5e4}},
      {Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.9, 0.00298, 1.2}), {1.0e-4, 1.0e2, 2.5e4}},
      {Isotherm(Isotherm::Type::Redlich_Peterson, {2.0e-5, 1.0e-5, 0.75}), {1.0e-2, 1.0e2, 1.0e6}},
      {Isotherm(Isotherm::Type::Toth, {1.5, 1.0e-5, 0.8}), {1.0e-2, 1.0e2, 1.0e6}},
      // value() of Unilan and Bingel-Walton cancels below p ~ 1 Pa (log(1 + x) / exp(-x) differences).
      {Isotherm(Isotherm::Type::Unilan, {3.0, 2.0e-5, 1.5}), {1.0e2, 1.0e4, 1.0e6}},
      {Isotherm(Isotherm::Type::OBrien_Myers, {2.0, 1.0e-5, 0.5}), {1.0e-2, 1.0e5, 1.0e6}},
      {Isotherm(Isotherm::Type::Quadratic, {2.0, 1.0e-5, 1.0e-10}), {1.0e-2, 1.0e5, 1.0e6}},
      {Isotherm(Isotherm::Type::Temkin, {2.0, 1.0e-5, 0.3}), {1.0e-2, 1.0e5, 1.0e6}},
      {Isotherm(Isotherm::Type::BingelWalton, {2.0, 1.0e-5, 2.0e-6}), {1.0e2, 1.0e5, 1.0e6}},
      {Isotherm(Isotherm::Type::Langmuir_pH, {2.0, 0.5, 7.0}), {1.0e-3, 0.5, 2.0}},
      {gab, {1.0, 500.0, 1500.0}},
  };

  for (const auto& [isotherm, pressures] : cases)
  {
    for (const double pressure : pressures)
    {
      SCOPED_TRACE(isotherm.repr() + " at p = " + std::to_string(pressure));
      const double step = 1.0e-5 * pressure;
      const double reference =
          (isotherm.value(pressure + step, 6.5) - isotherm.value(pressure - step, 6.5)) / (2.0 * step);
      EXPECT_NEAR(isotherm.slope(pressure, 6.5), reference, 1.0e-7 * std::abs(reference));
    }
  }
}

TEST(MixturePrediction, SensitivitiesMatchCentralDifferencesForEveryAnalyticalMethod)
{
  using Method = MixturePrediction::PredictionMethod;
  // JUC-77 xylene sites: Langmuir-Freundlich exponents below one have unbounded slopes at p -> 0.
  const std::array<std::vector<Isotherm>, 3> xyleneSites{
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.4, 0.0024, 0.8}),
                            Isotherm(Isotherm::Type::Langmuir, {0.9, 0.0201})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.53, 0.0124, 0.822}),
                            Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.77, 0.031, 1.1})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.4, 0.000702, 0.9}),
                            Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.9, 0.00298, 1.2})}};
  const std::array<std::vector<Isotherm>, 3> langmuirSites{
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}),
                            Isotherm(Isotherm::Type::Langmuir, {1.0, 3.0e-6})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {3.0, 2.0e-5}),
                            Isotherm(Isotherm::Type::Langmuir, {0.5, 1.0e-6})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {1.5, 5.0e-6}),
                            Isotherm(Isotherm::Type::Langmuir, {2.5, 2.0e-6})}};
  const std::array<std::vector<Isotherm>, 3> mixedSites{
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 0.7}),
                            Isotherm(Isotherm::Type::Henry, {2.0e-6})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir_Freundlich, {3.0, 2.0e-5, 0.5}),
                            Isotherm(Isotherm::Type::Toth, {1.5, 1.0e-5, 0.8})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Redlich_Peterson, {2.0e-5, 1.0e-5, 0.75}),
                            Isotherm(Isotherm::Type::Langmuir, {1.0, 4.0e-6})}};

  struct Case
  {
    std::string label;
    std::array<std::vector<Isotherm>, 3> sites;
    Method method;
    size_t iastMethod;
  };
  std::vector<Case> cases{
      {"IAST fast", xyleneSites, Method::IAST, 0},          {"IAST bisection", xyleneSites, Method::IAST, 1},
      {"SIAST fast", xyleneSites, Method::SIAST, 0},        {"SIAST bisection", xyleneSites, Method::SIAST, 1},
      {"EI", langmuirSites, Method::EI, 0},                 {"SEI", langmuirSites, Method::SEI, 0},
      {"SPI", mixedSites, Method::SPI, 0},                  {"SCI Langmuir", langmuirSites, Method::SCI, 0},
  };
  for (const Isotherm& model : {Isotherm(Isotherm::Type::Langmuir_Freundlich, {2.0, 1.0e-4, 0.8}),
                                Isotherm(Isotherm::Type::Sips, {2.0, 1.0e-5, 0.7}),
                                Isotherm(Isotherm::Type::Anti_Langmuir, {2.0e-5, 1.0e-6}),
                                Isotherm(Isotherm::Type::Toth, {1.5, 1.0e-5, 0.8}),
                                Isotherm(Isotherm::Type::Redlich_Peterson, {2.0e-5, 1.0e-5, 0.75})})
  {
    std::array<std::vector<Isotherm>, 3> sites;
    for (size_t i = 0; i < 3; ++i)
    {
      Isotherm scaled = model;
      scaled.parameters[0] *= 1.0 + 0.5 * static_cast<double>(i);
      scaled.parameters[1] *= 1.0 + 0.25 * static_cast<double>(i);
      sites[i] = {scaled};
    }
    cases.push_back({"SCI " + model.repr(), sites, Method::SCI, 0});
  }

  for (const Case& test : cases)
  {
    MixturePrediction prediction(test.label, carrierAndThreeAdsorbates(test.sites), 1, 0, 300.0, 1.0e3, 1.0e5, 2, 0,
                                 static_cast<size_t>(test.method), test.iastMethod);
    for (const std::vector<double>& fractions : sensitivityCompositions)
    {
      SCOPED_TRACE(test.label + ", adsorbate fractions " + std::to_string(fractions[1]) + ", " +
                   std::to_string(fractions[3]));
      EXPECT_LT(relativeSensitivityError(prediction, fractions, 1.0e5), 2.0e-5);
    }
  }
}

TEST(MixturePrediction, SensitivitiesReduceToFlooredPureSlopesWithoutAdsorbingGas)
{
  const Isotherm divergingSite(Isotherm::Type::Langmuir_Freundlich, {0.4, 0.0024, 0.8});
  const Isotherm langmuirSite(Isotherm::Type::Langmuir, {0.9, 0.0201});
  for (const size_t method : {size_t{0}, size_t{1}, size_t{4}, size_t{5}})
  {
    SCOPED_TRACE(method);
    MixturePrediction prediction("dilute", carrierAndThreeAdsorbates({std::vector<Isotherm>{divergingSite, langmuirSite},
                                                                      std::vector<Isotherm>{langmuirSite},
                                                                      std::vector<Isotherm>{divergingSite}}),
                                 1, 0, 300.0, 1.0e3, 1.0e5, 2, 0, method, 0);
    const std::vector<double> fractions{1.0, 0.0, 0.0, 0.0};
    std::vector<double> adsorbed(4), loadings(4), cachedPressure(8), cachedPotential(2), sensitivities(16);
    double temperature = 300.0;
    prediction.predictMixture(fractions, 1.0e5, adsorbed, loadings, cachedPressure, cachedPotential, temperature);
    prediction.predictMixtureSensitivities(fractions, 1.0e5, cachedPressure, temperature, 7.0, sensitivities);

    const double floor = 1.0e-10 * 1.0e5;
    EXPECT_DOUBLE_EQ(sensitivities[1 + 1 * 4], divergingSite.slope(floor) + langmuirSite.slope(floor));
    EXPECT_DOUBLE_EQ(sensitivities[2 + 2 * 4], langmuirSite.slope(floor));
    EXPECT_DOUBLE_EQ(sensitivities[3 + 3 * 4], divergingSite.slope(floor));
    EXPECT_TRUE(std::isfinite(sensitivities[3 + 3 * 4]));
    EXPECT_DOUBLE_EQ(sensitivities[1 + 2 * 4], 0.0);
    EXPECT_DOUBLE_EQ(sensitivities[0], 0.0);
  }
}

TEST(PhysisorptionEquilibrium, SensitivitiesFollowConcentrationsAtFrozenTotalPressure)
{
  const std::array<std::vector<Isotherm>, 3> sites{
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.4, 0.0024, 0.8}),
                            Isotherm(Isotherm::Type::Langmuir, {0.9, 0.0201})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir_Freundlich, {0.53, 0.0124, 0.822})},
      std::vector<Isotherm>{Isotherm(Isotherm::Type::Langmuir, {1.5, 5.0e-4})}};
  using Input = MixturePrediction::DrivingForceInput;
  for (const Input input : {Input::MoleFraction, Input::Concentration})
  {
    SCOPED_TRACE(input == Input::MoleFraction ? "gas" : "liquid");
    const size_t method = static_cast<size_t>(input == Input::MoleFraction ? MixturePrediction::PredictionMethod::IAST
                                                                           : MixturePrediction::PredictionMethod::SPI);
    std::vector<MixturePrediction> mixtures{
        MixturePrediction("column", carrierAndThreeAdsorbates(sites), 1, 0, 300.0, 1.0e3, 1.0e5, 2, 0, method, 1)};
    const size_t n = 4;
    const size_t terms = mixtures[0].maxIsothermTerms;
    const std::vector<double> fractionOfAdsorbent{1.0};
    const std::vector<bool> hasAdsorbentOfType{true};
    const std::vector<double> totalPressure{1.0e5};
    std::vector<double> temperature{300.0};
    const std::vector<double> pH{7.0};
    std::vector<double> idealFractions(n), adsorbed(n), molecules(n), cachedPressure(n * terms), cachedPotential(terms);
    std::vector<double> moleFraction(n), equilibrium(n), mixtureSensitivity(n * n), sensitivity(n * n);
    std::pair<size_t, size_t> performance{};

    const auto equilibriumAt = [&](const std::vector<double>& concentration)
    {
      double total = 0.0;
      for (double c : concentration) total += c;
      for (size_t comp = 0; comp < n; ++comp) moleFraction[comp] = concentration[comp] / total;
      std::fill(cachedPotential.begin(), cachedPotential.end(), 0.0);
      computePhysisorptionEquilibriumLoadings(mixtures, 0, n, 1, fractionOfAdsorbent, hasAdsorbentOfType, terms,
                                              performance, idealFractions, adsorbed, molecules, totalPressure,
                                              equilibrium, cachedPressure, cachedPotential, moleFraction, temperature,
                                              input, concentration, pH);
      return equilibrium;
    };

    const std::vector<double> concentration{37.0, 2.0, 0.5, 1.0};
    equilibriumAt(concentration);
    computePhysisorptionEquilibriumSensitivities(mixtures, 0, n, 1, fractionOfAdsorbent, hasAdsorbentOfType, terms,
                                                 idealFractions, mixtureSensitivity, totalPressure, cachedPressure,
                                                 moleFraction, temperature, input, concentration, pH, sensitivity);

    double scale = 0.0;
    double error = 0.0;
    for (size_t l = 0; l < n; ++l)
    {
      const double step = 1.0e-5 * concentration[l];
      std::vector<double> plus = concentration;
      std::vector<double> minus = concentration;
      plus[l] += step;
      minus[l] -= step;
      const std::vector<double> upper = equilibriumAt(plus);
      const std::vector<double> lower = equilibriumAt(minus);
      for (size_t i = 0; i < n; ++i)
      {
        const double reference = (upper[i] - lower[i]) / (2.0 * step);
        scale = std::max(scale, std::abs(reference));
        error = std::max(error, std::abs(sensitivity[i + l * n] - reference));
      }
    }
    EXPECT_GT(scale, 0.0);
    EXPECT_LT(error / scale, 1.0e-5);
  }
}

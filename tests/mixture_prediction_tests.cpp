#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "component.h"
#include "isotherm.h"
#include "mixture_prediction.h"
#include "multi_site_isotherm.h"

TEST(MultiSiteIsotherm, RetainsButIgnoresZeroLoadingSites)
{
  const Isotherm inactive(Isotherm::Type::Langmuir, {0.0, 0.0}, false);
  const Isotherm active(Isotherm::Type::Langmuir, {1.5, 2.0e-5}, false);
  MultiSiteIsotherm isotherm({inactive, active});

  ASSERT_EQ(isotherm.sites.size(), 2U);
  EXPECT_FALSE(isotherm.sites[0].enabled());
  EXPECT_TRUE(isotherm.sites[1].enabled());
  EXPECT_DOUBLE_EQ(isotherm.value(0, 1.0e5, 1.0), 0.0);
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

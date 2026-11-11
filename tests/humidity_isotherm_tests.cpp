#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "component.h"
#include "inputreader.h"
#include "isotherm.h"
#include "json.h"
#include "mixture_prediction.h"
#include "temporary_directory.h"
#include "utils.h"
#include "water_saturation_pressure.h"

namespace
{
double antoine(double temperature)
{
  return std::pow(10.0, 8.07131 - 1730.63 / (temperature - 39.724)) * 133.32236842105263;
}

double gabLoading(double relativeHumidity, double monolayer, double c, double k)
{
  const double x = k * relativeHumidity;
  return monolayer * c * x / ((1.0 - x) * (1.0 + x * (c - 1.0)));
}

double toth(double saturationLoading, double affinity, double exponent, double pressure)
{
  const double bp = affinity * pressure;
  return saturationLoading * bp / std::pow(1.0 + std::pow(bp, exponent), 1.0 / exponent);
}

/**
 * \brief Checks Isotherm::slope against a central difference of Isotherm::value at fixed humidity.
 */
void expectSlopeMatchesValue(const Isotherm& isotherm, double pressure, double relativeHumidity, double waterLoading)
{
  const double step = 1.0e-6 * pressure;
  const double difference = (isotherm.value(pressure + step, 7.0, relativeHumidity, waterLoading) -
                             isotherm.value(pressure - step, 7.0, relativeHumidity, waterLoading)) /
                            (2.0 * step);
  EXPECT_NEAR(isotherm.slope(pressure, 7.0, relativeHumidity, waterLoading), difference,
              1.0e-6 * std::abs(difference));
}

/**
 * \brief Liu et al. (2024) Figure 4 at 0.20 relative humidity: N2 carrier, CO2 on a GAB-modified Toth site, and H2O
 *        on its GAB site, 298.15 K and 1 bar.
 */
nlohmann::json humidInput()
{
  return {{"SimulationType", "Breakthrough"},
          {"BreakthroughIntegrator", "ROS34PW2"},
          {"MixturePredictionMethod", "SPI"},
          {"Temperature", 298.15},
          {"ColumnVoidFraction", 0.5},
          {"ParticleDensity", 1000.0},
          {"ParticleDiameter", 0.0005},
          {"DynamicViscosity", 1.78e-5},
          {"InletPressure", 1.0e5},
          {"PressureGradient", 0.0},
          {"ColumnEntranceVelocity", 0.016},
          {"ColumnLength", 0.02},
          {"BoundaryCondition", "FixedVelocity"},
          {"SimulationTime", 10.0},
          {"WriteInterval", 10.0},
          {"NumberOfGridPoints", 4},
          {"energyBalance", false},
          {"Geometry",
           {{"Type", "PackedBed"},
            {"ColumnVoidFraction", 0.5},
            {"ParticleDiameter", 0.0005},
            {"InternalDiameter", 0.01},
            {"OuterDiameter", 0.012}}},
          {"Components",
           {{{"Name", "N2"}, {"GasPhaseMolFraction", 0.9932875093666653}, {"CarrierGas", true}},
            {{"Name", "CO2"},
             {"GasPhaseMolFraction", 0.00039663312473638216},
             {"MassTransferCoefficient", 0.02},
             {"PhysisorptionSites",
              {{{"Type", "GAB-modified-Toth"}, {"Parameters", {1.0, 0.9869232667160128, 0.25, 0.006, 25.0}}}}}},
            {{"Name", "H2O"},
             {"Water", true},
             {"GasPhaseMolFraction", 0.00631585750859835},
             {"MassTransferCoefficient", 0.03},
             {"PhysisorptionSites", {{{"Type", "GAB"}, {"Parameters", {35.0, 0.15, 0.51}}}}}}}}};
}

InputReader readInput(const nlohmann::json& input)
{
  const std::filesystem::path path =
      testTemporaryDirectory() /
      (std::string{"ruptura_humidity_"} + ::testing::UnitTest::GetInstance()->current_test_info()->name() + ".json");
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
    FAIL() << "the input was accepted; expected: " << expectedMessage;
  }
  catch (const std::runtime_error& error)
  {
    EXPECT_NE(std::string{error.what()}.find(expectedMessage), std::string::npos) << error.what();
  }
}
}  // namespace

TEST(WaterSaturationPressure, WagnerReproducesTheIAPWSValues)
{
  const WaterSaturationPressure wagner{WaterSaturationPressure::Model::Wagner};
  // Triple point, IAPWS-IF97 verification values at 300, 500, and 600 K, and the normal boiling point.
  EXPECT_NEAR(wagner(273.16), 611.657, 1.0e-3);
  EXPECT_NEAR(wagner(300.0) / 3536.58941 - 1.0, 0.0, 1.0e-4);
  EXPECT_NEAR(wagner(500.0) / 2.63889776e6 - 1.0, 0.0, 2.0e-4);
  EXPECT_NEAR(wagner(600.0) / 12.3443146e6 - 1.0, 0.0, 1.0e-4);
  EXPECT_NEAR(wagner(373.124) / 101325.0 - 1.0, 0.0, 2.0e-5);
  EXPECT_THROW(wagner(647.096), std::domain_error);

  EXPECT_EQ(WaterSaturationPressure{}(298.15), antoine(298.15));
  EXPECT_EQ((WaterSaturationPressure{WaterSaturationPressure::Model::Constant, 3169.0}(350.0)), 3169.0);
}

TEST(Isotherm, GABUsesTheChosenSaturationPressure)
{
  Isotherm gab(Isotherm::Type::GAB, {35.0, 0.15, 0.51});
  gab.setTemperature(298.15);
  const double antoineLoading = gab.value(1000.0);
  EXPECT_NEAR(antoineLoading, gabLoading(1000.0 / antoine(298.15), 35.0, 0.15, 0.51), 1.0e-13);

  gab.waterSaturationPressure = {WaterSaturationPressure::Model::Wagner};
  const double wagnerPressure = WaterSaturationPressure{WaterSaturationPressure::Model::Wagner}(298.15);
  EXPECT_NEAR(gab.value(1000.0), gabLoading(1000.0 / wagnerPressure, 35.0, 0.15, 0.51), 1.0e-13);
  EXPECT_LT(gab.value(1000.0), antoineLoading);
}

TEST(Isotherm, HumidityEnhancedLangmuirSaturatesTheAffinityEnhancement)
{
  const Isotherm site(Isotherm::Type::Langmuir_humidity, {2.0, 1.0e-3, 3.0, 4.0});
  const Isotherm dry(Isotherm::Type::Langmuir, {2.0, 1.0e-3});
  EXPECT_DOUBLE_EQ(site.value(500.0, 7.0, 0.0), dry.value(500.0));

  const double affinity = 1.0e-3 * (1.0 + 3.0 * 4.0 * 0.5 / (1.0 + 4.0 * 0.5));
  EXPECT_NEAR(site.value(500.0, 7.0, 0.5), 2.0 * affinity * 500.0 / (1.0 + affinity * 500.0), 1.0e-14);
  expectSlopeMatchesValue(site, 500.0, 0.5, 0.0);
}

TEST(Isotherm, MoistureSwingLangmuirBindsMoreStronglyWhenDry)
{
  const Isotherm site(Isotherm::Type::Langmuir_moisture_swing, {1.5, 2.0e-4, 1.0});
  EXPECT_DOUBLE_EQ(site.value(800.0, 7.0, 1.0), Isotherm(Isotherm::Type::Langmuir, {1.5, 2.0e-4}).value(800.0));
  EXPECT_DOUBLE_EQ(site.value(800.0, 7.0, 0.25), Isotherm(Isotherm::Type::Langmuir, {1.5, 8.0e-4}).value(800.0));
  expectSlopeMatchesValue(site, 800.0, 0.25, 0.0);

  // RH^(-n_w) diverges in a dry gas, which saturates the site at any positive pressure.
  EXPECT_DOUBLE_EQ(site.value(800.0, 7.0, 0.0), 1.5);
  EXPECT_DOUBLE_EQ(site.value(0.0, 7.0, 0.0), 0.0);
  EXPECT_DOUBLE_EQ(site.slope(800.0, 7.0, 0.0), 0.0);
}

TEST(Isotherm, GABCoadsorptionEnhancesTheAffinityAndBlocksTheCapacity)
{
  constexpr double waterLoading = 2.0;
  const Isotherm langmuir(Isotherm::Type::GAB_Langmuir, {3.0, 1.0e-4, 0.5, 0.1});
  const double affinity = 1.0e-4 * (1.0 + 0.5 * waterLoading);
  const double capacity = 3.0 * (1.0 - 0.1 * waterLoading);
  EXPECT_NEAR(langmuir.value(400.0, 7.0, 0.3, waterLoading),
              capacity * affinity * 400.0 / (1.0 + affinity * 400.0), 1.0e-14);
  EXPECT_DOUBLE_EQ(langmuir.value(400.0, 7.0, 0.3, 20.0), 0.0);
  expectSlopeMatchesValue(langmuir, 400.0, 0.3, waterLoading);

  // Without enhancement and blocking the base isotherms remain.
  const Isotherm drySips(Isotherm::Type::GAB_Sips, {3.0, 1.0e-4, 1.4, 0.0, 0.0});
  const Isotherm dryToth(Isotherm::Type::GAB_Toth, {3.0, 1.0e-4, 0.4, 0.0, 0.0});
  EXPECT_DOUBLE_EQ(drySips.value(400.0, 7.0, 0.3, waterLoading),
                   Isotherm(Isotherm::Type::Sips, {3.0, 1.0e-4, 1.4}).value(400.0));
  EXPECT_DOUBLE_EQ(dryToth.value(400.0, 7.0, 0.3, waterLoading),
                   Isotherm(Isotherm::Type::Toth, {3.0, 1.0e-4, 0.4}).value(400.0));

  const Isotherm humidSips(Isotherm::Type::GAB_Sips, {3.0, 1.0e-4, 1.4, 0.5, 0.1});
  const Isotherm humidToth(Isotherm::Type::GAB_Toth, {3.0, 1.0e-4, 0.4, 0.5, 0.1});
  EXPECT_NEAR(humidSips.value(400.0, 7.0, 0.3, waterLoading),
              Isotherm(Isotherm::Type::Sips, {capacity, affinity, 1.4}).value(400.0), 1.0e-14);
  EXPECT_NEAR(humidToth.value(400.0, 7.0, 0.3, waterLoading), toth(capacity, affinity, 0.4, 400.0), 1.0e-14);
  expectSlopeMatchesValue(humidSips, 400.0, 0.3, waterLoading);
  expectSlopeMatchesValue(humidToth, 400.0, 0.3, waterLoading);
}

TEST(Isotherm, GABModifiedTothReproducesTheFrozenLiuFigureFourParameters)
{
  // liu/generate.py froze n_s = 1/(1 - 0.006 q_w) and b = b_0 (1 + 25 q_w) at the GAB loading of 0.20 RH.
  const Isotherm site(Isotherm::Type::GAB_modified_Toth, {1.0, 0.9869232667160128, 0.25, 0.006, 25.0});
  const double waterLoading = gabLoading(0.2, 35.0, 0.15, 0.51);
  const double frozen = toth(1.0039330155910466, 17.096831891428014, 0.25, 40.0);
  EXPECT_NEAR(site.value(40.0, 7.0, 0.2, waterLoading) / frozen - 1.0, 0.0, 1.0e-14);
  expectSlopeMatchesValue(site, 40.0, 0.2, waterLoading);
  EXPECT_THROW(site.value(40.0, 7.0, 0.2, 200.0), std::domain_error);
}

TEST(Isotherm, ExponentialHumidityTothScalesTheCapacity)
{
  const Isotherm site(Isotherm::Type::Toth_exponential_humidity, {0.81, 6.0e-2, 0.4, 1.11});
  EXPECT_DOUBLE_EQ(site.value(40.0, 7.0, 0.0), Isotherm(Isotherm::Type::Toth, {0.81, 6.0e-2, 0.4}).value(40.0));
  EXPECT_NEAR(site.value(40.0, 7.0, 0.85) / site.value(40.0, 7.0, 0.0), std::exp(1.11 * 0.85), 1.0e-14);
  expectSlopeMatchesValue(site, 40.0, 0.85, 0.0);
}

TEST(Isotherm, HumiditySitesDefineNoGrandPotential)
{
  for (const Isotherm::Type type : {Isotherm::Type::Langmuir_humidity, Isotherm::Type::GAB_modified_Toth,
                                    Isotherm::Type::Toth_exponential_humidity})
  {
    const Isotherm site(type, {1.0, 1.0e-3, 0.5, 0.1, 0.1});
    EXPECT_FALSE(site.definesGrandPotential());
    EXPECT_THROW(static_cast<void>(site.psiForPressure(100.0)), std::runtime_error);
  }
  EXPECT_FALSE(Isotherm(Isotherm::Type::GAB, {35.0, 0.15, 0.51}).definesGrandPotential());
  EXPECT_TRUE(Isotherm(Isotherm::Type::Toth, {1.0, 1.0e-3, 0.5}).definesGrandPotential());
}

TEST(MixturePrediction, SPITakesTheHumidityFromTheWaterComponent)
{
  MixturePrediction prediction(readInput(humidInput()));
  const std::vector<double> gasFractions{0.9932875093666653, 0.00039663312473638216, 0.00631585750859835};
  constexpr double pressure = 1.0e5;
  std::vector<double> adsorbedFractions(3, 0.0);
  std::vector<double> loadings(3, 0.0);
  std::vector<double> cachedPressure(3, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  double temperature = 298.15;
  prediction.predictMixture(gasFractions, pressure, adsorbedFractions, loadings, cachedPressure,
                            cachedGrandPotential, temperature);

  // RH = 0.20 at the feed, so CO2 has the frozen Figure 4 parameters, and its water loading is the H2O loading.
  const double rh = gasFractions[2] * pressure / antoine(298.15);
  EXPECT_NEAR(rh, 0.2, 1.0e-14);
  EXPECT_NEAR(loadings[2], gabLoading(rh, 35.0, 0.15, 0.51), 1.0e-13);
  EXPECT_NEAR(loadings[1] / toth(1.0039330155910466, 17.096831891428014, 0.25, gasFractions[1] * pressure) - 1.0, 0.0,
              1.0e-12);

  // A pure component carries no water: its CO2 isotherm is dry.
  std::vector<double> pureLoadings(3, 0.0);
  prediction.predictPureComponentLoadings(40.0, pureLoadings, 298.15);
  EXPECT_NEAR(pureLoadings[1], toth(1.0, 0.9869232667160128, 0.25, 40.0), 1.0e-14);
}

TEST(MixturePrediction, FixedRelativeHumidityOverridesTheWaterPartialPressure)
{
  nlohmann::json input = humidInput();
  input["RelativeHumidity"] = 0.5;
  MixturePrediction prediction(readInput(input));
  const std::vector<double> gasFractions{0.9932875093666653, 0.00039663312473638216, 0.00631585750859835};
  std::vector<double> adsorbedFractions(3, 0.0);
  std::vector<double> loadings(3, 0.0);
  std::vector<double> cachedPressure(3, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  double temperature = 298.15;
  prediction.predictMixture(gasFractions, 1.0e5, adsorbedFractions, loadings, cachedPressure, cachedGrandPotential,
                            temperature);
  EXPECT_NEAR(loadings[1] / toth(1.0139577011713956, 57.59341031907539, 0.25, gasFractions[1] * 1.0e5) - 1.0, 0.0,
              1.0e-12);
  // The water component keeps its own partial pressure.
  EXPECT_NEAR(loadings[2], gabLoading(0.2, 35.0, 0.15, 0.51), 1.0e-12);

  std::vector<double> sensitivities(9, 0.0);
  prediction.predictMixtureSensitivities(gasFractions, 1.0e5, cachedPressure, temperature, 7.0, sensitivities);
  EXPECT_EQ(sensitivities[1 + 2 * 3], 0.0);
}

TEST(MixturePrediction, SensitivityToTheWaterPartialPressureMatchesTheLoadings)
{
  MixturePrediction prediction(readInput(humidInput()));
  constexpr double pressure = 1.0e5;
  double temperature = 298.15;
  std::vector<double> adsorbedFractions(3, 0.0);
  std::vector<double> cachedPressure(3, 0.0);
  std::vector<double> cachedGrandPotential(1, 0.0);
  const auto loadingsAt = [&](double waterFraction)
  {
    const std::vector<double> gasFractions{1.0 - 0.00039663312473638216 - waterFraction, 0.00039663312473638216,
                                           waterFraction};
    std::vector<double> loadings(3, 0.0);
    prediction.predictMixture(gasFractions, pressure, adsorbedFractions, loadings, cachedPressure,
                              cachedGrandPotential, temperature);
    return loadings;
  };

  constexpr double waterFraction = 0.00631585750859835;
  constexpr double step = 1.0e-7;
  const double difference =
      (loadingsAt(waterFraction + step)[1] - loadingsAt(waterFraction - step)[1]) / (2.0 * step * pressure);

  const std::vector<double> gasFractions{1.0 - 0.00039663312473638216 - waterFraction, 0.00039663312473638216,
                                         waterFraction};
  loadingsAt(waterFraction);
  std::vector<double> sensitivities(9, 0.0);
  prediction.predictMixtureSensitivities(gasFractions, pressure, cachedPressure, temperature, 7.0, sensitivities);
  EXPECT_GT(difference, 0.0);
  EXPECT_NEAR(sensitivities[1 + 2 * 3], difference, 1.0e-5 * difference);
  EXPECT_NEAR(sensitivities[2 + 2 * 3], (loadingsAt(waterFraction + step)[2] - loadingsAt(waterFraction - step)[2]) /
                                            (2.0 * step * pressure),
              1.0e-5 * sensitivities[2 + 2 * 3]);
}

TEST(HumidityInput, ReadsTheWaterComponentAndTheSaturationPressure)
{
  nlohmann::json input = humidInput();
  input["WaterSaturationPressureModel"] = "Wagner";
  const InputReader reader = readInput(input);
  EXPECT_TRUE(reader.components[2].isWater);
  EXPECT_FALSE(reader.components[1].isWater);
  EXPECT_FALSE(reader.relativeHumidity.has_value());
  EXPECT_EQ(reader.waterSaturationPressure.model, WaterSaturationPressure::Model::Wagner);
  EXPECT_EQ(reader.components[2].isotherm.sites[0].waterSaturationPressure.model,
            WaterSaturationPressure::Model::Wagner);

  input.erase("WaterSaturationPressureModel");
  input["WaterSaturationPressure"] = 3169.0;
  const InputReader constant = readInput(input);
  EXPECT_EQ(constant.waterSaturationPressure.model, WaterSaturationPressure::Model::Constant);
  EXPECT_EQ(constant.waterSaturationPressure.constantPressure, 3169.0);
}

TEST(HumidityInput, RejectsInconsistentHumidityInput)
{
  nlohmann::json twoWaters = humidInput();
  twoWaters["Components"][1]["Water"] = true;
  expectRejected(twoWaters, "at most one component can be marked Water");

  nlohmann::json carrierWater = humidInput();
  carrierWater["Components"][0]["Water"] = true;
  carrierWater["Components"][2]["Water"] = false;
  expectRejected(carrierWater, "cannot be the CarrierGas");

  nlohmann::json noWater = humidInput();
  noWater["Components"][2]["Water"] = false;
  expectRejected(noWater, "needs a component with \"Water\": true or a RelativeHumidity");
  noWater["RelativeHumidity"] = 0.2;
  expectRejected(noWater, "take the water loading from a component with \"Water\": true");

  nlohmann::json iast = humidInput();
  iast["MixturePredictionMethod"] = "IAST";
  iast["Components"][2]["PhysisorptionSites"] = {{{"Type", "Langmuir"}, {"Parameters", {1.0, 1.0e-3}}}};
  expectRejected(iast, "humidity-dependent isotherms require MixturePredictionMethod SPI");

  nlohmann::json humidWater = humidInput();
  humidWater["Components"][2]["PhysisorptionSites"] = {
      {{"Type", "Humidity-enhanced-Langmuir"}, {"Parameters", {1.0, 1.0e-3, 1.0, 1.0}}}};
  expectRejected(humidWater, "cannot have a humidity-dependent isotherm");

  nlohmann::json constant = humidInput();
  constant["WaterSaturationPressureModel"] = "Constant";
  expectRejected(constant, "needs a positive WaterSaturationPressure");
  constant["WaterSaturationPressureModel"] = "Wagner";
  constant["WaterSaturationPressure"] = 3169.0;
  expectRejected(constant, "WaterSaturationPressure is the pressure of the Constant");

  nlohmann::json negative = humidInput();
  negative["RelativeHumidity"] = -0.1;
  expectRejected(negative, "RelativeHumidity must be a finite non-negative number");
}

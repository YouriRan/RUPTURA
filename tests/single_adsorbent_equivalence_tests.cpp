#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <span>
#include <string>

#include "column.h"
#include "inputreader.h"
#include "integrators/compute.h"
#include "integrators/rk3.h"
#include "mixture_prediction.h"
#include "utils.h"

namespace
{
struct Setting
{
  bool energyBalance;
  size_t gridPoints;
  MixturePrediction::PredictionMethod method;
};

constexpr std::array settings{
    Setting{false, 50, MixturePrediction::PredictionMethod::IAST},
    Setting{false, 50, MixturePrediction::PredictionMethod::SIAST},
    Setting{false, 100, MixturePrediction::PredictionMethod::IAST},
    Setting{false, 100, MixturePrediction::PredictionMethod::SIAST},
    Setting{true, 50, MixturePrediction::PredictionMethod::IAST},
    Setting{true, 50, MixturePrediction::PredictionMethod::SIAST},
    Setting{true, 100, MixturePrediction::PredictionMethod::IAST},
    Setting{true, 100, MixturePrediction::PredictionMethod::SIAST},
};

InputReader makeReader(const Setting& setting)
{
  const auto path = std::filesystem::path{RUPTURA_SOURCE_DIR} /
                    "examples/MOR-CO2-C3H8/breakthrough/simulation.json";
  InputReader reader(path.string());
  reader.energyBalance = setting.energyBalance;
  reader.numberOfGridPoints = setting.gridPoints;
  reader.columnDistances = makeUniformColumnDistances(reader.numberOfGridPoints, reader.columnLength);
  reader.mixturePredictionMethod = static_cast<size_t>(setting.method);
  reader.internalDiameter = 1.0e-2;
  reader.outerDiameter = 1.2e-2;
  reader.wallDensity = 7800.0;
  reader.gasThermalConductivity = 0.025;
  reader.wallThermalConductivity = 15.0;
  reader.heatTransferGasSolid = 20.0;
  reader.heatTransferGasWall = 10.0;
  reader.heatTransferWallExternal = 5.0;
  reader.heatCapacityGas = 1000.0;
  reader.heatCapacitySolid = 900.0;
  reader.heatCapacityWall = 500.0;
  reader.geometry = makeGeometry(PackedBedTubeSpec{.voidFraction = reader.columnVoidFraction,
                                                   .particleDiameter = reader.particleDiameter,
                                                   .internalDiameter = reader.internalDiameter,
                                                   .outerDiameter = reader.outerDiameter});
  return reader;
}

InputReader makeTwoIdenticalAdsorbentReader(const Setting& setting)
{
  InputReader reader = makeReader(setting);
  if (reader.adsorbentComponents.size() != 1U)
  {
    throw std::runtime_error("equivalence fixture requires one source adsorbent");
  }
  reader.adsorbentComponents.push_back(reader.adsorbentComponents.front());
  reader.adsorbentLengths = {0.5 * reader.columnLength, 0.5 * reader.columnLength};
  reader.adsorbentInterfaceLengths = {0.0};
  const size_t firstGridPoints = reader.numberOfGridPoints / 2;
  reader.adsorbentGridPoints = {firstGridPoints, reader.numberOfGridPoints - firstGridPoints};
  reader.adsorbentInterfaceGridPoints = {0};
  reader.adsorbentVoidFractions = {reader.columnVoidFraction, reader.columnVoidFraction};
  reader.adsorbentParticleDensities = {reader.particleDensity, reader.particleDensity};
  reader.adsorbentParticleDiameters = {reader.particleDiameter, reader.particleDiameter};
  reader.adsorbentGeometries = {reader.geometry, reader.geometry};
  return reader;
}

template <typename Left, typename Right>
void expectNumericallyEqual(const Left& left, const Right& right, const char* quantity,
                            double relativeTolerance = 1.0e-12, double absoluteTolerance = 1.0e-12)
{
  SCOPED_TRACE(quantity);
  const std::span<const double> leftValues{left};
  const std::span<const double> rightValues{right};
  ASSERT_EQ(leftValues.size(), rightValues.size());
  // Relative to the largest entry of the field: near-zero entries (e.g. derivatives) result from cancelling
  // terms of the field's magnitude, so their rounding error under -ffast-math scales with that magnitude.
  double scale = 0.0;
  for (size_t index = 0; index < leftValues.size(); ++index)
  {
    scale = std::max({scale, std::abs(leftValues[index]), std::abs(rightValues[index])});
  }
  for (size_t index = 0; index < leftValues.size(); ++index)
  {
    ASSERT_NEAR(leftValues[index], rightValues[index], absoluteTolerance + relativeTolerance * scale)
        << "index " << index;
  }
}

void synchronizeState(const Column& single, Column& multi)
{
  ASSERT_EQ(single.state.size(), multi.state.size());
  std::ranges::copy(single.state, multi.state.begin());
  std::ranges::copy(single.totalPressure, multi.totalPressure.begin());
  std::ranges::copy(single.totalConcentration, multi.totalConcentration.begin());
  std::ranges::copy(single.moleFraction, multi.moleFraction.begin());
  std::ranges::copy(single.partialPressure, multi.partialPressure.begin());
  std::ranges::copy(single.interstitialGasVelocity, multi.interstitialGasVelocity.begin());
  std::ranges::copy(single.gasDensity, multi.gasDensity.begin());
}

class SingleAdsorbentEquivalence : public ::testing::TestWithParam<size_t>
{
 protected:
  InputReader reader() const { return makeReader(settings.at(GetParam())); }
};

TEST_P(SingleAdsorbentEquivalence, TwoIdenticalAdsorbentsMatchOne)
{
  InputReader singleInput = reader();
  InputReader multiInput = makeTwoIdenticalAdsorbentReader(settings.at(GetParam()));
  Column single(singleInput);
  Column multi(multiInput);
  ASSERT_EQ(multi.numberOfAdsorbents, 2U);
  single.initialize();
  multi.initialize();
  synchronizeState(single, multi);

  computeDerivatives(single);
  computeDerivatives(multi);

  expectNumericallyEqual(single.equilibriumPhysisorption, multi.equilibriumPhysisorption,
                         "equilibrium physisorption");
  expectNumericallyEqual(single.equilibriumChemisorption, multi.equilibriumChemisorption,
                         "equilibrium chemisorption");
  expectNumericallyEqual(single.physisorptionDot, multi.physisorptionDot, "physisorption derivatives");
  expectNumericallyEqual(single.chemisorptionDot, multi.chemisorptionDot, "chemisorption derivatives");
  expectNumericallyEqual(single.bulkSpeciesSink, multi.bulkSpeciesSink, "bulk species sink");
  expectNumericallyEqual(single.interstitialGasVelocity, multi.interstitialGasVelocity, "interstitial velocity");
  expectNumericallyEqual(single.totalPressure, multi.totalPressure, "total pressure");
  expectNumericallyEqual(single.gasDensity, multi.gasDensity, "gas density");
  expectNumericallyEqual(single.partialPressure, multi.partialPressure, "partial pressure");
  expectNumericallyEqual(single.concentrationDot, multi.concentrationDot, "mass derivatives");
  expectNumericallyEqual(single.gasTemperatureDot, multi.gasTemperatureDot, "gas energy derivatives");
  expectNumericallyEqual(single.solidTemperatureDot, multi.solidTemperatureDot, "solid energy derivatives");
  expectNumericallyEqual(single.wallTemperatureDot, multi.wallTemperatureDot, "wall energy derivatives");
}

INSTANTIATE_TEST_SUITE_P(MultibedSettings, SingleAdsorbentEquivalence,
                         ::testing::Range<size_t>(0, settings.size()));
}  // namespace

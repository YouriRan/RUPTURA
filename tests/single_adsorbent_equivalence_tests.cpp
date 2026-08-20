#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <span>
#include <string>

#include "column.h"
#include "column_multibed.h"
#include "inputreader.h"
#include "integrators/compute.h"
#include "integrators/compute_multibed.h"
#include "integrators/rk3.h"
#include "integrators/sorption.h"
#include "integrators/transport.h"
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
  reader.adsorbentVoidFractions = {reader.columnVoidFraction, reader.columnVoidFraction};
  reader.adsorbentParticleDensities = {reader.particleDensity, reader.particleDensity};
  reader.adsorbentParticleDiameters = {reader.particleDiameter, reader.particleDiameter};
  reader.adsorbentGeometries = {reader.geometry, reader.geometry};
  return reader;
}

template <typename Left, typename Right>
void expectIdentical(const Left& left, const Right& right, const char* quantity)
{
  SCOPED_TRACE(quantity);
  const std::span<const double> leftValues{left};
  const std::span<const double> rightValues{right};
  ASSERT_EQ(leftValues.size(), rightValues.size());
  for (size_t index = 0; index < leftValues.size(); ++index)
  {
    ASSERT_DOUBLE_EQ(leftValues[index], rightValues[index]) << "index " << index;
  }
}

template <typename Left, typename Right>
void expectNumericallyEqual(const Left& left, const Right& right, const char* quantity,
                            double relativeTolerance = 1.0e-12, double absoluteTolerance = 1.0e-12)
{
  SCOPED_TRACE(quantity);
  const std::span<const double> leftValues{left};
  const std::span<const double> rightValues{right};
  ASSERT_EQ(leftValues.size(), rightValues.size());
  for (size_t index = 0; index < leftValues.size(); ++index)
  {
    const double scale = std::max(std::abs(leftValues[index]), std::abs(rightValues[index]));
    ASSERT_NEAR(leftValues[index], rightValues[index], absoluteTolerance + relativeTolerance * scale)
        << "index " << index;
  }
}

void synchronizeState(const Column& single, MultibedColumn& multi)
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

class SingleAdsorbentComputeEquivalence : public ::testing::TestWithParam<size_t>
{
 protected:
  InputReader reader() const { return makeReader(settings.at(GetParam())); }
};

TEST_P(SingleAdsorbentComputeEquivalence, PhysisorptionEquilibriumAndKineticsAreIdentical)
{
  InputReader input = reader();
  Column single(input);
  MultibedColumn multi(input);
  single.initialize();
  multi.initialize();
  synchronizeState(single, multi);

  computePhysisorptionEquilibriumLoadings(
      single.physisorptionMixture, single.numberOfGridPoints, single.numberOfComponents, single.maxIsothermTerms,
      single.iastPerformance, single.idealGasMolFractions, single.adsorbedMolFractions, single.numberOfMolecules,
      single.totalPressure, single.equilibriumPhysisorption, single.cachedPressure, single.cachedGrandPotential,
      single.moleFraction, single.gasTemperature, MixturePrediction::DrivingForceInput::MoleFraction,
      single.concentration, single.pH);
  computePhysisorptionEquilibriumLoadings(
      multi.physisorptionMixtures, multi.numberOfGridPoints, multi.numberOfComponents, multi.numberOfAdsorbents,
      multi.fractionOfAdsorbent, multi.hasAdsorbentOfType, multi.maxIsothermTerms, multi.iastPerformance,
      multi.idealGasMolFractions, multi.adsorbedMolFractions, multi.numberOfMolecules, multi.totalPressure,
      multi.equilibriumPhysisorption, multi.cachedPressure, multi.cachedGrandPotential, multi.moleFraction,
      multi.gasTemperature, MixturePrediction::DrivingForceInput::MoleFraction, multi.concentration, multi.pH);
  expectIdentical(single.equilibriumPhysisorption, multi.equilibriumPhysisorption, "equilibrium physisorption");

  computePhysisorption(single.components, single.numberOfGridPoints, single.numberOfComponents,
                       single.equilibriumPhysisorption, single.physisorption, single.physisorptionDot);
  computePhysisorption(multi.physisorptionMixtures, multi.numberOfGridPoints, multi.numberOfComponents,
                       multi.numberOfAdsorbents, multi.fractionOfAdsorbent, multi.equilibriumPhysisorption,
                       multi.physisorption, multi.physisorptionDot);
  expectIdentical(single.physisorptionDot, multi.physisorptionDot, "physisorption derivatives");
}

TEST_P(SingleAdsorbentComputeEquivalence, ChemisorptionRoutinesAreIdentical)
{
  InputReader input = reader();
  Chemisorption site;
  site.type = Chemisorption::Type::General;
  site.maximumLoading = 2.0;
  site.adsorptionRateCoefficient = 1.0e-4;
  input.components[1].chemisorption.add(site);
  input.adsorbentComponents[0][1].chemisorption.add(site);
  Column single(input);
  MultibedColumn multi(input);
  single.initialize();
  multi.initialize();
  synchronizeState(single, multi);
  computeDerivatives(single);
  computeDerivatives(multi);
  expectIdentical(single.equilibriumChemisorption, multi.equilibriumChemisorption, "equilibrium chemisorption");
  expectIdentical(single.chemisorptionDot, multi.chemisorptionDot, "chemisorption derivatives");
}

TEST_P(SingleAdsorbentComputeEquivalence, ReactionRoutinesAreIdentical)
{
  InputReader input = reader();
  Reaction reaction;
  reaction.phase = Reaction::Phase::Physisorbed;
  reaction.style = Reaction::Style::GeneralPowerLaw;
  reaction.reactants = {1};
  reaction.products = {2};
  reaction.reactantStoichiometry = {1.0};
  reaction.productStoichiometry = {1.0};
  reaction.forwardOrders = {1.0};
  reaction.backwardOrders = {1.0};
  reaction.forwardRateCoefficient = 1.0;
  reaction.equilibriumConstant = 1.0e9;
  input.reactions = {reaction};
  Column single(input);
  MultibedColumn multi(input);
  single.initialize();
  multi.initialize();
  single.physisorption[1] = 0.4;
  synchronizeState(single, multi);
  computeDerivatives(single);
  computeDerivatives(multi);
  expectIdentical(single.reactionPhysisorptionSource, multi.reactionPhysisorptionSource, "reaction source");
  expectIdentical(single.reactionHeat, multi.reactionHeat, "reaction heat");
}

TEST_P(SingleAdsorbentComputeEquivalence, TransportRoutinesAreIdentical)
{
  InputReader input = reader();
  Column single(input);
  MultibedColumn multi(input);
  single.initialize();
  multi.initialize();
  synchronizeState(single, multi);
  std::fill(single.bulkSpeciesSink.begin(), single.bulkSpeciesSink.end(), 0.0);
  std::fill(multi.bulkSpeciesSink.begin(), multi.bulkSpeciesSink.end(), 0.0);
  updateVelocityAndPressure(
      single.components, single.boundaryCondition, single.numberOfGridPoints, single.numberOfComponents,
      single.inletPressure, single.outletPressure, single.pressureGradient, single.columnLength, single.geometry,
      single.columnEntranceVelocity, single.dynamicViscosity, single.resolution, single.interstitialGasVelocity,
      single.gasDensity, single.totalConcentration, single.totalPressure, single.concentration, single.partialPressure,
      single.moleFraction, single.bulkSpeciesSink, single.gasTemperature, single.fluidPhase, single.liquidDensity,
      single.pHMode, single.pHValue, single.pKw, single.pHComponent, single.pH);
  updateVelocityAndPressure(
      multi.components, multi.boundaryCondition, multi.numberOfGridPoints, multi.numberOfComponents,
      multi.inletPressure, multi.outletPressure, multi.pressureGradient, multi.columnLength, multi.numberOfAdsorbents,
      multi.columnEntranceVelocity, multi.dynamicViscosity, multi.columnDistances, multi.fractionOfAdsorbent,
      multi.geometries, multi.adsorbentScaledVoidFraction, multi.totalVoidFraction,
      multi.interstitialGasVelocity, multi.gasDensity,
      multi.totalConcentration,
      multi.totalPressure, multi.concentration, multi.partialPressure, multi.moleFraction, multi.bulkSpeciesSink,
      multi.gasTemperature, multi.fluidPhase, multi.liquidDensity, multi.pHMode, multi.pHValue, multi.pKw,
      multi.pHComponent, multi.pH);
  EXPECT_DOUBLE_EQ(single.columnEntranceVelocity, multi.columnEntranceVelocity);
  expectIdentical(single.interstitialGasVelocity, multi.interstitialGasVelocity, "interstitial velocity");
  expectIdentical(single.totalPressure, multi.totalPressure, "pressure transport");
  expectIdentical(single.gasDensity, multi.gasDensity, "gas density");
  expectIdentical(single.partialPressure, multi.partialPressure, "partial pressure");
}

TEST_P(SingleAdsorbentComputeEquivalence, MassBalanceRoutinesAreIdentical)
{
  InputReader input = reader();
  Column single(input);
  MultibedColumn multi(input);
  single.initialize();
  multi.initialize();
  synchronizeState(single, multi);
  computeDerivatives(single);
  computeDerivatives(multi);
  expectIdentical(single.bulkSpeciesSink, multi.bulkSpeciesSink, "bulk species sink");
  expectNumericallyEqual(single.concentrationDot, multi.concentrationDot, "mass derivatives");
}

TEST_P(SingleAdsorbentComputeEquivalence, EnergyBalanceRoutinesAreIdentical)
{
  InputReader input = reader();
  input.energyBalance = true;
  Column single(input);
  MultibedColumn multi(input);
  single.initialize();
  multi.initialize();
  synchronizeState(single, multi);
  computeDerivatives(single);
  computeDerivatives(multi);
  expectNumericallyEqual(single.gasTemperatureDot, multi.gasTemperatureDot, "gas energy derivatives");
  expectNumericallyEqual(single.solidTemperatureDot, multi.solidTemperatureDot, "solid energy derivatives");
  expectNumericallyEqual(single.wallTemperatureDot, multi.wallTemperatureDot, "wall energy derivatives");
}

TEST_P(SingleAdsorbentComputeEquivalence, TwoIdenticalAdsorbentsMatchSingleBedCalculations)
{
  InputReader singleInput = reader();
  InputReader multiInput = makeTwoIdenticalAdsorbentReader(settings.at(GetParam()));
  Column single(singleInput);
  MultibedColumn multi(multiInput);
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

INSTANTIATE_TEST_SUITE_P(MultibedSettings, SingleAdsorbentComputeEquivalence,
                         ::testing::Range<size_t>(0, settings.size()));
}  // namespace

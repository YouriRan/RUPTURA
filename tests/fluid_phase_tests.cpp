#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "column.h"
#include "column_multibed.h"
#include "inputreader.h"
#include "integrators/cvode.h"
#include "integrators/rk3.h"
#include "isotherm.h"
#include "json.h"
#include "mixture_prediction.h"
#include "timing.h"
#include "temporary_directory.h"

namespace
{
struct TransportedPHInput
{
  std::filesystem::path path =
      testTemporaryDirectory() /
      (std::string{"ruptura_"} + ::testing::UnitTest::GetInstance()->current_test_info()->name() + ".json");
  nlohmann::json data;

  TransportedPHInput()
  {
    const auto source = std::filesystem::path{RUPTURA_SOURCE_DIR} /
                        "examples/Liquid/multibed-pH-front/simulation.json";
    std::ifstream input(source);
    data = nlohmann::json::parse(input);
  }

  void write() const
  {
    std::ofstream output(path);
    output << data;
  }

  ~TransportedPHInput()
  {
    std::error_code ignored;
    std::filesystem::remove(path, ignored);
  }
};
}  // namespace

TEST(LiquidPhaseInput, TransportedPHUsesZeroBasedComponentIdentifier)
{
  TransportedPHInput input;
  std::swap(input.data["Components"][0], input.data["Components"][1]);
  input.data["Components"][1]["Name"] = "transported acid";
  for (auto& adsorbent : input.data["Adsorbents"])
  {
    auto& parameters = adsorbent["ComponentParameters"];
    parameters["transported acid"] = parameters["H+"];
    parameters.erase("H+");
  }
  input.data["pHComponent"] = 1;
  input.write();

  InputReader reader(input.path.string());
  EXPECT_EQ(reader.pHComponent, 1U);
  EXPECT_EQ(reader.components[reader.pHComponent].name, "transported acid");
  MultibedColumn column(reader);
  column.initialize();
  EXPECT_EQ(column.pHComponent, 1U);
  EXPECT_NEAR(column.pH[0], 8.0, 1.0e-14);
  EXPECT_NEAR(column.pH[1], 7.0, 1.0e-14);

  input.data["pHMode"] = "OHMinus";
  input.data["pKw"] = 13.5;
  input.data.erase("pHComponent");
  input.data["phcomponent"] = 1;
  input.write();
  InputReader hydroxideReader(input.path.string());
  Column hydroxideColumn(hydroxideReader);
  hydroxideColumn.initialize();
  EXPECT_EQ(hydroxideColumn.pHComponent, 1U);
  EXPECT_NEAR(hydroxideColumn.pH[0], 5.5, 1.0e-14);
  EXPECT_NEAR(hydroxideColumn.pH[1], 6.5, 1.0e-14);
}

TEST(LiquidPhaseInput, TransportedPHRejectsInvalidComponentIdentifiers)
{
  TransportedPHInput input;
  for (const auto& identifier : nlohmann::json::array({"H+", "0", -1, 0.5, 0.0, true, nullptr, 2,
                                                     18446744073709551615ULL}))
  {
    SCOPED_TRACE(identifier.dump());
    input.data["pHComponent"] = identifier;
    input.write();
    EXPECT_THROW(InputReader reader(input.path.string()), std::runtime_error);
  }
  input.data.erase("pHComponent");
  input.write();
  EXPECT_THROW(InputReader reader(input.path.string()), std::runtime_error);

  input.data["pHMode"] = "Fixed";
  input.write();
  EXPECT_NO_THROW(InputReader reader(input.path.string()));
}

TEST(GasMixturePredictionRegression, MoleFractionsRemainTheDrivingForceForSPI)
{
  Component first(0, "first");
  first.isotherm.add(Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}, false));
  Component second(1, "second");
  second.isotherm.add(Isotherm(Isotherm::Type::Langmuir, {4.0, 2.0e-5}, false));

  MixturePrediction prediction("gas SPI regression", {first, second}, 0, 0, 298.15, 1.0e3, 1.0e5, 2, 1,
                               static_cast<size_t>(MixturePrediction::PredictionMethod::SPI),
                               static_cast<size_t>(MixturePrediction::IASTMethod::FastIAST));
  std::vector<double> gasMoleFractions{0.25, 0.75};
  std::vector<double> adsorbedMoleFractions(2);
  std::vector<double> loadings(2);
  std::vector<double> cachedPressure(2);
  std::vector<double> cachedGrandPotential(1);
  double temperature = 298.15;

  prediction.predictMixture(gasMoleFractions, 1.0e5, adsorbedMoleFractions, loadings, cachedPressure,
                            cachedGrandPotential, temperature);

  EXPECT_NEAR(loadings[0], 0.4, 1.0e-14);
  EXPECT_NEAR(loadings[1], 2.4, 1.0e-14);
  EXPECT_NEAR(adsorbedMoleFractions[0], 1.0 / 7.0, 1.0e-14);
  EXPECT_NEAR(adsorbedMoleFractions[1], 6.0 / 7.0, 1.0e-14);
}

TEST(GasMixturePredictionRegression, MoleFractionsRemainTheDrivingForceForSCI)
{
  Component first(0, "first");
  first.isotherm.add(Isotherm(Isotherm::Type::Langmuir, {2.0, 1.0e-5}, false));
  Component second(1, "second");
  second.isotherm.add(Isotherm(Isotherm::Type::Langmuir, {2.0, 2.0e-5}, false));

  MixturePrediction prediction("gas SCI regression", {first, second}, 0, 0, 298.15, 1.0e3, 1.0e5, 2, 1,
                               static_cast<size_t>(MixturePrediction::PredictionMethod::SCI),
                               static_cast<size_t>(MixturePrediction::IASTMethod::FastIAST));
  std::vector<double> gasMoleFractions{0.25, 0.75};
  std::vector<double> adsorbedMoleFractions(2);
  std::vector<double> loadings(2);
  std::vector<double> cachedPressure(2);
  std::vector<double> cachedGrandPotential(1);
  double temperature = 298.15;

  prediction.predictMixture(gasMoleFractions, 1.0e5, adsorbedMoleFractions, loadings, cachedPressure,
                            cachedGrandPotential, temperature);

  EXPECT_NEAR(loadings[0], 2.0 / 11.0, 1.0e-14);
  EXPECT_NEAR(loadings[1], 12.0 / 11.0, 1.0e-14);
}

TEST(GasBreakthroughRegression, ExistingColumnCvodePathRetainsItsOutcome)
{
  const std::filesystem::path path = std::filesystem::path{RUPTURA_SOURCE_DIR} /
                                     "examples/CoBDP-alkanes-C6/breakthrough/simulation.json";
  InputReader reader(path.string());
  Column column(reader);
  column.initialize();
  ASSERT_EQ(column.fluidPhase, Column::FluidPhase::Gas);

  CVODE integrator(reader.timeStep, false, 1);
  // The tolerances of the reference values (the CVODE defaults when they were recorded), so that later changes of
  // the defaults do not move them.
#if BUILD_SUNDIALS
  integrator.relativeTolerance = 1.0e-3;
  integrator.absoluteToleranceConcentration = 1.0e-6;
  integrator.absoluteToleranceLoading = 1.0e-6;
  integrator.absoluteToleranceTemperature = 1.0e-6;
#endif
  Timing timing;
  integrator.initialize(column);
  ASSERT_TRUE(integrator.propagate(column, 0, timing));

  const size_t index = column.numberOfComponents + 1;
  EXPECT_NEAR(column.concentration[index], 1.245870563589704, 1.0e-9);
  EXPECT_NEAR(column.moleFraction[index], 0.0022944618465073376, 1.0e-12);
}

TEST(GasBreakthroughRegression, ExistingMultibedRk3PathRetainsItsOutcome)
{
  const std::filesystem::path path =
      std::filesystem::path{RUPTURA_SOURCE_DIR} / "examples/Multibed/CO2-N2-polishing/simulation.json";
  InputReader reader(path.string());
  MultibedColumn column(reader);
  column.initialize();
  ASSERT_EQ(column.fluidPhase, MultibedColumn::FluidPhase::Gas);

  RungeKutta3 integrator(reader.timeStep, false, 1);
  Timing timing;
  ASSERT_TRUE(integrator.propagate(column, 0, timing));

  const size_t index = column.numberOfComponents + 1;
  EXPECT_NEAR(column.concentration[index], 0.29182638429337593, 1.0e-13);
  EXPECT_NEAR(column.moleFraction[index], 0.0072778674138634803, 1.0e-15);
}

TEST(LiquidPhase, FixedPHColumnUsesConcentrationsAndExpectedIsotherm)
{
  const std::filesystem::path path =
      std::filesystem::path{RUPTURA_SOURCE_DIR} / "examples/Liquid/pH-dependent/simulation.json";
  InputReader reader(path.string());
  Column column(reader);
  column.initialize();

  ASSERT_EQ(column.fluidPhase, Column::FluidPhase::Liquid);
  EXPECT_DOUBLE_EQ(column.concentration[0], 2.0);
  EXPECT_DOUBLE_EQ(column.concentration[1], 0.0);
  EXPECT_DOUBLE_EQ(column.gasDensity[0], 997.0);
  EXPECT_DOUBLE_EQ(column.pH[0], 6.5);
  const double pHFactor = std::pow(10.0, 0.5);
  const double activity = 0.5 * 2.0 / (1.0 + pHFactor);
  EXPECT_NEAR(column.equilibriumPhysisorption[0], 2.0 * activity / (1.0 + activity), 1.0e-14);
  EXPECT_DOUBLE_EQ(column.equilibriumPhysisorption[1], 0.0);
  EXPECT_NEAR(column.totalPressure[1],
              150000.0 - column.geometry.pressureDrop.gradient(0.00089, 997.0, 0.001) * column.resolution,
              1.0e-10);

  RungeKutta3 integrator(reader.timeStep, false, 1);
  Timing timing;
  ASSERT_TRUE(integrator.propagate(column, 0, timing));
  EXPECT_NEAR(column.concentration[1], 0.012766632220756373, 1.0e-14);
  EXPECT_NEAR(column.physisorption[0], 0.0015466082723401815, 1.0e-15);
}

TEST(LiquidPhase, TransportedPHAndSuperficialFluxWorkInMultibedColumn)
{
  const std::filesystem::path path =
      std::filesystem::path{RUPTURA_SOURCE_DIR} / "examples/Liquid/multibed-pH-front/simulation.json";
  InputReader reader(path.string());
  MultibedColumn column(reader);
  column.initialize();

  ASSERT_EQ(column.fluidPhase, MultibedColumn::FluidPhase::Liquid);
  ASSERT_EQ(column.pHMode, MultibedColumn::PHMode::HPlus);
  EXPECT_EQ(reader.pHComponent, 0U);
  EXPECT_EQ(column.pHComponent, 0U);
  EXPECT_NEAR(column.pH[0], 8.0, 1.0e-14);
  EXPECT_NEAR(column.pH[1], 7.0, 1.0e-14);
  EXPECT_NEAR(column.equilibriumPhysisorption[1], 0.625, 1.0e-14);
  for (size_t grid = 1; grid < column.numberOfGridPoints + 1; ++grid)
  {
    EXPECT_NEAR(column.totalVoidFraction[grid] * column.interstitialGasVelocity[grid],
                column.totalVoidFraction[0] * column.interstitialGasVelocity[0], 1.0e-15);
  }

  CVODE integrator(reader.timeStep, false, 1);
  // The tolerances of the reference value (the CVODE defaults when it was recorded).
#if BUILD_SUNDIALS
  integrator.relativeTolerance = 1.0e-3;
  integrator.absoluteToleranceConcentration = 1.0e-6;
  integrator.absoluteToleranceLoading = 1.0e-6;
  integrator.absoluteToleranceTemperature = 1.0e-6;
#endif
  Timing timing;
  integrator.initialize(column);
  ASSERT_TRUE(integrator.propagate(column, 0, timing));
  const size_t soluteAtFirstInteriorNode = column.numberOfComponents + 1;
  EXPECT_NEAR(column.concentration[soluteAtFirstInteriorNode], 0.010039044104904167, 1.0e-10);
  EXPECT_NEAR(column.pH[1], -std::log10(column.concentration[column.numberOfComponents] / 1000.0), 1.0e-14);
}

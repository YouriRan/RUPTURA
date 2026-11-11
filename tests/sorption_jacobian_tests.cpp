#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "column.h"
#include "inputreader.h"
#include "integrators/block_banded.h"
#include "integrators/compute.h"
#include "integrators/cvode.h"
#include "integrators/jacobian.h"
#include "integrators/rk3.h"
#include "integrators/sorption.h"
#include "integrators/transport.h"
#include "json.h"
#include "temporary_directory.h"
#include "timing.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

// The example on 12 grid cells with the given convection scheme, with every component given the axial dispersion
// coefficient and, if given, the boundary condition.
InputReader smallInput(const std::string& example, const std::string& name, double dispersion,
                       const std::string& boundaryCondition = "", const std::string& convectionScheme = "Upwind")
{
  std::ifstream input(examplePath(example));
  nlohmann::json data = nlohmann::json::parse(input);
  data["ConvectionScheme"] = convectionScheme;
  data["NumberOfGridPoints"] = 12;
  data["BreakthroughIntegrator"] = "RungeKutta3";
  if (!boundaryCondition.empty()) data["BoundaryCondition"] = boundaryCondition;
  for (auto& component : data["Components"]) component["AxialDispersionCoefficient"] = dispersion;
  const std::filesystem::path path = testTemporaryDirectory() / ("ruptura_sorption_jacobian_" + name + ".json");
  std::ofstream(path) << data;
  return InputReader(path.string());
}

// A front entering the column: the feed decays as exp(-0.6 i) along the nodes from the feed end (z = L in a reverse
// step), a carrier fills up the total concentration of the initial column, and every loading is 60% of its
// equilibrium value.
void setFrontState(Column& column)
{
  const size_t n = column.numberOfComponents;
  const bool reverse = column.flowDirection == Column::FlowDirection::Reverse;
  const std::vector<double> feed =
      reverse ? column.productEndConcentration
              : std::vector<double>(column.concentration.begin(),
                                    column.concentration.begin() + static_cast<std::ptrdiff_t>(n));
  double feedTotal = 0.0;
  for (double value : feed) feedTotal += value;
  const bool hasCarrier = std::any_of(column.components.begin(), column.components.end(),
                                      [](const Component& component) { return component.isCarrierGas; });
  for (size_t node = 1; node <= column.numberOfGridPoints; ++node)
  {
    const double front =
        std::exp(-0.6 * static_cast<double>(reverse ? column.numberOfGridPoints + 1 - node : node));
    double total = 0.0;
    for (size_t comp = 0; comp < n; ++comp) total += column.concentration[node * n + comp];
    double adsorbing = 0.0;
    for (size_t comp = 0; comp < n; ++comp)
    {
      if (column.components[comp].isCarrierGas) continue;
      column.concentration[node * n + comp] = (hasCarrier ? total / feedTotal : 1.0) * feed[comp] * front;
      adsorbing += column.concentration[node * n + comp];
    }
    for (size_t comp = 0; comp < n; ++comp)
    {
      if (column.components[comp].isCarrierGas) column.concentration[node * n + comp] = total - adsorbing;
    }
  }
  // In a reverse step node 0 takes the state of the first cell, as at the start of a reverse swing phase.
  if (reverse) std::copy_n(column.concentration.begin() + static_cast<std::ptrdiff_t>(n), n, column.concentration.begin());
  Timing timings;
  precompute(column, timings);
  for (size_t index = 0; index < column.physisorption.size(); ++index)
  {
    column.physisorption[index] = 0.6 * column.equilibriumPhysisorption[index];
  }
  precompute(column, timings);
}

// The concentration and loading right-hand side with the velocity and total pressure of the current state held
// fixed: the production kernels, without updateVelocityAndPressure.
std::vector<double> frozenClosureRightHandSide(Column& column, const SorptionJacobian& jacobian)
{
  const size_t n = column.numberOfComponents;
  for (size_t node = 0; node <= column.numberOfGridPoints; ++node)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < n; ++comp) sum += std::max(0.0, column.concentration[node * n + comp]);
    for (size_t comp = 0; comp < n; ++comp)
    {
      column.moleFraction[node * n + comp] = std::max(0.0, column.concentration[node * n + comp]) / sum;
    }
  }
  computePhysisorptionEquilibriumLoadings(
      column.physisorptionMixtures, column.numberOfGridPoints, n, column.numberOfAdsorbents, column.fractionOfAdsorbent,
      column.hasAdsorbentOfType, column.maxIsothermTerms, column.iastPerformance, column.idealGasMolFractions,
      column.adsorbedMolFractions, column.numberOfMolecules, column.totalPressure, column.equilibriumPhysisorption,
      column.cachedPressure, column.cachedGrandPotential, column.moleFraction, column.gasTemperature,
      MixturePrediction::DrivingForceInput::MoleFraction, column.concentration, column.pH);
  computePhysisorption(column.numberOfGridPoints, n, column.massTransferCoefficient, column.equilibriumPhysisorption,
                       column.physisorption, column.physisorptionDot);
  computeBulkSpeciesSink(column.physisorptionMixtures, column.numberOfGridPoints, n, column.numberOfAdsorbents,
                         column.maxChemisorptionSites, column.geometries, column.fractionOfAdsorbent,
                         column.loadingPrefactor, column.concentration, column.physisorptionDot,
                         column.chemisorptionDot, column.surfaceConcentration, column.bulkSpeciesSink,
                         column.reactionPhysisorptionSource, column.reactionChemisorptionSource,
                         column.reactionBulkSource);
  computeMassDerivatives(column.numberOfGridPoints, n, column.columnDistances, column.axialDispersionCoefficient,
                         column.interstitialGasVelocity, column.concentration, column.concentrationDot,
                         column.bulkSpeciesSink, column.fluidPhase, column.totalVoidFraction, column.convectionScheme,
                         column.inletCondition, column.wenoScratch, column.flowEnds());
  std::vector<double> rightHandSide(jacobian.size);
  jacobian.gather(column.concentrationDot, column.physisorptionDot, rightHandSide);
  return rightHandSide;
}

// Largest difference of the assembled Jacobian from central differences of the right-hand side, relative to the
// largest central-difference entry of its row. Entries outside the band must vanish.
double relativeJacobianError(Column& column, const std::function<std::vector<double>(Column&)>& rightHandSide)
{
  SorptionJacobian jacobian(column);
  assembleSorptionJacobian(column, column.concentration, column.gasTemperature, jacobian);

  std::vector<double> state(jacobian.size);
  jacobian.gather(column.concentration, column.physisorption, state);
  double concentrationScale = 0.0;
  double loadingScale = 0.0;
  for (size_t node = 0; node < jacobian.numberOfNodes; ++node)
  {
    for (size_t comp = 0; comp < jacobian.numberOfComponents; ++comp)
    {
      concentrationScale = std::max(concentrationScale, state[jacobian.concentrationIndex(node, comp)]);
      loadingScale = std::max(loadingScale, state[jacobian.loadingIndex(node, comp)]);
    }
  }

  std::vector<double> reference(jacobian.size * jacobian.size);
  for (size_t col = 0; col < jacobian.size; ++col)
  {
    const bool isLoading = (col / jacobian.numberOfComponents) % 2 == 1;
    // A zero state (the carrier loading) gets a step of 1e-6 of its scale, well above the noise of the warm-started
    // mixture solver.
    const double step =
        state[col] != 0.0 ? 1.0e-5 * std::abs(state[col]) : 1.0e-6 * (isLoading ? loadingScale : concentrationScale);
    std::vector<double> perturbed = state;
    perturbed[col] = state[col] + step;
    jacobian.scatter(perturbed, column.concentration, column.physisorption);
    const std::vector<double> upper = rightHandSide(column);
    perturbed[col] = state[col] - step;
    jacobian.scatter(perturbed, column.concentration, column.physisorption);
    const std::vector<double> lower = rightHandSide(column);
    for (size_t row = 0; row < jacobian.size; ++row)
    {
      reference[row + col * jacobian.size] = (upper[row] - lower[row]) / (2.0 * step);
    }
  }
  jacobian.scatter(state, column.concentration, column.physisorption);

  double error = 0.0;
  for (size_t row = 0; row < jacobian.size; ++row)
  {
    double rowScale = 0.0;
    for (size_t col = 0; col < jacobian.size; ++col)
    {
      rowScale = std::max(rowScale, std::abs(reference[row + col * jacobian.size]));
    }
    if (rowScale == 0.0) continue;
    for (size_t col = 0; col < jacobian.size; ++col)
    {
      const bool inBand = (row > col ? row - col : col - row) <= jacobian.bandwidth;
      const double assembled = inBand ? jacobian(row, col) : 0.0;
      error = std::max(error, std::abs(assembled - reference[row + col * jacobian.size]) / rowScale);
    }
  }
  return error;
}
}  // namespace

TEST(SorptionJacobian, MatchesCentralDifferencesOfTheGasRightHandSideAtFrozenClosure)
{
  // BEA C7: eight components with IAST on two-site Langmuir isotherms, and a helium carrier.
  InputReader reader = smallInput("examples/BEA-alkanes-C7/breakthrough/simulation.json", "bea", 1.0e-5);
  Column column(reader);
  column.initialize();
  setFrontState(column);

  SorptionJacobian jacobian(column);
  EXPECT_EQ(jacobian.size, 2U * 8U * 13U);
  EXPECT_EQ(jacobian.bandwidth, 16U);
  EXPECT_LT(relativeJacobianError(column, [&](Column& state) { return frozenClosureRightHandSide(state, jacobian); }),
            1.0e-5);
}

TEST(SorptionJacobian, MatchesCentralDifferencesWithHumidityDependentCO2)
{
  // Liu et al. (2024) Figure 4 at 0.20 relative humidity: the GAB-modified Toth affinity of CO2 rises 17-fold with the
  // local water loading, so dq*_CO2/dc_H2O is the largest off-diagonal entry.
  const nlohmann::json data = {
      {"SimulationType", "Breakthrough"},
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
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_sorption_jacobian_humid_source.json";
  std::ofstream(path) << data;

  InputReader reader = smallInput(path.string(), "humid", 1.0e-5);
  Column column(reader);
  column.initialize();
  setFrontState(column);

  SorptionJacobian jacobian(column);
  EXPECT_LT(relativeJacobianError(column, [&](Column& state) { return frozenClosureRightHandSide(state, jacobian); }),
            1.0e-5);
}

TEST(SorptionJacobian, MatchesCentralDifferencesOfTheLiquidRightHandSide)
{
  // The liquid velocity does not depend on the concentrations, so the production right-hand side is the reference.
  InputReader reader = smallInput("examples/Liquid/pH-dependent/simulation.json", "liquid", 1.0e-5);
  Column column(reader);
  column.initialize();
  setFrontState(column);

  SorptionJacobian jacobian(column);
  Timing timings;
  const auto productionRightHandSide = [&](Column& state)
  {
    precompute(state, timings);
    computeDerivatives(state);
    std::vector<double> rightHandSide(jacobian.size);
    jacobian.gather(state.concentrationDot, state.physisorptionDot, rightHandSide);
    return rightHandSide;
  };
  EXPECT_LT(relativeJacobianError(column, productionRightHandSide), 1.0e-6);
}

TEST(SorptionJacobian, MatchesCentralDifferencesOfTheWeno3GasRightHandSideAtFrozenClosure)
{
  // WENO3 reconstructs the mole fractions and carries them with the upwind total molar flux, so every component of a
  // node couples to every component of the stencil, two nodes upstream to one downstream.
  InputReader reader =
      smallInput("examples/BEA-alkanes-C7/breakthrough/simulation.json", "bea_weno", 1.0e-5, "", "WENO3");
  Column column(reader);
  column.initialize();
  setFrontState(column);

  SorptionJacobian jacobian(column);
  EXPECT_EQ(jacobian.bandwidth, 40U);
  EXPECT_LT(relativeJacobianError(column, [&](Column& state) { return frozenClosureRightHandSide(state, jacobian); }),
            1.0e-5);
}

TEST(SorptionJacobian, MatchesCentralDifferencesWithADanckwertsInlet)
{
  // The inlet-face value depends on the concentrations of the first cell, with either convection scheme.
  for (const std::string scheme : {"Upwind", "WENO3"})
  {
    SCOPED_TRACE(scheme);
    std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
    nlohmann::json data = nlohmann::json::parse(input);
    data["ConvectionScheme"] = scheme;
    data["InletCondition"] = "Danckwerts";
    data["NumberOfGridPoints"] = 12;
    data["BreakthroughIntegrator"] = "RungeKutta3";
    for (auto& component : data["Components"]) component["AxialDispersionCoefficient"] = 2.0e-4;
    const std::filesystem::path path = testTemporaryDirectory() / ("ruptura_sorption_jacobian_danckwerts.json");
    std::ofstream(path) << data;
    InputReader reader(path.string());
    Column column(reader);
    column.initialize();
    setFrontState(column);

    SorptionJacobian jacobian(column);
    EXPECT_LT(
        relativeJacobianError(column, [&](Column& state) { return frozenClosureRightHandSide(state, jacobian); }),
        1.0e-5);
  }
}

TEST(SorptionJacobian, MatchesCentralDifferencesInReverseFlowAndWithBackFlowAtTheProductEnd)
{
  // A reverse step feeds at z = L, so that every stencil is mirrored; a forward step whose velocity turns negative
  // near the outlet takes gas in from the product end there, that of the last cell or a held product-end gas.
  enum class Flow
  {
    Reverse,
    BackFlowFromLastCell,
    BackFlowFromProductEndGas
  };
  for (const std::string scheme : {"Upwind", "WENO3"})
  {
    for (const std::string inlet : {"Danckwerts", "Dirichlet"})
    {
      for (const Flow flow : {Flow::Reverse, Flow::BackFlowFromLastCell, Flow::BackFlowFromProductEndGas})
      {
        SCOPED_TRACE(scheme + " " + inlet + " " + std::to_string(static_cast<int>(flow)));
        std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
        nlohmann::json data = nlohmann::json::parse(input);
        data["ConvectionScheme"] = scheme;
        data["InletCondition"] = inlet;
        data["NumberOfGridPoints"] = 12;
        data["BreakthroughIntegrator"] = "RungeKutta3";
        data["FlowDirection"] = flow == Flow::Reverse ? "Reverse" : "Forward";
        for (auto& component : data["Components"]) component["AxialDispersionCoefficient"] = 2.0e-4;
        const std::filesystem::path path = testTemporaryDirectory() / "ruptura_sorption_jacobian_reverse.json";
        std::ofstream(path) << data;
        InputReader reader(path.string());
        Column column(reader);
        column.initialize();
        setFrontState(column);
        if (flow != Flow::Reverse)
        {
          // The last three faces carry gas back into the column.
          const size_t outlet = column.numberOfGridPoints;
          for (size_t face = outlet - 2; face <= outlet; ++face) column.interstitialGasVelocity[face] *= -0.7;
        }
        if (flow == Flow::BackFlowFromProductEndGas)
        {
          const size_t n = column.numberOfComponents;
          column.productEndConcentration.assign(column.concentration.begin(),
                                                column.concentration.begin() + static_cast<std::ptrdiff_t>(n));
          for (double& value : column.productEndConcentration) value *= 0.8;
          column.productEndTemperature = column.gasTemperature[0];
        }

        SorptionJacobian jacobian(column);
        EXPECT_LT(
            relativeJacobianError(column, [&](Column& state) { return frozenClosureRightHandSide(state, jacobian); }),
            1.0e-5);
      }
    }
  }
}

TEST(SorptionJacobian, MatchesCentralDifferencesInClosedEndPhases)
{
  // A closed end in both directions, for a pressurization (gas enters at the feed end, whose concentrations follow the
  // ramp) and a blowdown (gas leaves there).
  for (const std::string scheme : {"Upwind", "WENO3"})
  {
    for (const std::string direction : {"Forward", "Reverse"})
    {
      for (const double target : {2.0e6, 1.0e5})
      {
        SCOPED_TRACE(scheme + " " + direction + " " + std::to_string(target));
        std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
        nlohmann::json data = nlohmann::json::parse(input);
        data["ConvectionScheme"] = scheme;
        data["InletCondition"] = "Danckwerts";
        data["NumberOfGridPoints"] = 12;
        data["BreakthroughIntegrator"] = "RungeKutta3";
        data["BoundaryCondition"] = "InletPressureInletVelocity";
        data["FlowDirection"] = direction;
        for (auto& component : data["Components"]) component["AxialDispersionCoefficient"] = 2.0e-4;
        const std::filesystem::path path = testTemporaryDirectory() / "ruptura_sorption_jacobian_closed.json";
        std::ofstream(path) << data;
        InputReader reader(path.string());
        Column column(reader);
        column.initialize();
        const double start = column.inletPressure;
        column.closedEnd = true;
        column.pressureRamp = {Column::PressureRamp::Shape::Exponential, 5.0, 0.0, start, target, 0.0};
        column.advancePressureRamp(0.0);
        setFrontState(column);
        const size_t closedFace = direction == "Forward" ? column.numberOfGridPoints : 0;
        ASSERT_EQ(column.interstitialGasVelocity[closedFace], 0.0);

        SorptionJacobian jacobian(column);
        EXPECT_LT(
            relativeJacobianError(column, [&](Column& state) { return frozenClosureRightHandSide(state, jacobian); }),
            1.0e-5);
      }
    }
  }
}

TEST(SorptionJacobian, MatchesCentralDifferencesOfTheWeno3LiquidRightHandSide)
{
  InputReader reader = smallInput("examples/Liquid/pH-dependent/simulation.json", "liquid_weno", 1.0e-5, "", "WENO3");
  Column column(reader);
  column.initialize();
  setFrontState(column);

  SorptionJacobian jacobian(column);
  Timing timings;
  const auto productionRightHandSide = [&](Column& state)
  {
    precompute(state, timings);
    computeDerivatives(state);
    std::vector<double> rightHandSide(jacobian.size);
    jacobian.gather(state.concentrationDot, state.physisorptionDot, rightHandSide);
    return rightHandSide;
  };
  EXPECT_LT(relativeJacobianError(column, productionRightHandSide), 1.0e-6);
}

TEST(ThermalJacobian, MatchesCentralDifferencesOfTheEnergyBalanceAtFrozenFlowAndUptake)
{
  // Experiment T15 with a temperature profile; the gas density, velocity, conductivity, and uptake rates are held at
  // their values of the state, as in the Jacobian. Also in a reverse step, in a forward step with gas flowing back in
  // at the product end, with either inlet condition, and with a closed end in a blowdown in both directions (3, 4).
  for (const std::string inlet : {"Danckwerts", "Dirichlet"})
  {
    for (const int flow : {0, 1, 2, 3, 4})
    {
      if (flow > 2 && inlet == "Dirichlet") continue;
      SCOPED_TRACE(inlet + " " + std::to_string(flow));
      std::ifstream input(examplePath("marx2/ternary/T15/simulation.json"));
      nlohmann::json data = nlohmann::json::parse(input);
      data["InletCondition"] = inlet;
      data["FlowDirection"] = flow == 1 || flow == 4 ? "Reverse" : "Forward";
      const std::filesystem::path source = testTemporaryDirectory() / "ruptura_thermal_jacobian_source.json";
      std::ofstream(source) << data;
      InputReader reader = smallInput(source.string(), "t15_thermal", 4.3e-6, "", "WENO3");
      Column column(reader);
      column.initialize();
      if (flow > 2)
      {
        column.closedEnd = true;
        column.pressureRamp = {Column::PressureRamp::Shape::Exponential, 5.0, 0.0, column.inletPressure, 1.0e5, 0.0};
        column.advancePressureRamp(0.0);
        if (flow == 3) column.productEndTemperature = 0.0;
      }
      setFrontState(column);
      for (size_t node = 0; node <= column.numberOfGridPoints; ++node)
      {
        column.gasTemperature[node] += 8.0 * std::exp(-0.4 * static_cast<double>(node));
        column.solidTemperature[node] += 5.0 * std::exp(-0.3 * static_cast<double>(node));
        column.wallTemperature[node] += 2.0 * std::exp(-0.2 * static_cast<double>(node));
      }
      Timing timings;
      precompute(column, timings);
      computeDerivatives(column);
      if (flow == 2)
      {
        // The last three faces carry gas back into the column from a held product-end gas.
        for (size_t face = column.numberOfGridPoints - 2; face <= column.numberOfGridPoints; ++face)
        {
          column.interstitialGasVelocity[face] *= -0.7;
        }
        column.productEndTemperature = 330.0;
      }

      const size_t nodes = column.numberOfGridPoints + 1;
      BlockBandedMatrix jacobian(nodes, 3, 1);
      assembleThermalJacobian(column, jacobian);

      const auto energy = [&]()
      {
        computeEnergyDerivatives(
            column.numberOfGridPoints, column.numberOfComponents, column.externalTemperature, column.geometries,
            column.particleDensity, column.solidFluidContactArea, column.fluidSolidContactArea,
            column.heatOfAdsorption, column.heatOfChemisorption, column.wallDensity, column.gasThermalConductivity,
            column.wallThermalConductivity, column.heatTransferGasSolid, column.heatTransferGasWall,
            column.heatTransferWallExternal, column.heatCapacityGas, column.heatCapacitySolid,
            column.heatCapacityWall, column.columnDistances, column.interstitialGasVelocity, column.totalVoidFraction,
            column.gasDensity,
            column.coeffDiffusion, column.maxChemisorptionSites, column.physisorptionDot, column.chemisorptionDot,
            column.gasTemperature, column.gasTemperatureDot, column.solidTemperature, column.solidTemperatureDot,
            column.wallTemperature, column.wallTemperatureDot, column.reactionPhysisorptionSource,
            column.reactionChemisorptionSource, column.reactionHeat, column.inletCondition, column.flowEnds());
        std::vector<double> rates(3 * nodes);
        for (size_t node = 0; node < nodes; ++node)
        {
          rates[3 * node] = column.gasTemperatureDot[node];
          rates[3 * node + 1] = column.solidTemperatureDot[node];
          rates[3 * node + 2] = column.wallTemperatureDot[node];
        }
        return rates;
      };
      const std::array<std::span<double>, 3> temperatures{column.gasTemperature, column.solidTemperature,
                                                          column.wallTemperature};
      double error = 0.0;
      double scale = 0.0;
      for (size_t node = 0; node < nodes; ++node)
      {
        for (size_t field = 0; field < 3; ++field)
        {
          double& value = temperatures[field][node];
          const double original = value;
          value = original + 1.0e-4;
          const std::vector<double> upper = energy();
          value = original - 1.0e-4;
          const std::vector<double> lower = energy();
          value = original;
          for (size_t row = 0; row < 3 * nodes; ++row)
          {
            const double reference = (upper[row] - lower[row]) / 2.0e-4;
            const size_t rowNode = row / 3;
            const bool inBand = (rowNode > node ? rowNode - node : node - rowNode) <= 1;
            const double assembled = inBand ? jacobian(rowNode, node, row % 3, field) : 0.0;
            scale = std::max(scale, std::abs(reference));
            error = std::max(error, std::abs(assembled - reference));
          }
        }
      }
      EXPECT_GT(scale, 0.0);
      EXPECT_LT(error / scale, 1.0e-7);
    }
  }
}

#if BUILD_SUNDIALS
TEST(CvodeKlu, MatchesSpgmrThroughAnUpwindFrontWithFewJacobians)
{
  // With upwind convection the sorption Jacobian is CVODE's Jacobian up to the Ergun density term, so the sparse
  // direct solve follows the matrix-free one.
  InputReader reader = smallInput("examples/BEA-alkanes-C7/breakthrough/simulation.json", "bea_klu", 1.0e-5);
  const auto run = [&](CVODE::LinearSolverType solver, long& jacobians)
  {
    Column column(reader);
    column.initialize();
    CVODE integrator;
    integrator.relativeTolerance = 1.0e-7;
    integrator.absoluteToleranceConcentration = 1.0e-9;
    integrator.absoluteToleranceLoading = 1.0e-10;
    integrator.maximumOrder = 5;
    integrator.linearSolverType = solver;
    Timing timing;
    integrator.initialize(column);
    integrator.advance(column, 40.0, timing);
    CVodeGetNumJacEvals(integrator.cvodeMem, &jacobians);
    return std::vector<double>(column.concentration.begin(), column.concentration.end());
  };
  long kluJacobians = 0;
  long spgmrJacobians = 0;
  const std::vector<double> klu = run(CVODE::LinearSolverType::KLU, kluJacobians);
  const std::vector<double> spgmr = run(CVODE::LinearSolverType::SPGMR, spgmrJacobians);

  double feed = 0.0;
  double difference = 0.0;
  for (size_t index = 0; index < klu.size(); ++index)
  {
    feed = std::max(feed, spgmr[index]);
    difference = std::max(difference, std::abs(klu[index] - spgmr[index]));
  }
  EXPECT_GT(kluJacobians, 0);
  EXPECT_LT(difference / feed, 1.0e-5);
}
#endif

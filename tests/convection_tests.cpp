#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "column.h"
#include "component.h"
#include "inputreader.h"
#include "integrators/compute.h"
#include "integrators/jacobian.h"
#include "integrators/rk3.h"
#include "integrators/ros34pw2.h"
#include "json.h"
#include "swing_adsorption.h"
#include "temporary_directory.h"
#include "timing.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

// values[0] at the inlet face z_0, values[i] the average of f over the cell (z_{i-1}, z_i], for f(z) = a + b z.
std::vector<double> linearCellValues(const std::vector<double>& distances, double a, double b)
{
  std::vector<double> values(distances.size());
  values[0] = a + b * distances[0];
  for (size_t i = 1; i < distances.size(); ++i) values[i] = a + b * 0.5 * (distances[i - 1] + distances[i]);
  return values;
}

std::vector<double> uniformDistances(size_t intervals, double length)
{
  std::vector<double> distances(intervals + 1);
  for (size_t i = 0; i <= intervals; ++i) distances[i] = length * static_cast<double>(i) / static_cast<double>(intervals);
  return distances;
}

// Largest face error of the uniform WENO3 reconstruction of exp(z) on [0, 1] over the faces z_0 ... z_{N-1}.
double exponentialFaceError(size_t intervals)
{
  const std::vector<double> distances = uniformDistances(intervals, 1.0);
  const double h = distances[1];
  std::vector<double> values(intervals + 1);
  values[0] = 1.0;
  for (size_t i = 1; i <= intervals; ++i) values[i] = (std::exp(distances[i]) - std::exp(distances[i - 1])) / h;
  std::vector<double> faces(intervals + 1);
  computeWENO(values, distances, faces);
  double error = 0.0;
  for (size_t i = 0; i < intervals; ++i) error = std::max(error, std::abs(faces[i] - std::exp(distances[i])));
  return error;
}

nlohmann::json beaInput()
{
  std::ifstream input(examplePath("examples/BEA-alkanes-C7/breakthrough/simulation.json"));
  return nlohmann::json::parse(input);
}
}  // namespace

TEST(WENO3Reconstruction, IsExactForLinearProfilesOnUniformAndNonUniformGrids)
{
  const std::vector<double> uniform = uniformDistances(8, 0.4);
  const std::vector<double> stretched{0.0, 0.01, 0.025, 0.05, 0.06, 0.1, 0.16, 0.2, 0.21};
  for (const std::vector<double>& distances : {uniform, stretched})
  {
    const std::vector<double> values = linearCellValues(distances, 0.3, -1.7);
    std::vector<double> faces(values.size());
    computeWENO(values, distances, faces);
    // The outlet face is first order (last cell value); every other face is exact.
    for (size_t i = 0; i + 1 < values.size(); ++i) EXPECT_NEAR(faces[i], 0.3 - 1.7 * distances[i], 1.0e-14) << i;
  }
}

TEST(WENO3Reconstruction, IsThirdOrderForASmoothProfile)
{
  const double coarse = exponentialFaceError(20);
  const double medium = exponentialFaceError(40);
  const double fine = exponentialFaceError(80);
  EXPECT_GT(coarse / medium, 6.0);
  EXPECT_GT(medium / fine, 6.0);
}

TEST(WENO3Reconstruction, KeepsAStepWithinItsBounds)
{
  std::vector<double> values(21, 0.0);
  std::fill(values.begin(), values.begin() + 8, 1.0);
  std::vector<double> faces(values.size());
  computeWENO(values, uniformDistances(20, 1.0), faces);
  for (double face : faces)
  {
    EXPECT_GE(face, -1.0e-12);
    EXPECT_LE(face, 1.0 + 1.0e-12);
  }
  EXPECT_NEAR(faces[7], 1.0, 1.0e-12);  // last cell before the step
  EXPECT_NEAR(faces[8], 0.0, 1.0e-12);  // first cell after the step
}

TEST(MassBalanceWENO3, GasSpeciesFluxesSumToTheUpwindTotalFlux)
{
  // The velocity closures conserve the upwind total molar flux v_i c_tot,i, so the WENO3 species fluxes must sum to
  // it: sum_j dc_j/dt = -(v_i c_tot,i - v_{i-1} c_tot,i-1) / dz without dispersion and sinks.
  constexpr size_t gridPoints = 12;
  constexpr size_t componentCount = 3;
  constexpr double resolution = 0.01;
  std::vector<double> velocity(gridPoints + 1);
  std::vector<double> concentration((gridPoints + 1) * componentCount);
  for (size_t grid = 0; grid <= gridPoints; ++grid)
  {
    const double z = static_cast<double>(grid) / static_cast<double>(gridPoints);
    velocity[grid] = 0.1 - 0.03 * z;
    const double front = 0.5 + 0.5 * std::tanh(15.0 * (0.4 - z));
    concentration[grid * componentCount + 0] = 30.0 + 2.0 * z;
    concentration[grid * componentCount + 1] = 5.0 * front;
    concentration[grid * componentCount + 2] = 3.0 * (1.0 - front) + 0.2;
  }
  std::vector<double> concentrationDot(concentration.size());
  const std::vector<double> sink(concentration.size(), 0.0);
  const std::vector<double> ones(gridPoints + 1, 1.0);
  const std::vector<double> noDispersion(concentration.size(), 0.0);
  std::vector<double> wenoScratch((gridPoints + 1) * (componentCount + 4));
  computeMassDerivatives(gridPoints, componentCount, uniformDistances(gridPoints, gridPoints * resolution),
                         noDispersion, velocity, concentration, concentrationDot, sink, Column::FluidPhase::Gas, ones,
                         Column::ConvectionScheme::WENO3, Column::InletCondition::Danckwerts, wenoScratch);

  auto total = [&](size_t grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < componentCount; ++comp) sum += concentration[grid * componentCount + comp];
    return sum;
  };
  for (size_t grid = 1; grid <= gridPoints; ++grid)
  {
    double sumDot = 0.0;
    for (size_t comp = 0; comp < componentCount; ++comp) sumDot += concentrationDot[grid * componentCount + comp];
    const double upwindTotal = -(velocity[grid] * total(grid) - velocity[grid - 1] * total(grid - 1)) / resolution;
    EXPECT_NEAR(sumDot, upwindTotal, 1.0e-10 * std::abs(velocity[0] * total(0) / resolution)) << grid;
  }
}

TEST(MassBalanceWENO3, IsExactForLinearProfilesInGasAndLiquid)
{
  // Linear mole fractions at constant c_tot and v (gas), and linear concentrations at constant v (liquid): the
  // convective term is -v dc/dz at every node except the first-order outlet node.
  constexpr size_t gridPoints = 10;
  constexpr double resolution = 0.02;
  constexpr double v = 0.05;
  const std::vector<double> distances = uniformDistances(gridPoints, gridPoints * resolution);
  const std::vector<double> fraction = linearCellValues(distances, 0.8, -2.0);
  std::vector<double> concentration(2 * (gridPoints + 1));
  for (size_t grid = 0; grid <= gridPoints; ++grid)
  {
    concentration[2 * grid] = 40.0 * fraction[grid];
    concentration[2 * grid + 1] = 40.0 * (1.0 - fraction[grid]);
  }
  const std::vector<double> velocity(gridPoints + 1, v);
  const std::vector<double> sink(concentration.size(), 0.0);
  const std::vector<double> ones(gridPoints + 1, 1.0);
  const std::vector<double> noDispersion(concentration.size(), 0.0);
  std::vector<double> wenoScratch((gridPoints + 1) * (2 + 4));
  for (Column::FluidPhase phase : {Column::FluidPhase::Gas, Column::FluidPhase::Liquid})
  {
    std::vector<double> concentrationDot(concentration.size());
    computeMassDerivatives(gridPoints, 2, distances, noDispersion, velocity, concentration, concentrationDot, sink,
                           phase, ones, Column::ConvectionScheme::WENO3, Column::InletCondition::Danckwerts,
                           wenoScratch);
    for (size_t grid = 1; grid < gridPoints; ++grid)
    {
      EXPECT_NEAR(concentrationDot[2 * grid], -v * 40.0 * -2.0, 1.0e-10) << grid;
      EXPECT_NEAR(concentrationDot[2 * grid + 1], -v * 40.0 * 2.0, 1.0e-10) << grid;
    }
  }
}

TEST(ConvectionSchemeInput, DefaultsToWENO3ForEveryIntegrator)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_convection_scheme.json";
  nlohmann::json data = beaInput();
  std::ofstream(path) << data;
  {
    const InputReader reader(path.string());
    EXPECT_EQ(reader.convectionScheme, 1U);
    EXPECT_EQ(Column(reader).convectionScheme, Column::ConvectionScheme::WENO3);
  }

  data["ConvectionScheme"] = "Upwind";
  std::ofstream(path) << data;
  EXPECT_EQ(Column(InputReader(path.string())).convectionScheme, Column::ConvectionScheme::Upwind);

  data["ConvectionScheme"] = "Central";
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  data.erase("ConvectionScheme");
  data["BreakthroughIntegrator"] = "ROS34PW2";
  std::ofstream(path) << data;
  EXPECT_EQ(InputReader(path.string()).convectionScheme, 1U);
}

TEST(InletConditionInput, DefaultsToDanckwertsAndAcceptsDanckwertsForEveryIntegrator)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_inlet_condition.json";
  nlohmann::json data = beaInput();
  std::ofstream(path) << data;
  {
    const InputReader reader(path.string());
    EXPECT_EQ(reader.inletCondition, 1U);
    EXPECT_EQ(Column(reader).inletCondition, Column::InletCondition::Danckwerts);
  }

  data["InletCondition"] = "Dirichlet";
  std::ofstream(path) << data;
  EXPECT_EQ(Column(InputReader(path.string())).inletCondition, Column::InletCondition::Dirichlet);

  data["InletCondition"] = "Robin";
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  data["InletCondition"] = "Danckwerts";
  data["BreakthroughIntegrator"] = "ROS34PW2";
  std::ofstream(path) << data;
  EXPECT_EQ(InputReader(path.string()).inletCondition, 1U);
}

TEST(DanckwertsInlet, AdmitsExactlyTheFeedFluxWithBothSchemesInGasAndLiquid)
{
  // On a uniform grid with uniform dispersion, the dispersive face fluxes telescope, and sum_i h dc_i/dt is the inflow
  // across the inlet face minus the outflow v c_N. The Danckwerts inlet admits the feed flux v c_in; the Dirichlet
  // inlet adds the dispersive flux D (c_in - c_1)/delta, with delta = h/2 (WENO3) or h (upwind).
  constexpr size_t gridPoints = 16;
  constexpr size_t componentCount = 2;
  constexpr double length = 0.32;
  constexpr double v = 0.05;
  constexpr double dispersion = 2.0e-3;  // v h / D = 0.5
  const double h = length / static_cast<double>(gridPoints);
  const std::vector<double> distances = uniformDistances(gridPoints, length);
  const std::vector<double> feed{4.0, 32.0};
  std::vector<double> concentration((gridPoints + 1) * componentCount);
  concentration[0] = feed[0];
  concentration[1] = feed[1];
  for (size_t grid = 1; grid <= gridPoints; ++grid)
  {
    const double z = 0.5 * (distances[grid - 1] + distances[grid]) / length;
    concentration[grid * componentCount] = 3.0 * (0.5 + 0.5 * std::tanh(12.0 * (0.4 - z))) + 0.1;
    concentration[grid * componentCount + 1] = 36.0 - concentration[grid * componentCount];
  }
  const std::vector<double> velocity(gridPoints + 1, v);
  const std::vector<double> sink(concentration.size(), 0.0);
  const std::vector<double> ones(gridPoints + 1, 1.0);
  const std::vector<double> axialDispersion(concentration.size(), dispersion);
  std::vector<double> wenoScratch((gridPoints + 1) * (componentCount + 4));

  for (Column::ConvectionScheme scheme : {Column::ConvectionScheme::WENO3, Column::ConvectionScheme::Upwind})
  {
    const double delta = scheme == Column::ConvectionScheme::WENO3 ? 0.5 * h : h;
    for (Column::FluidPhase phase : {Column::FluidPhase::Gas, Column::FluidPhase::Liquid})
    {
      for (Column::InletCondition inlet : {Column::InletCondition::Danckwerts, Column::InletCondition::Dirichlet})
      {
        std::vector<double> concentrationDot(concentration.size());
        computeMassDerivatives(gridPoints, componentCount, distances, axialDispersion, velocity, concentration,
                               concentrationDot, sink, phase, ones, scheme, inlet, wenoScratch);
        for (size_t comp = 0; comp < componentCount; ++comp)
        {
          double inflow = v * concentration[gridPoints * componentCount + comp];
          for (size_t grid = 1; grid <= gridPoints; ++grid) inflow += h * concentrationDot[grid * componentCount + comp];
          const double dispersiveInflow =
              inlet == Column::InletCondition::Dirichlet
                  ? dispersion * (feed[comp] - concentration[componentCount + comp]) / delta
                  : 0.0;
          EXPECT_NEAR(inflow, v * feed[comp] + dispersiveInflow, 1.0e-12 * v * 36.0)
              << static_cast<int>(scheme) << " " << static_cast<int>(phase) << " " << static_cast<int>(inlet) << " "
              << comp;
        }
      }
    }
  }
}

TEST(DanckwertsInlet, ConductsNoHeatAcrossTheInletFace)
{
  // Gas temperature on a uniform grid without heat exchange: sum_i h dT_i/dt is v (T_in - T_N) with a Danckwerts
  // inlet, plus the conduction a (T_in - T_1)/h, a = lambda / (rho C_p), with a Dirichlet inlet.
  const Geometry geometry = makeGeometry(PackedBedTubeSpec{0.4, 3.0e-3, 0.025, 0.04});
  constexpr size_t gridPoints = 10;
  constexpr size_t nodes = gridPoints + 1;
  constexpr double v = 0.02;
  constexpr double gasDensityValue = 1.2;
  constexpr double heatCapacityGas = 1000.0;
  constexpr double gasThermalConductivity = 12.0;  // a = 1e-2 m^2/s
  const double a = gasThermalConductivity / (gasDensityValue * heatCapacityGas);
  const std::vector<double> distances = uniformDistances(gridPoints, 0.2);
  const double h = distances[1];
  std::vector<double> gasTemperature(nodes);
  gasTemperature[0] = 350.0;
  for (size_t grid = 1; grid < nodes; ++grid) gasTemperature[grid] = 300.0 + 40.0 * std::exp(-4.0 * distances[grid]);
  const std::vector<double> velocity(nodes, v);
  const std::vector<double> gasDensity(nodes, gasDensityValue);
  const std::vector<double> zeros(nodes, 0.0);
  const std::vector<double> particleDensity(nodes, 1000.0);
  const std::vector<double> wallTemperature(nodes, 300.0);
  const std::vector<Geometry> geometries{geometry};

  for (Column::InletCondition inlet : {Column::InletCondition::Danckwerts, Column::InletCondition::Dirichlet})
  {
    std::vector<double> coeffDiffusion(nodes, 0.0);
    std::vector<double> gasTemperatureDot(nodes, 0.0);
    std::vector<double> solidTemperatureDot(nodes, 0.0);
    std::vector<double> wallTemperatureDot(nodes, 0.0);
    computeEnergyDerivatives(gridPoints, 1, 300.0, geometries, particleDensity, zeros, zeros, zeros, {}, 4000.0,
                             gasThermalConductivity, 0.0, 0.0, 0.0, 0.0, heatCapacityGas, 900.0, 500.0, distances,
                             velocity, std::vector<double>(nodes, 1.0), gasDensity, coeffDiffusion, 0, zeros, {},
                             gasTemperature, gasTemperatureDot,
                             gasTemperature, solidTemperatureDot, wallTemperature, wallTemperatureDot, {}, {}, {},
                             inlet);
    double inflow = v * gasTemperature[gridPoints];
    for (size_t grid = 1; grid < nodes; ++grid) inflow += h * gasTemperatureDot[grid];
    const double conduction =
        inlet == Column::InletCondition::Dirichlet ? a * (gasTemperature[0] - gasTemperature[1]) / h : 0.0;
    EXPECT_NEAR(inflow, v * gasTemperature[0] + conduction, 1.0e-12 * v * 350.0) << static_cast<int>(inlet);
  }
}

TEST(DanckwertsInlet, SorptionJacobianCouplesTheFirstCellToTheInletFaceValue)
{
  // At a fixed velocity the upwind mass derivatives are linear in the concentrations, so differences in c_0 and c_1
  // are exact. Between the Dirichlet and the Danckwerts inlet, the Jacobian rows of the first cell must change by the
  // same amounts as these derivatives.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_danckwerts_jacobian.json";
  nlohmann::json data = beaInput();
  data["ConvectionScheme"] = "Upwind";
  data["NumberOfGridPoints"] = 12;
  data["BreakthroughIntegrator"] = "RungeKutta3";
  data["InletCondition"] = "Danckwerts";
  for (auto& component : data["Components"]) component["AxialDispersionCoefficient"] = 1.0e-3;
  std::ofstream(path) << data;
  Column column{InputReader(path.string())};
  column.initialize();
  Timing timings;
  precompute(column, timings);
  const size_t n = column.numberOfComponents;

  auto firstCellDerivative = [&](Column::InletCondition inlet, size_t node, size_t comp, size_t other)
  {
    std::vector<double> concentration(column.concentration.begin(), column.concentration.end());
    std::vector<double> plus(concentration.size());
    std::vector<double> minus(concentration.size());
    constexpr double step = 1.0e-2;
    concentration[node * n + other] += step;
    computeMassDerivatives(column.numberOfGridPoints, n, column.columnDistances, column.axialDispersionCoefficient,
                           column.interstitialGasVelocity, concentration, plus, column.bulkSpeciesSink,
                           column.fluidPhase, column.totalVoidFraction, Column::ConvectionScheme::Upwind, inlet);
    concentration[node * n + other] -= 2.0 * step;
    computeMassDerivatives(column.numberOfGridPoints, n, column.columnDistances, column.axialDispersionCoefficient,
                           column.interstitialGasVelocity, concentration, minus, column.bulkSpeciesSink,
                           column.fluidPhase, column.totalVoidFraction, Column::ConvectionScheme::Upwind, inlet);
    return (plus[n + comp] - minus[n + comp]) / (2.0 * step);
  };

  SorptionJacobian danckwerts(column);
  assembleSorptionJacobian(column, column.concentration, column.gasTemperature, danckwerts);
  column.inletCondition = Column::InletCondition::Dirichlet;
  SorptionJacobian dirichlet(column);
  assembleSorptionJacobian(column, column.concentration, column.gasTemperature, dirichlet);

  for (size_t comp = 0; comp < n; ++comp)
  {
    const size_t row = danckwerts.concentrationIndex(1, comp);
    for (size_t node : {size_t{0}, size_t{1}})
    {
      const size_t col = danckwerts.concentrationIndex(node, comp);
      const double expected = firstCellDerivative(Column::InletCondition::Danckwerts, node, comp, comp) -
                              firstCellDerivative(Column::InletCondition::Dirichlet, node, comp, comp);
      EXPECT_NEAR(danckwerts(row, col) - dirichlet(row, col), expected, 1.0e-8 * std::abs(dirichlet(row, col)))
          << comp << " " << node;
      EXPECT_GT(std::abs(expected), 0.0);
    }
  }
}

TEST(ContinuityClosure, KeepsTheTotalConcentrationWhenTheComponentsDisperseDifferently)
{
  // The continuity closures conserve the total molar flux of the species balances, uptake and net axial dispersion
  // included (C-BC-1), so sum_j dc_j/dt = 0 at every node also when the carrier does not disperse and the
  // adsorbing components disperse each with their own D_j.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_closure_dispersion.json";
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 16;
  data["OutletPressure"] = 1.2e5;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["AxialDispersionCoefficient"] = 2.0e-4 * (1.0 + 0.5 * static_cast<double>(comp));
  }
  for (const std::string boundary : {"InletPressureInletVelocity", "InletVelocityOutletPressure"})
  {
    for (const std::string scheme : {"WENO3", "Upwind"})
    {
      for (const std::string inlet : {"Danckwerts", "Dirichlet"})
      {
        data["BoundaryCondition"] = boundary;
        data["ConvectionScheme"] = scheme;
        data["InletCondition"] = inlet;
        std::ofstream(path) << data;
        Column column(InputReader(path.string()));
        column.initialize();

        // A steep front of every adsorbing component in the first half of the column, carrier to the balance.
        const size_t n = column.numberOfComponents;
        const double total = 1.2e5 / (R * column.gasTemperature[0]);
        for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
        {
          const double z = column.columnDistances[grid] / column.columnLength;
          const double front = 0.5 + 0.5 * std::tanh(20.0 * (0.4 - z));
          double adsorbing = 0.0;
          for (size_t comp = 1; comp < n; ++comp)
          {
            column.concentration[grid * n + comp] = column.concentration[comp] * front + 1.0e-3;
            adsorbing += column.concentration[grid * n + comp];
          }
          column.concentration[grid * n] = total - adsorbing;
        }

        computeDerivatives(column);

        double scale = 0.0;
        for (double rate : column.concentrationDot) scale = std::max(scale, std::abs(rate));
        ASSERT_GT(scale, 0.0);
        for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
        {
          double sumDot = 0.0;
          for (size_t comp = 0; comp < n; ++comp) sumDot += column.concentrationDot[grid * n + comp];
          EXPECT_NEAR(sumDot, 0.0, 1.0e-10 * scale) << boundary << " " << scheme << " " << inlet << " " << grid;
        }
      }
    }
  }
}

TEST(MultibedColumn, ConservesTheSuperficialMolarFlowWhereTheVoidFractionChanges)
{
  // The fluxes are superficial, eps v c (C-MBD-1): the helium inventory sum_i eps_i h_i c_i of a column whose void
  // fraction steps from 0.38 to 0.42 changes by the feed flux eps_0 v_0 c_0 less the outflow eps_N v_N c_N, while
  // CO2 and N2 adsorb in a front that crosses the interface.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_multibed_superficial_flux.json";
  std::ifstream input(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  for (const std::string boundary : {"InletPressureInletVelocity", "FixedVelocity"})
  {
    for (const std::string scheme : {"WENO3", "Upwind"})
    {
      data["BoundaryCondition"] = boundary;
      data["ConvectionScheme"] = scheme;
      std::ofstream(path) << data;
      Column column(InputReader(path.string()));
      column.initialize();
      ASSERT_NE(column.totalVoidFraction[1], column.totalVoidFraction[column.numberOfGridPoints]);

      const size_t n = column.numberOfComponents;
      const double total = column.totalConcentration[0];
      for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
      {
        const double z = column.columnDistances[grid] / column.columnLength;
        const double front = 0.5 + 0.5 * std::tanh(20.0 * (0.6 - z));
        double adsorbing = 0.0;
        for (size_t comp = 1; comp < n; ++comp)
        {
          column.concentration[grid * n + comp] = column.concentration[comp] * front + 1.0e-3;
          adsorbing += column.concentration[grid * n + comp];
        }
        column.concentration[grid * n] = total - adsorbing;
      }

      computeDerivatives(column);

      const size_t outlet = column.numberOfGridPoints;
      double accumulation = 0.0;
      for (size_t grid = 1; grid <= outlet; ++grid)
      {
        accumulation += column.totalVoidFraction[grid] * gridSpacing(column.columnDistances, grid) *
                        column.concentrationDot[grid * n];
      }
      const double inflow =
          column.totalVoidFraction[0] * column.interstitialGasVelocity[0] * column.concentration[0];
      const double outflow =
          column.totalVoidFraction[outlet] * column.interstitialGasVelocity[outlet] * column.concentration[outlet * n];
      EXPECT_NEAR(accumulation, inflow - outflow, 1.0e-12 * inflow) << boundary << " " << scheme;
      // The superficial velocity, not the interstitial one, is continuous at the interface.
      if (boundary == "FixedVelocity")
      {
        EXPECT_NEAR(column.totalVoidFraction[outlet] * column.interstitialGasVelocity[outlet],
                    column.totalVoidFraction[0] * column.interstitialGasVelocity[0], 1.0e-15);
      }
    }
  }
}

TEST(ContinuityClosure, CarriesTheThermalExpansionOfTheGas)
{
  // With the energy balance the closures carry the thermal expansion of the gas (B-BC-2): the total molar flux
  // leaving node i removes -(c_i/T_i) dT_i/dt, so that sum_j dc_ij/dt = -(sum_j c_ij / T_i) dT_i/dt at every node
  // and sum_j c_j T stays constant at constant pressure.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_closure_expansion.json";
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 16;
  data["OutletPressure"] = 1.2e5;
  data["energyBalance"] = true;
  data["wallDensity"] = 7800.0;
  data["heatCapacityWall"] = 500.0;
  data["heatCapacityGas"] = 1000.0;
  data["heatCapacitySolid"] = 900.0;
  data["gasThermalConductivity"] = 0.05;
  data["heatTransferGasSolid"] = 80.0;
  data["heatTransferGasWall"] = 40.0;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["AxialDispersionCoefficient"] = 2.0e-4 * (1.0 + 0.5 * static_cast<double>(comp));
  }
  for (const std::string boundary : {"InletPressureInletVelocity", "InletVelocityOutletPressure"})
  {
    for (const std::string scheme : {"WENO3", "Upwind"})
    {
      for (const std::string inlet : {"Danckwerts", "Dirichlet"})
      {
        data["BoundaryCondition"] = boundary;
        data["ConvectionScheme"] = scheme;
        data["InletCondition"] = inlet;
        std::ofstream(path) << data;
        Column column(InputReader(path.string()));
        column.initialize();

        // An adsorption front and a heat front: the gas is hotter behind the front, the solid hotter still, and the
        // wall at the feed temperature, so that every node heats or cools.
        const size_t n = column.numberOfComponents;
        const double feedTemperature = column.gasTemperature[0];
        const double total = 1.2e5 / (R * feedTemperature);
        for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
        {
          const double z = column.columnDistances[grid] / column.columnLength;
          const double front = 0.5 + 0.5 * std::tanh(20.0 * (0.4 - z));
          double adsorbing = 0.0;
          for (size_t comp = 1; comp < n; ++comp)
          {
            column.concentration[grid * n + comp] = column.concentration[comp] * front + 1.0e-3;
            adsorbing += column.concentration[grid * n + comp];
          }
          column.concentration[grid * n] = total - adsorbing;
          column.gasTemperature[grid] = feedTemperature + 30.0 * front + 5.0 * z;
          column.solidTemperature[grid] = column.gasTemperature[grid] + 12.0 * front;
          column.wallTemperature[grid] = feedTemperature;
        }

        computeDerivatives(column);

        double scale = 0.0;
        for (double rate : column.concentrationDot) scale = std::max(scale, std::abs(rate));
        ASSERT_GT(scale, 0.0);
        double expansionScale = 0.0;
        for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
        {
          double sumDot = 0.0;
          double sum = 0.0;
          for (size_t comp = 0; comp < n; ++comp)
          {
            sumDot += column.concentrationDot[grid * n + comp];
            sum += column.concentration[grid * n + comp];
          }
          const double expansion = -sum / column.gasTemperature[grid] * column.gasTemperatureDot[grid];
          expansionScale = std::max(expansionScale, std::abs(expansion));
          EXPECT_NEAR(sumDot, expansion, 1.0e-10 * scale) << boundary << " " << scheme << " " << inlet << " " << grid;
        }
        EXPECT_GT(expansionScale, 1.0e-3 * scale) << boundary << " " << scheme << " " << inlet;
      }
    }
  }
}

TEST(ColumnInitialization, FillsTheGasAtTheFeedTemperature)
{
  // The column starts at InfluxTemperature, so the initial and feed concentrations are P/(R T_in), also when the
  // ambient Temperature differs (C-BC-2; before, they took the ambient temperature).
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_initial_feed_temperature.json";
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 8;
  data["Temperature"] = 298.0;
  data["InfluxTemperature"] = 318.0;
  std::ofstream(path) << data;
  Column column(InputReader(path.string()));
  column.initialize();

  const size_t n = column.numberOfComponents;
  for (size_t grid = 0; grid <= column.numberOfGridPoints; ++grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < n; ++comp) sum += column.concentration[grid * n + comp];
    EXPECT_DOUBLE_EQ(column.gasTemperature[grid], 318.0) << grid;
    EXPECT_NEAR(R * column.gasTemperature[grid] * sum / column.totalPressure[grid], 1.0, 1.0e-6) << grid;
  }
}

TEST(AxialDispersionLaw, ScalesTheMolecularPartAndAddsTheMechanicalPart)
{
  // D_j = D_0,j (T/T_ref)^1.75 (P_ref/P) + lambda |v|, with the gas pressure P = R T sum_j c_j of the node.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_dispersion_law.json";
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 8;
  data["BoundaryCondition"] = "FixedVelocity";
  for (size_t comp = 0; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["AxialDispersionCoefficient"] = 1.0e-5 * static_cast<double>(comp + 1);
  }

  // Without a dispersivity and a reference state, the dispersion is the input value itself.
  std::ofstream(path) << data;
  {
    Column column(InputReader(path.string()));
    column.initialize();
    computeDerivatives(column);
    for (size_t index = 0; index < column.axialDispersionCoefficient.size(); ++index)
    {
      EXPECT_EQ(column.axialDispersionCoefficient[index], column.molecularDispersionCoefficient[index]);
    }
  }

  data["Geometry"]["AxialDispersivity"] = 1.5e-3;
  data["AxialDispersionReferenceTemperature"] = 300.0;
  data["AxialDispersionReferencePressure"] = 1.0e5;
  std::ofstream(path) << data;
  Column column(InputReader(path.string()));
  column.initialize();
  computeDerivatives(column);
  const size_t n = column.numberOfComponents;
  for (size_t grid = 0; grid <= column.numberOfGridPoints; ++grid)
  {
    const double temperature = column.gasTemperature[grid];
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < n; ++comp) concentrationSum += column.concentration[grid * n + comp];
    const double scaling = std::pow(temperature / 300.0, 1.75) * 1.0e5 / (R * temperature * concentrationSum);
    for (size_t comp = 0; comp < n; ++comp)
    {
      const double expected = 1.0e-5 * static_cast<double>(comp + 1) * scaling +
                              1.5e-3 * std::abs(column.interstitialGasVelocity[grid]);
      EXPECT_NEAR(column.axialDispersionCoefficient[grid * n + comp], expected, 1.0e-12 * expected) << grid << comp;
    }
  }
}

TEST(AxialDispersionLaw, AddsTheTaylorArisDispersionOfTheMonolithChannels)
{
  // D_j = D_0,j + v^2 d_h^2 / (192 D_0,j) in laminar channels; it divides by D_0, which every component needs.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_taylor_aris.json";
  std::ifstream input(examplePath("examples/Chemisorption-CO2-N2/monolith/simulation.json"));
  nlohmann::json data = nlohmann::json::parse(input);
  data["NumberOfGridPoints"] = 8;
  data["BoundaryCondition"] = "FixedVelocity";
  data["Geometry"]["TaylorArisDispersion"] = true;
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  const std::vector<double> molecular{7.0e-5, 1.6e-5, 2.0e-5};  // He, CO2, N2
  for (size_t comp = 0; comp < molecular.size(); ++comp)
  {
    data["Components"][comp]["AxialDispersionCoefficient"] = molecular[comp];
  }
  std::ofstream(path) << data;
  Column column(InputReader(path.string()));
  column.initialize();
  computeDerivatives(column);
  const double hydraulicDiameter = column.geometries[0].dimensions.hydraulicDiameter;
  EXPECT_DOUBLE_EQ(hydraulicDiameter, 1.0e-3);  // square channel: d_h is the side
  const size_t n = column.numberOfComponents;
  for (size_t grid = 0; grid <= column.numberOfGridPoints; ++grid)
  {
    const double v = column.interstitialGasVelocity[grid];
    for (size_t comp = 0; comp < n; ++comp)
    {
      const double expected = molecular[comp] + v * v * hydraulicDiameter * hydraulicDiameter / (192.0 * molecular[comp]);
      EXPECT_NEAR(column.axialDispersionCoefficient[grid * n + comp], expected, 1.0e-12 * expected) << grid << comp;
    }
  }
}

TEST(AxialDispersionLawInput, RejectsInconsistentSettingsAndOverridesTheDispersivityPerAdsorbent)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_dispersion_input.json";
  nlohmann::json data = beaInput();
  auto rejected = [&](const nlohmann::json& input)
  {
    std::ofstream(path) << input;
    EXPECT_THROW(InputReader{path.string()}, std::runtime_error) << input.dump();
  };

  nlohmann::json onlyTemperature = data;
  onlyTemperature["AxialDispersionReferenceTemperature"] = 300.0;
  rejected(onlyTemperature);
  nlohmann::json negativePressure = onlyTemperature;
  negativePressure["AxialDispersionReferencePressure"] = -1.0e5;
  rejected(negativePressure);
  nlohmann::json negativeDispersivity = data;
  negativeDispersivity["Geometry"]["AxialDispersivity"] = -1.0e-3;
  rejected(negativeDispersivity);
  nlohmann::json taylorArisInPackedBed = data;
  taylorArisInPackedBed["Geometry"]["TaylorArisDispersion"] = true;
  rejected(taylorArisInPackedBed);

  std::ifstream input(examplePath("examples/Multibed/CO2-N2-polishing/simulation.json"));
  nlohmann::json multibed = nlohmann::json::parse(input);
  multibed["Geometry"]["AxialDispersivity"] = 1.0e-3;
  multibed["Adsorbents"][1]["AxialDispersivity"] = 2.0e-3;
  std::ofstream(path) << multibed;
  const InputReader reader(path.string());
  ASSERT_EQ(reader.adsorbentGeometries.size(), 2U);
  EXPECT_DOUBLE_EQ(reader.adsorbentGeometries[0].dispersion.dispersivity, 1.0e-3);
  EXPECT_DOUBLE_EQ(reader.adsorbentGeometries[1].dispersion.dispersivity, 2.0e-3);
}

namespace
{
// The mirror image of node-major values on a grid of gridPoints cells: cell i becomes cell N + 1 - i; node 0 keeps the
// value given (node 0 is not used in a reverse step).
std::vector<double> mirroredCells(const std::vector<double>& values, size_t gridPoints, size_t width, double nodeZero)
{
  std::vector<double> mirrored(values.size(), nodeZero);
  for (size_t cell = 1; cell <= gridPoints; ++cell)
  {
    for (size_t k = 0; k < width; ++k) mirrored[cell * width + k] = values[(gridPoints + 1 - cell) * width + k];
  }
  return mirrored;
}

// Face j becomes face N - j with the opposite velocity.
std::vector<double> mirroredFaces(const std::vector<double>& velocity)
{
  std::vector<double> mirrored(velocity.size());
  for (size_t face = 0; face < velocity.size(); ++face) mirrored[face] = -velocity[velocity.size() - 1 - face];
  return mirrored;
}

// The input of ContinuityClosure.CarriesTheThermalExpansionOfTheGas with the given direction and settings, without
// pressure drop, so that the pressure is uniform and a reverse step is the exact mirror image of a forward one, and
// without wall conduction, which couples the wall to node 0 at z = 0 only.
Column mirrorTestColumn(nlohmann::json data, const std::string& direction)
{
  data["FlowDirection"] = direction;
  const std::filesystem::path path = testTemporaryDirectory() / ("ruptura_mirror_" + direction + ".json");
  std::ofstream(path) << data;
  Column column(InputReader(path.string()));
  for (Geometry& geometry : column.geometries) geometry.pressureDrop = PressureDropLaw{};
  column.initialize();
  return column;
}
}  // namespace

TEST(FlowReversal, MassDerivativesAreTheMirrorImageOfTheForwardOnes)
{
  // On a uniform grid a reverse step is a forward step mirrored: cell i of the forward column is cell N + 1 - i of the
  // reverse one, face i is face N - i with the opposite velocity, and the feed of node 0 is the product-end gas.
  constexpr size_t gridPoints = 14;
  constexpr size_t componentCount = 2;
  constexpr size_t nodes = gridPoints + 1;
  const std::vector<double> distances = uniformDistances(gridPoints, 0.28);
  const std::vector<double> feed{4.0, 32.0};
  std::vector<double> concentration(nodes * componentCount);
  std::vector<double> velocity(nodes);
  std::vector<double> dispersion(nodes * componentCount);
  std::vector<double> sink(nodes * componentCount);
  concentration[0] = feed[0];
  concentration[1] = feed[1];
  for (size_t grid = 0; grid < nodes; ++grid)
  {
    const double z = distances[grid] / 0.28;
    velocity[grid] = 0.05 * (1.0 - 0.2 * z);
    for (size_t comp = 0; comp < componentCount; ++comp)
    {
      dispersion[grid * componentCount + comp] = 2.0e-3 * (1.0 + 0.1 * z) * (1.0 + static_cast<double>(comp));
    }
    if (grid == 0) continue;
    concentration[grid * componentCount] = 3.0 * (0.5 + 0.5 * std::tanh(12.0 * (0.4 - z))) + 0.1;
    concentration[grid * componentCount + 1] = 36.0 - concentration[grid * componentCount];
    sink[grid * componentCount] = 0.02 * concentration[grid * componentCount];
    sink[grid * componentCount + 1] = -0.01 * concentration[grid * componentCount];
  }
  const std::vector<double> ones(nodes, 1.0);
  std::vector<double> wenoScratch(nodes * (componentCount + 4));
  // Node 0 of the reverse step holds a value that must not enter.
  const std::vector<double> reverseConcentration = mirroredCells(concentration, gridPoints, componentCount, 1.0e3);
  const std::vector<double> reverseDispersion = mirroredCells(dispersion, gridPoints, componentCount, 7.0);
  const std::vector<double> reverseSink = mirroredCells(sink, gridPoints, componentCount, 0.0);
  const std::vector<double> reverseVelocity = mirroredFaces(velocity);
  const Column::FlowEnds reverse{Column::FlowDirection::Reverse, feed, 0.0};

  for (Column::ConvectionScheme scheme : {Column::ConvectionScheme::WENO3, Column::ConvectionScheme::Upwind})
  {
    for (Column::FluidPhase phase : {Column::FluidPhase::Gas, Column::FluidPhase::Liquid})
    {
      for (Column::InletCondition inlet : {Column::InletCondition::Danckwerts, Column::InletCondition::Dirichlet})
      {
        std::vector<double> forwardDot(concentration.size());
        std::vector<double> reverseDot(concentration.size());
        computeMassDerivatives(gridPoints, componentCount, distances, dispersion, velocity, concentration, forwardDot,
                               sink, phase, ones, scheme, inlet, wenoScratch);
        computeMassDerivatives(gridPoints, componentCount, distances, reverseDispersion, reverseVelocity,
                               reverseConcentration, reverseDot, reverseSink, phase, ones, scheme, inlet, wenoScratch,
                               reverse);
        double scale = 0.0;
        for (double rate : forwardDot) scale = std::max(scale, std::abs(rate));
        ASSERT_GT(scale, 0.0);
        for (size_t cell = 1; cell <= gridPoints; ++cell)
        {
          for (size_t comp = 0; comp < componentCount; ++comp)
          {
            EXPECT_NEAR(reverseDot[cell * componentCount + comp],
                        forwardDot[(gridPoints + 1 - cell) * componentCount + comp], 1.0e-12 * scale)
                << static_cast<int>(scheme) << " " << static_cast<int>(phase) << " " << static_cast<int>(inlet) << " "
                << cell << " " << comp;
          }
        }
      }
    }
  }
}

TEST(FlowReversal, TemperatureDerivativesAreTheMirrorImageOfTheForwardOnes)
{
  // The gas, solid, and wall temperatures of a reverse step with the feed temperature at z = L, without wall
  // conduction (the wall node 0 conducts into the column at z = 0 only).
  const Geometry geometry = makeGeometry(PackedBedTubeSpec{0.4, 3.0e-3, 0.025, 0.04});
  constexpr size_t gridPoints = 10;
  constexpr size_t nodes = gridPoints + 1;
  const std::vector<double> distances = uniformDistances(gridPoints, 0.2);
  std::vector<double> gasTemperature(nodes), solidTemperature(nodes), wallTemperature(nodes);
  std::vector<double> velocity(nodes), gasDensity(nodes);
  for (size_t grid = 0; grid < nodes; ++grid)
  {
    const double z = distances[grid] / 0.2;
    velocity[grid] = 0.02 * (1.0 + 0.3 * z);
    gasDensity[grid] = 1.2 - 0.2 * z;
    gasTemperature[grid] = grid == 0 ? 350.0 : 300.0 + 40.0 * std::exp(-4.0 * z);
    solidTemperature[grid] = 305.0 + 30.0 * std::exp(-3.0 * z);
    wallTemperature[grid] = 300.0 + 5.0 * z;
  }
  const std::vector<double> zeros(nodes, 0.0);
  const std::vector<double> area(nodes, 600.0);
  const std::vector<double> particleDensity(nodes, 1000.0);
  const std::vector<Geometry> geometries{geometry};
  const std::vector<double> reverseGas = mirroredCells(gasTemperature, gridPoints, 1, 900.0);
  const std::vector<double> reverseSolid = mirroredCells(solidTemperature, gridPoints, 1, 900.0);
  const std::vector<double> reverseWall = mirroredCells(wallTemperature, gridPoints, 1, 900.0);
  const std::vector<double> reverseDensity = mirroredCells(gasDensity, gridPoints, 1, 1.2);
  const std::vector<double> reverseVelocity = mirroredFaces(velocity);
  const Column::FlowEnds reverse{Column::FlowDirection::Reverse, {}, 350.0};

  for (Column::InletCondition inlet : {Column::InletCondition::Danckwerts, Column::InletCondition::Dirichlet})
  {
    std::vector<double> rates[2][3];
    for (int direction = 0; direction < 2; ++direction)
    {
      for (std::vector<double>& rate : rates[direction]) rate.assign(nodes, 0.0);
      std::vector<double> coeffDiffusion(nodes, 0.0);
      const bool backward = direction == 1;
      computeEnergyDerivatives(gridPoints, 1, 300.0, geometries, particleDensity, area, area, zeros, {}, 4000.0, 6.0,
                               0.0, 20.0, 30.0, 5.0, 1000.0, 900.0, 500.0, distances,
                               backward ? reverseVelocity : velocity, std::vector<double>(nodes, 1.0),
                               backward ? reverseDensity : gasDensity,
                               coeffDiffusion, 0, zeros, {}, backward ? reverseGas : gasTemperature, rates[direction][0],
                               backward ? reverseSolid : solidTemperature, rates[direction][1],
                               backward ? reverseWall : wallTemperature, rates[direction][2], {}, {}, {}, inlet,
                               backward ? reverse : Column::FlowEnds{});
    }
    for (size_t field = 0; field < 3; ++field)
    {
      double scale = 0.0;
      for (double rate : rates[0][field]) scale = std::max(scale, std::abs(rate));
      ASSERT_GT(scale, 0.0);
      for (size_t cell = 1; cell <= gridPoints; ++cell)
      {
        EXPECT_NEAR(rates[1][field][cell], rates[0][field][gridPoints + 1 - cell], 1.0e-12 * scale)
            << static_cast<int>(inlet) << " " << field << " " << cell;
      }
    }
  }
}

TEST(FlowReversal, ClosuresMirrorTheForwardStepAndCarryTheThermalExpansion)
{
  // The full right-hand side, the continuity closures included: a reverse step marches from z = L with the inlet
  // velocity there, gives the mirrored velocities and rates, and keeps sum_j c_ij at P/(R T_i).
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 16;
  data["OutletPressure"] = 1.2e5;
  data["energyBalance"] = true;
  data["wallDensity"] = 7800.0;
  data["heatCapacityWall"] = 500.0;
  data["heatCapacityGas"] = 1000.0;
  data["heatCapacitySolid"] = 900.0;
  data["gasThermalConductivity"] = 0.05;
  data["wallThermalConductivity"] = 0.0;
  data["heatTransferGasSolid"] = 80.0;
  data["heatTransferGasWall"] = 40.0;
  data["ColumnEntranceVelocity"] = 0.098;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["AxialDispersionCoefficient"] = 2.0e-4 * (1.0 + 0.5 * static_cast<double>(comp));
  }
  for (const std::string boundary : {"InletPressureInletVelocity", "InletVelocityOutletPressure", "FixedVelocity"})
  {
    for (const std::string scheme : {"WENO3", "Upwind"})
    {
      for (const std::string inlet : {"Danckwerts", "Dirichlet"})
      {
        SCOPED_TRACE(boundary + " " + scheme + " " + inlet);
        data["BoundaryCondition"] = boundary;
        data["ConvectionScheme"] = scheme;
        data["InletCondition"] = inlet;
        Column forward = mirrorTestColumn(data, "Forward");
        Column reverse = mirrorTestColumn(data, "Reverse");
        const size_t n = forward.numberOfComponents;
        const size_t gridPoints = forward.numberOfGridPoints;

        // The front and heat front of CarriesTheThermalExpansionOfTheGas, entering at z = 0, and loadings below
        // equilibrium; the reverse column holds its mirror image.
        const double feedTemperature = forward.gasTemperature[0];
        const double total = 1.2e5 / (R * feedTemperature);
        for (size_t grid = 1; grid <= gridPoints; ++grid)
        {
          const double z = forward.columnDistances[grid] / forward.columnLength;
          const double front = 0.5 + 0.5 * std::tanh(20.0 * (0.4 - z));
          double adsorbing = 0.0;
          for (size_t comp = 1; comp < n; ++comp)
          {
            forward.concentration[grid * n + comp] = forward.concentration[comp] * front + 1.0e-3;
            forward.physisorption[grid * n + comp] = 0.05 * front;
            adsorbing += forward.concentration[grid * n + comp];
          }
          forward.concentration[grid * n] = total - adsorbing;
          forward.gasTemperature[grid] = feedTemperature + 30.0 * front + 5.0 * z;
          forward.solidTemperature[grid] = forward.gasTemperature[grid] + 12.0 * front;
          forward.wallTemperature[grid] = feedTemperature;
        }
        ASSERT_EQ(reverse.productEndConcentration.size(), n);
        for (size_t comp = 0; comp < n; ++comp)
        {
          EXPECT_NEAR(reverse.productEndConcentration[comp], forward.concentration[comp], 1.0e-12 * total);
        }
        for (size_t cell = 1; cell <= gridPoints; ++cell)
        {
          const size_t mirror = gridPoints + 1 - cell;
          for (size_t comp = 0; comp < n; ++comp)
          {
            reverse.concentration[cell * n + comp] = forward.concentration[mirror * n + comp];
            reverse.physisorption[cell * n + comp] = forward.physisorption[mirror * n + comp];
          }
          reverse.gasTemperature[cell] = forward.gasTemperature[mirror];
          reverse.solidTemperature[cell] = forward.solidTemperature[mirror];
          reverse.wallTemperature[cell] = forward.wallTemperature[mirror];
        }

        Timing timings;
        for (Column* column : {&forward, &reverse})
        {
          precompute(*column, timings);
          computeDerivatives(*column);
        }

        double velocityScale = 0.0;
        for (double velocity : forward.interstitialGasVelocity) velocityScale = std::max(velocityScale, velocity);
        for (size_t face = 0; face <= gridPoints; ++face)
        {
          EXPECT_LT(reverse.interstitialGasVelocity[face], 0.0) << face;
          EXPECT_NEAR(reverse.interstitialGasVelocity[face], -forward.interstitialGasVelocity[gridPoints - face],
                      1.0e-10 * velocityScale)
              << face;
        }
        const auto expectMirrored = [&](std::span<const double> expected, std::span<const double> actual,
                                        size_t width, const char* field)
        {
          double scale = 0.0;
          for (double value : expected) scale = std::max(scale, std::abs(value));
          for (size_t cell = 1; cell <= gridPoints; ++cell)
          {
            for (size_t k = 0; k < width; ++k)
            {
              EXPECT_NEAR(actual[cell * width + k], expected[(gridPoints + 1 - cell) * width + k], 1.0e-10 * scale)
                  << field << " " << cell << " " << k;
            }
          }
        };
        // The column has no wall area, so the wall temperatures do not change.
        expectMirrored(forward.concentrationDot, reverse.concentrationDot, n, "concentration");
        expectMirrored(forward.physisorptionDot, reverse.physisorptionDot, n, "physisorption");
        expectMirrored(forward.gasTemperatureDot, reverse.gasTemperatureDot, 1, "gas temperature");
        expectMirrored(forward.solidTemperatureDot, reverse.solidTemperatureDot, 1, "solid temperature");
        expectMirrored(forward.wallTemperatureDot, reverse.wallTemperatureDot, 1, "wall temperature");

        if (boundary == "FixedVelocity") continue;
        double scale = 0.0;
        for (double rate : reverse.concentrationDot) scale = std::max(scale, std::abs(rate));
        for (size_t grid = 1; grid <= gridPoints; ++grid)
        {
          double sumDot = 0.0;
          double sum = 0.0;
          for (size_t comp = 0; comp < n; ++comp)
          {
            sumDot += reverse.concentrationDot[grid * n + comp];
            sum += reverse.concentration[grid * n + comp];
          }
          EXPECT_NEAR(sumDot, -sum / reverse.gasTemperature[grid] * reverse.gasTemperatureDot[grid], 1.0e-10 * scale)
              << grid;
        }
      }
    }
  }
}

TEST(FlowReversal, BackFlowAtTheProductEndKeepsTheTotalConcentration)
{
  // B-CYC-1: a column full of feed with empty adsorbent takes up more than the feed brings, so that the march from
  // z = 0 turns the velocity negative before the outlet. The gas that flows back in at z = L is that of the last cell,
  // or the product-end gas if one is set, and the closure and the species balances keep sum_j c_j at P/RT.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_back_flow.json";
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 16;
  data["OutletPressure"] = 1.2e5;
  data["ColumnEntranceVelocity"] = 0.02;
  for (const std::string boundary : {"InletPressureInletVelocity", "InletVelocityOutletPressure"})
  {
    for (const std::string scheme : {"WENO3", "Upwind"})
    {
      for (const std::string inlet : {"Danckwerts", "Dirichlet"})
      {
        double lastCellRate = 0.0;  // d c_1/dt of the last cell with the gas of the last cell flowing back in
        for (const bool productEndGas : {false, true})
        {
          SCOPED_TRACE(boundary + " " + scheme + " " + inlet + (productEndGas ? " product-end gas" : ""));
          data["BoundaryCondition"] = boundary;
          data["ConvectionScheme"] = scheme;
          data["InletCondition"] = inlet;
          std::ofstream(path) << data;
          Column column(InputReader(path.string()));
          column.initialize();
          const size_t n = column.numberOfComponents;
          for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
          {
            for (size_t comp = 0; comp < n; ++comp) column.concentration[grid * n + comp] = column.concentration[comp];
          }
          if (productEndGas)
          {
            // Helium at the total concentration of the feed.
            double total = 0.0;
            for (size_t comp = 0; comp < n; ++comp) total += column.concentration[comp];
            column.productEndConcentration.assign(n, 0.0);
            column.productEndConcentration[0] = total;
            column.productEndTemperature = column.gasTemperature[0];
          }
          Timing timings;
          precompute(column, timings);
          computeDerivatives(column);

          ASSERT_LT(column.interstitialGasVelocity[column.numberOfGridPoints], 0.0);
          EXPECT_GT(column.interstitialGasVelocity[0], 0.0);
          double scale = 0.0;
          for (double rate : column.concentrationDot) scale = std::max(scale, std::abs(rate));
          ASSERT_GT(scale, 0.0);
          for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
          {
            double sumDot = 0.0;
            for (size_t comp = 0; comp < n; ++comp) sumDot += column.concentrationDot[grid * n + comp];
            EXPECT_NEAR(sumDot, 0.0, 1.0e-10 * scale) << grid;
          }
          // Helium that flows back in at z = L dilutes the adsorbing components of the last cell.
          const double rate = column.concentrationDot[column.numberOfGridPoints * n + 1];
          if (productEndGas) EXPECT_LT(rate, lastCellRate);
          lastCellRate = rate;
        }
      }
    }
  }
}

TEST(FlowReversal, ReverseBreakthroughIsTheMirrorImageOfTheForwardOne)
{
  // RK3 at a fixed step follows the mirror image to rounding; ROS34PW2, whose stage matrix has the sorption Jacobian of
  // the reverse step, to its tolerance.
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 12;
  data["BoundaryCondition"] = "InletPressureInletVelocity";
  data["ColumnEntranceVelocity"] = 0.098;
  data["BreakthroughIntegrator"] = "RungeKutta3";
  for (auto& component : data["Components"]) component["AxialDispersionCoefficient"] = 2.0e-4;
  for (const std::string scheme : {"WENO3", "Upwind"})
  {
    SCOPED_TRACE(scheme);
    data["ConvectionScheme"] = scheme;
    Column forward = mirrorTestColumn(data, "Forward");
    Column reverse = mirrorTestColumn(data, "Reverse");
    Column forwardRosenbrock = forward;
    Column reverseRosenbrock = reverse;
    Timing timings;
    RungeKutta3 forwardRungeKutta(5.0e-3);
    RungeKutta3 reverseRungeKutta(5.0e-3);
    forwardRungeKutta.advance(forward, 0.0, 30.0, timings);
    reverseRungeKutta.advance(reverse, 0.0, 30.0, timings);
    ROS34PW2 forwardIntegrator(1.0e-7, 1.0e-9, 1.0e-10, 1.0e-7);
    ROS34PW2 reverseIntegrator(1.0e-7, 1.0e-9, 1.0e-10, 1.0e-7);
    forwardIntegrator.advance(forwardRosenbrock, 0.0, 30.0, timings);
    reverseIntegrator.advance(reverseRosenbrock, 0.0, 30.0, timings);

    const size_t n = forward.numberOfComponents;
    const size_t gridPoints = forward.numberOfGridPoints;
    double scale = 0.0;
    for (double value : forward.concentration) scale = std::max(scale, value);
    EXPECT_GT(forward.concentration[2 * n + 1], 1.0e-3 * forward.concentration[1]);  // the front has entered
    for (size_t cell = 1; cell <= gridPoints; ++cell)
    {
      for (size_t comp = 0; comp < n; ++comp)
      {
        const size_t index = cell * n + comp;
        const size_t mirror = (gridPoints + 1 - cell) * n + comp;
        EXPECT_NEAR(reverse.concentration[index], forward.concentration[mirror], 1.0e-9 * scale) << cell << comp;
        EXPECT_NEAR(reverseRosenbrock.concentration[index], forwardRosenbrock.concentration[mirror], 1.0e-5 * scale)
            << cell << comp;
      }
    }
  }
}

TEST(FlowDirectionInput, ReadsTheColumnAndPhaseDirectionsAndTheProductEndGas)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_flow_direction.json";
  nlohmann::json data = beaInput();
  std::ofstream(path) << data;
  EXPECT_EQ(Column(InputReader(path.string())).flowDirection, Column::FlowDirection::Forward);
  data["FlowDirection"] = "Reverse";
  std::ofstream(path) << data;
  {
    const InputReader reader(path.string());
    EXPECT_EQ(reader.flowDirection, 1U);
    Column column(reader);
    EXPECT_EQ(column.flowDirection, Column::FlowDirection::Reverse);
    column.initialize();
    EXPECT_LT(column.interstitialGasVelocity[0], 0.0);
    ASSERT_EQ(column.productEndConcentration.size(), column.numberOfComponents);
    EXPECT_GT(column.productEndConcentration[1], 0.0);  // the feed enters at z = L
    EXPECT_EQ(column.concentration[1], 0.0);            // node 0 holds the initial gas
  }
  data["FlowDirection"] = "Sideways";
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);

  data["FlowDirection"] = "Reverse";
  data["SimulationType"] = "SwingAdsorption";
  data.erase("SimulationTime");
  data["SwingAdsorptionPhases"] = {
      {{"Name", "Inherited"}, {"Duration", 10.0}},
      {{"Name", "Feed"},
       {"Duration", 10.0},
       {"FlowDirection", "Forward"},
       {"ProductEndGasPhaseMolFractions", {{"Helium", 2.0}}}},
      {{"Name", "Purge"}, {"Duration", 10.0}, {"FlowDirection", "Reverse"}, {"GasPhaseMolFractions", {{"Helium", 1.0}}}}};
  std::ofstream(path) << data;
  {
    const InputReader reader(path.string());
    ASSERT_EQ(reader.swingAdsorptionPhases.size(), 3U);
    EXPECT_FALSE(reader.swingAdsorptionPhases[0].flowDirection.has_value());
    EXPECT_EQ(reader.swingAdsorptionPhases[1].flowDirection, std::optional<size_t>{0});
    EXPECT_EQ(reader.swingAdsorptionPhases[2].flowDirection, std::optional<size_t>{1});
    ASSERT_EQ(reader.swingAdsorptionPhases[1].productEndGasPhaseMolFractions.size(), reader.components.size());
    EXPECT_DOUBLE_EQ(reader.swingAdsorptionPhases[1].productEndGasPhaseMolFractions[0], 1.0);
    const SwingAdsorption swing(reader);
    EXPECT_EQ(swing.subStages[0].flowDirection, Column::FlowDirection::Reverse);
    EXPECT_EQ(swing.subStages[1].flowDirection, Column::FlowDirection::Forward);
    EXPECT_EQ(swing.subStages[2].flowDirection, Column::FlowDirection::Reverse);
  }
  data["SwingAdsorptionPhases"][2]["ProductEndGasPhaseMolFractions"] = {{"Helium", 1.0}};
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);  // the feed of a reverse phase is GasPhaseMolFractions
  data["SwingAdsorptionPhases"][2].erase("ProductEndGasPhaseMolFractions");
  data["SwingAdsorptionPhases"][1]["ProductEndGasPhaseMolFractions"] = {{"Argon", 1.0}};
  std::ofstream(path) << data;
  EXPECT_THROW(InputReader{path.string()}, std::runtime_error);
}

TEST(FlowReversal, SwingPhaseFeedsAtTheProductEnd)
{
  // A feed step and a helium purge from z = L: in the purge every face carries gas towards z = 0, and the column is
  // cleaned from the product end first.
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_reverse_swing.json";
  nlohmann::json data = beaInput();
  data["NumberOfGridPoints"] = 12;
  data["BoundaryCondition"] = "InletPressureInletVelocity";
  data["ColumnEntranceVelocity"] = 0.098;
  data["BreakthroughIntegrator"] = "ROS34PW2";
  data["SimulationType"] = "SwingAdsorption";
  data["WriteInterval"] = 50.0;
  data.erase("SimulationTime");
  data["SwingAdsorptionPhases"] = {
      {{"Name", "Feed"}, {"Duration", 200.0}},
      {{"Name", "Purge"}, {"Duration", 100.0}, {"FlowDirection", "Reverse"}, {"GasPhaseMolFractions", {{"Helium", 1.0}}}}};
  std::ofstream(path) << data;
  SwingAdsorption swing{InputReader(path.string())};
  swing.run();

  const Column& column = swing.breakthrough.column;
  const size_t n = column.numberOfComponents;
  const size_t gridPoints = column.numberOfGridPoints;
  EXPECT_EQ(column.flowDirection, Column::FlowDirection::Reverse);
  for (size_t face = 0; face <= gridPoints; ++face) EXPECT_LT(column.interstitialGasVelocity[face], 0.0) << face;
  ASSERT_EQ(column.productEndConcentration.size(), n);
  EXPECT_EQ(column.productEndConcentration[1], 0.0);
  auto adsorbing = [&](size_t grid)
  {
    double sum = 0.0;
    for (size_t comp = 1; comp < n; ++comp) sum += column.concentration[grid * n + comp];
    return sum;
  };
  EXPECT_LT(adsorbing(gridPoints), 1.0e-2 * adsorbing(1));
}

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "column.h"
#include "inputreader.h"
#include "integrators/block_banded.h"
#include "integrators/rk3.h"
#include "integrators/ros34pw2.h"
#include "json.h"
#include "temporary_directory.h"
#include "timing.h"

namespace
{
std::string examplePath(const std::string& relativePath)
{
  return (std::filesystem::path{RUPTURA_SOURCE_DIR} / relativePath).string();
}

// An input file on 12 grid cells, with the given keys replaced.
InputReader smallInput(const std::string& path, const std::string& name, const nlohmann::json& overrides)
{
  std::ifstream input(examplePath(path));
  nlohmann::json data = nlohmann::json::parse(input);
  data["NumberOfGridPoints"] = 12;
  data["BreakthroughIntegrator"] = "RungeKutta3";
  data.erase("TimeStep");
  for (const auto& [key, value] : overrides.items()) data[key] = value;
  const std::filesystem::path file = testTemporaryDirectory() / ("ruptura_ros34pw2_" + name + ".json");
  std::ofstream(file) << data;
  return InputReader(file.string());
}

// Largest difference of the concentrations and loadings (and temperatures, if integrated) of two columns, each
// relative to the largest value of its block in the first column.
double relativeStateDifference(const Column& reference, const Column& other)
{
  double difference = 0.0;
  const auto compare = [&](std::span<const double> expected, std::span<const double> actual)
  {
    double scale = 0.0;
    for (double value : expected) scale = std::max(scale, std::abs(value));
    if (scale == 0.0) return;
    for (size_t index = 0; index < expected.size(); ++index)
    {
      difference = std::max(difference, std::abs(actual[index] - expected[index]) / scale);
    }
  };
  compare(reference.concentration, other.concentration);
  compare(reference.physisorption, other.physisorption);
  compare(reference.chemisorption, other.chemisorption);
  if (reference.energyBalance)
  {
    compare(reference.gasTemperature, other.gasTemperature);
    compare(reference.solidTemperature, other.solidTemperature);
  }
  return difference;
}

// RK3 at a small step and ROS34PW2 at tight tolerances over the same interval; returns their relative difference.
double differenceFromRungeKutta3(const InputReader& reader, double rungeKuttaStep, double endTime,
                                 ROS34PW2& rosenbrock)
{
  Column explicitColumn(reader);
  Column rosenbrockColumn(reader);
  explicitColumn.initialize();
  rosenbrockColumn.initialize();
  Timing timings;
  RungeKutta3 rungeKutta(rungeKuttaStep);
  rungeKutta.advance(explicitColumn, 0.0, endTime, timings);
  rosenbrock.advance(rosenbrockColumn, 0.0, endTime, timings);
  return relativeStateDifference(explicitColumn, rosenbrockColumn);
}
}  // namespace

TEST(ROS34PW2, StageCoefficientsFollowFromTheTablesOfRangAndAngermann)
{
  // A, Gamma, b, and the embedded b of PETSc's TSROSWRA34PW2.
  const double g = 4.3586652150845900e-01;
  const std::array<std::array<double, 4>, 4> tableA{
      {{0.0, 0.0, 0.0, 0.0},
       {8.7173304301691801e-01, 0.0, 0.0, 0.0},
       {8.4457060015369423e-01, -1.1299064236484185e-01, 0.0, 0.0},
       {0.0, 0.0, 1.0, 0.0}}};
  const std::array<std::array<double, 4>, 4> tableGamma{
      {{g, 0.0, 0.0, 0.0},
       {-8.7173304301691801e-01, g, 0.0, 0.0},
       {-9.0338057013044082e-01, 5.4180672388095326e-02, g, 0.0},
       {2.4212380706095346e-01, -1.2232505839045147e+00, 5.4526025533510214e-01, g}}};
  const std::array<double, 4> b{2.4212380706095346e-01, -1.2232505839045147e+00, 1.5452602553351020e+00, g};
  const std::array<double, 4> embedded{3.7810903145819369e-01, -9.6042292212423178e-02, 0.5, 2.1793326075422950e-01};

  // Gamma^-1 by forward substitution, column by column.
  std::array<std::array<double, 4>, 4> inverse{};
  for (size_t col = 0; col < 4; ++col)
  {
    for (size_t row = col; row < 4; ++row)
    {
      double sum = row == col ? 1.0 : 0.0;
      for (size_t k = col; k < row; ++k) sum -= tableGamma[row][k] * inverse[k][col];
      inverse[row][col] = sum / tableGamma[row][row];
    }
  }

  EXPECT_DOUBLE_EQ(ROS34PW2::gamma, g);
  size_t packed = 0;
  for (size_t row = 1; row < 4; ++row)
  {
    for (size_t col = 0; col < row; ++col, ++packed)
    {
      double a = 0.0;
      for (size_t k = 0; k < 4; ++k) a += tableA[row][k] * inverse[k][col];
      EXPECT_NEAR(ROS34PW2::a[packed], a, 1.0e-13);
      EXPECT_NEAR(ROS34PW2::c[packed], -inverse[row][col], 1.0e-13);
    }
  }
  for (size_t col = 0; col < 4; ++col)
  {
    double m = 0.0;
    double e = 0.0;
    double alpha = 0.0;
    for (size_t k = 0; k < 4; ++k)
    {
      m += b[k] * inverse[k][col];
      e += (b[k] - embedded[k]) * inverse[k][col];
      alpha += tableA[col][k];
    }
    EXPECT_NEAR(ROS34PW2::m[col], m, 1.0e-13);
    EXPECT_NEAR(ROS34PW2::e[col], e, 1.0e-13);
    EXPECT_NEAR(ROS34PW2::alpha[col], alpha, 1.0e-15);
  }
}

TEST(ROS34PW2, KeepsThirdOrderWhenTheStageMatrixOnlyApproximatesTheJacobian)
{
  // u1' = -u1^2 + u2 + cos(t), u2' = -u2, u(0) = (1, 1), with W = J, J/2, 0, and J without its coupling term. The time
  // dependence enters through the stage times only: dF/dt is part of the Jacobian that a W-method may approximate.
  const auto rightHandSide = [](double t, std::span<const double> u, std::span<double> f)
  {
    f[0] = -u[0] * u[0] + u[1] + std::cos(t);
    f[1] = -u[1];
  };
  // Reference at t = 1 from 200000 classical RK4 steps.
  std::array<double, 2> reference{1.0, 1.0};
  {
    const double h = 1.0 / 200000.0;
    std::array<double, 2> k1, k2, k3, k4, stage;
    for (size_t step = 0; step < 200000; ++step)
    {
      const double t = static_cast<double>(step) * h;
      rightHandSide(t, reference, k1);
      for (size_t i = 0; i < 2; ++i) stage[i] = reference[i] + 0.5 * h * k1[i];
      rightHandSide(t + 0.5 * h, stage, k2);
      for (size_t i = 0; i < 2; ++i) stage[i] = reference[i] + 0.5 * h * k2[i];
      rightHandSide(t + 0.5 * h, stage, k3);
      for (size_t i = 0; i < 2; ++i) stage[i] = reference[i] + h * k3[i];
      rightHandSide(t + h, stage, k4);
      for (size_t i = 0; i < 2; ++i) reference[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
  }

  using Approximation = std::array<double, 4>;  // W row-major, from the state.
  const std::vector<std::pair<std::string, std::function<Approximation(std::span<const double>)>>> approximations{
      {"exact", [](std::span<const double> u) { return Approximation{-2.0 * u[0], 1.0, 0.0, -1.0}; }},
      {"half", [](std::span<const double> u) { return Approximation{-u[0], 0.5, 0.0, -0.5}; }},
      {"zero", [](std::span<const double>) { return Approximation{0.0, 0.0, 0.0, 0.0}; }},
      {"uncoupled", [](std::span<const double> u) { return Approximation{-2.0 * u[0], 0.0, 0.0, -1.0}; }}};

  for (const auto& [label, approximation] : approximations)
  {
    SCOPED_TRACE(label);
    std::array<double, 3> errors{};
    for (size_t refinement = 0; refinement < errors.size(); ++refinement)
    {
      const size_t steps = 40 * (size_t{1} << refinement);
      const double h = 1.0 / static_cast<double>(steps);
      std::vector<double> u{1.0, 1.0};
      std::vector<double> function(2), stageState(2);
      std::array<std::vector<double>, 4> increments{std::vector<double>(2), std::vector<double>(2),
                                                    std::vector<double>(2), std::vector<double>(2)};
      for (size_t step = 0; step < steps; ++step)
      {
        // I - gamma h W, inverted directly.
        const double t = static_cast<double>(step) * h;
        const Approximation w = approximation(u);
        const double gammaStep = ROS34PW2::gamma * h;
        const double m00 = 1.0 - gammaStep * w[0], m01 = -gammaStep * w[1];
        const double m10 = -gammaStep * w[2], m11 = 1.0 - gammaStep * w[3];
        const double determinant = m00 * m11 - m01 * m10;
        rightHandSide(t, u, function);
        ROS34PW2::stageIncrements(
            t, h, u, function, stageState, {increments[0], increments[1], increments[2], increments[3]},
            [&](double stageTime, std::span<const double> stage, std::span<double> stageFunction)
            { rightHandSide(stageTime, stage, stageFunction); },
            [&](std::span<double> x)
            {
              const double x0 = (m11 * x[0] - m01 * x[1]) / determinant;
              const double x1 = (m00 * x[1] - m10 * x[0]) / determinant;
              x[0] = x0;
              x[1] = x1;
            });
        for (size_t i = 0; i < 2; ++i)
        {
          for (size_t stage = 0; stage < 4; ++stage) u[i] += ROS34PW2::m[stage] * increments[stage][i];
        }
      }
      errors[refinement] = std::max(std::abs(u[0] - reference[0]), std::abs(u[1] - reference[1]));
    }
    for (size_t refinement = 1; refinement < errors.size(); ++refinement)
    {
      const double order = std::log2(errors[refinement - 1] / errors[refinement]);
      EXPECT_GT(order, 2.85);
      EXPECT_LT(order, 3.15);
    }
  }
}

TEST(BlockBandedMatrix, SolvesLikeADenseEliminationWithTwoBlocksBelowTheDiagonal)
{
  // Seven 3 x 3 blocks, two below and one above the diagonal, with deterministic entries and a dominant diagonal.
  constexpr size_t blocks = 7;
  constexpr size_t size = 3;
  constexpr size_t n = blocks * size;
  BlockBandedMatrix matrix(blocks, size, 2);
  std::vector<double> dense(n * n, 0.0);
  for (size_t blockRow = 0; blockRow < blocks; ++blockRow)
  {
    for (size_t blockColumn = blockRow >= 2 ? blockRow - 2 : 0; blockColumn < std::min(blocks, blockRow + 2);
         ++blockColumn)
    {
      for (size_t row = 0; row < size; ++row)
      {
        for (size_t col = 0; col < size; ++col)
        {
          const size_t i = blockRow * size + row;
          const size_t j = blockColumn * size + col;
          const double value = i == j ? 4.0 + 0.1 * static_cast<double>(i)
                                      : std::sin(static_cast<double>(3 * i + 7 * j + 1)) * 0.5;
          matrix(blockRow, blockColumn, row, col) = value;
          dense[i + j * n] = value;
        }
      }
    }
  }
  std::vector<double> expected(n);
  for (size_t i = 0; i < n; ++i) expected[i] = std::cos(static_cast<double>(i));
  std::vector<double> rightHandSide(n, 0.0);
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j) rightHandSide[i] += dense[i + j * n] * expected[j];
  }

  ASSERT_TRUE(matrix.factorize());
  matrix.solve(rightHandSide);
  for (size_t i = 0; i < n; ++i) EXPECT_NEAR(rightHandSide[i], expected[i], 1.0e-12) << i;
}

TEST(BlockBandedMatrix, SolvesLikeADenseEliminationWithTwoBlocksOnEitherSide)
{
  // As above with two blocks above the diagonal too, the band of WENO3 convection in either direction.
  constexpr size_t blocks = 8;
  constexpr size_t size = 3;
  constexpr size_t n = blocks * size;
  BlockBandedMatrix matrix(blocks, size, 2, 2);
  std::vector<double> dense(n * n, 0.0);
  for (size_t blockRow = 0; blockRow < blocks; ++blockRow)
  {
    for (size_t blockColumn = blockRow >= 2 ? blockRow - 2 : 0; blockColumn < std::min(blocks, blockRow + 3);
         ++blockColumn)
    {
      for (size_t row = 0; row < size; ++row)
      {
        for (size_t col = 0; col < size; ++col)
        {
          const size_t i = blockRow * size + row;
          const size_t j = blockColumn * size + col;
          const double value = i == j ? 5.0 + 0.1 * static_cast<double>(i)
                                      : std::sin(static_cast<double>(3 * i + 7 * j + 1)) * 0.5;
          matrix(blockRow, blockColumn, row, col) = value;
          dense[i + j * n] = value;
        }
      }
    }
  }
  std::vector<double> expected(n);
  for (size_t i = 0; i < n; ++i) expected[i] = std::cos(static_cast<double>(i));
  std::vector<double> rightHandSide(n, 0.0);
  for (size_t i = 0; i < n; ++i)
  {
    for (size_t j = 0; j < n; ++j) rightHandSide[i] += dense[i + j * n] * expected[j];
  }

  ASSERT_TRUE(matrix.factorize());
  matrix.solve(rightHandSide);
  for (size_t i = 0; i < n; ++i) EXPECT_NEAR(rightHandSide[i], expected[i], 1.0e-12) << i;
}

TEST(BlockBandedMatrix, ReportsASingularDiagonalBlock)
{
  BlockBandedMatrix matrix(3, 2, 1);
  matrix.scaleAndShift(1.0, 1.0);
  matrix(1, 1, 0, 0) = 0.0;
  matrix(1, 1, 0, 1) = 0.0;
  EXPECT_FALSE(matrix.factorize());
}

TEST(ROS34PW2Column, MatchesRungeKutta3ThroughAWeno3AdsorptionFrontWithLargeSteps)
{
  // BEA C7 with its pressure-driven velocity and WENO3 convection; the exact transport Jacobian lets ROS34PW2 take
  // steps of several times the cell crossing time.
  InputReader reader = smallInput("examples/BEA-alkanes-C7/breakthrough/simulation.json", "bea", {});
  ROS34PW2 tight(1.0e-7, 1.0e-9, 1.0e-10, 1.0e-7);
  EXPECT_LT(differenceFromRungeKutta3(reader, 2.0e-3, 60.0, tight), 2.0e-5);
  // One Jacobian per step: a rejected step is repeated with the Jacobian of its start.
  EXPECT_EQ(tight.numberOfJacobians, tight.numberOfSteps);

  // At a relative tolerance of 1e-5 the mean step exceeds the 0.44 s in which the gas crosses a cell.
  ROS34PW2 loose(1.0e-5, 1.0e-6, 1.0e-7, 1.0e-4);
  EXPECT_LT(differenceFromRungeKutta3(reader, 2.0e-3, 60.0, loose), 1.0e-3);
  EXPECT_LT(loose.numberOfSteps, 120U);
}

TEST(ROS34PW2Column, IntegratesTheEnergyBalanceAndContinuityLikeRungeKutta3)
{
  // Experiment T15: Sips isotherms, an energy balance, and the velocity from continuity.
  InputReader reader = smallInput("marx2/ternary/T15/simulation.json", "t15", {});
  ROS34PW2 rosenbrock(1.0e-7, 1.0e-9, 1.0e-10, 1.0e-7);
  EXPECT_LT(differenceFromRungeKutta3(reader, 1.0e-2, 120.0, rosenbrock), 2.0e-5);
}

TEST(ROS34PW2Column, IntegratesChemisorptionLikeRungeKutta3)
{
  // The monolith of Stampi-Bombelli et al.: Toth physisorption and chemisorption, dominated by axial dispersion. RK3
  // needs steps far below its stability limit (3e-2 s on this grid) to converge here (TODO.md B-ITG-2).
  InputReader reader =
      smallInput("stampi-bombelli/dl_from_stampibombelli/monolith_400ppm_1mmol/simulation.json", "monolith", {});
  ROS34PW2 rosenbrock(1.0e-7, 1.0e-9, 1.0e-10, 1.0e-7);
  EXPECT_LT(differenceFromRungeKutta3(reader, 2.5e-4, 40.0, rosenbrock), 2.0e-5);
}

TEST(ROS34PW2Column, IntegratesALiquidWithADanckwertsInletLikeRungeKutta3)
{
  InputReader reader = smallInput("examples/Liquid/pH-dependent/simulation.json", "liquid",
                                  {{"InletCondition", "Danckwerts"}});
  ROS34PW2 rosenbrock(1.0e-7, 1.0e-9, 1.0e-10, 1.0e-7);
  EXPECT_LT(differenceFromRungeKutta3(reader, 1.0e-2, 2000.0, rosenbrock), 2.0e-5);
}

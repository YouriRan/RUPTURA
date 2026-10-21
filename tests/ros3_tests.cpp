#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iomanip>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "integrators/ros3.h"

namespace
{
struct Problem
{
  Rosenbrock3Rhs rhs;
  Rosenbrock3Jacobian jacobian;
  Rosenbrock3Rhs partialTimeDerivative;
};

Problem scalarLinear(double rate, double forcing = 0.0)
{
  return {[=](double, std::span<const double> state, std::span<double> derivative)
          { derivative[0] = rate * state[0] + forcing; },
          [=](double, std::span<const double>, std::span<double> jacobian) { jacobian[0] = rate; },
          {}};
}

std::vector<double> integrate(const Problem& problem, std::vector<double> state, size_t steps, double start = 0.0,
                              double end = 1.0)
{
  const double h = (end - start) / static_cast<double>(steps);
  for (size_t step = 0; step < steps; ++step)
  {
    state = rosenbrock3Step(start + static_cast<double>(step) * h, h, state, problem.rhs, problem.jacobian,
                            problem.partialTimeDerivative)
                .state;
  }
  return state;
}

double maximumError(std::span<const double> state, std::span<const double> exact)
{
  double error = 0.0;
  for (size_t row = 0; row < state.size(); ++row) error = std::max(error, std::abs(state[row] - exact[row]));
  return error;
}

std::string preciseNumber(double value)
{
  std::ostringstream output;
  output << std::setprecision(17) << value;
  return output.str();
}

void expectThirdOrder(const Problem& problem, const std::vector<double>& initial, const std::vector<double>& exact,
                      double start = 0.0, double end = 1.0)
{
  std::array<double, 3> errors;
  for (size_t refinement = 0; refinement < errors.size(); ++refinement)
  {
    const size_t steps = 40 * (size_t{1} << refinement);
    errors[refinement] = maximumError(integrate(problem, initial, steps, start, end), exact);
    ASSERT_GT(errors[refinement], 1.0e-14);
    ::testing::Test::RecordProperty("error_" + std::to_string(steps), preciseNumber(errors[refinement]));
  }
  for (size_t refinement = 1; refinement < errors.size(); ++refinement)
  {
    const double order = std::log2(errors[refinement - 1] / errors[refinement]);
    SCOPED_TRACE("refinement " + std::to_string(refinement));
    ::testing::Test::RecordProperty("order_" + std::to_string(refinement), preciseNumber(order));
    EXPECT_GT(order, 2.90);
    EXPECT_LT(order, 3.10);
  }
}

// Independent stability function derived by imposing third order and L-stability
// on a rational function with three identical poles, rather than using stage coefficients.
template <typename T>
T stabilityFunction(T z)
{
  const double gamma = Rosenbrock3Coefficients::gamma;
  const T denominator = 1.0 - gamma * z;
  return (1.0 + (1.0 - 3.0 * gamma) * z + (0.5 - 3.0 * gamma + 3.0 * gamma * gamma) * z * z) /
         (denominator * denominator * denominator);
}
}  // namespace

TEST(Rosenbrock3, CoefficientsGiveThirdOrderAndAnAStableRationalFunction)
{
  const double gamma = Rosenbrock3Coefficients::gamma;
  // Vanishing cubic numerator coefficient gives the L-stable third-order root.
  EXPECT_NEAR(6.0 * gamma * gamma * gamma - 18.0 * gamma * gamma + 9.0 * gamma - 1.0, 0.0, 4.0e-15);
  EXPECT_GT(gamma, 0.0);  // All poles are on the positive real axis.
  const double quadratic = 0.5 - 3.0 * gamma + 3.0 * gamma * gamma;
  // On the imaginary axis, |denominator|^2 - |numerator|^2 equals
  // (3*gamma^4 - quadratic^2)*y^4 + gamma^6*y^6, which is nonnegative.
  EXPECT_GT(3.0 * std::pow(gamma, 4) - quadratic * quadratic, 0.0);
}

TEST(Rosenbrock3, ScalarDecayConvergesAtThirdOrder) { expectThirdOrder(scalarLinear(-2.0), {1.0}, {std::exp(-2.0)}); }

TEST(Rosenbrock3, AffineRelaxationConvergesAtThirdOrder)
{
  constexpr double rate = 3.0;
  constexpr double equilibrium = 2.5;
  expectThirdOrder(scalarLinear(-rate, rate * equilibrium), {0.3},
                   {equilibrium + (0.3 - equilibrium) * std::exp(-rate)});
}

TEST(Rosenbrock3, NonsymmetricCoupledLinearSystemConvergesAtThirdOrder)
{
  const Problem problem{[](double, std::span<const double> state, std::span<double> derivative)
                        {
                          derivative[0] = -state[0] + 2.0 * state[1];
                          derivative[1] = -3.0 * state[1];
                        },
                        [](double, std::span<const double>, std::span<double> jacobian)
                        {
                          jacobian[0] = -1.0;
                          jacobian[2] = 2.0;  // Entry (row 0, column 1), not a transposed Jacobian.
                          jacobian[3] = -3.0;
                        },
                        {}};
  expectThirdOrder(problem, {1.0, 1.0}, {2.0 * std::exp(-1.0) - std::exp(-3.0), std::exp(-3.0)});
}

TEST(Rosenbrock3, NonlinearStateDependentJacobianConvergesAtThirdOrder)
{
  const Problem problem{
      [](double, std::span<const double> state, std::span<double> derivative) { derivative[0] = -state[0] * state[0]; },
      [](double, std::span<const double> state, std::span<double> jacobian) { jacobian[0] = -2.0 * state[0]; },
      {}};
  expectThirdOrder(problem, {1.0}, {0.5});
}

TEST(Rosenbrock3, TimeDependentForcingConvergesAtThirdOrder)
{
  const Problem problem{[](double time, std::span<const double> state, std::span<double> derivative)
                        { derivative[0] = -2.0 * state[0] + std::sin(time); },
                        [](double, std::span<const double>, std::span<double> jacobian) { jacobian[0] = -2.0; },
                        [](double time, std::span<const double>, std::span<double> derivative)
                        { derivative[0] = std::cos(time); }};
  expectThirdOrder(problem, {1.0}, {(2.0 * std::sin(1.0) - std::cos(1.0)) / 5.0 + 1.2 * std::exp(-2.0)});
}

TEST(Rosenbrock3, TimeAndStateDependentProblemConvergesAtThirdOrder)
{
  const Problem problem{[](double time, std::span<const double> state, std::span<double> derivative)
                        { derivative[0] = std::sin(time) * state[0]; },
                        [](double time, std::span<const double>, std::span<double> jacobian)
                        { jacobian[0] = std::sin(time); },
                        [](double time, std::span<const double> state, std::span<double> derivative)
                        { derivative[0] = std::cos(time) * state[0]; }};
  expectThirdOrder(problem, {1.0}, {std::exp(std::cos(0.7) - std::cos(1.7))}, 0.7, 1.7);
}

TEST(Rosenbrock3, MainLocalErrorIsFourthOrderAndEmbeddedErrorIsThirdOrder)
{
  const auto problem = scalarLinear(-1.0);
  const std::vector<double> initial{1.0};
  std::array<double, 3> mainErrors, estimatedErrors, embeddedErrors;
  for (size_t refinement = 0; refinement < mainErrors.size(); ++refinement)
  {
    const double h = 0.1 / static_cast<double>(size_t{1} << refinement);
    const auto result = rosenbrock3Step(0.0, h, initial, problem.rhs, problem.jacobian);
    mainErrors[refinement] = std::abs(result.state[0] - std::exp(-h));
    estimatedErrors[refinement] = std::abs(result.error[0]);
    embeddedErrors[refinement] = std::abs(result.state[0] - result.error[0] - std::exp(-h));
  }
  EXPECT_GT(mainErrors[1] / mainErrors[2], 14.0);
  EXPECT_LT(mainErrors[1] / mainErrors[2], 17.0);
  EXPECT_GT(estimatedErrors[1] / estimatedErrors[2], 7.0);
  EXPECT_LT(estimatedErrors[1] / estimatedErrors[2], 9.0);
  EXPECT_GT(embeddedErrors[1] / embeddedErrors[2], 7.0);
  EXPECT_LT(embeddedErrors[1] / embeddedErrors[2], 9.0);
}

TEST(Rosenbrock3, StiffDecayMatchesIndependentStabilityFunctionAndTendsToZero)
{
  const auto problem = scalarLinear(-1.0);
  const std::vector<double> initial{1.0};
  for (const double h : {1.0e-6, 1.0e-3, 0.1, 1.0, 10.0, 1.0e3, 1.0e6, 1.0e12})
  {
    SCOPED_TRACE(h);
    const auto result = rosenbrock3Step(0.0, h, initial, problem.rhs, problem.jacobian);
    EXPECT_NEAR(result.state[0], stabilityFunction(-h), 2.0e-14);
    EXPECT_LE(std::abs(result.state[0]), 1.0 + 2.0e-14);
    if (h == 1.0e12) EXPECT_LT(std::abs(result.state[0]), 1.0e-10);
  }
}

TEST(Rosenbrock3, CoupledOscillationsRemainStableAcrossTheLeftHalfPlane)
{
  for (const double radius : {1.0e-4, 0.01, 0.1, 1.0, 10.0, 1.0e2, 1.0e4, 1.0e8})
  {
    for (size_t angleIndex = 0; angleIndex <= 24; ++angleIndex)
    {
      const double angle = std::numbers::pi * static_cast<double>(angleIndex) / 24.0;
      const double real = -radius * std::sin(angle);
      const double imaginary = radius * std::cos(angle);
      const Problem problem{[=](double, std::span<const double> state, std::span<double> derivative)
                            {
                              derivative[0] = real * state[0] - imaginary * state[1];
                              derivative[1] = imaginary * state[0] + real * state[1];
                            },
                            [=](double, std::span<const double>, std::span<double> jacobian)
                            {
                              jacobian[0] = jacobian[3] = real;
                              jacobian[1] = imaginary;
                              jacobian[2] = -imaginary;
                            },
                            {}};
      const std::vector<double> initial{1.0, 0.0};
      const auto result = rosenbrock3Step(0.0, 1.0, initial, problem.rhs, problem.jacobian);
      const auto exact = stabilityFunction(std::complex<double>{real, imaginary});
      SCOPED_TRACE("radius " + std::to_string(radius) + ", angle " + std::to_string(angleIndex));
      EXPECT_NEAR(result.state[0], exact.real(), 4.0e-14);
      EXPECT_NEAR(result.state[1], exact.imag(), 4.0e-14);
      EXPECT_LE(std::hypot(result.state[0], result.state[1]), 1.0 + 4.0e-14);
    }
  }
}

TEST(Rosenbrock3, PivotedLinearSolveHandlesAZeroFirstDiagonal)
{
  const double first = 1.0 / Rosenbrock3Coefficients::gamma;
  const double eigenvalue = (first - 1.0) / 2.0 - std::sqrt(std::pow((first + 1.0) / 2.0, 2) + 2.0);
  const std::vector<double> eigenvector{1.0, eigenvalue - first};
  const Problem problem{[=](double, std::span<const double> state, std::span<double> derivative)
                        {
                          derivative[0] = first * state[0] + state[1];
                          derivative[1] = 2.0 * state[0] - state[1];
                        },
                        [=](double, std::span<const double>, std::span<double> jacobian)
                        {
                          jacobian[0] = first;
                          jacobian[1] = 2.0;
                          jacobian[2] = 1.0;
                          jacobian[3] = -1.0;
                        },
                        {}};
  const auto result = rosenbrock3Step(0.0, 1.0, eigenvector, problem.rhs, problem.jacobian);
  for (size_t row = 0; row < eigenvector.size(); ++row)
  {
    EXPECT_NEAR(result.state[row], stabilityFunction(eigenvalue) * eigenvector[row], 5.0e-14);
  }
}

TEST(Rosenbrock3, CoupledAdsorptionPreservesTheTotalInventory)
{
  constexpr double capacityFactor = 3.0;
  constexpr double rate = 1000.0;
  constexpr double equilibrium = 2.0;
  const Problem problem{[](double, std::span<const double> state, std::span<double> derivative)
                        {
                          derivative[1] = rate * (equilibrium - state[1]);
                          derivative[0] = -capacityFactor * derivative[1];
                        },
                        [](double, std::span<const double>, std::span<double> jacobian)
                        {
                          jacobian[2] = capacityFactor * rate;
                          jacobian[3] = -rate;
                        },
                        {}};
  const std::vector<double> initial{10.0, 0.0};
  const auto result = rosenbrock3Step(0.0, 1.0, initial, problem.rhs, problem.jacobian);
  // Dense pivoting at h*rate = 1000 introduces roundoff in the conserved mode.
  EXPECT_NEAR(result.state[0] + capacityFactor * result.state[1], 10.0, 1.0e-12);
  EXPECT_NEAR(result.state[1], equilibrium * (1.0 - stabilityFunction(-rate)), 4.0e-14);
}

TEST(Rosenbrock3, ZeroRatesAndConstantSourcesHaveTheCorrectLimits)
{
  const std::vector<double> initial{-4.0};
  const auto stationary = scalarLinear(0.0);
  const auto unchanged = rosenbrock3Step(0.0, 0.3, initial, stationary.rhs, stationary.jacobian);
  EXPECT_EQ(unchanged.state, initial);
  EXPECT_DOUBLE_EQ(unchanged.error[0], 0.0);
  const auto source = scalarLinear(0.0, 2.0);
  const auto result = rosenbrock3Step(0.0, 0.3, initial, source.rhs, source.jacobian);
  EXPECT_NEAR(result.state[0], -3.4, 2.0e-15);
  EXPECT_NEAR(result.error[0], 0.0, 2.0e-15);
  EXPECT_EQ(initial, std::vector<double>{-4.0});
}

TEST(Rosenbrock3, EvaluatesTheAnalyticalCallbacksAtTheRequiredStatesAndTimes)
{
  size_t rhsCalls = 0, jacobianCalls = 0, timeDerivativeCalls = 0;
  const std::vector<double> initial{1.0};
  const auto rhs = [&](double time, std::span<const double> state, std::span<double> derivative)
  {
    EXPECT_DOUBLE_EQ(time, rhsCalls == 0 ? 0.7 : 0.7 + Rosenbrock3Coefficients::gamma * 0.1);
    derivative[0] = time * state[0];
    ++rhsCalls;
  };
  const auto jacobian = [&](double time, std::span<const double> state, std::span<double> matrix)
  {
    EXPECT_DOUBLE_EQ(time, 0.7);
    EXPECT_DOUBLE_EQ(state[0], initial[0]);
    matrix[0] = time;
    ++jacobianCalls;
  };
  const auto partialTimeDerivative = [&](double time, std::span<const double> state, std::span<double> derivative)
  {
    EXPECT_DOUBLE_EQ(time, 0.7);
    EXPECT_DOUBLE_EQ(state[0], initial[0]);
    derivative[0] = state[0];  // Partial time derivative, not dF/dt along the solution.
    ++timeDerivativeCalls;
  };
  rosenbrock3Step(0.7, 0.1, initial, rhs, jacobian, partialTimeDerivative);
  EXPECT_EQ(rhsCalls, 2U);
  EXPECT_EQ(jacobianCalls, 1U);
  EXPECT_EQ(timeDerivativeCalls, 1U);
  EXPECT_EQ(initial, std::vector<double>{1.0});
}

TEST(Rosenbrock3, RejectsInvalidInputsAndNonfiniteEvaluations)
{
  const auto problem = scalarLinear(-1.0);
  const std::vector<double> initial{1.0};
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  for (const double h : {0.0, -1.0, nan, infinity})
  {
    EXPECT_THROW(rosenbrock3Step(0.0, h, initial, problem.rhs, problem.jacobian), std::invalid_argument);
  }
  EXPECT_THROW(rosenbrock3Step(infinity, 1.0, initial, problem.rhs, problem.jacobian), std::invalid_argument);
  EXPECT_THROW(rosenbrock3Step(1.0e100, 1.0, initial, problem.rhs, problem.jacobian), std::invalid_argument);
  EXPECT_THROW(rosenbrock3Step(1.0e308, 1.0e308, initial, problem.rhs, problem.jacobian), std::invalid_argument);
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, {}, problem.rhs, problem.jacobian), std::invalid_argument);
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, initial, {}, problem.jacobian), std::invalid_argument);
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, initial, problem.rhs, {}), std::invalid_argument);
  const std::vector<double> invalidState{nan};
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, invalidState, problem.rhs, problem.jacobian), std::runtime_error);

  const auto invalid = [=](double, std::span<const double>, std::span<double> output) { output[0] = nan; };
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, initial, invalid, problem.jacobian), std::runtime_error);
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, initial, problem.rhs, invalid), std::runtime_error);
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0, initial, problem.rhs, problem.jacobian, invalid), std::runtime_error);
  EXPECT_EQ(initial, std::vector<double>{1.0});
}

TEST(Rosenbrock3, FailedTrialsPreserveTheInitialState)
{
  const std::vector<double> initial{1.0};
  const auto singular = scalarLinear(1.0);
  EXPECT_THROW(rosenbrock3Step(0.0, 1.0 / Rosenbrock3Coefficients::gamma, initial, singular.rhs, singular.jacobian),
               std::runtime_error);
  size_t evaluations = 0;
  const auto failingRhs = [&](double, std::span<const double> state, std::span<double> derivative)
  {
    if (++evaluations == 2) throw std::runtime_error("stage evaluation failed");
    derivative[0] = -state[0];
  };
  const auto problem = scalarLinear(-1.0);
  EXPECT_THROW(rosenbrock3Step(0.0, 0.1, initial, failingRhs, problem.jacobian), std::runtime_error);
  EXPECT_EQ(initial, std::vector<double>{1.0});
  EXPECT_NO_THROW(rosenbrock3Step(0.0, 0.1, initial, problem.rhs, problem.jacobian));
}

#include "ros3.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

extern "C"
{
  void dgetrf_(int* m, int* n, double* a, int* lda, int* ipiv, int* info);
  void dgetrs_(char* trans, int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
}

namespace
{
void requireFinite(std::span<const double> values, const char* context)
{
  if (!std::all_of(values.begin(), values.end(), [](double value) { return std::isfinite(value); }))
  {
    throw std::runtime_error(std::string{"ROS3: nonfinite "} + context);
  }
}
}  // namespace

Rosenbrock3Step rosenbrock3Step(double time, double timeStep, std::span<const double> state, const Rosenbrock3Rhs& rhs,
                                const Rosenbrock3Jacobian& jacobian, const Rosenbrock3Rhs& partialTimeDerivative)
{
  using Coefficients = Rosenbrock3Coefficients;
  const size_t size = state.size();
  if (size == 0 || size > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      size > std::vector<double>{}.max_size() / size)
  {
    throw std::invalid_argument("ROS3: invalid state dimension");
  }
  if (!rhs || !jacobian)
  {
    throw std::invalid_argument("ROS3: RHS and analytical Jacobian callbacks are required");
  }
  if (!std::isfinite(time) || !std::isfinite(timeStep) || timeStep <= 0.0 || !std::isfinite(time + timeStep) ||
      time + timeStep == time || time + Coefficients::gamma * timeStep == time)
  {
    throw std::invalid_argument("ROS3: time and a positive, representable time step are required");
  }
  requireFinite(state, "initial state");

  std::vector<double> function(size, 0.0);
  rhs(time, state, function);
  requireFinite(function, "RHS");

  std::vector<double> matrix(size * size, 0.0);
  jacobian(time, state, matrix);
  requireFinite(matrix, "Jacobian");

  std::vector<double> timeDerivative(size, 0.0);
  if (partialTimeDerivative)
  {
    partialTimeDerivative(time, state, timeDerivative);
    requireFinite(timeDerivative, "partial time derivative");
  }

  const double gammaStep = Coefficients::gamma * timeStep;
  for (double& value : matrix) value *= -gammaStep;
  for (size_t row = 0; row < size; ++row) matrix[row + row * size] += 1.0;
  requireFinite(matrix, "stage matrix");

  int dimension = static_cast<int>(size);
  int info = 0;
  std::vector<int> pivots(size);
  dgetrf_(&dimension, &dimension, matrix.data(), &dimension, pivots.data(), &info);
  if (info != 0)
  {
    throw std::runtime_error("ROS3: LU factorization failed (LAPACK info " + std::to_string(info) + ")");
  }
  requireFinite(matrix, "LU factors");

  auto solveStage = [&](std::span<double> increment)
  {
    requireFinite(increment, "stage right-hand side");
    char transpose = 'N';
    int rightHandSides = 1;
    dgetrs_(&transpose, &dimension, &rightHandSides, matrix.data(), &dimension, pivots.data(), increment.data(),
            &dimension, &info);
    if (info != 0)
    {
      throw std::runtime_error("ROS3: stage solve failed (LAPACK info " + std::to_string(info) + ")");
    }
    requireFinite(increment, "stage increment");
  };

  std::array<std::vector<double>, 3> increments{std::vector<double>(size), std::vector<double>(size),
                                                std::vector<double>(size)};
  for (size_t row = 0; row < size; ++row)
  {
    increments[0][row] =
        gammaStep * (function[row] + timeStep * Coefficients::timeDerivativeWeights[0] * timeDerivative[row]);
  }
  solveStage(increments[0]);

  std::vector<double> stageState(size);
  for (size_t row = 0; row < size; ++row)
  {
    stageState[row] = state[row] + Coefficients::a[0] * increments[0][row];
  }
  requireFinite(stageState, "stage state");
  std::fill(function.begin(), function.end(), 0.0);
  rhs(time + Coefficients::alpha[1] * timeStep, stageState, function);
  requireFinite(function, "stage RHS");

  for (size_t row = 0; row < size; ++row)
  {
    increments[1][row] =
        gammaStep * (function[row] + timeStep * Coefficients::timeDerivativeWeights[1] * timeDerivative[row]) +
        Coefficients::gamma * Coefficients::c[0] * increments[0][row];
  }
  solveStage(increments[1]);

  // a31 = a21, a32 = 0 and alpha3 = alpha2: the third stage reuses the second RHS evaluation.
  for (size_t row = 0; row < size; ++row)
  {
    increments[2][row] =
        gammaStep * (function[row] + timeStep * Coefficients::timeDerivativeWeights[2] * timeDerivative[row]) +
        Coefficients::gamma * (Coefficients::c[1] * increments[0][row] + Coefficients::c[2] * increments[1][row]);
  }
  solveStage(increments[2]);

  Rosenbrock3Step result{std::vector<double>(state.begin(), state.end()), std::vector<double>(size, 0.0)};
  for (size_t stage = 0; stage < increments.size(); ++stage)
  {
    for (size_t row = 0; row < size; ++row)
    {
      result.state[row] += Coefficients::m[stage] * increments[stage][row];
      result.error[row] += Coefficients::e[stage] * increments[stage][row];
    }
  }
  requireFinite(result.state, "trial solution");
  requireFinite(result.error, "error estimate");
  return result;
}

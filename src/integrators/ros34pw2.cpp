#include "ros34pw2.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include <print>
#include <stdexcept>
#include <string>

#include "rk3.h"
#include "utils.h"

ROS34PW2::ROS34PW2(const InputReader& inputReader)
    : ROS34PW2(inputReader.ros34pw2RelativeTolerance, inputReader.ros34pw2AbsoluteToleranceConcentration,
               inputReader.ros34pw2AbsoluteToleranceLoading, inputReader.ros34pw2AbsoluteToleranceTemperature,
               inputReader.ros34pw2MaximumTimeStep)
{
}

ROS34PW2::ROS34PW2(double relativeTolerance, double absoluteToleranceConcentration, double absoluteToleranceLoading,
                   double absoluteToleranceTemperature, double maximumTimeStep)
    : relativeTolerance(relativeTolerance),
      absoluteToleranceConcentration(absoluteToleranceConcentration),
      absoluteToleranceLoading(absoluteToleranceLoading),
      absoluteToleranceTemperature(absoluteToleranceTemperature),
      maximumTimeStep(maximumTimeStep)
{
}

void ROS34PW2::allocate(const Column& column)
{
  const size_t n = column.numberOfComponents;
  const size_t nodes = column.numberOfGridPoints + 1;
  layout = column.stateLayout();
  energyBalance = column.energyBalance;
  keepTotalConcentration = column.fluidPhase == Column::FluidPhase::Gas &&
                           (column.boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity ||
                            column.boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure);
  totalIncrement.assign(keepTotalConcentration ? nodes : 0, 0.0);
  stepExpansion.assign(keepTotalConcentration && energyBalance ? nodes : 0, 0.0);
  sorptionJacobian.emplace(column);
  loadingDiagonal.assign(nodes * n, 1.0);
  // A WENO3 face flux reaches two nodes upstream, an upwind one one node; upstream lies on either side with the sign of
  // the velocity.
  const size_t reach = column.convectionScheme == Column::ConvectionScheme::WENO3 ? 2 : 1;
  concentrationMatrix = BlockBandedMatrix(nodes, n, reach, reach);
  thermalJacobian = energyBalance ? BlockBandedMatrix(nodes, 3, 1) : BlockBandedMatrix();
  thermalMatrix = thermalJacobian;
  thermalVector.assign(energyBalance ? 3 * nodes : 0, 0.0);
  for (std::vector<double>* vector : {&state, &rates, &function, &stageState})
  {
    vector->assign(column.state.size(), 0.0);
  }
  for (std::vector<double>& increment : increments) increment.assign(column.state.size(), 0.0);

  // The tolerance blocks of CVODE: loadings and temperatures have their own absolute tolerance, every concentration
  // (bulk, surface, pore) that of the concentrations.
  absoluteTolerance.assign(column.state.size(), absoluteToleranceConcentration);
  double* base = absoluteTolerance.data();
  const auto fill = [](std::span<double> block, double value) { std::fill(block.begin(), block.end(), value); };
  fill(layout.physisorption(base), absoluteToleranceLoading);
  fill(layout.chemisorption(base), absoluteToleranceLoading);
  fill(layout.gasTemperature(base), absoluteToleranceTemperature);
  fill(layout.solidTemperature(base), absoluteToleranceTemperature);
  fill(layout.wallTemperature(base), absoluteToleranceTemperature);
}

bool ROS34PW2::factorize(double timeStep)
{
  const SorptionJacobian& jacobian = *sorptionJacobian;
  const size_t n = jacobian.numberOfComponents;
  const size_t nodes = jacobian.numberOfNodes;
  const double gammaStep = gamma * timeStep;

  // The loading rows of J~ are -k_i on the diagonal and k_i dq*_i/dc on the concentrations of their node, and a
  // concentration row couples to the loading of its own component and node only. Eliminating the loadings leaves
  // W_cc - W_cq W_qq^-1 W_qc on the concentrations, with the transport blocks of the neighbouring nodes.
  for (size_t node = 0; node < nodes; ++node)
  {
    for (size_t comp = 0; comp < n; ++comp)
    {
      const size_t q = jacobian.loadingIndex(node, comp);
      loadingDiagonal[node * n + comp] = 1.0 - gammaStep * jacobian(q, q);
    }
  }
  concentrationMatrix.setZero();
  for (size_t node = 0; node < nodes; ++node)
  {
    const size_t firstNode = node > concentrationMatrix.lowerBlocks ? node - concentrationMatrix.lowerBlocks : 0;
    for (size_t other = firstNode; other < std::min(nodes, node + concentrationUpperBlocks + 1); ++other)
    {
      for (size_t i = 0; i < n; ++i)
      {
        const size_t row = jacobian.concentrationIndex(node, i);
        for (size_t k = 0; k < n; ++k)
        {
          double value = (other == node && i == k ? 1.0 : 0.0) -
                         gammaStep * jacobian.entry(row, jacobian.concentrationIndex(other, k));
          if (other == node)
          {
            const size_t q = jacobian.loadingIndex(node, i);
            value -= gammaStep * gammaStep * jacobian(row, q) * jacobian(q, jacobian.concentrationIndex(node, k)) /
                     loadingDiagonal[node * n + i];
          }
          concentrationMatrix(node, other, i, k) = value;
        }
      }
    }
  }
  if (!concentrationMatrix.factorize(concentrationUpperBlocks)) return false;

  if (energyBalance)
  {
    thermalMatrix = thermalJacobian;
    thermalMatrix.scaleAndShift(-gammaStep, 1.0);
    if (!thermalMatrix.factorize()) return false;
  }
  stageGammaStep = gammaStep;
  return true;
}

void ROS34PW2::solve(std::span<double> x) const
{
  const SorptionJacobian& jacobian = *sorptionJacobian;
  const size_t n = jacobian.numberOfComponents;
  const size_t nodes = jacobian.numberOfNodes;
  const std::span<double> concentration = layout.concentration(x.data());
  const std::span<double> loading = layout.physisorption(x.data());
  const std::span<double> gasTemperature = layout.gasTemperature(x.data());
  for (size_t node = 1; node < totalIncrement.size(); ++node)
  {
    totalIncrement[node] = std::reduce(concentration.begin() + static_cast<std::ptrdiff_t>(node * n),
                                       concentration.begin() + static_cast<std::ptrdiff_t>((node + 1) * n));
  }

  if (energyBalance)
  {
    const std::span<double> solidTemperature = layout.solidTemperature(x.data());
    const std::span<double> wallTemperature = layout.wallTemperature(x.data());
    for (size_t node = 0; node < nodes; ++node)
    {
      thermalVector[3 * node] = gasTemperature[node];
      thermalVector[3 * node + 1] = solidTemperature[node];
      thermalVector[3 * node + 2] = wallTemperature[node];
    }
    thermalMatrix.solve(thermalVector);
    for (size_t node = 1; node < stepExpansion.size(); ++node)
    {
      totalIncrement[node] -= stepExpansion[node] * (thermalVector[3 * node] - gasTemperature[node]);
    }
    for (size_t node = 0; node < nodes; ++node)
    {
      gasTemperature[node] = thermalVector[3 * node];
      solidTemperature[node] = thermalVector[3 * node + 1];
      wallTemperature[node] = thermalVector[3 * node + 2];
    }
  }

  // Concentrations from the system with the loadings eliminated, then the loadings of every node.
  for (size_t node = 0; node < nodes; ++node)
  {
    for (size_t comp = 0; comp < n; ++comp)
    {
      const size_t index = node * n + comp;
      concentration[index] += stageGammaStep *
                              jacobian(jacobian.concentrationIndex(node, comp), jacobian.loadingIndex(node, comp)) *
                              loading[index] / loadingDiagonal[index];
    }
  }
  concentrationMatrix.solve(concentration);
  // The total concentration of each cell follows the closure (see the class description).
  for (size_t node = 1; node < totalIncrement.size(); ++node)
  {
    const auto first = concentration.begin() + static_cast<std::ptrdiff_t>(node * n);
    const double correction = totalIncrement[node] - std::reduce(first, first + static_cast<std::ptrdiff_t>(n));
    for (size_t comp = 0; comp < n; ++comp) first[comp] += correction * stepMoleFractions[node * n + comp];
  }
  for (size_t node = 0; node < nodes; ++node)
  {
    for (size_t comp = 0; comp < n; ++comp)
    {
      const size_t q = jacobian.loadingIndex(node, comp);
      double sum = loading[node * n + comp];
      for (size_t other = 0; other < n; ++other)
      {
        sum += stageGammaStep * jacobian(q, jacobian.concentrationIndex(node, other)) * concentration[node * n + other];
      }
      loading[node * n + comp] = sum / loadingDiagonal[node * n + comp];
    }
  }
}

void ROS34PW2::advance(Column& column, double startTime, double endTime, Timing& timings)
{
  auto totalTimer = timings.scoped(timings.total);
  if (endTime <= startTime) return;
  if (!sorptionJacobian || state.size() != column.state.size() ||
      sorptionJacobian->numberOfNodes != column.numberOfGridPoints + 1)
  {
    allocate(column);
  }

  // Only a pressure ramp makes the right-hand side depend on the time.
  const auto evaluate = [&](double stageTime, std::span<const double> stage, std::span<double> stageRates)
  {
    std::copy(stage.begin(), stage.end(), column.state.begin());
    clampNonnegative(column.state);
    column.advancePressureRamp(stageTime);
    precompute(column, timings);
    ++numberOfMixturePredictions;
    timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column); });
    ++numberOfEvaluations;
    std::copy(column.stateDot.begin(), column.stateDot.end(), stageRates.begin());
  };
  const auto solveStage = [&](std::span<double> x) { solve(x); };

  // A first step of a thousandth of the time in which the fluid crosses the shortest grid cell.
  if (timeStep <= 0.0)
  {
    double cellTime = endTime - startTime;
    for (size_t node = 1; node <= column.numberOfGridPoints; ++node)
    {
      const double velocity = std::abs(column.interstitialGasVelocity[node]);
      if (velocity > 0.0) cellTime = std::min(cellTime, gridSpacing(column.columnDistances, node) / velocity);
    }
    timeStep = 1.0e-3 * cellTime;
  }

  std::array<std::span<double>, stages> stageIncrementSpans;
  for (size_t stage = 0; stage < stages; ++stage) stageIncrementSpans[stage] = increments[stage];
  const size_t size = state.size();

  double time = startTime;
  while (time < endTime)
  {
    // The column holds the fields of the precompute at this state; J~ is assembled once per step.
    timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column); });
    ++numberOfEvaluations;
    std::copy(column.state.begin(), column.state.end(), state.begin());
    std::copy(column.stateDot.begin(), column.stateDot.end(), rates.begin());
    assembleSorptionJacobian(column, column.concentration, column.gasTemperature, *sorptionJacobian);
    if (keepTotalConcentration) stepMoleFractions.assign(column.moleFraction.begin(), column.moleFraction.end());
    for (size_t node = 1; node < stepExpansion.size(); ++node)
    {
      stepExpansion[node] = column.totalConcentration[node] / std::max(1e-10, column.gasTemperature[node]);
    }
    if (energyBalance) assembleThermalJacobian(column, thermalJacobian);
    concentrationUpperBlocks =
        std::any_of(column.interstitialGasVelocity.begin(), column.interstitialGasVelocity.end(),
                    [](double velocity) { return velocity < 0.0; })
            ? concentrationMatrix.upperBlocks
            : 1;
    ++numberOfJacobians;

    bool rejected = false;
    std::string stageError;  // the last error of a stage evaluation that rejected the step
    while (true)
    {
      // Equal steps of at most timeStep up to endTime, so that no sliver remains before an output.
      const double remaining = endTime - time;
      const double longest = maximumTimeStep > 0.0 ? std::min(timeStep, maximumTimeStep) : timeStep;
      const double stepsToEnd = std::max(1.0, std::ceil(remaining / longest * (1.0 - 1.0e-12)));
      const double h = remaining / stepsToEnd;
      const bool reachesEnd = stepsToEnd == 1.0;
      if (rejected && h <= 1.0e-12 * std::max(1.0, std::abs(time)))
      {
        throw std::runtime_error(std::format("ROS34PW2: the step size fell to {:.3g} s at t = {:.6g} s{}", h, time,
                                             stageError.empty() ? "" : "; last stage error: " + stageError));
      }

      // A stage state for which the closures or the mixture prediction fail (negative outlet pressure, no velocity
      // bracket, failed IAST inversion) rejects the step, as a failed factorization does (B-ITG-18).
      double errorNorm = std::numeric_limits<double>::infinity();
      bool stagesEvaluated = false;
      if (factorize(h))
      {
        std::copy(rates.begin(), rates.end(), function.begin());
        try
        {
          stageIncrements(time, h, state, function, stageState, stageIncrementSpans, evaluate, solveStage);
          stagesEvaluated = true;
        }
        catch (const std::runtime_error& error)
        {
          stageError = error.what();
        }
      }
      if (stagesEvaluated)
      {
        // The new state, and the weighted RMS norm of the error estimate, as in CVODE.
        double sum = 0.0;
        for (size_t row = 0; row < size; ++row)
        {
          double next = state[row];
          double error = 0.0;
          for (size_t stage = 0; stage < stages; ++stage)
          {
            next += m[stage] * increments[stage][row];
            error += e[stage] * increments[stage][row];
          }
          stageState[row] = next;
          const double scale =
              absoluteTolerance[row] + relativeTolerance * std::max(std::abs(state[row]), std::abs(next));
          sum += (error / scale) * (error / scale);
        }
        errorNorm = std::sqrt(sum / static_cast<double>(size));
      }

      if (std::isfinite(errorNorm) && errorNorm <= 1.0)
      {
        std::copy(stageState.begin(), stageState.end(), column.state.begin());
        clampNonnegative(column.state);
        time = reachesEnd ? endTime : time + h;
        column.advancePressureRamp(time);
        precompute(column, timings);
        ++numberOfMixturePredictions;
        ++numberOfSteps;
        // The embedded solution is of second order, so the error estimate scales with h^3; no growth right after a
        // rejected step, and a step shortened to land on endTime keeps the longer proposal when its error allows.
        const double factor =
            std::min(rejected ? 1.0 : 5.0, std::max(0.2, 0.9 * std::pow(std::max(errorNorm, 1.0e-10), -1.0 / 3.0)));
        timeStep = reachesEnd && factor >= 1.0 ? std::max(timeStep, h * factor) : h * factor;
        break;
      }

      // A rejected step is repeated from the step start, with the fields of the column recomputed there.
      ++numberOfRejectedSteps;
      rejected = true;
      timeStep = h * (std::isfinite(errorNorm) ? std::max(0.2, 0.9 * std::pow(errorNorm, -1.0 / 3.0)) : 0.25);
      std::copy(state.begin(), state.end(), column.state.begin());
      std::copy(rates.begin(), rates.end(), column.stateDot.begin());
      column.advancePressureRamp(time);
      precompute(column, timings);
      ++numberOfMixturePredictions;
      timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column); });
      ++numberOfEvaluations;
      std::copy(column.stateDot.begin(), column.stateDot.end(), rates.begin());
    }
  }
}

void ROS34PW2::printStatistics() const
{
  std::print(
      "ROS34PW2 statistics:\n"
      "  steps:                     {}\n"
      "  rejected steps:            {}\n"
      "  rhs evaluations:           {}\n"
      "  jacobian evaluations:      {}\n"
      "  mixture predictions:       {}\n"
      "  next step:                 {:e} s\n",
      numberOfSteps, numberOfRejectedSteps, numberOfEvaluations, numberOfJacobians, numberOfMixturePredictions,
      timeStep);
}

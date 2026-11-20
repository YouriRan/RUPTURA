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
  pressureWork = energyBalance && column.boundaryCondition == Column::BoundaryCondition::MomentumBalance;
  totalIncrement.assign(keepTotalConcentration || pressureWork ? nodes : 0, 0.0);
  stepExpansion.assign((keepTotalConcentration && energyBalance) || pressureWork ? nodes : 0, 0.0);
  coupledVector.assign(pressureWork ? (n + 3) * nodes : 0, 0.0);
  gasRateByVelocity.assign(pressureWork ? 2 * nodes : 0, 0.0);
  sorptionJacobian.emplace(column);
  loadingDiagonal.assign(nodes * n, 1.0);
  // A WENO3 face flux reaches two nodes upstream, an upwind one one node; upstream lies on either side with the sign of
  // the velocity.
  const size_t reach = column.convectionScheme != Column::ConvectionScheme::Upwind ? 2 : 1;
  concentrationMatrix = BlockBandedMatrix(nodes, n, reach, reach);
  coupledMatrix = pressureWork ? BlockBandedMatrix(nodes, n + 3, reach, reach) : BlockBandedMatrix();
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
  fill(layout.coolantTemperature(base), absoluteToleranceTemperature);
  fill(layout.shellLoading(base), absoluteToleranceLoading);
  // The stream amounts are integrals of the solution: they follow it (solve) and take no part in the step control.
  fill(layout.streamAmount(base), 1.0e30);
}

bool ROS34PW2::factorize(double h)
{
  const SorptionJacobian& jacobian = *sorptionJacobian;
  const size_t n = jacobian.numberOfComponents;
  const size_t nodes = jacobian.numberOfNodes;
  const double gammaStep = gamma * h;

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
  // With pressureWork the concentrations are the first entries of the node blocks of the coupled matrix.
  BlockBandedMatrix& matrix = pressureWork ? coupledMatrix : concentrationMatrix;
  matrix.setZero();
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
          matrix(node, other, i, k) = value;
        }
      }
    }
  }
  if (pressureWork)
  {
    const size_t gasT = n;
    for (size_t node = 0; node < nodes; ++node)
    {
      for (size_t other = node > 0 ? node - 1 : 0; other < std::min(nodes, node + 2); ++other)
      {
        for (size_t row = 0; row < 3; ++row)
        {
          for (size_t col = 0; col < 3; ++col)
          {
            coupledMatrix(node, other, n + row, n + col) =
                (other == node && row == col ? 1.0 : 0.0) - gammaStep * thermalJacobian(node, other, row, col);
          }
        }
        // The concentration rates by the gas temperature of the node, through the face velocities.
        for (size_t k = 0; node > 0 && k < n; ++k)
        {
          coupledMatrix(node, other, k, gasT) =
              -gammaStep * jacobian.temperatureCoupling[(3 * node + other + 1 - node) * n + k];
        }
      }
    }
    for (size_t node = 1; node < nodes; ++node)
    {
      // The convection of the gas temperature by its two face velocities, which follow the gas on either side.
      for (size_t face = node - 1; face <= node; ++face)
      {
        const double rateByVelocity = gasRateByVelocity[2 * node + face + 1 - node];
        for (size_t side = 0; side < 2 && face + side < nodes && rateByVelocity != 0.0; ++side)
        {
          for (size_t comp = 0; comp < n; ++comp)
          {
            coupledMatrix(node, face + side, gasT, comp) -=
                gammaStep * rateByVelocity * jacobian.velocityByConcentration[(2 * face + side) * n + comp];
          }
          coupledMatrix(node, face + side, gasT, gasT) -=
              gammaStep * rateByVelocity * jacobian.velocityByTemperature[2 * face + side];
        }
      }
      // The pressure work: stepExpansion times the rate of the total concentration, the sum of the concentration rows.
      const size_t firstNode = node > coupledMatrix.lowerBlocks ? node - coupledMatrix.lowerBlocks : 0;
      for (size_t other = firstNode; other < std::min(nodes, node + concentrationUpperBlocks + 1); ++other)
      {
        for (size_t col = 0; col <= gasT; ++col)
        {
          double sum = other == node && col < n ? -1.0 : 0.0;
          for (size_t k = 0; k < n; ++k) sum += coupledMatrix(node, other, k, col);
          coupledMatrix(node, other, gasT, col) += stepExpansion[node] * sum;
        }
      }
    }
    stageGammaStep = gammaStep;
    return coupledMatrix.factorize(concentrationUpperBlocks);
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

  if (energyBalance && !pressureWork)
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
  if (pressureWork)
  {
    // The coupled system; the loadings that the lines above eliminated from the concentrations also leave the pressure
    // work of the gas temperature, by what they added to the total concentration of the cell.
    const std::span<double> solidTemperature = layout.solidTemperature(x.data());
    const std::span<double> wallTemperature = layout.wallTemperature(x.data());
    const size_t blockSize = n + 3;
    for (size_t node = 0; node < nodes; ++node)
    {
      const auto first = concentration.begin() + static_cast<std::ptrdiff_t>(node * n);
      std::copy_n(first, n, coupledVector.begin() + static_cast<std::ptrdiff_t>(node * blockSize));
      coupledVector[node * blockSize + n] =
          gasTemperature[node] +
          stepExpansion[node] * (std::reduce(first, first + static_cast<std::ptrdiff_t>(n)) - totalIncrement[node]);
      coupledVector[node * blockSize + n + 1] = solidTemperature[node];
      coupledVector[node * blockSize + n + 2] = wallTemperature[node];
    }
    coupledMatrix.solve(coupledVector);
    for (size_t node = 0; node < nodes; ++node)
    {
      std::copy_n(coupledVector.begin() + static_cast<std::ptrdiff_t>(node * blockSize), n,
                  concentration.begin() + static_cast<std::ptrdiff_t>(node * n));
      gasTemperature[node] = coupledVector[node * blockSize + n];
      solidTemperature[node] = coupledVector[node * blockSize + n + 1];
      wallTemperature[node] = coupledVector[node * blockSize + n + 2];
    }
  }
  else
  {
    concentrationMatrix.solve(concentration);
  }
  // The total concentration of each cell follows the closure (see the class description).
  for (size_t node = 1; keepTotalConcentration && node < totalIncrement.size(); ++node)
  {
    const auto first = concentration.begin() + static_cast<std::ptrdiff_t>(node * n);
    const double correction = totalIncrement[node] - std::reduce(first, first + static_cast<std::ptrdiff_t>(n));
    for (size_t comp = 0; comp < n; ++comp)
      first[static_cast<std::ptrdiff_t>(comp)] += correction * stepMoleFractions[node * n + comp];
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

  // MomentumBalance: the flux through an end follows the velocity of its face, and that the change that W made to the
  // gas on either side of it; the stream amounts take it, so that they stay those of the implicit solution.
  const std::span<double> streamAmount = layout.streamAmount(x.data());
  for (size_t index = 0; index < 2 * n && index < streamAmount.size() && !jacobian.streamByVelocity.empty(); ++index)
  {
    const size_t face = index < n ? 0 : nodes - 1;
    double velocityChange = 0.0;
    for (size_t side = 0; side < 2 && face + side < nodes; ++side)
    {
      const size_t source = face + side;
      for (size_t comp = 0; comp < n; ++comp)
      {
        velocityChange +=
            jacobian.velocityByConcentration[(2 * face + side) * n + comp] * concentration[source * n + comp];
      }
      if (energyBalance) velocityChange += jacobian.velocityByTemperature[2 * face + side] * gasTemperature[source];
    }
    streamAmount[index] += stageGammaStep * jacobian.streamByVelocity[index] * velocityChange;
  }
}

double ROS34PW2::advance(Column& column, double startTime, double endTime, Timing& timings, bool singleStep)
{
  auto totalTimer = timings.scoped(timings.total);
  if (endTime <= startTime) return startTime;
  if (!sorptionJacobian || state.size() != column.state.size() ||
      sorptionJacobian->numberOfNodes != column.numberOfGridPoints + 1)
  {
    allocate(column);
  }

  // Only a pressure ramp makes the right-hand side depend on the time.
  const auto evaluate = [&](double stageTime, std::span<const double> stage, std::span<double> stageRates)
  {
    std::copy(stage.begin(), stage.end(), column.state.begin());
    clampNonnegative(std::span<double>(column.state).first(column.state.size() - column.streamAmount.size()));
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
    // The column holds the accepted state clamped to max(0, c), which would add the undershoots of WENO3 ahead of a
    // front at every step (B-ITG-17); while it still holds that clamp, the step starts from the unclamped state.
    const bool clampedOnly = std::ranges::equal(column.state, stageState, [](double held, double accepted)
                                                { return held == accepted || (accepted < 0.0 && held == 0.0); });
    std::ranges::copy(clampedOnly ? stageState : column.state, state.begin());
    std::copy(column.stateDot.begin(), column.stateDot.end(), rates.begin());
    stepCachedPressure = column.cachedPressure;
    stepCachedGrandPotential = column.cachedGrandPotential;
    assembleSorptionJacobian(column, column.concentration, column.gasTemperature, column.solidTemperature,
                             *sorptionJacobian);
    if (keepTotalConcentration) stepMoleFractions.assign(column.moleFraction.begin(), column.moleFraction.end());
    for (size_t node = 1; node < stepExpansion.size(); ++node)
    {
      stepExpansion[node] = column.totalConcentration[node] / std::max(1e-10, column.gasTemperature[node]);
      if (!pressureWork) continue;
      const double kappa = R * column.totalConcentration[node] /
                           (std::max(1e-10, column.gasDensity[node]) *
                            (column.gasHeatCapacity.empty() ? column.heatCapacityGas : column.gasHeatCapacity[node]));
      stepExpansion[node] = kappa / ((1.0 - kappa) * stepExpansion[node]);
    }
    if (energyBalance) assembleThermalJacobian(column, thermalJacobian, gasRateByVelocity);
    concentrationUpperBlocks =
        std::any_of(column.interstitialGasVelocity.begin(), column.interstitialGasVelocity.end(),
                    [](double velocity) { return velocity < 0.0; })
            ? concentrationMatrix.upperBlocks
            : 1;
    ++numberOfJacobians;

    // A ramped port pressure of MomentumBalance drives the stiff relaxation of the pressure through the time: dF/dt
    // from the right-hand side a little later, at the same state.
    timeDerivative.clear();
    if (column.boundaryCondition == Column::BoundaryCondition::MomentumBalance &&
        std::ranges::any_of(column.ports, [](const Column::Port& port)
                            { return port.ramp.shape != Column::PressureRamp::Shape::None; }))
    {
      // A short time on, as the clock resolves it: the pressure difference that drives the gas is small against the
      // pressure, and the difference quotient must stay within it.
      const double later = time + 1.0e-6 * std::min(timeStep, endTime - time);
      const double delta = later - time;
      timeDerivative.resize(size);
      evaluate(later, state, timeDerivative);
      for (size_t row = 0; row < size; ++row) timeDerivative[row] = (timeDerivative[row] - rates[row]) / delta;
    }

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
          stageIncrements(time, h, state, function, stageState, stageIncrementSpans, evaluate, solveStage,
                          timeDerivative);
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
        clampNonnegative(std::span<double>(column.state).first(column.state.size() - column.streamAmount.size()));
        time = reachesEnd ? endTime : time + h;
        column.advancePressureRamp(time);
        precompute(column, timings);
        ++numberOfMixturePredictions;
        ++numberOfSteps;
        // As CVODE does (B-ITG-19, B-ITG-25): stops when 100 000 steps in a row before endTime advanced the time by
        // less than 1e-3 of the time in which the fluid crosses the shortest grid cell, or of the step cap if that is
        // smaller. Without flow there is no crossing time.
        constexpr size_t stepsPerCrawlCheck = 100000;
        timeSinceEndTime = reachesEnd ? 0.0 : timeSinceEndTime + h;
        stepsSinceEndTime = reachesEnd ? 0 : stepsSinceEndTime + 1;
        if (stepsSinceEndTime == stepsPerCrawlCheck)
        {
          double cellTime = std::numeric_limits<double>::infinity();
          for (size_t node = 1; node <= column.numberOfGridPoints; ++node)
          {
            const double velocity = std::abs(column.interstitialGasVelocity[node]);
            if (velocity > 0.0) cellTime = std::min(cellTime, gridSpacing(column.columnDistances, node) / velocity);
          }
          if (std::isfinite(cellTime) && maximumTimeStep > 0.0) cellTime = std::min(cellTime, maximumTimeStep);
          const double meanStep = timeSinceEndTime / static_cast<double>(stepsPerCrawlCheck);
          if (std::isfinite(cellTime) && meanStep < 1.0e-3 * cellTime)
          {
            throw std::runtime_error(std::format(
                "ROS34PW2 crawls at t = {} s: {} steps advanced the time by {:.3e} s (mean step {:.3e} s, {:.1e} of "
                "the {:.3e} s in which the fluid crosses the smallest cell or of the step cap).",
                time, stepsPerCrawlCheck, timeSinceEndTime, meanStep, meanStep / cellTime, cellTime));
          }
          timeSinceEndTime = 0.0;
          stepsSinceEndTime = 0;
        }
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
      clampNonnegative(std::span<double>(column.state).first(column.state.size() - column.streamAmount.size()));
      std::copy(rates.begin(), rates.end(), column.stateDot.begin());
      if (!stagesEvaluated)
      {
        // A stage that failed left the mixture prediction with starting values from which that of the step start
        // need not converge again.
        column.cachedPressure = stepCachedPressure;
        column.cachedGrandPotential = stepCachedGrandPotential;
      }
      column.advancePressureRamp(time);
      try
      {
        precompute(column, timings);
      }
      catch (const std::runtime_error&)
      {
        // Stages that converged but failed the error test can leave starting values from which the mixture prediction
        // of the step start does not converge again (B-ITG-23): those of the step start do.
        column.cachedPressure = stepCachedPressure;
        column.cachedGrandPotential = stepCachedGrandPotential;
        precompute(column, timings);
      }
      ++numberOfMixturePredictions;
      timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column); });
      ++numberOfEvaluations;
      std::copy(column.stateDot.begin(), column.stateDot.end(), rates.begin());
    }
    if (singleStep) break;
  }
  return time;
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

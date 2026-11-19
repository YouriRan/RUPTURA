#include "swing_adsorption.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "rk3.h"
#include "utils.h"

namespace
{
constexpr std::array<std::string_view, 3> portModeNames{"Closed", "Pressure", "Velocity"};

/**
 * \brief Value of the metric of a stop criterion for a column state; the pressures are those of the gas, R T sum_j c_j.
 */
double stopMetric(const Column& column, const ColumnStateViews& views, const InputReader::StopCriterion& criterion)
{
  using Metric = InputReader::StopCriterion::Metric;
  const size_t n = column.numberOfComponents;
  const size_t cell = criterion.end == 0 ? 1 : column.numberOfGridPoints;
  const auto totalConcentration = [&](size_t grid)
  {
    const auto first = views.concentration.begin() + static_cast<std::ptrdiff_t>(grid * n);
    return std::reduce(first, first + static_cast<std::ptrdiff_t>(n));
  };
  switch (criterion.metric)
  {
    case Metric::Pressure:
      return R * views.gasTemperature[cell] * totalConcentration(cell);
    case Metric::AveragePressure:
    {
      double fluidVolume = 0.0;
      double sum = 0.0;
      for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
      {
        const double volume = column.totalVoidFraction[grid] * gridSpacing(column.columnDistances, grid);
        fluidVolume += volume;
        sum += volume * R * views.gasTemperature[grid] * totalConcentration(grid);
      }
      return sum / fluidVolume;
    }
    case Metric::MaximumTemperature:
      return *std::max_element(views.solidTemperature.begin() + 1, views.solidTemperature.end());
    case Metric::MinimumTemperature:
      return *std::min_element(views.solidTemperature.begin() + 1, views.solidTemperature.end());
    case Metric::MoleFraction:
      return views.concentration[cell * n + *criterion.component] / std::max(1e-30, totalConcentration(cell));
    case Metric::AmountIn:
    case Metric::AmountOut:
    {
      // The net amount that entered at the end; what left is its negative.
      const auto first = views.streamAmount.begin() + static_cast<std::ptrdiff_t>(criterion.end * n);
      const double entered = criterion.component.has_value()
                                 ? first[static_cast<std::ptrdiff_t>(*criterion.component)]
                                 : std::reduce(first, first + static_cast<std::ptrdiff_t>(n));
      return criterion.metric == Metric::AmountIn ? entered : -entered;
    }
  }
  return 0.0;
}
}  // namespace

SwingAdsorption::SwingAdsorption(const InputReader& inputReader)
    : breakthrough(inputReader),
      columnNames(inputReader.columnNames),
      buffers(inputReader.buffers),
      phasesPerCycle(inputReader.swingAdsorptionPhases.size() / inputReader.numberOfCycles /
                     std::max<size_t>(1, inputReader.columnNames.size())),
      steadyStateTolerance(inputReader.cyclicSteadyStateTolerance),
      steadyStateCycles(inputReader.cyclicSteadyStateCycles)
{
  if (inputReader.swingAdsorptionPhases.empty())
  {
    throw std::invalid_argument("SwingAdsorption: at least one swing adsorption phase must be provided.");
  }
  for (size_t bed = 1; bed < columnNames.size(); ++bed) otherColumns.emplace_back(inputReader);
  for (const InputReader::Buffer& buffer : buffers)
  {
    std::vector<double>& amounts = bufferAmounts.emplace_back(buffer.moleFractions);
    for (double& amount : amounts) amount *= buffer.pressure * buffer.volume / (R * buffer.temperature);
  }

  // With MomentumBalance the ports set the ends; the inlet pressure and velocity are then those of the port at z = 0.
  const bool momentumBalance = breakthrough.column.boundaryCondition == Column::BoundaryCondition::MomentumBalance;
  subStages.reserve(inputReader.swingAdsorptionPhases.size());
  for (const InputReader::SwingAdsorptionPhase& phase : inputReader.swingAdsorptionPhases)
  {
    const double temperature = phase.temperature.value_or(inputReader.temperature);
    const double pressure = momentumBalance ? phase.ports[0].pressure.value_or(inputReader.outletPressure)
                                            : phase.inletPressure.value_or(inputReader.inletPressure);
    if (!(phase.duration > 0.0))
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has a non-positive duration.");
    }
    if (temperature < 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid temperature.");
    }
    if (pressure <= 0.0)
    {
      throw std::invalid_argument("SwingAdsorption: phase '" + phase.name + "' has invalid inlet pressure.");
    }
    const double velocity = momentumBalance ? phase.ports[0].velocity.value_or(0.0)
                                            : phase.columnEntranceVelocity.value_or(inputReader.columnEntranceVelocity);
    subStages.emplace_back(
        SwingAdsorptionSubStage{phase.name,
                                temperature,
                                pressure,
                                phase.duration,
                                phase.gasPhaseMolFractions,
                                velocity,
                                Column::FlowDirection(phase.flowDirection.value_or(inputReader.flowDirection)),
                                phase.productEndGasPhaseMolFractions,
                                phase.closedEnd,
                                Column::PressureRamp::Shape(phase.pressureRampShape),
                                phase.pressureRampTimeConstant,
                                phase.feedFromPhase,
                                phase.productEndFromPhase,
                                phase.stopCriteria,
                                phase.stopOnAny,
                                phase.minimumDuration,
                                phase.errorOnTimeout,
                                phase.ports,
                                phase.jacketTemperature.value_or(temperature),
                                phase.jacketHeatTransfer.value_or(inputReader.heatTransferWallExternal),
                                phase.coolantMassFlowRate.value_or(inputReader.coolantMassFlowRate),
                                phase.connections});
  }
}

void SwingAdsorption::run()
{
  // The columns: one, or those of Columns, each with its output files in a directory of its name.
  std::vector<Breakthrough*> beds{&breakthrough};
  for (Breakthrough& other : otherColumns) beds.push_back(&other);
  const size_t numberOfBeds = beds.size();
  const size_t n = breakthrough.column.numberOfComponents;
  const size_t lastNode = breakthrough.column.numberOfGridPoints;
  const double crossSection = breakthrough.column.crossSection();
  if (!columnNames.empty() && !(crossSection > 0.0))
  {
    throw std::runtime_error(
        "Error: Columns needs the cross-section of the column, through which the connections "
        "carry their gas: give the Geometry its InternalDiameter\n");
  }

  // phases.data is continued with the other output files (Breakthrough::writeOutput), which the first output starts.
  std::vector<std::ofstream> phaseStreams;
  for (size_t bed = 0; bed < numberOfBeds; ++bed)
  {
    if (!columnNames.empty())
    {
      std::filesystem::create_directories(columnNames[bed]);
      beds[bed]->outputPrefix = columnNames[bed] + "/";
    }
    const std::string fileName = beds[bed]->outputPrefix + "phases.data";
    const bool appendPhases = breakthrough.appendOutput && std::filesystem::exists(fileName);
    std::ofstream& phaseStream =
        phaseStreams.emplace_back(fileName, appendPhases ? std::ios_base::app : std::ios_base::trunc);
    if (appendPhases) continue;
    std::print(phaseStream,
               "# column 1: Phase start time [s]\n# column 2: Phase end time [s]\n"
               "# column 3: Pressure at z = 0 at the end of the phase [Pa]\n"
               "# column 4: Pressure at z = L at the end of the phase [Pa]\n"
               "# column 5: Max |R T sum(c) / P - 1| at the end of the phase [-]\n");
    for (size_t comp = 0; comp < n; ++comp)
    {
      std::print(phaseStream,
                 "# columns {} to {}: {}: gas and adsorbed at the phase start, gas and adsorbed at its end [mol/m^2]\n",
                 6 + 4 * comp, 9 + 4 * comp, breakthrough.column.components[comp].name);
    }
    if (breakthrough.column.streamAmountsEnabled)
    {
      std::print(phaseStream,
                 "# columns {} to {}: amount that entered at z = 0, left at z = 0, entered at z = L, and left at z = L "
                 "in the phase, each for every component in turn [mol/m^2]\n"
                 "# column {}: Heat that entered through the outer wall in the phase [J/m^2]\n"
                 "# column {}: Work that brings the gas that left below ExhaustPressure to it [J/m^2]\n",
                 6 + 4 * n, 5 + 8 * n, 6 + 8 * n, 7 + 8 * n);
    }
  }
  phaseRecords.clear();

  // The gas that the connections of each phase carried, and the buffers at its end.
  std::ofstream connectionStream;
  std::ofstream bufferStream;
  if (!columnNames.empty())
  {
    connectionStream.open("connections.data");
    bufferStream.open("buffers.data");
    std::print(connectionStream,
               "# One line per phase and connection\n# column 1: Phase start time [s]\n"
               "# column 2: Phase end time [s]\n"
               "# columns 3 to {}: Amount of each component that passed from the first endpoint to the second in "
               "the phase [mol]\n# column {}: Work of the compressor in the phase [J]\n",
               2 + n, 3 + n);
    std::print(bufferStream,
               "# One line per phase and buffer\n# column 1: Phase end time [s]\n"
               "# column 2: Pressure at the end of the phase [Pa]\n"
               "# columns 3 to {}: Amount of each component at the end of the phase [mol]\n"
               "# columns {} to {}: Amount of each component that left above MaximumPressure in the phase [mol]\n",
               2 + n, 3 + n, 2 + 2 * n);
  }

  // Cyclic steady state: the end of each cycle against the end of the one before (see the class description).
  std::ofstream cycleStream;
  if (steadyStateTolerance > 0.0)
  {
    const bool appendCycles = breakthrough.appendOutput && std::filesystem::exists("cycles.data");
    cycleStream.open("cycles.data", appendCycles ? std::ios_base::app : std::ios_base::trunc);
    if (!appendCycles)
    {
      std::print(cycleStream,
                 "# column 1: Cycle\n# column 2: Cycle end time [s]\n"
                 "# columns 3 to 5: Largest relative change from the cycle before of a state value, of the amount of "
                 "a component in the column, and of a stream amount of a phase [-]\n"
                 "# column 6: Largest amount of a component that the streams and the column of the cycle do not "
                 "account for, relative to the largest stream amount [-]\n");
    }
  }
  std::vector<std::vector<double>> previousState(numberOfBeds);
  std::vector<double> previousAmount;
  size_t steadyCycles = 0;
  // The mean gas that left each column in each phase of the cycle the last time it ran, for phases that take it.
  std::vector<std::vector<std::vector<double>>> products(numberOfBeds,
                                                         std::vector<std::vector<double>>(phasesPerCycle));

  auto applyStageConditions = [&](Breakthrough& bed, const SwingAdsorptionSubStage& stage,
                                  const std::vector<std::vector<double>>& phaseProducts)
  {
    Column& column = bed.column;
    const bool reverse = stage.flowDirection == Column::FlowDirection::Reverse;
    using Boundary = Column::BoundaryCondition;
    const bool momentumBalance = column.boundaryCondition == Boundary::MomentumBalance;
    const std::vector<double>& feedFractions =
        stage.feedFromPhase.has_value() && !phaseProducts[*stage.feedFromPhase].empty()
            ? phaseProducts[*stage.feedFromPhase]
            : stage.moleFractions;
    const std::vector<double>& productEndFractions =
        stage.productEndFromPhase.has_value() && !phaseProducts[*stage.productEndFromPhase].empty()
            ? phaseProducts[*stage.productEndFromPhase]
            : stage.productEndMoleFractions;
    std::ranges::fill(column.streamAmount, 0.0);
    // The gas takes the pressure profile of the new phase at its composition, where a continuity closure would keep the
    // old amount: after a pressure step without a ramp, a temperature step without the energy balance, and when the end
    // where the pressure is set changes, which moves the profile by the pressure drop. Left alone, that shift would add
    // up over the cycles. FixedVelocity and a fixed pressure profile let the feed displace the gas,
    // InletVelocityOutletPressure does not follow the phase pressure, and MomentumBalance keeps the gas it holds.
    const bool ramped = stage.rampShape != Column::PressureRamp::Shape::None;
    const bool rescaleGas =
        column.fluidPhase == Column::FluidPhase::Gas && column.boundaryCondition != Boundary::FixedVelocity &&
        column.boundaryCondition != Boundary::FixedPressureInletVelocity &&
        column.boundaryCondition != Boundary::InletVelocityOutletPressure && !momentumBalance &&
        ((!ramped && stage.pressure != column.inletPressure) || stage.flowDirection != column.flowDirection ||
         (!column.energyBalance && stage.temperature != column.gasTemperature[1]));
    column.pressureRamp = {stage.rampShape, stage.rampTimeConstant,
                           bed.time,        column.totalPressure[reverse ? lastNode : 0],
                           stage.pressure,  bed.time};
    column.pressureRate = 0.0;
    column.relativePressureRate = 0.0;
    column.flowDirection = stage.flowDirection;
    column.closedEnd = stage.closedEnd;
    column.setTemperature(stage.temperature);
    column.externalTemperature = stage.jacketTemperature;
    column.heatTransferWallExternal = stage.jacketHeatTransfer;
    column.coolantMassFlow = stage.coolantMassFlow;
    column.inletPressure = stage.pressure;
    // The ports of MomentumBalance. A ramp starts from the pressure of the end cell, a port without a pressure holds
    // gas at that pressure, and a molar flux is the velocity of the port gas at the phase temperature. The gas of the
    // port at z = 0 is the feed below, at its pressure.
    for (size_t end = 0; momentumBalance && end < column.ports.size(); ++end)
    {
      const InputReader::Port& port = stage.ports[end];
      const size_t cell = end == 0 ? 1 : lastNode;
      const double cellPressure = column.totalPressure[cell];
      const double pressure = port.pressureRampShape != 0 ? cellPressure : port.pressure.value_or(cellPressure);
      const double fluxVelocity = port.molarFlux.value_or(0.0) * R * stage.temperature /
                                  (pressure * column.totalVoidFraction[end == 0 ? 0 : lastNode]);
      column.ports[end] = {Column::Port::Mode(port.mode),
                           pressure,
                           port.velocity.value_or(fluxVelocity),
                           {Column::PressureRamp::Shape(port.pressureRampShape), port.pressureRampTimeConstant,
                            bed.time, cellPressure, port.pressure.value_or(cellPressure), bed.time}};
      if (end == 0) column.inletPressure = pressure;
    }
    column.advancePressureRamp(bed.time);
    // Boundary conditions that derive the velocity from the pressures recompute it in precompute.
    column.columnEntranceVelocity = stage.velocity;

    // The feed enters from node 0 in a forward phase and from the product-end gas in a reverse phase. There node 0,
    // which gas that flows back in at z = 0 does not use, takes the state of the first cell, so that the output at
    // z = 0 starts from the gas that leaves there.
    std::vector<double> feed(column.concentration.begin(), column.concentration.begin() + static_cast<std::ptrdiff_t>(n));
    if (column.fluidPhase == decltype(column.fluidPhase)::Liquid)
    {
      for (size_t comp = 0; comp < n; ++comp) feed[comp] = column.components[comp].inletLiquidConcentration;
    }
    else if (column.inletPressure > 0.0)
    {
      const double inletConcentration = column.inletPressure / (R * std::max(1e-10, stage.temperature));
      for (size_t comp = 0; comp < n; ++comp)
      {
        const double feedFraction =
            feedFractions.empty() ? column.components[comp].initialGasMoleFraction : feedFractions[comp];
        feed[comp] = feedFraction * inletConcentration;
      }
    }
    if (reverse)
    {
      column.productEndConcentration = feed;
      column.productEndTemperature = stage.temperature;
      std::copy_n(column.concentration.begin() + static_cast<std::ptrdiff_t>(n), n, column.concentration.begin());
      if (!column.gasTemperature.empty()) column.gasTemperature[0] = column.gasTemperature[1];
    }
    else
    {
      column.productEndConcentration.clear();
      column.productEndTemperature = 0.0;
      std::copy(feed.begin(), feed.end(), column.concentration.begin());
      if (!column.gasTemperature.empty()) column.gasTemperature[0] = stage.temperature;
    }
    if (momentumBalance && column.ports[0].mode == Column::Port::Mode::Closed)
    {
      // No gas passes a closed z = 0: node 0 shows the gas of the first cell in the output, as in a reverse phase.
      std::copy_n(column.concentration.begin() + static_cast<std::ptrdiff_t>(n), n, column.concentration.begin());
      column.gasTemperature[0] = column.gasTemperature[1];
    }
    // Derived fields (velocity, pressure, equilibrium), then the rates and the sink they feed back into.
    if (rescaleGas) std::fill(column.stateDot.begin(), column.stateDot.end(), 0.0);  // the rates of the old state
    precompute(column, timings);
    // The second pass takes the pressure drop of the rescaled gas.
    for (size_t pass = 0; rescaleGas && pass < 2; ++pass)
    {
      for (size_t grid = reverse ? 0 : 1; grid <= lastNode; ++grid)
      {
        const double factor = column.totalPressure[grid] /
                              std::max(1e-30, R * column.gasTemperature[grid] * column.totalConcentration[grid]);
        for (size_t comp = 0; comp < n; ++comp) column.concentration[grid * n + comp] *= factor;
      }
      precompute(column, timings);
    }
    if (!reverse && !productEndFractions.empty())
    {
      // The gas that flows back in at z = L, at the pressure there (that of its port with MomentumBalance) and the
      // phase temperature.
      const double productEndTotal =
          (momentumBalance ? column.ports[1].pressure : column.totalPressure[column.numberOfGridPoints]) /
          (R * std::max(1e-10, stage.temperature));
      column.productEndConcentration.resize(n);
      for (size_t comp = 0; comp < n; ++comp)
      {
        column.productEndConcentration[comp] = productEndFractions[comp] * productEndTotal;
      }
      column.productEndTemperature = stage.temperature;
      precompute(column, timings);
    }
    computeDerivatives(column);
    precompute(column, timings);
  };

  // Several columns advance together in coupling steps (see the class description). The gas at an endpoint of a
  // connection is that of the end cell of a column, or of a buffer.
  struct EndpointGas
  {
    double pressure{0.0};      // In Pa.
    double temperature{0.0};   // In K.
    double capacity{0.0};      // Gas per pressure of the column or the buffer, in mol/Pa.
    double heatCapacity{0.0};  // In J/(mol K).
    double molarMass{0.0};     // In kg/mol.
    std::vector<double> fractions{};
  };
  std::vector<EndpointGas> gas(2 * numberOfBeds + buffers.size());
  std::vector<double> previousFlows;  // Those of the coupling step before, in mol/s; empty at the start of a phase.
  double previousStep = 0.0;
  // Per connection of the phase the amount of each component that passed from its first endpoint to its second in
  // mol, then the work of its compressor in J; per buffer what left above its maximum pressure in the phase, in mol.
  std::vector<std::vector<double>> connectionAmounts;
  std::vector<std::vector<double>> exported(buffers.size(), std::vector<double>(n, 0.0));

  const auto advanceCoupled = [&](double target, const SwingAdsorptionSubStage& stage, const StopResidual* stop)
  {
    using FlowModel = InputReader::Connection::FlowModel;
    const std::vector<InputReader::Connection>& connections = stage.connections;
    std::vector<double> flows(connections.size());    // From the first endpoint to the second, in mol/s.
    std::vector<double> evenOut(connections.size());  // Gas that evens out the pressures of the endpoints, in mol.
    std::vector<double> carried(connections.size());  // From the donor to the receiver in this step, in mol/s.
    std::vector<std::array<size_t, 2>> route(connections.size());  // The donor and the receiver endpoint.
    std::vector<std::vector<double>> passed(connections.size(), std::vector<double>(n));  // In mol.
    std::vector<std::vector<double>> states(numberOfBeds);
    std::vector<std::vector<double>> rates(numberOfBeds);
    while (breakthrough.time < target)
    {
      const double time = breakthrough.time;
      for (size_t endpoint = 0; endpoint < gas.size(); ++endpoint)
      {
        EndpointGas& point = gas[endpoint];
        if (endpoint < 2 * numberOfBeds)
        {
          const Column& column = beds[endpoint / 2]->column;
          const size_t cell = endpoint % 2 == 0 ? 1 : lastNode;
          const std::vector<double> inventory = column.gasInventory();
          point.fractions.assign(column.concentration.begin() + static_cast<std::ptrdiff_t>(cell * n),
                                 column.concentration.begin() + static_cast<std::ptrdiff_t>((cell + 1) * n));
          point.temperature = column.gasTemperature[cell];
          point.pressure = R * point.temperature * std::reduce(point.fractions.begin(), point.fractions.end());
          point.capacity = crossSection * std::reduce(inventory.begin(), inventory.end()) / point.pressure;
        }
        else
        {
          const InputReader::Buffer& buffer = buffers[endpoint - 2 * numberOfBeds];
          point.fractions = bufferAmounts[endpoint - 2 * numberOfBeds];
          point.temperature = buffer.temperature;
          point.capacity = buffer.volume / (R * buffer.temperature);
          point.pressure = std::reduce(point.fractions.begin(), point.fractions.end()) / point.capacity;
        }
        const double total = std::reduce(point.fractions.begin(), point.fractions.end());
        point.molarMass = 0.0;
        point.heatCapacity = 0.0;
        for (size_t comp = 0; comp < n; ++comp)
        {
          point.fractions[comp] /= total;
          point.molarMass += point.fractions[comp] * breakthrough.column.components[comp].molecularWeight;
          point.heatCapacity += point.fractions[comp] * breakthrough.column.components[comp].gasHeatCapacity;
        }
        if (breakthrough.column.gasHeatCapacity.empty())
        {
          point.heatCapacity = breakthrough.column.heatCapacityGas * point.molarMass;
        }
      }

      // The flow that each connection carries at these pressures, and the coupling step: a twentieth of the shortest
      // time in which such a flow evens out the pressures of its endpoints (a compressor: moves that much gas), from
      // the gas that the columns and buffers hold per pressure.
      // shortcut: no shorter than 1e-4 of the phase, where a flow is cut to what half evens out the pressures; an
      // implicit coupling is needed for valves that are much wider than that.
      double step = target - time;
      for (size_t index = 0; index < connections.size(); ++index)
      {
        const InputReader::Connection& connection = connections[index];
        const EndpointGas& from = gas[connection.endpoints[0]];
        const EndpointGas& to = gas[connection.endpoints[1]];
        const double difference = from.pressure - to.pressure;
        const EndpointGas& upstream = difference >= 0.0 ? from : to;
        const bool compressor = connection.flowModel == FlowModel::Compressor;
        double flow = connection.conductance;
        if (connection.flowModel == FlowModel::Linear) flow *= std::abs(difference);
        if (connection.flowModel == FlowModel::SquareRoot) flow *= std::sqrt(std::abs(difference));
        if ((compressor || connection.flowModel == FlowModel::Orifice) && !(upstream.heatCapacity > R))
        {
          throw std::runtime_error(std::format(
              "Error: connection '{}' needs a gas heat capacity above R per mole of gas: give every component its "
              "GasHeatCapacity in J/(mol K)\n",
              connection.name));
        }
        if (connection.flowModel == FlowModel::Orifice)
        {
          // Isentropic flow of an ideal gas through the area; below the critical pressure ratio it is choked.
          const double gamma = upstream.heatCapacity / (upstream.heatCapacity - R);
          const double ratio = std::max(std::min(from.pressure, to.pressure) / upstream.pressure,
                                        std::pow(2.0 / (gamma + 1.0), gamma / (gamma - 1.0)));
          flow *= upstream.pressure / std::sqrt(upstream.molarMass * R * upstream.temperature) *
                  std::sqrt(2.0 * gamma / (gamma - 1.0) *
                            (std::pow(ratio, 2.0 / gamma) - std::pow(ratio, (gamma + 1.0) / gamma)));
        }
        flow = std::min(connection.opening * flow, connection.maximumMolarFlow);
        flows[index] = connection.checkValve && difference < 0.0 ? 0.0 : compressor || difference >= 0.0 ? flow : -flow;
        evenOut[index] = (compressor ? std::min(from.pressure, to.pressure) : std::abs(difference)) /
                         (1.0 / from.capacity + 1.0 / to.capacity);
        if (flows[index] != 0.0) step = std::min(step, 0.05 * evenOut[index] / flow);
      }
      step = std::max(step, std::min(1.0e-4 * stage.duration, target - time));
      for (size_t index = 0; index < connections.size(); ++index)
      {
        // The flow at the middle of the step, from this one and the one before, while that keeps its direction.
        double flow = flows[index];
        if (!previousFlows.empty())
        {
          const double ahead = flow + 0.5 * step / previousStep * (flow - previousFlows[index]);
          if (ahead * flow > 0.0) flow = ahead;
        }
        carried[index] = std::min(std::abs(flow), 0.5 * evenOut[index] / step);
        route[index] = {connections[index].endpoints[flow >= 0.0 ? 0 : 1],
                        connections[index].endpoints[flow >= 0.0 ? 1 : 0]};
      }

      for (size_t bed = 0; bed < numberOfBeds; ++bed)
      {
        states[bed] = beds[bed]->column.state;
        rates[bed] = beds[bed]->column.stateDot;
      }
      // The step in which the stop criteria turn is repeated from its start, with its end bisected between the last
      // time before they turned and the first after it, as Breakthrough::advanceTo does for one column.
      double stepEnd = step >= target - time ? target : time + step;
      const double resolution = 1.0e-6 * (stepEnd - time);
      double low = time;
      double high = std::numeric_limits<double>::infinity();
      for (bool repeat = true; repeat;)
      {
        // A buffer gives gas of its composition at the start of the step.
        for (size_t index = 0; index < connections.size(); ++index)
        {
          for (size_t comp = 0; comp < n; ++comp)
          {
            passed[index][comp] = route[index][0] < 2 * numberOfBeds
                                      ? 0.0
                                      : carried[index] * (stepEnd - time) * gas[route[index][0]].fractions[comp];
          }
        }
        // A column takes its step once the columns that give it gas have taken theirs: what left those, by their
        // stream amounts, is what enters it, so that every component is conserved.
        std::vector<bool> done(numberOfBeds, false);
        for (size_t pass = 0; pass < numberOfBeds; ++pass)
        {
          for (size_t bed = 0; bed < numberOfBeds; ++bed)
          {
            bool ready = !done[bed];
            for (size_t index = 0; index < connections.size(); ++index)
            {
              const size_t donor = route[index][0];
              ready = ready && !(carried[index] > 0.0 && route[index][1] / 2 == bed && donor < 2 * numberOfBeds &&
                                 !done[donor / 2]);
            }
            if (!ready) continue;

            Column& column = beds[bed]->column;
            std::array<double, 2> leaving{};  // Through each end, in mol/s.
            for (size_t end = 0; end < 2; ++end)
            {
              const size_t endpoint = 2 * bed + end;
              bool connected = false;
              std::vector<double> entering(n, 0.0);  // In mol.
              double warmth = 0.0;                   // Their sum times the temperature of the gas.
              for (size_t index = 0; index < connections.size(); ++index)
              {
                connected = connected || route[index][0] == endpoint || route[index][1] == endpoint;
                if (route[index][0] == endpoint) leaving[end] += carried[index];
                for (size_t comp = 0; route[index][1] == endpoint && comp < n; ++comp)
                {
                  entering[comp] += passed[index][comp];
                  warmth += passed[index][comp] * gas[route[index][0]].temperature;
                }
              }
              if (!connected) continue;
              const double entered = std::reduce(entering.begin(), entering.end());
              if (entered > 0.0 && leaving[end] > 0.0)
              {
                // shortcut: one end takes gas or gives it; a tee that does both at once needs a junction of its own.
                throw std::runtime_error(std::format(
                    "Error: in swing phase '{}' gas enters and leaves column {} at z = {} at once through its "
                    "connections; give the two flows a buffer between them\n",
                    stage.name, columnNames[bed], end == 0 ? "0" : "L"));
              }
              column.ports[end].mode = Column::Port::Mode::Flux;
              column.ports[end].flux = (entered / (stepEnd - time) - leaving[end]) / crossSection;
              if (!(entered > 0.0)) continue;
              // The gas that enters, at the pressure of the end cell. An ideal gas keeps its temperature through a
              // valve; that of a compressor is cooled back to it.
              const double temperature = column.energyBalance ? warmth / entered : gas[endpoint].temperature;
              const double total = gas[endpoint].pressure / (R * temperature);
              if (end == 1) column.productEndConcentration.resize(n);
              for (size_t comp = 0; comp < n; ++comp)
              {
                (end == 0 ? column.concentration[comp] : column.productEndConcentration[comp]) =
                    entering[comp] / entered * total;
              }
              if (end == 0 && column.energyBalance) column.gasTemperature[0] = temperature;
              if (end == 1) column.productEndTemperature = temperature;
            }
            precompute(column, timings);
            const std::vector<double> before(column.streamAmount.begin(), column.streamAmount.end());
            beds[bed]->advanceTo(stepEnd, timings);
            for (size_t index = 0; index < connections.size(); ++index)
            {
              const size_t end = route[index][0] % 2;
              for (size_t comp = 0; route[index][0] / 2 == bed && carried[index] > 0.0 && comp < n; ++comp)
              {
                passed[index][comp] = carried[index] / leaving[end] * crossSection *
                                      (before[end * n + comp] - column.streamAmount[end * n + comp]);
              }
            }
            done[bed] = true;
          }
        }
        if (std::ranges::find(done, false) != done.end())
        {
          // shortcut: no ring of columns that each give gas to the next; that needs their steps solved together.
          throw std::runtime_error(std::format(
              "Error: in swing phase '{}' the connections carry gas around a ring of columns, which is not supported\n",
              stage.name));
        }

        if (stop != nullptr) ((*stop)(stepEnd, breakthrough.column) < 0.0 ? low : high) = stepEnd;
        if (!std::isfinite(high)) break;
        repeat = high - low > resolution || stepEnd != high;
        for (size_t bed = 0; repeat && bed < numberOfBeds; ++bed)
        {
          Column& column = beds[bed]->column;
          std::ranges::copy(states[bed], column.state.begin());
          std::ranges::copy(rates[bed], column.stateDot.begin());
          beds[bed]->time = time;
          precompute(column, timings);
        }
        if (repeat) stepEnd = high - low > resolution ? 0.5 * (low + high) : high;
      }

      for (size_t index = 0; index < connections.size(); ++index)
      {
        const double forward = route[index][0] == connections[index].endpoints[0] ? 1.0 : -1.0;
        for (size_t comp = 0; comp < n; ++comp)
        {
          for (size_t side = 0; side < 2; ++side)
          {
            if (route[index][side] < 2 * numberOfBeds) continue;
            bufferAmounts[route[index][side] - 2 * numberOfBeds][comp] +=
                (side == 0 ? -1.0 : 1.0) * passed[index][comp];
          }
          connectionAmounts[index][comp] += forward * passed[index][comp];
        }
        // The work of a compressor: the isentropic compression of its gas to the pressure behind it, over its
        // efficiency; the pressures are those of the middle of the step, by the gas that the endpoints hold per
        // pressure.
        const EndpointGas& from = gas[connections[index].endpoints[0]];
        const EndpointGas& to = gas[connections[index].endpoints[1]];
        const double moved = std::reduce(passed[index].begin(), passed[index].end());
        const double ratio = (to.pressure + 0.5 * moved / to.capacity) / (from.pressure - 0.5 * moved / from.capacity);
        if (connections[index].flowModel == FlowModel::Compressor && ratio > 1.0)
        {
          connectionAmounts[index][n] += moved * from.heatCapacity * from.temperature *
                                         std::expm1(R / from.heatCapacity * std::log(ratio)) /
                                         connections[index].efficiency;
        }
      }
      for (size_t index = 0; index < buffers.size(); ++index)
      {
        // Gas above the maximum pressure leaves the buffer: its export.
        const double amount = std::reduce(bufferAmounts[index].begin(), bufferAmounts[index].end());
        const double most = buffers[index].maximumPressure * buffers[index].volume / (R * buffers[index].temperature);
        for (size_t comp = 0; buffers[index].maximumPressure > 0.0 && amount > most && comp < n; ++comp)
        {
          exported[index][comp] += bufferAmounts[index][comp] * (1.0 - most / amount);
          bufferAmounts[index][comp] *= most / amount;
        }
      }
      previousFlows = flows;
      previousStep = stepEnd - time;
      if (std::isfinite(high)) return true;
    }
    return false;
  };

  const double writeInterval = breakthrough.writeInterval;
  try
  {
    // An output that falls on the start of a phase is written once that phase's conditions are applied, as is the
    // initial state, unless the files hold it already (a continued run); the state at the end of the last phase is
    // always written.
    bool outputPending = !breakthrough.appendOutput;
    double phaseEndTime = breakthrough.time;
    for (size_t stageIndex = 0; stageIndex < subStages.size() / numberOfBeds; ++stageIndex)
    {
      // The sub-stage of the first column holds what the columns share: the times, the stop criteria, the connections.
      const SwingAdsorptionSubStage& stage = subStages[stageIndex * numberOfBeds];
      std::vector<SwingAdsorptionPhaseRecord> records;
      for (size_t bed = 0; bed < numberOfBeds; ++bed)
      {
        Breakthrough& simulation = *beds[bed];
        Column& column = simulation.column;
        const SwingAdsorptionSubStage& bedStage = subStages[stageIndex * numberOfBeds + bed];
        try
        {
          applyStageConditions(simulation, bedStage, products[bed]);
        }
        catch (const std::runtime_error& error)
        {
          throw std::runtime_error(
              std::format("{}at the start of swing phase '{}'; a large pressure step can need a "
                          "PressureRamp or a ClosedEnd\n",
                          error.what(), bedStage.name));
        }
        if (simulation.integrationScheme == BreakthroughIntegrationScheme::CVODE)
        {
          // Restart the multistep history and re-derive the 4 dz/v step cap from this phase's velocity.
          simulation.cvode.reinitialize(column);
        }
        // The switched boundary conditions start a new transient, which ROS34PW2 starts with a small step.
        if (simulation.integrationScheme == BreakthroughIntegrationScheme::ROS34PW2) simulation.ros34pw2.reinitialize();
        const bool reverse = bedStage.flowDirection == Column::FlowDirection::Reverse;
        std::print(
            "\nSwing phase {}: {} (T={} K, inlet pressure={} Pa, inlet velocity={} m/s, duration={} s{}{}{}{}{})\n",
            stageIndex + 1, bedStage.name, bedStage.temperature, bedStage.pressure, column.columnEntranceVelocity,
            bedStage.duration, reverse && !bedStage.closedEnd ? ", reverse flow" : "",
            !bedStage.closedEnd ? ""
            : reverse           ? ", closed at z = 0"
                                : ", closed at z = L",
            bedStage.rampShape == Column::PressureRamp::Shape::None ? "" : ", pressure ramp",
            column.boundaryCondition != Column::BoundaryCondition::MomentumBalance
                ? ""
                : std::format(", ports {} at z = 0 and {} at z = L", portModeNames[bedStage.ports[0].mode],
                              portModeNames[bedStage.ports[1].mode]),
            columnNames.empty() ? "" : ", column " + columnNames[bed]);
        records.push_back({.name = bedStage.name,
                           .startTime = simulation.time,
                           .gasAtStart = column.gasInventory(),
                           .adsorbedAtStart = column.adsorbedInventory()});
      }
      connectionAmounts.assign(stage.connections.size(), std::vector<double>(n + 1, 0.0));
      for (std::vector<double>& amounts : exported) std::ranges::fill(amounts, 0.0);
      previousFlows.clear();
      if (outputPending)
      {
        for (Breakthrough* bed : beds) bed->writeOutput();
        outputPending = false;
      }

      // The stop criteria as one residual: negative until, after the minimum duration, all (or any) of them hold.
      const double phaseStartTime = breakthrough.time;
      const StopResidual stopResidual = [&](double time, const ColumnStateViews& views)
      {
        double residual =
            stage.stopOnAny ? -std::numeric_limits<double>::infinity() : std::numeric_limits<double>::infinity();
        for (const InputReader::StopCriterion& criterion : stage.stopCriteria)
        {
          // One column: the state that its integrator passes; several: that of the column of the criterion.
          const Column& column = beds[criterion.column]->column;
          const double excess =
              stopMetric(column, columnNames.empty() ? views : static_cast<const ColumnStateViews&>(column),
                         criterion) -
              criterion.value;
          residual = stage.stopOnAny ? std::max(residual, criterion.atLeast ? excess : -excess)
                                     : std::min(residual, criterion.atLeast ? excess : -excess);
        }
        return std::min(residual, time - phaseStartTime - stage.minimumDuration);
      };
      const StopResidual* stop = stage.stopCriteria.empty() ? nullptr : &stopResidual;

      // Integrate to each output time inside the phase and then to the phase end, where the conditions switch, or to
      // the time at which the stop criteria turn true.
      phaseEndTime = phaseStartTime + stage.duration;
      const double tolerance = 1.0e-9 * writeInterval;
      bool stopped = stop != nullptr && stopResidual(phaseStartTime, breakthrough.column) >= 0.0;
      while (!stopped && breakthrough.time < phaseEndTime - tolerance)
      {
        // Output times are computed from their index, so that they do not drift over many phases.
        const double nextOutputTime = static_cast<double>(breakthrough.numberOfOutputs + 1) * writeInterval;
        const bool outputBeforePhaseEnd = nextOutputTime < phaseEndTime - tolerance;
        const bool outputAtPhaseEnd = !outputBeforePhaseEnd && nextOutputTime <= phaseEndTime + tolerance;
        const double target = outputBeforePhaseEnd ? nextOutputTime : phaseEndTime;
        stopped =
            columnNames.empty() ? breakthrough.advanceTo(target, timings, stop) : advanceCoupled(target, stage, stop);
        if (stopped) break;
        for (Breakthrough* bed : beds)
        {
          if (outputBeforePhaseEnd || outputAtPhaseEnd) ++bed->numberOfOutputs;
          if (outputBeforePhaseEnd) bed->writeOutput();
        }
        outputPending = outputAtPhaseEnd;
      }
      if (stop != nullptr)
      {
        if (!stopped && stage.errorOnTimeout)
        {
          throw std::runtime_error(
              std::format("Error: swing phase '{}' reached its Duration of {} s before its "
                          "StopCriterion held (TimeoutPolicy Error)\n",
                          stage.name, stage.duration));
        }
        std::print("    Stop criterion {} at t = {:6.5f} [s]\n", stopped ? "met" : "not met within the Duration",
                   breakthrough.time);
      }

      for (size_t bed = 0; bed < numberOfBeds; ++bed)
      {
        const Column& column = beds[bed]->column;
        SwingAdsorptionPhaseRecord& record = records[bed];
        record.endTime = breakthrough.time;
        record.pressureAtZero = column.totalPressure[0];
        record.pressureAtLength = column.totalPressure[lastNode];
        record.gasDeviation = column.idealGasDeviation();
        record.gasAtEnd = column.gasInventory();
        record.adsorbedAtEnd = column.adsorbedInventory();
        for (size_t index = 0; index < 4 * n && column.streamAmountsEnabled; ++index)
        {
          // Entered and left at z = 0, then at z = L, from the net amount that entered at each end.
          const double entered = column.streamAmount[(index / (2 * n)) * n + index % n];
          record.streamAmounts.push_back(std::max(0.0, (index / n) % 2 == 0 ? entered : -entered));
        }
        std::ofstream& phaseStream = phaseStreams[bed];
        std::print(phaseStream, "{} {} {} {} {}", record.startTime, record.endTime, record.pressureAtZero,
                   record.pressureAtLength, record.gasDeviation);
        for (size_t comp = 0; comp < n; ++comp)
        {
          std::print(phaseStream, " {} {} {} {}", record.gasAtStart[comp], record.adsorbedAtStart[comp],
                     record.gasAtEnd[comp], record.adsorbedAtEnd[comp]);
        }
        for (const double amount : record.streamAmounts) std::print(phaseStream, " {}", amount);
        if (column.streamAmountsEnabled)
        {
          record.wallHeat = column.streamAmount[2 * n];
          record.exhaustWork = column.streamAmount[2 * n + 1];
          std::print(phaseStream, " {} {}", record.wallHeat, record.exhaustWork);
        }
        std::print(phaseStream, "  # {}\n", record.name);
        phaseStream.flush();
        if (!phaseStream) throw std::runtime_error("Error: cannot write the output file 'phases.data'");

        // The mean gas that left at either end, for the phases that take it as their gas.
        std::vector<double> left(n, 0.0);
        for (size_t comp = 0; comp < n && column.streamAmountsEnabled; ++comp)
        {
          left[comp] = record.streamAmounts[n + comp] + record.streamAmounts[3 * n + comp];
        }
        if (const double total = std::reduce(left.begin(), left.end()); total > 0.0)
        {
          for (double& amount : left) amount /= total;
          products[bed][stageIndex % phasesPerCycle] = left;
        }
        phaseRecords.push_back(std::move(record));
      }
      for (size_t index = 0; index < stage.connections.size(); ++index)
      {
        std::print(connectionStream, "{} {}", phaseStartTime, breakthrough.time);
        for (const double amount : connectionAmounts[index]) std::print(connectionStream, " {}", amount);
        std::print(connectionStream, "  # {}: {}\n", stage.name, stage.connections[index].name);
      }
      for (size_t index = 0; index < buffers.size(); ++index)
      {
        std::print(bufferStream, "{} {}", breakthrough.time,
                   R * buffers[index].temperature *
                       std::reduce(bufferAmounts[index].begin(), bufferAmounts[index].end()) / buffers[index].volume);
        for (const double amount : bufferAmounts[index]) std::print(bufferStream, " {}", amount);
        for (const double amount : exported[index]) std::print(bufferStream, " {}", amount);
        std::print(bufferStream, "  # {}: {}\n", stage.name, buffers[index].name);
      }
      connectionStream.flush();
      bufferStream.flush();
      if (!columnNames.empty() && (!connectionStream || !bufferStream))
      {
        throw std::runtime_error("Error: cannot write the output file 'connections.data' or 'buffers.data'");
      }

      if (steadyStateTolerance <= 0.0 || (stageIndex + 1) % phasesPerCycle != 0) continue;
      // The end of a cycle: its changes from the end of the cycle before, each relative to the size of the value plus
      // the largest of its kind.
      const auto relativeChange = [](double value, double previous, double largest)
      {
        return value == previous
                   ? 0.0
                   : std::abs(value - previous) / (largest + std::max(std::abs(value), std::abs(previous)));
      };
      const size_t cycle = stageIndex / phasesPerCycle;
      // The records of the cycle, one per phase and column. The amounts are those of every column, then of every
      // buffer per unit column cross-section; what the streams do not account for is that of each column.
      const size_t recordsPerCycle = phasesPerCycle * numberOfBeds;
      const std::span<const SwingAdsorptionPhaseRecord> cycleRecords(
          phaseRecords.end() - static_cast<std::ptrdiff_t>(recordsPerCycle), recordsPerCycle);
      std::vector<double> amount(numberOfBeds * n);
      std::vector<double> unaccounted(numberOfBeds * n, 0.0);
      double largestStream = 0.0;
      for (size_t index = 0; index < recordsPerCycle; ++index)
      {
        const SwingAdsorptionPhaseRecord& phase = cycleRecords[index];
        const size_t first = (index % numberOfBeds) * n;
        largestStream =
            std::max(largestStream, *std::max_element(phase.streamAmounts.begin(), phase.streamAmounts.end()));
        for (size_t comp = 0; comp < n; ++comp)
        {
          unaccounted[first + comp] +=
              phase.streamAmounts[comp] - phase.streamAmounts[n + comp] + phase.streamAmounts[2 * n + comp] -
              phase.streamAmounts[3 * n + comp] -
              (phase.gasAtEnd[comp] + phase.adsorbedAtEnd[comp] - phase.gasAtStart[comp] - phase.adsorbedAtStart[comp]);
          amount[first + comp] = phase.gasAtEnd[comp] + phase.adsorbedAtEnd[comp];
        }
      }
      for (double& value : unaccounted) value = std::abs(value) / std::max(largestStream, 1e-300);
      for (const std::vector<double>& amounts : bufferAmounts)
      {
        for (const double value : amounts) amount.push_back(value / crossSection);
      }
      double stateChange = std::numeric_limits<double>::infinity();
      double amountChange = stateChange;
      double streamChange = stateChange;
      if (cycle > 0)
      {
        // The state field by field: the concentrations, each kind of loading, and the three temperatures.
        stateChange = 0.0;
        for (size_t bed = 0; bed < numberOfBeds; ++bed)
        {
          const Column& column = beds[bed]->column;
          const ColumnStateLayout layout = column.stateLayout();
          for (size_t field = 0; field < layout.componentBlockCount() + 3; ++field)
          {
            const bool thermal = field >= layout.componentBlockCount();
            const size_t first =
                thermal ? layout.temperatureOffset() + (field - layout.componentBlockCount()) * layout.scalarBlockSize()
                        : field * layout.componentBlockSize();
            const std::span<const double> values(column.state.data() + first,
                                                 thermal ? layout.scalarBlockSize() : layout.componentBlockSize());
            const double largest = std::ranges::max(values, {}, [](double value) { return std::abs(value); });
            for (size_t index = 0; index < values.size(); ++index)
            {
              stateChange = std::max(
                  stateChange, relativeChange(values[index], previousState[bed][first + index], std::abs(largest)));
            }
          }
        }
        amountChange = 0.0;
        const double largestAmount = std::max(*std::max_element(amount.begin(), amount.end()),
                                              *std::max_element(previousAmount.begin(), previousAmount.end()));
        for (size_t comp = 0; comp < amount.size(); ++comp)
        {
          amountChange = std::max(amountChange, relativeChange(amount[comp], previousAmount[comp], largestAmount));
        }
        streamChange = 0.0;
        for (size_t phase = 0; phase < recordsPerCycle; ++phase)
        {
          const std::vector<double>& before = (cycleRecords.data() - recordsPerCycle + phase)->streamAmounts;
          for (size_t index = 0; index < before.size(); ++index)
          {
            streamChange = std::max(
                streamChange, relativeChange(cycleRecords[phase].streamAmounts[index], before[index], largestStream));
          }
        }
      }
      for (size_t bed = 0; bed < numberOfBeds; ++bed) previousState[bed] = beds[bed]->column.state;
      previousAmount = amount;
      std::print(cycleStream, "{} {} {} {} {} {}\n", cycle + 1, breakthrough.time, stateChange, amountChange,
                 streamChange, *std::max_element(unaccounted.begin(), unaccounted.end()));
      cycleStream.flush();
      if (!cycleStream) throw std::runtime_error("Error: cannot write the output file 'cycles.data'");
      steadyCycles = std::max({stateChange, amountChange, streamChange}) <= steadyStateTolerance ? steadyCycles + 1 : 0;
      if (steadyCycles >= steadyStateCycles)
      {
        std::print("\nCyclic steady state after {} cycles (tolerance {})\n", cycle + 1, steadyStateTolerance);
        break;
      }
    }
    for (Breakthrough* bed : beds) bed->writeOutput();
  }
  catch (const std::exception&)
  {
    for (const Breakthrough* bed : beds) bed->column.writeJSON(bed->outputPrefix + "failed_state.json", bed->time);
    throw;
  }

  std::print("Final time: {:6.5f} [s]\n\n", breakthrough.time);
  timings.print();
  for (Breakthrough* bed : beds)
  {
    if (bed->integrationScheme == BreakthroughIntegrationScheme::CVODE) bed->cvode.printStatistics();
    if (bed->integrationScheme == BreakthroughIntegrationScheme::ROS34PW2) bed->ros34pw2.printStatistics();
  }
}

void SwingAdsorption::print() const
{
  std::print("{}", repr());
}

std::string SwingAdsorption::repr() const
{
  std::string result = breakthrough.repr();
  result += std::format(
      "Swing adsorption phases\n"
      "=======================================================\n");
  for (size_t i = 0; i < subStages.size(); ++i)
  {
    const SwingAdsorptionSubStage& stage = subStages[i];
    const bool reverse = stage.flowDirection == Column::FlowDirection::Reverse;
    const size_t numberOfBeds = otherColumns.size() + 1;
    result +=
        std::format("{}: {}{}  T={} [K], inlet pressure={} [Pa], inlet velocity={} [m/s], duration={} [s]{}{}{}{}{}\n",
                    i / numberOfBeds + 1, stage.name,
                    columnNames.empty() ? "" : ", column " + columnNames[i % numberOfBeds], stage.temperature,
                    stage.pressure, stage.velocity, stage.duration, reverse ? ", reverse flow (feed at z = L)" : "",
                    !stage.closedEnd ? ""
                    : reverse        ? ", closed at z = 0"
                                     : ", closed at z = L",
                    stage.rampShape == Column::PressureRamp::Shape::None
                        ? ""
                        : std::format(", pressure ramp over {} [s]", stage.rampTimeConstant),
                    breakthrough.column.boundaryCondition != Column::BoundaryCondition::MomentumBalance
                        ? ""
                        : std::format(", ports {} at z = 0 and {} at z = L", portModeNames[stage.ports[0].mode],
                                      portModeNames[stage.ports[1].mode]),
                    stage.stopCriteria.empty() ? "" : ", ends early on its stop criterion");
    for (const auto& [label, fractions] : {std::pair{"feed", &stage.moleFractions},
                                            std::pair{"product-end gas", &stage.productEndMoleFractions}})
    {
      if (fractions->empty()) continue;
      result += std::format("    {} mole fractions:", label);
      for (size_t comp = 0; comp < fractions->size(); ++comp)
      {
        result += std::format(" {}={}", breakthrough.column.components[comp].name, (*fractions)[comp]);
      }
      result += "\n";
    }
    for (const InputReader::Connection& connection : stage.connections)
    {
      result += std::format("    connection: {}\n", connection.name);
    }
  }
  result += "\n\n";
  return result;
}

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "column.h"
#include "inputreader.h"
#include "integrators/rk3.h"
#include "json.h"
#include "swing_adsorption.h"
#include "temporary_directory.h"
#include "timing.h"
#include "utils.h"

namespace
{
// Experiment T15 of Marx et al.: CO2/N2/H2 on activated carbon with the energy balance, in a helium-filled column.
nlohmann::json t15Input(size_t gridPoints)
{
  std::ifstream input(std::filesystem::path{RUPTURA_SOURCE_DIR} / "marx2/ternary/T15/simulation.json");
  nlohmann::json data = nlohmann::json::parse(input);
  data["NumberOfGridPoints"] = gridPoints;
  data["BreakthroughIntegrator"] = "CVODE";
  data["SimulationType"] = "SwingAdsorption";
  data.erase("SimulationTime");
  return data;
}

InputReader readInput(const nlohmann::json& data)
{
  const std::filesystem::path path = testTemporaryDirectory() / "ruptura_pressure_swing.json";
  std::ofstream(path) << data;
  return InputReader(path.string());
}

nlohmann::json closedPhase(const std::string& closedEnd, double pressure, double duration, const std::string& shape,
                           double timeConstant)
{
  return {{"Name", closedEnd + " closed"},
          {"ClosedEnd", closedEnd},
          {"InletPressure", pressure},
          {"Duration", duration},
          {"PressureRamp", {{"Shape", shape}, {"TimeConstant", timeConstant}}}};
}

// Whitespace-separated values of one column of a data file, over all rows.
std::vector<double> readColumn(const std::string& fileName, size_t column)
{
  std::ifstream file(fileName);
  std::vector<double> values;
  for (std::string line; std::getline(file, line);)
  {
    if (line.empty() || line.front() == '#') continue;
    std::istringstream row(line);
    double value = 0.0;
    for (size_t index = 0; index <= column; ++index) row >> value;
    values.push_back(value);
  }
  return values;
}

// Amount in the gas and on the solid after a phase minus before it, per component.
std::vector<double> amountChange(const SwingAdsorptionPhaseRecord& record)
{
  std::vector<double> change(record.gasAtEnd.size());
  for (size_t comp = 0; comp < change.size(); ++comp)
  {
    change[comp] = record.gasAtEnd[comp] + record.adsorbedAtEnd[comp] - record.gasAtStart[comp] -
                   record.adsorbedAtStart[comp];
  }
  return change;
}
}  // namespace

TEST(PressureRamp, FollowsItsShapeAndScalesTheProductEndGas)
{
  SwingAdsorption swing(readInput([]
                                  {
                                    nlohmann::json data = t15Input(4);
                                    data["SwingAdsorptionPhases"] = {{{"Duration", 1.0}}};
                                    return data;
                                  }()));
  Column& column = swing.breakthrough.column;
  column.productEndConcentration = {2.0};

  column.inletPressure = 1.0e5;
  column.pressureRamp = {Column::PressureRamp::Shape::Exponential, 20.0, 100.0, 1.0e5, 1.0e6, 100.0};
  column.advancePressureRamp(120.0);
  const double expected = 1.0e6 - 9.0e5 * std::exp(-1.0);
  EXPECT_DOUBLE_EQ(column.inletPressure, expected);
  EXPECT_DOUBLE_EQ(column.pressureRate, 9.0e5 * std::exp(-1.0) / 20.0);
  EXPECT_DOUBLE_EQ(column.relativePressureRate, column.pressureRate / expected);
  EXPECT_DOUBLE_EQ(column.productEndConcentration[0], 2.0 * expected / 1.0e5);
  EXPECT_EQ(column.pressureRamp.time, 120.0);

  column.pressureRamp = {Column::PressureRamp::Shape::Linear, 10.0, 0.0, 1.0e6, 2.0e5, 0.0};
  column.advancePressureRamp(5.0);
  EXPECT_DOUBLE_EQ(column.inletPressure, 6.0e5);
  EXPECT_DOUBLE_EQ(column.pressureRate, -8.0e4);
  column.advancePressureRamp(12.0);
  EXPECT_DOUBLE_EQ(column.inletPressure, 2.0e5);
  EXPECT_EQ(column.pressureRate, 0.0);
  EXPECT_EQ(column.relativePressureRate, 0.0);
}

TEST(ClosedEnd, InertColumnFollowsTheAnalyticVelocity)
{
  // Inert helium: v(z) = (L - z) (dP/dt) / P with z = L closed, and -z (dP/dt) / P with z = 0 closed, for pressurization
  // and blowdown, whatever the scheme; the pressure and the total concentration follow the ramp.
  nlohmann::json data = t15Input(16);
  data["energyBalance"] = false;
  data["Components"][0]["GasPhaseMolFraction"] = 1.0;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["GasPhaseMolFraction"] = 0.0;
  }
  data["WriteInterval"] = 1000.0;
  for (const std::string scheme : {"WENO3", "Upwind"})
  {
    for (const std::string closedEnd : {"ProductEnd", "FeedEnd"})
    {
      for (const std::string shape : {"Linear", "Exponential"})
      {
        for (const double target : {1.0e5, 2.0e6})
        {
          SCOPED_TRACE(scheme + " " + closedEnd + " " + shape + " " + std::to_string(target));
          data["ConvectionScheme"] = scheme;
          data["SwingAdsorptionPhases"] = nlohmann::json::array({closedPhase(closedEnd, target, 25.0, shape, 50.0)});
          SwingAdsorption swing(readInput(data));
          swing.run();

          const Column& column = swing.breakthrough.column;
          // The inlet pressure starts from the pressure of its end in the initial state, lower at z = L.
          const double start = column.pressureRamp.startPressure;
          EXPECT_NEAR(start, 1.0e6, 100.0);
          const double pressure = shape == "Linear"
                                      ? start + (target - start) * 0.5
                                      : target + (start - target) * std::exp(-0.5);
          const double rate = shape == "Linear" ? (target - start) / 50.0
                                                : (target - start) * std::exp(-0.5) / 50.0;
          EXPECT_NEAR(column.inletPressure, pressure, 1.0e-9 * pressure);
          const size_t gridPoints = column.numberOfGridPoints;
          const double velocityScale = column.columnLength * std::abs(rate) / pressure;
          for (size_t face = closedEnd == "FeedEnd" ? 1 : 0; face <= gridPoints; ++face)
          {
            const double z = column.columnDistances[face];
            const double analytic = (closedEnd == "ProductEnd" ? column.columnLength - z : -z) * rate / pressure;
            EXPECT_NEAR(column.interstitialGasVelocity[face], analytic, 1.0e-4 * velocityScale) << face;
          }
          EXPECT_EQ(column.interstitialGasVelocity[closedEnd == "ProductEnd" ? gridPoints : 0], 0.0);
          EXPECT_LT(column.idealGasDeviation(), 1.0e-4);
          EXPECT_NEAR(column.pressureRate, rate, 1.0e-9 * std::abs(rate));
        }
      }
    }
  }
}

TEST(ClosedEnd, ClosureCarriesTheAccumulationThePressureWorkAndConservesEverySpecies)
{
  // A front and a temperature profile with the energy balance: in every cell sum_j dc_j/dt = C ((dP/dt)/P - (dT/dt)/T),
  // and the species balances over the column close on the flux through the open end only.
  for (const std::string scheme : {"WENO3", "Upwind"})
  {
    for (const bool reverse : {false, true})
    {
      for (const double target : {1.0e5, 3.0e6})
      {
        SCOPED_TRACE(scheme + (reverse ? " closed at z = 0 " : " closed at z = L ") + std::to_string(target));
        nlohmann::json data = t15Input(16);
        data["ConvectionScheme"] = scheme;
        data["heatCapacityGas"] = 5193.0;  // of helium, the carrier here; T15 gives the value of the feed
        data["SwingAdsorptionPhases"] = {{{"Duration", 1.0}}};
        SwingAdsorption swing(readInput(data));
        Column& column = swing.breakthrough.column;
        const size_t n = column.numberOfComponents;
        const size_t gridPoints = column.numberOfGridPoints;
        const double feedTemperature = column.gasTemperature[0];
        const std::vector<double> feed(column.concentration.begin(),
                                       column.concentration.begin() + static_cast<std::ptrdiff_t>(n));
        column.flowDirection = reverse ? Column::FlowDirection::Reverse : Column::FlowDirection::Forward;
        column.closedEnd = true;
        column.pressureRamp = {Column::PressureRamp::Shape::Exponential, 20.0, 0.0, 1.0e6, target, 0.0};
        column.advancePressureRamp(0.0);
        if (reverse)
        {
          column.productEndConcentration = feed;
          column.productEndTemperature = feedTemperature;
        }
        for (size_t grid = 1; grid <= gridPoints; ++grid)
        {
          const double z = column.columnDistances[grid] / column.columnLength;
          const double front = 0.5 + 0.5 * std::tanh(15.0 * (0.4 - (reverse ? 1.0 - z : z)));
          // Small differences in the blowdown, whose flow would otherwise follow the cooling of the gas.
          const double thermal = target > 1.0e6 ? 1.0 : 0.02;
          column.gasTemperature[grid] = feedTemperature + thermal * (20.0 * front + 4.0 * z);
          column.solidTemperature[grid] = column.gasTemperature[grid] + thermal * 6.0 * front;
          column.wallTemperature[grid] = feedTemperature + thermal;
          const double total = 1.0e6 / (R * column.gasTemperature[grid]);
          double adsorbing = 0.0;
          for (size_t comp = 1; comp < n; ++comp)
          {
            column.concentration[grid * n + comp] = feed[comp] * front * 0.8 + 1.0e-3;
            adsorbing += column.concentration[grid * n + comp];
          }
          column.concentration[grid * n] = total - adsorbing;
        }
        if (reverse)
        {
          std::copy_n(column.concentration.begin() + static_cast<std::ptrdiff_t>(n), n, column.concentration.begin());
        }

        // Near equilibrium, so that the ramp sets the sign of the flow.
        Timing timings;
        precompute(column, timings);
        for (size_t index = 0; index < column.physisorption.size(); ++index)
        {
          column.physisorption[index] = 0.99 * column.equilibriumPhysisorption[index];
        }
        precompute(column, timings);
        computeDerivatives(column);
        precompute(column, timings);
        computeDerivatives(column);

        EXPECT_EQ(column.interstitialGasVelocity[reverse ? 0 : gridPoints], 0.0);
        double scale = 0.0;
        for (double rate : column.concentrationDot) scale = std::max(scale, std::abs(rate));
        ASSERT_GT(scale, 0.0);
        std::vector<double> accumulation(n, 0.0);
        for (size_t grid = 1; grid <= gridPoints; ++grid)
        {
          double sum = 0.0;
          double sumDot = 0.0;
          for (size_t comp = 0; comp < n; ++comp)
          {
            const size_t index = grid * n + comp;
            sum += column.concentration[index];
            sumDot += column.concentrationDot[index];
            accumulation[comp] += column.totalVoidFraction[grid] * gridSpacing(column.columnDistances, grid) *
                                  (column.concentrationDot[index] +
                                   column.loadingPrefactor[grid] * column.physisorptionDot[index]);
          }
          const double expected =
              sum * (column.relativePressureRate - column.gasTemperatureDot[grid] / column.gasTemperature[grid]);
          EXPECT_NEAR(sumDot, expected, 1.0e-10 * scale) << grid;
        }

        // Pressurization draws the feed in through the open end; a blowdown lets out what is in the cell next to it.
        const double openFlux = reverse
                                    ? -column.totalVoidFraction[gridPoints] * column.interstitialGasVelocity[gridPoints]
                                    : column.totalVoidFraction[0] * column.interstitialGasVelocity[0];
        EXPECT_EQ(openFlux > 0.0, target > 1.0e6);
        const double fluxScale = std::abs(openFlux) * column.totalConcentration[1];
        const size_t cell = reverse ? gridPoints : 1;
        for (size_t comp = 0; comp < n; ++comp)
        {
          const double gas = openFlux > 0.0 ? feed[comp] * column.inletPressure / 1.0e6
                                            : column.concentration[cell * n + comp];
          EXPECT_NEAR(accumulation[comp], openFlux * gas, 1.0e-9 * fluxScale) << comp;
        }
      }
    }
  }
}

TEST(PressureRamp, AddsThePressureWorkToTheGasTemperature)
{
  // At uniform temperatures only the pressure work changes the gas temperature: dT/dt = (dP/dt) / (rho c_p).
  nlohmann::json data = t15Input(8);
  data["SwingAdsorptionPhases"] = {{{"Duration", 1.0}}};
  SwingAdsorption swing(readInput(data));
  Column& column = swing.breakthrough.column;
  column.closedEnd = true;
  column.pressureRamp = {Column::PressureRamp::Shape::Linear, 100.0, 0.0, 1.0e6, 2.0e6, 0.0};
  column.advancePressureRamp(0.0);
  Timing timings;
  precompute(column, timings);
  computeDerivatives(column);
  ASSERT_EQ(column.pressureRate, 1.0e4);
  for (size_t grid = 1; grid <= column.numberOfGridPoints; ++grid)
  {
    const double expected = 1.0e4 / (column.gasDensity[grid] * column.heatCapacityGas);
    EXPECT_NEAR(column.gasTemperatureDot[grid], expected, 1.0e-10 * expected) << grid;
  }
}

TEST(SwingStep, RescalesTheGasAtAStepOrWhenThePressureMovesToTheOtherEnd)
{
  nlohmann::json data = t15Input(12);
  data["energyBalance"] = false;
  data["WriteInterval"] = 1000.0;
  const auto run = [&](const nlohmann::json& second, const std::string& boundary = "InletPressureInletVelocity")
  {
    data["BoundaryCondition"] = boundary;
    data["SwingAdsorptionPhases"] = {{{"Name", "Feed"}, {"Duration", 50.0}}, second};
    auto swing = std::make_unique<SwingAdsorption>(readInput(data));
    swing->run();
    return swing;
  };

  {
    SCOPED_TRACE("pressure step");
    const auto swing = run({{"Name", "Low"}, {"Duration", 50.0}, {"InletPressure", 1.0e5}});
    const std::vector<SwingAdsorptionPhaseRecord>& records = swing->phaseRecords;
    ASSERT_EQ(records.size(), 2U);
    EXPECT_EQ(records[1].adsorbedAtStart, records[0].adsorbedAtEnd);
    for (size_t comp = 0; comp < records[0].gasAtEnd.size(); ++comp)
    {
      EXPECT_NEAR(records[1].gasAtStart[comp], 0.1 * records[0].gasAtEnd[comp], 2.0e-3 * records[0].gasAtEnd[comp]);
    }
    // What remains is the change of the pressure drop, about 100 Pa here.
    EXPECT_LT(swing->breakthrough.column.idealGasDeviation(), 3.0e-3);
    EXPECT_NEAR(records[1].pressureAtZero, 1.0e5, 1.0e-9 * 1.0e5);
    EXPECT_LT(records[1].pressureAtLength, 1.0e5);
  }
  {
    SCOPED_TRACE("temperature step");
    const auto swing = run({{"Name", "Hot"}, {"Duration", 50.0}, {"Temperature", 373.0}});
    const std::vector<SwingAdsorptionPhaseRecord>& records = swing->phaseRecords;
    for (size_t comp = 0; comp < records[0].gasAtEnd.size(); ++comp)
    {
      EXPECT_NEAR(records[1].gasAtStart[comp], 298.0 / 373.0 * records[0].gasAtEnd[comp],
                  2.0e-3 * records[0].gasAtEnd[comp]);
    }
    EXPECT_LT(swing->breakthrough.column.idealGasDeviation(), 3.0e-3);
  }
  {
    SCOPED_TRACE("no step");
    const auto swing = run({{"Name", "Same"}, {"Duration", 50.0}, {"ColumnEntranceVelocity", 0.05}});
    EXPECT_EQ(swing->phaseRecords[1].gasAtStart, swing->phaseRecords[0].gasAtEnd);
  }
  {
    SCOPED_TRACE("pressure set at the other end");
    // The inlet pressure moves to z = L at the same value, so that the profile drops by the pressure drop, about 40 Pa.
    const auto swing =
        run({{"Name", "Reverse"}, {"Duration", 50.0}, {"FlowDirection", "Reverse"}, {"GasPhaseMolFractions", {{"H2", 1.0}}}});
    const std::vector<SwingAdsorptionPhaseRecord>& records = swing->phaseRecords;
    double before = 0.0;
    double after = 0.0;
    for (size_t comp = 0; comp < records[0].gasAtEnd.size(); ++comp)
    {
      before += records[0].gasAtEnd[comp];
      after += records[1].gasAtStart[comp];
    }
    EXPECT_LT(after, before);
    EXPECT_GT(after, (1.0 - 1.0e-4) * before);
  }
  {
    SCOPED_TRACE("FixedVelocity keeps the gas");
    const auto swing = run({{"Name", "Low"}, {"Duration", 50.0}, {"InletPressure", 1.0e5}}, "FixedVelocity");
    EXPECT_EQ(swing->phaseRecords[1].gasAtStart, swing->phaseRecords[0].gasAtEnd);
  }
}

TEST(SwingStep, AStepThatTheFlowCannotCarryNamesTheRemedy)
{
  // BEA at 10 mbar: the C7 that desorbs at once after a step from 1 bar cannot leave through the bed at that pressure.
  // A ramp, or an open end on the low-pressure side, lets it.
  std::ifstream input(std::filesystem::path{RUPTURA_SOURCE_DIR} / "examples/BEA-alkanes-C7/psa/simulation.json");
  nlohmann::json data = nlohmann::json::parse(input);
  data["NumberOfGridPoints"] = 12;
  data["WriteInterval"] = 1000.0;
  data["SwingAdsorptionPhases"][0]["Duration"] = 500.0;
  data["SwingAdsorptionPhases"][1]["Duration"] = 200.0;
  data["SwingAdsorptionPhases"][1].erase("PressureRamp");
  try
  {
    SwingAdsorption swing(readInput(data));
    swing.run();
    ADD_FAILURE() << "the step should not be possible";
  }
  catch (const std::runtime_error& error)
  {
    EXPECT_NE(std::string(error.what()).find("swing phase 'Blowdown'"), std::string::npos) << error.what();
    EXPECT_NE(std::string(error.what()).find("PressureRamp"), std::string::npos) << error.what();
  }

  data["SwingAdsorptionPhases"][1]["PressureRamp"] = {{"Shape", "Exponential"}, {"TimeConstant", 50.0}};
  EXPECT_NO_THROW({
    SwingAdsorption swing(readInput(data));
    swing.run();
  });
  data["SwingAdsorptionPhases"][1]["ClosedEnd"] = "ProductEnd";
  SwingAdsorption swing(readInput(data));
  swing.run();
  EXPECT_LT(swing.breakthrough.column.interstitialGasVelocity[0], 0.0);
  EXPECT_LT(swing.breakthrough.column.idealGasDeviation(), 1.0e-2);
}

TEST(ClosedEnd, IntegratorsAgreeOnAPressurizationAndABlowdown)
{
  for (const std::string shape : {"Exponential", "Linear"})
  {
    nlohmann::json data = t15Input(12);
    data["WriteInterval"] = 1000.0;
    data["InletPressure"] = 1.0e5;
    data["SwingAdsorptionPhases"] = {closedPhase("ProductEnd", 1.0e6, 25.0, shape, 10.0),
                                     closedPhase("ProductEnd", 1.0e5, 25.0, shape, 10.0)};
    std::vector<std::vector<double>> results;
    for (const std::string integrator : {"RungeKutta3", "CVODE", "ROS34PW2"})
    {
      SCOPED_TRACE(shape + " " + integrator);
      data["BreakthroughIntegrator"] = integrator;
      data.erase("TimeStep");
      data.erase("CVODERelativeTolerance");
      data.erase("ROS34PW2RelativeTolerance");
      if (integrator == "RungeKutta3") data["TimeStep"] = 4.0e-3;
      if (integrator == "CVODE") data["CVODERelativeTolerance"] = 1.0e-7;
      if (integrator == "ROS34PW2") data["ROS34PW2RelativeTolerance"] = 1.0e-6;
      SwingAdsorption swing(readInput(data));
      swing.run();
      std::vector<double> result;
      for (const SwingAdsorptionPhaseRecord& record : swing.phaseRecords)
      {
        result.insert(result.end(), record.gasAtEnd.begin(), record.gasAtEnd.end());
        result.insert(result.end(), record.adsorbedAtEnd.begin(), record.adsorbedAtEnd.end());
      }
      const Column& column = swing.breakthrough.column;
      for (double temperature : column.solidTemperature) result.push_back(temperature);
      results.push_back(std::move(result));
    }
    double scale = 0.0;
    for (double value : results[0]) scale = std::max(scale, std::abs(value));
    for (size_t index = 0; index < results[0].size(); ++index)
    {
      EXPECT_NEAR(results[1][index], results[0][index], 3.0e-4 * scale) << shape << " CVODE " << index;
      EXPECT_NEAR(results[2][index], results[0][index], 3.0e-4 * scale) << shape << " ROS34PW2 " << index;
    }
  }
}

TEST(TemperatureSwing, Ros34pw2KeepsTheGasAtTheClosurePressureAsCvodeDoes)
{
  // Two cycles of adsorption at 298 K and regeneration at 373 K with the energy balance (the T15 TSA of 3.4.2 on 12
  // cells). W leaves out the closure; without its projection (B-ITG-16) max |R T sum(c) / P - 1| reached 2.1e-3 and
  // the gas amounts were 1.6e-4 off those of CVODE.
  nlohmann::json data = t15Input(12);
  data["WriteInterval"] = 1000.0;
  data["NumberOfCycles"] = 2;
  data["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"Temperature", 298.0}, {"Duration", 1000.0}},
                                   {{"Name", "Regeneration"}, {"Temperature", 373.0}, {"Duration", 1000.0}}};
  std::vector<std::vector<double>> gas;
  for (const std::string integrator : {"CVODE", "ROS34PW2"})
  {
    SCOPED_TRACE(integrator);
    data["BreakthroughIntegrator"] = integrator;
    data.erase("CVODERelativeTolerance");
    if (integrator == "CVODE") data["CVODERelativeTolerance"] = 1.0e-7;
    SwingAdsorption swing(readInput(data));
    swing.run();
    ASSERT_EQ(swing.phaseRecords.size(), 4U);
    std::vector<double> amounts;
    for (const SwingAdsorptionPhaseRecord& record : swing.phaseRecords)
    {
      EXPECT_LT(record.gasDeviation, 2.0e-4) << record.name;
      amounts.insert(amounts.end(), record.gasAtEnd.begin(), record.gasAtEnd.end());
    }
    gas.push_back(std::move(amounts));
  }
  const double scale = *std::max_element(gas[0].begin(), gas[0].end());
  for (size_t index = 0; index < gas[0].size(); ++index)
  {
    EXPECT_NEAR(gas[1][index], gas[0][index], 2.0e-5 * scale) << index;
  }
}

TEST(SkarstromCycle, MassBalanceClosesOnThePhaseRecordsAndTheOutletFluxes)
{
  // Feed pressurization, adsorption, counter-current blowdown, and a counter-current H2 purge, twice: what the outlet
  // files show leaving an open step is its known inflow minus the change of the recorded amounts.
  const double high = 1.0e6;
  const double low = 1.0e5;
  const double feedVelocity = 0.03;
  const double purgeVelocity = 0.05;
  const double interval = 0.25;
  nlohmann::json data = t15Input(30);
  data["ConvectionScheme"] = "Upwind";  // the face value at the outlet is that of the last cell
  data["InletPressure"] = low;
  data["WriteInterval"] = interval;
  data["NumberOfCycles"] = 2;
  data["SwingAdsorptionPhases"] = {
      closedPhase("ProductEnd", high, 30.0, "Exponential", 5.0),
      {{"Name", "Adsorption"}, {"InletPressure", high}, {"ColumnEntranceVelocity", feedVelocity}, {"Duration", 60.0}},
      closedPhase("ProductEnd", low, 30.0, "Exponential", 5.0),
      {{"Name", "Purge"},
       {"InletPressure", low},
       {"ColumnEntranceVelocity", purgeVelocity},
       {"FlowDirection", "Reverse"},
       {"GasPhaseMolFractions", {{"H2", 1.0}}},
       {"Duration", 60.0}}};
  const InputReader reader = readInput(data);
  SwingAdsorption swing(reader);
  swing.run();

  const Column& column = swing.breakthrough.column;
  const size_t n = column.numberOfComponents;
  const size_t rows = column.numberOfGridPoints + 1;
  ASSERT_EQ(swing.phaseRecords.size(), 8U);

  const std::vector<double> velocity = readColumn("column.data", 3);
  std::vector<std::vector<double>> concentration;
  for (size_t comp = 0; comp < n; ++comp)
  {
    concentration.push_back(
        readColumn("component_" + std::to_string(comp) + "_" + column.components[comp].name + ".data", 3));
    ASSERT_EQ(concentration.back().size(), velocity.size());
  }

  // The amount of a component leaving through a node in [start, end) of the output blocks: the output at the end of the
  // phase holds the velocity of the next one, so that the last interval takes its previous value.
  const auto leaving = [&](size_t comp, double start, double end, size_t velocityNode, size_t concentrationNode)
  {
    const size_t first = static_cast<size_t>(std::round(start / interval));
    const size_t last = static_cast<size_t>(std::round(end / interval)) - 1;
    const double fraction = column.totalVoidFraction[velocityNode];
    const auto flux = [&](size_t block)
    {
      return fraction * std::abs(velocity[block * rows + velocityNode]) *
             concentration[comp][block * rows + concentrationNode];
    };
    double amount = interval * flux(last);
    for (size_t block = first; block < last; ++block) amount += 0.5 * interval * (flux(block) + flux(block + 1));
    return amount;
  };

  const std::vector<double> feedFractions = [&]
  {
    std::vector<double> fractions;
    for (const Component& component : reader.components) fractions.push_back(component.initialGasMoleFraction);
    return fractions;
  }();
  const double inletFraction = column.totalVoidFraction[0];
  const double temperature = 298.0;
  for (size_t cycle = 0; cycle < 2; ++cycle)
  {
    SCOPED_TRACE(cycle);
    const SwingAdsorptionPhaseRecord& pressurization = swing.phaseRecords[4 * cycle];
    const SwingAdsorptionPhaseRecord& adsorption = swing.phaseRecords[4 * cycle + 1];
    const SwingAdsorptionPhaseRecord& blowdown = swing.phaseRecords[4 * cycle + 2];
    const SwingAdsorptionPhaseRecord& purge = swing.phaseRecords[4 * cycle + 3];

    // Only the feed enters a pressurization.
    const std::vector<double> entered = amountChange(pressurization);
    double total = 0.0;
    for (double value : entered) total += value;
    EXPECT_GT(total, 0.0);
    // Only gas leaves a blowdown.
    const std::vector<double> released = amountChange(blowdown);
    for (size_t comp = 1; comp < n; ++comp)
    {
      EXPECT_NEAR(entered[comp], feedFractions[comp] * total, 1.0e-4 * total) << comp;
      EXPECT_LE(released[comp], 1.0e-9) << comp;
    }

    const std::vector<double> adsorbed = amountChange(adsorption);
    const std::vector<double> purged = amountChange(purge);
    const double feedAmount =
        inletFraction * feedVelocity * high / (R * temperature) * (adsorption.endTime - adsorption.startTime);
    const double purgeAmount =
        inletFraction * purgeVelocity * low / (R * temperature) * (purge.endTime - purge.startTime);
    double productTotal = 0.0;
    for (size_t comp = 0; comp < n; ++comp)
    {
      const double product = leaving(comp, adsorption.startTime, adsorption.endTime, rows - 1, rows - 1);
      productTotal += product;
      EXPECT_NEAR(feedAmount * feedFractions[comp] - adsorbed[comp], product, 2.0e-3 * feedAmount) << comp;
      const double purgeIn = column.components[comp].name == "H2" ? purgeAmount : 0.0;
      EXPECT_NEAR(purgeIn - purged[comp], leaving(comp, purge.startTime, purge.endTime, 0, 1), 3.0e-3 * purgeAmount)
          << comp;
    }
    EXPECT_GT(productTotal, 0.1 * feedAmount);
  }
}

TEST(PressureSwingInput, ReadsClosedEndsAndRampsInAnyCaseAndRejectsWhatTheyCannotUse)
{
  nlohmann::json data = t15Input(8);
  data["NumberOfCycles"] = 2;
  const nlohmann::json pressurization = {{"Name", "Pressurization"},
                                         {"Duration", 10.0},
                                         {"InletPressure", 1.0e6},
                                         {"closedend", "productend"},
                                         {"PRESSURERAMP", {{"shape", "exponential"}, {"TIMECONSTANT", 2.0}}}};
  const nlohmann::json blowdown = closedPhase("FeedEnd", 1.0e5, 10.0, "Linear", 5.0);
  data["SwingAdsorptionPhases"] = {pressurization, blowdown};
  {
    const InputReader reader = readInput(data);
    ASSERT_EQ(reader.swingAdsorptionPhases.size(), 4U);
    const SwingAdsorption swing(reader);
    EXPECT_TRUE(swing.subStages[2].closedEnd);
    EXPECT_EQ(swing.subStages[2].flowDirection, Column::FlowDirection::Forward);
    EXPECT_EQ(swing.subStages[2].rampShape, Column::PressureRamp::Shape::Exponential);
    EXPECT_EQ(swing.subStages[2].rampTimeConstant, 2.0);
    EXPECT_TRUE(swing.subStages[3].closedEnd);
    EXPECT_EQ(swing.subStages[3].flowDirection, Column::FlowDirection::Reverse);
    EXPECT_EQ(swing.subStages[3].rampShape, Column::PressureRamp::Shape::Linear);
  }

  const auto expectRejected = [&](nlohmann::json input, const std::string& expected)
  {
    try
    {
      static_cast<void>(readInput(input));
      ADD_FAILURE() << "accepted; expected: " << expected;
    }
    catch (const std::runtime_error& error)
    {
      EXPECT_NE(std::string(error.what()).find(expected), std::string::npos) << error.what();
    }
  };
  const auto withPhase = [&](const std::string& key, const nlohmann::json& value)
  {
    nlohmann::json input = data;
    input["SwingAdsorptionPhases"][1][key] = value;
    return input;
  };
  expectRejected(withPhase("ClosedEnd", "Middle"), "invalid ClosedEnd");
  expectRejected(withPhase("FlowDirection", "Forward"), "remove FlowDirection");
  expectRejected(withPhase("ColumnEntranceVelocity", 0.1), "closed-end phase");
  expectRejected(withPhase("ProductEndGasPhaseMolFractions", {{"H2", 1.0}}), "closed-end phase");
  expectRejected(withPhase("PressureRamp", 5.0), "PressureRamp must be an object");
  expectRejected(withPhase("PressureRamp", {{"Shape", "Linear"}}), "TimeConstant");
  expectRejected(withPhase("PressureRamp", {{"TimeConstant", 5.0}}), "Shape");
  expectRejected(withPhase("PressureRamp", {{"Shape", "Cubic"}, {"TimeConstant", 5.0}}), "invalid Shape");
  expectRejected(withPhase("PressureRamp", {{"Shape", "Linear"}, {"TimeConstant", 0.0}}), "positive time");
  expectRejected(withPhase("PressureRamp", {{"Shape", "Linear"}, {"TimeConstant", 5.0}, {"Tau", 1.0}}), "Tau");
  nlohmann::json input = data;
  input["InletCondition"] = "Dirichlet";
  expectRejected(input, "needs InletCondition Danckwerts");
  input = data;
  input["BoundaryCondition"] = "FixedVelocity";
  expectRejected(input, "InletPressureInletVelocity");
}

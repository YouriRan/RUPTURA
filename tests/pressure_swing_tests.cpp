#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
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
  // The test reads back the output files, which a test running at the same time would overwrite in a shared directory.
  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_skarstrom_mass_balance");
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

namespace
{
// A phase that feeds the column until the criterion holds, for at most the duration.
nlohmann::json stoppedPhase(const nlohmann::json& criterion, double duration)
{
  return {{"Name", "Adsorption"}, {"Duration", duration}, {"StopCriterion", criterion}};
}

const nlohmann::json hydrogenAtTheOutlet = {
    {"Metric", "MoleFraction"}, {"Component", "H2"}, {"End", "ProductEnd"}, {"Operator", ">="}, {"Value", 0.15}};
}  // namespace

TEST(StopCriterion, EndsThePhaseWhereTheOutletMoleFractionCrossesWithEveryIntegrator)
{
  // T15 on 8 cells: the phase ends once the hydrogen mole fraction at z = L reaches 0.15, long before its 1000 s. The
  // state at the end holds the value, and the three integrators end at the same time.
  nlohmann::json data = t15Input(8);
  data["WriteInterval"] = 50.0;
  data["SwingAdsorptionPhases"] = nlohmann::json::array({stoppedPhase(hydrogenAtTheOutlet, 1000.0)});
  std::vector<double> endTimes;
  for (const std::string integrator : {"CVODE", "ROS34PW2", "RungeKutta3"})
  {
    SCOPED_TRACE(integrator);
    data["BreakthroughIntegrator"] = integrator;
    data.erase("TimeStep");
    if (integrator == "RungeKutta3") data["TimeStep"] = 4.0e-3;
    SwingAdsorption swing(readInput(data));
    swing.run();
    const Column& column = swing.breakthrough.column;
    ASSERT_EQ(swing.phaseRecords.size(), 1U);
    endTimes.push_back(swing.phaseRecords[0].endTime);
    EXPECT_LT(endTimes.back(), 500.0);
    EXPECT_NEAR(column.moleFraction[column.numberOfGridPoints * column.numberOfComponents + 3], 0.15, 1.0e-6);
  }
  EXPECT_NEAR(endTimes[1], endTimes[0], 1.0e-3 * endTimes[0]);
  EXPECT_NEAR(endTimes[2], endTimes[0], 1.0e-3 * endTimes[0]);

  // The same run without the criterion, written every second: its outlet file crosses 0.15 at that time, and the phase
  // after a stopped one starts where it ended.
  data["BreakthroughIntegrator"] = "CVODE";
  data.erase("TimeStep");
  data["WriteInterval"] = 1.0;
  data["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"Duration", 200.0}}};
  const InputReader reference = readInput(data);
  data["SwingAdsorptionPhases"] = {stoppedPhase(hydrogenAtTheOutlet, 1000.0), {{"Name", "Hold"}, {"Duration", 20.0}}};
  const InputReader twoPhases = readInput(data);
  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_stop_criterion");
  SwingAdsorption stopped(twoPhases);
  stopped.run();
  ASSERT_EQ(stopped.phaseRecords.size(), 2U);
  EXPECT_NEAR(stopped.phaseRecords[0].endTime, endTimes[0], 1.0e-5 * endTimes[0]);  // other output times
  EXPECT_EQ(stopped.phaseRecords[1].startTime, stopped.phaseRecords[0].endTime);
  EXPECT_DOUBLE_EQ(stopped.phaseRecords[1].endTime, stopped.phaseRecords[0].endTime + 20.0);

  SwingAdsorption unstopped(reference);
  unstopped.run();
  const std::vector<double> fraction = readColumn("component_3_H2.data", 5);
  const size_t rows = unstopped.breakthrough.column.numberOfGridPoints + 1;
  double crossing = 0.0;
  for (size_t output = 1; output < fraction.size() / rows && crossing == 0.0; ++output)
  {
    const double before = fraction[output * rows - 1];
    const double after = fraction[(output + 1) * rows - 1];
    if (after >= 0.15) crossing = static_cast<double>(output - 1) + (0.15 - before) / (after - before);
  }
  EXPECT_NEAR(crossing, endTimes[0], 0.2);
}

TEST(StopCriterion, FollowsTheMinimumDurationTheAmountsAndTheTimeoutPolicy)
{
  nlohmann::json data = t15Input(8);
  data["WriteInterval"] = 1000.0;
  const auto endTime = [&](const nlohmann::json& phase)
  {
    data["SwingAdsorptionPhases"] = nlohmann::json::array({phase});
    SwingAdsorption swing(readInput(data));
    swing.run();
    return swing.phaseRecords[0].endTime;
  };

  // A criterion that holds from the start ends the phase at its minimum duration.
  const nlohmann::json always = {{"Metric", "Pressure"}, {"End", "FeedEnd"}, {"Operator", ">="}, {"Value", 1.0}};
  for (const std::string integrator : {"CVODE", "ROS34PW2", "RungeKutta3"})
  {
    SCOPED_TRACE(integrator);
    data["BreakthroughIntegrator"] = integrator;
    data.erase("TimeStep");
    if (integrator == "RungeKutta3") data["TimeStep"] = 4.0e-3;
    nlohmann::json phase = stoppedPhase(always, 30.0);
    phase["MinimumDuration"] = 12.5;
    EXPECT_NEAR(endTime(phase), 12.5, 1.0e-5);
  }
  data["BreakthroughIntegrator"] = "CVODE";
  data.erase("TimeStep");

  // The feed enters at eps v P / (R T) per column cross-section; the phase ends once 100 mol/m^2 have entered.
  Column column(readInput(data));
  column.initialize();
  const double feedFlux = column.totalVoidFraction[0] * column.columnEntranceVelocity * 1.0e6 / (R * 298.0);
  const nlohmann::json fed = {{"Metric", "AmountIn"}, {"End", "FeedEnd"}, {"Operator", ">="}, {"Value", 100.0}};
  const double fedTime = endTime(stoppedPhase(fed, 500.0));
  EXPECT_NEAR(fedTime, 100.0 / feedFlux, 1.0e-4 * fedTime);

  // All criteria, or any of them.
  const double hydrogenTime = endTime(stoppedPhase(hydrogenAtTheOutlet, 500.0));
  EXPECT_GT(hydrogenTime, 2.0 * fedTime);
  EXPECT_NEAR(endTime(stoppedPhase({{"All", {fed, hydrogenAtTheOutlet}}}, 500.0)), hydrogenTime, 1.0e-6 * hydrogenTime);
  EXPECT_NEAR(endTime(stoppedPhase({{"Any", {hydrogenAtTheOutlet, fed}}}, 500.0)), fedTime, 1.0e-6 * fedTime);

  // A criterion that never holds: the phase lasts its duration, or the run stops with TimeoutPolicy Error.
  const nlohmann::json never = {{"Metric", "MaximumTemperature"}, {"Operator", ">="}, {"Value", 1.0e4}};
  EXPECT_DOUBLE_EQ(endTime(stoppedPhase(never, 5.0)), 5.0);
  nlohmann::json phase = stoppedPhase(never, 5.0);
  phase["TimeoutPolicy"] = "Error";
  data["SwingAdsorptionPhases"] = nlohmann::json::array({phase});
  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_stop_timeout");
  SwingAdsorption swing(readInput(data));
  EXPECT_THROW(swing.run(), std::runtime_error);
}

TEST(CyclicSteadyState, EndsTheRunAtTheFirstCyclesWithinTheToleranceWithStreamsThatCloseTheMassBalance)
{
  // The isothermal Skarstrom cycle of T15 on 10 cells, for at most 60 cycles.
  nlohmann::json data = t15Input(10);
  data["energyBalance"] = false;
  data["InletPressure"] = 1.0e5;
  data["WriteInterval"] = 1000.0;
  data["NumberOfCycles"] = 60;
  data["CyclicSteadyState"] = {{"Tolerance", 2.0e-3}, {"ConsecutiveCycles", 2}};
  data["SwingAdsorptionPhases"] = {
      closedPhase("ProductEnd", 1.0e6, 30.0, "Linear", 20.0),
      {{"Name", "Adsorption"}, {"InletPressure", 1.0e6}, {"ColumnEntranceVelocity", 0.03}, {"Duration", 60.0}},
      closedPhase("ProductEnd", 1.0e5, 30.0, "Linear", 20.0),
      {{"Name", "Purge"},
       {"InletPressure", 1.0e5},
       {"ColumnEntranceVelocity", 0.05},
       {"FlowDirection", "Reverse"},
       {"GasPhaseMolFractions", "Adsorption"},
       {"Duration", 60.0}}};
  const InputReader reader = readInput(data);
  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_cyclic_steady_state");
  SwingAdsorption swing(reader);
  swing.run();

  const size_t cycles = readColumn("cycles.data", 0).size();
  ASSERT_GT(cycles, 3U);
  EXPECT_LT(cycles, 60U);
  EXPECT_EQ(swing.phaseRecords.size(), 4 * cycles);
  // The largest of the three changes of each cycle: within the tolerance in the last two cycles, not in the one before.
  std::vector<double> change(cycles, 0.0);
  for (size_t column = 2; column < 5; ++column)
  {
    const std::vector<double> values = readColumn("cycles.data", column);
    for (size_t cycle = 0; cycle < cycles; ++cycle) change[cycle] = std::max(change[cycle], values[cycle]);
  }
  EXPECT_LE(change[cycles - 1], 2.0e-3);
  EXPECT_LE(change[cycles - 2], 2.0e-3);
  EXPECT_GT(change[cycles - 3], 2.0e-3);
  for (const double unaccounted : readColumn("cycles.data", 5)) EXPECT_LT(unaccounted, 1.0e-3);

  // The streams of the last cycle: the feed of the adsorption is known, nothing passes a closed end, and the purge
  // takes the gas that left the adsorption of its cycle.
  const size_t n = swing.breakthrough.column.numberOfComponents;
  const std::span<const SwingAdsorptionPhaseRecord> last(swing.phaseRecords.end() - 4, 4);
  const std::vector<double>& adsorption = last[1].streamAmounts;
  const std::vector<double>& purge = last[3].streamAmounts;
  ASSERT_EQ(adsorption.size(), 4 * n);
  const double fed = swing.breakthrough.column.totalVoidFraction[0] * 0.03 * 1.0e6 / (R * 298.0) * 60.0;
  double product = 0.0;
  double purgeGas = 0.0;
  for (size_t comp = 0; comp < n; ++comp)
  {
    EXPECT_NEAR(adsorption[comp], fed * reader.components[comp].initialGasMoleFraction, 1.0e-4 * fed) << comp;
    EXPECT_EQ(last[0].streamAmounts[2 * n + comp] + last[0].streamAmounts[3 * n + comp], 0.0) << comp;
    product += adsorption[3 * n + comp];
    purgeGas += purge[2 * n + comp];
  }
  ASSERT_GT(product, 0.1 * fed);
  for (size_t comp = 0; comp < n; ++comp)
  {
    EXPECT_NEAR(purge[2 * n + comp] / purgeGas, adsorption[3 * n + comp] / product, 1.0e-6) << comp;
  }
}

namespace
{
// T15 with the pressure of the gas and the velocity of the pressure-drop law, ROS34PW2, and the molar heat capacities
// of helium, CO2, N2, and H2; the column starts at the outlet pressure.
nlohmann::json momentumBalanceInput(size_t gridPoints, double outletPressure)
{
  nlohmann::json data = t15Input(gridPoints);
  data["BreakthroughIntegrator"] = "ROS34PW2";
  data["BoundaryCondition"] = "MomentumBalance";
  data["OutletPressure"] = outletPressure;
  data["WriteInterval"] = 1000.0;
  const std::vector<double> heatCapacities{20.8, 37.1, 29.1, 28.8};
  for (size_t comp = 0; comp < heatCapacities.size(); ++comp)
  {
    data["Components"][comp]["GasHeatCapacity"] = heatCapacities[comp];
  }
  return data;
}
}  // namespace

TEST(MomentumBalance, MatchesTheContinuityClosureInABreakthrough)
{
  // 400 s of the T15 feed at its velocity into a column at 1 MPa, with gas leaving to 1 MPa at z = L: the closure
  // InletVelocityOutletPressure and the momentum balance give the same fronts, temperatures, and pressure drop.
  for (const bool energyBalance : {false, true})
  {
    SCOPED_TRACE(energyBalance);
    nlohmann::json data = momentumBalanceInput(12, 1.0e6);
    data["energyBalance"] = energyBalance;
    data["SwingAdsorptionPhases"] = {
        {{"Name", "Adsorption"},
         {"Duration", 400.0},
         {"FeedEnd", {{"Mode", "Velocity"}, {"Velocity", data["ColumnEntranceVelocity"]}, {"Pressure", 1.0e6}}},
         {"ProductEnd", {{"Mode", "Pressure"}, {"Pressure", 1.0e6}}}}};
    SwingAdsorption momentum(readInput(data));
    momentum.run();

    data["BoundaryCondition"] = "InletVelocityOutletPressure";
    data["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"Duration", 400.0}}};
    SwingAdsorption closure(readInput(data));
    closure.run();

    const Column& column = momentum.breakthrough.column;
    const Column& reference = closure.breakthrough.column;
    const size_t n = column.numberOfComponents;
    const size_t last = column.numberOfGridPoints;
    const double feedConcentration = 1.0e6 / (R * 298.0);
    for (size_t index = n; index < column.concentration.size(); ++index)
    {
      EXPECT_NEAR(column.concentration[index], reference.concentration[index], 2.0e-3 * feedConcentration) << index;
      EXPECT_NEAR(column.physisorption[index], reference.physisorption[index], 5.0e-3) << index;
    }
    for (size_t node = 1; node <= last; ++node)
    {
      EXPECT_NEAR(column.solidTemperature[node], reference.solidTemperature[node], 0.05) << node;
      EXPECT_NEAR(column.interstitialGasVelocity[node], reference.interstitialGasVelocity[node], 2.0e-4) << node;
    }
    const double pressureDrop = reference.totalPressure[1] - reference.totalPressure[last];
    EXPECT_GT(pressureDrop, 10.0);
    EXPECT_NEAR(column.totalPressure[1] - column.totalPressure[last], pressureDrop, 0.05 * pressureDrop);
    EXPECT_NEAR(column.totalPressure[last], 1.0e6, 5.0);
  }
}

TEST(MomentumBalance, APressurePortFillsAndEmptiesAnInertColumnAndClosedPortsKeepItsGas)
{
  // Helium only: a port at z = 0 takes the column from 0.1 to 0.5 MPa, the column is closed, and a port at z = L takes
  // it back. What passes a port is the change of the gas, P V / (R T) of the fluid volume, and a closed column keeps
  // its gas and levels its pressure.
  nlohmann::json data = momentumBalanceInput(10, 1.0e5);
  data["energyBalance"] = false;
  data["Components"][0]["GasPhaseMolFraction"] = 1.0;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["GasPhaseMolFraction"] = 0.0;
  }
  data["SwingAdsorptionPhases"] = {
      {{"Name", "Fill"},
       {"Duration", 30.0},
       {"FeedEnd",
        {{"Mode", "Pressure"}, {"Pressure", 5.0e5}, {"PressureRamp", {{"Shape", "Linear"}, {"TimeConstant", 20.0}}}}}},
      {{"Name", "Hold"}, {"Duration", 10.0}},
      {{"Name", "Empty"},
       {"Duration", 40.0},
       {"ProductEnd",
        {{"Mode", "Pressure"},
         {"Pressure", 1.0e5},
         {"PressureRamp", {{"Shape", "Exponential"}, {"TimeConstant", 5.0}}}}}}};
  for (const std::string integrator : {"ROS34PW2", "CVODE"})
  {
    SCOPED_TRACE(integrator);
    data["BreakthroughIntegrator"] = integrator;
    data["CVODELinearSolver"] = "KLU";
    SwingAdsorption swing(readInput(data));
    swing.run();
    ASSERT_EQ(swing.phaseRecords.size(), 3U);
    const Column& column = swing.breakthrough.column;
    const size_t n = column.numberOfComponents;
    double fluidVolume = 0.0;
    for (size_t node = 1; node <= column.numberOfGridPoints; ++node)
    {
      fluidVolume += column.totalVoidFraction[node] * gridSpacing(column.columnDistances, node);
    }
    const double amount = fluidVolume * 4.0e5 / (R * 298.0);
    // ROS34PW2 keeps the stream amounts with the column amounts. With CVODE they hold its tolerance on a flux that
    // follows a small pressure difference, which is less.
    const double streamTolerance = (integrator == "CVODE" ? 1.0e-2 : 1.0e-5) * amount;

    const SwingAdsorptionPhaseRecord& fill = swing.phaseRecords[0];
    EXPECT_NEAR(amountChange(fill)[0], amount, 1.0e-3 * amount);
    EXPECT_NEAR(fill.streamAmounts[0], amountChange(fill)[0], streamTolerance);
    EXPECT_NEAR(fill.pressureAtLength, 5.0e5, 500.0);
    for (size_t index = 1; index < 4 * n; ++index) EXPECT_EQ(fill.streamAmounts[index], 0.0) << index;

    const SwingAdsorptionPhaseRecord& hold = swing.phaseRecords[1];
    for (const double passed : hold.streamAmounts) EXPECT_NEAR(passed, 0.0, 1.0e-12 * amount);
    EXPECT_NEAR(hold.gasAtEnd[0], hold.gasAtStart[0], 1.0e-9 * amount);

    const SwingAdsorptionPhaseRecord& empty = swing.phaseRecords[2];
    EXPECT_NEAR(-amountChange(empty)[0], amount, 2.0e-3 * amount);
    EXPECT_NEAR(empty.streamAmounts[3 * n], -amountChange(empty)[0], streamTolerance);
    EXPECT_NEAR(empty.pressureAtLength, 1.0e5, 500.0);
    EXPECT_NEAR(column.totalPressure[1], column.totalPressure[column.numberOfGridPoints], 50.0);
  }
}

TEST(MomentumBalance, RunsTheExampleCycleWithTheProductAsPurgeGasAndAClosedMassBalance)
{
  // examples/PSA-H2-activated-carbon-ports on 12 cells, for three cycles.
  std::ifstream input(std::filesystem::path{RUPTURA_SOURCE_DIR} /
                      "examples/PSA-H2-activated-carbon-ports/simulation.json");
  nlohmann::json data = nlohmann::json::parse(input);
  data["NumberOfGridPoints"] = 12;
  data["NumberOfCycles"] = 3;
  data["WriteInterval"] = 1000.0;
  const InputReader reader = readInput(data);
  WorkingDirectoryGuard workingDirectory(testTemporaryDirectory() / "ruptura_momentum_example");
  SwingAdsorption swing(reader);
  swing.run();
  ASSERT_EQ(swing.phaseRecords.size(), 12U);
  for (const double unaccounted : readColumn("cycles.data", 5)) EXPECT_LT(unaccounted, 1.0e-4);

  const size_t n = swing.breakthrough.column.numberOfComponents;
  const std::span<const SwingAdsorptionPhaseRecord> cycle(swing.phaseRecords.end() - 4, 4);
  EXPECT_NEAR(cycle[0].pressureAtLength, 1.0e6, 2.0e3);  // pressurized
  EXPECT_NEAR(cycle[2].pressureAtLength, 1.0e5, 2.0e3);  // blown down
  double product = 0.0;
  double purgeGas = 0.0;
  for (size_t comp = 0; comp < n; ++comp)
  {
    product += cycle[1].streamAmounts[3 * n + comp];
    purgeGas += cycle[3].streamAmounts[2 * n + comp];
  }
  ASSERT_GT(product, 0.0);
  EXPECT_GT(cycle[1].streamAmounts[3 * n + 3] / product, 0.9);  // hydrogen
  for (size_t comp = 0; comp < n; ++comp)
  {
    EXPECT_NEAR(cycle[3].streamAmounts[2 * n + comp] / purgeGas, cycle[1].streamAmounts[3 * n + comp] / product, 1.0e-6)
        << comp;
  }
}

TEST(GasHeatCapacity, OfTheComponentsGivesTheHeatCapacityPerMassItAddsUpTo)
{
  // Every component with M_j c_p per mole is the gas with c_p per kg, whatever its composition.
  nlohmann::json data = t15Input(8);
  data["BreakthroughIntegrator"] = "ROS34PW2";
  data["WriteInterval"] = 1000.0;
  data["SwingAdsorptionPhases"] = {{{"Name", "Adsorption"}, {"Duration", 200.0}}};
  SwingAdsorption perMass(readInput(data));
  perMass.run();
  for (auto& component : data["Components"])
  {
    component["GasHeatCapacity"] = component["MolecularWeight"].get<double>() * data["heatCapacityGas"].get<double>();
  }
  SwingAdsorption perMole(readInput(data));
  perMole.run();
  const Column& column = perMole.breakthrough.column;
  const Column& reference = perMass.breakthrough.column;
  ASSERT_EQ(column.gasHeatCapacity.size(), column.numberOfGridPoints + 1);
  EXPECT_TRUE(reference.gasHeatCapacity.empty());
  for (size_t node = 1; node <= column.numberOfGridPoints; ++node)
  {
    EXPECT_NEAR(column.gasHeatCapacity[node], reference.heatCapacityGas, 1.0e-9 * reference.heatCapacityGas) << node;
    EXPECT_NEAR(column.gasTemperature[node], reference.gasTemperature[node], 1.0e-6) << node;
    EXPECT_NEAR(column.solidTemperature[node], reference.solidTemperature[node], 1.0e-6) << node;
  }
}

TEST(CycleInput, ReadsPortsStopCriteriaAndSteadyStateInAnyCaseAndRejectsWhatDoesNotFit)
{
  nlohmann::json data = momentumBalanceInput(8, 1.0e5);
  data["NumberOfCycles"] = 3;
  data["cyclicsteadystate"] = {{"tolerance", 1.0e-4}, {"CONSECUTIVECYCLES", 2}};
  data["SwingAdsorptionPhases"] = {
      {{"Name", "Fill"},
       {"Duration", 30.0},
       {"minimumduration", 5.0},
       {"TIMEOUTPOLICY", "error"},
       {"stopcriterion",
        {{"any",
          {{{"metric", "pressure"}, {"end", "productend"}, {"operator", ">="}, {"value", 9.0e5}},
           {{"Metric", "AmountIn"}, {"End", "FeedEnd"}, {"Component", "co2"}, {"Operator", ">="}, {"Value", 50.0}}}}}},
       {"feedend",
        {{"mode", "pressure"}, {"PRESSURE", 1.0e6}, {"pressureramp", {{"shape", "linear"}, {"timeconstant", 20.0}}}}}},
      {{"Name", "Purge"},
       {"Duration", 10.0},
       {"ProductEndGasPhaseMolFractions", "Fill"},
       {"FeedEnd", {{"Mode", "Pressure"}, {"Pressure", 1.0e5}}},
       {"ProductEnd", {{"Mode", "Velocity"}, {"MolarFlux", 2.0}}}}};
  {
    const InputReader reader = readInput(data);
    EXPECT_TRUE(reader.streamAmounts);
    EXPECT_EQ(reader.numberOfCycles, 3U);
    EXPECT_EQ(reader.cyclicSteadyStateTolerance, 1.0e-4);
    EXPECT_EQ(reader.cyclicSteadyStateCycles, 2U);
    ASSERT_EQ(reader.swingAdsorptionPhases.size(), 6U);
    const InputReader::SwingAdsorptionPhase& fill = reader.swingAdsorptionPhases[2];
    EXPECT_TRUE(fill.stopOnAny);
    EXPECT_TRUE(fill.errorOnTimeout);
    EXPECT_EQ(fill.minimumDuration, 5.0);
    ASSERT_EQ(fill.stopCriteria.size(), 2U);
    EXPECT_EQ(fill.stopCriteria[0].metric, InputReader::StopCriterion::Metric::Pressure);
    EXPECT_EQ(fill.stopCriteria[0].end, 1U);
    EXPECT_EQ(fill.stopCriteria[1].metric, InputReader::StopCriterion::Metric::AmountIn);
    EXPECT_EQ(fill.stopCriteria[1].end, 0U);
    EXPECT_EQ(fill.stopCriteria[1].component, 1U);
    EXPECT_EQ(fill.ports[0].mode, 1U);
    EXPECT_EQ(fill.ports[0].pressureRampShape, 2U);
    EXPECT_EQ(fill.ports[1].mode, 0U);
    const InputReader::SwingAdsorptionPhase& purge = reader.swingAdsorptionPhases[3];
    EXPECT_EQ(purge.productEndFromPhase, 0U);
    EXPECT_EQ(purge.ports[1].molarFlux, 2.0);

    // A molar flux is a velocity of the port gas, here of the gas of the last cell.
    SwingAdsorption swing(reader);
    EXPECT_EQ(swing.phasesPerCycle, 2U);
    EXPECT_EQ(swing.steadyStateCycles, 2U);
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
  // The first phase with its keys as documented, and one of them changed.
  data["SwingAdsorptionPhases"][0] = {{"Name", "Fill"},
                                      {"Duration", 30.0},
                                      {"StopCriterion", {{"Metric", "Pressure"}, {"Operator", ">="}, {"Value", 9.0e5}}},
                                      {"FeedEnd", {{"Mode", "Pressure"}, {"Pressure", 1.0e6}}}};
  const auto withPhase = [&](const std::string& key, const nlohmann::json& value)
  {
    nlohmann::json input = data;
    input["SwingAdsorptionPhases"][0][key] = value;
    return input;
  };
  expectRejected(withPhase("InletPressure", 1.0e6), "the ports FeedEnd and ProductEnd set the ends");
  expectRejected(withPhase("ClosedEnd", "ProductEnd"), "the ports FeedEnd and ProductEnd set the ends");
  expectRejected(withPhase("FeedEnd", {{"Mode", "Pressure"}}), "a Pressure port a positive Pressure");
  expectRejected(withPhase("FeedEnd", {{"Mode", "Closed"}, {"Pressure", 1.0e5}}), "a Closed port takes no other key");
  expectRejected(withPhase("FeedEnd", {{"Mode", "Velocity"}, {"Velocity", 0.1}, {"MolarFlux", 2.0}}), "or a MolarFlux");
  expectRejected(withPhase("FeedEnd", {{"Mode", "Velocity"}, {"Velocity", 0.1}}), "the Pressure of its gas");
  expectRejected(withPhase("FeedEnd", {{"Mode", "Open"}}), "invalid Mode");
  expectRejected(withPhase("FeedEnd", {{"Pressure", 1.0e5}}), "Mode");
  expectRejected(withPhase("StopCriterion", {{"Metric", "Pressure"}, {"Value", 1.0}}), "Operator");
  expectRejected(withPhase("StopCriterion", {{"Metric", "Speed"}, {"Operator", ">="}, {"Value", 1.0}}),
                 "invalid Metric");
  expectRejected(withPhase("StopCriterion", {{"Metric", "MoleFraction"}, {"Operator", ">="}, {"Value", 0.1}}),
                 "needs a Component");
  expectRejected(withPhase("StopCriterion",
                           {{"Metric", "MoleFraction"}, {"Component", "Argon"}, {"Operator", ">="}, {"Value", 0.1}}),
                 "unknown component 'Argon'");
  expectRejected(withPhase("StopCriterion", {{"All", nlohmann::json::array()}}), "StopCriterion must be");
  expectRejected(withPhase("MinimumDuration", 40.0), "MinimumDuration lies between 0 and the Duration");
  expectRejected(withPhase("GasPhaseMolFractions", "Blowdown"), "names no phase 'Blowdown'");

  nlohmann::json input = data;
  input["SwingAdsorptionPhases"][1].erase("StopCriterion");
  input["SwingAdsorptionPhases"][1]["MinimumDuration"] = 2.0;
  expectRejected(input, "belong to a StopCriterion");
  input = data;
  input.erase("NumberOfCycles");
  expectRejected(input, "CyclicSteadyState needs");
  input = data;
  input["BoundaryCondition"] = "InletPressureInletVelocity";
  input["InletPressure"] = 1.0e5;
  expectRejected(input, "FeedEnd and ProductEnd need BoundaryCondition MomentumBalance");
  input = data;
  input["InletCondition"] = "Dirichlet";
  expectRejected(input, "MomentumBalance needs a gas, InletCondition Danckwerts");
  for (const auto& [integrator, solver] : {std::pair{"RungeKutta3", "Dense"}, std::pair{"CVODE", "Dense"},
                                           std::pair{"CVODE", "SPGMR"}})
  {
    input = data;
    input["BreakthroughIntegrator"] = integrator;
    input["CVODELinearSolver"] = solver;
    expectRejected(input, "MomentumBalance is stiff in the pressure");
  }
  input = data;
  input["Components"][1].erase("GasHeatCapacity");
  expectRejected(input, "to every component or to none");
}

namespace
{
// Columns of helium without the energy balance, named A and B: a port fills A to 0.5 MPa in 30 s, the others stay
// closed at 0.1 MPa.
nlohmann::json filledColumnInput(const nlohmann::json& columns)
{
  nlohmann::json data = momentumBalanceInput(10, 1.0e5);
  data["energyBalance"] = false;
  data["Components"][0]["GasPhaseMolFraction"] = 1.0;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["GasPhaseMolFraction"] = 0.0;
  }
  data["Columns"] = columns;
  data["SwingAdsorptionPhases"] =
      nlohmann::json::array({{{"Name", "Fill A"},
                              {"Duration", 30.0},
                              {"Columns",
                               {{"A",
                                 {{"FeedEnd",
                                   {{"Mode", "Pressure"},
                                    {"Pressure", 5.0e5},
                                    {"PressureRamp", {{"Shape", "Linear"}, {"TimeConstant", 20.0}}}}}}}}}}});
  return data;
}
}  // namespace

TEST(Connections, AValveEvensOutTwoInertColumnsAtItsRateAndWhatLeavesOneEntersTheOther)
{
  // A linear valve of conductance k between the product ends. Each column holds C = V / (R T) of gas per pressure, so
  // that their pressure difference falls as exp(-2 k t / C): B reaches 0.2 MPa, where the first phase stops, after
  // ln(2) C / (2 k), and both end at the mean, 0.3 MPa.
  const double conductance = 1.0e-8;
  nlohmann::json data = filledColumnInput(nlohmann::json::array({"A", "B"}));
  const nlohmann::json valve = nlohmann::json::array(
      {{{"Name", "Valve"}, {"From", "A.ProductEnd"}, {"To", "B.ProductEnd"}, {"Conductance", conductance}}});
  data["SwingAdsorptionPhases"].push_back(
      {{"Name", "Half"},
       {"Duration", 60.0},
       {"Connections", valve},
       {"StopCriterion", {{"Column", "B"}, {"Metric", "Pressure"}, {"Operator", ">="}, {"Value", 2.0e5}}}});
  data["SwingAdsorptionPhases"].push_back({{"Name", "Rest"}, {"Duration", 60.0}, {"Connections", valve}});
  SwingAdsorption swing(readInput(data));
  swing.run();
  ASSERT_EQ(swing.phaseRecords.size(), 6U);
  const Column& column = swing.breakthrough.column;
  const size_t n = column.numberOfComponents;
  double capacity = 0.0;
  for (size_t node = 1; node <= column.numberOfGridPoints; ++node)
  {
    capacity += column.crossSection() * column.totalVoidFraction[node] * gridSpacing(column.columnDistances, node) /
                (R * 298.0);
  }

  const SwingAdsorptionPhaseRecord& half = swing.phaseRecords[3];
  EXPECT_NEAR(half.endTime - half.startTime, std::log(2.0) * capacity / (2.0 * conductance), 0.02);
  EXPECT_NEAR(half.pressureAtLength, 2.0e5, 5.0);
  EXPECT_EQ(swing.phaseRecords[2].endTime, half.endTime);

  double passed = 0.0;
  for (size_t phase = 1; phase < 3; ++phase)
  {
    const SwingAdsorptionPhaseRecord& a = swing.phaseRecords[2 * phase];
    const SwingAdsorptionPhaseRecord& b = swing.phaseRecords[2 * phase + 1];
    const double left = a.streamAmounts[3 * n];
    passed += left;
    EXPECT_GT(left, 0.0);
    EXPECT_NEAR(b.streamAmounts[2 * n], left, 1.0e-12 * left);
    EXPECT_NEAR(amountChange(b)[0], left, 1.0e-9 * left);
    EXPECT_NEAR(-amountChange(a)[0], left, 1.0e-5 * left);
  }
  EXPECT_NEAR(passed * column.crossSection(), 2.0e5 * capacity, 1.0e-3 * 2.0e5 * capacity);
  EXPECT_NEAR(swing.phaseRecords[4].pressureAtLength, 3.0e5, 100.0);
  EXPECT_NEAR(swing.phaseRecords[5].pressureAtLength, 3.0e5, 100.0);

  // Each column writes to its own directory, and the connection its amounts in mol.
  EXPECT_EQ(readColumn("A/phases.data", 0).size(), 3U);
  EXPECT_TRUE(std::filesystem::exists("B/column.data"));
  const std::vector<double> carried = readColumn("connections.data", 2);
  ASSERT_EQ(carried.size(), 2U);
  EXPECT_NEAR(carried[0] + carried[1], passed * column.crossSection(), 1.0e-12 * carried[0]);
}

TEST(Connections, ABufferKeepsWhatAColumnGivesUpToItsPressureAndACompressorCostsTheWorkOfItsGas)
{
  // A square-root valve empties A into a tank that lets go of gas above 0.15 MPa, where both end. A compressor then
  // brings gas back at a set flow F: the pressures follow from the gas that column and tank hold per pressure, and the
  // work is the integral of F c_p T ((P_column / P_tank)^(R / c_p) - 1) / efficiency.
  const double volume = 1.0e-3;
  const double flow = 2.0e-4;
  nlohmann::json data = filledColumnInput(nlohmann::json::array({"A"}));
  data["Buffers"] = {{"Tank",
                      {{"Volume", volume},
                       {"Pressure", 1.0e5},
                       {"GasPhaseMolFractions", {{"Helium", 1.0}}},
                       {"MaximumPressure", 1.5e5}}}};
  data["SwingAdsorptionPhases"].push_back(
      {{"Name", "Vent"},
       {"Duration", 60.0},
       {"Connections",
        nlohmann::json::array(
            {{{"From", "A.ProductEnd"}, {"To", "Tank"}, {"FlowModel", "SquareRoot"}, {"Conductance", 8.0e-6}}})}});
  data["SwingAdsorptionPhases"].push_back({{"Name", "Compress"},
                                           {"Duration", 20.0},
                                           {"Connections", nlohmann::json::array({{{"From", "Tank"},
                                                                                   {"To", "A.ProductEnd"},
                                                                                   {"FlowModel", "Compressor"},
                                                                                   {"MolarFlow", flow},
                                                                                   {"Efficiency", 0.75}}})}});
  SwingAdsorption swing(readInput(data));
  swing.run();
  ASSERT_EQ(swing.phaseRecords.size(), 3U);
  const Column& column = swing.breakthrough.column;
  const size_t n = column.numberOfComponents;
  const double area = column.crossSection();
  const double tankCapacity = volume / (R * 298.0);
  const double columnCapacity = area * swing.phaseRecords[0].gasAtEnd[0] / 5.0e5;

  // The vent: the tank holds what its pressure allows, the rest is its export, and no gas is lost.
  const SwingAdsorptionPhaseRecord& vent = swing.phaseRecords[1];
  const std::vector<double> tankPressure = readColumn("buffers.data", 1);
  const std::vector<double> tankAmount = readColumn("buffers.data", 2);
  const std::vector<double> exported = readColumn("buffers.data", 2 + n);
  ASSERT_EQ(tankPressure.size(), 3U);
  EXPECT_NEAR(tankPressure[1], 1.5e5, 1.0);
  EXPECT_NEAR(vent.pressureAtLength, 1.5e5, 200.0);
  EXPECT_GT(exported[1], 0.0);
  EXPECT_NEAR(tankAmount[1] - tankAmount[0] + exported[1], area * vent.streamAmounts[3 * n], 1.0e-12 * tankAmount[1]);
  EXPECT_NEAR(tankAmount[1] - tankAmount[0] + exported[1], -area * amountChange(vent)[0], 1.0e-5 * tankAmount[1]);

  // The compressor.
  const SwingAdsorptionPhaseRecord& compress = swing.phaseRecords[2];
  EXPECT_NEAR(area * compress.streamAmounts[2 * n], 20.0 * flow, 1.0e-9 * flow);
  EXPECT_NEAR(tankAmount[1] - tankAmount[2], 20.0 * flow, 1.0e-9 * flow);
  EXPECT_NEAR(swing.bufferAmounts[0][0], tankAmount[2], 1.0e-12 * tankAmount[2]);
  double work = 0.0;
  for (size_t step = 0; step < 2000; ++step)
  {
    const double moved = (static_cast<double>(step) + 0.5) * 0.01 * flow;
    const double ratio = (vent.pressureAtLength + moved / columnCapacity) / (tankPressure[1] - moved / tankCapacity);
    work += 0.01 * flow * 20.8 * 298.0 * (std::pow(ratio, R / 20.8) - 1.0) / 0.75;
  }
  EXPECT_GT(work, 1.0);
  EXPECT_NEAR(readColumn("connections.data", 2 + n).back(), work, 0.01 * work);
}

TEST(Connections, TheTwoBedExampleKeepsEveryMoleBetweenItsColumnsItsTankAndItsProduct)
{
  // examples/PSA-H2-activated-carbon-two-beds on 12 cells, for two cycles of two halves: adsorption in one bed, whose
  // product end feeds the tank and the purge of the other bed, and an equalization.
  std::ifstream input(std::filesystem::path{RUPTURA_SOURCE_DIR} /
                      "examples/PSA-H2-activated-carbon-two-beds/simulation.json");
  nlohmann::json data = nlohmann::json::parse(input);
  data["NumberOfGridPoints"] = 12;
  data["NumberOfCycles"] = 2;
  data["WriteInterval"] = 1000.0;
  SwingAdsorption swing(readInput(data));
  const std::vector<double> tankAtStart = swing.bufferAmounts[0];
  swing.run();
  ASSERT_EQ(swing.phaseRecords.size(), 24U);
  const size_t n = swing.breakthrough.column.numberOfComponents;
  const double area = swing.breakthrough.column.crossSection();
  std::vector<double> product(n, 0.0);
  for (size_t comp = 0; comp < n; ++comp)
  {
    // Per half: the product valve, the purge valve, and the equalization valve.
    const std::vector<double> carried = readColumn("connections.data", 2 + comp);
    const std::vector<double> exported = readColumn("buffers.data", 2 + n + comp);
    ASSERT_EQ(carried.size(), 12U);
    ASSERT_EQ(exported.size(), 12U);
    const double tolerance = 1.0e-10 * area * swing.phaseRecords[2].streamAmounts[3 * n + n - 1];
    double toTank = 0.0;
    for (size_t half = 0; half < 4; ++half)
    {
      const SwingAdsorptionPhaseRecord& giving = swing.phaseRecords[2 * (3 * half + 1) + half % 2];
      const SwingAdsorptionPhaseRecord& taking = swing.phaseRecords[2 * (3 * half + 1) + 1 - half % 2];
      EXPECT_NEAR(area * giving.streamAmounts[3 * n + comp], carried[3 * half] + carried[3 * half + 1], tolerance);
      EXPECT_NEAR(area * taking.streamAmounts[2 * n + comp], carried[3 * half + 1], tolerance);
      toTank += carried[3 * half];
    }
    product[comp] = std::reduce(exported.begin() + 6, exported.end());
    EXPECT_NEAR(swing.bufferAmounts[0][comp] - tankAtStart[comp] + std::reduce(exported.begin(), exported.end()),
                toTank, tolerance);
  }
  EXPECT_GT(product[n - 1], 0.0);
  EXPECT_GT(product[n - 1] / std::reduce(product.begin(), product.end()), 0.9);
  for (const double unaccounted : readColumn("cycles.data", 5)) EXPECT_LT(unaccounted, 1.0e-4);
}

TEST(ThermalControl, AJacketHeatsAClosedColumnWhichGivesTheSameHeatBack)
{
  // Helium in a closed column with the energy balance: a jacket at 350 K heats wall and solid by 52 K, and one at the
  // phase temperature takes that heat out again.
  nlohmann::json data = momentumBalanceInput(8, 1.0e5);
  data["Components"][0]["GasPhaseMolFraction"] = 1.0;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["GasPhaseMolFraction"] = 0.0;
  }
  data["ExhaustPressure"] = 1.0e5;
  data["CompressionEfficiency"] = 0.7;
  data["SwingAdsorptionPhases"] = {
      {{"Name", "Heat"},
       {"Duration", 3000.0},
       {"ThermalControl", {{"Temperature", 350.0}, {"HeatTransferCoefficient", 200.0}}}},
      {{"Name", "Cool"}, {"Duration", 3000.0}, {"ThermalControl", {{"HeatTransferCoefficient", 200.0}}}}};
  SwingAdsorption swing(readInput(data));
  swing.run();
  ASSERT_EQ(swing.phaseRecords.size(), 2U);
  const Column& column = swing.breakthrough.column;
  const size_t n = column.numberOfComponents;

  // Wall and solid per unit column cross-section; the gas adds 1e-4 of that.
  const double wallArea = (0.04 * 0.04 - 0.025 * 0.025) / (0.025 * 0.025);
  const double heat = 52.0 * 1.2 * (4000.0 * 1000.0 * wallArea + (1.0 - 0.403) * 456.6 * 1000.0);
  EXPECT_NEAR(swing.phaseRecords[0].wallHeat, heat, 0.01 * heat);
  EXPECT_NEAR(swing.phaseRecords[1].wallHeat, -swing.phaseRecords[0].wallHeat, 1.0e-3 * heat);
  EXPECT_EQ(swing.phaseRecords[0].exhaustWork, 0.0);
  EXPECT_NEAR(readColumn("phases.data", 5 + 8 * n).front(), swing.phaseRecords[0].wallHeat, 1.0e-9 * heat);
}

TEST(ExhaustPressure, GasThatLeavesThroughAPortBelowItCostsTheWorkOfItsCompression)
{
  // A column of helium at 0.5 MPa and 298 K is opened to 0.02 MPa: every mole that leaves costs the isentropic
  // compression to the exhaust pressure of 0.1 MPa, over the efficiency.
  nlohmann::json data = filledColumnInput(nlohmann::json::array({"A"}));
  data["ExhaustPressure"] = 1.0e5;
  data["CompressionEfficiency"] = 0.7;
  data["SwingAdsorptionPhases"].push_back(
      {{"Name", "Evacuate"}, {"Duration", 60.0}, {"ProductEnd", {{"Mode", "Pressure"}, {"Pressure", 2.0e4}}}});
  SwingAdsorption swing(readInput(data));
  swing.run();
  ASSERT_EQ(swing.phaseRecords.size(), 2U);
  const SwingAdsorptionPhaseRecord& evacuate = swing.phaseRecords[1];
  const double left = evacuate.streamAmounts[3 * swing.breakthrough.column.numberOfComponents];
  const double work = left * 20.8 * 298.0 * (std::pow(5.0, R / 20.8) - 1.0) / 0.7;
  EXPECT_NEAR(left, 0.96 * evacuate.gasAtStart[0], 1.0e-3 * left);
  EXPECT_EQ(swing.phaseRecords[0].exhaustWork, 0.0);
  EXPECT_NEAR(evacuate.exhaustWork, work, 1.0e-3 * work);
}

TEST(CycleInput, ReadsColumnsBuffersAndConnectionsInAnyCaseAndRejectsWhatDoesNotFit)
{
  nlohmann::json data = filledColumnInput(nlohmann::json::array({"A", "B"}));
  data["buffers"] = {{"Tank", {{"volume", 1.0e-3}, {"PRESSURE", 2.0e5}, {"gasphasemolfractions", {{"helium", 2.0}}}}}};
  data["SwingAdsorptionPhases"].push_back(
      {{"Name", "Share"},
       {"Duration", 10.0},
       {"Temperature", 300.0},
       {"columns", {{"b", {{"temperature", 310.0}, {"feedend", {{"Mode", "Pressure"}, {"Pressure", 1.0e5}}}}}}},
       {"connections",
        nlohmann::json::array(
            {{{"from", "a.productend"}, {"to", "B.PRODUCTEND"}, {"flowmodel", "orifice"}, {"area", 1.0e-7}},
             {{"From", "tank"}, {"To", "A.FeedEnd"}, {"MolarFlow", 1.0e-4}, {"FlowModel", "Compressor"}}})}});
  {
    const InputReader reader = readInput(data);
    ASSERT_EQ(reader.columnNames.size(), 2U);
    ASSERT_EQ(reader.buffers.size(), 1U);
    EXPECT_EQ(reader.buffers[0].moleFractions[0], 1.0);
    EXPECT_EQ(reader.buffers[0].temperature, 298.0);
    EXPECT_EQ(reader.simulationTime, 40.0);
    // One entry per phase and column: the keys of the phase, with those of the column in their place.
    ASSERT_EQ(reader.swingAdsorptionPhases.size(), 4U);
    EXPECT_EQ(reader.swingAdsorptionPhases[0].ports[0].mode, 1U);
    EXPECT_EQ(reader.swingAdsorptionPhases[1].ports[0].mode, 0U);
    EXPECT_EQ(reader.swingAdsorptionPhases[2].temperature, 300.0);
    EXPECT_EQ(reader.swingAdsorptionPhases[3].temperature, 310.0);
    EXPECT_EQ(reader.swingAdsorptionPhases[3].ports[0].mode, 1U);
    EXPECT_TRUE(reader.swingAdsorptionPhases[3].connections.empty());
    const std::vector<InputReader::Connection>& connections = reader.swingAdsorptionPhases[2].connections;
    ASSERT_EQ(connections.size(), 2U);
    EXPECT_EQ(connections[0].name, "a.productend to B.PRODUCTEND");
    EXPECT_EQ(connections[0].endpoints, (std::array<size_t, 2>{1, 3}));
    EXPECT_EQ(connections[0].flowModel, InputReader::Connection::FlowModel::Orifice);
    EXPECT_EQ(connections[0].conductance, 1.0e-7);
    EXPECT_EQ(connections[1].endpoints, (std::array<size_t, 2>{4, 0}));
    EXPECT_EQ(SwingAdsorption(reader).phasesPerCycle, 2U);
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
  const auto withConnection = [&](const std::string& key, const nlohmann::json& value)
  {
    nlohmann::json input = data;
    input["SwingAdsorptionPhases"][1]["connections"][0][key] = value;
    return input;
  };
  expectRejected(withConnection("to", "B.FeedEnd"), "which then takes no port");
  expectRejected(withConnection("to", "C.ProductEnd"), "is no buffer and no <column>.FeedEnd");
  expectRejected(withConnection("to", "B"), "is no buffer and no <column>.FeedEnd");
  expectRejected(withConnection("to", "A.ProductEnd"), "two different endpoints");
  expectRejected(withConnection("Opening", 1.5), "an Opening from 0 to 1");
  expectRejected(withConnection("Efficiency", 0.8), "for a Compressor only");
  expectRejected(withConnection("Conductance", 1.0e-8), "a positive Area as its only size");
  expectRejected(withConnection("flowmodel", "Pipe"), "invalid FlowModel");

  nlohmann::json input = data;
  input["SwingAdsorptionPhases"][1]["columns"]["b"]["Duration"] = 5.0;
  expectRejected(input, "Duration");
  input = data;
  input["SwingAdsorptionPhases"][1]["columns"]["C"] = nlohmann::json::object();
  expectRejected(input, "names unknown column 'C'");
  input = data;
  input["Columns"] = nlohmann::json::array({"A", "a"});
  expectRejected(input, "distinct, non-empty names");
  input = data;
  input["buffers"]["Tank"]["MaximumPressure"] = 1.0e5;
  expectRejected(input, "a MaximumPressure of at least its Pressure");
  input = data;
  input["buffers"]["A"] = input["buffers"]["Tank"];
  expectRejected(input, "a name of its own");
  input = data;
  input["BreakthroughIntegrator"] = "CVODE";
  expectRejected(input, "Columns needs BoundaryCondition MomentumBalance and BreakthroughIntegrator ROS34PW2");
  input = data;
  input.erase("Columns");
  expectRejected(input, "Columns");
  input = data;
  input["SwingAdsorptionPhases"][0]["ThermalControl"] = {{"Temperature", 350.0}};
  expectRejected(input, "ThermalControl needs the energy balance");
  input = data;
  input["CompressionEfficiency"] = 0.0;
  expectRejected(input, "CompressionEfficiency lies in (0, 1]");
}

TEST(ParticleTransport, HsdmFillsASphereAtAConstantSurfaceLoadingAsItsSeriesDoes)
{
  // T15 with so little sorbent that the gas stays the feed, fed fast: every particle then sees a constant surface
  // loading q*, for which its loading follows q / q* = 1 - (6 / pi^2) sum_n exp(-n^2 pi^2 D_s t / R^2) / n^2. With the
  // linear driving force of the same coefficient, k = 15 D_s / R^2, it follows 1 - exp(-k t), which is slower at first.
  std::ifstream file(std::filesystem::path{RUPTURA_SOURCE_DIR} / "marx2/ternary/T15/simulation.json");
  nlohmann::json data = nlohmann::json::parse(file);
  data["NumberOfGridPoints"] = 6;
  data["energyBalance"] = false;
  data["ParticleDensity"] = 1.0e-6;
  data["ColumnEntranceVelocity"] = 0.5;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["MassTransferCoefficient"] = 0.015;
  }
  const double pi = 3.14159265358979323846;
  const auto loaded = [&](const std::string& integrator, double time, bool diffusion)
  {
    data["BreakthroughIntegrator"] = integrator;
    data.erase("ParticleTransport");
    if (diffusion) data["ParticleTransport"] = {{"model", "hsdm"}, {"RADIALCELLS", 24}};
    Breakthrough breakthrough(readInput(data));
    Timing timings;
    breakthrough.advanceTo(time, timings);
    const Column& column = breakthrough.column;
    const size_t n = column.numberOfComponents;
    const size_t index = 3 * n + 1;  // CO2 in the third cell
    // The loading of the particle is the mean of its radial cells over their volumes.
    double mean = 0.0;
    for (size_t shell = 0; shell < column.numberOfShells; ++shell)
    {
      const double lower = static_cast<double>(shell) / 24.0;
      const double upper = static_cast<double>(shell + 1) / 24.0;
      mean += (upper * upper * upper - lower * lower * lower) *
              column.shellLoading[shell * column.physisorption.size() + index];
    }
    if (diffusion) EXPECT_NEAR(mean, column.physisorption[index], 1.0e-3 * column.physisorption[index]);
    return column.physisorption[index] / column.equilibriumPhysisorption[index];
  };
  for (const double time : {50.0, 200.0})
  {
    double series = 0.0;
    for (int term = 1; term < 40; ++term) series += std::exp(-term * term * pi * pi * 1.0e-3 * time) / (term * term);
    EXPECT_NEAR(loaded("ROS34PW2", time, true), 1.0 - 6.0 / (pi * pi) * series, 0.01) << time;
    EXPECT_NEAR(loaded("CVODE", time, true), 1.0 - 6.0 / (pi * pi) * series, 0.01) << time;
    EXPECT_NEAR(loaded("ROS34PW2", time, false), 1.0 - std::exp(-0.015 * time), 0.01) << time;
  }

  data["BreakthroughIntegrator"] = "CVODE";
  data["CVODELinearSolver"] = "KLU";
  data["ParticleTransport"] = {{"Model", "HSDM"}};
  EXPECT_THROW(readInput(data), std::runtime_error);
  data["ParticleTransport"] = {{"Model", "HSDM"}, {"RadialCells", 1}};
  EXPECT_THROW(readInput(data), std::runtime_error);
}

TEST(Coolant, AFluidAlongTheWallCoolsOnItsWayAndHeatsTheColumnToItsInletTemperature)
{
  // Helium in a closed column, as in the jacket test, with water along the outer wall that enters at 350 K. After 20 s
  // the water has given heat to the wall on its way; counter-current it does so from the other end. In the end the
  // column is at 350 K, and the wall has taken the heat that its 52 K and those of the solid hold.
  nlohmann::json data = momentumBalanceInput(8, 1.0e5);
  data["Components"][0]["GasPhaseMolFraction"] = 1.0;
  for (size_t comp = 1; comp < data["Components"].size(); ++comp)
  {
    data["Components"][comp]["GasPhaseMolFraction"] = 0.0;
  }
  data["Coolant"] = {{"massflowrate", 0.002}, {"VOLUME", 2.0e-4}};
  const nlohmann::json control = {{"Temperature", 350.0}, {"HeatTransferCoefficient", 200.0}};
  data["SwingAdsorptionPhases"] = {{{"Name", "Start"}, {"Duration", 20.0}, {"ThermalControl", control}}};
  std::vector<std::vector<double>> profiles;
  for (const std::string direction : {"CoCurrent", "CounterCurrent"})
  {
    data["Coolant"]["FlowDirection"] = direction;
    SwingAdsorption swing(readInput(data));
    swing.run();
    const Column& column = swing.breakthrough.column;
    profiles.emplace_back(column.coolantTemperature.begin(), column.coolantTemperature.end());
    EXPECT_GT(swing.phaseRecords[0].wallHeat, 0.0);
  }
  const size_t last = profiles[0].size() - 1;
  for (size_t grid = 1; grid <= last; ++grid)
  {
    EXPECT_LT(profiles[0][grid], grid == 1 ? 350.0 : profiles[0][grid - 1]) << grid;
    EXPECT_NEAR(profiles[1][last + 1 - grid], profiles[0][grid], 0.5) << grid;
  }
  EXPECT_LT(profiles[0][last], 330.0);

  // A phase without flow leaves the water where it is; with flow again the column reaches the inlet temperature.
  data["SwingAdsorptionPhases"].push_back(
      {{"Name", "Hold"},
       {"Duration", 1.0},
       {"ThermalControl", {{"Temperature", 350.0}, {"HeatTransferCoefficient", 0.0}, {"MassFlowRate", 0.0}}}});
  data["SwingAdsorptionPhases"].push_back({{"Name", "Heat"}, {"Duration", 6000.0}, {"ThermalControl", control}});
  SwingAdsorption swing(readInput(data));
  swing.run();
  const Column& column = swing.breakthrough.column;
  for (size_t grid = 1; grid <= last; ++grid)
  {
    EXPECT_NEAR(column.coolantTemperature[grid], 350.0, 0.05) << grid;
    EXPECT_NEAR(column.solidTemperature[grid], 350.0, 0.05) << grid;
  }
  EXPECT_EQ(swing.phaseRecords[1].wallHeat, 0.0);
  const double wallArea = (0.04 * 0.04 - 0.025 * 0.025) / (0.025 * 0.025);
  const double heat = 52.0 * 1.2 * (4000.0 * 1000.0 * wallArea + (1.0 - 0.403) * 456.6 * 1000.0);
  EXPECT_NEAR(swing.phaseRecords[0].wallHeat + swing.phaseRecords[2].wallHeat, heat, 0.01 * heat);

  nlohmann::json input = data;
  input["energyBalance"] = false;
  EXPECT_THROW(readInput(input), std::runtime_error);
  input = data;
  input.erase("Coolant");
  EXPECT_THROW(readInput(input), std::runtime_error);  // the MassFlowRate of a phase without a Coolant
}

TEST(ROS34PW2, RecomputesTheStartOfARejectedStepFromItsOwnMixturePrediction)
{
  // B-ITG-23: the blowdown column of an experiment of the two-bed campaign (round 0, exp. 5), a column of helium at
  // rest at 1.68 bar with its feed end open at that pressure, into which traces of the feed disperse. A step whose
  // stages converged but failed the error test left the IAST with starting values from which that of the step start
  // did not converge, and the run ended with the IAST error.
  std::ifstream input(std::filesystem::path{RUPTURA_SOURCE_DIR} /
                      "examples/PSA-H2-activated-carbon-two-beds/simulation.json");
  nlohmann::json data = nlohmann::json::parse(input);
  data.erase("Columns");
  data.erase("Buffers");
  data.erase("CyclicSteadyState");
  const double temperature = 307.5483651012182;
  const double pressure = 167800.59710144997;
  const double duration = 20.71618204919959;
  data["Temperature"] = temperature;
  data["InfluxTemperature"] = temperature;
  data["OutletPressure"] = pressure;
  data["ColumnLength"] = 1.436670905082086;
  data["Geometry"]["InternalDiameter"] = 0.02284819956215412;
  data["Geometry"]["OuterDiameter"] = 0.03784819956215412;
  data["NumberOfGridPoints"] = 80;
  data["NumberOfCycles"] = 1;
  data["WriteInterval"] = 1.0e7;
  data["ROS34PW2RelativeTolerance"] = 1.0e-5;
  data["SwingAdsorptionPhases"] = {{{"Name", "Blowdown"},
                                    {"Temperature", temperature},
                                    {"Duration", duration},
                                    {"FeedEnd",
                                     {{"Mode", "Pressure"},
                                      {"Pressure", pressure},
                                      {"PressureRamp", {{"Shape", "Linear"}, {"TimeConstant", duration}}}}}}};
  SwingAdsorption swing(readInput(data));
  EXPECT_NO_THROW(swing.run());
  EXPECT_NEAR(swing.breakthrough.time, duration, 1.0e-9);
}

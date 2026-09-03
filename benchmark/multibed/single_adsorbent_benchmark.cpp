#include <benchmark/benchmark.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <string>
#include <utility>

#include "benchmark_cases.h"
#include "column.h"
#include "column_multibed.h"
#include "inputreader.h"
#include "integrators/rk3.h"
#include "mixture_prediction.h"
#include "timing.h"
#include "utils.h"

namespace
{
using multibed_benchmark::absoluteTolerance;
using multibed_benchmark::relativeTolerance;
using multibed_benchmark::Setting;
using multibed_benchmark::settings;
using multibed_benchmark::simulationSteps;

const Setting& currentSetting(const benchmark::State& state)
{
  return settings.at(static_cast<std::size_t>(state.range(0)));
}

template <typename ColumnType>
ColumnType makeInitializedColumn(const InputReader& input)
{
  ColumnType column(input);
  column.initialize();
  return column;
}

template <typename ColumnType>
void runSimulation(ColumnType& column, const InputReader& input)
{
  RungeKutta3 integrator(input.timeStep, false, simulationSteps);
  Timing timings;
  for (std::size_t step = 0; step < simulationSteps; ++step)
  {
    static_cast<void>(integrator.propagate(column, step, timings));
  }
}

struct Comparison
{
  double maximumAbsoluteError{0.0};
  double maximumRelativeError{0.0};
  bool equal{true};
};

Comparison compare(std::span<const double> singleBed, std::span<const double> multiBed)
{
  Comparison result;
  if (singleBed.size() != multiBed.size())
  {
    result.equal = false;
    result.maximumAbsoluteError = std::numeric_limits<double>::max();
    result.maximumRelativeError = std::numeric_limits<double>::max();
    return result;
  }

  for (std::size_t index = 0; index < singleBed.size(); ++index)
  {
    const double absoluteError = std::abs(singleBed[index] - multiBed[index]);
    const double scale = std::max({1.0, std::abs(singleBed[index]), std::abs(multiBed[index])});
    const double relativeError = absoluteError / scale;
    result.maximumAbsoluteError = std::max(result.maximumAbsoluteError, absoluteError);
    result.maximumRelativeError = std::max(result.maximumRelativeError, relativeError);
    result.equal = result.equal && (absoluteError <= absoluteTolerance || relativeError <= relativeTolerance);
  }
  return result;
}

void BM_SingleBed(benchmark::State& state)
{
  const Setting& setting = currentSetting(state);
  const InputReader input = multibed_benchmark::makeInput(setting);
  const Column initial = makeInitializedColumn<Column>(input);
  for (auto _ : state)
  {
    state.PauseTiming();
    Column column = initial;
    state.ResumeTiming();
    runSimulation(column, input);
    benchmark::DoNotOptimize(column.state.data());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(simulationSteps));
  state.SetLabel(multibed_benchmark::settingName(setting));
}

void BM_MultibedWithOneAdsorbent(benchmark::State& state)
{
  const Setting& setting = currentSetting(state);
  const InputReader input = multibed_benchmark::makeInput(setting);
  if (input.adsorbentComponents.size() != 1)
  {
    state.SkipWithError("The benchmark input must contain exactly one adsorbent");
    return;
  }
  const MultibedColumn initial = makeInitializedColumn<MultibedColumn>(input);
  for (auto _ : state)
  {
    state.PauseTiming();
    MultibedColumn column = initial;
    state.ResumeTiming();
    runSimulation(column, input);
    benchmark::DoNotOptimize(column.state.data());
    benchmark::ClobberMemory();
  }
  state.SetItemsProcessed(state.iterations() * static_cast<std::int64_t>(simulationSteps));
  state.SetLabel(multibed_benchmark::settingName(setting));
}

void BM_SingleBedVsMultibedEquality(benchmark::State& state)
{
  const Setting& setting = currentSetting(state);
  const InputReader input = multibed_benchmark::makeInput(setting);
  for (auto _ : state)
  {
    state.PauseTiming();
    Column singleBed = makeInitializedColumn<Column>(input);
    MultibedColumn multiBed = makeInitializedColumn<MultibedColumn>(input);
    state.ResumeTiming();

    runSimulation(singleBed, input);
    runSimulation(multiBed, input);

    state.PauseTiming();
    const Comparison stateComparison = compare(singleBed.state, multiBed.state);
    const Comparison pressureComparison = compare(singleBed.totalPressure, multiBed.totalPressure);
    const Comparison velocityComparison = compare(singleBed.interstitialGasVelocity, multiBed.interstitialGasVelocity);
    const double maximumAbsoluteError =
        std::max({stateComparison.maximumAbsoluteError, pressureComparison.maximumAbsoluteError,
                  velocityComparison.maximumAbsoluteError});
    const double maximumRelativeError =
        std::max({stateComparison.maximumRelativeError, pressureComparison.maximumRelativeError,
                  velocityComparison.maximumRelativeError});
    state.counters["max_abs_error"] = maximumAbsoluteError;
    state.counters["max_rel_error"] = maximumRelativeError;
    const bool resultsEqual = stateComparison.equal && pressureComparison.equal && velocityComparison.equal;
    state.counters["results_equal"] = resultsEqual ? 1.0 : 0.0;
    state.SetLabel(multibed_benchmark::settingName(setting) +
                   (resultsEqual ? "/within tolerance" : "/outside tolerance"));
    state.ResumeTiming();
  }
}
}  // namespace

BENCHMARK(BM_SingleBed)
    ->DenseRange(0, static_cast<std::int64_t>(settings.size() - 1))
    ->ArgName("setting")
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_MultibedWithOneAdsorbent)
    ->DenseRange(0, static_cast<std::int64_t>(settings.size() - 1))
    ->ArgName("setting")
    ->Unit(benchmark::kMillisecond);
BENCHMARK(BM_SingleBedVsMultibedEquality)
    ->DenseRange(0, static_cast<std::int64_t>(settings.size() - 1))
    ->ArgName("setting")
    ->Unit(benchmark::kMillisecond)
    ->Iterations(1);

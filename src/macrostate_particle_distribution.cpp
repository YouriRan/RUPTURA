#include "macrostate_particle_distribution.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
constexpr long double boltzmannConstant = 1.380649e-23L;  // J K^-1, exact SI definition
constexpr const char* noStateError = "Error: MPD reweighting produced no finite-probability macrostates";

std::string rowContext(const std::string& fileName, std::size_t lineNumber)
{
  return " (MPD file '" + fileName + "', line " + std::to_string(lineNumber) + ")";
}
}  // namespace

MacrostateParticleDistribution::MacrostateParticleDistribution(MPDSettings settings) : settings_(std::move(settings))
{
  if (settings_.fileName.empty()) throw std::runtime_error("Error: MPD FileName must not be empty");
  if (!std::isfinite(settings_.referenceTemperature) || settings_.referenceTemperature <= 0.0)
  {
    throw std::runtime_error("Error: MPD ReferenceTemperature must be a positive finite number");
  }
  if (!std::isfinite(settings_.referenceFugacity) || settings_.referenceFugacity <= 0.0)
  {
    throw std::runtime_error("Error: MPD ReferenceFugacity must be a positive finite number");
  }
  if (!std::isfinite(settings_.referenceFrameworkMass) || settings_.referenceFrameworkMass <= 0.0)
  {
    throw std::runtime_error("Error: MPD ReferenceFrameworkMass must be a positive finite number");
  }
  if (settings_.componentBounds.empty())
  {
    throw std::runtime_error("Error: MPD ComponentBounds must contain at least one component");
  }

  extents_.reserve(settings_.componentBounds.size());
  std::size_t numberOfStates = 1;
  for (const MPDComponentBounds& bounds : settings_.componentBounds)
  {
    if (bounds.deltaN == 0) throw std::runtime_error("Error: MPD DeltaN must be positive");
    if (bounds.nMax < bounds.nMin) throw std::runtime_error("Error: MPD NMax must be at least NMin");

    const std::size_t width = bounds.nMax - bounds.nMin;
    if (width % bounds.deltaN != 0)
    {
      throw std::runtime_error("Error: MPD [NMin, NMax] must be exactly divisible by DeltaN");
    }
    const std::size_t intervals = width / bounds.deltaN;
    if (intervals == std::numeric_limits<std::size_t>::max())
    {
      throw std::runtime_error("Error: MPD component extent overflows the addressable macrostate count");
    }
    const std::size_t extent = intervals + 1;
    if (numberOfStates > std::numeric_limits<std::size_t>::max() / extent)
    {
      throw std::runtime_error("Error: MPD dimensions overflow the addressable macrostate count");
    }
    numberOfStates *= extent;
    extents_.push_back(extent);
  }

  strides_.assign(extents_.size(), 1);
  for (std::size_t component = extents_.size(); component-- > 1;)
  {
    strides_[component - 1] = strides_[component] * extents_[component];
  }

  std::ifstream input(settings_.fileName);
  if (!input)
  {
    throw std::runtime_error("Error: MPD distribution file '" + settings_.fileName + "' could not be opened");
  }

  std::vector<double> probabilities;
  probabilities.reserve(numberOfStates);
  std::vector<double> parsedEnergies;
  parsedEnergies.reserve(numberOfStates);
  std::size_t expectedColumns = 0;
  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(input, line))
  {
    ++lineNumber;
    const std::size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;

    std::istringstream row(line);
    std::vector<double> values;
    double value = 0.0;
    while (row >> value) values.push_back(value);
    if (!row.eof())
    {
      throw std::runtime_error("Error: MPD rows must contain only numeric values" +
                               rowContext(settings_.fileName, lineNumber));
    }
    if (values.size() != 1 && values.size() != 2)
    {
      throw std::runtime_error("Error: MPD rows must contain Pi or Pi and <H>" +
                               rowContext(settings_.fileName, lineNumber));
    }
    if (expectedColumns == 0) expectedColumns = values.size();
    if (values.size() != expectedColumns)
    {
      throw std::runtime_error("Error: every MPD row must have the same number of columns" +
                               rowContext(settings_.fileName, lineNumber));
    }
    if (!std::isfinite(values[0]) || values[0] < 0.0 || values[0] > 1.0)
    {
      throw std::runtime_error("Error: MPD probabilities must be finite and in [0, 1]" +
                               rowContext(settings_.fileName, lineNumber));
    }
    if (values.size() == 2 && !std::isfinite(values[1]))
    {
      throw std::runtime_error("Error: MPD mean Hamiltonians must be finite" +
                               rowContext(settings_.fileName, lineNumber));
    }

    probabilities.push_back(values[0]);
    if (values.size() == 2) parsedEnergies.push_back(values[1]);
  }

  if (probabilities.size() != numberOfStates)
  {
    throw std::runtime_error("Error: MPD distribution file contains " + std::to_string(probabilities.size()) +
                             " entries, but ComponentBounds require " + std::to_string(numberOfStates));
  }
  bool hasPositiveProbability = false;
  for (double probability : probabilities) hasPositiveProbability = hasPositiveProbability || probability > 0.0;
  if (!hasPositiveProbability) throw std::runtime_error("Error: MPD distribution has zero total probability");

  // Reweighting only ever visits populated macrostates, so keep just those, with their particle numbers decoded
  // from the linear index once here instead of at every evaluation.
  numberOfMacrostates_ = numberOfStates;
  states_.particleNumbers.resize(rank());
  for (std::size_t state = 0; state < numberOfStates; ++state)
  {
    if (probabilities[state] == 0.0) continue;

    states_.logProbabilities.push_back(std::log(static_cast<long double>(probabilities[state])));
    for (std::size_t component = 0; component < rank(); ++component)
    {
      states_.particleNumbers[component].push_back(static_cast<double>(particleNumber(state, component)));
    }
    if (expectedColumns == 2) states_.meanEnergies.push_back(parsedEnergies[state]);
  }

  // A gas without component j reaches only the states with n_j = 0; keep them apart so it reweights just those.
  emptyComponentStates_.resize(rank());
  for (std::size_t empty = 0; empty < rank(); ++empty)
  {
    States& slice = emptyComponentStates_[empty];
    slice.particleNumbers.resize(rank());
    for (std::size_t state = 0; state < states_.logProbabilities.size(); ++state)
    {
      if (states_.particleNumbers[empty][state] != 0.0) continue;
      slice.logProbabilities.push_back(states_.logProbabilities[state]);
      for (std::size_t component = 0; component < rank(); ++component)
      {
        slice.particleNumbers[component].push_back(states_.particleNumbers[component][state]);
      }
      if (hasMeanEnergies()) slice.meanEnergies.push_back(states_.meanEnergies[state]);
    }
  }
}

std::size_t MacrostateParticleDistribution::particleNumber(std::size_t linearIndex,
                                                           std::size_t component) const noexcept
{
  const std::size_t coordinate = (linearIndex / strides_[component]) % extents_[component];
  const MPDComponentBounds& bounds = settings_.componentBounds[component];
  return bounds.nMin + coordinate * bounds.deltaN;
}

std::vector<double> MacrostateParticleDistribution::meanParticleNumbers(std::span<const double> gasMoleFractions,
                                                                        double fugacity, double temperature,
                                                                        std::span<double> covariance) const
{
  if (gasMoleFractions.size() != rank())
  {
    throw std::runtime_error("Error: MPD gas mole-fraction count does not match its rank");
  }
  if (!covariance.empty() && covariance.size() != rank() * rank())
  {
    throw std::runtime_error("Error: MPD covariance size does not match its rank squared");
  }
  std::fill(covariance.begin(), covariance.end(), 0.0);
  if (!std::isfinite(fugacity) || fugacity < 0.0)
  {
    throw std::runtime_error("Error: MPD target fugacity must be a non-negative finite number");
  }
  if (!std::isfinite(temperature) || temperature <= 0.0)
  {
    throw std::runtime_error("Error: MPD target temperature must be a positive finite number");
  }
  for (double fraction : gasMoleFractions)
  {
    if (!std::isfinite(fraction) || fraction < 0.0)
    {
      throw std::runtime_error("Error: MPD gas mole fractions must be non-negative finite numbers");
    }
  }

  std::vector<double> means(rank(), 0.0);
  // Without fugacity or without any adsorbing component in the gas, only the empty macrostate (file index 0) counts.
  if (fugacity == 0.0 || std::ranges::all_of(gasMoleFractions, [](double fraction) { return fraction == 0.0; }))
  {
    if (std::ranges::any_of(states_.particleNumbers, [](const auto& n) { return n.front() != 0.0; }))
    {
      throw std::runtime_error(noStateError);
    }
    return means;
  }

  const long double logFugacityRatio = std::log(static_cast<long double>(fugacity) / settings_.referenceFugacity);
  const long double logBetaRatio = std::log(static_cast<long double>(settings_.referenceTemperature) / temperature);
  const long double betaDifference = (1.0L / temperature - 1.0L / settings_.referenceTemperature) / boltzmannConstant;

  // An absent component excludes every state that holds it: reweight the smallest slice of an absent component.
  const States* states = &states_;
  for (std::size_t component = 0; component < rank(); ++component)
  {
    const States& slice = emptyComponentStates_[component];
    if (gasMoleFractions[component] == 0.0 && slice.logProbabilities.size() < states->logProbabilities.size())
    {
      states = &slice;
    }
  }
  const std::size_t count = states->logProbabilities.size();
  if (count == 0) throw std::runtime_error(noStateError);

  // Pass 1: log Pi_ref + sum_j n_j (log y_j + log f/f_ref + log T_ref/T) - <H> dbeta. One buffer per thread: this
  // runs ~1e8 times in a breakthrough.
  thread_local std::vector<long double> logWeights;
  logWeights.assign(states->logProbabilities.begin(), states->logProbabilities.end());
  for (std::size_t component = 0; component < rank(); ++component)
  {
    const std::vector<double>& n = states->particleNumbers[component];
    if (gasMoleFractions[component] > 0.0)
    {
      const long double slope =
          std::log(static_cast<long double>(gasMoleFractions[component])) + logFugacityRatio + logBetaRatio;
      for (std::size_t state = 0; state < count; ++state) logWeights[state] += n[state] * slope;
    }
    else
    {
      for (std::size_t state = 0; state < count; ++state)
      {
        if (n[state] != 0.0) logWeights[state] = -std::numeric_limits<long double>::infinity();
      }
    }
  }
  if (hasMeanEnergies() && betaDifference != 0.0L)
  {
    for (std::size_t state = 0; state < count; ++state) logWeights[state] -= states->meanEnergies[state] * betaDifference;
  }

  // Four running maxima: a single one is a serial compare-and-select chain, four times slower.
  long double maxima[4] = {-std::numeric_limits<long double>::infinity(), -std::numeric_limits<long double>::infinity(),
                           -std::numeric_limits<long double>::infinity(), -std::numeric_limits<long double>::infinity()};
  std::size_t state = 0;
  for (; state + 4 <= count; state += 4)
  {
    for (std::size_t lane = 0; lane < 4; ++lane) maxima[lane] = std::max(maxima[lane], logWeights[state + lane]);
  }
  for (; state < count; ++state) maxima[0] = std::max(maxima[0], logWeights[state]);
  const long double maximumLogWeight = std::max({maxima[0], maxima[1], maxima[2], maxima[3]});
  if (!std::isfinite(maximumLogWeight)) throw std::runtime_error(noStateError);

  // Pass 2: states more than 50 below the maximum add < e^-50 each to a sum of at least 1, below double rounding.
  long double normalization = 0.0L;
  std::vector<long double> weightedParticleNumbers(rank(), 0.0L);
  std::vector<long double> weightedProducts(covariance.size(), 0.0L);
  for (state = 0; state < count; ++state)
  {
    const long double difference = logWeights[state] - maximumLogWeight;
    if (difference < -50.0L) continue;
    const long double weight = std::exp(difference);
    normalization += weight;
    for (std::size_t component = 0; component < rank(); ++component)
    {
      weightedParticleNumbers[component] += weight * states->particleNumbers[component][state];
    }
    for (std::size_t a = 0; a < rank() && !covariance.empty(); ++a)
    {
      for (std::size_t b = 0; b <= a; ++b)
      {
        weightedProducts[a * rank() + b] +=
            weight * states->particleNumbers[a][state] * states->particleNumbers[b][state];
      }
    }
  }

  for (std::size_t component = 0; component < rank(); ++component)
  {
    means[component] = static_cast<double>(weightedParticleNumbers[component] / normalization);
  }
  for (std::size_t a = 0; a < rank() && !covariance.empty(); ++a)
  {
    for (std::size_t b = 0; b <= a; ++b)
    {
      const long double value = weightedProducts[a * rank() + b] / normalization -
                                weightedParticleNumbers[a] / normalization * weightedParticleNumbers[b] / normalization;
      covariance[a * rank() + b] = covariance[b * rank() + a] = static_cast<double>(value);
    }
  }
  return means;
}

#include "macrostate_particle_distribution.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace
{
constexpr long double boltzmannConstant = 1.380649e-23L;  // J K^-1, exact SI definition

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

  probabilities_.reserve(numberOfStates);
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

    probabilities_.push_back(values[0]);
    if (values.size() == 2) parsedEnergies.push_back(values[1]);
  }

  if (probabilities_.size() != numberOfStates)
  {
    throw std::runtime_error("Error: MPD distribution file contains " + std::to_string(probabilities_.size()) +
                             " entries, but ComponentBounds require " + std::to_string(numberOfStates));
  }
  bool hasPositiveProbability = false;
  for (double probability : probabilities_) hasPositiveProbability = hasPositiveProbability || probability > 0.0;
  if (!hasPositiveProbability) throw std::runtime_error("Error: MPD distribution has zero total probability");

  if (expectedColumns == 2) meanEnergies_ = std::move(parsedEnergies);
}

std::size_t MacrostateParticleDistribution::particleNumber(std::size_t linearIndex,
                                                           std::size_t component) const noexcept
{
  const std::size_t coordinate = (linearIndex / strides_[component]) % extents_[component];
  const MPDComponentBounds& bounds = settings_.componentBounds[component];
  return bounds.nMin + coordinate * bounds.deltaN;
}

std::vector<double> MacrostateParticleDistribution::meanParticleNumbers(std::span<const double> gasMoleFractions,
                                                                        double fugacity, double temperature) const
{
  if (gasMoleFractions.size() != rank())
  {
    throw std::runtime_error("Error: MPD gas mole-fraction count does not match its rank");
  }
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

  const long double logFugacityRatio =
      fugacity > 0.0 ? std::log(static_cast<long double>(fugacity) / settings_.referenceFugacity) : 0.0L;
  const long double logBetaRatio = std::log(static_cast<long double>(settings_.referenceTemperature) / temperature);
  const long double betaDifference = (1.0L / temperature - 1.0L / settings_.referenceTemperature) / boltzmannConstant;

  bool hasAccumulatedWeight = false;
  long double maximumLogWeight = 0.0L;
  long double normalization = 0.0L;
  std::vector<long double> weightedParticleNumbers(rank(), 0.0L);

  for (std::size_t state = 0; state < probabilities_.size(); ++state)
  {
    const double probability = probabilities_[state];
    if (probability == 0.0) continue;

    long double logWeight = std::log(static_cast<long double>(probability));
    long double totalParticleNumber = 0.0L;
    bool compatible = true;
    for (std::size_t component = 0; component < rank(); ++component)
    {
      const std::size_t n = particleNumber(state, component);
      totalParticleNumber += static_cast<long double>(n);
      if (n == 0) continue;
      if (gasMoleFractions[component] == 0.0)
      {
        compatible = false;
        break;
      }
      logWeight += static_cast<long double>(n) * std::log(static_cast<long double>(gasMoleFractions[component]));
    }
    if (!compatible || (fugacity == 0.0 && totalParticleNumber != 0.0L)) continue;

    logWeight += totalParticleNumber * (logFugacityRatio + logBetaRatio);
    if (hasMeanEnergies()) logWeight -= static_cast<long double>(meanEnergies_[state]) * betaDifference;

    if (!hasAccumulatedWeight || logWeight > maximumLogWeight)
    {
      const long double scale = hasAccumulatedWeight ? std::exp(maximumLogWeight - logWeight) : 0.0L;
      normalization = normalization * scale + 1.0L;
      for (std::size_t component = 0; component < rank(); ++component)
      {
        weightedParticleNumbers[component] =
            weightedParticleNumbers[component] * scale + particleNumber(state, component);
      }
      maximumLogWeight = logWeight;
      hasAccumulatedWeight = true;
    }
    else
    {
      const long double weight = std::exp(logWeight - maximumLogWeight);
      normalization += weight;
      for (std::size_t component = 0; component < rank(); ++component)
      {
        weightedParticleNumbers[component] += weight * particleNumber(state, component);
      }
    }
  }

  if (!hasAccumulatedWeight || !(normalization > 0.0L) || !std::isfinite(normalization))
  {
    throw std::runtime_error("Error: MPD reweighting produced no finite-probability macrostates");
  }

  std::vector<double> means(rank(), 0.0);
  for (std::size_t component = 0; component < rank(); ++component)
  {
    means[component] = static_cast<double>(weightedParticleNumbers[component] / normalization);
  }
  return means;
}

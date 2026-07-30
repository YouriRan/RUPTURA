#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

/**
 * \brief Universal gas constant in J/(mol K).
 */
const double R = 8.31446261815324;

/**
 * \brief Evaluates an Arrhenius law, optionally relative to a reference temperature.
 *
 * Returns k0 exp[-Ea/R (1/T - 1/Tref)]. Without Tref, the reference-temperature
 * term is zero and this reduces to the conventional k0 exp[-Ea/(RT)] form.
 */
inline double arrhenius(double preExponentialFactor, double activationEnergy, double temperature,
                        std::optional<double> referenceTemperature = std::nullopt)
{
  const double inverseTemperature = 1.0 / (R * std::max(temperature, 1.0e-10));
  const double inverseReferenceTemperature =
      referenceTemperature.has_value() ? 1.0 / (R * std::max(referenceTemperature.value(), 1.0e-10)) : 0.0;
  return preExponentialFactor *
         std::exp(-activationEnergy * (inverseTemperature - inverseReferenceTemperature));
}

/**
 * \brief Returns the maximum absolute element-wise difference between two vectors.
 */
inline double maxVectorDifference(const std::vector<double>& v, const std::vector<double>& w)
{
  if (v.empty() || w.empty()) return 0.0;
  if (v.size() != w.size()) throw std::runtime_error("Error: unequal vector size\n");

  double max = std::abs(v[0] - w[0]);
  for (size_t i = 1; i < v.size(); ++i)
  {
    double temp = std::abs(v[i] - w[i]);
    if (temp > max) max = temp;
  }
  return max;
}

/**
 * \brief Builds uniformly spaced node positions from 0 to columnLength.
 */
inline std::vector<double> makeUniformColumnDistances(size_t numberOfGridPoints, double columnLength)
{
  std::vector<double> distances(numberOfGridPoints + 1, 0.0);
  if (numberOfGridPoints == 0) return distances;

  const double dz = columnLength / static_cast<double>(numberOfGridPoints);
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    distances[grid] = static_cast<double>(grid) * dz;
  }
  distances.back() = columnLength;
  return distances;
}

/**
 * \brief Validates strictly increasing node positions from 0 to columnLength.
 */
inline void validateColumnDistances(const std::vector<double>& distances, size_t numberOfGridPoints,
                                    double columnLength, const std::string& name = "ColumnDistances")
{
  if (distances.size() != numberOfGridPoints + 1)
  {
    throw std::runtime_error("Error: " + name + " size must equal NumberOfGridPoints + 1");
  }
  if (distances.empty() || std::abs(distances.front()) > 1e-12)
  {
    throw std::runtime_error("Error: " + name + " must start at 0");
  }
  if (std::abs(distances.back() - columnLength) > std::max(1e-12, std::abs(columnLength) * 1e-10))
  {
    throw std::runtime_error("Error: " + name + " must end at ColumnLength");
  }
  for (size_t grid = 1; grid < distances.size(); ++grid)
  {
    if (distances[grid] <= distances[grid - 1])
    {
      throw std::runtime_error("Error: " + name + " must be strictly increasing");
    }
  }
}

/**
 * \brief Returns local node spacing between grid and grid - 1.
 */
inline double gridSpacing(std::span<const double> distances, size_t grid)
{
  if (distances.size() < 2) return 0.0;
  if (grid == 0) return distances[1] - distances[0];
  if (grid >= distances.size()) grid = distances.size() - 1;
  return distances[grid] - distances[grid - 1];
}

/**
 * \brief Returns the element-wise sum of two pairs.
 */
template <typename T, typename U>
std::pair<T, U> operator+(const std::pair<T, U>& l, const std::pair<T, U>& r)
{
  return {l.first + r.first, l.second + r.second};
}

/**
 * \brief Adds the second pair into the first pair element by element.
 */
template <typename T, typename U>
std::pair<T, U>& operator+=(std::pair<T, U>& l, const std::pair<T, U>& r)
{
  l.first += r.first;
  l.second += r.second;
  return l;
}

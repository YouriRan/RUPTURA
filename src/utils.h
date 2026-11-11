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
 * \brief Clamps values to the physically meaningful non-negative range.
 */
inline void clampNonnegative(std::span<double> values)
{
  for (double& value : values)
  {
    value = std::max(0.0, value);
  }
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

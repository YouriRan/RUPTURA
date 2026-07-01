#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

/**
 * \brief Universal gas constant in J/(mol K).
 */
const double R = 8.31446261815324;

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

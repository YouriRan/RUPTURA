#pragma once

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

/**
 * \brief Evaluates the dilogarithm function Li_2(x).
 */
extern double li2(double x);

/**
 * \brief Evaluates the Gauss hypergeometric function 2F1(a, b; c; z).
 */
extern double hypergeometric2F1(double a, double b, double c, double z);

/**
 * \brief Returns indices that sort a vector in ascending order without modifying the vector.
 */
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T>& v)
{
  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

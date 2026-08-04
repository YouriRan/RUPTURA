#pragma once

#include <cstddef>
#include <string>

#include "column_multibed.h"
#include "inputreader.h"
#include "rk3.h"
#include "timing.h"

/**
 * \brief Runs breakthrough simulations containing two or more adsorbent beds.
 */
struct BreakthroughMultibed
{
  explicit BreakthroughMultibed(const InputReader& inputReader);

  void print() const;
  [[nodiscard]] std::string repr() const;
  void run();

  const std::string displayName;
  size_t numberOfComponents;
  size_t numberOfGridPoints;
  size_t printEvery;
  size_t writeEvery;
  double timeStep;
  size_t numberOfInitTimeSteps;
  size_t numberOfSteps;
  bool autoNumberOfSteps;
  size_t maxIsothermTerms;

  ColumnMultibed column;
  RungeKutta3 rk3;
  Timing timings;
};

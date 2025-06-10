#pragma once

#include "breakthrough_state.h"
#include "compute.h"

struct RungeKutta3
{
  RungeKutta3(double timeStep, bool autoSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoSteps(autoSteps), numberOfSteps(numberOfSteps) {};

  double timeStep;
  bool autoSteps;
  size_t numberOfSteps;

  bool propagate(BreakthroughState& state, std::vector<Component>& components, MixturePrediction& mixture, size_t step);
};
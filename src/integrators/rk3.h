#pragma once

#include "column.h"
#include "compute.h"

struct RungeKutta3
{
  RungeKutta3(const InputReader& inputReader)
      : timeStep(inputReader.timeStep),
        autoSteps(inputReader.autoNumberOfTimeSteps),
        numberOfSteps(inputReader.numberOfTimeSteps) {};

  RungeKutta3(double timeStep, bool autoSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoSteps(autoSteps), numberOfSteps(numberOfSteps) {};

  double timeStep;
  bool autoSteps;
  size_t numberOfSteps;

  bool propagate(Column& column, size_t step);
};
#pragma once

#include "column.h"
#include "compute.h"

struct SemiImplicitRungeKutta3
{
  SemiImplicitRungeKutta3(double timeStep, bool autoSteps, size_t numberOfSteps)
      : timeStep(timeStep), autoSteps(autoSteps), numberOfSteps(numberOfSteps) {};

  double timeStep;
  bool autoSteps;
  size_t numberOfSteps;

  bool propagate(Column& state, size_t step);
};

void computePressureUpdateMatrix(Column& state, double timeStep, std::vector<double>& solved);
void computePressureUpdateMatrixFinal(Column& state, double timeStep, std::vector<double>& solved);
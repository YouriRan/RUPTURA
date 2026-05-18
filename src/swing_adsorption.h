#include <cstddef>
#include <tuple>
#include <vector>

#include "breakthrough.h"
#include "column.h"
#include "component.h"
#include "cvode.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "rk3.h"
#include "rk3_si.h"

struct SwingAdsorption
{
  struct SubStage
  {
    double temperature;
    double pressure;
    size_t numberOfSteps;
  };

  SwingAdsorption(const InputReader& inputReader);

  void run();
  void print() const;
  std::string repr() const;

  Breakthrough breakthrough;
  std::vector<SubStage> subStages;
};
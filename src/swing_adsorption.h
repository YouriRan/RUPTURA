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

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

struct SwingAdsorption
{
  SwingAdsorption(const InputReader& inputReader);

  void run();

  std::vector<std::pair<double, double>> pressTempPairs;
  std::vector<double> temperatures;
  std::vector<double> pressures;
  std::vector<size_t> steps;
}
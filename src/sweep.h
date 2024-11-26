#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "mixture_prediction.h"
namespace py = pybind11;

struct Sweep
{
  Sweep(std::string _displayName, std::vector<Component> _components, size_t _numberOfCarrierGases,
        size_t _carrierGasComponent, double _temperature, double _pressureStart, double _pressureEnd,
        size_t _numberOfPressurePoints, size_t _pressureScale, size_t _predictionMethod, size_t _iastMethod);

  MixturePrediction mix;
  py::array_t<double> phi;

  py::array_t<double> compute(py::array_t<double> phi);
}
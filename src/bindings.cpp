#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

#include "test.h"

PYBIND11_MODULE(ruptura, m) {
  py::class_<Test>(m, "Test")
    .def(py::init<size_t>())
    .def(py::init<std::vector<double>>())
    .def("getVector", &Test::getVector);
}


/*
#include "mixture_prediction.h"
#include "isotherm.h"
#include "component.h"
#include "multi_site_isotherm.h"
#include "inputreader.h"
#include "breakthrough.h"
PYBIND11_MODULE(ruptura, m) {

  py::class_<MultiSiteIsotherm>(m, "MultiSiteIsotherm");
  py::class_<Component>(m, "Component")
    .def(py::init<size_t, std::string>());
  py::class_<MixturePrediction>(m, "MixturePrediction")
    .def(py::init<std::string, std::vector<Component>, double, double, double, size_t, size_t, size_t, size_t>())
    .def("compute", &MixturePrediction::compute);
  py::class_<Isotherm>(m, "Isotherm")
    .def(py::init<size_t, std::vector<double>, size_t>());
  }
  py::class_<Breakthrough>(m, "Breakthrough")
    .def(py::init<std::string, std::vector<Component>, size_t, size_t, size_t, double, double, double, double, double, double, double, double, size_t, size_t, double, double, const MixturePrediction)
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

#include "mixture_prediction.h"
#include "isotherm.h"
#include "component.h"
#include "breakthrough.h"
#include "multi_site_isotherm.h"
#include "fitting.h"

PYBIND11_MODULE(_ruptura, m)
{
    py::class_<Component>(m, "Component")
        .def(py::init<size_t, std::string, std::vector<Isotherm>, double, double, double, bool>())
        .def_readonly("gasPhaseMolFraction", &Component::Yi0)
        .def_readonly("name", &Component::name)
        .def("__repr__", &Component::repr);
    py::class_<Isotherm>(m, "Isotherm")
        .def(py::init<size_t, std::vector<double>, size_t>())
        .def("__repr__", &Isotherm::repr);
    py::class_<MixturePrediction>(m, "MixturePrediction")
        .def(py::init<std::string, std::vector<Component>, size_t, size_t, double, double, double, size_t, size_t, size_t, size_t>())
        .def("getComponentsParameters", &MixturePrediction::getComponentsParameters)
        .def("setComponentsParameters", &MixturePrediction::setComponentsParameters)
        .def("setPressure", &MixturePrediction::setPressure)
        .def("compute", &MixturePrediction::compute)
        .def("__repr__", &MixturePrediction::repr);
    py::class_<Breakthrough>(m, "Breakthrough")
        .def(py::init<std::string, std::vector<Component>, size_t, size_t, size_t, size_t, double, double, double, double, double, double, double, double, size_t, bool, bool, double, const MixturePrediction>())
        .def("getComponentsParameters", &Breakthrough::getComponentsParameters)
        .def("setComponentsParameters", &Breakthrough::setComponentsParameters)
        .def("compute", &Breakthrough::compute)
        .def("__repr__", &Breakthrough::repr);
    py::class_<Fitting>(m, "Fitting")
        .def(py::init<std::string, std::vector<Component>, std::vector<std::vector<double>>, size_t>())
        .def("evaluate", &Fitting::evaluate)
        .def("compute", &Fitting::compute);
}

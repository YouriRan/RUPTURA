#pragma once

#include <vector>
#include <tuple>

#include "inputreader.h"
#include "component.h"

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>

// namespace py = pybind11;

class MixturePrediction
{
  public:
    enum class PredictionMethod
    {
      IAST = 0,
      SIAST = 1,
      EI = 2,
      SEI = 3
    };

    enum class IASTMethod
    {
      FastIAST = 0,
      NestedLoopBisection = 1
    };

    MixturePrediction(const InputReader &inputreader);
    MixturePrediction(
      std::string _displayName,
      std::vector <Component> _components,
      double _temperature,
      double _pressureStart,
      double _pressureEnd,
      size_t _numberOfPressurePoints,
      size_t _pressureScale,
      size_t _predictionMethod = 0,
      size_t _iastMethod = 0
    );

    void print() const;
    void run();
    void createPureComponentsPlotScript();
    void createMixturePlotScript();
    void createMixtureAdsorbedMolFractionPlotScript();
    void createPlotScript();
    std::vector<std::vector<std::vector<double>>> compute();

    // Yi  = gas phase mol-fraction
    // P   = total pressure
    // Xi  = adsorbed phase mol-fraction
    // Ni  = number of adsorbed molecules of component i
    std::pair<size_t, size_t> predictMixture(const std::vector<double> &Yi,
                                             const double &P,
                                             std::vector<double> &Xi,
                                             std::vector<double> &Ni,
                                             double *cachedP0,
                                             double *cachedPsi);

  private:
    std::string displayName;
    const std::vector<Component> components;
    std::vector<Component> sortedComponents;
    const size_t Ncomp;
    const size_t Nsorted;
    size_t numberOfCarrierGases;
    size_t carrierGasComponent;
    PredictionMethod predictionMethod;
    IASTMethod iastMethod;
    size_t maxIsothermTerms;
    std::vector<std::vector<Component>> segregatedSortedComponents;

    std::vector<double> alpha1;
    std::vector<double> alpha2;
    std::vector<double> alpha_prod;
    std::vector<double> x;

    std::vector<double> pstar;
    std::vector<double> psi;
    std::vector<double> G;
    std::vector<double> delta;
    std::vector<double> Phi;

    enum class PressureScale
    {
      Log = 0,
      Normal = 1
    };
    double temperature{ 300.0 };
    double pressureStart{ 1e3 };
    double pressureEnd{ 1e8 };
    size_t numberOfPressurePoints{ 100 };
    PressureScale pressureScale{ PressureScale::Log };

    std::pair<size_t, size_t> computeFastIAST(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeFastSIAST(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeFastSIAST(size_t term,
                                      const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);

    std::pair<size_t, size_t> computeIASTNestedLoopBisection(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeSIASTNestedLoopBisection(const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeSIASTNestedLoopBisection(size_t term,
                                      const std::vector<double> &Yi,
                                      const double &P,
                                      std::vector<double> &Xi,
                                      std::vector<double> &Ni,
                                      double *cachedP0,
                                      double *cachedPsi);
    std::pair<size_t, size_t> computeExplicitIsotherm(const std::vector<double> &Yi,
                                                      const double &P,
                                                      std::vector<double> &Xi,
                                                      std::vector<double> &Ni);
    std::pair<size_t, size_t> computeSegratedExplicitIsotherm(const std::vector<double> &Yi,
                                                      const double &P,
                                                      std::vector<double> &Xi,
                                                      std::vector<double> &Ni);
    std::pair<size_t, size_t> computeSegratedExplicitIsotherm(size_t site, const std::vector<double> &Yi,
                                                      const double &P,
                                                      std::vector<double> &Xi,
                                                      std::vector<double> &Ni);

    void printErrorStatus(double psi, double sum, double P, const std::vector<double> Yi, double cachedP0[]);
};

/* 
PYBIND11_MODULE(ruptura, m) {
  py::class_<MixturePrediction>(m, "MixturePrediction")
    .def(py::init<std::string, std::vector<Component>, double, double, double, size_t, size_t, size_t, size_t>())
    .def("compute", &MixturePrediction::compute)
}
*/
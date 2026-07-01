#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "breakthrough.h"
#include "column.h"
#include "component.h"
#include "compute.h"
#include "cvode.h"
#include "fitting.h"
#include "inputreader.h"
#include "isotherm.h"
#include "mixture_prediction.h"
#include "multi_site_isotherm.h"
#include "rk3.h"
#include "rk3_si.h"
#include "swing_adsorption.h"

namespace nb = nanobind;

#define EXPAND_MODULE(name) NB_MODULE(name, m)

namespace
{
template <typename T>
std::vector<T> spanToVector(std::span<T> values)
{
  return std::vector<T>(values.begin(), values.end());
}

template <typename T>
std::vector<T> spanToVector(std::span<const T> values)
{
  return std::vector<T>(values.begin(), values.end());
}
}  // namespace

EXPAND_MODULE(MODULE_NAME)
{
  nb::class_<InputReader> inputReader(m, "InputReader");

  nb::enum_<InputReader::SimulationType>(inputReader, "SimulationType", nb::is_arithmetic())
      .value("BREAKTHROUGH", InputReader::SimulationType::Breakthrough)
      .value("MIXTURE_PREDICTION", InputReader::SimulationType::MixturePrediction)
      .value("FITTING", InputReader::SimulationType::Fitting)
      .value("SWING_ADSORPTION", InputReader::SimulationType::SwingAdsorption)
      .value("TEST", InputReader::SimulationType::Test);

  inputReader.def(nb::init<const std::string>(), nb::arg("file_name"))
      .def_rw("components", &InputReader::components)
      .def_rw("number_of_carrier_gases", &InputReader::numberOfCarrierGases)
      .def_rw("carrier_gas_component", &InputReader::carrierGasComponent)
      .def_rw("max_isotherm_terms", &InputReader::maxIsothermTerms)
      .def_rw("simulation_type", &InputReader::simulationType)
      .def_rw("mixture_prediction_method", &InputReader::mixturePredictionMethod)
      .def_rw("iast_method", &InputReader::IASTMethod)
      .def_rw("breakthrough_integrator", &InputReader::breakthroughIntegrator)
      .def_rw("velocity_profile", &InputReader::velocityProfile)
      .def_rw("boundary_condition", &InputReader::boundaryCondition)
      .def_rw("display_name", &InputReader::displayName)
      .def_rw("temperature", &InputReader::temperature)
      .def_rw("column_void_fraction", &InputReader::columnVoidFraction)
      .def_rw("dynamic_viscosity", &InputReader::dynamicViscosity)
      .def_rw("particle_diameter", &InputReader::particleDiameter)
      .def_rw("particle_density", &InputReader::particleDensity)
      .def_rw("inlet_pressure", &InputReader::inletPressure)
      .def_rw("outlet_pressure", &InputReader::outletPressure)
      .def_rw("pressure_gradient", &InputReader::pressureGradient)
      .def_rw("column_entrance_velocity", &InputReader::columnEntranceVelocity)
      .def_rw("column_length", &InputReader::columnLength)
      .def_rw("influx_temperature", &InputReader::influxTemperature)
      .def_rw("internal_diameter", &InputReader::internalDiameter)
      .def_rw("outer_diameter", &InputReader::outerDiameter)
      .def_rw("wall_density", &InputReader::wallDensity)
      .def_rw("gas_thermal_conductivity", &InputReader::gasThermalConductivity)
      .def_rw("wall_thermal_conductivity", &InputReader::wallThermalConductivity)
      .def_rw("heat_transfer_gas_solid", &InputReader::heatTransferGasSolid)
      .def_rw("heat_transfer_gas_wall", &InputReader::heatTransferGasWall)
      .def_rw("heat_transfer_wall_external", &InputReader::heatTransferWallExternal)
      .def_rw("heat_capacity_gas", &InputReader::heatCapacityGas)
      .def_rw("heat_capacity_solid", &InputReader::heatCapacitySolid)
      .def_rw("heat_capacity_wall", &InputReader::heatCapacityWall)
      .def_rw("energy_balance", &InputReader::energyBalance)
      .def_rw("number_of_time_steps", &InputReader::numberOfTimeSteps)
      .def_rw("number_of_init_time_steps", &InputReader::numberOfInitTimeSteps)
      .def_rw("auto_number_of_time_steps", &InputReader::autoNumberOfTimeSteps)
      .def_rw("time_step", &InputReader::timeStep)
      .def_rw("print_every", &InputReader::printEvery)
      .def_rw("write_every", &InputReader::writeEvery)
      .def_rw("number_of_grid_points", &InputReader::numberOfGridPoints)
      .def_rw("pressure_start", &InputReader::pressureStart)
      .def_rw("pressure_end", &InputReader::pressureEnd)
      .def_rw("number_of_pressure_points", &InputReader::numberOfPressurePoints)
      .def_rw("pressure_scale", &InputReader::pressureScale)
      .def_rw("column_pressure", &InputReader::columnPressure)
      .def_rw("column_loading", &InputReader::columnLoading)
      .def_rw("column_error", &InputReader::columnError)
      .def_rw("swing_temperatures", &InputReader::swingTemperatures)
      .def_rw("swing_pressures", &InputReader::swingPressures)
      .def_rw("swing_steps", &InputReader::swingSteps)
      .def_rw("read_column_file", &InputReader::readColumnFile);

  nb::class_<Isotherm> isotherm(m, "Isotherm");

  nb::enum_<Isotherm::Type>(isotherm, "Type", nb::is_arithmetic())
      .value("LANGMUIR", Isotherm::Type::Langmuir)
      .value("ANTI_LANGMUIR", Isotherm::Type::Anti_Langmuir)
      .value("BET", Isotherm::Type::BET)
      .value("HENRY", Isotherm::Type::Henry)
      .value("FREUNDLICH", Isotherm::Type::Freundlich)
      .value("SIPS", Isotherm::Type::Sips)
      .value("LANGMUIR_FREUNDLICH", Isotherm::Type::Langmuir_Freundlich)
      .value("REDLICH_PETERSON", Isotherm::Type::Redlich_Peterson)
      .value("TOTH", Isotherm::Type::Toth)
      .value("UNILAN", Isotherm::Type::Unilan)
      .value("OBRIEN_MYERS", Isotherm::Type::OBrien_Myers)
      .value("QUADRATIC", Isotherm::Type::Quadratic)
      .value("TEMKIN", Isotherm::Type::Temkin)
      .value("BINGEL_WALTON", Isotherm::Type::BingelWalton);

  isotherm
      .def(nb::init<Isotherm::Type, std::vector<double>, bool>(), nb::arg("type"), nb::arg("parameters"),
           nb::arg("non_isothermal") = false)
      .def(nb::init<Isotherm::Type, const std::vector<double>&, std::size_t, bool>(), nb::arg("type"),
           nb::arg("parameters"), nb::arg("number_of_parameters"), nb::arg("non_isothermal") = false)
      .def(nb::init<std::size_t, const std::vector<double>&, std::size_t, bool>(), nb::arg("type"),
           nb::arg("parameters"), nb::arg("number_of_parameters"), nb::arg("non_isothermal") = false)
      .def_rw("type", &Isotherm::type)
      .def_rw("parameters", &Isotherm::parameters)
      .def_rw("number_of_parameters", &Isotherm::numberOfParameters)
      .def_rw("non_isothermal", &Isotherm::nonIsothermal)
      .def("print", &Isotherm::print)
      .def("__repr__", &Isotherm::repr)
      .def("value", &Isotherm::value, nb::arg("pressure"), nb::arg("scale") = 1.0)
      .def("psi_for_pressure", &Isotherm::psiForPressure, nb::arg("pressure"), nb::arg("scale") = 1.0)
      .def(
          "inverse_pressure_for_psi",
          [](const Isotherm& self, double reducedGrandPotential, double cachedPressure, double scale)
          {
            double cachedPressureValue = cachedPressure;
            double inversePressure = self.inversePressureForPsi(reducedGrandPotential, cachedPressureValue, scale);
            return nb::make_tuple(inversePressure, cachedPressureValue);
          },
          nb::arg("reduced_grand_potential"), nb::arg("cached_pressure") = 0.0, nb::arg("scale") = 1.0)
      .def("randomize", &Isotherm::randomize, nb::arg("maximum_loading"))
      .def("is_unphysical", &Isotherm::isUnphysical);

  nb::class_<MultiSiteIsotherm>(m, "MultiSiteIsotherm")
      .def(nb::init<>())
      .def(nb::init<std::vector<Isotherm>>(), nb::arg("sites"))
      .def_rw("sites", &MultiSiteIsotherm::sites)
      .def_rw("number_of_sites", &MultiSiteIsotherm::numberOfSites)
      .def_rw("number_of_parameters", &MultiSiteIsotherm::numberOfParameters)
      .def("add", &MultiSiteIsotherm::add, nb::arg("site"))
      .def("print", &MultiSiteIsotherm::print)
      .def("__repr__", &MultiSiteIsotherm::repr)
      .def(
          "value", [](const MultiSiteIsotherm& self, double pressure, double scale)
          { return self.value(pressure, scale); }, nb::arg("pressure"), nb::arg("scale") = 1.0)
      .def(
          "value_site", [](const MultiSiteIsotherm& self, std::size_t site, double pressure, double scale)
          { return self.value(site, pressure, scale); }, nb::arg("site"), nb::arg("pressure"), nb::arg("scale") = 1.0)
      .def(
          "psi_for_pressure", [](const MultiSiteIsotherm& self, double pressure, double scale)
          { return self.psiForPressure(pressure, scale); }, nb::arg("pressure"), nb::arg("scale") = 1.0)
      .def(
          "psi_for_pressure_site", [](const MultiSiteIsotherm& self, std::size_t site, double pressure, double scale)
          { return self.psiForPressure(site, pressure, scale); }, nb::arg("site"), nb::arg("pressure"),
          nb::arg("scale") = 1.0)
      .def(
          "inverse_pressure_for_psi",
          [](const MultiSiteIsotherm& self, double reducedGrandPotential, double cachedPressure, double scale)
          {
            double cachedPressureValue = cachedPressure;
            double inversePressure = self.inversePressureForPsi(reducedGrandPotential, cachedPressureValue, scale);
            return nb::make_tuple(inversePressure, cachedPressureValue);
          },
          nb::arg("reduced_grand_potential"), nb::arg("cached_pressure") = 0.0, nb::arg("scale") = 1.0)
      .def(
          "inverse_pressure_for_psi_site",
          [](const MultiSiteIsotherm& self, std::size_t site, double reducedGrandPotential, double cachedPressure,
             double scale)
          {
            double cachedPressureValue = cachedPressure;
            double inversePressure =
                self.inversePressureForPsi(site, reducedGrandPotential, cachedPressureValue, scale);
            return nb::make_tuple(inversePressure, cachedPressureValue);
          },
          nb::arg("site"), nb::arg("reduced_grand_potential"), nb::arg("cached_pressure") = 0.0, nb::arg("scale") = 1.0)
      .def("randomized", &MultiSiteIsotherm::randomized, nb::arg("maximum_loading"))
      .def("fitness", &MultiSiteIsotherm::fitness)
      .def("get_parameters", &MultiSiteIsotherm::getParameters)
      .def("set_parameters", &MultiSiteIsotherm::setParameters, nb::arg("parameters"));

  nb::class_<Component>(m, "Component")
      .def(nb::init<std::size_t, std::string>(), nb::arg("id"), nb::arg("name"))
      .def(nb::init<std::size_t, std::string, std::vector<Isotherm>, double, double, double, bool, double, double>(),
           nb::arg("id"), nb::arg("name"), nb::arg("isotherms"), nb::arg("yi0"), nb::arg("kl"), nb::arg("diffusion"),
           nb::arg("is_carrier_gas") = false, nb::arg("molecular_weight") = 1.0, nb::arg("heat_of_adsorption") = 0.0)
      .def(nb::init<std::size_t, std::string, MultiSiteIsotherm, double, double, double, bool, double, double>(),
           nb::arg("id"), nb::arg("name"), nb::arg("isotherm"), nb::arg("yi0"), nb::arg("kl"), nb::arg("diffusion"),
           nb::arg("is_carrier_gas") = false, nb::arg("molecular_weight") = 1.0, nb::arg("heat_of_adsorption") = 0.0)
      .def_rw("id", &Component::id)
      .def_rw("name", &Component::name)
      .def_rw("filename", &Component::filename)
      .def_rw("isotherm", &Component::isotherm)
      .def_rw("yi0", &Component::initialGasMoleFraction)
      .def_rw("kl", &Component::massTransferCoefficient)
      .def_rw("diffusion", &Component::axialDispersionCoefficient)
      .def_rw("heat_of_adsorption", &Component::heatOfAdsorption)
      .def_rw("is_carrier_gas", &Component::isCarrierGas)
      .def_rw("molecular_weight", &Component::molecularWeight)
      .def_rw("non_isothermal", &Component::nonIsothermal)
      .def_rw("reference_temperature", &Component::referenceTemperature)
      .def(
          "scale",
          [](const Component& self, double temperature)
          {
            double temperatureValue = temperature;
            return self.scale(temperatureValue);
          },
          nb::arg("temperature"))
      .def("print", &Component::print)
      .def("__repr__", &Component::repr);

  nb::class_<MixturePrediction> mixturePrediction(m, "MixturePrediction");

  nb::enum_<MixturePrediction::PredictionMethod>(mixturePrediction, "PredictionMethod", nb::is_arithmetic())
      .value("IAST", MixturePrediction::PredictionMethod::IAST)
      .value("SIAST", MixturePrediction::PredictionMethod::SIAST)
      .value("EI", MixturePrediction::PredictionMethod::EI)
      .value("SEI", MixturePrediction::PredictionMethod::SEI);

  nb::enum_<MixturePrediction::IASTMethod>(mixturePrediction, "IASTMethod", nb::is_arithmetic())
      .value("FAST_IAST", MixturePrediction::IASTMethod::FastIAST)
      .value("NESTED_LOOP_BISECTION", MixturePrediction::IASTMethod::NestedLoopBisection);

  nb::enum_<MixturePrediction::PressureScale>(mixturePrediction, "PressureScale", nb::is_arithmetic())
      .value("LOG", MixturePrediction::PressureScale::Log)
      .value("NORMAL", MixturePrediction::PressureScale::Normal);

  mixturePrediction.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(nb::init<std::string, std::vector<Component>, std::size_t, std::size_t, MixturePrediction::PredictionMethod,
                    MixturePrediction::IASTMethod, std::size_t, double, double, double, std::size_t,
                    MixturePrediction::PressureScale>(),
           nb::arg("display_name"), nb::arg("components"), nb::arg("number_of_carrier_gases"),
           nb::arg("carrier_gas_component"), nb::arg("prediction_method") = MixturePrediction::PredictionMethod::IAST,
           nb::arg("iast_method") = MixturePrediction::IASTMethod::FastIAST, nb::arg("max_isotherm_terms"),
           nb::arg("temperature") = 300.0, nb::arg("pressure_start") = 1.0e3, nb::arg("pressure_end") = 1.0e8,
           nb::arg("number_of_pressure_points") = 100,
           nb::arg("pressure_scale") = MixturePrediction::PressureScale::Log)
      .def("print", &MixturePrediction::print)
      .def("__repr__", &MixturePrediction::repr)
      .def_rw("display_name", &MixturePrediction::displayName)
      .def_rw("components", &MixturePrediction::components)
      .def_rw("sorted_components", &MixturePrediction::sortedComponents)
      .def_rw("number_of_components", &MixturePrediction::numberOfComponents)
      .def_rw("number_of_sorted_components", &MixturePrediction::numberOfSortedComponents)
      .def_rw("number_of_carrier_gases", &MixturePrediction::numberOfCarrierGases)
      .def_rw("carrier_gas_component", &MixturePrediction::carrierGasComponent)
      .def_rw("prediction_method", &MixturePrediction::predictionMethod)
      .def_rw("iast_method", &MixturePrediction::iastMethod)
      .def_rw("max_isotherm_terms", &MixturePrediction::maxIsothermTerms)
      .def_rw("segregated_sorted_components", &MixturePrediction::segregatedSortedComponents)
      .def_rw("firstExplicitIsothermAlpha", &MixturePrediction::firstExplicitIsothermAlpha)
      .def_rw("secondExplicitIsothermAlpha", &MixturePrediction::secondExplicitIsothermAlpha)
      .def_rw("explicitIsothermAlphaProduct", &MixturePrediction::explicitIsothermAlphaProduct)
      .def_rw("adsorbedMoleFractionsScratch", &MixturePrediction::adsorbedMoleFractionsScratch)
      .def_rw("hypotheticalPressure", &MixturePrediction::hypotheticalPressure)
      .def_rw("reducedGrandPotential", &MixturePrediction::reducedGrandPotential)
      .def_rw("g", &MixturePrediction::residualVector)
      .def_rw("correctionVector", &MixturePrediction::correctionVector)
      .def_rw("phi", &MixturePrediction::jacobianMatrix)
      .def_rw("temperature", &MixturePrediction::temperature)
      .def_rw("pressure_start", &MixturePrediction::pressureStart)
      .def_rw("pressure_end", &MixturePrediction::pressureEnd)
      .def_rw("number_of_pressure_points", &MixturePrediction::numberOfPressurePoints)
      .def_rw("pressure_scale", &MixturePrediction::pressureScale)
      .def("run", &MixturePrediction::run)
      .def("init_pressures", &MixturePrediction::initPressures)
      .def("sort_components", &MixturePrediction::sortComponents)
      .def("compute",
           [](MixturePrediction& self) -> nb::ndarray<nb::numpy, double, nb::ndim<3>>
           {
             // Based on run(), but returns a NumPy array with shape:
             // [numberOfPressurePoints, numberOfComponents, 6].
             std::vector<double> idealGasMolFractions(self.numberOfComponents);
             std::vector<double> adsorbedMolFractions(self.numberOfComponents);
             std::vector<double> numberOfMolecules(self.numberOfComponents);
             std::vector<double> cachedPressure(self.numberOfComponents * self.maxIsothermTerms);
             std::vector<double> cachedGrandPotential(self.maxIsothermTerms);

             for (size_t i = 0; i < self.numberOfComponents; ++i)
             {
               idealGasMolFractions[i] = self.components[i].initialGasMoleFraction;
             }

             std::vector<double> pressures = self.initPressures();

             auto* buffer = new std::vector<double>(self.numberOfPressurePoints * self.numberOfComponents * 6, 0.0);
             nb::capsule owner(buffer,
                               [](void* pointer) noexcept { delete static_cast<std::vector<double>*>(pointer); });

             double* data = buffer->data();
             for (size_t i = 0; i < self.numberOfPressurePoints; ++i)
             {
               // Check for an interrupt/error from the Python side.
               if (PyErr_CheckSignals() != 0)
               {
                 throw nb::python_error();
               }

               double gasTemperature = self.temperature;
               self.predictMixture(idealGasMolFractions, pressures[i], adsorbedMolFractions, numberOfMolecules,
                                   cachedPressure.data(), cachedGrandPotential.data(), gasTemperature);

               for (size_t j = 0; j < self.numberOfComponents; j++)
               {
                 double pStar = idealGasMolFractions[j] * pressures[i] / adsorbedMolFractions[j];
                 size_t index = (i * self.numberOfComponents + j) * 6;
                 data[index] = pressures[i];
                 data[index + 1] = self.components[j].isotherm.value(pressures[i], gasTemperature);
                 data[index + 2] = numberOfMolecules[j];
                 data[index + 3] = idealGasMolFractions[j];
                 data[index + 4] = adsorbedMolFractions[j];
                 data[index + 5] = self.components[j].isotherm.psiForPressure(pStar, gasTemperature);
               }
             }

             return nb::ndarray<nb::numpy, double, nb::ndim<3>>(
                 buffer->data(), {self.numberOfPressurePoints, self.numberOfComponents, static_cast<size_t>(6)}, owner);
           })
      .def("get_max_isotherm_terms", &MixturePrediction::getMaxIsothermTerms)
      .def(
          "predict_mixture",
          [](MixturePrediction& self, const std::vector<double>& idealGasMolFractions, double externalPressure,
             double gasTemperature, std::vector<double> cachedPressure, std::vector<double> cachedGrandPotential)
          {
            if (idealGasMolFractions.empty())
            {
              throw std::invalid_argument("ideal_gas_mol_fractions must not be empty");
            }

            const std::size_t numberOfComponents = idealGasMolFractions.size();
            const std::size_t maxIsothermTerms = self.getMaxIsothermTerms();

            if (cachedPressure.empty())
            {
              cachedPressure.assign(numberOfComponents * maxIsothermTerms, 0.0);
            }
            if (cachedGrandPotential.empty())
            {
              cachedGrandPotential.assign(maxIsothermTerms, 0.0);
            }

            if (cachedPressure.size() < numberOfComponents * maxIsothermTerms)
            {
              throw std::invalid_argument("cached_pressure is smaller than number_of_components * max_isotherm_terms");
            }
            if (cachedGrandPotential.size() < maxIsothermTerms)
            {
              throw std::invalid_argument("cached_grand_potential is smaller than max_isotherm_terms");
            }

            std::vector<double> adsorbedMolFractions(numberOfComponents, 0.0);
            std::vector<double> numberOfMolecules(numberOfComponents, 0.0);
            double gasTemperatureValue = gasTemperature;

            std::pair<std::size_t, std::size_t> performance = self.predictMixture(
                std::span<const double>(idealGasMolFractions.data(), idealGasMolFractions.size()), externalPressure,
                std::span<double>(adsorbedMolFractions.data(), adsorbedMolFractions.size()),
                std::span<double>(numberOfMolecules.data(), numberOfMolecules.size()), cachedPressure.data(),
                cachedGrandPotential.data(), gasTemperatureValue);

            return nb::make_tuple(performance, adsorbedMolFractions, numberOfMolecules, cachedPressure,
                                  cachedGrandPotential, gasTemperatureValue);
          },
          nb::arg("ideal_gas_mol_fractions"), nb::arg("external_pressure"), nb::arg("gas_temperature"),
          nb::arg("cached_pressure") = std::vector<double>(), nb::arg("cached_grand_potential") = std::vector<double>())
      .def(
          "set_pressure",
          [](MixturePrediction& self, double pressureStart, double pressureEnd)
          {
            self.pressureStart = pressureStart;
            self.pressureEnd = pressureEnd;
          },
          nb::arg("pressure_start"), nb::arg("pressure_end"))
      .def(
          "set_components_parameters",
          [](MixturePrediction& self, std::vector<double> molfracs, std::vector<double> params)
          {
            size_t index = 0;
            for (size_t i = 0; i < self.numberOfComponents; ++i)
            {
              self.components[i].initialGasMoleFraction = molfracs[i];
              size_t numberOfParameters = self.components[i].isotherm.parameters.size();
              std::vector<double> slicedVec(params.begin() + index, params.begin() + index + numberOfParameters);
              index += numberOfParameters;
              self.components[i].isotherm.setParameters(slicedVec);
            }
            self.sortedComponents = self.components;
            self.segregatedSortedComponents =
                std::vector<std::vector<Component>>(self.maxIsothermTerms, std::vector<Component>(self.components));
            self.sortComponents();
          },
          nb::arg("molfracs"), nb::arg("params"))
      .def("get_components_parameters",
           [](const MixturePrediction& self)
           {
             std::vector<double> params;
             for (size_t i = 0; i < self.numberOfComponents; ++i)
             {
               std::vector<double> compParams = self.components[i].isotherm.getParameters();
               params.insert(params.end(), compParams.begin(), compParams.end());
             }
             return params;
           });

  nb::class_<Column> column(m, "Column");

  nb::enum_<Column::VelocityProfile>(column, "VelocityProfile", nb::is_arithmetic())
      .value("FIXED_PRESSURE_GRADIENT", Column::VelocityProfile::FixedPressureGradient)
      .value("ERGUN", Column::VelocityProfile::Ergun)
      .value("FIXED_VELOCITY", Column::VelocityProfile::FixedVelocity);

  nb::enum_<Column::BoundaryCondition>(column, "BoundaryCondition", nb::is_arithmetic())
      .value("INLET_PRESSURE_INLET_VELOCITY", Column::BoundaryCondition::InletPressureInletVelocity)
      .value("INLET_PRESSURE_OUTLET_PRESSURE", Column::BoundaryCondition::InletPressureOutletPressure)
      .value("INLET_VELOCITY_OUTLET_PRESSURE", Column::BoundaryCondition::InletVelocityOutletPressure);

  column.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(nb::init<MixturePrediction, std::vector<Component>, Column::VelocityProfile, Column::BoundaryCondition, bool,
                    std::size_t, std::size_t, std::size_t, double, double, double, double, double, double, double,
                    double, double, double, double, double, double, double, double, double, double, double, double,
                    double, double, double>(),
           nb::arg("mixture"), nb::arg("components"), nb::arg("velocity_profile"), nb::arg("boundary_condition"),
           nb::arg("energy_balance"), nb::arg("number_of_grid_points"), nb::arg("max_isotherm_terms"),
           nb::arg("carrier_gas_component"), nb::arg("temperature"), nb::arg("inlet_pressure"),
           nb::arg("outlet_pressure"), nb::arg("pressure_gradient"), nb::arg("column_void_fraction"),
           nb::arg("particle_density"), nb::arg("column_entrance_velocity"), nb::arg("column_length"),
           nb::arg("dynamic_viscosity"), nb::arg("particle_diameter"), nb::arg("influx_temperature"),
           nb::arg("internal_diameter"), nb::arg("outer_diameter"), nb::arg("wall_density"),
           nb::arg("gas_thermal_conductivity"), nb::arg("wall_thermal_conductivity"),
           nb::arg("heat_transfer_gas_solid"), nb::arg("heat_transfer_gas_wall"),
           nb::arg("heat_transfer_wall_external"), nb::arg("heat_capacity_gas"), nb::arg("heat_capacity_solid"),
           nb::arg("heat_capacity_wall"))
      .def_rw("mixture", &Column::mixture)
      .def_rw("components", &Column::components)
      .def_rw("velocity_profile", &Column::velocityProfile)
      .def_rw("boundary_condition", &Column::boundaryCondition)
      .def_rw("energy_balance", &Column::energyBalance)
      .def_rw("number_of_grid_points", &Column::numberOfGridPoints)
      .def_rw("number_of_components", &Column::numberOfComponents)
      .def_rw("max_isotherm_terms", &Column::maxIsothermTerms)
      .def_rw("num_calls", &Column::numberOfCalls)
      .def_rw("carrier_gas_component", &Column::carrierGasComponent)
      .def_rw("external_temperature", &Column::externalTemperature)
      .def_rw("inlet_pressure", &Column::inletPressure)
      .def_rw("outlet_pressure", &Column::outletPressure)
      .def_rw("pressure_gradient", &Column::pressureGradient)
      .def_rw("void_fraction", &Column::voidFraction)
      .def_rw("particle_density", &Column::particleDensity)
      .def_rw("column_entrance_velocity", &Column::columnEntranceVelocity)
      .def_rw("column_length", &Column::columnLength)
      .def_rw("dynamic_viscosity", &Column::dynamicViscosity)
      .def_rw("particle_diameter", &Column::particleDiameter)
      .def_rw("influx_temperature", &Column::influxTemperature)
      .def_rw("internal_diameter", &Column::internalDiameter)
      .def_rw("outer_diameter", &Column::outerDiameter)
      .def_rw("wall_density", &Column::wallDensity)
      .def_rw("gas_thermal_conductivity", &Column::gasThermalConductivity)
      .def_rw("wall_thermal_conductivity", &Column::wallThermalConductivity)
      .def_rw("heat_transfer_gas_solid", &Column::heatTransferGasSolid)
      .def_rw("heat_transfer_gas_wall", &Column::heatTransferGasWall)
      .def_rw("heat_transfer_wall_external", &Column::heatTransferWallExternal)
      .def_rw("heat_capacity_gas", &Column::heatCapacityGas)
      .def_rw("heat_capacity_solid", &Column::heatCapacitySolid)
      .def_rw("heat_capacity_wall", &Column::heatCapacityWall)
      .def_rw("resolution", &Column::resolution)
      .def_rw("time_normalization_factor", &Column::timeNormalizationFactor)
      .def_rw("iast_performance", &Column::iastPerformance)
      .def_rw("prefactor_mass_transfer", &Column::prefactorMassTransfer)
      .def_rw("ideal_gas_mol_fractions", &Column::idealGasMolFractions)
      .def_rw("adsorbed_mol_fractions", &Column::adsorbedMolFractions)
      .def_rw("number_of_molecules", &Column::numberOfMolecules)
      .def_rw("interstitial_gas_velocity", &Column::interstitialGasVelocity)
      .def_rw("gas_density", &Column::gasDensity)
      .def_rw("total_concentration", &Column::totalConcentration)
      .def_rw("total_pressure", &Column::totalPressure)
      .def_rw("partial_pressure", &Column::partialPressure)
      .def_rw("equilibrium_adsorption", &Column::equilibriumAdsorption)
      .def_prop_ro("mole_fraction",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.moleFraction)); })
      .def_prop_ro("mole_fraction_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.moleFractionDot)); })
      .def_rw("cached_pressure", &Column::cachedPressure)
      .def_rw("cached_grand_potential", &Column::cachedGrandPotential)
      .def_rw("coeff_gas_gas", &Column::coeffGasGas)
      .def_rw("coeff_gas_solid", &Column::coeffGasSolid)
      .def_rw("coeff_gas_wall", &Column::coeffGasWall)
      .def_rw("coeff_diffusion", &Column::coeffDiffusion)
      .def_rw("face_pressures", &Column::facePressures)
      .def_rw("mass_flux", &Column::massFlux)
      .def_rw("state", &Column::state)
      .def_rw("state_dot", &Column::stateDot)
      .def_prop_ro("concentration",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.concentration)); })
      .def_prop_ro("adsorption",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.adsorption)); })
      .def_prop_ro("adsorption_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.adsorptionDot)); })
      .def_prop_ro("gas_temperature",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.gasTemperature)); })
      .def_prop_ro("gas_temperature_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.gasTemperatureDot)); })
      .def_prop_ro("solid_temperature",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.solidTemperature)); })
      .def_prop_ro("solid_temperature_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.solidTemperatureDot)); })
      .def_prop_ro("wall_temperature",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.wallTemperature)); })
      .def_prop_ro("wall_temperature_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.wallTemperatureDot)); })
      .def("state_size", &Column::stateSize)
      .def("bind_state_views", &Column::bindStateViews)
      .def("initialize", &Column::initialize)
      .def("set_temperature", &Column::setTemperature, nb::arg("temperature"))
      .def("__repr__", &Column::repr)
      .def("write_json", &Column::writeJSON, nb::arg("file_name"))
      .def("read_json", &Column::readJSON, nb::arg("file_name"));

  m.def("compute_pressure", static_cast<void (*)(Column&)>(&computePressure), nb::arg("column"));

  m.def("compute_equilibrium_loadings", static_cast<void (*)(Column&)>(&computeEquilibriumLoadings), nb::arg("column"));

  m.def("compute_velocity", static_cast<void (*)(Column&)>(&computeVelocity), nb::arg("column"));

  m.def("compute_first_derivatives", static_cast<void (*)(Column&)>(&computeDerivatives), nb::arg("column"));

  m.def("enforce_boundary_condition", &enforceBoundaryCondition, nb::arg("column"));

  m.def(
      "compute_weno",
      [](const std::vector<double>& input)
      {
        std::vector<double> output(input.size(), 0.0);
        computeWENO(std::span<const double>(input.data(), input.size()),
                    std::span<double>(output.data(), output.size()));
        return output;
      },
      nb::arg("input"));

  nb::class_<RungeKutta3>(m, "RungeKutta3")
      .def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(nb::init<double, bool, std::size_t>(), nb::arg("time_step"), nb::arg("auto_steps"),
           nb::arg("number_of_steps"))
      .def_rw("time_step", &RungeKutta3::timeStep)
      .def_rw("auto_steps", &RungeKutta3::autoNumberOfSteps)
      .def_rw("number_of_steps", &RungeKutta3::numberOfSteps)
      .def(
          "propagate",
          [](RungeKutta3& self, Column& column, std::size_t step)
          {
            Timing timings;
            return self.propagate(column, step, timings);
          },
          nb::arg("column"), nb::arg("step"));

  nb::class_<SemiImplicitRungeKutta3>(m, "SemiImplicitRungeKutta3")
      .def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(nb::init<double, bool, std::size_t>(), nb::arg("time_step"), nb::arg("auto_steps"),
           nb::arg("number_of_steps"))
      .def_rw("time_step", &SemiImplicitRungeKutta3::timeStep)
      .def_rw("auto_steps", &SemiImplicitRungeKutta3::autoNumberOfSteps)
      .def_rw("number_of_steps", &SemiImplicitRungeKutta3::numberOfSteps)
      .def(
          "propagate",
          [](SemiImplicitRungeKutta3& self, Column& column, std::size_t step)
          {
            Timing timings;
            return self.propagate(column, step, timings);
          },
          nb::arg("column"), nb::arg("step"));

  nb::class_<CVODE>(m, "CVODE")
      .def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(nb::init<double, bool, std::size_t>(), nb::arg("time_step"), nb::arg("auto_steps"),
           nb::arg("number_of_steps"))
      .def_rw("time_step", &CVODE::timeStep)
      .def_rw("auto_steps", &CVODE::autoNumberOfSteps)
      .def_rw("number_of_steps", &CVODE::numberOfSteps)
      .def("initialize", &CVODE::initialize, nb::arg("column"))
      .def(
          "propagate",
          [](CVODE& self, Column& column, std::size_t step)
          {
            Timing timings;
            return self.propagate(column, step, timings);
          },
          nb::arg("column"), nb::arg("step"));

  nb::class_<Breakthrough> breakthrough(m, "Breakthrough");

  nb::enum_<Breakthrough::IntegrationScheme>(breakthrough, "IntegrationScheme", nb::is_arithmetic())
      .value("SSP_RK", Breakthrough::IntegrationScheme::SSP_RK)
      .value("CVODE", Breakthrough::IntegrationScheme::CVODE)
      .value("ITERATIVE", Breakthrough::IntegrationScheme::Iterative)
      .value("SIRK3", Breakthrough::IntegrationScheme::SIRK3);

  breakthrough.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(
          "__init__",
          [](Breakthrough* self, std::string displayName, std::size_t carrierGasComponent,
             std::size_t numberOfComponents, std::size_t numberOfGridPoints, std::size_t printEvery,
             std::size_t writeEvery, double timeStep, std::size_t numberOfInitTimeSteps, std::size_t numberOfTimeSteps,
             bool autoNumberOfTimeSteps, std::size_t maxIsothermTerms, Column column,
             Breakthrough::IntegrationScheme integrationScheme, std::optional<std::string> readColumnFile)
          {
            new (self) Breakthrough(std::move(displayName), carrierGasComponent, numberOfComponents, numberOfGridPoints,
                                    printEvery, writeEvery, timeStep, numberOfInitTimeSteps, numberOfTimeSteps,
                                    autoNumberOfTimeSteps, maxIsothermTerms, std::move(column),
                                    RungeKutta3(timeStep, autoNumberOfTimeSteps, numberOfTimeSteps),
                                    SemiImplicitRungeKutta3(timeStep, autoNumberOfTimeSteps, numberOfTimeSteps),
                                    CVODE(timeStep, autoNumberOfTimeSteps, numberOfTimeSteps), integrationScheme,
                                    std::move(readColumnFile));
          },
          nb::arg("display_name"), nb::arg("carrier_gas_component"), nb::arg("number_of_components"),
          nb::arg("number_of_grid_points"), nb::arg("print_every"), nb::arg("write_every"), nb::arg("time_step"),
          nb::arg("number_of_init_time_steps"), nb::arg("number_of_time_steps"), nb::arg("auto_number_of_time_steps"),
          nb::arg("max_isotherm_terms"), nb::arg("column"), nb::arg("integration_scheme"),
          nb::arg("read_column_file") = nb::none())
      .def_ro("display_name", &Breakthrough::displayName)
      .def_rw("carrier_gas_component", &Breakthrough::carrierGasComponent)
      .def_rw("number_of_components", &Breakthrough::numberOfComponents)
      .def_rw("number_of_grid_points", &Breakthrough::numberOfGridPoints)
      .def_rw("print_every", &Breakthrough::printEvery)
      .def_rw("write_every", &Breakthrough::writeEvery)
      .def_rw("time_step", &Breakthrough::timeStep)
      .def_rw("number_of_init_time_steps", &Breakthrough::numberOfInitTimeSteps)
      .def_rw("number_of_time_steps", &Breakthrough::numberOfSteps)
      .def_rw("auto_number_of_time_steps", &Breakthrough::autoNumberOfSteps)
      .def_rw("max_isotherm_terms", &Breakthrough::maxIsothermTerms)
      .def_rw("column", &Breakthrough::column)
      .def_rw("rk3", &Breakthrough::rk3)
      .def_rw("sirk3", &Breakthrough::sirk3)
      .def_prop_ro(
          "cvode", [](Breakthrough& self) -> CVODE& { return self.cvode; }, nb::rv_policy::reference_internal)
      .def_rw("integration_scheme", &Breakthrough::integrationScheme)
      .def("print", &Breakthrough::print)
      .def("__repr__", &Breakthrough::repr)
      .def("run", &Breakthrough::run)
      .def("compute",
           [](Breakthrough& self) -> nb::ndarray<nb::numpy, double, nb::ndim<3>>
           {
             const size_t columnSize = 5 * self.numberOfComponents + 5;
             auto* buffer = new std::vector<double>();
             buffer->reserve(((self.numberOfSteps / std::max<size_t>(self.writeEvery, 1)) + 1) *
                             (self.numberOfGridPoints + 1) * columnSize);

             // Loop can quit early if autoNumberOfSteps.
             for (size_t step = 0; (step < self.numberOfSteps || self.autoNumberOfSteps); ++step)
             {
               // Check for an interrupt/error from the Python side.
               if (PyErr_CheckSignals() != 0)
               {
                 delete buffer;
                 throw nb::python_error();
               }

               self.computeStep(step);
               double time = static_cast<double>(step) * self.timeStep;

               if (step % self.writeEvery == 0)
               {
                 for (size_t grid = 0; grid < self.numberOfGridPoints + 1; ++grid)
                 {
                   buffer->push_back(time * self.column.columnEntranceVelocity / self.column.columnLength);
                   buffer->push_back(time / 60.0);
                   buffer->push_back(static_cast<double>(grid) * self.column.resolution);
                   buffer->push_back(self.column.interstitialGasVelocity[grid]);
                   buffer->push_back(self.column.totalPressure[grid]);

                   for (size_t comp = 0; comp < self.numberOfComponents; ++comp)
                   {
                     const size_t index = grid * self.numberOfComponents + comp;
                     const double normalizedPartialPressureDenominator =
                         self.column.totalPressure[grid] * self.column.components[comp].initialGasMoleFraction;

                     buffer->push_back(self.column.adsorption[index]);
                     buffer->push_back(self.column.equilibriumAdsorption[index]);
                     buffer->push_back(self.column.partialPressure[index]);
                     buffer->push_back(normalizedPartialPressureDenominator != 0.0
                                           ? self.column.partialPressure[index] / normalizedPartialPressureDenominator
                                           : 0.0);
                     buffer->push_back(self.column.adsorptionDot[index]);
                   }
                 }
               }

               if (step % self.printEvery == 0)
               {
                 const double averageNumberOfMixturePredictionSteps =
                     self.column.iastPerformance.second > 0
                         ? static_cast<double>(self.column.iastPerformance.first) /
                               static_cast<double>(self.column.iastPerformance.second)
                         : 0.0;
                 std::cout << "Timestep " + std::to_string(step) + ", time: " + std::to_string(time) + " [s]"
                           << std::endl;
                 std::cout << "    Average number of mixture-prediction steps: " +
                                  std::to_string(averageNumberOfMixturePredictionSteps)
                           << std::endl;
               }
             }

             std::cout << "Final timestep " + std::to_string(self.numberOfSteps) +
                              ", time: " + std::to_string(self.timeStep * static_cast<double>(self.numberOfSteps)) +
                              " [s]"
                       << std::endl;

             const size_t numberOfRows = self.numberOfGridPoints + 1;
             const size_t numberOfSnapshots = buffer->size() / (numberOfRows * columnSize);
             nb::capsule owner(buffer,
                               [](void* pointer) noexcept { delete static_cast<std::vector<double>*>(pointer); });

             return nb::ndarray<nb::numpy, double, nb::ndim<3>>(buffer->data(),
                                                                {numberOfSnapshots, numberOfRows, columnSize}, owner);
           })
      .def("compute_step", &Breakthrough::computeStep, nb::arg("step"))
      .def(
          "set_components_parameters",
          [](Breakthrough& self, std::vector<double> molfracs, std::vector<double> params)
          {
            size_t index = 0;
            for (size_t i = 0; i < self.numberOfComponents; ++i)
            {
              self.column.components[i].initialGasMoleFraction = molfracs[i];
              size_t numberOfParameters = self.column.components[i].isotherm.numberOfParameters;
              std::vector<double> slicedVec(params.begin() + index, params.begin() + index + numberOfParameters);
              index += numberOfParameters;
              self.column.components[i].isotherm.setParameters(slicedVec);
            }

            // Also set for mixture prediction.
            MixturePrediction& mixture = self.column.mixture;
            index = 0;
            for (size_t i = 0; i < mixture.numberOfComponents; ++i)
            {
              mixture.components[i].initialGasMoleFraction = molfracs[i];
              size_t numberOfParameters = mixture.components[i].isotherm.parameters.size();
              std::vector<double> slicedVec(params.begin() + index, params.begin() + index + numberOfParameters);
              index += numberOfParameters;
              mixture.components[i].isotherm.setParameters(slicedVec);
            }
            mixture.sortedComponents = mixture.components;
            mixture.segregatedSortedComponents = std::vector<std::vector<Component>>(
                mixture.maxIsothermTerms, std::vector<Component>(mixture.components));
            mixture.sortComponents();
          },
          nb::arg("molfracs"), nb::arg("params"))
      .def("get_components_parameters",
           [](const Breakthrough& self)
           {
             std::vector<double> params;
             for (size_t i = 0; i < self.numberOfComponents; ++i)
             {
               std::vector<double> compParams = self.column.components[i].isotherm.getParameters();
               params.insert(params.end(), compParams.begin(), compParams.end());
             }
             return params;
           });

  nb::class_<Fitting> fitting(m, "Fitting");

  nb::enum_<Fitting::PressureScale>(fitting, "PressureScale", nb::is_arithmetic())
      .value("LOG", Fitting::PressureScale::Log)
      .value("NORMAL", Fitting::PressureScale::Normal);

  fitting.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def_rw("number_of_components", &Fitting::numberOfComponents)
      .def_rw("components", &Fitting::components)
      .def_rw("display_name", &Fitting::displayName)
      .def_rw("component_name", &Fitting::componentName)
      .def_rw("filename", &Fitting::filename)
      .def_rw("isotherms", &Fitting::isotherms)
      .def_rw("column_pressure", &Fitting::columnPressure)
      .def_rw("column_loading", &Fitting::columnLoading)
      .def_rw("column_error", &Fitting::columnError)
      .def_rw("maximum_loading", &Fitting::maximumLoading)
      .def_rw("external_temperature", &Fitting::externalTemperature)
      .def_rw("pressure_scale", &Fitting::pressureScale)
      .def_rw("raw_data", &Fitting::rawData)
      .def("read_data", &Fitting::readData, nb::arg("component"))
      .def("print_solution", &Fitting::printSolution, nb::arg("component"))
      .def("run", &Fitting::run)
      .def("write_components_json", &Fitting::writeComponentsJson, nb::arg("path"));

  nb::class_<SwingAdsorption::SubStage>(m, "SubStage")
      .def_rw("temperature", &SwingAdsorption::SubStage::temperature)
      .def_rw("pressure", &SwingAdsorption::SubStage::pressure)
      .def_rw("number_of_steps", &SwingAdsorption::SubStage::numberOfSteps);

  nb::class_<SwingAdsorption>(m, "SwingAdsorption")
      .def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def_rw("breakthrough", &SwingAdsorption::breakthrough)
      .def_rw("sub_stages", &SwingAdsorption::subStages)
      .def("run", &SwingAdsorption::run)
      .def("print", &SwingAdsorption::print)
      .def("__repr__", &SwingAdsorption::repr);

  m.def("read_input", [](const std::string& fileName) { return InputReader(fileName); }, nb::arg("file_name"));

  m.def(
      "simulation_from_file",
      [](const std::string& fileName) -> nb::object
      {
        InputReader inputReader(fileName);
        switch (inputReader.simulationType)
        {
          case InputReader::SimulationType::Breakthrough:
            return nb::cast(Breakthrough(inputReader));
          case InputReader::SimulationType::MixturePrediction:
            return nb::cast(MixturePrediction(inputReader));
          case InputReader::SimulationType::Fitting:
            return nb::cast(Fitting(inputReader));
          case InputReader::SimulationType::SwingAdsorption:
            return nb::cast(SwingAdsorption(inputReader));
          case InputReader::SimulationType::Test:
            throw std::invalid_argument("SimulationType 'Test' does not have a Python simulation object");
        }

        throw std::invalid_argument("Unknown SimulationType");
      },
      nb::arg("file_name"));
}

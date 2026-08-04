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
#include <type_traits>
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
#include "reaction.h"
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

template <typename ColumnType>
nb::ndarray<nb::numpy, double, nb::ndim<3>> computeBreakthrough(Breakthrough<ColumnType>& self)
{
  ColumnType& column = self.column;
  const size_t columnSize = 5 * self.numberOfComponents + 5;
  auto* buffer = new std::vector<double>();
  buffer->reserve(((self.numberOfSteps / std::max<size_t>(self.writeEvery, 1)) + 1) * (self.numberOfGridPoints + 1) *
                  columnSize);

  for (size_t step = 0; (step < self.numberOfSteps || self.autoNumberOfSteps); ++step)
  {
    if (PyErr_CheckSignals() != 0)
    {
      delete buffer;
      throw nb::python_error();
    }

    self.computeStep(step);
    const double time = static_cast<double>(step) * self.timeStep;

    if (step % self.writeEvery == 0)
    {
      for (size_t grid = 0; grid < self.numberOfGridPoints + 1; ++grid)
      {
        buffer->push_back(time * column.columnEntranceVelocity / column.columnLength);
        buffer->push_back(time / 60.0);
        if constexpr (std::is_same_v<ColumnType, Column>)
        {
          buffer->push_back(static_cast<double>(grid) * column.resolution);
        }
        else
        {
          buffer->push_back(column.columnDistances[grid]);
        }
        buffer->push_back(column.interstitialGasVelocity[grid]);
        buffer->push_back(column.totalPressure[grid]);

        for (size_t comp = 0; comp < self.numberOfComponents; ++comp)
        {
          const size_t index = grid * self.numberOfComponents + comp;
          const double normalizedPartialPressureDenominator =
              column.totalPressure[grid] * column.components[comp].initialGasMoleFraction;

          buffer->push_back(column.physisorption[index]);
          buffer->push_back(column.equilibriumPhysisorption[index]);
          buffer->push_back(column.partialPressure[index]);
          buffer->push_back(normalizedPartialPressureDenominator != 0.0
                                ? column.partialPressure[index] / normalizedPartialPressureDenominator
                                : 0.0);
          buffer->push_back(column.physisorptionDot[index]);
        }
      }
    }

    if (step % self.printEvery == 0)
    {
      const double averageNumberOfMixturePredictionSteps =
          column.iastPerformance.second > 0
              ? static_cast<double>(column.iastPerformance.first) / static_cast<double>(column.iastPerformance.second)
              : 0.0;
      std::cout << "Timestep " + std::to_string(step) + ", time: " + std::to_string(time) + " [s]" << std::endl;
      std::cout << "    Average number of mixture-prediction steps: " +
                       std::to_string(averageNumberOfMixturePredictionSteps)
                << std::endl;
    }
  }

  const size_t numberOfRows = self.numberOfGridPoints + 1;
  const size_t numberOfSnapshots = buffer->size() / (numberOfRows * columnSize);
  nb::capsule owner(buffer, [](void* pointer) noexcept { delete static_cast<std::vector<double>*>(pointer); });

  return nb::ndarray<nb::numpy, double, nb::ndim<3>>(buffer->data(), {numberOfSnapshots, numberOfRows, columnSize},
                                                     owner);
}

template <typename ColumnType>
void bindBreakthrough(nb::module_& module, const char* name)
{
  using Simulation = Breakthrough<ColumnType>;
  nb::class_<Simulation> simulation(module, name);
  simulation.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def_ro("display_name", &Simulation::displayName)
      .def_rw("carrier_gas_component", &Simulation::carrierGasComponent)
      .def_rw("number_of_components", &Simulation::numberOfComponents)
      .def_rw("number_of_grid_points", &Simulation::numberOfGridPoints)
      .def_rw("print_every", &Simulation::printEvery)
      .def_rw("write_every", &Simulation::writeEvery)
      .def_rw("time_step", &Simulation::timeStep)
      .def_rw("number_of_init_time_steps", &Simulation::numberOfInitTimeSteps)
      .def_rw("number_of_time_steps", &Simulation::numberOfSteps)
      .def_rw("auto_number_of_time_steps", &Simulation::autoNumberOfSteps)
      .def_rw("max_isotherm_terms", &Simulation::maxIsothermTerms)
      .def_rw("column", &Simulation::column)
      .def_rw("rk3", &Simulation::rk3)
      .def_rw("sirk3", &Simulation::sirk3)
      .def_prop_ro(
          "cvode", [](Simulation& self) -> CVODE& { return self.cvode; }, nb::rv_policy::reference_internal)
      .def_rw("integration_scheme", &Simulation::integrationScheme)
      .def("print", &Simulation::print)
      .def("__repr__", &Simulation::repr)
      .def("run", &Simulation::run)
      .def("compute", &computeBreakthrough<ColumnType>)
      .def("compute_step", &Simulation::computeStep, nb::arg("step"));

  if constexpr (std::is_same_v<ColumnType, Column>)
  {
    simulation
        .def(
            "set_components_parameters",
            [](Simulation& self, std::vector<double> molfracs, std::vector<double> params)
            {
              size_t index = 0;
              for (size_t i = 0; i < self.numberOfComponents; ++i)
              {
                self.column.components[i].initialGasMoleFraction = molfracs[i];
                const size_t numberOfParameters = self.column.components[i].isotherm.numberOfParameters;
                std::vector<double> slicedVec(params.begin() + static_cast<std::ptrdiff_t>(index),
                                              params.begin() + static_cast<std::ptrdiff_t>(index + numberOfParameters));
                index += numberOfParameters;
                self.column.components[i].isotherm.setParameters(slicedVec);
              }

              MixturePrediction& mixture = self.column.physisorptionMixture;
              index = 0;
              for (size_t i = 0; i < mixture.numberOfComponents; ++i)
              {
                mixture.components[i].initialGasMoleFraction = molfracs[i];
                const size_t numberOfParameters = mixture.components[i].isotherm.parameters.size();
                std::vector<double> slicedVec(params.begin() + static_cast<std::ptrdiff_t>(index),
                                              params.begin() + static_cast<std::ptrdiff_t>(index + numberOfParameters));
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
             [](const Simulation& self)
             {
               std::vector<double> params;
               for (size_t i = 0; i < self.numberOfComponents; ++i)
               {
                 std::vector<double> componentParameters = self.column.components[i].isotherm.getParameters();
                 params.insert(params.end(), componentParameters.begin(), componentParameters.end());
               }
               return params;
             });
  }
}

template <typename ColumnType>
void bindSwingAdsorption(nb::module_& module, const char* name)
{
  using Simulation = SwingAdsorption<ColumnType>;
  nb::class_<Simulation>(module, name)
      .def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def_rw("breakthrough", &Simulation::breakthrough)
      .def_rw("sub_stages", &Simulation::subStages)
      .def("run", &Simulation::run)
      .def("print", &Simulation::print)
      .def("__repr__", &Simulation::repr);
}
}  // namespace

EXPAND_MODULE(MODULE_NAME)
{
  nb::class_<Reaction> reaction(m, "Reaction");
  nb::enum_<Reaction::Phase>(reaction, "Phase", nb::is_arithmetic())
      .value("PHYSISORBED", Reaction::Phase::Physisorbed)
      .value("CHEMISORBED", Reaction::Phase::Chemisorbed)
      .value("PORE_CONCENTRATION", Reaction::Phase::PoreConcentration);
  nb::enum_<Reaction::Style>(reaction, "Style", nb::is_arithmetic())
      .value("GENERAL_POWER_LAW", Reaction::Style::GeneralPowerLaw)
      .value("LANGMUIR_HINSHELWOOD", Reaction::Style::LangmuirHinshelwood)
      .value("LANGMUIR_HINSHELWOOD_HOUGEN_WATSON", Reaction::Style::LangmuirHinshelwoodHougenWatson);
  reaction.def(nb::init<>())
      .def_rw("phase", &Reaction::phase)
      .def_rw("style", &Reaction::style)
      .def_rw("site", &Reaction::site)
      .def_rw("reactants", &Reaction::reactants)
      .def_rw("products", &Reaction::products)
      .def_rw("reactant_stoichiometry", &Reaction::reactantStoichiometry)
      .def_rw("product_stoichiometry", &Reaction::productStoichiometry)
      .def_rw("forward_orders", &Reaction::forwardOrders)
      .def_rw("backward_orders", &Reaction::backwardOrders)
      .def_rw("forward_rate_coefficient", &Reaction::forwardRateCoefficient)
      .def_rw("forward_activation_energy", &Reaction::forwardActivationEnergy)
      .def_rw("equilibrium_constant", &Reaction::equilibriumConstant)
      .def_rw("gibbs_free_energy", &Reaction::gibbsFreeEnergy)
      .def_rw("rate_limit_time", &Reaction::rateLimitTime)
      .def("equilibrium_constant_at", &Reaction::equilibriumConstantAt, nb::arg("temperature"))
      .def("forward_rate_constant", &Reaction::forwardRateConstant, nb::arg("temperature"))
      .def("__repr__", &Reaction::repr);

  nb::class_<InputReader> inputReader(m, "InputReader");

  nb::enum_<InputReader::SimulationType>(inputReader, "SimulationType", nb::is_arithmetic())
      .value("BREAKTHROUGH", InputReader::SimulationType::Breakthrough)
      .value("MIXTURE_PREDICTION", InputReader::SimulationType::MixturePrediction)
      .value("FITTING", InputReader::SimulationType::Fitting)
      .value("SWING_ADSORPTION", InputReader::SimulationType::SwingAdsorption)
      .value("TEST", InputReader::SimulationType::Test);

  nb::class_<InputReader::SwingAdsorptionPhase>(m, "InputSwingAdsorptionPhase")
      .def(nb::init<>())
      .def_rw("name", &InputReader::SwingAdsorptionPhase::name)
      .def_rw("temperature", &InputReader::SwingAdsorptionPhase::temperature)
      .def_rw("inlet_pressure", &InputReader::SwingAdsorptionPhase::inletPressure)
      .def_rw("number_of_steps", &InputReader::SwingAdsorptionPhase::numberOfSteps);

  inputReader.def(nb::init<const std::string>(), nb::arg("file_name"))
      .def_rw("components", &InputReader::components)
      .def_rw("reactions", &InputReader::reactions)
      .def_rw("adsorbent_components", &InputReader::adsorbentComponents)
      .def_rw("number_of_carrier_gases", &InputReader::numberOfCarrierGases)
      .def_rw("carrier_gas_component", &InputReader::carrierGasComponent)
      .def_rw("max_isotherm_terms", &InputReader::maxIsothermTerms)
      .def_rw("simulation_type", &InputReader::simulationType)
      .def_rw("mixture_prediction_method", &InputReader::mixturePredictionMethod)
      .def_rw("iast_method", &InputReader::IASTMethod)
      .def_rw("breakthrough_integrator", &InputReader::breakthroughIntegrator)
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
      .def_rw("column_distances", &InputReader::columnDistances)
      .def_rw("adsorbent_lengths", &InputReader::adsorbentLengths)
      .def_rw("adsorbent_interface_lengths", &InputReader::adsorbentInterfaceLengths)
      .def_rw("adsorbent_grid_points", &InputReader::adsorbentGridPoints)
      .def_rw("adsorbent_void_fractions", &InputReader::adsorbentVoidFractions)
      .def_rw("adsorbent_particle_densities", &InputReader::adsorbentParticleDensities)
      .def_rw("adsorbent_particle_diameters", &InputReader::adsorbentParticleDiameters)
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
      .def_rw("swing_adsorption_phases", &InputReader::swingAdsorptionPhases)
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
      .def_prop_ro("number_of_sites", [](const MultiSiteIsotherm& self) { return self.sites.size(); })
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

  nb::class_<Chemisorption> chemisorption(m, "Chemisorption");
  nb::enum_<Chemisorption::Type>(chemisorption, "Type", nb::is_arithmetic())
      .value("None_", Chemisorption::Type::None)
      .value("FirstOrder", Chemisorption::Type::FirstOrder)
      .value("PseudoNth", Chemisorption::Type::PseudoNth)
      .value("Avrami", Chemisorption::Type::Avrami)
      .value("General", Chemisorption::Type::General)
      .value("Elovich", Chemisorption::Type::Elovich);
  chemisorption.def(nb::init<>())
      .def_rw("type", &Chemisorption::type)
      .def_rw("rate_coefficient", &Chemisorption::rateCoefficient)
      .def_rw("order", &Chemisorption::order)
      .def_rw("maximum_loading", &Chemisorption::maximumLoading)
      .def_rw("heat_of_chemisorption", &Chemisorption::heatOfChemisorption)
      .def_rw("adsorption_rate_coefficient", &Chemisorption::adsorptionRateCoefficient)
      .def_rw("adsorption_activation_energy", &Chemisorption::adsorptionActivationEnergy)
      .def_rw("desorption_rate_coefficient", &Chemisorption::desorptionRateCoefficient)
      .def_rw("desorption_activation_energy", &Chemisorption::desorptionActivationEnergy)
      .def_rw("pore_concentration_order", &Chemisorption::poreConcentrationOrder)
      .def_rw("capacity_order", &Chemisorption::capacityOrder)
      .def_rw("desorption_order", &Chemisorption::desorptionOrder)
      .def_rw("elovich_alpha", &Chemisorption::elovichAlpha)
      .def_rw("elovich_beta", &Chemisorption::elovichBeta)
      .def_rw("film_mass_transfer_coefficient", &Chemisorption::filmMassTransferCoefficient)
      .def_rw("pore_diffusivity", &Chemisorption::poreDiffusivity)
      .def_rw("use_pore_surface_transport", &Chemisorption::usePoreSurfaceTransport)
      .def_rw("isotherm", &Chemisorption::isotherm)
      .def("rate", &Chemisorption::rate, nb::arg("equilibrium_loading"), nb::arg("loading"),
           nb::arg("concentration") = 0.0, nb::arg("temperature") = 298.15)
      .def("uses_surface_pore_transport", &Chemisorption::usesSurfacePoreTransport)
      .def("__repr__", &Chemisorption::repr);

  nb::class_<MultiSiteChemisorption>(m, "MultiSiteChemisorption")
      .def(nb::init<>())
      .def(nb::init<std::vector<Chemisorption>>(), nb::arg("sites"))
      .def_rw("sites", &MultiSiteChemisorption::sites)
      .def_rw("number_of_sites", &MultiSiteChemisorption::numberOfSites)
      .def("add", &MultiSiteChemisorption::add, nb::arg("site"))
      .def("enabled", &MultiSiteChemisorption::enabled)
      .def("uses_surface_pore_transport", &MultiSiteChemisorption::usesSurfacePoreTransport)
      .def("maximum_loading", &MultiSiteChemisorption::maximumLoading)
      .def("__repr__", &MultiSiteChemisorption::repr);

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
      .def_rw("chemisorption", &Component::chemisorption)
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
      .value("SEI", MixturePrediction::PredictionMethod::SEI)
      .value("SCI", MixturePrediction::PredictionMethod::SCI)
      .value("SPI", MixturePrediction::PredictionMethod::SPI);

  nb::enum_<MixturePrediction::IASTMethod>(mixturePrediction, "IASTMethod", nb::is_arithmetic())
      .value("FAST_IAST", MixturePrediction::IASTMethod::FastIAST)
      .value("NESTED_LOOP_BISECTION", MixturePrediction::IASTMethod::NestedLoopBisection);

  nb::enum_<MixturePrediction::PressureScale>(mixturePrediction, "PressureScale", nb::is_arithmetic())
      .value("LOG", MixturePrediction::PressureScale::Log)
      .value("LINEAR", MixturePrediction::PressureScale::Linear);

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
      .def_rw("segregated_number_of_sorted_components", &MixturePrediction::segregatedNumberOfSortedComponents)
      .def_rw("equilibrium_site_loadings", &MixturePrediction::equilibriumSiteLoadings)
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
                                   cachedPressure, cachedGrandPotential, gasTemperature);

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
                std::span<double>(numberOfMolecules.data(), numberOfMolecules.size()), cachedPressure,
                cachedGrandPotential, gasTemperatureValue);

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

  nb::enum_<Column::BoundaryCondition>(column, "BoundaryCondition", nb::is_arithmetic())
      .value("INLET_PRESSURE_INLET_VELOCITY", Column::BoundaryCondition::InletPressureInletVelocity)
      .value("INLET_PRESSURE_OUTLET_PRESSURE", Column::BoundaryCondition::InletPressureOutletPressure)
      .value("INLET_VELOCITY_OUTLET_PRESSURE", Column::BoundaryCondition::InletVelocityOutletPressure)
      .value("FIXED_VELOCITY", Column::BoundaryCondition::FixedVelocity)
      .value("FIXED_PRESSURE_INLET_VELOCITY", Column::BoundaryCondition::FixedPressureInletVelocity);

  column.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def(
          nb::init<MixturePrediction, std::vector<Component>, Column::BoundaryCondition, bool, std::size_t, std::size_t,
                   std::size_t, double, double, double, double, double, double, double, double, double, double, double,
                   double, double, double, double, double, double, double, double, double, double, double>(),
          nb::arg("mixture"), nb::arg("components"), nb::arg("boundary_condition"), nb::arg("energy_balance"),
          nb::arg("number_of_grid_points"), nb::arg("max_isotherm_terms"), nb::arg("carrier_gas_component"),
          nb::arg("temperature"), nb::arg("inlet_pressure"), nb::arg("outlet_pressure"), nb::arg("pressure_gradient"),
          nb::arg("column_void_fraction"), nb::arg("particle_density"), nb::arg("column_entrance_velocity"),
          nb::arg("column_length"), nb::arg("dynamic_viscosity"), nb::arg("particle_diameter"),
          nb::arg("influx_temperature"), nb::arg("internal_diameter"), nb::arg("outer_diameter"),
          nb::arg("wall_density"), nb::arg("gas_thermal_conductivity"), nb::arg("wall_thermal_conductivity"),
          nb::arg("heat_transfer_gas_solid"), nb::arg("heat_transfer_gas_wall"), nb::arg("heat_transfer_wall_external"),
          nb::arg("heat_capacity_gas"), nb::arg("heat_capacity_solid"), nb::arg("heat_capacity_wall"))
      .def_rw("physisorption_mixture", &Column::physisorptionMixture)
      .def_rw("chemisorption_mixture", &Column::chemisorptionMixture)
      .def_rw("components", &Column::components)
      .def_rw("reactions", &Column::reactions)
      .def_rw("boundary_condition", &Column::boundaryCondition)
      .def_rw("energy_balance", &Column::energyBalance)
      .def_rw("number_of_grid_points", &Column::numberOfGridPoints)
      .def_rw("number_of_components", &Column::numberOfComponents)
      .def_rw("max_isotherm_terms", &Column::maxIsothermTerms)
      .def_ro("max_chemisorption_sites", &Column::maxChemisorptionSites)
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
      .def_rw("equilibrium_physisorption", &Column::equilibriumPhysisorption)
      .def_rw("equilibrium_chemisorption", &Column::equilibriumChemisorption)
      .def_prop_ro("concentration",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.concentration)); })
      .def_prop_ro("concentration_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.concentrationDot)); })
      .def_prop_ro("mole_fraction",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.moleFraction)); })
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
      .def_prop_ro("physisorption",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.physisorption)); })
      .def_prop_ro("physisorption_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.physisorptionDot)); })
      .def_prop_ro("chemisorption",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.chemisorption)); })
      .def_prop_ro("chemisorption_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.chemisorptionDot)); })
      .def_prop_ro("surface_concentration",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.surfaceConcentration)); })
      .def_prop_ro("surface_concentration_dot", [](const Column& self)
                   { return spanToVector(std::span<const double>(self.surfaceConcentrationDot)); })
      .def_prop_ro("pore_concentration",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.poreConcentration)); })
      .def_prop_ro("pore_concentration_dot",
                   [](const Column& self) { return spanToVector(std::span<const double>(self.poreConcentrationDot)); })
      .def_rw("bulk_species_sink", &Column::bulkSpeciesSink)
      .def_rw("reaction_physisorption_source", &Column::reactionPhysisorptionSource)
      .def_rw("reaction_chemisorption_source", &Column::reactionChemisorptionSource)
      .def_rw("reaction_pore_concentration_source", &Column::reactionPoreConcentrationSource)
      .def_rw("reaction_heat", &Column::reactionHeat)
      .def_prop_ro("surface_pore_transport_enabled",
                   [](const Column& self) { return self.surfacePoreTransportEnabled; })
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

  nb::class_<MultibedColumn> multibedColumn(m, "MultibedColumn");
  nb::enum_<MultibedColumn::BoundaryCondition>(multibedColumn, "BoundaryCondition", nb::is_arithmetic())
      .value("INLET_PRESSURE_INLET_VELOCITY", MultibedColumn::BoundaryCondition::InletPressureInletVelocity)
      .value("INLET_PRESSURE_OUTLET_PRESSURE", MultibedColumn::BoundaryCondition::InletPressureOutletPressure)
      .value("INLET_VELOCITY_OUTLET_PRESSURE", MultibedColumn::BoundaryCondition::InletVelocityOutletPressure)
      .value("FIXED_VELOCITY", MultibedColumn::BoundaryCondition::FixedVelocity)
      .value("FIXED_PRESSURE_INLET_VELOCITY", MultibedColumn::BoundaryCondition::FixedPressureInletVelocity);

  multibedColumn.def(nb::init<const InputReader&>(), nb::arg("input_reader"))
      .def_rw("physisorption_mixtures", &MultibedColumn::physisorptionMixtures)
      .def_rw("chemisorption_mixtures", &MultibedColumn::chemisorptionMixtures)
      .def_rw("components", &MultibedColumn::components)
      .def_rw("reactions", &MultibedColumn::reactions)
      .def_rw("boundary_condition", &MultibedColumn::boundaryCondition)
      .def_rw("energy_balance", &MultibedColumn::energyBalance)
      .def_rw("number_of_grid_points", &MultibedColumn::numberOfGridPoints)
      .def_rw("number_of_components", &MultibedColumn::numberOfComponents)
      .def_rw("number_of_adsorbents", &MultibedColumn::numberOfAdsorbents)
      .def_rw("max_isotherm_terms", &MultibedColumn::maxIsothermTerms)
      .def_ro("max_chemisorption_sites", &MultibedColumn::maxChemisorptionSites)
      .def_rw("carrier_gas_component", &MultibedColumn::carrierGasComponent)
      .def_rw("adsorbent_lengths", &MultibedColumn::adsorbentLengths)
      .def_rw("adsorbent_interface_lengths", &MultibedColumn::adsorbentInterfaceLengths)
      .def_rw("adsorbent_void_fractions", &MultibedColumn::adsorbentVoidFractions)
      .def_rw("particle_densities", &MultibedColumn::particleDensities)
      .def_rw("particle_diameters", &MultibedColumn::particleDiameters)
      .def_rw("external_temperature", &MultibedColumn::externalTemperature)
      .def_rw("inlet_pressure", &MultibedColumn::inletPressure)
      .def_rw("outlet_pressure", &MultibedColumn::outletPressure)
      .def_rw("pressure_gradient", &MultibedColumn::pressureGradient)
      .def_rw("column_entrance_velocity", &MultibedColumn::columnEntranceVelocity)
      .def_rw("column_length", &MultibedColumn::columnLength)
      .def_rw("column_distances", &MultibedColumn::columnDistances)
      .def_rw("interstitial_gas_velocity", &MultibedColumn::interstitialGasVelocity)
      .def_rw("gas_density", &MultibedColumn::gasDensity)
      .def_rw("total_concentration", &MultibedColumn::totalConcentration)
      .def_rw("total_pressure", &MultibedColumn::totalPressure)
      .def_rw("mole_fraction", &MultibedColumn::moleFraction)
      .def_rw("partial_pressure", &MultibedColumn::partialPressure)
      .def_rw("equilibrium_physisorption", &MultibedColumn::equilibriumPhysisorption)
      .def_rw("equilibrium_chemisorption", &MultibedColumn::equilibriumChemisorption)
      .def_rw("fraction_of_adsorbent", &MultibedColumn::fractionOfAdsorbent)
      .def_rw("bulk_species_sink", &MultibedColumn::bulkSpeciesSink)
      .def_rw("reaction_physisorption_source", &MultibedColumn::reactionPhysisorptionSource)
      .def_rw("reaction_chemisorption_source", &MultibedColumn::reactionChemisorptionSource)
      .def_rw("reaction_pore_concentration_source", &MultibedColumn::reactionPoreConcentrationSource)
      .def_rw("reaction_heat", &MultibedColumn::reactionHeat)
      .def_rw("state", &MultibedColumn::state)
      .def_rw("state_dot", &MultibedColumn::stateDot)
      .def_prop_ro("concentration",
                   [](const MultibedColumn& self) { return spanToVector(std::span<const double>(self.concentration)); })
      .def_prop_ro("concentration_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.concentrationDot)); })
      .def_prop_ro("physisorption",
                   [](const MultibedColumn& self) { return spanToVector(std::span<const double>(self.physisorption)); })
      .def_prop_ro("physisorption_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.physisorptionDot)); })
      .def_prop_ro("chemisorption",
                   [](const MultibedColumn& self) { return spanToVector(std::span<const double>(self.chemisorption)); })
      .def_prop_ro("chemisorption_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.chemisorptionDot)); })
      .def_prop_ro("surface_concentration", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.surfaceConcentration)); })
      .def_prop_ro("surface_concentration_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.surfaceConcentrationDot)); })
      .def_prop_ro("pore_concentration", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.poreConcentration)); })
      .def_prop_ro("pore_concentration_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.poreConcentrationDot)); })
      .def_prop_ro("surface_pore_transport_enabled",
                   [](const MultibedColumn& self) { return self.surfacePoreTransportEnabled; })
      .def_prop_ro("gas_temperature", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.gasTemperature)); })
      .def_prop_ro("gas_temperature_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.gasTemperatureDot)); })
      .def_prop_ro("solid_temperature", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.solidTemperature)); })
      .def_prop_ro("solid_temperature_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.solidTemperatureDot)); })
      .def_prop_ro("wall_temperature", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.wallTemperature)); })
      .def_prop_ro("wall_temperature_dot", [](const MultibedColumn& self)
                   { return spanToVector(std::span<const double>(self.wallTemperatureDot)); })
      .def("state_size", &MultibedColumn::stateSize)
      .def("bind_state_views", &MultibedColumn::bindStateViews)
      .def("initialize", &MultibedColumn::initialize)
      .def("set_temperature", &MultibedColumn::setTemperature, nb::arg("temperature"))
      .def("__repr__", &MultibedColumn::repr)
      .def("write_json", &MultibedColumn::writeJSON, nb::arg("file_name"))
      .def("read_json", &MultibedColumn::readJSON, nb::arg("file_name"));

  m.def("update_velocity_and_pressure", &RK3Helpers::updateVelocityAndPressure, nb::arg("column"));

  m.def("compute_equilibrium_loadings", &RK3Helpers::computeEquilibriumLoadings, nb::arg("column"));

  m.def("compute_first_derivatives", &RK3Helpers::computeDerivatives, nb::arg("column"));

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
          nb::arg("column"), nb::arg("step"))
      .def(
          "propagate",
          [](RungeKutta3& self, MultibedColumn& column, std::size_t step)
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
      .def("reinitialize", &CVODE::reinitialize)
      .def(
          "initialize", [](CVODE& self, Column& column) { self.initialize(column); }, nb::arg("column"))
      .def(
          "initialize", [](CVODE& self, MultibedColumn& column) { self.initialize(column); }, nb::arg("column"))
      .def(
          "propagate",
          [](CVODE& self, Column& column, std::size_t step)
          {
            Timing timings;
            return self.propagate(column, step, timings);
          },
          nb::arg("column"), nb::arg("step"))
      .def(
          "propagate",
          [](CVODE& self, MultibedColumn& column, std::size_t step)
          {
            Timing timings;
            return self.propagate(column, step, timings);
          },
          nb::arg("column"), nb::arg("step"));

  nb::enum_<BreakthroughIntegrationScheme>(m, "BreakthroughIntegrationScheme", nb::is_arithmetic())
      .value("SSP_RK", BreakthroughIntegrationScheme::SSP_RK)
      .value("CVODE", BreakthroughIntegrationScheme::CVODE)
      .value("SIRK3", BreakthroughIntegrationScheme::SIRK3);

  bindBreakthrough<Column>(m, "Breakthrough");
  bindBreakthrough<MultibedColumn>(m, "MultibedBreakthrough");

  nb::class_<Fitting> fitting(m, "Fitting");

  nb::enum_<Fitting::PressureScale>(fitting, "PressureScale", nb::is_arithmetic())
      .value("LOG", Fitting::PressureScale::Log)
      .value("LINEAR", Fitting::PressureScale::Linear);

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

  nb::class_<SwingAdsorptionSubStage>(m, "SwingAdsorptionSubStage")
      .def(nb::init<>())
      .def_rw("name", &SwingAdsorptionSubStage::name)
      .def_rw("temperature", &SwingAdsorptionSubStage::temperature)
      .def_rw("pressure", &SwingAdsorptionSubStage::pressure)
      .def_rw("number_of_steps", &SwingAdsorptionSubStage::numberOfSteps);

  bindSwingAdsorption<Column>(m, "SwingAdsorption");
  bindSwingAdsorption<MultibedColumn>(m, "MultibedSwingAdsorption");

  m.def("read_input", [](const std::string& fileName) { return InputReader(fileName); }, nb::arg("file_name"));

  m.def(
      "simulation_from_file",
      [](const std::string& fileName) -> nb::object
      {
        InputReader inputReader(fileName);
        switch (inputReader.simulationType)
        {
          case InputReader::SimulationType::Breakthrough:
            if (inputReader.adsorbentComponents.size() > 1)
            {
              return nb::cast(Breakthrough<MultibedColumn>(inputReader));
            }
            return nb::cast(Breakthrough<Column>(inputReader));
          case InputReader::SimulationType::MixturePrediction:
            return nb::cast(MixturePrediction(inputReader));
          case InputReader::SimulationType::Fitting:
            return nb::cast(Fitting(inputReader));
          case InputReader::SimulationType::SwingAdsorption:
            if (inputReader.adsorbentComponents.size() > 1)
            {
              return nb::cast(SwingAdsorption<MultibedColumn>(inputReader));
            }
            return nb::cast(SwingAdsorption<Column>(inputReader));
          case InputReader::SimulationType::Test:
            throw std::invalid_argument("SimulationType 'Test' does not have a Python simulation object");
        }

        throw std::invalid_argument("Unknown SimulationType");
      },
      nb::arg("file_name"));
}

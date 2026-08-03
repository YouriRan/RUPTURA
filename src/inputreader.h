#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "component.h"
#include "geometry.h"

/**
 * \brief Returns true when a string starts with the supplied prefix.
 */
extern bool startsWith(const std::string& str, const std::string& prefix);

/**
 * \brief Removes leading and trailing whitespace from a string.
 */
extern std::string trim(const std::string& s);

/**
 * \brief Parses JSON input files and stores simulation parameters.
 *
 * Supports case-insensitive keys for user-facing input fields. The reader also
 * performs basic consistency checks for breakthrough simulations.
 */
struct InputReader
{
  /**
   * \brief Constructs an InputReader and parses the given input file.
   *
   * \param fileName The name of the input file to parse.
   */
  InputReader(const std::string fileName);

  /**
   * \brief Enumerates the types of simulations supported.
   */
  enum class SimulationType
  {
    Breakthrough = 0,       ///< Breakthrough simulation.
    MixturePrediction = 1,  ///< Mixture prediction simulation.
    Fitting = 2,            ///< Fitting simulation.
    SwingAdsorption = 3,    ///< Swing adsorption simulation.
    Test = 4                ///< Test simulation.
  };

  std::vector<Component> components;  ///< The list of components involved in the simulation.
  std::vector<std::vector<Component>>
      adsorbentComponents;  ///< Per-adsorbent component lists with adsorbent-specific isotherm parameters.
  size_t numberOfCarrierGases{0};     ///< The number of carrier gas components.
  size_t carrierGasComponent{0};      ///< The index of the carrier gas component.
  size_t maxIsothermTerms{0};         ///< The maximum number of isotherm terms among all components.

  SimulationType simulationType{SimulationType::Breakthrough};  ///< The type of simulation to perform.
  size_t mixturePredictionMethod{0};                            ///< The method used for mixture prediction.
  size_t IASTMethod{0};                                         ///< The method used for IAST calculations.
  size_t breakthroughIntegrator{0};                             ///< The integrator used for breakthrough calculations.
  size_t boundaryCondition{0};                                  ///< The breakthrough boundary condition.
  std::string displayName{"Column"};                            ///< The display name for the simulation.

  double temperature{433.0};            ///< The simulation temperature in K.
  double columnVoidFraction{0.4};       ///< The void fraction of the packed column.
  double dynamicViscosity{1e-5};        ///< Dynamic viscosity of the gas phase in Pa s.
  double particleDiameter{1e-3};        ///< Diameter of the packed particles in m.
  double particleDensity{1000.0};       ///< Particle density in kg/m^3.
  double inletPressure{-1.0};           ///< Inlet pressure, P_in, in Pa.
  double outletPressure{-1.0};          ///< Outlet pressure, P_out, in Pa.
  double pressureGradient{0.0};         ///< Pressure-gradient parameter used by fixed-pressure boundary conditions.
  double columnEntranceVelocity{-1.0};  ///< Inlet velocity, v_in, in m/s.
  double columnLength{0.3};             ///< Column length in m.
  std::vector<double> columnDistances;  ///< Spatial grid node positions in m.

  std::vector<double> adsorbentLengths;            ///< Pure adsorbent-region lengths in m.
  std::vector<double> adsorbentInterfaceLengths;   ///< Linear interface lengths between adsorbents in m.
  std::vector<size_t> adsorbentGridPoints;          ///< Spatial grid intervals per adsorbent section.
  std::vector<double> adsorbentVoidFractions;      ///< Packed-bed void fractions per adsorbent.
  std::vector<double> adsorbentParticleDensities;  ///< Particle densities per adsorbent in kg/m^3.
  std::vector<double> adsorbentParticleDiameters;  ///< Particle diameters per adsorbent in m.

  Geometry geometry;                 ///< Parsed column geometry with precomputed shape parameters.
  double influxTemperature{433.0};   ///< Feed/influx gas temperature in K.
  double internalDiameter{0.0};      ///< Column internal diameter in m.
  double outerDiameter{0.0};         ///< Column outer diameter in m.
  double wallDensity{0.0};           ///< Wall density in kg/m^3.
  double gasThermalConductivity{0.0};   ///< Gas thermal conductivity in W/(m K).
  double wallThermalConductivity{0.0};  ///< Wall thermal conductivity in W/(m K).
  double heatTransferGasSolid{0.0};     ///< Gas-solid heat-transfer coefficient in W/(m^2 K).
  double heatTransferGasWall{0.0};      ///< Gas-wall heat-transfer coefficient in W/(m^2 K).
  double heatTransferWallExternal{0.0};  ///< Wall-external heat-transfer coefficient in W/(m^2 K).
  double heatCapacityGas{1.0};           ///< Gas heat capacity in J/(kg K).
  double heatCapacitySolid{1.0};         ///< Solid heat capacity in J/(kg K).
  double heatCapacityWall{1.0};          ///< Wall heat capacity in J/(kg K).
  bool energyBalance{false};             ///< Enables gas/solid/wall temperature dynamics when true.

  size_t numberOfTimeSteps{0};       ///< The number of time steps in the simulation.
  size_t numberOfInitTimeSteps{0};   ///< The number of initialization time steps.
  bool autoNumberOfTimeSteps{true};  ///< Whether to automatically determine the number of time steps.
  double timeStep{0.0005};           ///< The time step size in s.
  size_t printEvery{10000};          ///< The interval at which to print output.
  size_t writeEvery{10000};          ///< The interval at which to write output.
  size_t numberOfGridPoints{100};    ///< The number of grid points in the column.

  double pressureStart{-1.0};          ///< The starting pressure for isotherm calculations in Pa.
  double pressureEnd{-1.0};            ///< The ending pressure for isotherm calculations in Pa.
  size_t numberOfPressurePoints{100};  ///< The number of pressure points to calculate.
  size_t pressureScale{0};             ///< The pressure scale: 0 = log, 1 = linear.

  size_t columnPressure{0};  ///< The index of the column for pressure data.
  size_t columnLoading{1};   ///< The index of the column for loading data.
  size_t columnError{2};     ///< The index of the column for error data.

  std::vector<double> swingTemperatures;  ///< Swing-stage temperatures in K.
  std::vector<double> swingPressures;     ///< Swing-stage pressures in Pa.
  std::vector<size_t> swingSteps;         ///< Swing-stage step counts.

  std::optional<std::string> readColumnFile;  ///< Optional JSON file used to initialize column state.
};

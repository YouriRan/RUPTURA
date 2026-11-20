#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "component.h"
#include "geometry.h"
#include "macrostate_particle_distribution.h"
#include "reaction.h"

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
    SwingAdsorption = 3     ///< Swing adsorption simulation.
  };

  /**
   * \brief A condition that ends a swing phase before its Duration: a metric of the column at or beyond a value.
   */
  struct StopCriterion
  {
    enum class Metric
    {
      Pressure = 0,            ///< Gas pressure of the cell at the end, in Pa.
      AveragePressure = 1,     ///< Gas pressure averaged over the fluid volume of the cells, in Pa.
      MaximumTemperature = 2,  ///< Highest solid temperature of the cells, in K.
      MinimumTemperature = 3,  ///< Lowest solid temperature of the cells, in K.
      MoleFraction = 4,        ///< Mole fraction of the component in the cell at the end.
      AmountIn = 5,  ///< Net amount that entered at the end in this phase, of the component or in total, in mol/m^2.
      AmountOut = 6  ///< Net amount that left at the end in this phase, of the component or in total, in mol/m^2.
    };
    Metric metric{Metric::Pressure};
    size_t column{0};                 ///< Index of the column in columnNames.
    size_t end{1};                    ///< 0 = z = 0 (FeedEnd), 1 = z = L (ProductEnd).
    std::optional<size_t> component;  ///< Component index; none: all components (the amounts).
    bool atLeast{true};               ///< The criterion holds for metric >= value; else for metric <= value.
    double value{0.0};
  };

  /**
   * \brief One end of the column in a swing phase with BoundaryCondition MomentumBalance (Column::Port).
   */
  struct Port
  {
    size_t mode{0};                   ///< 0 = Closed, 1 = Pressure, 2 = Velocity.
    std::optional<double> pressure;   ///< Pressure of the gas beyond the end in Pa.
    std::optional<double> velocity;   ///< Interstitial velocity into the column in m/s.
    std::optional<double> molarFlux;  ///< Molar flux into the column per column cross-section in mol/(m^2 s).
    size_t pressureRampShape{0};      ///< Pressure from that of the end cell: 0 = a step, 1 = exponential, 2 = linear.
    double pressureRampTimeConstant{0.0};  ///< Time constant (exponential) or duration (linear) of the ramp in s.
  };

  /**
   * \brief A valve or compressor of a swing phase between two column ends or buffers; gas flows from the higher to the
   *        lower pressure, a positive flow from the first endpoint to the second.
   *
   * An endpoint is 2 column + end for the end (0 = z = 0, 1 = z = L) of a column, or 2 numberOfColumns + buffer.
   */
  struct Connection
  {
    enum class FlowModel
    {
      Linear = 0,      ///< conductance (P_from - P_to), conductance in mol/(s Pa).
      SquareRoot = 1,  ///< conductance sqrt(|P_from - P_to|), conductance in mol/(s Pa^0.5).
      Orifice = 2,     ///< Compressible flow through the effective area conductance in m^2, choked or not.
      Compressor = 3   ///< The molar flow conductance in mol/s from the first endpoint to the second.
    };
    std::string name;
    std::array<size_t, 2> endpoints{};
    FlowModel flowModel{FlowModel::Linear};
    double conductance{0.0};
    double opening{1.0};     ///< Fraction of the flow of the open valve.
    bool checkValve{false};  ///< Gas flows from the first endpoint to the second only.
    double maximumMolarFlow{std::numeric_limits<double>::infinity()};  ///< In mol/s.
    double efficiency{1.0};                                            ///< Isentropic efficiency of a compressor.
  };

  /**
   * \brief A well-mixed gas vessel at a constant temperature, which connections fill and empty.
   */
  struct Buffer
  {
    std::string name;
    double volume{0.0};                   ///< In m^3.
    double pressure{0.0};                 ///< Initial pressure in Pa.
    double temperature{0.0};              ///< In K.
    std::vector<double> moleFractions{};  ///< Initial composition per component.
    double maximumPressure{0.0};          ///< Gas above this pressure leaves the buffer as its export, in Pa; 0: none.
  };

  /**
   * \brief User-provided operating conditions for one swing adsorption phase.
   *
   * Temperature and inlet pressure are optional at input level; omitted values
   * inherit the base breakthrough conditions when the swing simulation is built.
   */
  struct SwingAdsorptionPhase
  {
    std::string name;                         ///< Optional phase label for logs and frontends.
    std::optional<double> temperature;        ///< Phase temperature in K.
    std::optional<double> inletPressure;      ///< Phase inlet pressure in Pa.
    double duration{0.0};                     ///< Phase duration in s.
    /// Phase feed mole fractions per component (normalized); empty means the Components feed.
    std::vector<double> gasPhaseMolFractions{};
    std::optional<double> columnEntranceVelocity;  ///< Phase inlet (interstitial) velocity in m/s.
    std::optional<size_t> flowDirection;           ///< Phase flow direction (0 Forward, 1 Reverse); default the column's.
    /// Mole fractions of the gas that flows back in at z = L in a forward phase (normalized); empty: the gas of the
    /// last cell.
    std::vector<double> productEndGasPhaseMolFractions{};
    bool closedEnd{false};  ///< The end opposite the one where flowDirection feeds is closed.
    size_t pressureRampShape{0};  ///< Inlet pressure from the previous phase: 0 = a step, 1 = exponential, 2 = linear.
    double pressureRampTimeConstant{0.0};  ///< Time constant (exponential) or duration (linear) of the ramp in s.
    /// The feed, or the product-end gas, is the mean gas that left the column in this phase of the cycle (its index in
    /// the cycle) the last time it ran, in place of mole fractions; until it has run, the Components composition.
    std::optional<size_t> feedFromPhase;
    std::optional<size_t> productEndFromPhase;
    /// The phase ends before its Duration once the stop criteria hold: all of them, or with stopOnAny any of them.
    std::vector<StopCriterion> stopCriteria{};
    bool stopOnAny{false};
    double minimumDuration{0.0};  ///< The stop criteria apply from this time into the phase, in s.
    bool errorOnTimeout{false};   ///< Reaching the Duration with the stop criteria unmet is an error.
    std::array<Port, 2> ports{};  ///< MomentumBalance: the ports at z = 0 (FeedEnd) and z = L (ProductEnd).
    /// Temperature beyond the outer wall in K and its heat-transfer coefficient in W/(m^2 K) in this phase
    /// (ThermalControl); without them the phase temperature and heatTransferWallExternal.
    std::optional<double> jacketTemperature;
    std::optional<double> jacketHeatTransfer;
    std::optional<double> coolantMassFlowRate;  ///< Of the Coolant in this phase, in kg/s (ThermalControl).
    std::vector<Connection> connections{};      ///< Between column ends and buffers; on the entry of the first column.
  };

  std::vector<Component> components;  ///< The list of components involved in the simulation.
  std::vector<Reaction> reactions;    ///< Optional reactions coupled to pore or adsorbed component balances.
  std::vector<std::vector<Component>>
      adsorbentComponents;  ///< Per-adsorbent component lists with adsorbent-specific isotherm parameters.
  size_t numberOfCarrierGases{0};     ///< The number of carrier gas components.
  size_t carrierGasComponent{0};      ///< The index of the carrier gas component.
  size_t maxIsothermTerms{0};         ///< The maximum number of isotherm terms among all components.

  SimulationType simulationType{SimulationType::Breakthrough};  ///< The type of simulation to perform.
  size_t mixturePredictionMethod{0};                            ///< The method used for mixture prediction.
  size_t IASTMethod{0};                                         ///< The method used for IAST calculations.
  std::optional<MPDSettings> mpdSettings;                        ///< Settings used by MPD mixture prediction.
  size_t breakthroughIntegrator{0};  ///< Breakthrough integrator: 0 = RungeKutta3, 1 = CVODE, 2 = ROS34PW2.
  size_t convectionScheme{1};        ///< Convective face values: 0 = Upwind, 1 = WENO3, 2 = TVD, 3 = WENO5.
  size_t boundaryCondition{0};       ///< The breakthrough boundary condition (Column::BoundaryCondition).
  size_t inletCondition{1};  ///< Species and fluid-temperature inlet: 0 = Dirichlet (feed state), 1 = Danckwerts (feed flux).
  size_t flowDirection{0};   ///< End at which the feed enters: 0 = Forward (z = 0), 1 = Reverse (z = L).
  size_t fluidPhase{0};                                         ///< 0 = gas, 1 = incompressible liquid.
  size_t pHMode{0};                                             ///< 0 = fixed, 1 = H+, 2 = OH-.
  size_t pHComponent{0};                                        ///< Zero-based acid/base component identifier.
  std::string displayName{"Column"};                            ///< The display name for the simulation.

  double temperature{298.0};            ///< The simulation temperature in K.
  double columnVoidFraction{0.4};       ///< The void fraction of the packed column.
  double dynamicViscosity{1e-5};        ///< Dynamic viscosity of the gas phase in Pa s.
  /// Reference state of the AxialDispersionCoefficient values, which then scale as (T/T_ref)^1.75 (P_ref/P);
  /// 0 means no scaling.
  double axialDispersionReferenceTemperature{0.0};
  double axialDispersionReferencePressure{0.0};  ///< See axialDispersionReferenceTemperature; in Pa.
  double liquidDensity{1000.0};         ///< Constant liquid density in kg/m^3.
  double pHValue{7.0};                  ///< Fixed pH when pHMode is fixed.
  double pKw{14.0};                     ///< Water ion-product exponent for OH--derived pH.
  /// Fixed relative humidity [-] of the humidity-dependent isotherms; without it, RH = p_H2O / p_sat(T) from the
  /// component marked Water.
  std::optional<double> relativeHumidity{};
  /// Saturation pressure of water for the relative humidity and GAB (WaterSaturationPressureModel and
  /// WaterSaturationPressure).
  WaterSaturationPressure waterSaturationPressure{};
  double particleDiameter{1e-3};        ///< Diameter of the packed particles in m.
  double particleDensity{1000.0};       ///< Particle density in kg/m^3.
  double inletPressure{-1.0};           ///< Inlet pressure, P_in, in Pa.
  double outletPressure{-1.0};          ///< Outlet pressure, P_out, in Pa.
  double pressureGradient{0.0};         ///< Signed axial pressure gradient dP/dz in Pa/m.
  double columnEntranceVelocity{-1.0};  ///< Inlet velocity, v_in, in m/s.
  double columnLength{0.3};             ///< Column length in m.
  std::vector<double> columnDistances;  ///< Spatial grid node positions in m.

  std::vector<double> adsorbentLengths;            ///< Pure adsorbent-region lengths in m; empty for a uniform mix.
  std::vector<double> adsorbentInterfaceLengths;   ///< Linear interface lengths between adsorbents in m.
  std::vector<size_t> adsorbentInterfaceGridPoints;  ///< Intervals per interface; zero for sharp interfaces.
  std::vector<double> adsorbentMixFractions;       ///< Uniform adsorbent volume fractions over the whole column.
  std::vector<size_t> adsorbentGridPoints;          ///< Spatial grid intervals per adsorbent section.
  std::vector<double> adsorbentVoidFractions;      ///< Packed-bed void fractions per adsorbent.
  std::vector<double> adsorbentParticleDensities;  ///< Particle densities per adsorbent in kg/m^3.
  std::vector<double> adsorbentParticleDiameters;  ///< Particle diameters per adsorbent in m.
  std::vector<Geometry> adsorbentGeometries;        ///< Derived geometry for each adsorbent section.

  Geometry geometry{makeGeometry(PackedBedTubeSpec{})};  ///< Parsed geometry with precomputed terms.
  double influxTemperature{298.0};   ///< Feed/influx gas temperature in K.
  double internalDiameter{0.0};      ///< Column internal diameter in m.
  double outerDiameter{0.0};         ///< Column outer diameter in m.
  double wallDensity{0.0};           ///< Wall density in kg/m^3.
  double gasThermalConductivity{0.0};   ///< Gas thermal conductivity in W/(m K).
  double wallThermalConductivity{0.0};  ///< Wall thermal conductivity in W/(m K).
  double heatTransferGasSolid{0.0};     ///< Gas-solid heat-transfer coefficient in W/(m^2 K).
  double heatTransferGasWall{0.0};      ///< Gas-wall heat-transfer coefficient in W/(m^2 K).
  double heatTransferWallExternal{0.0};  ///< Wall-external heat-transfer coefficient in W/(m^2 K).
  /// A heat-transfer fluid along the outer wall (Coolant): it enters at the temperature beyond the wall, at z = 0 or,
  /// counter-current, at z = L. Mass flow in kg/s, heat capacity in J/(kg K), density in kg/m^3, volume in m^3.
  bool coolant{false};
  double coolantMassFlowRate{0.0};
  double coolantHeatCapacity{4180.0};
  double coolantDensity{1000.0};
  double coolantVolume{0.0};
  bool coolantCounterCurrent{false};
  size_t radialCells{0};  ///< Radial cells of the particle with ParticleTransport HSDM; 0: linear driving force.
  double heatCapacityGas{1.0};           ///< Gas heat capacity in J/(kg K).
  double liquidThermalConductivity{0.6};  ///< Liquid thermal conductivity in W/(m K).
  double heatTransferLiquidSolid{0.0};    ///< Liquid-solid heat-transfer coefficient in W/(m^2 K).
  double heatTransferLiquidWall{0.0};     ///< Liquid-wall heat-transfer coefficient in W/(m^2 K).
  double heatCapacityLiquid{4180.0};      ///< Liquid heat capacity in J/(kg K).
  double heatCapacitySolid{1.0};         ///< Solid heat capacity in J/(kg K).
  double heatCapacityWall{1.0};          ///< Wall heat capacity in J/(kg K).
  bool energyBalance{false};             ///< Enables gas/solid/wall temperature dynamics when true.

  double simulationTime{0.0};       ///< Simulated time in s; set from the phase durations for SwingAdsorption.
  bool autoSimulationTime{true};    ///< Run until the breakthrough or reaction steady-state criterion is met.
  double writeInterval{5.0};        ///< Interval between outputs (and progress lines) in s.
  bool writeOutletOnly{false};      ///< The output files get the outlet node only, one row per output.
  size_t numberOfInitTimeSteps{0};  ///< RK3 steps over which the step ramps up from 1e-4 TimeStep.
  double timeStep{0.0005};          ///< RK3 time step in s; not used by the adaptive CVODE and ROS34PW2.

  // CVODE integrator controls (only used when BreakthroughIntegrator is CVODE).
  double cvodeRelativeTolerance{1.0e-5};  ///< CVODE relative tolerance, applied to every state component.
  double cvodeAbsoluteToleranceConcentration{1.0e-6};   ///< CVODE absolute tolerance for concentrations in mol/m^3.
  double cvodeAbsoluteToleranceLoading{1.0e-7};        ///< CVODE absolute tolerance for loadings in mol/kg.
  double cvodeAbsoluteToleranceTemperature{1.0e-4};     ///< CVODE absolute tolerance for temperatures in K.
  double cvodeMaximumTimeStep{-1.0};  ///< CVODE maximum internal step in s: <0 = 4 dz/v, 0 = unlimited, >0 = explicit.
  size_t cvodeMaximumOrder{5};        ///< CVODE maximum BDF order (1-5); 1 for very steep fronts that stall at higher orders.
  size_t cvodeLinearSolver{1};      ///< CVODE linear solver: 0 = Dense, 1 = SPGMR (matrix-free), 2 = KLU (sparse).
  size_t cvodeKrylovDimension{30};  ///< Krylov subspace dimension used by the SPGMR linear solver.
  size_t cvodePreconditioner{0};    ///< SPGMR preconditioner: 0 = None, 1 = SorptionJacobian (banded I - gamma J).

  // ROS34PW2 integrator controls (only used when BreakthroughIntegrator is ROS34PW2). The tolerances are ten times
  // smaller than those of CVODE, whose step cap of 4 dz/v adds accuracy that ROS34PW2, without a cap, does not get.
  double ros34pw2RelativeTolerance{1.0e-6};               ///< ROS34PW2 relative tolerance, for every state component.
  double ros34pw2AbsoluteToleranceConcentration{1.0e-7};  ///< ROS34PW2 absolute tolerance for concentrations in mol/m^3.
  double ros34pw2AbsoluteToleranceLoading{1.0e-8};        ///< ROS34PW2 absolute tolerance for loadings in mol/kg.
  double ros34pw2AbsoluteToleranceTemperature{1.0e-5};    ///< ROS34PW2 absolute tolerance for temperatures in K.
  double ros34pw2MaximumTimeStep{0.0};                    ///< ROS34PW2 largest step in s; 0 = unlimited.
  size_t numberOfGridPoints{100};    ///< The number of grid points in the column.

  double pressureStart{-1.0};          ///< The starting pressure for isotherm calculations in Pa.
  double pressureEnd{-1.0};            ///< The ending pressure for isotherm calculations in Pa.
  size_t numberOfPressurePoints{100};  ///< The number of pressure points to calculate.
  size_t pressureScale{0};             ///< The pressure scale: 0 = log, 1 = linear.

  size_t columnPressure{1};  ///< One-based number of the data-file column for pressure.
  size_t columnLoading{2};   ///< One-based number of the data-file column for loading.
  size_t columnError{3};     ///< One-based number of the data-file column for the error.
  size_t fittingWeighting{0};  ///< The fitting residual: 0 = absolute, 1 = relative.

  /// Ordered swing adsorption phases. With columnNames every phase is there once per column, in the order of the names.
  std::vector<SwingAdsorptionPhase> swingAdsorptionPhases;
  /// Names of the identical columns of a swing run that connections couple (Columns); empty: one column.
  std::vector<std::string> columnNames;
  std::vector<Buffer> buffers;  ///< Gas vessels of a swing run (Buffers).
  /// Gas that leaves through a pressure port below this pressure costs the work of its compression to it, in Pa, at
  /// the isentropic efficiency below, which is also that of the vacuum pump; 0: no such work.
  double exhaustPressure{0.0};
  double compressionEfficiency{1.0};
  size_t numberOfCycles{1};  ///< Times the phase list is repeated in swingAdsorptionPhases (NumberOfCycles).
  /// Cyclic steady state: the run ends after cyclicSteadyStateCycles cycles in a row whose state, amounts, and streams
  /// differ from those of the cycle before by less than this relative tolerance; 0 = run NumberOfCycles cycles.
  double cyclicSteadyStateTolerance{0.0};
  size_t cyclicSteadyStateCycles{1};
  /// The column integrates the amounts through its ends (Column::streamAmount): set by stop criteria, cyclic steady
  /// state, gas taken from a phase, and MomentumBalance.
  bool streamAmounts{false};

  std::optional<std::string> readColumnFile;  ///< Optional JSON file used to initialize column state.
  /// A ReadColumnFile run starts at t = 0 with new output files, instead of going on from the time of the state.
  bool resetClock{false};
};

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "component.h"
#include "geometry.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "reaction.h"
#include "utils.h"

/**
 * \brief Centralized view of the flattened Column ODE state layout.
 *
 * The component blocks are grid-major arrays of size (numberOfGridPoints + 1) * numberOfComponents.
 * Optional surface/pore concentration blocks are present only for mechanistic chemisorption.
 */
struct ColumnStateLayout
{
  size_t numberOfGridPoints{0};
  size_t numberOfComponents{0};
  size_t numberOfChemisorptionSites{1};
  bool includeSurfacePoreTransport{false};

  [[nodiscard]] size_t nodeCount() const noexcept { return numberOfGridPoints + 1; }
  [[nodiscard]] size_t componentBlockSize() const noexcept { return nodeCount() * numberOfComponents; }
  [[nodiscard]] size_t scalarBlockSize() const noexcept { return nodeCount(); }
  [[nodiscard]] size_t chemisorptionSiteCount() const noexcept
  {
    return std::max<size_t>(1, numberOfChemisorptionSites);
  }
  [[nodiscard]] size_t componentBlockCount() const noexcept
  {
    return 2 + chemisorptionSiteCount() * (includeSurfacePoreTransport ? 3 : 1);
  }
  [[nodiscard]] size_t temperatureOffset() const noexcept { return componentBlockCount() * componentBlockSize(); }
  [[nodiscard]] size_t stateSize() const noexcept { return temperatureOffset() + 3 * scalarBlockSize(); }

  [[nodiscard]] std::span<double> componentBlock(double* base, size_t block) const noexcept
  {
    return {base + block * componentBlockSize(), componentBlockSize()};
  }

  [[nodiscard]] std::span<const double> componentBlock(const double* base, size_t block) const noexcept
  {
    return {base + block * componentBlockSize(), componentBlockSize()};
  }

  [[nodiscard]] std::span<double> concentration(double* base) const noexcept { return componentBlock(base, 0); }
  [[nodiscard]] std::span<double> physisorption(double* base) const noexcept { return componentBlock(base, 1); }
  [[nodiscard]] std::span<double> chemisorption(double* base) const noexcept
  {
    return {base + 2 * componentBlockSize(), chemisorptionSiteCount() * componentBlockSize()};
  }
  [[nodiscard]] std::span<double> surfaceConcentration(double* base) const noexcept
  {
    return includeSurfacePoreTransport
               ? std::span<double>{base + (2 + chemisorptionSiteCount()) * componentBlockSize(),
                                   chemisorptionSiteCount() * componentBlockSize()}
               : std::span<double>{};
  }
  [[nodiscard]] std::span<double> poreConcentration(double* base) const noexcept
  {
    return includeSurfacePoreTransport
               ? std::span<double>{base + (2 + 2 * chemisorptionSiteCount()) * componentBlockSize(),
                                   chemisorptionSiteCount() * componentBlockSize()}
               : std::span<double>{};
  }
  [[nodiscard]] std::span<double> gasTemperature(double* base) const noexcept
  {
    return {base + temperatureOffset(), scalarBlockSize()};
  }
  [[nodiscard]] std::span<double> solidTemperature(double* base) const noexcept
  {
    return {base + temperatureOffset() + scalarBlockSize(), scalarBlockSize()};
  }
  [[nodiscard]] std::span<double> wallTemperature(double* base) const noexcept
  {
    return {base + temperatureOffset() + 2 * scalarBlockSize(), scalarBlockSize()};
  }
};

/**
 * \brief Packed-column model state and configuration.
 *
 * Stores column parameters, component data, cache arrays, scratch arrays, and
 * canonical ODE state storage.
 */
struct Column
{
  /**
   * \brief Boundary-condition pair supplied by the input file.
   */
  enum class BoundaryCondition
  {
    InletPressureInletVelocity = 0,   ///< Boundary data: P_in and v_in.
    InletPressureOutletPressure = 1,  ///< Boundary data: P_in and P_out.
    InletVelocityOutletPressure = 2,  ///< Boundary data: v_in and P_out.
    FixedVelocity = 3,                ///< Boundary data: fixed v, with P_in or P_out.
    FixedPressureInletVelocity = 4    ///< Boundary data: fixed pressure profile and v_in.
  };

  /**
   * \brief Constructs a column from explicit model/configuration arguments.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  Column(MixturePrediction physisorptionMixture, std::vector<Component> components,
         BoundaryCondition boundaryCondition,
         bool energyBalance, size_t numberOfGridPoints, size_t maxIsothermTerms, size_t carrierGasComponent,
         double temperature, double inletPressure, double outletPressure, double pressureGradient,
         double columnVoidFraction, double particleDensity, double columnEntranceVelocity, double columnLength,
         double dynamicViscosity, double particleDiameter, double influxTemperature, double internalDiameter,
         double outerDiameter, double wallDensity, double gasThermalConductivity, double wallThermalConductivity,
         double heatTransferGasSolid, double heatTransferGasWall, double heatTransferWallExternal,
         double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
         std::vector<Reaction> reactions = {})
      : Column(std::move(physisorptionMixture), std::move(components), boundaryCondition, energyBalance,
               numberOfGridPoints, maxIsothermTerms, carrierGasComponent, temperature, inletPressure,
               outletPressure, pressureGradient, columnVoidFraction, particleDensity, columnEntranceVelocity,
               columnLength, dynamicViscosity, particleDiameter, influxTemperature, internalDiameter,
               outerDiameter, wallDensity, gasThermalConductivity, wallThermalConductivity,
               heatTransferGasSolid, heatTransferGasWall, heatTransferWallExternal, heatCapacityGas,
               heatCapacitySolid, heatCapacityWall,
               makeGeometry(PackedBedTubeSpec{.voidFraction = columnVoidFraction,
                                              .particleDiameter = particleDiameter,
                                              .internalDiameter = internalDiameter,
                                              .outerDiameter = outerDiameter}),
               std::move(reactions))
  {
  }

  Column(MixturePrediction physisorptionMixture, std::vector<Component> components,
         BoundaryCondition boundaryCondition,
         bool energyBalance, size_t numberOfGridPoints, size_t maxIsothermTerms, size_t carrierGasComponent,
         double temperature, double inletPressure, double outletPressure, double pressureGradient,
         double columnVoidFraction, double particleDensity, double columnEntranceVelocity, double columnLength,
         double dynamicViscosity, double particleDiameter, double influxTemperature, double internalDiameter,
         double outerDiameter, double wallDensity, double gasThermalConductivity, double wallThermalConductivity,
         double heatTransferGasSolid, double heatTransferGasWall, double heatTransferWallExternal,
         double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall, Geometry geometry,
         std::vector<Reaction> reactions = {})
      : physisorptionMixture(std::move(physisorptionMixture)),
        components(std::move(components)),
        boundaryCondition(boundaryCondition),
        energyBalance(energyBalance),
        geometry(std::move(geometry)),
        reactions(std::move(reactions)),
        numberOfGridPoints(numberOfGridPoints),
        numberOfComponents(this->components.size()),
        maxIsothermTerms(maxIsothermTerms),
        maxChemisorptionSites(maximumChemisorptionSites(this->components)),
        numberOfCalls(0),
        carrierGasComponent(carrierGasComponent),
        externalTemperature(temperature),
        inletPressure(inletPressure),
        outletPressure(outletPressure),
        pressureGradient(pressureGradient),
        voidFraction(columnVoidFraction),
        particleDensity(particleDensity),
        columnEntranceVelocity(columnEntranceVelocity),
        columnLength(columnLength),
        dynamicViscosity(dynamicViscosity),
        particleDiameter(particleDiameter),
        influxTemperature(influxTemperature),
        internalDiameter(internalDiameter),
        outerDiameter(outerDiameter),
        wallDensity(wallDensity),
        gasThermalConductivity(gasThermalConductivity),
        wallThermalConductivity(wallThermalConductivity),
        heatTransferGasSolid(heatTransferGasSolid),
        heatTransferGasWall(heatTransferGasWall),
        heatTransferWallExternal(heatTransferWallExternal),
        heatCapacityGas(heatCapacityGas),
        heatCapacitySolid(heatCapacitySolid),
        heatCapacityWall(heatCapacityWall),
        resolution(this->columnLength / static_cast<double>(this->numberOfGridPoints)),
        timeNormalizationFactor(this->columnEntranceVelocity / this->columnLength),
        surfacePoreTransportEnabled(requiresSurfacePoreTransport(this->components) ||
                                    reactionsRequirePoreConcentration(this->reactions)),
        prefactorMassTransfer(this->numberOfComponents),
        idealGasMolFractions(this->numberOfComponents),
        adsorbedMolFractions(this->numberOfComponents),
        numberOfMolecules(this->numberOfComponents),
        interstitialGasVelocity(this->numberOfGridPoints + 1),
        gasDensity(this->numberOfGridPoints + 1),
        totalConcentration(this->numberOfGridPoints + 1),
        totalPressure(this->numberOfGridPoints + 1),
        moleFraction((this->numberOfGridPoints + 1) * this->numberOfComponents),
        partialPressure((this->numberOfGridPoints + 1) * this->numberOfComponents),
        equilibriumPhysisorption((this->numberOfGridPoints + 1) * this->numberOfComponents),
        equilibriumChemisorption(this->maxChemisorptionSites *
                                 (this->numberOfGridPoints + 1) * this->numberOfComponents),
        cachedPressure((this->numberOfGridPoints + 1) * this->numberOfComponents * this->maxIsothermTerms),
        cachedGrandPotential((this->numberOfGridPoints + 1) * this->maxIsothermTerms),
        cachedChemisorptionPressure(this->maxChemisorptionSites *
                                    (this->numberOfGridPoints + 1) * this->numberOfComponents),
        cachedChemisorptionGrandPotential(this->maxChemisorptionSites *
                                          (this->numberOfGridPoints + 1)),
        coeffDiffusion(this->numberOfGridPoints + 1),
        facePressures(this->numberOfGridPoints),
        massFlux((this->numberOfGridPoints + 1) * this->numberOfComponents),
        bulkSpeciesSink((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionPhysisorptionSource((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionChemisorptionSource(this->maxChemisorptionSites *
                                    (this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionPoreConcentrationSource(this->maxChemisorptionSites *
                                        (this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionHeat(this->numberOfGridPoints + 1),
        state(ColumnStateLayout{this->numberOfGridPoints, this->numberOfComponents,
                                this->maxChemisorptionSites,
                                this->surfacePoreTransportEnabled}
                  .stateSize(),
              0.0),
        stateDot(ColumnStateLayout{this->numberOfGridPoints, this->numberOfComponents,
                                   this->maxChemisorptionSites,
                                   this->surfacePoreTransportEnabled}
                     .stateSize(),
                 0.0)
  {
    bindStateViews();
    chemisorptionMixture = makeChemisorptionMixture(this->physisorptionMixture, this->components);
  }

  /**
   * \brief Constructs a column from parsed input data.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  Column(const InputReader& inputReader)
      : Column(MixturePrediction(inputReader), inputReader.components, BoundaryCondition(inputReader.boundaryCondition),
               inputReader.energyBalance, inputReader.numberOfGridPoints, inputReader.maxIsothermTerms,
               inputReader.carrierGasComponent, inputReader.temperature, inputReader.inletPressure,
               inputReader.outletPressure, inputReader.pressureGradient, inputReader.columnVoidFraction,
               inputReader.particleDensity, inputReader.columnEntranceVelocity, inputReader.columnLength,
               inputReader.dynamicViscosity, inputReader.particleDiameter, inputReader.influxTemperature,
               inputReader.internalDiameter, inputReader.outerDiameter, inputReader.wallDensity,
               inputReader.gasThermalConductivity, inputReader.wallThermalConductivity,
               inputReader.heatTransferGasSolid, inputReader.heatTransferGasWall,
               inputReader.heatTransferWallExternal, inputReader.heatCapacityGas,
               inputReader.heatCapacitySolid, inputReader.heatCapacityWall, inputReader.geometry,
               inputReader.reactions)
  {
  }
  /**
   * \brief Copy constructor; rebinds spans to this object's state storage.
   */
  Column(const Column& other);

  /**
   * \brief Copy assignment; rebinds spans to this object's state storage.
   */
  Column& operator=(const Column& other);

  // Model configuration and component data.
  MixturePrediction physisorptionMixture;  ///< Competitive physisorption equilibrium model.
  MixturePrediction chemisorptionMixture;  ///< Multisite competitive chemisorption equilibrium model.
  std::vector<Component> components;    ///< Component definitions and isotherm parameters; size numberOfComponents.
  BoundaryCondition boundaryCondition;  ///< Selected breakthrough boundary-condition pair.
  bool energyBalance;                   ///< Enables gas/solid/wall temperature dynamics when true.
  Geometry geometry;                    ///< Column/tube/monolith geometry and precomputed shape terms.
  std::vector<Reaction> reactions;      ///< Optional reactions coupled to adsorbed or pore-phase variables.

  // Dimensions and counters.
  size_t numberOfGridPoints;   ///< Number of spatial grid intervals; node count is numberOfGridPoints + 1.
  size_t numberOfComponents;   ///< Number of gas components.
  size_t maxIsothermTerms;     ///< Maximum number of isotherm sites across all components.
  size_t maxChemisorptionSites;  ///< Maximum number of chemisorption sites across all components.
  size_t numberOfCalls;        ///< Counter for model/evaluation calls.
  size_t carrierGasComponent;  ///< Index of the carrier-gas component.

  // Column operating conditions and geometry.
  double externalTemperature;     ///< External/reference gas temperature, T, in K.
  double inletPressure;           ///< Inlet pressure, P_in, in Pa.
  double outletPressure;          ///< Outlet pressure, P_out, in Pa.
  double pressureGradient;        ///< Input pressure-gradient parameter in Pa/m.
  double voidFraction;            ///< Packed-bed void fraction, epsilon.
  double particleDensity;         ///< Particle density in kg/m^3.
  double columnEntranceVelocity;  ///< Inlet/interstitial velocity, v_in, in m/s.
  double columnLength;            ///< Column length, L, in m.
  double dynamicViscosity;        ///< Gas dynamic viscosity used in Ergun calculations in Pa s.
  double particleDiameter;        ///< Particle diameter used in Ergun calculations in m.

  // Energy-balance parameters.
  double influxTemperature;         ///< Feed/influx gas temperature in K.
  double internalDiameter;          ///< Column internal diameter in m.
  double outerDiameter;             ///< Column outer diameter in m.
  double wallDensity;               ///< Wall density in kg/m^3.
  double gasThermalConductivity;    ///< Gas thermal conductivity in W/(m K).
  double wallThermalConductivity;   ///< Wall thermal conductivity in W/(m K).
  double heatTransferGasSolid;      ///< Gas-solid heat-transfer coefficient in W/(m^2 K).
  double heatTransferGasWall;       ///< Gas-wall heat-transfer coefficient in W/(m^2 K).
  double heatTransferWallExternal;  ///< Wall-external heat-transfer coefficient in W/(m^2 K).
  double heatCapacityGas;           ///< Gas heat capacity in J/(kg K).
  double heatCapacitySolid;         ///< Solid heat capacity in J/(kg K).
  double heatCapacityWall;          ///< Wall heat capacity in J/(kg K).

  // Derived scalar quantities.
  double resolution;               ///< Spatial grid spacing, dz = L / numberOfGridPoints, in m.
  double timeNormalizationFactor;  ///< Dimensionless-time factor, v_in / L, in 1/s.
  bool surfacePoreTransportEnabled;  ///< Adds surface/pore concentration state blocks when true.

  std::pair<size_t, size_t> iastPerformance{0, 0};  ///< Accumulated mixture-prediction diagnostics.

  // Size numberOfComponents. Indexed as value[comp].
  std::vector<double> prefactorMassTransfer;  ///< Per-component mass-transfer prefactor in kg/m^3/s.
  std::vector<double> idealGasMolFractions;   ///< Temporary gas-phase mole fractions.
  std::vector<double> adsorbedMolFractions;   ///< Temporary adsorbed-phase mole fractions.
  std::vector<double> numberOfMolecules;      ///< Temporary equilibrium loading result per component.

  // Size numberOfGridPoints + 1. Indexed as value[grid].
  std::vector<double> interstitialGasVelocity;  ///< Interstitial gas velocity at each grid node in m/s.
  std::vector<double> gasDensity;               ///< Gas density at each grid node in kg/m^3.
  std::vector<double> totalConcentration;       ///< Total gas concentration at each grid node in mol/m^3.
  std::vector<double> totalPressure;            ///< Total pressure at each grid node in Pa.

  // Size (numberOfGridPoints + 1) * numberOfComponents. Grid-major index: grid * numberOfComponents + comp.
  std::vector<double> moleFraction;          ///< Derived gas-phase mole fraction y_i.
  std::vector<double> partialPressure;        ///< Component partial pressure at each grid node in Pa.
  std::vector<double> equilibriumPhysisorption;  ///< Equilibrium physisorbed loading in mol/kg.
  std::vector<double> equilibriumChemisorption;  ///< Site-major equilibrium chemisorbed loading in mol/kg.

  // Mixture-prediction cache arrays.
  // cachedPressure: (numberOfGridPoints + 1) * numberOfComponents * maxIsothermTerms; cachedGrandPotential:
  // (numberOfGridPoints + 1) * maxIsothermTerms.
  std::vector<double> cachedPressure;        ///< Cached hypothetical pressures for mixture prediction.
  std::vector<double> cachedGrandPotential;  ///< Cached reduced grand potentials for mixture prediction.
  std::vector<double> cachedChemisorptionPressure;  ///< Site-major chemical-phase pressure cache.
  std::vector<double> cachedChemisorptionGrandPotential;  ///< Site-major chemical-phase spreading-pressure cache.

  // Scratch/work arrays.
  // coeffDiffusion is size numberOfGridPoints + 1; facePressures is size numberOfGridPoints; massFlux is grid-major
  // size (numberOfGridPoints + 1) * numberOfComponents.
  std::vector<double> coeffDiffusion;  ///< Temporary diffusion coefficient field.
  std::vector<double> facePressures;   ///< Pressure values at cell faces in Pa.
  std::vector<double> massFlux;        ///< Component mass flux at each grid node in mol/(m^2 s).
  std::vector<double> bulkSpeciesSink;  ///< Species sink in the bulk concentration equation, mol/(m^3 s).
  std::vector<double> reactionPhysisorptionSource;  ///< Reaction-only physisorbed loading source, mol/(kg s).
  std::vector<double> reactionChemisorptionSource;  ///< Reaction-only chemisorbed loading source, mol/(kg s).
  std::vector<double> reactionPoreConcentrationSource;  ///< Reaction-only pore concentration source, mol/(m^3 s).
  std::vector<double> reactionHeat;  ///< Reaction heat release on a solid-mass basis, J/(kg s).

  // Canonical ODE storage.
  // Layout: concentration, physisorption, chemisorption, optional surface concentration, optional pore concentration,
  // gas T, solid T, wall T.
  std::vector<double> state;     ///< Canonical ODE state vector.
  std::vector<double> stateDot;  ///< Time derivative of the canonical ODE state vector.

  // Views into state/stateDot. Non-owning; rebind after copy/assignment.
  // Size (numberOfGridPoints + 1) * numberOfComponents. Grid-major index: grid * numberOfComponents + comp.
  std::span<double> concentration;  ///< Gas concentration c_i in mol/m^3; size (numberOfGridPoints + 1) * numberOfComponents.
  std::span<double> concentrationDot;  ///< Time derivative dc_i/dt in mol/(m^3 s).
  std::span<double> physisorption;  ///< Physical adsorbed loading in mol/kg.
  std::span<double> physisorptionDot;  ///< Time derivative of physical adsorbed loading in mol/(kg s).
  std::span<double> chemisorption;     ///< Chemical adsorbed loading in mol/kg.
  std::span<double> chemisorptionDot;  ///< Time derivative of chemical adsorbed loading in mol/(kg s).
  std::span<double> surfaceConcentration;     ///< Surface concentration c_s in mol/m^3 for mechanistic chemisorption.
  std::span<double> surfaceConcentrationDot;  ///< Time derivative dc_s/dt in mol/(m^3 s).
  std::span<double> poreConcentration;        ///< Pore concentration c_p in mol/m^3 for mechanistic chemisorption.
  std::span<double> poreConcentrationDot;     ///< Time derivative dc_p/dt in mol/(m^3 s).
  std::span<double> gasTemperature;       ///< Gas temperature in K; size numberOfGridPoints + 1.
  std::span<double> gasTemperatureDot;    ///< Time derivative of gas temperature in K/s; size numberOfGridPoints + 1.
  std::span<double> solidTemperature;     ///< Solid temperature in K; size numberOfGridPoints + 1.
  std::span<double> solidTemperatureDot;  ///< Time derivative of solid temperature in K/s; size numberOfGridPoints + 1.
  std::span<double> wallTemperature;      ///< Wall temperature in K; size numberOfGridPoints + 1.
  std::span<double> wallTemperatureDot;   ///< Time derivative of wall temperature in K/s; size numberOfGridPoints + 1.

  /**
   * \brief Returns the canonical ODE state size.
   */
  size_t stateSize() const noexcept;

  /**
   * \brief Returns true when any component uses mechanistic surface/pore chemisorption.
   */
  static bool requiresSurfacePoreTransport(const std::vector<Component>& components) noexcept;

  /**
   * \brief Returns the maximum number of chemisorption sites on any component.
   */
  static size_t maximumChemisorptionSites(const std::vector<Component>& components) noexcept;

  /**
   * \brief Builds the multisite competitive chemisorption equilibrium model.
   */
  static MixturePrediction makeChemisorptionMixture(
      const MixturePrediction& physisorptionMixture, const std::vector<Component>& components);

  /**
   * \brief Returns the centralized state layout for this column.
   */
  ColumnStateLayout stateLayout() const noexcept;

  /**
   * \brief Rebinds all state/stateDot spans to the current vector storage.
   */
  void bindStateViews() noexcept;

  /**
   * \brief Initializes pressure, velocity, mole fractions, loadings, and temperatures.
   */
  void initialize();

  /**
   * \brief Sets the operating temperature and fixed temperature fields when energy balance is disabled.
   */
  void setTemperature(double temperature);

  /**
   * \brief Writes headers for component output files and the column output file.
   */
  void writeOutputHeader(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream) const;

  /**
   * \brief Writes one output snapshot at the given simulation time.
   */
  void writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream, double time) const;

  /**
   * \brief Returns a compact text representation of key column settings.
   */
  std::string repr() const;

  /**
   * \brief Serializes selected column state and cache arrays to JSON.
   */
  void writeJSON(const std::string& fileName) const;

  /**
   * \brief Loads selected column state and cache arrays from JSON.
   */
  void readJSON(const std::string& fileName);
};

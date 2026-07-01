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
#include "inputreader.h"
#include "mixture_prediction.h"

/**
 * \brief Packed-column model state and configuration.
 *
 * Stores column parameters, component data, cache arrays, scratch arrays, and
 * canonical ODE state storage.
 */
struct Column
{
  /**
   * \brief Velocity-profile model used for pressure/velocity calculations.
   */
  enum class VelocityProfile
  {
    FixedPressureGradient = 0,  ///< Use the input pressure-gradient parameter.
    Ergun = 1,                  ///< Compute velocity from pressure gradient using Ergun relation.
    FixedVelocity = 2           ///< Use the input inlet velocity everywhere.
  };

  /**
   * \brief Boundary-condition pair supplied by the input file.
   */
  enum class BoundaryCondition
  {
    InletPressureInletVelocity = 0,   ///< Boundary data: P_in and v_in.
    InletPressureOutletPressure = 1,  ///< Boundary data: P_in and P_out.
    InletVelocityOutletPressure = 2   ///< Boundary data: v_in and P_out.
  };

  /**
   * \brief Constructs a column from explicit model/configuration arguments.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  Column(std::vector<MixturePrediction> mixture, std::vector<Component> components, VelocityProfile velocityProfile,
         BoundaryCondition boundaryCondition, bool energyBalance, size_t numberOfGridPoints, size_t maxIsothermTerms,
         size_t carrierGasComponent, double temperature, double inletPressure, double outletPressure,
         double pressureGradient, double columnVoidFraction, double particleDensity, double columnEntranceVelocity,
         double columnLength, double dynamicViscosity, double particleDiameter, double influxTemperature,
         double internalDiameter, double outerDiameter, double wallDensity, double gasThermalConductivity,
         double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
         double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall)
      : mixture(std::move(mixture)),
        components(std::move(components)),
        velocityProfile(velocityProfile),
        boundaryCondition(boundaryCondition),
        energyBalance(energyBalance),
        Ngrid(numberOfGridPoints),
        Ncomp(this->components.size()),
        Nads(this->mixture.size()),
        maxIsothermTerms(maxIsothermTerms),
        numCalls(0),
        carrierGasComponent(carrierGasComponent),
        externalTemperature(temperature),
        inletPressure(inletPressure),
        outletPressure(outletPressure),
        pressureGradient(pressureGradient),
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
        resolution(this->columnLength / static_cast<double>(this->Ngrid)),
        timeNormalizationFactor(this->columnEntranceVelocity / this->columnLength),
        prefactorMassTransfer(this->Ncomp),
        idealGasMolFractions(this->Ncomp),
        adsorbedMolFractions(this->Ncomp),
        numberOfMolecules(this->Ncomp),
        interstitialGasVelocity(this->Ngrid + 1),
        gasDensity(this->Ngrid + 1),
        totalConcentration(this->Ngrid + 1),
        totalPressure(this->Ngrid + 1),
        partialPressure((this->Ngrid + 1) * this->Ncomp),
        equilibriumAdsorption((this->Ngrid + 1) * this->Ncomp),
        moleFraction((this->Ngrid + 1) * this->Ncomp),
        fractionOfAdsorbent((this->Ngrid + 1) * this->Nads),
        hasAdsorbentOfType((this->Ngrid + 1) * this->Nads),
        cachedPressure((this->Ngrid + 1) * this->Ncomp * this->maxIsothermTerms),
        cachedGrandPotential((this->Ngrid + 1) * this->maxIsothermTerms),
        coeffGasGas(this->Ngrid + 1),
        coeffGasSolid(this->Ngrid + 1),
        coeffGasWall(this->Ngrid + 1),
        coeffDiffusion(this->Ngrid + 1),
        facePressures(this->Ngrid),
        massFlux((this->Ngrid + 1) * this->Ncomp),
        state((2 * this->Ncomp + 3) * (this->Ngrid + 1), 0.0),
        stateDot((2 * this->Ncomp + 3) * (this->Ngrid + 1), 0.0)
  {
    bindStateViews();
  }

  /**
   * \brief Constructs a column from parsed input data.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  Column(const InputReader& inputReader)
      : Column(MixturePrediction(inputReader), inputReader.components, VelocityProfile(inputReader.velocityProfile),
               BoundaryCondition(inputReader.boundaryCondition), inputReader.energyBalance,
               inputReader.numberOfGridPoints, inputReader.maxIsothermTerms, inputReader.carrierGasComponent,
               inputReader.temperature, inputReader.inletPressure, inputReader.outletPressure,
               inputReader.pressureGradient, inputReader.columnVoidFraction, inputReader.particleDensity,
               inputReader.columnEntranceVelocity, inputReader.columnLength, inputReader.dynamicViscosity,
               inputReader.particleDiameter, inputReader.influxTemperature, inputReader.internalDiameter,
               inputReader.outerDiameter, inputReader.wallDensity, inputReader.gasThermalConductivity,
               inputReader.wallThermalConductivity, inputReader.heatTransferGasSolid, inputReader.heatTransferGasWall,
               inputReader.heatTransferWallExternal, inputReader.heatCapacityGas, inputReader.heatCapacitySolid,
               inputReader.heatCapacityWall)
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
  std::vector<MixturePrediction> mixture;            ///< Mixture-prediction object used for adsorption equilibrium.
  std::vector<Component> components;    ///< Component definitions and isotherm parameters; size Ncomp.
  VelocityProfile velocityProfile;      ///< Selected velocity-profile model.
  BoundaryCondition boundaryCondition;  ///< Selected breakthrough boundary-condition pair.
  bool energyBalance;                   ///< Enables gas/solid/wall temperature dynamics when true.

  // Dimensions and counters.
  size_t Ngrid;                ///< Number of spatial grid intervals; node count is Ngrid + 1.
  size_t Ncomp;                ///< Number of gas components.
  size_t Nads;                 ///< Number of Adsorbents.
  size_t maxIsothermTerms;     ///< Maximum number of isotherm sites across all components.
  size_t numCalls;             ///< Counter for model/evaluation calls.
  size_t carrierGasComponent;  ///< Index of the carrier-gas component.

  // Column operating conditions and geometry.
  double externalTemperature;     ///< External/reference gas temperature, T, in K.
  double inletPressure;           ///< Inlet pressure, P_in, in Pa.
  double outletPressure;          ///< Outlet pressure, P_out, in Pa.
  double pressureGradient;        ///< Input pressure-gradient parameter.
  double columnEntranceVelocity;  ///< Inlet/interstitial velocity, v_in, in m/s.
  double columnLength;            ///< Column length, L, in m.
  double dynamicViscosity;        ///< Gas dynamic viscosity used in Ergun calculations.
  double particleDiameter;        ///< Particle diameter used in Ergun calculations.

  // Energy-balance parameters.
  double influxTemperature;         ///< Feed/influx gas temperature in K.
  double internalDiameter;          ///< Column internal diameter in m.
  double outerDiameter;             ///< Column outer diameter in m.
  double wallDensity;               ///< Wall density in kg/m^3.
  double gasThermalConductivity;    ///< Gas thermal conductivity.
  double wallThermalConductivity;   ///< Wall thermal conductivity.
  double heatTransferGasSolid;      ///< Gas-solid heat-transfer coefficient.
  double heatTransferGasWall;       ///< Gas-wall heat-transfer coefficient.
  double heatTransferWallExternal;  ///< Wall-external heat-transfer coefficient.
  double heatCapacityGas;           ///< Gas heat capacity.
  double heatCapacitySolid;         ///< Solid heat capacity.
  double heatCapacityWall;          ///< Wall heat capacity.

  // Derived scalar quantities.
  double resolution;               ///< Spatial grid spacing, dz = L / Ngrid.
  double timeNormalizationFactor;  ///< Dimensionless-time factor, v_in / L.

  std::pair<size_t, size_t> iastPerformance{0, 0};  ///< Accumulated mixture-prediction diagnostics.

  // Size Ncomp. Indexed as value[comp].
  std::vector<double> prefactorMassTransfer;  ///< Per-component mass-transfer prefactor.
  std::vector<double> idealGasMolFractions;   ///< Temporary gas-phase mole fractions.
  std::vector<double> adsorbedMolFractions;   ///< Temporary adsorbed-phase mole fractions.
  std::vector<double> numberOfMolecules;      ///< Temporary equilibrium loading result per component.

  // Size Ngrid + 1. Indexed as value[grid].
  std::vector<double> interstitialGasVelocity;  ///< Interstitial gas velocity at each grid node.
  std::vector<double> gasDensity;               ///< Gas density at each grid node.
  std::vector<double> totalConcentration;       ///< Total gas concentration at each grid node.
  std::vector<double> totalPressure;            ///< Total pressure at each grid node.
  std::vector<double> totalVoidFraction;            ///< Packed-bed void fraction, epsilon.
  std::vector<double> particleDensity;         ///< Particle density in kg/m^3.

  // Size Nads, Indexed as value[ads].
  std::vector<double> adsorbentSpecificVoidFraction;
  std::vector<double> particleDiameters;
  
  // Size (Ngrid + 1) * Ncomp. Grid-major index: grid * Ncomp + comp.
  std::vector<double> partialPressure;        ///< Component partial pressure at each grid node.
  std::vector<double> equilibriumAdsorption;  ///< Component equilibrium loading at each grid node.
  std::vector<double> moleFraction;           ///< Component gas-phase mole fraction at each grid node.

  // Size (Ngrid + 1) * Nads. Grid-major index: ads * Ncomp + comp.
  std::vector<double> fractionOfAdsorbent; ///< Fraction of the column that has this adsorbent
  std::vector<bool> hasAdsorbentOfType;  ///< Boolean switch to determine if this adsorbent is here.
  std::vector<double> adsorbentScaledVoidFractions; ///< fractionOfAdsorbent * adsorbentSpecificVoidFraction / totalVoidFraction

  // Mixture-prediction cache arrays.
  // cachedPressure: (Ngrid + 1) * Nads * Ncomp * maxIsothermTerms
  // Grid-major index: (grid * Nads + ads) * Ncomp * maxIsothermTerms
  std::vector<double> cachedPressure;        ///< Cached hypothetical pressures for mixture prediction.
  // cachedGrandPotential: (Ngrid + 1) * Nads * maxIsothermTerms. 
  // Grid-major index: (grid * Nads + ads) * maxIsothermTerms
  std::vector<double> cachedGrandPotential;  ///< Cached reduced grand potentials for mixture prediction.

  // Scratch/work arrays.
  // coeff* arrays are size Ngrid + 1; facePressures is size Ngrid; massFlux is grid-major size (Ngrid + 1) * Ncomp.
  std::vector<double> coeffGasGas;     ///< Temporary gas-gas heat-transfer coefficient field.
  std::vector<double> coeffGasSolid;   ///< Temporary gas-solid heat-transfer coefficient field.
  std::vector<double> coeffGasWall;    ///< Temporary gas-wall heat-transfer coefficient field.
  std::vector<double> coeffDiffusion;  ///< Temporary diffusion coefficient field.
  std::vector<double> facePressures;   ///< Pressure values at cell faces.
  std::vector<double> massFlux;        ///< Component mass flux at each grid node.

  // Canonical ODE storage.
  // Size: (2 * Ncomp + 3) * (Ngrid + 1). Layout: concentration, adsorption, gas T, solid T, wall T.
  std::vector<double> state;     ///< Canonical ODE state vector.
  std::vector<double> stateDot;  ///< Time derivative of the canonical ODE state vector.

  // Views into state/stateDot. Non-owning; rebind after copy/assignment.
  // Size (Ngrid + 1) * Ncomp. Grid-major index: grid * Ncomp + comp.
  std::span<double> concentration;        ///< Gas concentration c_i; size (Ngrid + 1) * Ncomp.
  std::span<double> concentrationDot;     ///< Time derivative dc_i/dt; size (Ngrid + 1) * Ncomp.
  std::span<double> adsorption;           ///< Adsorbed loading q_i; size (Ngrid + 1) * Ncomp.
  std::span<double> adsorptionDot;        ///< Time derivative dq_i/dt; size (Ngrid + 1) * Ncomp.
  std::span<double> gasTemperature;       ///< Gas temperature; size Ngrid + 1.
  std::span<double> gasTemperatureDot;    ///< Time derivative of gas temperature; size Ngrid + 1.
  std::span<double> solidTemperature;     ///< Solid temperature; size Ngrid + 1.
  std::span<double> solidTemperatureDot;  ///< Time derivative of solid temperature; size Ngrid + 1.
  std::span<double> wallTemperature;      ///< Wall temperature; size Ngrid + 1.
  std::span<double> wallTemperatureDot;   ///< Time derivative of wall temperature; size Ngrid + 1.

  /**
   * \brief Returns the canonical ODE state size.
   */
  size_t stateSize() const noexcept;

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
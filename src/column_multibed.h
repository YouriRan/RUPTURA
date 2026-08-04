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
#include "reaction.h"
#include "utils.h"

/**
 * \brief Centralized view of the flattened multibed column ODE state layout.
 *
 * The layout matches ColumnStateLayout: gas concentration, physisorption,
 * chemisorption, optional surface/pore concentration, and temperature fields.
 */
struct MultibedColumnStateLayout
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

  [[nodiscard]] std::span<double> concentration(double* base) const noexcept { return componentBlock(base, 0); }
  [[nodiscard]] std::span<double> physisorption(double* base) const noexcept { return componentBlock(base, 1); }
  [[nodiscard]] std::span<double> chemisorption(double* base) const noexcept
  {
    return {base + 2 * componentBlockSize(), chemisorptionSiteCount() * componentBlockSize()};
  }
  [[nodiscard]] std::span<double> surfaceConcentration(double* base) const noexcept
  {
    return includeSurfacePoreTransport ? std::span<double>{base + (2 + chemisorptionSiteCount()) * componentBlockSize(),
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
struct MultibedColumn
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
  MultibedColumn(std::vector<MixturePrediction> physisorptionMixtures, std::vector<Component> components,
                 BoundaryCondition boundaryCondition, bool energyBalance, size_t numberOfGridPoints,
                 size_t maxIsothermTerms, size_t carrierGasComponent, double temperature, double inletPressure,
                 double outletPressure, double pressureGradient, std::vector<double> adsorbentVoidFractions,
                 std::vector<double> particleDensities, double columnEntranceVelocity,
                 std::vector<double> adsorbentLengths, std::vector<double> adsorbentInterfaceLengths,
                 std::vector<size_t> adsorbentGridPoints, double dynamicViscosity,
                 std::vector<double> particleDiameters, double influxTemperature, double internalDiameter,
                 double outerDiameter, double wallDensity, double gasThermalConductivity,
                 double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
                 double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
                 double heatCapacityWall, std::vector<double> columnDistances = {},
                 std::vector<Reaction> reactions = {})
      : physisorptionMixtures(std::move(physisorptionMixtures)),
        components(std::move(components)),
        boundaryCondition(boundaryCondition),
        energyBalance(energyBalance),
        reactions(std::move(reactions)),
        numberOfGridPoints(numberOfGridPoints),
        numberOfComponents(this->components.size()),
        numberOfAdsorbents(this->physisorptionMixtures.size()),
        maxIsothermTerms(maxIsothermTerms),
        maxChemisorptionSites(maximumChemisorptionSites(this->physisorptionMixtures)),
        numberOfCalls(0),
        carrierGasComponent(carrierGasComponent),
        adsorbentLengths(std::move(adsorbentLengths)),
        adsorbentInterfaceLengths(std::move(adsorbentInterfaceLengths)),
        adsorbentGridPoints(std::move(adsorbentGridPoints)),
        adsorbentVoidFractions(std::move(adsorbentVoidFractions)),
        particleDensities(std::move(particleDensities)),
        particleDiameters(std::move(particleDiameters)),
        externalTemperature(temperature),
        inletPressure(inletPressure),
        outletPressure(outletPressure),
        pressureGradient(pressureGradient),
        columnEntranceVelocity(columnEntranceVelocity),
        columnLength(std::reduce(this->adsorbentLengths.begin(), this->adsorbentLengths.end(), 0.0)),
        dynamicViscosity(dynamicViscosity),
        columnDistances(columnDistances.empty()
                            ? makeUniformColumnDistances(this->numberOfGridPoints, this->columnLength)
                            : std::move(columnDistances)),
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
        surfacePoreTransportEnabled(requiresSurfacePoreTransport(this->physisorptionMixtures) ||
                                    reactionsRequirePoreConcentration(this->reactions)),
        prefactorMassTransfer(this->numberOfComponents),
        idealGasMolFractions(this->numberOfComponents),
        adsorbedMolFractions(this->numberOfComponents),
        numberOfMolecules(this->numberOfComponents),
        interstitialGasVelocity(this->numberOfGridPoints + 1),
        gasDensity(this->numberOfGridPoints + 1),
        totalConcentration(this->numberOfGridPoints + 1),
        totalPressure(this->numberOfGridPoints + 1),
        totalVoidFraction(this->numberOfGridPoints + 1),
        particleDensity(this->numberOfGridPoints + 1),
        moleFraction((this->numberOfGridPoints + 1) * this->numberOfComponents),
        partialPressure((this->numberOfGridPoints + 1) * this->numberOfComponents),
        equilibriumPhysisorption((this->numberOfGridPoints + 1) * this->numberOfComponents),
        equilibriumChemisorption(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                 this->numberOfComponents),
        fractionOfAdsorbent((this->numberOfGridPoints + 1) * this->numberOfAdsorbents),
        hasAdsorbentOfType((this->numberOfGridPoints + 1) * this->numberOfAdsorbents),
        adsorbentScaledVoidFraction((this->numberOfGridPoints + 1) * this->numberOfAdsorbents),
        cachedPressure((this->numberOfGridPoints + 1) * this->numberOfAdsorbents * this->numberOfComponents *
                       this->maxIsothermTerms),
        cachedGrandPotential((this->numberOfGridPoints + 1) * this->numberOfAdsorbents * this->maxIsothermTerms),
        cachedChemisorptionPressure(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                    this->numberOfAdsorbents * this->numberOfComponents),
        cachedChemisorptionGrandPotential(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                          this->numberOfAdsorbents),
        coeffDiffusion(this->numberOfGridPoints + 1),
        facePressures(this->numberOfGridPoints),
        massFlux((this->numberOfGridPoints + 1) * this->numberOfComponents),
        bulkSpeciesSink((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionPhysisorptionSource((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionChemisorptionSource(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                    this->numberOfComponents),
        reactionPoreConcentrationSource(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                        this->numberOfComponents),
        reactionHeat(this->numberOfGridPoints + 1),
        state(MultibedColumnStateLayout{this->numberOfGridPoints, this->numberOfComponents, this->maxChemisorptionSites,
                                        this->surfacePoreTransportEnabled}
                  .stateSize(),
              0.0),
        stateDot(MultibedColumnStateLayout{this->numberOfGridPoints, this->numberOfComponents,
                                           this->maxChemisorptionSites, this->surfacePoreTransportEnabled}
                     .stateSize(),
                 0.0)
  {
    validateColumnDistances(this->columnDistances, this->numberOfGridPoints, this->columnLength);
    bindStateViews();
    chemisorptionMixtures = makeChemisorptionMixtures(this->physisorptionMixtures);
  }

  /**
   * \brief Constructs a column from parsed input data.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  MultibedColumn(const InputReader& inputReader)
      : MultibedColumn(
            [&inputReader]()
            {
              std::vector<MixturePrediction> mixtures;
              mixtures.reserve(inputReader.adsorbentComponents.size());
              for (const std::vector<Component>& componentsForAdsorbent : inputReader.adsorbentComponents)
              {
                mixtures.emplace_back(inputReader.displayName, componentsForAdsorbent, inputReader.numberOfCarrierGases,
                                      inputReader.carrierGasComponent, inputReader.temperature,
                                      inputReader.pressureStart, inputReader.pressureEnd,
                                      inputReader.numberOfPressurePoints, inputReader.pressureScale,
                                      inputReader.mixturePredictionMethod, inputReader.IASTMethod);
              }
              return mixtures;
            }(),
            inputReader.components, BoundaryCondition(inputReader.boundaryCondition), inputReader.energyBalance,
            inputReader.numberOfGridPoints, inputReader.maxIsothermTerms, inputReader.carrierGasComponent,
            inputReader.temperature, inputReader.inletPressure, inputReader.outletPressure,
            inputReader.pressureGradient, inputReader.adsorbentVoidFractions, inputReader.adsorbentParticleDensities,
            inputReader.columnEntranceVelocity, inputReader.adsorbentLengths, inputReader.adsorbentInterfaceLengths,
            inputReader.adsorbentGridPoints, inputReader.dynamicViscosity, inputReader.adsorbentParticleDiameters,
            inputReader.influxTemperature, inputReader.internalDiameter, inputReader.outerDiameter,
            inputReader.wallDensity, inputReader.gasThermalConductivity, inputReader.wallThermalConductivity,
            inputReader.heatTransferGasSolid, inputReader.heatTransferGasWall, inputReader.heatTransferWallExternal,
            inputReader.heatCapacityGas, inputReader.heatCapacitySolid, inputReader.heatCapacityWall,
            inputReader.columnDistances, inputReader.reactions)
  {
  }
  /**
   * \brief Copy constructor; rebinds spans to this object's state storage.
   */
  MultibedColumn(const MultibedColumn& other);

  /**
   * \brief Copy assignment; rebinds spans to this object's state storage.
   */
  MultibedColumn& operator=(const MultibedColumn& other);

  // Model configuration and component data.
  std::vector<MixturePrediction> physisorptionMixtures;  ///< One physisorption mixture-prediction object per bed.
  std::vector<MixturePrediction> chemisorptionMixtures;  ///< One competitive chemisorption model per bed.
  std::vector<Component> components;                     ///< Feed component definitions; size numberOfComponents.
  BoundaryCondition boundaryCondition;                   ///< Selected breakthrough boundary-condition pair.
  bool energyBalance;                                    ///< Enables gas/solid/wall temperature dynamics when true.
  std::vector<Reaction> reactions;                       ///< Reactions coupled to adsorbed or pore-phase variables.

  // Dimensions and counters.
  size_t numberOfGridPoints;     ///< Number of spatial grid intervals; node count is numberOfGridPoints + 1.
  size_t numberOfComponents;     ///< Number of gas components.
  size_t numberOfAdsorbents;     ///< Number of adsorbents.
  size_t maxIsothermTerms;       ///< Maximum number of isotherm sites across all components.
  size_t maxChemisorptionSites;  ///< Maximum chemisorption-site count across all beds and components.
  size_t numberOfCalls;          ///< Counter for model/evaluation calls.
  size_t carrierGasComponent;    ///< Index of the carrier-gas component.

  // Size numberOfAdsorbents. Indexed as value[ads].
  std::vector<double> adsorbentLengths;           ///< Pure adsorbent-region lengths in m.
  std::vector<double> adsorbentInterfaceLengths;  ///< Linear interface lengths between adsorbents in m.
  std::vector<size_t> adsorbentGridPoints;        ///< Spatial grid intervals per adsorbent section.
  std::vector<double> adsorbentVoidFractions;     ///< Packed-bed void fraction for each adsorbent.
  std::vector<double> particleDensities;          ///< Particle density for each adsorbent in kg/m^3.
  std::vector<double> particleDiameters;          ///< Particle diameter for each adsorbent in m.

  // Column operating conditions and geometry.
  double externalTemperature;           ///< External/reference gas temperature, T, in K.
  double inletPressure;                 ///< Inlet pressure, P_in, in Pa.
  double outletPressure;                ///< Outlet pressure, P_out, in Pa.
  double pressureGradient;              ///< Input pressure-gradient parameter.
  double columnEntranceVelocity;        ///< Inlet/interstitial velocity, v_in, in m/s.
  double columnLength;                  ///< Column length, L, in m.
  double dynamicViscosity;              ///< Gas dynamic viscosity used in Ergun calculations.
  std::vector<double> columnDistances;  ///< Spatial grid node positions in m.

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
  double resolution;                 ///< Spatial grid spacing, dz = L / numberOfGridPoints.
  double timeNormalizationFactor;    ///< Dimensionless-time factor, v_in / L.
  bool surfacePoreTransportEnabled;  ///< Adds surface/pore state for mechanistic or pore reactions.

  std::pair<size_t, size_t> iastPerformance{0, 0};  ///< Accumulated mixture-prediction diagnostics.

  // Size numberOfComponents. Indexed as value[comp].
  std::vector<double> prefactorMassTransfer;  ///< Per-component mass-transfer prefactor.
  std::vector<double> idealGasMolFractions;   ///< Temporary gas-phase mole fractions.
  std::vector<double> adsorbedMolFractions;   ///< Temporary adsorbed-phase mole fractions.
  std::vector<double> numberOfMolecules;      ///< Temporary equilibrium loading result per component.

  // Size numberOfGridPoints + 1. Indexed as value[grid].
  std::vector<double> interstitialGasVelocity;  ///< Interstitial gas velocity at each grid node.
  std::vector<double> gasDensity;               ///< Gas density at each grid node.
  std::vector<double> totalConcentration;       ///< Total gas concentration at each grid node.
  std::vector<double> totalPressure;            ///< Total pressure at each grid node.
  std::vector<double> totalVoidFraction;        ///< Packed-bed void fraction, epsilon.
  std::vector<double> particleDensity;          ///< Particle density in kg/m^3.

  // Size (numberOfGridPoints + 1) * numberOfComponents. Grid-major index: grid * numberOfComponents + comp.
  std::vector<double> moleFraction;              ///< Derived gas-phase mole fraction y_i.
  std::vector<double> partialPressure;           ///< Component partial pressure at each grid node.
  std::vector<double> equilibriumPhysisorption;  ///< Component equilibrium loading at each grid node.
  std::vector<double> equilibriumChemisorption;  ///< Site-major equilibrium chemisorption loading.

  // Size (numberOfGridPoints + 1) * numberOfAdsorbents. Grid-major index: grid * numberOfAdsorbents + ads.
  std::vector<double> fractionOfAdsorbent;  ///< Fraction of the column that has this adsorbent
  std::vector<bool> hasAdsorbentOfType;     ///< Boolean switch to determine if this adsorbent is here.
  std::vector<double> adsorbentScaledVoidFraction;

  // Mixture-prediction cache arrays.
  // cachedPressure: (numberOfGridPoints + 1) * numberOfAdsorbents * numberOfComponents * maxIsothermTerms
  // Grid-major index: (grid * numberOfAdsorbents + ads) * numberOfComponents * maxIsothermTerms
  std::vector<double> cachedPressure;  ///< Cached hypothetical pressures for mixture prediction.
  // cachedGrandPotential: (numberOfGridPoints + 1) * numberOfAdsorbents * maxIsothermTerms.
  // Grid-major index: (grid * numberOfAdsorbents + ads) * maxIsothermTerms
  std::vector<double> cachedGrandPotential;               ///< Cached reduced grand potentials for mixture prediction.
  std::vector<double> cachedChemisorptionPressure;        ///< Per-bed competitive chemisorption pressure cache.
  std::vector<double> cachedChemisorptionGrandPotential;  ///< Per-bed chemisorption spreading-pressure cache.

  // Scratch/work arrays.
  // coeffDiffusion is size numberOfGridPoints + 1; facePressures is size numberOfGridPoints; massFlux is grid-major
  // size (numberOfGridPoints + 1) * numberOfComponents.
  std::vector<double> coeffDiffusion;                   ///< Temporary diffusion coefficient field.
  std::vector<double> facePressures;                    ///< Pressure values at cell faces.
  std::vector<double> massFlux;                         ///< Component mass flux at each grid node.
  std::vector<double> bulkSpeciesSink;                  ///< Solid-uptake sink in the bulk gas balance.
  std::vector<double> reactionPhysisorptionSource;      ///< Reaction-only physisorbed source.
  std::vector<double> reactionChemisorptionSource;      ///< Reaction-only site-major chemisorbed source.
  std::vector<double> reactionPoreConcentrationSource;  ///< Reaction-only site-major pore source.
  std::vector<double> reactionHeat;                     ///< Reaction heat release on a solid-mass basis.

  // Canonical ODE storage.
  // Layout: concentration, physisorption, chemisorption, optional surface/pore concentration, gas T, solid T, wall T.
  std::vector<double> state;     ///< Canonical ODE state vector.
  std::vector<double> stateDot;  ///< Time derivative of the canonical ODE state vector.

  // Views into state/stateDot. Non-owning; rebind after copy/assignment.
  // Size (numberOfGridPoints + 1) * numberOfComponents. Grid-major index: grid * numberOfComponents + comp.
  std::span<double>
      concentration;  ///< Gas concentration c_i in mol/m^3; size (numberOfGridPoints + 1) * numberOfComponents.
  std::span<double> concentrationDot;  ///< Time derivative dc_i/dt; size (numberOfGridPoints + 1) * numberOfComponents.
  std::span<double> physisorption;     ///< Adsorbed loading q_i; size (numberOfGridPoints + 1) * numberOfComponents.
  std::span<double> physisorptionDot;  ///< Time derivative dq_i/dt; size (numberOfGridPoints + 1) * numberOfComponents.
  std::span<double> chemisorption;     ///< Site-major chemisorbed loading.
  std::span<double> chemisorptionDot;  ///< Site-major chemisorbed loading derivative.
  std::span<double> surfaceConcentration;     ///< Site-major particle-surface concentration.
  std::span<double> surfaceConcentrationDot;  ///< Surface-concentration derivative.
  std::span<double> poreConcentration;        ///< Site-major pore concentration.
  std::span<double> poreConcentrationDot;     ///< Pore-concentration derivative.
  std::span<double> gasTemperature;           ///< Gas temperature; size numberOfGridPoints + 1.
  std::span<double> gasTemperatureDot;        ///< Time derivative of gas temperature; size numberOfGridPoints + 1.
  std::span<double> solidTemperature;         ///< Solid temperature; size numberOfGridPoints + 1.
  std::span<double> solidTemperatureDot;      ///< Time derivative of solid temperature; size numberOfGridPoints + 1.
  std::span<double> wallTemperature;          ///< Wall temperature; size numberOfGridPoints + 1.
  std::span<double> wallTemperatureDot;       ///< Time derivative of wall temperature; size numberOfGridPoints + 1.

  /**
   * \brief Returns the canonical ODE state size.
   */
  size_t stateSize() const noexcept;

  static bool requiresSurfacePoreTransport(const std::vector<MixturePrediction>& mixtures) noexcept;
  static size_t maximumChemisorptionSites(const std::vector<MixturePrediction>& mixtures) noexcept;
  static std::vector<MixturePrediction> makeChemisorptionMixtures(
      const std::vector<MixturePrediction>& physisorptionMixtures);

  /**
   * \brief Returns the canonical multibed ODE state layout.
   */
  [[nodiscard]] MultibedColumnStateLayout stateLayout() const noexcept;

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

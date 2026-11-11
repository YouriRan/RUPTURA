#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
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
 * \brief Non-owning views of the state fields and their time derivatives in a state vector and a derivative vector
 *        with the ColumnStateLayout layout: the column's own state, or a trial state of CVODE.
 */
struct ColumnStateViews
{
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
};

/**
 * \brief Centralized view of the flattened column ODE state layout: gas concentration, physisorption,
 * chemisorption, optional surface/pore concentration, and temperature fields.
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

  /**
   * \brief Views of the fields of the state vector at base and of their derivatives in the vector at baseDot.
   */
  [[nodiscard]] ColumnStateViews views(double* base, double* baseDot) const noexcept
  {
    return {.concentration = concentration(base),
            .concentrationDot = concentration(baseDot),
            .physisorption = physisorption(base),
            .physisorptionDot = physisorption(baseDot),
            .chemisorption = chemisorption(base),
            .chemisorptionDot = chemisorption(baseDot),
            .surfaceConcentration = surfaceConcentration(base),
            .surfaceConcentrationDot = surfaceConcentration(baseDot),
            .poreConcentration = poreConcentration(base),
            .poreConcentrationDot = poreConcentration(baseDot),
            .gasTemperature = gasTemperature(base),
            .gasTemperatureDot = gasTemperature(baseDot),
            .solidTemperature = solidTemperature(base),
            .solidTemperatureDot = solidTemperature(baseDot),
            .wallTemperature = wallTemperature(base),
            .wallTemperatureDot = wallTemperature(baseDot)};
  }
};

/**
 * \brief Owned data of a column: parameters, component data, cache arrays, scratch arrays, and the canonical ODE state
 * storage.
 *
 * Column adds the non-owning views into the state and the methods. Keeping every owned member here lets Column copy
 * them with the implicit copy operations of this base and only rebind its views.
 */
struct ColumnData
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

  enum class FluidPhase
  {
    Gas = 0,
    Liquid = 1
  };

  enum class PHMode
  {
    Fixed = 0,
    HPlus = 1,
    OHMinus = 2
  };

  /**
   * \brief Reconstruction of the convective face values (Supplementary "Spatial discretization").
   */
  enum class ConvectionScheme
  {
    Upwind = 0,  ///< First-order upwind: the face value is the upstream node value.
    WENO3 = 1    ///< Third-order WENO; gas phase: mole fractions at the faces, total molar flux upwind.
  };

  /**
   * \brief Species and fluid-temperature condition at the inlet (Supplementary "Boundary conditions").
   */
  enum class InletCondition
  {
    Dirichlet = 0,  ///< The inlet face holds the feed state; axial dispersion adds to the feed flux.
    Danckwerts = 1  ///< The total convective-dispersive flux into the first cell equals the feed flux.
  };

  /**
   * \brief End of the column at which the feed enters (Supplementary "Boundary conditions").
   */
  enum class FlowDirection
  {
    Forward = 0,  ///< The feed enters at z = 0 from node 0 and leaves at z = L.
    Reverse = 1   ///< The feed enters at z = L from the product-end gas and leaves at z = 0.
  };

  /**
   * \brief The flow direction and the gas at the product end z = L, as the transport and closure kernels take them.
   *
   * In a forward step the feed enters at z = 0 from node 0, and gas that flows back in at z = L is the product-end
   * gas, or that of the last cell when none is set. In a reverse step the product-end gas is the feed, the boundary
   * settings apply at z = L for the inlet and at z = 0 for the outlet, and gas that flows back in at z = 0 is that of
   * the first cell; node 0 is not used. A closed end is the end opposite the one where the feed enters: z = L in a
   * forward step, z = 0 in a reverse step. Gas then enters or leaves at the other end only, as the pressure rate and
   * the uptake require.
   */
  struct FlowEnds
  {
    FlowDirection direction{FlowDirection::Forward};
    std::span<const double> productEndConcentration{};  ///< In mol/m^3; empty: the gas of the last cell.
    double productEndTemperature{0.0};                  ///< In K; 0: that of the last cell.
    bool closedEnd{false};                              ///< The velocity is zero at the end opposite the feed.
    double pressureRate{0.0};                           ///< dP/dt of the column pressure in Pa/s.
    double relativePressureRate{0.0};                   ///< (dP/dt) / P in 1/s.
  };

  /**
   * \brief Inlet pressure of a swing phase: from startPressure to targetPressure, exponentially with the time constant,
   *        or linearly over the time constant and constant after it.
   */
  struct PressureRamp
  {
    enum class Shape
    {
      None = 0,  ///< No ramp: the pressure is that of the boundary settings.
      Exponential = 1,
      Linear = 2
    };
    Shape shape{Shape::None};
    double timeConstant{0.0};    ///< In s.
    double startTime{0.0};       ///< In s.
    double startPressure{0.0};   ///< In Pa.
    double targetPressure{0.0};  ///< In Pa.
    double time{0.0};            ///< Time of the current pressure in s.
  };

  /**
   * \brief Constructs a column from explicit model/configuration arguments.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  ColumnData(std::vector<MixturePrediction> physisorptionMixtures, std::vector<Component> components,
                 BoundaryCondition boundaryCondition, bool energyBalance, size_t numberOfGridPoints,
                 size_t maxIsothermTerms, size_t carrierGasComponent, double temperature, double inletPressure,
                 double outletPressure, double pressureGradient, std::vector<double> adsorbentVoidFractions,
                 std::vector<double> particleDensities, double columnEntranceVelocity,
                 std::vector<double> adsorbentLengths, std::vector<double> adsorbentInterfaceLengths,
                 std::vector<size_t> adsorbentInterfaceGridPoints, double dynamicViscosity, std::vector<double> particleDiameters, double influxTemperature,
                 double internalDiameter, double outerDiameter, double wallDensity, double gasThermalConductivity,
                 double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
                 double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
                 double heatCapacityWall, std::vector<double> columnDistances = {},
                 std::vector<Reaction> reactions = {}, std::vector<Geometry> geometries = {},
                 std::vector<double> adsorbentMixFractions = {}, double configuredColumnLength = -1.0,
                 FluidPhase fluidPhase = FluidPhase::Gas, double liquidDensity = 1000.0, PHMode pHMode = PHMode::Fixed,
                 double pHValue = 7.0, double pKw = 14.0, size_t pHComponent = 0)
      : physisorptionMixtures(std::move(physisorptionMixtures)),
        components(std::move(components)),
        boundaryCondition(boundaryCondition),
        fluidPhase(fluidPhase),
        pHMode(pHMode),
        energyBalance(energyBalance),
        reactions(std::move(reactions)),
        numberOfGridPoints(numberOfGridPoints),
        numberOfComponents(this->components.size()),
        numberOfAdsorbents(this->physisorptionMixtures.size()),
        maxIsothermTerms(maxIsothermTerms),
        maxChemisorptionSites(maximumChemisorptionSites(this->physisorptionMixtures)),
        carrierGasComponent(carrierGasComponent),
        pHComponent(pHComponent),
        adsorbentLengths(std::move(adsorbentLengths)),
        adsorbentInterfaceLengths(std::move(adsorbentInterfaceLengths)),
        adsorbentMixFractions(std::move(adsorbentMixFractions)),
        adsorbentInterfaceGridPoints(std::move(adsorbentInterfaceGridPoints)),
        adsorbentVoidFractions(std::move(adsorbentVoidFractions)),
        particleDensities(std::move(particleDensities)),
        particleDiameters(std::move(particleDiameters)),
        geometries(geometries.empty() ? makePackedBedGeometries(this->adsorbentVoidFractions, this->particleDiameters,
                                                                internalDiameter, outerDiameter)
                                      : std::move(geometries)),
        externalTemperature(temperature),
        inletPressure(inletPressure),
        outletPressure(outletPressure),
        pressureGradient(pressureGradient),
        columnEntranceVelocity(columnEntranceVelocity),
        columnLength(
            !this->adsorbentMixFractions.empty() && configuredColumnLength > 0.0
                ? configuredColumnLength
                : std::reduce(this->adsorbentLengths.begin(), this->adsorbentLengths.end(), 0.0) +
                      std::reduce(this->adsorbentInterfaceLengths.begin(), this->adsorbentInterfaceLengths.end(), 0.0)),
        dynamicViscosity(dynamicViscosity),
        liquidDensity(liquidDensity),
        pHValue(pHValue),
        pKw(pKw),
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
        timeNormalizationFactor(this->columnEntranceVelocity / this->columnLength),
        surfacePoreTransportEnabled(requiresSurfacePoreTransport(this->physisorptionMixtures) ||
                                    reactionsRequirePoreConcentration(this->reactions)),
        idealGasMolFractions(this->numberOfComponents),
        adsorbedMolFractions(this->numberOfComponents),
        numberOfMolecules(this->numberOfComponents),
        interstitialGasVelocity(this->numberOfGridPoints + 1),
        gasDensity(this->numberOfGridPoints + 1),
        totalConcentration(this->numberOfGridPoints + 1),
        totalPressure(this->numberOfGridPoints + 1),
        pH(this->numberOfGridPoints + 1, pHValue),
        totalVoidFraction(this->numberOfGridPoints + 1),
        particleDensity(this->numberOfGridPoints + 1),
        moleFraction((this->numberOfGridPoints + 1) * this->numberOfComponents),
        partialPressure((this->numberOfGridPoints + 1) * this->numberOfComponents),
        equilibriumPhysisorption((this->numberOfGridPoints + 1) * this->numberOfComponents),
        equilibriumChemisorption(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                 this->numberOfComponents),
        fractionOfAdsorbent((this->numberOfGridPoints + 1) * this->numberOfAdsorbents),
        hasAdsorbentOfType((this->numberOfGridPoints + 1) * this->numberOfAdsorbents),
        loadingPrefactor(this->numberOfGridPoints + 1),
        poreLoadingConversion(this->numberOfGridPoints + 1),
        solidFluidContactArea(this->numberOfGridPoints + 1),
        fluidSolidContactArea(this->numberOfGridPoints + 1),
        pressureDropLaw(this->numberOfGridPoints + 1),
        dispersionLaw(this->numberOfGridPoints + 1),
        massTransferCoefficient((this->numberOfGridPoints + 1) * this->numberOfComponents),
        molecularDispersionCoefficient((this->numberOfGridPoints + 1) * this->numberOfComponents),
        axialDispersionCoefficient((this->numberOfGridPoints + 1) * this->numberOfComponents),
        heatOfAdsorption((this->numberOfGridPoints + 1) * this->numberOfComponents),
        heatOfChemisorption(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) * this->numberOfComponents),
        cachedPressure((this->numberOfGridPoints + 1) * this->numberOfAdsorbents * this->numberOfComponents *
                       this->maxIsothermTerms),
        cachedGrandPotential((this->numberOfGridPoints + 1) * this->numberOfAdsorbents * this->maxIsothermTerms),
        cachedChemisorptionPressure(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                    this->numberOfAdsorbents * this->numberOfComponents),
        cachedChemisorptionGrandPotential(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                          this->numberOfAdsorbents),
        coeffDiffusion(this->numberOfGridPoints + 1),
        bulkSpeciesSink((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionPhysisorptionSource((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionChemisorptionSource(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                    this->numberOfComponents),
        reactionPoreConcentrationSource(this->maxChemisorptionSites * (this->numberOfGridPoints + 1) *
                                        this->numberOfComponents),
        reactionBulkSource((this->numberOfGridPoints + 1) * this->numberOfComponents),
        reactionHeat(this->numberOfGridPoints + 1),
        reactionCoverage(this->numberOfComponents),
        wenoScratch((this->numberOfGridPoints + 1) * (this->numberOfComponents + 4)),
        state(ColumnStateLayout{this->numberOfGridPoints, this->numberOfComponents, this->maxChemisorptionSites,
                                        this->surfacePoreTransportEnabled}
                  .stateSize(),
              0.0),
        stateDot(ColumnStateLayout{this->numberOfGridPoints, this->numberOfComponents,
                                           this->maxChemisorptionSites, this->surfacePoreTransportEnabled}
                     .stateSize(),
                 0.0)
  {
    validateColumnDistances(this->columnDistances, this->numberOfGridPoints, this->columnLength);
    chemisorptionMixtures = makeChemisorptionMixtures(this->physisorptionMixtures);
  }

  /**
   * \brief Constructs a column from parsed input data.
   *
   * Allocates grid, component, cache, scratch, and ODE-state arrays.
   */
  ColumnData(const InputReader& inputReader)
      : ColumnData(
            [&inputReader]()
            {
              std::vector<MixturePrediction> mixtures;
              mixtures.reserve(inputReader.adsorbentComponents.size());
              for (const std::vector<Component>& componentsForAdsorbent : inputReader.adsorbentComponents)
              {
                // MPD loads its particle distribution from MPDSettings, which describes the single adsorbent.
                if (inputReader.mixturePredictionMethod ==
                    static_cast<size_t>(MixturePrediction::PredictionMethod::MPD))
                {
                  mixtures.emplace_back(inputReader);
                  continue;
                }
                mixtures.emplace_back(inputReader.displayName, componentsForAdsorbent, inputReader.numberOfCarrierGases,
                                      inputReader.carrierGasComponent, inputReader.temperature,
                                      inputReader.pressureStart, inputReader.pressureEnd,
                                      inputReader.numberOfPressurePoints, inputReader.pressureScale,
                                      inputReader.mixturePredictionMethod, inputReader.IASTMethod);
                mixtures.back().relativeHumidity = inputReader.relativeHumidity;
                mixtures.back().waterSaturationPressure = inputReader.waterSaturationPressure;
              }
              return mixtures;
            }(),
            inputReader.components, BoundaryCondition(inputReader.boundaryCondition), inputReader.energyBalance,
            inputReader.numberOfGridPoints, inputReader.maxIsothermTerms, inputReader.carrierGasComponent,
            inputReader.temperature, inputReader.inletPressure, inputReader.outletPressure,
            inputReader.pressureGradient, inputReader.adsorbentVoidFractions, inputReader.adsorbentParticleDensities,
            inputReader.columnEntranceVelocity, inputReader.adsorbentLengths, inputReader.adsorbentInterfaceLengths,
            inputReader.adsorbentInterfaceGridPoints, inputReader.dynamicViscosity,
            inputReader.adsorbentParticleDiameters, inputReader.influxTemperature, inputReader.internalDiameter,
            inputReader.outerDiameter, inputReader.wallDensity,
            inputReader.fluidPhase == 1 ? inputReader.liquidThermalConductivity : inputReader.gasThermalConductivity,
            inputReader.wallThermalConductivity,
            inputReader.fluidPhase == 1 ? inputReader.heatTransferLiquidSolid : inputReader.heatTransferGasSolid,
            inputReader.fluidPhase == 1 ? inputReader.heatTransferLiquidWall : inputReader.heatTransferGasWall,
            inputReader.heatTransferWallExternal,
            inputReader.fluidPhase == 1 ? inputReader.heatCapacityLiquid : inputReader.heatCapacityGas,
            inputReader.heatCapacitySolid, inputReader.heatCapacityWall, inputReader.columnDistances,
            inputReader.reactions, inputReader.adsorbentGeometries, inputReader.adsorbentMixFractions,
            inputReader.columnLength, FluidPhase(inputReader.fluidPhase), inputReader.liquidDensity,
            PHMode(inputReader.pHMode), inputReader.pHValue, inputReader.pKw, inputReader.pHComponent)
  {
    convectionScheme = ConvectionScheme(inputReader.convectionScheme);
    inletCondition = InletCondition(inputReader.inletCondition);
    flowDirection = FlowDirection(inputReader.flowDirection);
    dispersionReferenceTemperature = inputReader.axialDispersionReferenceTemperature;
    dispersionReferencePressure = inputReader.axialDispersionReferencePressure;
  }

  // Model configuration and component data.
  std::vector<MixturePrediction> physisorptionMixtures;  ///< One physisorption mixture-prediction object per bed.
  std::vector<MixturePrediction> chemisorptionMixtures;  ///< One competitive chemisorption model per bed.
  std::vector<Component> components;                     ///< Feed component definitions; size numberOfComponents.
  BoundaryCondition boundaryCondition;                   ///< Selected breakthrough boundary-condition pair.
  FluidPhase fluidPhase;                                 ///< Gas or liquid mobile phase.
  PHMode pHMode;                                         ///< Fixed or transported pH representation.
  ConvectionScheme convectionScheme{ConvectionScheme::WENO3};  ///< Reconstruction of the convective face values.
  InletCondition inletCondition{InletCondition::Danckwerts};   ///< Species and fluid-temperature inlet condition.
  FlowDirection flowDirection{FlowDirection::Forward};         ///< End at which the feed enters; set per swing phase.
  /// Gas at the product end z = L in mol/m^3, per component: the feed of a reverse step, or the gas that flows back in
  /// at z = L in a forward step; empty: the gas of the last cell.
  std::vector<double> productEndConcentration;
  double productEndTemperature{0.0};  ///< Temperature of the product-end gas in K; 0: that of the last cell.
  bool closedEnd{false};             ///< The end opposite the feed end is closed; set per swing phase.
  PressureRamp pressureRamp{};       ///< Pressure ramp of the swing phase.
  double pressureRate{0.0};          ///< dP/dt of the inlet pressure in Pa/s; see advancePressureRamp.
  double relativePressureRate{0.0};  ///< (dP/dt) / P of the inlet pressure in 1/s; see advancePressureRamp.
  [[nodiscard]] FlowEnds flowEnds() const noexcept
  {
    return {flowDirection, productEndConcentration, productEndTemperature, closedEnd, pressureRate,
            relativePressureRate};
  }

  /**
   * \brief Sets the inlet pressure, its rates, and the product-end gas to their values at the given time; without a
   *        ramp it does nothing. The integrators call it before each right-hand side, which depends on the time through
   *        these fields only.
   */
  void advancePressureRamp(double time) noexcept
  {
    if (pressureRamp.shape == PressureRamp::Shape::None) return;
    pressureRamp.time = time;
    const double elapsed = time - pressureRamp.startTime;
    const double change = pressureRamp.targetPressure - pressureRamp.startPressure;
    double pressure = pressureRamp.targetPressure;
    pressureRate = 0.0;
    if (pressureRamp.shape == PressureRamp::Shape::Exponential)
    {
      const double decay = std::exp(-elapsed / pressureRamp.timeConstant);
      pressure -= change * decay;
      pressureRate = change * decay / pressureRamp.timeConstant;
    }
    else if (elapsed < pressureRamp.timeConstant)
    {
      pressure = pressureRamp.startPressure + change * elapsed / pressureRamp.timeConstant;
      pressureRate = change / pressureRamp.timeConstant;
    }
    relativePressureRate = pressureRate / pressure;
    for (double& value : productEndConcentration) value *= pressure / inletPressure;
    inletPressure = pressure;
  }
  /// Reference state of molecularDispersionCoefficient, which scales as (T/T_ref)^1.75 (P_ref/P); 0 = no scaling.
  double dispersionReferenceTemperature{0.0};
  double dispersionReferencePressure{0.0};  ///< See dispersionReferenceTemperature; in Pa.
  bool energyBalance;                                    ///< Enables gas/solid/wall temperature dynamics when true.
  std::vector<Reaction> reactions;                       ///< Reactions coupled to adsorbed or pore-phase variables.

  // Dimensions and counters.
  size_t numberOfGridPoints;     ///< Number of spatial grid intervals; node count is numberOfGridPoints + 1.
  size_t numberOfComponents;     ///< Number of gas components.
  size_t numberOfAdsorbents;     ///< Number of adsorbents.
  size_t maxIsothermTerms;       ///< Maximum number of isotherm sites across all components.
  size_t maxChemisorptionSites;  ///< Maximum chemisorption-site count across all beds and components.
  size_t carrierGasComponent;    ///< Index of the carrier-gas component.
  size_t pHComponent;            ///< Component used to derive transported pH.

  // Size numberOfAdsorbents. Indexed as value[ads].
  std::vector<double> adsorbentLengths;           ///< Pure adsorbent-region lengths in m; empty for a uniform mix.
  std::vector<double> adsorbentInterfaceLengths;  ///< Linear interface lengths between adsorbents in m.
  std::vector<double> adsorbentMixFractions;      ///< Uniform adsorbent volume fractions over the whole column.
  std::vector<size_t> adsorbentInterfaceGridPoints;  ///< Intervals per interface; zero for sharp interfaces.
  std::vector<double> adsorbentVoidFractions;     ///< Packed-bed void fraction for each adsorbent.
  std::vector<double> particleDensities;          ///< Particle density for each adsorbent in kg/m^3.
  std::vector<double> particleDiameters;          ///< Particle diameter for each adsorbent in m.
  std::vector<Geometry> geometries;                ///< Derived geometry for each adsorbent section.

  // Column operating conditions and geometry.
  double externalTemperature;           ///< External/reference gas temperature, T, in K.
  double inletPressure;                 ///< Inlet pressure, P_in, in Pa.
  double outletPressure;                ///< Outlet pressure, P_out, in Pa.
  double pressureGradient;              ///< Signed axial pressure gradient dP/dz in Pa/m.
  double columnEntranceVelocity;        ///< Inlet/interstitial velocity, v_in, in m/s.
  double columnLength;                  ///< Column length, L, in m.
  double dynamicViscosity;              ///< Gas dynamic viscosity used in Ergun calculations.
  double liquidDensity;                 ///< Constant liquid density in kg/m^3.
  double pHValue;                       ///< Fixed pH value when pHMode is Fixed.
  double pKw;                           ///< Water ion-product exponent used for OH- transport.
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
  double timeNormalizationFactor;    ///< Dimensionless-time factor, v_in / L.
  bool surfacePoreTransportEnabled;  ///< Adds surface/pore state for mechanistic or pore reactions.

  std::pair<size_t, size_t> iastPerformance{0, 0};  ///< Accumulated mixture-prediction diagnostics.

  // Size numberOfComponents. Indexed as value[comp].
  std::vector<double> idealGasMolFractions;   ///< Temporary gas-phase mole fractions.
  std::vector<double> adsorbedMolFractions;   ///< Temporary adsorbed-phase mole fractions.
  std::vector<double> numberOfMolecules;      ///< Temporary equilibrium loading result per component.

  // Size numberOfGridPoints + 1. Indexed as value[grid].
  std::vector<double> interstitialGasVelocity;  ///< Interstitial gas velocity at each grid node.
  std::vector<double> gasDensity;               ///< Gas density at each grid node.
  std::vector<double> totalConcentration;       ///< Total gas concentration at each grid node.
  std::vector<double> totalPressure;            ///< Total pressure at each grid node.
  std::vector<double> pH;                       ///< Local liquid pH; fixed at 7 for gas simulations.
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

  // Adsorbent blends, fixed after initialize() (computeAdsorbentBlends). Size numberOfGridPoints + 1.
  std::vector<double> loadingPrefactor;  ///< Bulk-sink factor sum_a f_a (V_s/V_f)_a eps_a rho_a / eps_total, kg/m^3.
  std::vector<double> poreLoadingConversion;     ///< Solid-volume-weighted particle density in kg/m^3.
  std::vector<double> solidFluidContactArea;     ///< Fluid-solid contact area per solid volume in 1/m.
  std::vector<double> fluidSolidContactArea;     ///< Fluid-solid contact area per fluid volume in 1/m.
  std::vector<PressureDropLaw> pressureDropLaw;  ///< Fraction-weighted pressure-drop law.
  std::vector<DispersionLaw> dispersionLaw;      ///< Fraction-weighted dispersivity and Taylor-Aris coefficient.
  // Size (numberOfGridPoints + 1) * numberOfComponents. Grid-major index: grid * numberOfComponents + comp.
  std::vector<double> massTransferCoefficient;  ///< LDF coefficient in 1/s.
  /// Velocity-independent part D_0 of the axial dispersion (AxialDispersionCoefficient) in m^2/s, at the reference
  /// state when one is set.
  std::vector<double> molecularDispersionCoefficient;
  /// Axial dispersion coefficient in m^2/s; set from the dispersion law in every updateVelocityAndPressure.
  std::vector<double> axialDispersionCoefficient;
  std::vector<double> heatOfAdsorption;            ///< Physisorption heat in J/mol.
  // Size maxChemisorptionSites * (numberOfGridPoints + 1) * numberOfComponents, site-major.
  std::vector<double> heatOfChemisorption;  ///< Chemisorption heat per site in J/mol.

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
  std::vector<double> coeffDiffusion;                   ///< Gas thermal diffusivity at each grid node.
  std::vector<double> bulkSpeciesSink;                  ///< Solid-uptake sink in the bulk gas balance.
  std::vector<double> reactionPhysisorptionSource;      ///< Reaction-only physisorbed source.
  std::vector<double> reactionChemisorptionSource;      ///< Reaction-only site-major chemisorbed source.
  std::vector<double> reactionPoreConcentrationSource;  ///< Reaction-only site-major pore source.
  std::vector<double> reactionBulkSource;  ///< Bulk-concentration reaction source per fluid volume, mol/(m^3 s).
  std::vector<double> reactionHeat;                     ///< Reaction heat release on a solid-mass basis.
  std::vector<double> reactionCoverage;  ///< Langmuir-Hinshelwood coverages of one node; size numberOfComponents.
  /// WENO3 node values, face values, face fluxes, total concentrations, and face-fraction sums of
  /// computeMassDerivatives; size (numberOfGridPoints + 1) * (numberOfComponents + 4).
  std::vector<double> wenoScratch;

  // Canonical ODE storage.
  // Layout: concentration, physisorption, chemisorption, optional surface/pore concentration, gas T, solid T, wall T.
  std::vector<double> state;     ///< Canonical ODE state vector.
  std::vector<double> stateDot;  ///< Time derivative of the canonical ODE state vector.

  static bool requiresSurfacePoreTransport(const std::vector<MixturePrediction>& mixtures) noexcept;
  static size_t maximumChemisorptionSites(const std::vector<MixturePrediction>& mixtures) noexcept;
  static std::vector<Geometry> makePackedBedGeometries(std::span<const double> voidFractions,
                                                        std::span<const double> particleDiameters,
                                                        double internalDiameter, double outerDiameter);
  static std::vector<MixturePrediction> makeChemisorptionMixtures(
      const std::vector<MixturePrediction>& physisorptionMixtures);
};

/**
 * \brief Packed-column model state and configuration.
 *
 * Holds one or more adsorbents, as layered sections (with optional linear interfaces) or as a uniform mixture; a
 * single adsorbent is the one-section case. The owned data are in ColumnData; Column adds the views into its own ODE
 * state (ColumnStateViews) and the methods.
 */
struct Column : ColumnData, ColumnStateViews
{
  /**
   * \brief Constructs the column data with a ColumnData constructor (explicit arguments or an InputReader) and binds
   * the state views.
   */
  template <typename... Arguments>
    requires std::constructible_from<ColumnData, Arguments...>
  explicit Column(Arguments&&... arguments) : ColumnData(std::forward<Arguments>(arguments)...)
  {
    bindStateViews();
  }

  /**
   * \brief Copy constructor; rebinds spans to this object's state storage.
   */
  Column(const Column& other) : ColumnData(other), ColumnStateViews{} { bindStateViews(); }

  /**
   * \brief Copy assignment; rebinds spans to this object's state storage.
   */
  Column& operator=(const Column& other)
  {
    ColumnData::operator=(other);
    bindStateViews();
    return *this;
  }

  /**
   * \brief Returns the canonical ODE state size.
   */
  size_t stateSize() const noexcept;

  /**
   * \brief Returns the canonical ODE state layout.
   */
  [[nodiscard]] ColumnStateLayout stateLayout() const noexcept;

  /**
   * \brief Rebinds all state/stateDot spans (the ColumnStateViews base) to the current vector storage; needed after
   *        copy/assignment.
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
   * \brief Fills the adsorbent blends (loadingPrefactor ... heatOfChemisorption) from fractionOfAdsorbent and the
   *        per-adsorbent properties; called by initialize() and readJSON(), and needed again after changing either.
   */
  void computeAdsorbentBlends();

  /**
   * \brief Largest |R T sum_j c_j / P - 1| over the nodes for the gas continuity closures, which keep sum_j c_j at
   *        P/RT; 0 for FixedVelocity, FixedPressureInletVelocity, and a liquid, which do not. A reverse step leaves out
   *        node 0, which it does not use.
   */
  [[nodiscard]] double idealGasDeviation() const;

  /**
   * \brief Amount of each component in the gas of the cells, per unit column cross-section, in mol/m^2.
   */
  [[nodiscard]] std::vector<double> gasInventory() const;

  /**
   * \brief Physisorbed and chemisorbed amount of each component in the cells, per unit column cross-section, in
   *        mol/m^2; without the surface and pore reservoirs of mechanistic chemisorption.
   */
  [[nodiscard]] std::vector<double> adsorbedInventory() const;

  /**
   * \brief Writes headers for component output files and the column output file.
   */
  void writeOutputHeader(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream) const;

  /**
   * \brief Writes one output snapshot at the given simulation time: per node, the values of appendColumnOutput in
   *        column.data and those of appendComponentOutput in the file of each component.
   */
  void writeOutput(std::vector<std::ofstream>& componentStreams, std::ofstream& columnStream, double time) const;

  /**
   * \brief The chemisorption sites of component comp with surface-pore transport, whose surface and pore
   *        concentrations the component output adds after its fixed columns.
   */
  [[nodiscard]] std::vector<size_t> surfacePoreTransportSites(size_t comp) const;

  /**
   * \brief Appends the column.data values of node grid at the given time to row (the columns of writeOutputHeader):
   *        tau, t in min, z, v, p_t, the gas, solid, and wall temperatures and their rates, the fluid density, and in
   *        the liquid phase the pH. The in-memory Python result (Breakthrough.compute) uses the same values.
   */
  void appendColumnOutput(std::vector<double>& row, size_t grid, double time) const;

  /**
   * \brief Appends the values of the output file of component comp at node grid to row: tau, t in min, z, c, dc/dt, y,
   *        q_phy, dq_phy/dt, q_chem, dq_chem/dt, p (c for a liquid), q_phy^*, the normalized driving force, q_chem^*,
   *        and the surface and pore concentrations of each site of surfacePoreTransportSites.
   */
  void appendComponentOutput(std::vector<double>& row, size_t comp, size_t grid, double time) const;

  /**
   * \brief Returns a compact text representation of key column settings.
   */
  std::string repr() const;

  /**
   * \brief Serializes the column state (concentrations, loadings, temperatures) to JSON.
   */
  void writeJSON(const std::string& fileName) const;

  /**
   * \brief Loads the column state from JSON, and the derived fields and caches that older files also hold.
   */
  void readJSON(const std::string& fileName);
};

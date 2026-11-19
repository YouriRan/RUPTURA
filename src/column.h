#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <fstream>
#include <functional>
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
  /// Temperature of the heat-transfer fluid along the outer wall in K, per node; empty without a Coolant.
  std::span<double> coolantTemperature;
  std::span<double> coolantTemperatureDot;
  /// Loading of the radial cells of the particle (ParticleTransport HSDM), from the centre outwards: one block of the
  /// size and order of physisorption per cell; empty with the linear driving force.
  std::span<double> shellLoading;
  std::span<double> shellLoadingDot;
  /// Net amount of each component that entered the column at z = 0, then at z = L, since the start of the swing phase,
  /// per unit column cross-section in mol/m^2 (negative: it left), and after these 2 numberOfComponents values the heat
  /// that entered through the outer wall and the work that brings the gas leaving through pressure ports to the
  /// exhaust pressure, in J/m^2; empty without them (ColumnStateLayout).
  std::span<double> streamAmount;
  std::span<double> streamAmountDot;  ///< The molar fluxes into the column at the ends in mol/(m^2 s), then W/m^2.
};

/**
 * \brief Residual of the stop criteria of a swing phase at a time in s and a column state: the integration of the
 *        phase ends at the time at which it turns from negative to non-negative.
 */
using StopResidual = std::function<double(double, const ColumnStateViews&)>;

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
  /// The state ends with the net amounts that entered through the two ends (ColumnStateViews::streamAmount), which
  /// the integrators then integrate with the column: for stop criteria, cyclic steady state, and MomentumBalance.
  bool includeStreamAmounts{false};
  bool includeCoolant{false};  ///< After the temperatures: that of the heat-transfer fluid (coolantTemperature).
  size_t numberOfShells{0};    ///< Then the loadings of that many radial cells of the particle (shellLoading).

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
  [[nodiscard]] size_t streamAmountSize() const noexcept
  {
    return includeStreamAmounts ? 2 * numberOfComponents + 2 : 0;
  }
  [[nodiscard]] size_t coolantSize() const noexcept { return includeCoolant ? scalarBlockSize() : 0; }
  [[nodiscard]] size_t shellOffset() const noexcept
  {
    return temperatureOffset() + 3 * scalarBlockSize() + coolantSize();
  }
  [[nodiscard]] size_t stateSize() const noexcept
  {
    return shellOffset() + numberOfShells * componentBlockSize() + streamAmountSize();
  }

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
  [[nodiscard]] std::span<double> coolantTemperature(double* base) const noexcept
  {
    return {base + temperatureOffset() + 3 * scalarBlockSize(), coolantSize()};
  }
  [[nodiscard]] std::span<double> shellLoading(double* base) const noexcept
  {
    return {base + shellOffset(), numberOfShells * componentBlockSize()};
  }
  [[nodiscard]] std::span<double> streamAmount(double* base) const noexcept
  {
    return {base + shellOffset() + numberOfShells * componentBlockSize(), streamAmountSize()};
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
            .wallTemperatureDot = wallTemperature(baseDot),
            .coolantTemperature = coolantTemperature(base),
            .coolantTemperatureDot = coolantTemperature(baseDot),
            .shellLoading = shellLoading(base),
            .shellLoadingDot = shellLoading(baseDot),
            .streamAmount = streamAmount(base),
            .streamAmountDot = streamAmount(baseDot)};
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
    FixedPressureInletVelocity = 4,   ///< Boundary data: fixed pressure profile and v_in.
    MomentumBalance = 5  ///< Gas: pressure of the gas in each cell, velocity from the pressure-drop law; ports.
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
    WENO3 = 1,   ///< Third-order WENO; gas phase: mole fractions at the faces, total molar flux upwind.
    TVD = 2,     ///< Second-order MUSCL with the minmod limiter, on the quantities and faces of WENO3.
    WENO5 = 3    ///< Fifth-order WENO on a grid of one cell width, WENO3 next to the ends and on other grids.
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

    /**
     * \brief The pressure in Pa and its rate in Pa/s at the given time in s, for a shape other than None.
     */
    [[nodiscard]] std::pair<double, double> at(double atTime) const noexcept
    {
      const double elapsed = atTime - startTime;
      const double change = targetPressure - startPressure;
      if (shape == Shape::Exponential)
      {
        const double decay = std::exp(-elapsed / timeConstant);
        return {targetPressure - change * decay, change * decay / timeConstant};
      }
      if (elapsed < timeConstant) return {startPressure + change * elapsed / timeConstant, change / timeConstant};
      return {targetPressure, 0.0};
    }
  };

  /**
   * \brief An end of the column with MomentumBalance: closed, open to gas at a pressure, or with a set velocity.
   *
   * Gas that flows in is that of node 0 at z = 0 and the product-end gas at z = L; gas that flows out is that of the
   * end cell. A pressure port takes the velocity that the pressure-drop law gives over the half cell between its gas
   * and the centre of the end cell, so that the gas flows either way; its pressure can follow a ramp, which stands for
   * the valve or the pump that sets how fast the column follows. A flux port carries a set molar flux, whatever the gas
   * that carries it: the connections of the swing driver to other columns and to buffers.
   */
  struct Port
  {
    enum class Mode
    {
      Closed = 0,
      Pressure = 1,
      Velocity = 2,
      Flux = 3
    };
    Mode mode{Mode::Closed};
    double pressure{0.0};  ///< Pressure of the gas beyond the end in Pa.
    double velocity{0.0};  ///< Interstitial velocity into the column in m/s (Velocity).
    PressureRamp ramp{};   ///< Ramp of the pressure (Pressure).
    double flux{0.0};      ///< Molar flux into the column per column cross-section in mol/(m^2 s) (Flux).
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
    std::span<const Port> ports{};                      ///< MomentumBalance: the ports at z = 0 and at z = L.
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
        fluidWallContactArea(this->numberOfGridPoints + 1),
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
        wenoScratch((this->numberOfGridPoints + 1) * (2 * this->numberOfComponents + 5) + this->numberOfComponents),
        temperatureRateScratch(3 * (this->numberOfGridPoints + 1)),
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
    for (MixturePrediction& mixture : this->physisorptionMixtures) mixture.storeSiteLoadings = false;
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
    // MomentumBalance without swing phases: the feed at the inlet velocity, and gas at the outlet pressure at z = L.
    exhaustPressure = inputReader.exhaustPressure;
    compressionEfficiency = inputReader.compressionEfficiency;
    ports = {Port{Port::Mode::Velocity, inputReader.outletPressure, inputReader.columnEntranceVelocity},
             Port{Port::Mode::Pressure, inputReader.outletPressure}};
    if (std::ranges::all_of(this->components,
                            [](const Component& component) { return component.gasHeatCapacity > 0.0; }))
    {
      gasHeatCapacity.assign(this->numberOfGridPoints + 1, heatCapacityGas);
    }
    // The optional blocks of the state (ColumnStateLayout): coolant, radial cells of the particle, stream amounts.
    streamAmountsEnabled = inputReader.streamAmounts;
    numberOfShells = inputReader.radialCells;
    coolantEnabled = inputReader.coolant;
    coolantMassFlow = inputReader.coolantMassFlowRate;
    coolantHeatCapacity = inputReader.coolantHeatCapacity;
    coolantHeatCapacityPerLength =
        inputReader.coolantDensity * inputReader.coolantHeatCapacity * inputReader.coolantVolume / columnLength;
    coolantCounterCurrent = inputReader.coolantCounterCurrent;
    const size_t nodes = this->numberOfGridPoints + 1;
    state.resize(state.size() + (coolantEnabled ? nodes : 0) + numberOfShells * nodes * numberOfComponents +
                     (streamAmountsEnabled ? 2 * numberOfComponents + 2 : 0),
                 0.0);
    stateDot.resize(state.size(), 0.0);
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
  std::array<Port, 2> ports{};       ///< MomentumBalance: the ports at z = 0 and at z = L; set per swing phase.
  /// Gas that leaves through a pressure port below this pressure costs the work of its compression to it; 0: none.
  double exhaustPressure{0.0};
  double compressionEfficiency{1.0};  ///< Isentropic efficiency of that compression.
  /// Radial cells of the particle with homogeneous surface diffusion (ParticleTransport HSDM); 0: the linear driving
  /// force. The surface diffusivity is that of the mass-transfer coefficient, D_s / R^2 = k / 15.
  size_t numberOfShells{0};
  /// A heat-transfer fluid flows along the outer wall (Coolant), from z = 0 or, counter-current, from z = L, where it
  /// enters at externalTemperature: its mass flow in kg/s, heat capacity in J/(kg K), and heat capacity per column
  /// length, rho c_p A, in J/(m K).
  bool coolantEnabled{false};
  double coolantMassFlow{0.0};
  double coolantHeatCapacity{0.0};
  double coolantHeatCapacityPerLength{0.0};
  bool coolantCounterCurrent{false};
  bool streamAmountsEnabled{false};  ///< The state holds the stream amounts (ColumnStateLayout::includeStreamAmounts).
  [[nodiscard]] FlowEnds flowEnds() const noexcept
  {
    return {
        flowDirection, productEndConcentration, productEndTemperature, closedEnd, pressureRate, relativePressureRate,
        ports};
  }

  /**
   * \brief Sets the inlet pressure, its rates, and the product-end gas to their values at the given time; without a
   *        ramp it does nothing. The integrators call it before each right-hand side, which depends on the time through
   *        these fields only.
   */
  void advancePressureRamp(double time) noexcept
  {
    if (boundaryCondition == BoundaryCondition::MomentumBalance)
    {
      // The gas of a ramped port follows its pressure: that of node 0 through its rate (computeMassDerivatives), the
      // product-end gas directly. The pressure work of the cells is in their own energy balance, not a pressure rate.
      pressureRamp.time = time;
      relativePressureRate = 0.0;
      for (size_t end = 0; end < ports.size(); ++end)
      {
        Port& port = ports[end];
        if (port.ramp.shape == PressureRamp::Shape::None) continue;
        const auto [pressure, rate] = port.ramp.at(time);
        if (end == 0) relativePressureRate = rate / pressure;
        for (size_t comp = 0; end == 1 && comp < productEndConcentration.size(); ++comp)
        {
          productEndConcentration[comp] *= pressure / port.pressure;
        }
        port.pressure = pressure;
      }
      return;
    }
    if (pressureRamp.shape == PressureRamp::Shape::None) return;
    pressureRamp.time = time;
    const auto [pressure, rate] = pressureRamp.at(time);
    pressureRate = rate;
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
  /// Heat capacity of the gas of each node in J/(kg K), from the GasHeatCapacity of its components (precompute); empty
  /// when they give none: the gas then has heatCapacityGas everywhere.
  std::vector<double> gasHeatCapacity;
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
  std::vector<double> fluidWallContactArea;      ///< Fluid-wall contact area per fluid volume in 1/m.
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
  /// Work array of computeMassDerivatives; size (numberOfGridPoints + 1) * (2 numberOfComponents + 5) +
  /// numberOfComponents.
  std::vector<double> wenoScratch;
  std::vector<double> temperatureRateScratch;  ///< Work array of the velocity closures with the energy balance.
  /// Chemisorption sites of each component, the largest count over the adsorbents; set by initialize().
  std::vector<size_t> chemisorptionSites;

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
   *        per-adsorbent properties; called by initialize(), and needed again after changing either.
   */
  void computeAdsorbentBlends();

  /**
   * \brief Cross-section of the column in m^2, to which the amounts per unit cross-section refer: that of the bed of
   *        a packed tube, and the outer one of a monolith; 0 without the diameter.
   */
  [[nodiscard]] double crossSection() const noexcept
  {
    const Geometry& geometry = geometries.front();
    return geometry.kind == GeometryKind::Monolith ? geometry.dimensions.outerArea : geometry.dimensions.openArea;
  }

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
   * \brief Serializes the column state (concentrations, loadings, temperatures) and its time in s to JSON.
   */
  void writeJSON(const std::string& fileName, double time = 0.0) const;

  /**
   * \brief Loads the column state from JSON and returns its time in s (0 for a file without one). Other keys, which
   *        older files hold, are not read: the derived fields follow from the state.
   */
  double readJSON(const std::string& fileName);
};

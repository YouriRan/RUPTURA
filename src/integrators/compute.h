#pragma once
#include <cstddef>
#include <span>
#include <vector>

#include "column.h"

/**
 * \brief The gas-temperature balance of computeEnergyDerivatives: dT_f/dt of a node for given velocities at its faces.
 *
 * computeEnergyDerivatives evaluates the gas temperature with it, and the velocity closures with the thermal
 * expansion of the gas (updateVelocityAndPressure).
 */
struct GasTemperatureBalance
{
  std::span<const double> columnDistances;
  std::span<const double> gasTemperature;
  std::span<const double> solidTemperature;
  std::span<const double> wallTemperature;
  std::span<const double> gasDensity;
  std::span<const double> fluidSolidContactArea;
  std::span<const double> totalVoidFraction;
  double gasThermalConductivity;
  double heatTransferGasSolid;
  double heatCapacityGas;
  double gasWallPrefactor;  ///< Gas-wall area per fluid volume times h_gw / c_p,g.
  Column::InletCondition inletCondition;
  Column::FlowEnds ends{};  ///< Flow direction and product-end gas temperature.

  /**
   * \brief Returns dT_f/dt of node grid (1 to numberOfGridPoints) with the velocities v_{i-1} (lowerVelocity) and v_i
   *        (upperVelocity) at its faces: gas-solid and gas-wall exchange, upwind convection, and conduction.
   *
   * The velocity v_k of face z_k is the interstitial velocity of node k there (computeMassDerivatives), so node i
   * flows through its lower face at w_{i-1} = eps_{i-1} v_{i-1} / eps_i, the same superficial velocity. Convection is
   * -max(v_i, 0) (T_i - T_{i-1}) / h_i + min(w_{i-1}, 0) (T_i - T_{i+1}) / h_i: the velocity of the downstream face
   * with the upstream neighbour, and nothing in a node whose faces both carry gas away. With a Danckwerts inlet, the
   * node next to an end with inflow takes the face temperature T_b, with u (T_b - T_res) = g (T_node - T_b) and u the
   * inflow speed of its inner face (v_1 at z = 0, -w_{N-1} at z = L), in place of the reservoir temperature, so that no
   * heat is conducted across that face; T_b depends on u. At an end with outflow, and at z = 0 in a reverse step, the
   * gradient is zero. A pressure ramp adds the pressure work (dP/dt) / (rho c_p).
   */
  [[nodiscard]] double rate(size_t grid, double lowerVelocity, double upperVelocity) const;
};

/**
 * \brief Updates the velocity, pressure, total concentration, mole fractions and pH fields from the boundary
 *        settings and the bulk species sink, and the axial dispersion from its law.
 *
 * The axial dispersion is D_j = D_0,j theta + lambda |v| + kappa v^2 / (D_0,j theta), with the molecular part D_0,j
 * (molecularDispersionCoefficient), the dispersivity lambda and Taylor-Aris coefficient kappa (dispersionLaw), and
 * theta = (T/T_ref)^1.75 (P_ref/P) when a reference temperature is set (1 otherwise). The gas continuity closures
 * conserve the total superficial molar flux eps v sum_j c_j of the species balances, uptake and net axial dispersion
 * included, so that sum_j c_j does not drift when the D_j differ; convectionScheme and inletCondition select the
 * dispersion stencils that computeMassDerivatives uses. With gasTemperatureBalance (the energy balance) they also carry
 * the thermal expansion of the gas, d(sum_j c_ij)/dt = -(c_i/T_i) dT_i/dt, so that sum_j c_j T stays constant at
 * constant pressure. The velocity of face z_i is the interstitial velocity of node i there: where the void fraction
 * eps changes, the superficial velocity eps v is continuous and v is not (C-MBD-1). FixedVelocity, and the liquid
 * phase, hold eps v at its value at the feed end, eps_in v_in.
 *
 * In a reverse step (ends) the feed enters at z = L: the inlet velocity and pressure apply at z = L, the outlet pressure
 * at z = 0, the closures march from z = L, and the velocities are negative. The molar flux through a face is its
 * velocity times the total concentration upstream of it, so that a closure that reverses the flow inside the column
 * (gas that flows back in at the outlet end) takes the gas of the node beyond the face, or the product-end gas.
 *
 * With a closed end (ends) the gas continuity closure marches from that end with v = 0, and the pressure follows from
 * inletPressure at the open end, whatever the boundary condition. A pressure ramp (ends) adds the accumulation
 * C_i (dP/dt) / P of every node, which neglects the time derivative of the pressure drop: the pressure drop must stay
 * small against the lowest pressure of the ramp.
 */
void updateVelocityAndPressure(
    const std::vector<Component>& components, const Column::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const PressureDropLaw> pressureDropLaw,
    std::span<const double> totalVoidFraction, std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature, Column::FluidPhase fluidPhase,
    double liquidDensity, Column::PHMode pHMode, double pHValue, double pKw, size_t pHComponent,
    std::span<double> pH, std::span<const double> molecularDispersionCoefficient,
    std::span<const DispersionLaw> dispersionLaw, double dispersionReferenceTemperature,
    double dispersionReferencePressure, std::span<double> axialDispersionCoefficient,
    Column::ConvectionScheme convectionScheme, Column::InletCondition inletCondition,
    const GasTemperatureBalance* gasTemperatureBalance = nullptr, const Column::FlowEnds& ends = {});

/**
 * \brief Adds reaction source terms to explicit adsorbed/pore state arrays.
 *
 * BulkConcentration reactions write their source, per unit fluid volume, into reactionBulkSource; the caller
 * subtracts it from the bulk species sink. Their heat enters reactionHeat (per kg of solid) divided by
 * loadingPrefactor, the solid mass per fluid volume. A reaction restricted to Adsorbents is weighted by their local
 * fraction (fractionOfAdsorbent, grid-major, numberOfAdsorbents per node). coverage is a work array of
 * numberOfComponents values for the Langmuir-Hinshelwood coverages.
 */
void computeReactionDerivatives(
    const std::vector<Component>& components, const std::vector<Reaction>& reactions,
    size_t numberOfGridPoints, size_t numberOfComponents, size_t maxChemisorptionSites,
    double externalTemperature, std::span<const double> physisorption,
    std::span<double> physisorptionDot, std::span<const double> chemisorption,
    std::span<double> chemisorptionDot, std::span<const double> poreConcentration,
    std::span<double> poreConcentrationDot, std::span<const double> solidTemperature,
    std::span<double> reactionPhysisorptionSource,
    std::span<double> reactionChemisorptionSource,
    std::span<double> reactionPoreConcentrationSource, std::span<double> reactionHeat,
    std::span<const double> concentration, std::span<double> reactionBulkSource,
    std::span<const double> fractionOfAdsorbent, size_t numberOfAdsorbents, std::span<const double> loadingPrefactor,
    std::span<double> coverage);

/**
 * \brief Computes the concentration derivatives: convection, axial dispersion and the bulk species sink.
 *
 * The convective fluxes are superficial, eps_i v_i c through face z_i, with v_i the interstitial velocity of node i,
 * and node i divides them by its fluid volume eps_i h_i. With WENO3 the convective face values are reconstructed with
 * computeWENO: in the gas phase the mole fractions, with the total molar flux eps v c_tot kept upwind (the flux that
 * the velocity closures conserve), in the liquid phase the concentrations. Node 0 holds the feed concentrations, which
 * follow a pressure ramp in a forward step; with a Danckwerts inlet, the stencils use the inlet-face value at which the
 * convective-dispersive flux into the first cell equals the feed flux. wenoScratch is the work array of WENO3,
 * (numberOfGridPoints + 1) * (numberOfComponents + 4) values; Upwind does not use it.
 *
 * The face values are upwinded by the sign of the velocity of each face (computeWENOReverse for v < 0). Gas that flows
 * in at z = L is the product-end gas of ends, or that of the last cell, with the same Danckwerts or Dirichlet condition
 * as at z = 0. In a reverse step (ends) the feed enters at z = L, and z = 0 has zero gradient: node 0 is replaced by
 * the first cell.
 */
void computeMassDerivatives(size_t numberOfGridPoints, size_t numberOfComponents,
                            std::span<const double> columnDistances, std::span<const double> axialDispersionCoefficient,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink,
                            Column::FluidPhase fluidPhase, std::span<const double> totalVoidFraction,
                            Column::ConvectionScheme convectionScheme,
                            Column::InletCondition inletCondition = Column::InletCondition::Danckwerts,
                            std::span<double> wenoScratch = {}, const Column::FlowEnds& ends = {});

/**
 * \brief Computes the gas, solid, and wall temperature derivatives.
 *
 * Node 0 holds the feed temperature; with a Danckwerts inlet, the upwind stencil of node 1 uses the inlet-face
 * temperature at which no heat is conducted across the inlet face. The gas temperature follows the flow direction and
 * the product-end gas of ends (GasTemperatureBalance::rate).
 */
void computeEnergyDerivatives(
    size_t numberOfGridPoints, size_t numberOfComponents, double externalTemperature,
    std::span<const Geometry> geometries, std::span<const double> particleDensity,
    std::span<const double> solidFluidContactArea, std::span<const double> fluidSolidContactArea,
    std::span<const double> heatOfAdsorption, std::span<const double> heatOfChemisorption, double wallDensity,
    double gasThermalConductivity, double wallThermalConductivity, double heatTransferGasSolid,
    double heatTransferGasWall, double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
    double heatCapacityWall, std::span<const double> columnDistances, std::span<const double> interstitialGasVelocity,
    std::span<const double> totalVoidFraction, std::span<const double> gasDensity, std::span<double> coeffDiffusion,
    size_t maxChemisorptionSites,
    std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot,
    std::span<const double> reactionPhysisorptionSource = {}, std::span<const double> reactionChemisorptionSource = {},
    std::span<const double> reactionHeat = {}, Column::InletCondition inletCondition = Column::InletCondition::Danckwerts,
    const Column::FlowEnds& ends = {});

/**
 * \brief Third-order WENO face values (v > 0) on the grid with node positions columnDistances (z_0 ... z_N).
 *
 * values[0] is the inlet value at the inlet face z_0, and values[i], i >= 1, the value of the cell (z_{i-1}, z_i].
 * faces[i] is the value at the face z_i reconstructed from the upstream side: faces[0] = values[0], and the outlet
 * face takes the last cell value (first order).
 */
void computeWENO(std::span<const double> values, std::span<const double> columnDistances, std::span<double> faces);

/**
 * \brief Derivatives of the computeWENO face values with respect to the node values, nonlinear weights included:
 *        face i changes by lower[i] du_{i-1} + centre[i] du_i + upper[i] du_{i+1}; the inlet and outlet faces, which
 *        take the node values, have centre 1.
 */
void computeWENOLinearization(std::span<const double> values, std::span<const double> columnDistances,
                              std::span<double> lower, std::span<double> centre, std::span<double> upper);

/**
 * \brief The mirror image of computeWENO, for flow towards z = 0 (v < 0).
 *
 * values[i], i >= 1, are the cell values as in computeWENO (values[0] is not used), and inflow is the value at the face
 * z_N. faces[i] is the value at the face z_i reconstructed from the side of z = L: faces[N] = inflow, and the face z_0
 * takes the first cell value (first order). The result is that of computeWENO on the mirrored grid.
 */
void computeWENOReverse(std::span<const double> values, double inflow, std::span<const double> columnDistances,
                        std::span<double> faces);

/**
 * \brief Derivatives of the computeWENOReverse face values: face i changes by downstream[i] du_i + centre[i] du_{i+1}
 *        + upstream[i] du_{i+2}, where u_{N+1} stands for the inflow value; faces[N] has centre 1 (the inflow), and
 *        faces[0] centre 1 (the first cell).
 */
void computeWENOReverseLinearization(std::span<const double> values, double inflow,
                                     std::span<const double> columnDistances, std::span<double> downstream,
                                     std::span<double> centre, std::span<double> upstream);

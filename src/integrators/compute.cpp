#include "compute.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mdspan>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

namespace
{
constexpr double unboundedLoading = std::numeric_limits<double>::max();

bool isUnbounded(double value) noexcept
{
  return value >= 0.5 * unboundedLoading;
}

double finiteNonNegative(double value) noexcept
{
  if (!std::isfinite(value)) return 0.0;
  return std::max(0.0, value);
}

double isothermMaximumLoading(const Isotherm& isotherm) noexcept
{
  if (isotherm.parameters.empty()) return unboundedLoading;

  switch (isotherm.type)
  {
    case Isotherm::Type::Langmuir:
    case Isotherm::Type::Langmuir_pH:
    case Isotherm::Type::Sips:
    case Isotherm::Type::Langmuir_Freundlich:
    case Isotherm::Type::Toth:
    case Isotherm::Type::Unilan:
    case Isotherm::Type::OBrien_Myers:
    case Isotherm::Type::Quadratic:
    case Isotherm::Type::Temkin:
    case Isotherm::Type::BingelWalton:
    // The humidity modifiers of these types change the affinity or reduce the capacity.
    case Isotherm::Type::Langmuir_humidity:
    case Isotherm::Type::Langmuir_moisture_swing:
    case Isotherm::Type::GAB_Langmuir:
    case Isotherm::Type::GAB_Sips:
    case Isotherm::Type::GAB_Toth:
    {
      const double maximum = finiteNonNegative(isotherm.parameters.front());
      return maximum > 0.0 ? maximum : unboundedLoading;
    }
    case Isotherm::Type::Anti_Langmuir:
    case Isotherm::Type::BET:
    case Isotherm::Type::Henry:
    case Isotherm::Type::Freundlich:
    case Isotherm::Type::Redlich_Peterson:
    case Isotherm::Type::GAB:
    // The capacity of these types grows with the water loading or the humidity.
    case Isotherm::Type::GAB_modified_Toth:
    case Isotherm::Type::Toth_exponential_humidity:
      return unboundedLoading;
  }
  return unboundedLoading;
}

double physisorptionMaximumLoading(const Component& component) noexcept
{
  double maximumLoading = 0.0;
  bool hasBoundedSite = false;
  for (const Isotherm& site : component.isotherm.sites)
  {
    const double siteMaximum = isothermMaximumLoading(site);
    if (isUnbounded(siteMaximum)) return unboundedLoading;
    maximumLoading += siteMaximum;
    hasBoundedSite = true;
  }
  return hasBoundedSite && maximumLoading > 0.0 ? maximumLoading : unboundedLoading;
}

double chemisorptionMaximumLoading(const Component& component, size_t site) noexcept
{
  if (site >= component.chemisorption.numberOfSites) return unboundedLoading;

  const Chemisorption& kinetics = component.chemisorption.sites[site];
  if (std::isfinite(kinetics.maximumLoading) && kinetics.maximumLoading > 0.0)
  {
    return kinetics.maximumLoading;
  }
  if (kinetics.isotherm.has_value())
  {
    return isothermMaximumLoading(*kinetics.isotherm);
  }
  return unboundedLoading;
}

void limitByInventory(double& limit, double amount, double stoichiometry, double limitTime) noexcept
{
  if (stoichiometry <= 0.0) return;

  amount = finiteNonNegative(amount);
  if (amount <= 0.0)
  {
    limit = 0.0;
    return;
  }
  if (limitTime > 0.0)
  {
    limit = std::min(limit, amount / (stoichiometry * limitTime));
  }
}

void limitByCapacity(double& limit, double maximumLoading, double loading, double stoichiometry,
                     double limitTime) noexcept
{
  if (stoichiometry <= 0.0 || isUnbounded(maximumLoading)) return;
  limitByInventory(limit, maximumLoading - loading, stoichiometry, limitTime);
}

double limitedAdsorbedRate(const Reaction& reaction, std::span<const double> activity,
                           std::span<const double> maximumLoadings, double rate) noexcept
{
  if (!std::isfinite(rate) || rate == 0.0) return 0.0;

  const double limitTime = reaction.rateLimitTime;
  double extentLimit = std::abs(rate);

  if (rate > 0.0)
  {
    for (size_t i = 0; i < reaction.reactants.size(); ++i)
    {
      const size_t comp = reaction.reactants[i];
      const double stoichiometry = i < reaction.reactantStoichiometry.size() ? reaction.reactantStoichiometry[i] : 1.0;
      limitByInventory(extentLimit, comp < activity.size() ? activity[comp] : 0.0, stoichiometry, limitTime);
    }
    for (size_t i = 0; i < reaction.products.size(); ++i)
    {
      const size_t comp = reaction.products[i];
      const double stoichiometry = i < reaction.productStoichiometry.size() ? reaction.productStoichiometry[i] : 1.0;
      const double loading = comp < activity.size() ? activity[comp] : 0.0;
      const double maximumLoading = comp < maximumLoadings.size() ? maximumLoadings[comp] : unboundedLoading;
      limitByCapacity(extentLimit, maximumLoading, loading, stoichiometry, limitTime);
    }
    return std::min(rate, extentLimit);
  }

  for (size_t i = 0; i < reaction.products.size(); ++i)
  {
    const size_t comp = reaction.products[i];
    const double stoichiometry = i < reaction.productStoichiometry.size() ? reaction.productStoichiometry[i] : 1.0;
    limitByInventory(extentLimit, comp < activity.size() ? activity[comp] : 0.0, stoichiometry, limitTime);
  }
  for (size_t i = 0; i < reaction.reactants.size(); ++i)
  {
    const size_t comp = reaction.reactants[i];
    const double stoichiometry = i < reaction.reactantStoichiometry.size() ? reaction.reactantStoichiometry[i] : 1.0;
    const double loading = comp < activity.size() ? activity[comp] : 0.0;
    const double maximumLoading = comp < maximumLoadings.size() ? maximumLoadings[comp] : unboundedLoading;
    limitByCapacity(extentLimit, maximumLoading, loading, stoichiometry, limitTime);
  }
  return -std::min(-rate, extentLimit);
}

double limitedPoreRate(const Reaction& reaction, std::span<const double> activity, double rate) noexcept
{
  if (!std::isfinite(rate) || rate == 0.0) return 0.0;

  const double limitTime = reaction.rateLimitTime;
  double extentLimit = std::abs(rate);
  const std::vector<size_t>& consumed = rate > 0.0 ? reaction.reactants : reaction.products;
  const std::vector<double>& stoichiometry =
      rate > 0.0 ? reaction.reactantStoichiometry : reaction.productStoichiometry;

  for (size_t i = 0; i < consumed.size(); ++i)
  {
    const size_t comp = consumed[i];
    const double coefficient = i < stoichiometry.size() ? stoichiometry[i] : 1.0;
    limitByInventory(extentLimit, comp < activity.size() ? activity[comp] : 0.0, coefficient, limitTime);
  }
  return std::copysign(std::min(std::abs(rate), extentLimit), rate);
}

}  // namespace

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
    std::span<double> coverage)
{
  std::fill(reactionPhysisorptionSource.begin(), reactionPhysisorptionSource.end(), 0.0);
  std::fill(reactionChemisorptionSource.begin(), reactionChemisorptionSource.end(), 0.0);
  std::fill(reactionPoreConcentrationSource.begin(), reactionPoreConcentrationSource.end(), 0.0);
  std::fill(reactionBulkSource.begin(), reactionBulkSource.end(), 0.0);
  std::fill(reactionHeat.begin(), reactionHeat.end(), 0.0);

  if (reactions.empty()) return;

  mdspan2d_const spanPhysisorption(physisorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanPhysisorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanReactionPhysisorptionSource(reactionPhysisorptionSource.data(),
                                               numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorption(chemisorption.data(), maxChemisorptionSites,
                                   numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites,
                                    numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanReactionChemisorptionSource(reactionChemisorptionSource.data(),
                                               maxChemisorptionSites, numberOfGridPoints + 1,
                                               numberOfComponents);
  mdspan3d_const spanPoreConcentration(poreConcentration.data(), maxChemisorptionSites,
                                       numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanPoreConcentrationDot(poreConcentrationDot.data(), maxChemisorptionSites,
                                        numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanReactionPoreConcentrationSource(reactionPoreConcentrationSource.data(),
                                                   maxChemisorptionSites, numberOfGridPoints + 1,
                                                   numberOfComponents);
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanReactionBulkSource(reactionBulkSource.data(), numberOfGridPoints + 1, numberOfComponents);

  std::vector<double> activity(numberOfComponents, 0.0);
  std::vector<double> maximumLoadings(numberOfComponents, unboundedLoading);

  for (const Reaction& reaction : reactions)
  {
    if (!reaction.enabled()) continue;
    // The reader checks the reaction sites; the column holds the pore state whenever a reaction needs it.

    if (reaction.phase == Reaction::Phase::Physisorbed)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        maximumLoadings[comp] = physisorptionMaximumLoading(components[comp]);
      }
    }
    else if (reaction.phase == Reaction::Phase::Chemisorbed)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        maximumLoadings[comp] = chemisorptionMaximumLoading(components[comp], reaction.site);
      }
    }

    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      // A reaction restricted to some adsorbents (e.g. a catalyst section) acts with their local fraction.
      double sectionWeight = 1.0;
      if (!reaction.adsorbents.empty())
      {
        sectionWeight = 0.0;
        for (size_t ads : reaction.adsorbents)
        {
          sectionWeight += fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
        }
        if (sectionWeight <= 0.0) continue;
      }

      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        switch (reaction.phase)
        {
          case Reaction::Phase::Physisorbed:
            activity[comp] = spanPhysisorption[grid, comp];
            break;
          case Reaction::Phase::Chemisorbed:
            activity[comp] = spanChemisorption[reaction.site, grid, comp];
            break;
          case Reaction::Phase::PoreConcentration:
            activity[comp] = spanPoreConcentration[reaction.site, grid, comp];
            break;
          case Reaction::Phase::BulkConcentration:
            activity[comp] = spanConcentration[grid, comp];
            break;
        }
      }

      const double temperature = solidTemperature.empty() ? externalTemperature : solidTemperature[grid];
      double rate = sectionWeight * reaction.rate(activity, temperature, coverage);
      if (reaction.phase == Reaction::Phase::PoreConcentration ||
          reaction.phase == Reaction::Phase::BulkConcentration)
      {
        rate = limitedPoreRate(reaction, activity, rate);
      }
      else
      {
        rate = limitedAdsorbedRate(reaction, activity, maximumLoadings, rate);
      }
      if (rate == 0.0) continue;

      // Reactants are consumed (sign -1), products formed (+1); a missing stoichiometric coefficient is 1.
      for (const auto& [species, stoichiometries, sign] :
           {std::tuple{&reaction.reactants, &reaction.reactantStoichiometry, -1.0},
            std::tuple{&reaction.products, &reaction.productStoichiometry, 1.0}})
      {
        for (size_t i = 0; i < species->size(); ++i)
        {
          const size_t comp = (*species)[i];
          const double stoichiometry = i < stoichiometries->size() ? (*stoichiometries)[i] : 1.0;
          const double source = sign * stoichiometry * rate;
          switch (reaction.phase)
          {
            case Reaction::Phase::Physisorbed:
              spanPhysisorptionDot[grid, comp] += source;
              spanReactionPhysisorptionSource[grid, comp] += source;
              break;
            case Reaction::Phase::Chemisorbed:
              spanChemisorptionDot[reaction.site, grid, comp] += source;
              spanReactionChemisorptionSource[reaction.site, grid, comp] += source;
              break;
            case Reaction::Phase::PoreConcentration:
              spanPoreConcentrationDot[reaction.site, grid, comp] += source;
              spanReactionPoreConcentrationSource[reaction.site, grid, comp] += source;
              break;
            case Reaction::Phase::BulkConcentration:
              // Applied to the bulk balance through computeBulkSpeciesSink, per unit fluid volume.
              spanReactionBulkSource[grid, comp] += source;
              break;
          }
        }
      }

      // reactionHeat is per kg of solid. A bulk reaction's rate is per fluid volume; its heat is released into the
      // solid (the catalyst), so it is divided by the solid mass per fluid volume, as in the bulk species sink.
      const double heat = reaction.heatOfReaction * rate;
      reactionHeat[grid] += reaction.phase == Reaction::Phase::BulkConcentration ? heat / loadingPrefactor[grid] : heat;
    }
  }
}

void updateVelocityAndPressure(
    const std::vector<Component>& components, const Column::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const PressureDropLaw> pressureDropLaw,
    std::span<const double> totalVoidFraction, std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature,
    Column::FluidPhase fluidPhase, double liquidDensity, Column::PHMode pHMode,
    double pHValue, double pKw, size_t pHComponent, std::span<double> pH,
    std::span<const double> molecularDispersionCoefficient, std::span<const DispersionLaw> dispersionLaw,
    double dispersionReferenceTemperature, double dispersionReferencePressure,
    std::span<double> axialDispersionCoefficient, Column::ConvectionScheme convectionScheme,
    Column::InletCondition inletCondition, const GasTemperatureBalance* gasTemperatureBalance,
    const Column::FlowEnds& ends)
{
  const bool reverse = ends.direction == Column::FlowDirection::Reverse;
  const bool danckwerts = inletCondition == Column::InletCondition::Danckwerts;
  const size_t outlet = numberOfGridPoints;
  // The continuity closures march from the end where the feed enters, or from the closed end, which is the other one.
  const bool marchFromLength = ends.closedEnd != reverse;

  auto refreshNode = [&](size_t grid)
  {
    totalConcentration[grid] = totalPressure[grid] / (R * std::max(1e-10, gasTemperature[grid]));
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }

    gasDensity[grid] = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = moleFraction[index] * totalPressure[grid];
      gasDensity[grid] += moleFraction[index] * totalConcentration[grid] * components[comp].molecularWeight;
    }
  };

  auto ergunGrad = [&](size_t grid)
  {
    double grad = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ads++)
    {
      grad += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
              geometries[ads].pressureDrop.gradient(dynamicViscosity, gasDensity[grid],
                                                    interstitialGasVelocity[grid]);
    }
    return grad;
  };

  // Uptake eps_i h_i sum_j S_ij of node i per column cross-section: the closures balance superficial molar fluxes, and
  // the sink is per fluid volume. A pressure ramp adds the accumulation eps_i h_i C_i (dP/dt) / P, the (dP/dt) / (R T_i)
  // of a node at the column pressure, so that C_i / P keeps its value and an offset from the pressure drop stays small
  // when the pressure falls.
  auto sinkTerm = [&](size_t grid)
  {
    const auto begin = static_cast<std::ptrdiff_t>(grid * numberOfComponents);
    const auto end = static_cast<std::ptrdiff_t>((grid + 1) * numberOfComponents);
    return totalVoidFraction[grid] * gridSpacing(columnDistances, grid) *
           (std::reduce(bulkSpeciesSink.begin() + begin, bulkSpeciesSink.begin() + end) +
            ends.relativePressureRate * totalConcentration[grid]);
  };

  // Axial dispersion D_j = D_0,j theta + lambda |v| + kappa v^2 / (D_0,j theta). The Fuller scaling
  // theta = (T/T_ref)^1.75 (P_ref/P) of the molecular part takes the gas temperature and the gas pressure
  // R T sum_j c_j of the node, so that the dispersion is a function of the state.
  auto setAxialDispersion = [&](size_t grid, double velocity)
  {
    double scaling = 1.0;
    if (dispersionReferenceTemperature > 0.0)
    {
      const double temperature = std::max(1e-10, gasTemperature[grid]);
      double concentrationSum = 0.0;
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
      }
      scaling = std::pow(temperature / dispersionReferenceTemperature, 1.75) * dispersionReferencePressure /
                std::max(1e-10, R * temperature * concentrationSum);
    }
    const DispersionLaw& law = dispersionLaw[grid];
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      const double molecular = scaling * molecularDispersionCoefficient[index];
      double dispersion = molecular + law.dispersivity * std::abs(velocity);
      if (law.taylorArisCoefficient > 0.0)
      {
        dispersion += law.taylorArisCoefficient * velocity * velocity / molecular;
      }
      axialDispersionCoefficient[index] = dispersion;
    }
  };

  // Dispersive inflow eps_i h_i sum_j D_ij d2c_ij/dz2 of node i on the stencils of computeMassDerivatives. The
  // continuity closures add it to the uptake, so that they balance the same total flux as the species balances and
  // sum_j c_j stays constant also when the D_j differ (C-BC-1; RUPTURA 1 had this term). With a Danckwerts inlet the
  // convective and dispersive fluxes through the inlet face add up to the feed flux eps_0 v_0 c_0 that the closure
  // starts from, so the inlet face is left out (inlet value c_1); the same holds for gas that flows in at z = L, whose
  // total flux is eps_N v_N times the product-end gas. In a reverse step z = 0 has zero gradient, and a Dirichlet inlet
  // at z = L conducts towards the feed, as at z = 0 in a forward step.
  const double lastSpacing = gridSpacing(columnDistances, outlet);
  const double outletDistance = convectionScheme == Column::ConvectionScheme::WENO3 ? 0.5 * lastSpacing
                                : outlet > 1 ? 0.5 * (lastSpacing + gridSpacing(columnDistances, outlet - 1))
                                             : lastSpacing;
  auto dispersionTerm = [&](size_t grid)
  {
    const double hm = std::max(1e-30, gridSpacing(columnDistances, grid));
    double term = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double centre = concentration[grid * numberOfComponents + comp];
      const double upstream =
          grid == 1 && (danckwerts || reverse) ? centre : concentration[(grid - 1) * numberOfComponents + comp];
      double d2cDz2 = 0.0;
      if (grid == numberOfGridPoints)
      {
        d2cDz2 = (upstream - centre) / (hm * hm);
        if (reverse && !danckwerts)
        {
          const double feed = ends.productEndConcentration.empty() ? centre : ends.productEndConcentration[comp];
          d2cDz2 += (feed - centre) / std::max(1e-30, lastSpacing * outletDistance);
        }
      }
      else
      {
        const double downstream = concentration[(grid + 1) * numberOfComponents + comp];
        const double hp = gridSpacing(columnDistances, grid + 1);
        d2cDz2 = grid == 1 && convectionScheme == Column::ConvectionScheme::WENO3
                     ? ((downstream - centre) / (0.5 * (hm + hp)) - (centre - upstream) / (0.5 * hm)) / hm
                     : 2.0 * (hm * downstream - (hm + hp) * centre + hp * upstream) /
                           std::max(1e-30, hm * hp * (hm + hp));
      }
      term += axialDispersionCoefficient[grid * numberOfComponents + comp] * d2cDz2;
    }
    return totalVoidFraction[grid] * gridSpacing(columnDistances, grid) * term;
  };

  // Velocity of a face of node i in a march through node i: of its upper face z_i in a march from z = 0, from
  // eps_i v_i c_up = flux + eps_i h_i (c_i/T_i) dT_i/dt(v_{i-1}, v_i), and of its lower face z_{i-1} in a march from
  // z = L (fromLength), from eps_{i-1} v_{i-1} c_up = flux - eps_i h_i (c_i/T_i) dT_i/dt(v_{i-1}, v_i). The velocity
  // of face z_k is the interstitial velocity of node k there, so that eps_k v_k is the superficial velocity of the face
  // (C-MBD-1). flux is the superficial molar flux of the species balance of node i through that face without the
  // expansion, and c_up the total concentration upstream of the face: that of node i (nodeConcentration) for flow away
  // from it, of the node beyond the face (otherConcentration) otherwise. With the energy balance the gas of node i
  // expands as it heats, d(sum_j c_ij)/dt = -(c_i/T_i) dT_i/dt (B-BC-2), which the flux leaving node i carries.
  // dT_i/dt = a_i + b_i max(v_i, 0) + r_i min(v_{i-1}, 0) is piecewise linear, and its coefficients do not depend on
  // the march, except at the node next to a Danckwerts inlet, whose face temperature depends on the velocity of its
  // other face: there the root is found by secant steps from the velocity without expansion.
  std::vector<double> temperatureRateAtRest;   // a_i, filled once the gas densities are updated
  std::vector<double> temperatureRateSlope;    // b_i
  std::vector<double> temperatureRateReverse;  // r_i, in a reverse step
  const bool nonlinearFirstNode = danckwerts && !reverse;
  const bool nonlinearLastNode = danckwerts && reverse;
  auto faceVelocity = [&](size_t grid, double superficialFlux, double nodeConcentration, double otherConcentration,
                          bool fromLength)
  {
    const double direction = fromLength ? -1.0 : 1.0;
    const double faceFraction = totalVoidFraction[fromLength ? grid - 1 : grid];
    const double flux = superficialFlux / faceFraction;
    const double total = std::max(1e-10, nodeConcentration);
    const double other = std::max(1e-10, otherConcentration);
    const double velocity = flux / total;
    if (gasTemperatureBalance == nullptr) return direction * velocity >= 0.0 ? velocity : flux / other;

    const double lowerVelocity = interstitialGasVelocity[grid - 1];
    const double upperVelocity = interstitialGasVelocity[grid];
    const double scale = totalVoidFraction[grid] / faceFraction * gridSpacing(columnDistances, grid) * total /
                         std::max(1e-10, gasTemperature[grid]);
    if (fromLength ? grid == outlet && nonlinearLastNode : grid == 1 && nonlinearFirstNode)
    {
      auto residual = [&](double v)
      {
        const double rate = fromLength ? gasTemperatureBalance->rate(grid, v, upperVelocity)
                                       : gasTemperatureBalance->rate(grid, lowerVelocity, v);
        return v * (direction * v >= 0.0 ? total : other) - flux - direction * scale * rate;
      };
      double previous = velocity;
      double previousResidual = residual(previous);
      double current = velocity + std::max(1.0e-7 * std::abs(velocity), 1.0e-12);
      double currentResidual = residual(current);
      for (size_t iteration = 0; iteration < 50 && currentResidual != previousResidual; ++iteration)
      {
        const double next = current - currentResidual * (current - previous) / (currentResidual - previousResidual);
        previous = current;
        previousResidual = currentResidual;
        current = next;
        currentResidual = residual(current);
        if (std::abs(current - previous) <= 1.0e-14 * std::max(std::abs(current), 1.0e-12)) break;
      }
      return current;
    }

    // dT_i/dt = atRest + slope v for flow away from node i through the face.
    double atRest = 0.0;
    double slope = 0.0;
    if (fromLength)
    {
      // Marching from a closed z = L, the face temperature of a Danckwerts node 1 depends on v_1, which is known here.
      atRest = grid == 1 && nonlinearFirstNode
                   ? gasTemperatureBalance->rate(1, 0.0, upperVelocity)
                   : temperatureRateAtRest[grid] + std::max(0.0, upperVelocity) * temperatureRateSlope[grid];
      slope = temperatureRateReverse[grid];
    }
    else if (grid == outlet)
    {
      atRest = gasTemperatureBalance->rate(grid, lowerVelocity, 0.0);
      slope = gasTemperatureBalance->rate(grid, lowerVelocity, 1.0) - atRest;
    }
    else
    {
      atRest = temperatureRateAtRest[grid];
      if (lowerVelocity < 0.0)
      {
        atRest += lowerVelocity * (temperatureRateAtRest[grid] - gasTemperatureBalance->rate(grid, -1.0, 0.0));
      }
      slope = temperatureRateSlope[grid];
    }
    const double away = (flux + direction * scale * atRest) / (total - direction * scale * slope);
    return direction * away >= 0.0 ? away : (flux + direction * scale * atRest) / other;
  };

  auto bisection = [&](auto&& func)
  {
    // The bracket width is tested relative to the velocity scale, so the criterion no longer depends on the
    // magnitude of the flow: an absolute 1e-6 m/s was ~2e-5 relative at 0.05 m/s but a 1e-2 relative error at
    // 1e-4 m/s. The residual staircase this leaves in the right-hand side has to stay well below the
    // difference-quotient increment CVODE uses (~1e-8 relative), otherwise every J*v product is noise.
    // The absolute floor guarantees termination as the bracket approaches zero; the iteration cap is a
    // safety net, since bisection reaches the relative tolerance in roughly 40 halvings.
    constexpr double relativeTolerance = 1.0e-12;
    constexpr double absoluteFloor = 1.0e-15;
    constexpr size_t maximumIterations = 200;

    double a = 1e-7;
    double b = std::max(10.0 * std::abs(columnEntranceVelocity), 1e-6);
    double fa = func(a);
    double fb = func(b);

    for (size_t expansion = 0; fa * fb > 0.0 && expansion < 20; expansion++)
    {
      b *= 2.0;
      fb = func(b);
    }

    if (fa * fb > 0.0)
    {
      throw std::runtime_error("Bounds for bisection method to solve for velocity improperly set.\n");
    }

    double c = 0.5 * (a + b);
    for (size_t iteration = 0; iteration < maximumIterations; ++iteration)
    {
      if ((b - a) <= std::max(relativeTolerance * std::max(std::abs(a), std::abs(b)), absoluteFloor)) break;

      c = 0.5 * (a + b);
      const double fc = func(c);

      if (fc == 0.0) break;

      if (fa * fc < 0.0)
      {
        b = c;
        fb = fc;
      }
      else
      {
        a = c;
        fa = fc;
      }
    }

    return c;
  };

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }

    totalConcentration[grid] = concentrationSum;
    if (fluidPhase == Column::FluidPhase::Gas)
    {
      totalPressure[grid] = totalConcentration[grid] * R * std::max(1e-10, gasTemperature[grid]);
      gasDensity[grid] = 0.0;
    }
    else
    {
      gasDensity[grid] = liquidDensity;
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = fluidPhase == Column::FluidPhase::Gas
                                   ? concentration[index] * R * std::max(1e-10, gasTemperature[grid])
                                   : concentration[index];
      if (fluidPhase == Column::FluidPhase::Gas)
      {
        gasDensity[grid] += concentration[index] * components[comp].molecularWeight;
      }
    }

    if (pHMode == Column::PHMode::HPlus)
    {
      pH[grid] = -std::log10(std::max(1.0e-300, concentration[grid * numberOfComponents + pHComponent] / 1000.0));
    }
    else if (pHMode == Column::PHMode::OHMinus)
    {
      pH[grid] = pKw +
                 std::log10(std::max(1.0e-300, concentration[grid * numberOfComponents + pHComponent] / 1000.0));
    }
    else
    {
      pH[grid] = pHValue;
    }
  }

  if (gasTemperatureBalance != nullptr)
  {
    temperatureRateAtRest.assign(numberOfGridPoints + 1, 0.0);
    temperatureRateSlope.assign(numberOfGridPoints + 1, 0.0);
    temperatureRateReverse.assign(numberOfGridPoints + 1, 0.0);
    for (size_t grid = nonlinearFirstNode ? 2 : 1; grid < numberOfGridPoints + 1; ++grid)
    {
      temperatureRateAtRest[grid] = gasTemperatureBalance->rate(grid, 0.0, 0.0);
      temperatureRateSlope[grid] = gasTemperatureBalance->rate(grid, 0.0, 1.0) - temperatureRateAtRest[grid];
      // r_i, needed by the march from z = L; a march from z = 0 evaluates it where v_{i-1} < 0.
      if (marchFromLength && (grid < outlet || !nonlinearLastNode))
      {
        temperatureRateReverse[grid] = temperatureRateAtRest[grid] - gasTemperatureBalance->rate(grid, -1.0, 0.0);
      }
    }
    if (marchFromLength && nonlinearFirstNode)
    {
      temperatureRateReverse[1] = gasTemperatureBalance->rate(1, 0.0, 0.0) - gasTemperatureBalance->rate(1, -1.0, 0.0);
    }
  }

  // Gas at the ends: the total concentration upstream of face i for either sign of v_i, and the density of the
  // product-end gas, which the pressure drop of the last interval takes in a reverse step.
  double productEndTotal = totalConcentration[outlet];
  double productEndDensity = gasDensity[outlet];
  if (!ends.productEndConcentration.empty())
  {
    productEndTotal = 0.0;
    productEndDensity = fluidPhase == Column::FluidPhase::Gas ? 0.0 : liquidDensity;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      productEndTotal += std::max(0.0, ends.productEndConcentration[comp]);
      if (fluidPhase == Column::FluidPhase::Gas)
      {
        productEndDensity += ends.productEndConcentration[comp] * components[comp].molecularWeight;
      }
    }
  }
  auto upwindTotal = [&](size_t face)
  {
    if (interstitialGasVelocity[face] >= 0.0) return face == 0 && reverse ? totalConcentration[1] : totalConcentration[face];
    return face < outlet ? totalConcentration[face + 1] : productEndTotal;
  };
  // In a march from z = L the interval (z_{i-1}, z_i] takes the velocity of its upstream face z_i and the gas of the
  // node upstream of that face, the mirror image of ergunGrad(i - 1) in a march from z = 0.
  auto reverseErgunGrad = [&](size_t grid)
  {
    const size_t node = grid < outlet ? grid + 1 : outlet;
    const double density = grid < outlet ? gasDensity[grid + 1] : productEndDensity;
    double grad = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ads++)
    {
      grad += fractionOfAdsorbent[node * numberOfAdsorbents + ads] *
              geometries[ads].pressureDrop.gradient(dynamicViscosity, density, interstitialGasVelocity[grid]);
    }
    return grad;
  };

  // Pressure profile from a prescribed pressure at z = 0 (fromZero) or at z = L, with the interval rules of the march
  // direction. The inlet pressure applies at the feed end, the outlet pressure at the other end.
  auto integratePressure = [&](double pressure, bool fromZero)
  {
    if (fromZero)
    {
      totalPressure[0] = pressure;
      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        totalPressure[grid] = totalPressure[grid - 1] - (reverse ? reverseErgunGrad(grid) : ergunGrad(grid - 1)) *
                                                            gridSpacing(columnDistances, grid);
      }
      return;
    }
    totalPressure[numberOfGridPoints] = pressure;
    for (size_t grid = numberOfGridPoints; grid > 0; --grid)
    {
      const size_t current = grid - 1;
      totalPressure[current] = totalPressure[current + 1] + (reverse ? reverseErgunGrad(grid) : ergunGrad(current)) *
                                                                gridSpacing(columnDistances, current + 1);
    }
  };

  // Velocities of a constant superficial velocity, eps_in v_in, for an inlet velocity v_in in the bed at the feed end.
  const double inletVoidFraction = reverse ? totalVoidFraction.back() : totalVoidFraction.front();
  auto setIncompressibleVelocities = [&](double inletVelocity)
  {
    const double direction = reverse ? -1.0 : 1.0;
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      interstitialGasVelocity[grid] = direction * inletVelocity * inletVoidFraction / totalVoidFraction[grid];
    }
  };

  if (fluidPhase == Column::FluidPhase::Liquid)
  {
    if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure ||
        boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
    {
      const double targetPressure =
          boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure
              ? outletPressure
              : inletPressure + pressureGradient * columnLength;
      auto shoot = [&](double inletVelocity)
      {
        setIncompressibleVelocities(inletVelocity);
        integratePressure(inletPressure, !reverse);
        return totalPressure[reverse ? 0 : numberOfGridPoints] - targetPressure;
      };
      columnEntranceVelocity = bisection(shoot);
      shoot(columnEntranceVelocity);
    }
    else
    {
      setIncompressibleVelocities(columnEntranceVelocity);
      if (boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure ||
          (boundaryCondition == Column::BoundaryCondition::FixedVelocity && inletPressure <= 0.0))
      {
        integratePressure(outletPressure, reverse);
      }
      else
      {
        integratePressure(inletPressure, !reverse);
      }
    }

    if (totalPressure[reverse ? 0 : numberOfGridPoints] <= 0.0)
    {
      throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
    }
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      setAxialDispersion(grid, interstitialGasVelocity[grid]);
    }
    return;
  }

  // The continuity closures march from the feed, or with v = 0 from a closed end; the dispersion of a node takes the
  // velocity of the face the march reaches it through, so that the march stays explicit. A march from z = L gives the
  // velocity of the lower face of each node from the balance of that node. The inlet velocity v_in points from the feed
  // end to the other.
  const double startVelocity = ends.closedEnd ? 0.0 : reverse ? -columnEntranceVelocity : columnEntranceVelocity;
  if (ends.closedEnd || boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity ||
      boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure)
  {
    if (marchFromLength)
    {
      interstitialGasVelocity[numberOfGridPoints] = startVelocity;
      for (size_t grid = numberOfGridPoints; grid > 0; --grid)
      {
        setAxialDispersion(grid, interstitialGasVelocity[grid]);
        interstitialGasVelocity[grid - 1] = faceVelocity(
            grid,
            totalVoidFraction[grid] * interstitialGasVelocity[grid] * upwindTotal(grid) + sinkTerm(grid) -
                dispersionTerm(grid),
            totalConcentration[grid], totalConcentration[grid > 1 || !reverse ? grid - 1 : 1], true);
      }
      setAxialDispersion(0, interstitialGasVelocity[0]);
    }
    else
    {
      interstitialGasVelocity[0] = startVelocity;
      setAxialDispersion(0, interstitialGasVelocity[0]);
      for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
      {
        setAxialDispersion(grid, interstitialGasVelocity[grid - 1]);
        interstitialGasVelocity[grid] = faceVelocity(
            grid,
            totalVoidFraction[grid - 1] * interstitialGasVelocity[grid - 1] * upwindTotal(grid - 1) - sinkTerm(grid) +
                dispersionTerm(grid),
            totalConcentration[grid], grid < outlet ? totalConcentration[grid + 1] : productEndTotal, false);
      }
    }
    if (ends.closedEnd || boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity)
    {
      integratePressure(inletPressure, !reverse);
    }
    else
    {
      integratePressure(outletPressure, reverse);
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
  {
    // The march carries the molar flux with the total concentrations of the pressure profile, P / (R T), at the
    // pressure of the face it starts from.
    auto shoot = [&](double v0)
    {
      if (reverse)
      {
        interstitialGasVelocity[outlet] = -v0;
        totalPressure[outlet] = inletPressure;
        for (size_t grid = numberOfGridPoints; grid > 0; --grid)
        {
          const double beyondTemperature =
              grid < outlet ? gasTemperature[grid + 1]
                            : (ends.productEndTemperature > 0.0 ? ends.productEndTemperature : gasTemperature[grid]);
          const double cbeyond = totalPressure[grid] / (R * std::max(1e-10, beyondTemperature));
          const double cnow = totalPressure[grid] / (R * std::max(1e-10, gasTemperature[grid]));
          const double cbelow = totalPressure[grid] / (R * std::max(1e-10, gasTemperature[grid > 1 ? grid - 1 : 1]));
          setAxialDispersion(grid, interstitialGasVelocity[grid]);
          const double velocity = interstitialGasVelocity[grid];
          interstitialGasVelocity[grid - 1] = faceVelocity(
              grid,
              totalVoidFraction[grid] * velocity * (velocity < 0.0 ? cbeyond : cnow) + sinkTerm(grid) -
                  dispersionTerm(grid),
              cnow, cbelow, true);
          totalPressure[grid - 1] = totalPressure[grid] + reverseErgunGrad(grid) * gridSpacing(columnDistances, grid);
          if (totalPressure[grid - 1] <= 0.0)
          {
            return -outletPressure;
          }
        }
        setAxialDispersion(0, interstitialGasVelocity[0]);
        return totalPressure[0] - outletPressure;
      }

      interstitialGasVelocity[0] = v0;
      totalPressure[0] = inletPressure;
      setAxialDispersion(0, interstitialGasVelocity[0]);

      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        const double cprev = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid - 1]));
        const double cnow = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid]));
        const double beyondTemperature =
            grid < outlet ? gasTemperature[grid + 1]
                          : (ends.productEndTemperature > 0.0 ? ends.productEndTemperature : gasTemperature[grid]);
        const double cbeyond = totalPressure[grid - 1] / (R * std::max(1e-10, beyondTemperature));
        setAxialDispersion(grid, interstitialGasVelocity[grid - 1]);
        const double velocity = interstitialGasVelocity[grid - 1];
        interstitialGasVelocity[grid] = faceVelocity(
            grid,
            totalVoidFraction[grid - 1] * velocity * (velocity >= 0.0 ? cprev : cnow) - sinkTerm(grid) +
                dispersionTerm(grid),
            cnow, cbeyond, false);
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
        if (totalPressure[grid] <= 0.0)
        {
          return -outletPressure;
        }
      }

      return totalPressure[numberOfGridPoints] - outletPressure;
    };

    columnEntranceVelocity = bisection(shoot);
    shoot(columnEntranceVelocity);
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedVelocity)
  {
    setIncompressibleVelocities(columnEntranceVelocity);

    if (inletPressure > 0.0)
    {
      integratePressure(inletPressure, !reverse);
    }
    else
    {
      integratePressure(outletPressure, reverse);
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
  {
    // The prescribed gradient points from the feed end to the other end.
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      totalPressure[grid] =
          inletPressure + pressureGradient * (reverse ? columnLength - columnDistances[grid] : columnDistances[grid]);
      refreshNode(grid);
      interstitialGasVelocity[grid] = pressureDropLaw[grid].velocity(dynamicViscosity, gasDensity[grid],
                                                                     reverse ? -pressureGradient : pressureGradient);
    }
    columnEntranceVelocity = reverse ? -interstitialGasVelocity[outlet] : interstitialGasVelocity[0];
  }

  if (boundaryCondition == Column::BoundaryCondition::FixedVelocity ||
      boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      setAxialDispersion(grid, interstitialGasVelocity[grid]);
    }
  }

  if (totalPressure[reverse ? 0 : numberOfGridPoints] <= 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeMassDerivatives(size_t numberOfGridPoints, size_t numberOfComponents,
                            std::span<const double> columnDistances, std::span<const double> axialDispersionCoefficient,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink,
                            Column::FluidPhase fluidPhase, std::span<const double> totalVoidFraction,
                            Column::ConvectionScheme convectionScheme, Column::InletCondition inletCondition,
                            std::span<double> wenoScratch, const Column::FlowEnds& ends)
{
  const bool weno = convectionScheme == Column::ConvectionScheme::WENO3;
  const bool danckwerts = inletCondition == Column::InletCondition::Danckwerts;
  const bool reverse = ends.direction == Column::FlowDirection::Reverse;
  const size_t outlet = numberOfGridPoints;
  const std::span<const double> feed = concentration.first(numberOfComponents);

  // Danckwerts inlet: the inlet-face value c_b, with u (c_b - c_in) = g (c_1 - c_b), replaces the feed value c_0 in
  // the stencils, so that the convective and dispersive inflow u c_b + g (c_b - c_1) equals the feed flux u c_in. g is
  // D over the distance that the first dispersion stencil assigns to the inlet: h_1/2 to the WENO3 cell centre,
  // (h_1 + h_2)/2 in the upwind three-point difference (h_1 with one grid interval). In a reverse step z = 0 is the
  // outlet end: the stencils take the first cell for node 0 (zero gradient).
  std::vector<double> boundaryConcentration;
  if (danckwerts || reverse)
  {
    boundaryConcentration.assign(concentration.begin(), concentration.end());
    const double firstSpacing = gridSpacing(columnDistances, 1);
    const double inletDistance = weno                    ? 0.5 * firstSpacing
                                 : numberOfGridPoints > 1 ? 0.5 * (firstSpacing + gridSpacing(columnDistances, 2))
                                                          : firstSpacing;
    const double inflowVelocity =
        std::max(0.0, totalVoidFraction[0] * interstitialGasVelocity[0] / totalVoidFraction[1]);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      if (reverse)
      {
        boundaryConcentration[comp] = concentration[numberOfComponents + comp];
        continue;
      }
      const double conductance =
          axialDispersionCoefficient[numberOfComponents + comp] / std::max(1e-30, inletDistance);
      if (inflowVelocity + conductance > 0.0)
      {
        boundaryConcentration[comp] =
            (inflowVelocity * concentration[comp] + conductance * concentration[numberOfComponents + comp]) /
            (inflowVelocity + conductance);
      }
    }
    concentration = boundaryConcentration;
  }

  // Gas that flows in at z = L (v_N < 0) is the product-end gas, or that of the last cell; with a Danckwerts condition
  // the face value c_b, with u (c_b - c_res) = g (c_N - c_b), u = -v_N, and g = D_N over the mirror image of the inlet
  // distance. The dispersion of the last cell then takes the inflow g (c_b - c_N) across the face (outletDispersion),
  // so that the total inflow is u c_res; with a Dirichlet inlet in a reverse step it conducts towards the feed c_res.
  const double lastSpacing = gridSpacing(columnDistances, outlet);
  const double outletDistance = weno         ? 0.5 * lastSpacing
                                : outlet > 1 ? 0.5 * (lastSpacing + gridSpacing(columnDistances, outlet - 1))
                                             : lastSpacing;
  const double outletInflowVelocity = std::max(0.0, -interstitialGasVelocity[outlet]);
  const bool outletDispersion = danckwerts ? outletInflowVelocity > 0.0 : reverse;
  std::vector<double> outletFace(concentration.begin() + static_cast<std::ptrdiff_t>(outlet * numberOfComponents),
                                 concentration.end());
  if (outletInflowVelocity > 0.0 || reverse)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double last = concentration[outlet * numberOfComponents + comp];
      const double reservoir = ends.productEndConcentration.empty() ? last : ends.productEndConcentration[comp];
      if (!danckwerts)
      {
        outletFace[comp] = reservoir;
        continue;
      }
      const double conductance =
          axialDispersionCoefficient[outlet * numberOfComponents + comp] / std::max(1e-30, outletDistance);
      if (outletInflowVelocity > 0.0)
      {
        outletFace[comp] =
            (outletInflowVelocity * reservoir + conductance * last) / (outletInflowVelocity + conductance);
      }
    }
  }

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  auto axialDispersion = [&](size_t grid, size_t comp)
  { return axialDispersionCoefficient[grid * numberOfComponents + comp]; };

  auto secondDerivative = [&](size_t grid, size_t comp)
  {
    const double hm = gridSpacing(columnDistances, grid);
    const double hp = gridSpacing(columnDistances, grid + 1);
    return 2.0 *
           (hm * spanConcentration[grid + 1, comp] - (hm + hp) * spanConcentration[grid, comp] +
            hp * spanConcentration[grid - 1, comp]) /
           std::max(1e-30, hm * hp * (hm + hp));
  };

  // Last cell: zero gradient, plus the dispersion across the face z = L when gas flows in there.
  auto outletInflowDispersion = [&](size_t comp)
  {
    return (outletFace[comp] - spanConcentration[outlet, comp]) / std::max(1e-30, lastSpacing * outletDistance);
  };

  // Upwind face value of face i: node i for v_i >= 0, node i + 1, or the value at z = L, for v_i < 0.
  auto upwindFace = [&](size_t face, size_t comp)
  {
    if (interstitialGasVelocity[face] >= 0.0) return spanConcentration[face, comp];
    return face < outlet ? spanConcentration[face + 1, comp] : outletFace[comp];
  };

  // Node 0 holds the feed of a forward step, which follows a pressure ramp at its composition.
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    spanConcentrationDot[0, comp] = reverse ? 0.0 : ends.relativePressureRate * feed[comp];
  }

  if (weno)
  {
    // Finite volumes with the local spacings h_i = z_i - z_{i-1}: node i holds the cell (z_{i-1}, z_i], node 0 the
    // inlet face value. Fluxes are superficial, eps_t v c, so that they are continuous where eps_t changes (C-MBD-1).
    // Faces with v < 0 take the reconstruction from the side of z = L.
    const size_t nodes = numberOfGridPoints + 1;
    if (wenoScratch.size() < (numberOfComponents + 4) * nodes)
    {
      throw std::runtime_error("Error: the WENO3 work array of computeMassDerivatives is too small");
    }
    const std::span<double> values = wenoScratch.subspan(0, nodes);
    const std::span<double> faces = wenoScratch.subspan(nodes, nodes);
    const std::span<double> totalConcentration = wenoScratch.subspan(2 * nodes, nodes);
    const std::span<double> faceFractionSum = wenoScratch.subspan(3 * nodes, nodes);
    const std::span<double> faceFlux = wenoScratch.subspan(4 * nodes, numberOfComponents * nodes);
    const bool anyReverseFace = std::any_of(interstitialGasVelocity.begin(), interstitialGasVelocity.end(),
                                            [](double velocity) { return velocity < 0.0; });
    std::vector<double> reverseFaces(anyReverseFace ? nodes : 0);
    std::ranges::fill(totalConcentration, 0.0);
    std::ranges::fill(faceFractionSum, 0.0);
    const bool gasPhase = fluidPhase == Column::FluidPhase::Gas;
    double outletFaceTotal = 0.0;
    if (gasPhase)
    {
      for (size_t grid = 0; grid < nodes; ++grid)
      {
        for (size_t comp = 0; comp < numberOfComponents; ++comp)
        {
          totalConcentration[grid] += std::max(0.0, spanConcentration[grid, comp]);
        }
      }
      for (double value : outletFace) outletFaceTotal += std::max(0.0, value);
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      for (size_t grid = 0; grid < nodes; ++grid)
      {
        if (gasPhase)
        {
          values[grid] = totalConcentration[grid] > 0.0
                             ? std::max(0.0, spanConcentration[grid, comp]) / totalConcentration[grid]
                             : 0.0;
        }
        else
        {
          values[grid] = spanConcentration[grid, comp];
        }
      }
      computeWENO(values, columnDistances, faces);
      if (anyReverseFace)
      {
        const double inflow = !gasPhase               ? outletFace[comp]
                              : outletFaceTotal > 0.0 ? std::max(0.0, outletFace[comp]) / outletFaceTotal
                                                      : 0.0;
        computeWENOReverse(values, inflow, columnDistances, reverseFaces);
        for (size_t grid = 0; grid < nodes; ++grid)
        {
          if (interstitialGasVelocity[grid] < 0.0) faces[grid] = reverseFaces[grid];
        }
      }
      for (size_t grid = 0; grid < nodes; ++grid)
      {
        faceFlux[comp * nodes + grid] =
            gasPhase ? faces[grid] : totalVoidFraction[grid] * interstitialGasVelocity[grid] * faces[grid];
        faceFractionSum[grid] += faces[grid];
      }
    }
    if (gasPhase)
    {
      // Normalised mole fractions times the upwind total molar flux, which the velocity closures conserve.
      for (size_t grid = 0; grid < nodes; ++grid)
      {
        const double velocity = interstitialGasVelocity[grid];
        const double upwindTotal = velocity >= 0.0     ? totalConcentration[grid]
                                   : grid < outlet ? totalConcentration[grid + 1]
                                                   : outletFaceTotal;
        const double totalFlux = totalVoidFraction[grid] * velocity * upwindTotal;
        for (size_t comp = 0; comp < numberOfComponents; ++comp)
        {
          faceFlux[comp * nodes + grid] *= faceFractionSum[grid] > 0.0 ? totalFlux / faceFractionSum[grid] : 0.0;
        }
      }
    }

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double* flux = faceFlux.data() + comp * nodes;
      for (size_t grid = 1; grid <= numberOfGridPoints; ++grid)
      {
        const double hm = std::max(1e-30, gridSpacing(columnDistances, grid));
        const double fluidFraction = totalVoidFraction[grid];
        double d2cDz2 = 0.0;
        if (grid == numberOfGridPoints)
        {
          d2cDz2 = (spanConcentration[grid - 1, comp] - spanConcentration[grid, comp]) / (hm * hm);
          if (outletDispersion) d2cDz2 += outletInflowDispersion(comp);
        }
        else if (grid == 1)
        {
          // The inlet value lies on the inlet face, half a cell from the centre of the first cell.
          const double hp = gridSpacing(columnDistances, grid + 1);
          d2cDz2 = ((spanConcentration[grid + 1, comp] - spanConcentration[grid, comp]) / (0.5 * (hm + hp)) -
                    (spanConcentration[grid, comp] - spanConcentration[grid - 1, comp]) / (0.5 * hm)) /
                   hm;
        }
        else
        {
          d2cDz2 = secondDerivative(grid, comp);
        }
        spanConcentrationDot[grid, comp] = -(flux[grid] - flux[grid - 1]) / (hm * fluidFraction) +
                                           axialDispersion(grid, comp) * d2cDz2 - spanBulkSpeciesSink[grid, comp];
      }
    }
    return;
  }

  // Upwind finite volumes with superficial fluxes, as with WENO3.
  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    const double hm = gridSpacing(columnDistances, grid);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double superficialFluxGradient =
          (totalVoidFraction[grid] * interstitialGasVelocity[grid] * upwindFace(grid, comp) -
           totalVoidFraction[grid - 1] * interstitialGasVelocity[grid - 1] * upwindFace(grid - 1, comp)) /
          std::max(1e-30, hm * totalVoidFraction[grid]);
      spanConcentrationDot[grid, comp] = -superficialFluxGradient +
                                         axialDispersion(grid, comp) * secondDerivative(grid, comp) -
                                         spanBulkSpeciesSink[grid, comp];
    }
  }

  const double hm = gridSpacing(columnDistances, outlet);
  const double invHm = 1.0 / std::max(1e-30, hm);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double superficialFluxGradient =
        (totalVoidFraction[outlet] * interstitialGasVelocity[outlet] * upwindFace(outlet, comp) -
         totalVoidFraction[outlet - 1] * interstitialGasVelocity[outlet - 1] * upwindFace(outlet - 1, comp)) *
        invHm / totalVoidFraction[outlet];
    double d2cDz2 = (spanConcentration[outlet - 1, comp] - spanConcentration[outlet, comp]) * invHm * invHm;
    if (outletDispersion) d2cDz2 += outletInflowDispersion(comp);
    spanConcentrationDot[outlet, comp] =
        -superficialFluxGradient + axialDispersion(outlet, comp) * d2cDz2 - spanBulkSpeciesSink[outlet, comp];
  }
}

double GasTemperatureBalance::rate(size_t grid, double lowerVelocity, double upperVelocity) const
{
  lowerVelocity *= totalVoidFraction[grid - 1] / totalVoidFraction[grid];
  const size_t outlet = gasTemperature.size() - 1;
  const double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
  const double diffusivity = gasThermalConductivity * invGasDensity / heatCapacityGas;
  const bool danckwerts = inletCondition == Column::InletCondition::Danckwerts;
  const bool reverse = ends.direction == Column::FlowDirection::Reverse;

  // The value the stencils of node 1 take for node 0: with a Danckwerts inlet, the inlet-face temperature T_b with
  // v (T_b - T_in) = g (T_1 - T_b), g the diffusivity of node 1 over (h_1 + h_2)/2 (h_1 with one grid interval), as
  // in computeMassDerivatives. In a reverse step z = 0 is the outlet end, with zero gradient.
  double upstream = gasTemperature[grid - 1];
  if (grid == 1 && reverse)
  {
    upstream = gasTemperature[1];
  }
  else if (grid == 1 && danckwerts)
  {
    const double firstSpacing = gridSpacing(columnDistances, 1);
    const double inletDistance = outlet > 1 ? 0.5 * (firstSpacing + gridSpacing(columnDistances, 2)) : firstSpacing;
    const double inflowVelocity = std::max(0.0, upperVelocity);
    const double conductance = diffusivity / std::max(1e-30, inletDistance);
    if (inflowVelocity + conductance > 0.0)
    {
      upstream = (inflowVelocity * gasTemperature[0] + conductance * gasTemperature[1]) / (inflowVelocity + conductance);
    }
  }

  // The value beyond the last node, for gas that flows in at z = L: the product-end gas (that of the last node when
  // none is set), with a Danckwerts end its face temperature with u (T_b - T_res) = g (T_N - T_b), u = -v_{N-1} and g
  // over (h_N + h_{N-1})/2, the mirror image of node 1. The conduction across that face (outletGhost) is that of the
  // stencil of node 1: T_b at a Danckwerts end with inflow, the product-end gas at a Dirichlet inlet in a reverse step.
  double downstream = grid < outlet ? gasTemperature[grid + 1] : gasTemperature[grid];
  bool conductsAtOutlet = false;
  double outletDistance = 0.0;
  if (grid == outlet)
  {
    const double reservoir = ends.productEndTemperature > 0.0 ? ends.productEndTemperature : gasTemperature[grid];
    const double lastSpacing = gridSpacing(columnDistances, outlet);
    outletDistance = outlet > 1 ? 0.5 * (lastSpacing + gridSpacing(columnDistances, outlet - 1)) : lastSpacing;
    const double inflowVelocity = std::max(0.0, -lowerVelocity);
    if (danckwerts && inflowVelocity > 0.0)
    {
      const double conductance = diffusivity / std::max(1e-30, outletDistance);
      downstream = (inflowVelocity * reservoir + conductance * gasTemperature[grid]) / (inflowVelocity + conductance);
      conductsAtOutlet = true;
    }
    else if (!danckwerts)
    {
      downstream = reservoir;
      conductsAtOutlet = reverse;
    }
  }

  const double gasSolidExchange = fluidSolidContactArea[grid] * heatTransferGasSolid * invGasDensity / heatCapacityGas;
  const double gasWallExchange = gasWallPrefactor * invGasDensity;
  double result = gasSolidExchange * (solidTemperature[grid] - gasTemperature[grid]) +
                  gasWallExchange * (wallTemperature[grid] - gasTemperature[grid]);
  result += ends.pressureRate * invGasDensity / heatCapacityGas;  // the pressure work of a ramp

  const double hm = gridSpacing(columnDistances, grid);
  const double inverseSpacing = 1.0 / std::max(1e-30, hm);
  result -= std::max(0.0, upperVelocity) * (gasTemperature[grid] - upstream) * inverseSpacing;
  if (lowerVelocity < 0.0) result += lowerVelocity * (gasTemperature[grid] - downstream) * inverseSpacing;
  if (grid < outlet)
  {
    const double hp = gridSpacing(columnDistances, grid + 1);
    result += diffusivity * (2.0 * (hm * gasTemperature[grid + 1] - (hm + hp) * gasTemperature[grid] + hp * upstream) /
                             std::max(1e-30, hm * hp * (hm + hp)));
  }
  else
  {
    result += diffusivity * (upstream - gasTemperature[grid]) * (inverseSpacing * inverseSpacing);
    if (conductsAtOutlet)
    {
      result += diffusivity * (downstream - gasTemperature[grid]) / std::max(1e-30, hm * outletDistance);
    }
  }
  return result;
}

void computeEnergyDerivatives(size_t numberOfGridPoints, size_t numberOfComponents, double externalTemperature,
                              std::span<const Geometry> geometries, std::span<const double> particleDensity,
                              std::span<const double> solidFluidContactArea,
                              std::span<const double> fluidSolidContactArea, std::span<const double> heatOfAdsorption,
                              std::span<const double> heatOfChemisorption, double wallDensity,
                              double gasThermalConductivity, double wallThermalConductivity,
                              double heatTransferGasSolid, double heatTransferGasWall, double heatTransferWallExternal,
                              double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
                              std::span<const double> columnDistances, std::span<const double> interstitialGasVelocity,
                              std::span<const double> totalVoidFraction, std::span<const double> gasDensity,
                              std::span<double> coeffDiffusion,
                              size_t maxChemisorptionSites, std::span<const double> physisorptionDot,
                              std::span<const double> chemisorptionDot, std::span<const double> gasTemperature,
                              std::span<double> gasTemperatureDot, std::span<const double> solidTemperature,
                              std::span<double> solidTemperatureDot, std::span<const double> wallTemperature,
                              std::span<double> wallTemperatureDot, std::span<const double> reactionPhysisorptionSource,
                              std::span<const double> reactionChemisorptionSource, std::span<const double> reactionHeat,
                              Column::InletCondition inletCondition, const Column::FlowEnds& ends)
{
  mdspan2d_const spanAdsorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                      numberOfComponents);
  auto idx = [&](size_t grid) { return 1.0 / std::max(1e-30, gridSpacing(columnDistances, grid)); };
  auto d2Temperature = [&](std::span<const double> values, size_t grid)
  {
    const double hm = gridSpacing(columnDistances, grid);
    const double hp = gridSpacing(columnDistances, grid + 1);
    return 2.0 * (hm * values[grid + 1] - (hm + hp) * values[grid] + hp * values[grid - 1]) /
           std::max(1e-30, hm * hp * (hm + hp));
  };

  // fluidWallPerFluidVolume is per fluid volume (4/(eps d_in) for a packed bed), consistent with the wall balance.
  double prefactorGasWall = geometries.front().contactAreas.fluidWallPerFluidVolume *
                            heatTransferGasWall / heatCapacityGas;

  auto coeffSolidGasAt = [&](size_t grid)
  {
    return solidFluidContactArea[grid] * heatTransferGasSolid /
           (heatCapacitySolid * std::max(1e-10, particleDensity[grid]));
  };

  double internalArea = geometries.front().contactAreas.wallInnerPerWallVolume;
  double externalArea = geometries.front().contactAreas.wallOuterPerWallVolume;
  double invHeatDensityWall = 1.0 / (heatCapacityWall * wallDensity);
  double coeffWallGas = heatTransferGasWall * internalArea * invHeatDensityWall;
  double coeffWallExternal = heatTransferWallExternal * externalArea * invHeatDensityWall;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    coeffDiffusion[grid] = gasThermalConductivity * invGasDensity / heatCapacityGas;
  }

  const GasTemperatureBalance gasTemperatureBalance{
      columnDistances,       gasTemperature,         solidTemperature,     wallTemperature,
      gasDensity,            fluidSolidContactArea,  totalVoidFraction,    gasThermalConductivity,
      heatTransferGasSolid,  heatCapacityGas,        prefactorGasWall,     inletCondition,
      ends};

  auto solidHeatExchange = [&](size_t grid)
  { return coeffSolidGasAt(grid) * (gasTemperature[grid] - solidTemperature[grid]); };

  auto wallHeatExchange = [&](size_t grid)
  {
    return coeffWallGas * (gasTemperature[grid] - wallTemperature[grid]) +
           coeffWallExternal * (externalTemperature - wallTemperature[grid]);
  };

  auto physisorptionSource = [&](size_t grid, size_t comp)
  {
    const size_t index = grid * numberOfComponents + comp;
    return spanAdsorptionDot[grid, comp] -
           (reactionPhysisorptionSource.empty() ? 0.0 : reactionPhysisorptionSource[index]);
  };

  auto chemisorptionHeat = [&](size_t grid, size_t comp)
  {
    double heat = 0.0;
    const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
    for (size_t site = 0; site < maxChemisorptionSites; ++site)
    {
      const size_t index = site * componentBlockSize + grid * numberOfComponents + comp;
      const double source = spanChemisorptionDot[site, grid, comp] -
                            (reactionChemisorptionSource.empty() ? 0.0 : reactionChemisorptionSource[index]);
      heat += heatOfChemisorption[index] * source;
    }
    return heat;
  };

  auto reactionHeatAt = [&](size_t grid) { return reactionHeat.empty() ? 0.0 : reactionHeat[grid]; };

  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = solidHeatExchange(0);
  wallTemperatureDot[0] = wallHeatExchange(0);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[0] +=
        (heatOfAdsorption[0 * numberOfComponents + comp] * physisorptionSource(0, comp) + chemisorptionHeat(0, comp)) /
        heatCapacitySolid;
  }
  solidTemperatureDot[0] += reactionHeatAt(0) / heatCapacitySolid;

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    gasTemperatureDot[grid] =
        gasTemperatureBalance.rate(grid, interstitialGasVelocity[grid - 1], interstitialGasVelocity[grid]);
    solidTemperatureDot[grid] = solidHeatExchange(grid);
    wallTemperatureDot[grid] = wallHeatExchange(grid);

    wallTemperatureDot[grid] += (wallThermalConductivity * invHeatDensityWall) * d2Temperature(wallTemperature, grid);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      solidTemperatureDot[grid] +=
          (heatOfAdsorption[grid * numberOfComponents + comp] * physisorptionSource(grid, comp) +
           chemisorptionHeat(grid, comp)) /
          heatCapacitySolid;
    }
    solidTemperatureDot[grid] += reactionHeatAt(grid) / heatCapacitySolid;
  }

  gasTemperatureDot[numberOfGridPoints] = gasTemperatureBalance.rate(
      numberOfGridPoints, interstitialGasVelocity[numberOfGridPoints - 1], interstitialGasVelocity[numberOfGridPoints]);
  solidTemperatureDot[numberOfGridPoints] = solidHeatExchange(numberOfGridPoints);
  wallTemperatureDot[numberOfGridPoints] = wallHeatExchange(numberOfGridPoints);

  const double outletTemperatureIdx2 = idx(numberOfGridPoints) * idx(numberOfGridPoints);
  wallTemperatureDot[numberOfGridPoints] +=
      outletTemperatureIdx2 * (wallThermalConductivity * invHeatDensityWall) *
      (wallTemperature[numberOfGridPoints - 1] - wallTemperature[numberOfGridPoints]);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[numberOfGridPoints] += (heatOfAdsorption[numberOfGridPoints * numberOfComponents + comp] *
                                                    physisorptionSource(numberOfGridPoints, comp) +
                                                chemisorptionHeat(numberOfGridPoints, comp)) /
                                               heatCapacitySolid;
  }
  solidTemperatureDot[numberOfGridPoints] += reactionHeatAt(numberOfGridPoints) / heatCapacitySolid;
}

namespace
{
// Smoothness offset of the WENO weights; the reconstructed quantities are mole fractions (gas) or concentrations
// in mol/m^3 (liquid).
constexpr double wenoEpsilon = 1.0e-10;

/**
 * \brief WENO3 value at the downstream face of cell i from the candidates of the downstream stencil {i, i+1}
 *        (q0, linear weight d0) and the upstream stencil {i-1, i} (q1, d1); beta are the smoothness indicators.
 */
double wenoFace(double q0, double q1, double beta0, double beta1, double d0, double d1)
{
  const double alpha0 = d0 / ((wenoEpsilon + beta0) * (wenoEpsilon + beta0));
  const double alpha1 = d1 / ((wenoEpsilon + beta1) * (wenoEpsilon + beta1));
  return (alpha0 * q0 + alpha1 * q1) / (alpha0 + alpha1);
}

/**
 * \brief Geometry of the WENO3 face of cell i, 0 < i < N: the factors h_i / distance that turn the differences to the
 *        downstream and the upstream node into slopes over cell i, and the linear weights of the two stencils.
 */
struct WenoStencil
{
  double downstreamFactor;  ///< h_i over the distance from the centre of cell i to that of cell i + 1.
  double upstreamFactor;    ///< h_i over the distance to the centre of cell i - 1, or to the inlet face for i = 1.
  double downstreamWeight;  ///< Linear weight d0 of the downstream stencil {i, i+1}.
  double upstreamWeight;    ///< Linear weight d1 of the upstream stencil {i-1, i}.
};

WenoStencil wenoStencil(std::span<const double> columnDistances, size_t i, bool reverse)
{
  // Local coordinate with the centre of the cell at zero; the face is at h/2. The downstream cell has width b, the
  // upstream cell width a. For flow towards z = L the face z_i is that of cell i, with the upstream value at the inlet
  // face (-h/2) for i = 1; for flow towards z = 0 it is that of cell i + 1, the mirror image, with the upstream value
  // at the face z_N for i + 1 = N.
  const size_t last = columnDistances.size() - 1;
  const double h = reverse ? columnDistances[i + 1] - columnDistances[i] : columnDistances[i] - columnDistances[i - 1];
  const double b = reverse ? columnDistances[i] - columnDistances[i - 1] : columnDistances[i + 1] - columnDistances[i];
  const double downstreamDistance = 0.5 * (h + b);
  const double downstreamMoment = downstreamDistance * downstreamDistance + b * b / 12.0 - h * h / 12.0;
  double upstreamDistance = 0.5 * h;
  double upstreamMoment = 0.25 * h * h - h * h / 12.0;
  if (reverse ? i + 1 < last : i > 1)
  {
    const double a =
        reverse ? columnDistances[i + 2] - columnDistances[i + 1] : columnDistances[i - 1] - columnDistances[i - 2];
    upstreamDistance = 0.5 * (h + a);
    upstreamMoment = upstreamDistance * upstreamDistance + a * a / 12.0 - h * h / 12.0;
  }

  // The quadratic that reproduces the three values gives u_face = u_i + p (u_{i+1} - u_i) + m (u_i - u_{i-1});
  // p and m follow from exactness for the linear and quadratic terms.
  const double determinant = -downstreamDistance * upstreamMoment - upstreamDistance * downstreamMoment;
  const double p = (-0.5 * h * upstreamMoment - upstreamDistance * h * h / 6.0) / determinant;
  const double m = (downstreamDistance * h * h / 6.0 - 0.5 * h * downstreamMoment) / determinant;
  return {h / downstreamDistance, h / upstreamDistance, p * downstreamDistance / (0.5 * h),
          m * upstreamDistance / (0.5 * h)};
}
}  // namespace

void computeWENO(std::span<const double> values, std::span<const double> columnDistances, std::span<double> faces)
{
  const size_t last = values.size() - 1;
  if (values.size() < 3 || columnDistances.size() != values.size())
  {
    throw std::runtime_error("Unable to call WENO with numberOfGridPoints smaller than 2 or without node positions.");
  }

  faces[0] = values[0];
  for (size_t i = 1; i < last; ++i)
  {
    const WenoStencil stencil = wenoStencil(columnDistances, i, false);
    const double downstreamSlope = stencil.downstreamFactor * (values[i + 1] - values[i]);
    const double upstreamSlope = stencil.upstreamFactor * (values[i] - values[i - 1]);
    faces[i] = wenoFace(values[i] + 0.5 * downstreamSlope, values[i] + 0.5 * upstreamSlope,
                        downstreamSlope * downstreamSlope, upstreamSlope * upstreamSlope, stencil.downstreamWeight,
                        stencil.upstreamWeight);
  }
  faces[last] = values[last];
}

void computeWENOLinearization(std::span<const double> values, std::span<const double> columnDistances,
                              std::span<double> lower, std::span<double> centre, std::span<double> upper)
{
  const size_t last = values.size() - 1;
  if (values.size() < 3 || columnDistances.size() != values.size())
  {
    throw std::runtime_error("Unable to call WENO with numberOfGridPoints smaller than 2 or without node positions.");
  }

  // The inlet and outlet faces take the node values.
  for (size_t face : {size_t{0}, last})
  {
    lower[face] = 0.0;
    centre[face] = 1.0;
    upper[face] = 0.0;
  }
  for (size_t i = 1; i < last; ++i)
  {
    // face = u_i + (omega s_d + (1 - omega) s_u) / 2 with the slopes s_d = f_d (u_{i+1} - u_i) and
    // s_u = f_u (u_i - u_{i-1}), and the weight omega = alpha0 / (alpha0 + alpha1) of wenoFace, alpha_k = d_k /
    // (eps + s_k^2)^2.
    const WenoStencil stencil = wenoStencil(columnDistances, i, false);
    const double downstreamSlope = stencil.downstreamFactor * (values[i + 1] - values[i]);
    const double upstreamSlope = stencil.upstreamFactor * (values[i] - values[i - 1]);
    const double downstreamIndicator = wenoEpsilon + downstreamSlope * downstreamSlope;
    const double upstreamIndicator = wenoEpsilon + upstreamSlope * upstreamSlope;
    const double alpha0 = stencil.downstreamWeight / (downstreamIndicator * downstreamIndicator);
    const double alpha1 = stencil.upstreamWeight / (upstreamIndicator * upstreamIndicator);
    const double sum = alpha0 + alpha1;
    const double omega = alpha0 / sum;
    // d omega / d s_d and d omega / d s_u
    const double weightProduct = 2.0 * alpha0 * alpha1 / (sum * sum);
    const double downstreamWeightSlope = -2.0 * downstreamSlope * weightProduct / downstreamIndicator;
    const double upstreamWeightSlope = 2.0 * upstreamSlope * weightProduct / upstreamIndicator;
    const double halfGap = 0.5 * (downstreamSlope - upstreamSlope);
    const double faceByDownstreamSlope = 0.5 * omega + halfGap * downstreamWeightSlope;
    const double faceByUpstreamSlope = 0.5 * (1.0 - omega) + halfGap * upstreamWeightSlope;
    upper[i] = stencil.downstreamFactor * faceByDownstreamSlope;
    lower[i] = -stencil.upstreamFactor * faceByUpstreamSlope;
    centre[i] = 1.0 - upper[i] - lower[i];
  }
}

void computeWENOReverse(std::span<const double> values, double inflow, std::span<const double> columnDistances,
                        std::span<double> faces)
{
  const size_t last = values.size() - 1;
  if (values.size() < 3 || columnDistances.size() != values.size())
  {
    throw std::runtime_error("Unable to call WENO with numberOfGridPoints smaller than 2 or without node positions.");
  }

  faces[last] = inflow;
  for (size_t i = 1; i < last; ++i)
  {
    // The face z_i of cell i + 1, reconstructed as computeWENO does on the mirrored grid.
    const WenoStencil stencil = wenoStencil(columnDistances, i, true);
    const double centre = values[i + 1];
    const double upstream = i + 1 < last ? values[i + 2] : inflow;
    const double downstreamSlope = stencil.downstreamFactor * (values[i] - centre);
    const double upstreamSlope = stencil.upstreamFactor * (centre - upstream);
    faces[i] = wenoFace(centre + 0.5 * downstreamSlope, centre + 0.5 * upstreamSlope, downstreamSlope * downstreamSlope,
                        upstreamSlope * upstreamSlope, stencil.downstreamWeight, stencil.upstreamWeight);
  }
  faces[0] = values[1];
}

void computeWENOReverseLinearization(std::span<const double> values, double inflow,
                                     std::span<const double> columnDistances, std::span<double> downstream,
                                     std::span<double> centre, std::span<double> upstream)
{
  const size_t last = values.size() - 1;
  if (values.size() < 3 || columnDistances.size() != values.size())
  {
    throw std::runtime_error("Unable to call WENO with numberOfGridPoints smaller than 2 or without node positions.");
  }

  downstream[last] = 0.0;
  centre[last] = 1.0;
  upstream[last] = 0.0;
  downstream[0] = 0.0;
  centre[0] = 1.0;
  upstream[0] = 0.0;
  for (size_t i = 1; i < last; ++i)
  {
    // As computeWENOLinearization, with the downstream value u_i and the upstream value u_{i+2} of the face of cell
    // i + 1.
    const WenoStencil stencil = wenoStencil(columnDistances, i, true);
    const double downstreamSlope = stencil.downstreamFactor * (values[i] - values[i + 1]);
    const double upstreamSlope = stencil.upstreamFactor * (values[i + 1] - (i + 1 < last ? values[i + 2] : inflow));
    const double downstreamIndicator = wenoEpsilon + downstreamSlope * downstreamSlope;
    const double upstreamIndicator = wenoEpsilon + upstreamSlope * upstreamSlope;
    const double alpha0 = stencil.downstreamWeight / (downstreamIndicator * downstreamIndicator);
    const double alpha1 = stencil.upstreamWeight / (upstreamIndicator * upstreamIndicator);
    const double sum = alpha0 + alpha1;
    const double omega = alpha0 / sum;
    const double weightProduct = 2.0 * alpha0 * alpha1 / (sum * sum);
    const double downstreamWeightSlope = -2.0 * downstreamSlope * weightProduct / downstreamIndicator;
    const double upstreamWeightSlope = 2.0 * upstreamSlope * weightProduct / upstreamIndicator;
    const double halfGap = 0.5 * (downstreamSlope - upstreamSlope);
    const double faceByDownstreamSlope = 0.5 * omega + halfGap * downstreamWeightSlope;
    const double faceByUpstreamSlope = 0.5 * (1.0 - omega) + halfGap * upstreamWeightSlope;
    downstream[i] = stencil.downstreamFactor * faceByDownstreamSlope;
    upstream[i] = -stencil.upstreamFactor * faceByUpstreamSlope;
    centre[i] = 1.0 - downstream[i] - upstream[i];
  }
}

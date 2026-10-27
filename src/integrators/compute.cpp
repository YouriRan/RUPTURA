#include "compute.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mdspan>
#include <stdexcept>
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
    std::span<const double> fractionOfAdsorbent, size_t numberOfAdsorbents)
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
    if ((reaction.phase == Reaction::Phase::Chemisorbed ||
         reaction.phase == Reaction::Phase::PoreConcentration) &&
        reaction.site >= maxChemisorptionSites)
    {
      throw std::runtime_error("Error: Reaction site index exceeds available chemisorption/pore sites");
    }
    if (reaction.phase == Reaction::Phase::PoreConcentration && poreConcentration.empty())
    {
      throw std::runtime_error("Error: Pore concentration reaction requires pore concentration state");
    }
    if (reaction.phase == Reaction::Phase::BulkConcentration && (concentration.empty() || reactionBulkSource.empty()))
    {
      throw std::runtime_error("Error: Bulk concentration reaction requires the bulk concentration state");
    }

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
      double rate = sectionWeight * reaction.rate(activity, temperature);
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

      switch (reaction.phase)
      {
        case Reaction::Phase::Physisorbed:
          for (size_t i = 0; i < reaction.reactants.size(); ++i)
          {
            const size_t comp = reaction.reactants[i];
            const double stoichiometry =
                i < reaction.reactantStoichiometry.size() ? reaction.reactantStoichiometry[i] : 1.0;
            const double source = -stoichiometry * rate;
            spanPhysisorptionDot[grid, comp] += source;
            spanReactionPhysisorptionSource[grid, comp] += source;
          }
          for (size_t i = 0; i < reaction.products.size(); ++i)
          {
            const size_t comp = reaction.products[i];
            const double stoichiometry =
                i < reaction.productStoichiometry.size() ? reaction.productStoichiometry[i] : 1.0;
            const double source = stoichiometry * rate;
            spanPhysisorptionDot[grid, comp] += source;
            spanReactionPhysisorptionSource[grid, comp] += source;
          }
          break;
        case Reaction::Phase::Chemisorbed:
          for (size_t i = 0; i < reaction.reactants.size(); ++i)
          {
            const size_t comp = reaction.reactants[i];
            const double stoichiometry =
                i < reaction.reactantStoichiometry.size() ? reaction.reactantStoichiometry[i] : 1.0;
            const double source = -stoichiometry * rate;
            spanChemisorptionDot[reaction.site, grid, comp] += source;
            spanReactionChemisorptionSource[reaction.site, grid, comp] += source;
          }
          for (size_t i = 0; i < reaction.products.size(); ++i)
          {
            const size_t comp = reaction.products[i];
            const double stoichiometry =
                i < reaction.productStoichiometry.size() ? reaction.productStoichiometry[i] : 1.0;
            const double source = stoichiometry * rate;
            spanChemisorptionDot[reaction.site, grid, comp] += source;
            spanReactionChemisorptionSource[reaction.site, grid, comp] += source;
          }
          break;
        case Reaction::Phase::PoreConcentration:
          for (size_t i = 0; i < reaction.reactants.size(); ++i)
          {
            const size_t comp = reaction.reactants[i];
            const double stoichiometry =
                i < reaction.reactantStoichiometry.size() ? reaction.reactantStoichiometry[i] : 1.0;
            const double source = -stoichiometry * rate;
            spanPoreConcentrationDot[reaction.site, grid, comp] += source;
            spanReactionPoreConcentrationSource[reaction.site, grid, comp] += source;
          }
          for (size_t i = 0; i < reaction.products.size(); ++i)
          {
            const size_t comp = reaction.products[i];
            const double stoichiometry =
                i < reaction.productStoichiometry.size() ? reaction.productStoichiometry[i] : 1.0;
            const double source = stoichiometry * rate;
            spanPoreConcentrationDot[reaction.site, grid, comp] += source;
            spanReactionPoreConcentrationSource[reaction.site, grid, comp] += source;
          }
          break;
        case Reaction::Phase::BulkConcentration:
          // Applied to the bulk balance through computeBulkSpeciesSink, per unit fluid volume.
          for (size_t i = 0; i < reaction.reactants.size(); ++i)
          {
            const size_t comp = reaction.reactants[i];
            const double stoichiometry =
                i < reaction.reactantStoichiometry.size() ? reaction.reactantStoichiometry[i] : 1.0;
            spanReactionBulkSource[grid, comp] -= stoichiometry * rate;
          }
          for (size_t i = 0; i < reaction.products.size(); ++i)
          {
            const size_t comp = reaction.products[i];
            const double stoichiometry =
                i < reaction.productStoichiometry.size() ? reaction.productStoichiometry[i] : 1.0;
            spanReactionBulkSource[grid, comp] += stoichiometry * rate;
          }
          break;
      }

      // reactionHeat is per kg of solid; a bulk reaction's heat per fluid volume has no conversion here, which the
      // input reader enforces by rejecting BulkConcentration reactions together with the energy balance.
      if (reaction.phase != Reaction::Phase::BulkConcentration)
      {
        reactionHeat[grid] += -reaction.gibbsFreeEnergy * rate;
      }
    }
  }
}

void updateVelocityAndPressure(
    const std::vector<Component>& components, const Column::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const double> totalVoidFraction,
    std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature,
    Column::FluidPhase fluidPhase, double liquidDensity, Column::PHMode pHMode,
    double pHValue, double pKw, size_t pHComponent, std::span<double> pH)
{
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

  auto sinkTerm = [&](size_t grid)
  {
    const auto begin = static_cast<std::ptrdiff_t>(grid * numberOfComponents);
    const auto end = static_cast<std::ptrdiff_t>((grid + 1) * numberOfComponents);
    return gridSpacing(columnDistances, grid) *
           std::reduce(bulkSpeciesSink.begin() + begin, bulkSpeciesSink.begin() + end);
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

  if (fluidPhase == Column::FluidPhase::Liquid)
  {
    const double inletVoidFraction = totalVoidFraction.front();
    auto setLiquidVelocities = [&](double inletVelocity)
    {
      for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
      {
        interstitialGasVelocity[grid] = inletVelocity * inletVoidFraction / totalVoidFraction[grid];
      }
    };

    if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure ||
        boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
    {
      const double targetPressure =
          boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure
              ? outletPressure
              : inletPressure + pressureGradient * columnLength;
      auto shoot = [&](double inletVelocity)
      {
        setLiquidVelocities(inletVelocity);
        totalPressure[0] = inletPressure;
        for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
        {
          totalPressure[grid] = totalPressure[grid - 1] -
                                ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
        }
        return totalPressure[numberOfGridPoints] - targetPressure;
      };
      columnEntranceVelocity = bisection(shoot);
      shoot(columnEntranceVelocity);
    }
    else
    {
      setLiquidVelocities(columnEntranceVelocity);
      if (boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure ||
          (boundaryCondition == Column::BoundaryCondition::FixedVelocity && inletPressure <= 0.0))
      {
        totalPressure[numberOfGridPoints] = outletPressure;
        for (size_t grid = numberOfGridPoints; grid > 0; --grid)
        {
          const size_t current = grid - 1;
          totalPressure[current] = totalPressure[current + 1] +
                                   ergunGrad(current) * gridSpacing(columnDistances, current + 1);
        }
      }
      else
      {
        totalPressure[0] = inletPressure;
        for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
        {
          totalPressure[grid] = totalPressure[grid - 1] -
                                ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
        }
      }
    }

    if (totalPressure[numberOfGridPoints] <= 0.0)
    {
      throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
    }
    return;
  }

  if (boundaryCondition == Column::BoundaryCondition::InletPressureInletVelocity)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    totalPressure[0] = inletPressure;

    for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
    {
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
      totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
  {
    auto shoot = [&](double v0)
    {
      interstitialGasVelocity[0] = v0;
      totalPressure[0] = inletPressure;

      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        const double cprev = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid - 1]));
        const double cnow = totalPressure[grid - 1] / (R * std::max(1e-10, gasTemperature[grid]));
        interstitialGasVelocity[grid] =
            (interstitialGasVelocity[grid - 1] * cprev - sinkTerm(grid)) / std::max(1e-10, cnow);
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
  else if (boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure)
  {
    interstitialGasVelocity[0] = columnEntranceVelocity;
    totalPressure[numberOfGridPoints] = outletPressure;

    for (size_t grid = 1; grid < numberOfGridPoints + 1; grid++)
    {
      interstitialGasVelocity[grid] =
          (interstitialGasVelocity[grid - 1] * totalConcentration[grid - 1] - sinkTerm(grid)) /
          std::max(1e-10, totalConcentration[grid]);
    }
    for (size_t grid = numberOfGridPoints; grid > 0; --grid)
    {
      const size_t current = grid - 1;
      totalPressure[current] =
          totalPressure[current + 1] + ergunGrad(current) * gridSpacing(columnDistances, current + 1);
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedVelocity)
  {
    std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), columnEntranceVelocity);

    if (inletPressure > 0.0)
    {
      totalPressure[0] = inletPressure;
      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * gridSpacing(columnDistances, grid);
      }
    }
    else
    {
      totalPressure[numberOfGridPoints] = outletPressure;
      for (size_t grid = numberOfGridPoints; grid > 0; --grid)
      {
        const size_t current = grid - 1;
        totalPressure[current] =
            totalPressure[current + 1] + ergunGrad(current) * gridSpacing(columnDistances, current + 1);
      }
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      totalPressure[grid] = inletPressure + pressureGradient * columnDistances[grid];
      refreshNode(grid);

      PressureDropLaw localLaw{};
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        const size_t index = grid * numberOfAdsorbents + ads;
        localLaw.viscousCoefficient +=
            fractionOfAdsorbent[index] * geometries[ads].pressureDrop.viscousCoefficient;
        localLaw.inertialCoefficient +=
            fractionOfAdsorbent[index] * geometries[ads].pressureDrop.inertialCoefficient;
      }
      interstitialGasVelocity[grid] = localLaw.velocity(dynamicViscosity, gasDensity[grid], pressureGradient);
    }
    columnEntranceVelocity = interstitialGasVelocity[0];
  }

  if (totalPressure[numberOfGridPoints] <= 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeMassDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink,
                            Column::FluidPhase fluidPhase, std::span<const double> totalVoidFraction,
                            Column::ConvectionScheme convectionScheme)
{
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  std::vector<double> dispersion(numberOfAdsorbents * numberOfComponents);
  for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      dispersion[ads * numberOfComponents + comp] = physisorptionMixtures[ads].components[comp].axialDispersionCoefficient;
    }
  }
  auto axialDispersion = [&](size_t grid, size_t comp)
  {
    double coefficient = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      coefficient += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * dispersion[ads * numberOfComponents + comp];
    }
    return coefficient;
  };

  auto secondDerivative = [&](size_t grid, size_t comp)
  {
    const double hm = gridSpacing(columnDistances, grid);
    const double hp = gridSpacing(columnDistances, grid + 1);
    return 2.0 *
           (hm * spanConcentration[grid + 1, comp] - (hm + hp) * spanConcentration[grid, comp] +
            hp * spanConcentration[grid - 1, comp]) /
           std::max(1e-30, hm * hp * (hm + hp));
  };

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    spanConcentrationDot[0, comp] = 0.0;
  }

  if (convectionScheme == Column::ConvectionScheme::WENO3)
  {
    // Finite volumes with the local spacings h_i = z_i - z_{i-1}: node i holds the cell (z_{i-1}, z_i], node 0 the
    // inlet face value. Liquid fluxes are superficial, eps_t v c.
    const size_t nodes = numberOfGridPoints + 1;
    std::vector<double> values(nodes);
    std::vector<double> faces(nodes);
    std::vector<double> faceFlux(numberOfComponents * nodes);
    std::vector<double> totalConcentration(nodes, 0.0);
    std::vector<double> faceFractionSum(nodes, 0.0);
    const bool gasPhase = fluidPhase == Column::FluidPhase::Gas;
    if (gasPhase)
    {
      for (size_t grid = 0; grid < nodes; ++grid)
      {
        for (size_t comp = 0; comp < numberOfComponents; ++comp)
        {
          totalConcentration[grid] += std::max(0.0, spanConcentration[grid, comp]);
        }
      }
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
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        for (size_t grid = 0; grid < nodes; ++grid)
        {
          const double totalFlux = interstitialGasVelocity[grid] * totalConcentration[grid];
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
        const double fluidFraction = gasPhase ? 1.0 : totalVoidFraction[grid];
        double d2cDz2 = 0.0;
        if (grid == numberOfGridPoints)
        {
          d2cDz2 = (spanConcentration[grid - 1, comp] - spanConcentration[grid, comp]) / (hm * hm);
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

  if (fluidPhase == Column::FluidPhase::Liquid)
  {
    for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
    {
      const double hm = gridSpacing(columnDistances, grid);
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        const double superficialFluxGradient =
            (totalVoidFraction[grid] * interstitialGasVelocity[grid] * spanConcentration[grid, comp] -
             totalVoidFraction[grid - 1] * interstitialGasVelocity[grid - 1] *
                 spanConcentration[grid - 1, comp]) /
            std::max(1e-30, hm * totalVoidFraction[grid]);
        spanConcentrationDot[grid, comp] =
            -superficialFluxGradient + axialDispersion(grid, comp) * secondDerivative(grid, comp) -
            spanBulkSpeciesSink[grid, comp];
      }
    }

    const size_t outlet = numberOfGridPoints;
    const double hm = gridSpacing(columnDistances, outlet);
    const double invHm = 1.0 / std::max(1e-30, hm);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double superficialFluxGradient =
          (totalVoidFraction[outlet] * interstitialGasVelocity[outlet] * spanConcentration[outlet, comp] -
           totalVoidFraction[outlet - 1] * interstitialGasVelocity[outlet - 1] *
               spanConcentration[outlet - 1, comp]) *
          invHm / totalVoidFraction[outlet];
      const double d2cDz2 =
          (spanConcentration[outlet - 1, comp] - spanConcentration[outlet, comp]) * invHm * invHm;
      spanConcentrationDot[outlet, comp] =
          -superficialFluxGradient + axialDispersion(outlet, comp) * d2cDz2 -
          spanBulkSpeciesSink[outlet, comp];
    }
    return;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    const double hm = gridSpacing(columnDistances, grid);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double dvcDz = (interstitialGasVelocity[grid] * spanConcentration[grid, comp] -
                            interstitialGasVelocity[grid - 1] * spanConcentration[grid - 1, comp]) /
                           std::max(1e-30, hm);
      spanConcentrationDot[grid, comp] =
          -dvcDz + axialDispersion(grid, comp) * secondDerivative(grid, comp) - spanBulkSpeciesSink[grid, comp];
    }
  }

  const size_t outlet = numberOfGridPoints;
  const double hm = gridSpacing(columnDistances, outlet);
  const double invHm = 1.0 / std::max(1e-30, hm);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double dvcDz = (interstitialGasVelocity[outlet] * spanConcentration[outlet, comp] -
                          interstitialGasVelocity[outlet - 1] * spanConcentration[outlet - 1, comp]) *
                         invHm;
    const double d2cDz2 = (spanConcentration[outlet - 1, comp] - spanConcentration[outlet, comp]) * invHm * invHm;
    spanConcentrationDot[outlet, comp] =
        -dvcDz + axialDispersion(outlet, comp) * d2cDz2 - spanBulkSpeciesSink[outlet, comp];
  }
}

void computeEnergyDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                              size_t numberOfComponents, size_t numberOfAdsorbents, double externalTemperature,
                              std::span<const Geometry> geometries, std::span<const double> particleDensities,
                              std::span<const double> fractionOfAdsorbent, double wallDensity,
                              double gasThermalConductivity, double wallThermalConductivity,
                              double heatTransferGasSolid, double heatTransferGasWall, double heatTransferWallExternal,
                              double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
                              std::span<const double> columnDistances, std::span<const double> interstitialGasVelocity,
                              std::span<const double> gasDensity, std::span<double> coeffDiffusion,
                              size_t maxChemisorptionSites, std::span<const double> physisorptionDot,
                              std::span<const double> chemisorptionDot, std::span<const double> gasTemperature,
                              std::span<double> gasTemperatureDot, std::span<const double> solidTemperature,
                              std::span<double> solidTemperatureDot, std::span<const double> wallTemperature,
                              std::span<double> wallTemperatureDot, std::span<const double> reactionPhysisorptionSource,
                              std::span<const double> reactionChemisorptionSource, std::span<const double> reactionHeat)
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

  auto accessibleSurface = [&](size_t grid)
  {
    double value = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      value += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
               geometries[ads].contactAreas.solidFluidPerSolidVolume;
    }
    return value;
  };

  auto solidDensity = [&](size_t grid)
  {
    double value = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      value += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * particleDensities[ads];
    }
    return value;
  };

  auto coeffSolidGasAt = [&](size_t grid)
  {
    return accessibleSurface(grid) * heatTransferGasSolid / (heatCapacitySolid * std::max(1e-10, solidDensity(grid)));
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

  auto gasHeatExchange = [&](size_t grid)
  {
    const double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    double fluidSolidArea = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      fluidSolidArea += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
                        geometries[ads].contactAreas.fluidSolidPerFluidVolume;
    const double gasSolidExchange = fluidSolidArea * heatTransferGasSolid * invGasDensity / heatCapacityGas;
    const double gasWallExchange = prefactorGasWall * invGasDensity;
    return gasSolidExchange * (solidTemperature[grid] - gasTemperature[grid]) +
           gasWallExchange * (wallTemperature[grid] - gasTemperature[grid]);
  };

  auto solidHeatExchange = [&](size_t grid)
  { return coeffSolidGasAt(grid) * (gasTemperature[grid] - solidTemperature[grid]); };

  auto wallHeatExchange = [&](size_t grid)
  {
    return coeffWallGas * (gasTemperature[grid] - wallTemperature[grid]) +
           coeffWallExternal * (externalTemperature - wallTemperature[grid]);
  };

  auto heatOfAdsorption = [&](size_t grid, size_t comp)
  {
    double value = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      value += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
               physisorptionMixtures[ads].components[comp].heatOfAdsorption;
    }
    return value;
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
      double heatOfChemisorption = 0.0;
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        const MultiSiteChemisorption& multisite = physisorptionMixtures[ads].components[comp].chemisorption;
        if (site < multisite.numberOfSites)
        {
          heatOfChemisorption +=
              fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * multisite.sites[site].heatOfChemisorption;
        }
      }
      const size_t index = site * componentBlockSize + grid * numberOfComponents + comp;
      const double source = spanChemisorptionDot[site, grid, comp] -
                            (reactionChemisorptionSource.empty() ? 0.0 : reactionChemisorptionSource[index]);
      heat += heatOfChemisorption * source;
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
        (heatOfAdsorption(0, comp) * physisorptionSource(0, comp) + chemisorptionHeat(0, comp)) / heatCapacitySolid;
  }
  solidTemperatureDot[0] += reactionHeatAt(0) / heatCapacitySolid;

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    gasTemperatureDot[grid] = gasHeatExchange(grid);
    solidTemperatureDot[grid] = solidHeatExchange(grid);
    wallTemperatureDot[grid] = wallHeatExchange(grid);

    gasTemperatureDot[grid] -=
        interstitialGasVelocity[grid] * (gasTemperature[grid] - gasTemperature[grid - 1]) * idx(grid);
    gasTemperatureDot[grid] += coeffDiffusion[grid] * d2Temperature(gasTemperature, grid);

    wallTemperatureDot[grid] += (wallThermalConductivity * invHeatDensityWall) * d2Temperature(wallTemperature, grid);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      solidTemperatureDot[grid] +=
          (heatOfAdsorption(grid, comp) * physisorptionSource(grid, comp) + chemisorptionHeat(grid, comp)) /
          heatCapacitySolid;
    }
    solidTemperatureDot[grid] += reactionHeatAt(grid) / heatCapacitySolid;
  }

  gasTemperatureDot[numberOfGridPoints] = gasHeatExchange(numberOfGridPoints);
  solidTemperatureDot[numberOfGridPoints] = solidHeatExchange(numberOfGridPoints);
  wallTemperatureDot[numberOfGridPoints] = wallHeatExchange(numberOfGridPoints);

  gasTemperatureDot[numberOfGridPoints] -=
      interstitialGasVelocity[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints] - gasTemperature[numberOfGridPoints - 1]) * idx(numberOfGridPoints);
  const double outletTemperatureIdx2 = idx(numberOfGridPoints) * idx(numberOfGridPoints);
  gasTemperatureDot[numberOfGridPoints] +=
      coeffDiffusion[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints - 1] - gasTemperature[numberOfGridPoints]) * outletTemperatureIdx2;

  wallTemperatureDot[numberOfGridPoints] +=
      outletTemperatureIdx2 * (wallThermalConductivity * invHeatDensityWall) *
      (wallTemperature[numberOfGridPoints - 1] - wallTemperature[numberOfGridPoints]);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[numberOfGridPoints] +=
        (heatOfAdsorption(numberOfGridPoints, comp) * physisorptionSource(numberOfGridPoints, comp) +
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
    // Local coordinate with the centre of cell i at zero; the face is at h/2. The downstream cell has width b, the
    // upstream cell width a; for i = 1 the upstream value is the inlet face value at -h/2.
    const double h = columnDistances[i] - columnDistances[i - 1];
    const double b = columnDistances[i + 1] - columnDistances[i];
    const double downstreamDistance = 0.5 * (h + b);
    const double downstreamMoment = downstreamDistance * downstreamDistance + b * b / 12.0 - h * h / 12.0;
    double upstreamDistance = 0.5 * h;
    double upstreamMoment = 0.25 * h * h - h * h / 12.0;
    if (i > 1)
    {
      const double a = columnDistances[i - 1] - columnDistances[i - 2];
      upstreamDistance = 0.5 * (h + a);
      upstreamMoment = upstreamDistance * upstreamDistance + a * a / 12.0 - h * h / 12.0;
    }

    // The quadratic that reproduces the three values gives u_face = u_i + p (u_{i+1} - u_i) + m (u_i - u_{i-1});
    // p and m follow from exactness for the linear and quadratic terms.
    const double determinant = -downstreamDistance * upstreamMoment - upstreamDistance * downstreamMoment;
    const double p = (-0.5 * h * upstreamMoment - upstreamDistance * h * h / 6.0) / determinant;
    const double m = (downstreamDistance * h * h / 6.0 - 0.5 * h * downstreamMoment) / determinant;

    const double downstream = values[i + 1] - values[i];
    const double upstream = values[i] - values[i - 1];
    const double downstreamSlope = h / downstreamDistance * downstream;
    const double upstreamSlope = h / upstreamDistance * upstream;
    faces[i] = wenoFace(values[i] + 0.5 * downstreamSlope, values[i] + 0.5 * upstreamSlope,
                        downstreamSlope * downstreamSlope, upstreamSlope * upstreamSlope,
                        p * downstreamDistance / (0.5 * h), m * upstreamDistance / (0.5 * h));
  }
  faces[last] = values[last];
}

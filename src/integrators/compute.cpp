#include "compute.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mdspan>
#include <stdexcept>
#include <vector>

#include "transport.h"
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
    std::span<double> reactionPoreConcentrationSource, std::span<double> reactionHeat)
{
  std::fill(reactionPhysisorptionSource.begin(), reactionPhysisorptionSource.end(), 0.0);
  std::fill(reactionChemisorptionSource.begin(), reactionChemisorptionSource.end(), 0.0);
  std::fill(reactionPoreConcentrationSource.begin(), reactionPoreConcentrationSource.end(), 0.0);
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
        }
      }

      const double temperature = solidTemperature.empty() ? externalTemperature : solidTemperature[grid];
      double rate = reaction.rate(activity, temperature);
      if (reaction.phase == Reaction::Phase::PoreConcentration)
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
      }

      reactionHeat[grid] += -reaction.gibbsFreeEnergy * rate;
    }
  }
}

void updateVelocityAndPressure(const std::vector<Component>& components,
                               const Column::BoundaryCondition& boundaryCondition, size_t numberOfGridPoints,
                               size_t numberOfComponents, double inletPressure, double outletPressure,
                               double pressureGradient, double columnLength, const Geometry& geometry,
                               double& columnEntranceVelocity, double dynamicViscosity, double resolution,
                               std::span<double> interstitialGasVelocity, std::span<double> gasDensity,
                               std::span<double> totalConcentration, std::span<double> totalPressure,
                               std::span<const double> concentration, std::span<double> partialPressure,
                               std::span<double> moleFraction, std::span<const double> bulkSpeciesSink,
                               std::span<const double> gasTemperature, Column::FluidPhase fluidPhase,
                               double liquidDensity, Column::PHMode pHMode, double pHValue, double pKw,
                               size_t pHComponent, std::span<double> pH)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; grid++)
  {
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; comp++)
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
    for (size_t comp = 0; comp < numberOfComponents; comp++)
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

  auto ergunGrad = [&](size_t grid)
  {
    return geometry.pressureDrop.gradient(dynamicViscosity, gasDensity[grid], interstitialGasVelocity[grid]);
  };

  auto sinkTerm = [&](size_t grid)
  {
    const auto begin = static_cast<std::ptrdiff_t>(grid * numberOfComponents);
    const auto end = static_cast<std::ptrdiff_t>((grid + 1) * numberOfComponents);
    return resolution * std::reduce(bulkSpeciesSink.begin() + begin, bulkSpeciesSink.begin() + end);
  };

  auto gridRatio = [&](size_t grid) -> double
  { return numberOfGridPoints == 0 ? 0.0 : static_cast<double>(grid) / static_cast<double>(numberOfGridPoints); };

  auto refreshNode = [&](size_t grid)
  {
    totalConcentration[grid] = totalPressure[grid] / (R * std::max(1e-10, gasTemperature[grid]));
    gasDensity[grid] = 0.0;
    double concentrationSum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      concentrationSum += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      moleFraction[index] = concentrationSum > 0.0 ? std::max(0.0, concentration[index]) / concentrationSum : 0.0;
      partialPressure[index] = moleFraction[index] * totalPressure[grid];
      gasDensity[grid] += moleFraction[index] * totalConcentration[grid] * components[comp].molecularWeight;
    }
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

  if (fluidPhase == Column::FluidPhase::Liquid)
  {
    if (boundaryCondition == Column::BoundaryCondition::InletPressureOutletPressure)
    {
      columnEntranceVelocity = geometry.pressureDrop.velocity(
          dynamicViscosity, liquidDensity, (outletPressure - inletPressure) / columnLength);
    }

    std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), columnEntranceVelocity);

    if (boundaryCondition == Column::BoundaryCondition::InletVelocityOutletPressure ||
        (boundaryCondition == Column::BoundaryCondition::FixedVelocity && inletPressure <= 0.0))
    {
      totalPressure[numberOfGridPoints] = outletPressure;
      for (size_t grid = numberOfGridPoints; grid > 0; --grid)
      {
        const size_t current = grid - 1;
        totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * resolution;
      }
    }
    else if (boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
    {
      for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
      {
        totalPressure[grid] = inletPressure + pressureGradient * columnLength * gridRatio(grid);
      }
      columnEntranceVelocity = geometry.pressureDrop.velocity(dynamicViscosity, liquidDensity, pressureGradient);
      std::fill(interstitialGasVelocity.begin(), interstitialGasVelocity.end(), columnEntranceVelocity);
    }
    else
    {
      totalPressure[0] = inletPressure;
      for (size_t grid = 1; grid < numberOfGridPoints + 1; ++grid)
      {
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
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
      totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
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
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
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
      totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * resolution;
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
        totalPressure[grid] = totalPressure[grid - 1] - ergunGrad(grid - 1) * resolution;
      }
    }
    else
    {
      totalPressure[numberOfGridPoints] = outletPressure;
      for (size_t grid = numberOfGridPoints; grid > 0; --grid)
      {
        const size_t current = grid - 1;
        totalPressure[current] = totalPressure[current + 1] + ergunGrad(current) * resolution;
      }
    }
  }
  else if (boundaryCondition == Column::BoundaryCondition::FixedPressureInletVelocity)
  {
    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      const double position = columnLength * gridRatio(grid);
      totalPressure[grid] = inletPressure + pressureGradient * position;
      refreshNode(grid);
      interstitialGasVelocity[grid] =
          geometry.pressureDrop.velocity(dynamicViscosity, gasDensity[grid], pressureGradient);
    }
    columnEntranceVelocity = interstitialGasVelocity[0];
  }

  if (totalPressure[numberOfGridPoints] <= 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeMassDerivatives(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, double resolution,
                            std::span<const double> interstitialGasVelocity,
                            std::span<const double> concentration, std::span<double> concentrationDot,
                            std::span<const double> bulkSpeciesSink)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    spanConcentrationDot[0, comp] = 0.0;
  }

  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const double dvcDz = (interstitialGasVelocity[grid] * spanConcentration[grid, comp] -
                           interstitialGasVelocity[grid - 1] * spanConcentration[grid - 1, comp]) *
                           idx;
      const double d2cDz2 =
          (spanConcentration[grid + 1, comp] - 2.0 * spanConcentration[grid, comp] +
           spanConcentration[grid - 1, comp]) *
          idx2;

      spanConcentrationDot[grid, comp] = -dvcDz + components[comp].axialDispersionCoefficient * d2cDz2 -
                                         spanBulkSpeciesSink[grid, comp];
    }
  }

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double dvcDz = (interstitialGasVelocity[numberOfGridPoints] *
                              spanConcentration[numberOfGridPoints, comp] -
                         interstitialGasVelocity[numberOfGridPoints - 1] *
                             spanConcentration[numberOfGridPoints - 1, comp]) *
                         idx;
    const double d2cDz2 =
        (spanConcentration[numberOfGridPoints - 1, comp] - spanConcentration[numberOfGridPoints, comp]) * idx2;

    spanConcentrationDot[numberOfGridPoints, comp] =
        -dvcDz + components[comp].axialDispersionCoefficient * d2cDz2 -
        spanBulkSpeciesSink[numberOfGridPoints, comp];
  }
}

void computeEnergyDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, const Geometry& geometry,
    double particleDensity, double wallDensity, double gasThermalConductivity,
    double wallThermalConductivity, double heatTransferGasSolid, double heatTransferGasWall,
    double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid, double heatCapacityWall,
    double resolution, std::span<const double> interstitialGasVelocity, std::span<const double> gasDensity,
    std::span<double> coeffDiffusion, std::span<const double> physisorptionDot,
    std::span<const double> chemisorptionDot,
    std::span<const double> gasTemperature, std::span<double> gasTemperatureDot,
    std::span<const double> solidTemperature, std::span<double> solidTemperatureDot,
    std::span<const double> wallTemperature, std::span<double> wallTemperatureDot,
    std::span<const double> reactionPhysisorptionSource,
    std::span<const double> reactionChemisorptionSource,
    std::span<const double> reactionHeat)
{
  mdspan2d_const spanPhysisorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  std::mdspan<const double, std::dextents<size_t, 3>> spanChemisorptionDot(
      chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1, numberOfComponents);
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  // commented out the parts for weno, seems to be unstable
  // std::vector<double> gasTemperatureFlux(numberOfGridPoints + 1);
  // computeWENO(gasTemperature, gasTemperatureFlux);

  const double prefactorGasSolid =
      geometry.contactAreas.fluidSolidPerFluidVolume * heatTransferGasSolid / heatCapacityGas;

  // is this extra 1/eps necessary? It's in python not in eqs
  const double prefactorGasWall =
      geometry.contactAreas.fluidWallPerFluidVolume * heatTransferGasWall / heatCapacityGas;

  const double coeffSolidGas =
      geometry.contactAreas.solidFluidPerSolidVolume * heatTransferGasSolid /
      (heatCapacitySolid * particleDensity);

  // prefactor 4 in python, 2 in eqs
  double invHeatDensityWall = 1.0 / (heatCapacityWall * wallDensity);
  const double coeffWallGas =
      heatTransferGasWall * geometry.contactAreas.wallInnerPerWallVolume * invHeatDensityWall;
  const double coeffWallExternal =
      heatTransferWallExternal * geometry.contactAreas.wallOuterPerWallVolume * invHeatDensityWall;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    coeffDiffusion[grid] = gasThermalConductivity * invGasDensity / heatCapacityGas;
  }

  auto gasHeatExchange = [&](size_t grid)
  {
    const double invGasDensity = 1.0 / std::max(1e-10, gasDensity[grid]);
    const double gasSolidExchange = prefactorGasSolid * invGasDensity;
    const double gasWallExchange = prefactorGasWall * invGasDensity;
    return gasSolidExchange * (solidTemperature[grid] - gasTemperature[grid]) +
           gasWallExchange * (wallTemperature[grid] - gasTemperature[grid]);
  };

  auto solidHeatExchange = [&](size_t grid)
  {
    return coeffSolidGas * (gasTemperature[grid] - solidTemperature[grid]);
  };

  auto wallHeatExchange = [&](size_t grid)
  {
    return coeffWallGas * (gasTemperature[grid] - wallTemperature[grid]) +
           coeffWallExternal * (externalTemperature - wallTemperature[grid]);
  };

  auto chemisorptionHeat = [&](size_t grid, size_t comp)
  {
    double heat = 0.0;
    const MultiSiteChemisorption& multisite = components[comp].chemisorption;
    for (size_t site = 0; site < multisite.numberOfSites; ++site)
    {
      double source = spanChemisorptionDot[site, grid, comp];
      if (!reactionChemisorptionSource.empty())
      {
        const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
        source -= reactionChemisorptionSource[site * componentBlockSize + grid * numberOfComponents + comp];
      }
      heat += multisite.sites[site].heatOfChemisorption * source;
    }
    return heat;
  };

  auto physisorptionSource = [&](size_t grid, size_t comp)
  {
    double source = spanPhysisorptionDot[grid, comp];
    if (!reactionPhysisorptionSource.empty())
    {
      source -= reactionPhysisorptionSource[grid * numberOfComponents + comp];
    }
    return source;
  };

  auto reactionHeatAt = [&](size_t grid)
  {
    return reactionHeat.empty() ? 0.0 : reactionHeat[grid];
  };

  // first grid point
  gasTemperatureDot[0] = 0.0;
  solidTemperatureDot[0] = solidHeatExchange(0);
  wallTemperatureDot[0] = wallHeatExchange(0);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[0] +=
        (components[comp].heatOfAdsorption * physisorptionSource(0, comp) +
         chemisorptionHeat(0, comp)) /
        heatCapacitySolid;
  }
  solidTemperatureDot[0] += reactionHeatAt(0) / heatCapacitySolid;

  // middle grid points
  for (size_t grid = 1; grid < numberOfGridPoints; ++grid)
  {
    // heat flux transfer
    gasTemperatureDot[grid] = gasHeatExchange(grid);
    solidTemperatureDot[grid] = solidHeatExchange(grid);
    wallTemperatureDot[grid] = wallHeatExchange(grid);

    // flux from gas diffusion and advection
    gasTemperatureDot[grid] -= interstitialGasVelocity[grid] * (gasTemperature[grid] - gasTemperature[grid - 1]) * idx;
    gasTemperatureDot[grid] += coeffDiffusion[grid] *
                               (gasTemperature[grid - 1] - 2.0 * gasTemperature[grid] + gasTemperature[grid + 1]) *
                               idx2;

    // flux from heat diffusion in wall
    wallTemperatureDot[grid] += idx2 * (wallThermalConductivity * invHeatDensityWall) *
                                (wallTemperature[grid - 1] - 2.0 * wallTemperature[grid] + wallTemperature[grid + 1]);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      solidTemperatureDot[grid] +=
          (components[comp].heatOfAdsorption * physisorptionSource(grid, comp) +
           chemisorptionHeat(grid, comp)) /
          heatCapacitySolid;
    }
    solidTemperatureDot[grid] += reactionHeatAt(grid) / heatCapacitySolid;
  }

  // last gridpoint
  // heat flux transfer
  gasTemperatureDot[numberOfGridPoints] = gasHeatExchange(numberOfGridPoints);
  solidTemperatureDot[numberOfGridPoints] = solidHeatExchange(numberOfGridPoints);
  wallTemperatureDot[numberOfGridPoints] = wallHeatExchange(numberOfGridPoints);

  // flux from gas diffusion and advection
  gasTemperatureDot[numberOfGridPoints] -=
      interstitialGasVelocity[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints] - gasTemperature[numberOfGridPoints - 1]) * idx;
  gasTemperatureDot[numberOfGridPoints] +=
      coeffDiffusion[numberOfGridPoints] *
      (gasTemperature[numberOfGridPoints - 1] - gasTemperature[numberOfGridPoints]) * idx2;

  // flux from heat diffusion in wall
  wallTemperatureDot[numberOfGridPoints] +=
      idx2 * (wallThermalConductivity * invHeatDensityWall) *
      (wallTemperature[numberOfGridPoints - 1] - wallTemperature[numberOfGridPoints]);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    solidTemperatureDot[numberOfGridPoints] +=
        (components[comp].heatOfAdsorption * physisorptionSource(numberOfGridPoints, comp) +
         chemisorptionHeat(numberOfGridPoints, comp)) /
        heatCapacitySolid;
  }
  solidTemperatureDot[numberOfGridPoints] += reactionHeatAt(numberOfGridPoints) / heatCapacitySolid;
}

void computeWENO(std::span<const double> input, std::span<double> output)
{
  double tol = 1e-10;
  double df0, df1, alpha_0, alpha_1, beta_0, beta_1, first_term, second_term;

  size_t size = input.size();
  if (size < 3)
  {
    throw std::runtime_error("Unable to call WENO with numberOfGridPoints smaller than 3.");
  }

  // inlet boundary flux: prescribed from Dirichlet inflow state
  output[0] = input[0];

  // first interior interface, special one-sided closure
  df0 = input[2] - input[1];
  df1 = input[1] - input[0];
  beta_0 = df0 * df0;
  beta_1 = df1 * df1;

  alpha_0 = (2.0 / 3.0) / ((beta_0 + tol) * (beta_0 + tol));
  alpha_1 = (1.0 / 3.0) / (16.0 * (beta_1 + tol) * (beta_1 + tol));

  first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (input[2] + input[1]);
  second_term = (alpha_1 / (alpha_0 + alpha_1)) * (2.0 * input[1] - input[0]);
  output[1] = first_term + second_term;

  // interior interfaces
  for (size_t i = 2; i < size - 1; ++i)
  {
    df0 = input[i + 1] - input[i];
    df1 = input[i] - input[i - 1];
    beta_0 = df0 * df0;
    beta_1 = df1 * df1;

    alpha_0 = (2.0 / 3.0) / ((beta_0 + tol) * (beta_0 + tol));
    alpha_1 = (1.0 / 3.0) / ((beta_1 + tol) * (beta_1 + tol));

    first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (input[i + 1] + input[i]);
    second_term = (alpha_1 / (alpha_0 + alpha_1)) * ((3.0 / 2.0) * input[i] - (1.0 / 2.0) * input[i - 1]);
    output[i] = first_term + second_term;
  }

  // outlet boundary flux: outflow closure
  output[size - 1] = input[size - 1];  // simplest option
}

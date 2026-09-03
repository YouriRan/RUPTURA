#include "compute_multibed.h"

#include <algorithm>
#include <cmath>
#include <mdspan>

#include "utils.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

void updateVelocityAndPressure(
    const std::vector<Component>& components, const MultibedColumn::BoundaryCondition& boundaryCondition,
    size_t numberOfGridPoints, size_t numberOfComponents, double inletPressure, double outletPressure,
    double pressureGradient, double columnLength, size_t numberOfAdsorbents, double& columnEntranceVelocity,
    double dynamicViscosity, std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const double> adsorbentScaledVoidFraction,
    std::span<const double> totalVoidFraction,
    std::span<double> interstitialGasVelocity,
    std::span<double> gasDensity, std::span<double> totalConcentration, std::span<double> totalPressure,
    std::span<const double> concentration, std::span<double> partialPressure, std::span<double> moleFraction,
    std::span<const double> bulkSpeciesSink, std::span<const double> gasTemperature,
    MultibedColumn::FluidPhase fluidPhase, double liquidDensity, MultibedColumn::PHMode pHMode,
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
    if (fluidPhase == MultibedColumn::FluidPhase::Gas)
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
      partialPressure[index] = fluidPhase == MultibedColumn::FluidPhase::Gas
                                   ? concentration[index] * R * std::max(1e-10, gasTemperature[grid])
                                   : concentration[index];
      if (fluidPhase == MultibedColumn::FluidPhase::Gas)
      {
        gasDensity[grid] += concentration[index] * components[comp].molecularWeight;
      }
    }

    if (pHMode == MultibedColumn::PHMode::HPlus)
    {
      pH[grid] = -std::log10(std::max(1.0e-300, concentration[grid * numberOfComponents + pHComponent] / 1000.0));
    }
    else if (pHMode == MultibedColumn::PHMode::OHMinus)
    {
      pH[grid] = pKw +
                 std::log10(std::max(1.0e-300, concentration[grid * numberOfComponents + pHComponent] / 1000.0));
    }
    else
    {
      pH[grid] = pHValue;
    }
  }

  if (fluidPhase == MultibedColumn::FluidPhase::Liquid)
  {
    const double inletVoidFraction = totalVoidFraction.front();
    auto setLiquidVelocities = [&](double inletVelocity)
    {
      for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
      {
        interstitialGasVelocity[grid] = inletVelocity * inletVoidFraction / totalVoidFraction[grid];
      }
    };

    if (boundaryCondition == MultibedColumn::BoundaryCondition::InletPressureOutletPressure ||
        boundaryCondition == MultibedColumn::BoundaryCondition::FixedPressureInletVelocity)
    {
      const double targetPressure =
          boundaryCondition == MultibedColumn::BoundaryCondition::InletPressureOutletPressure
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
      if (boundaryCondition == MultibedColumn::BoundaryCondition::InletVelocityOutletPressure ||
          (boundaryCondition == MultibedColumn::BoundaryCondition::FixedVelocity && inletPressure <= 0.0))
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

  if (boundaryCondition == MultibedColumn::BoundaryCondition::InletPressureInletVelocity)
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
  else if (boundaryCondition == MultibedColumn::BoundaryCondition::InletPressureOutletPressure)
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
  else if (boundaryCondition == MultibedColumn::BoundaryCondition::InletVelocityOutletPressure)
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
  else if (boundaryCondition == MultibedColumn::BoundaryCondition::FixedVelocity)
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
  else if (boundaryCondition == MultibedColumn::BoundaryCondition::FixedPressureInletVelocity)
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

void computePhysisorptionEquilibriumLoadings(
    std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature,
    MixturePrediction::DrivingForceInput input, std::span<const double> concentration,
    std::span<const double> pH)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(
          0.0, input == MixturePrediction::DrivingForceInput::MoleFraction
                   ? moleFraction[grid * numberOfComponents + comp]
                   : concentration[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
    }
    if (input == MixturePrediction::DrivingForceInput::MoleFraction && sum > 0.0)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        idealGasMolFractions[comp] /= sum;
      }
    }
    else if (input == MixturePrediction::DrivingForceInput::MoleFraction)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        idealGasMolFractions[comp] = 1.0 / static_cast<double>(numberOfComponents);
      }
    }

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumPhysisorption[grid * numberOfComponents + comp] = 0.0;
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (!hasAdsorbentOfType[grid * numberOfAdsorbents + ads]) continue;

      // use Yi and Pt[i] to compute the loadings in the physisorption mixture via mixture prediction
      std::span<double> spanCachedPressure =
          cachedPressure.subspan((grid * numberOfAdsorbents + ads) * numberOfComponents * maxIsothermTerms,
                                 numberOfComponents * maxIsothermTerms);
      std::span<double> spanCachedGrandPotential =
          cachedGrandPotential.subspan((grid * numberOfAdsorbents + ads) * maxIsothermTerms, maxIsothermTerms);
      iastPerformance += physisorptionMixtures[ads].predictMixture(
          idealGasMolFractions,
          input == MixturePrediction::DrivingForceInput::MoleFraction ? totalPressure[grid] : 1.0,
          adsorbedMolFractions, numberOfMolecules, spanCachedPressure, spanCachedGrandPotential,
          gasTemperature[grid], pH[grid], input);

      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        equilibriumPhysisorption[grid * numberOfComponents + comp] +=
            fractionOfAdsorbent[grid * numberOfAdsorbents + ads] * numberOfMolecules[comp];
      }
    }
  }
}

void computeChemisorptionEquilibriumLoadings(
    std::vector<MixturePrediction>& chemisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxChemisorptionSites, std::pair<size_t, size_t>& iastPerformance, std::span<double> idealGasMolFractions,
    std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumChemisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction, std::span<double> gasTemperature,
    MixturePrediction::DrivingForceInput input, std::span<const double> concentration,
    std::span<const double> pH)
{
  std::fill(equilibriumChemisorption.begin(), equilibriumChemisorption.end(), 0.0);
  const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(
          0.0, input == MixturePrediction::DrivingForceInput::MoleFraction
                   ? moleFraction[grid * numberOfComponents + comp]
                   : concentration[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      if (input == MixturePrediction::DrivingForceInput::MoleFraction)
      {
        idealGasMolFractions[comp] =
            sum > 0.0 ? idealGasMolFractions[comp] / sum : 1.0 / static_cast<double>(numberOfComponents);
      }
    }

    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (!hasAdsorbentOfType[grid * numberOfAdsorbents + ads]) continue;

      MixturePrediction& mixture = chemisorptionMixtures[ads];
      std::span<double> spanCachedPressure =
          cachedPressure.subspan((grid * numberOfAdsorbents + ads) * numberOfComponents * maxChemisorptionSites,
                                 numberOfComponents * maxChemisorptionSites);
      std::span<double> spanCachedGrandPotential = cachedGrandPotential.subspan(
          (grid * numberOfAdsorbents + ads) * maxChemisorptionSites, maxChemisorptionSites);
      iastPerformance += mixture.predictMixture(
          idealGasMolFractions,
          input == MixturePrediction::DrivingForceInput::MoleFraction ? totalPressure[grid] : 1.0,
          adsorbedMolFractions, numberOfMolecules, spanCachedPressure, spanCachedGrandPotential,
          gasTemperature[grid], pH[grid], input);

      const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      for (size_t site = 0; site < mixture.maxIsothermTerms; ++site)
      {
        for (size_t comp = 0; comp < numberOfComponents; ++comp)
        {
          equilibriumChemisorption[site * componentBlockSize + grid * numberOfComponents + comp] +=
              fraction * mixture.equilibriumSiteLoadings[site * numberOfComponents + comp];
        }
      }
    }
  }
}

void computePhysisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents,
                          std::span<const double> fractionOfAdsorbent, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      double massTransferCoefficient = 0.0;
      for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
      {
        massTransferCoefficient += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
                                   physisorptionMixtures[ads].components[comp].massTransferCoefficient;
      }

      const size_t index = grid * numberOfComponents + comp;
      physisorptionDot[index] = massTransferCoefficient * (equilibriumPhysisorption[index] - physisorption[index]);
    }
  }
}

void computeChemisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                          double externalTemperature, double elapsedTime, std::span<const double> fractionOfAdsorbent,
                          std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
                          std::span<const double> equilibriumChemisorption, std::span<const double> concentration,
                          std::span<const double> chemisorption, std::span<double> chemisorptionDot,
                          std::span<const double> poreConcentration, std::span<const double> solidTemperature)
{
  std::fill(chemisorptionDot.begin(), chemisorptionDot.end(), 0.0);
  mdspan3d_const spanEquilibrium(equilibriumChemisorption.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                 numberOfComponents);
  mdspan3d_const spanChemisorption(chemisorption.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                   numberOfComponents);
  mdspan3d_mut spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                    numberOfComponents);
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanPoreConcentration(poreConcentration.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                       numberOfComponents);
  constexpr double poreInventoryLimitTime = 1.0e-4;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    const double temperature = solidTemperature.empty() ? externalTemperature : solidTemperature[grid];
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      for (size_t site = 0; site < maxChemisorptionSites; ++site)
      {
        double blendedRate = 0.0;
        for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
        {
          const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
          if (fraction <= 0.0) continue;

          const MultiSiteChemisorption& multisite = physisorptionMixtures[ads].components[comp].chemisorption;
          if (site >= multisite.numberOfSites) continue;
          const Chemisorption& kinetics = multisite.sites[site];
          if (!kinetics.enabled()) continue;

          const double drivingConcentration = kinetics.usesSurfacePoreTransport()
                                                  ? spanPoreConcentration[site, grid, comp]
                                                  : spanConcentration[grid, comp];
          double equilibriumLoading = std::max(0.0, spanEquilibrium[site, grid, comp]);
          if (kinetics.maximumLoading > 0.0)
          {
            equilibriumLoading = std::min(equilibriumLoading, kinetics.maximumLoading);
          }
          double rate = kinetics.rate(equilibriumLoading, spanChemisorption[site, grid, comp], drivingConcentration,
                                      temperature, elapsedTime);
          if (kinetics.usesSurfacePoreTransport())
          {
            const double gamma = (1.0 - adsorbentVoidFractions[ads]) * particleDensities[ads] /
                                 std::max(1.0e-10, adsorbentVoidFractions[ads]);
            if (rate > 0.0 && gamma > 0.0)
            {
              rate = std::min(rate, drivingConcentration / (gamma * poreInventoryLimitTime));
            }
          }

          if (spanChemisorption[site, grid, comp] >= equilibriumLoading)
          {
            rate = std::min(0.0, rate);
          }
          else if (spanChemisorption[site, grid, comp] <= 0.0)
          {
            rate = std::max(0.0, rate);
          }
          blendedRate += fraction * rate;
        }
        spanChemisorptionDot[site, grid, comp] = blendedRate;
      }
    }
  }
}

void computeChemisorptionTransportDerivatives(
    const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, size_t maxChemisorptionSites, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries,
    std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
    std::span<const double> particleDiameters, std::span<const double> totalVoidFraction,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot)
{
  std::fill(surfaceConcentrationDot.begin(), surfaceConcentrationDot.end(), 0.0);
  std::fill(poreConcentrationDot.begin(), poreConcentrationDot.end(), 0.0);
  if (surfaceConcentration.empty() || poreConcentration.empty()) return;

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                      numberOfComponents);
  mdspan3d_const spanSurface(surfaceConcentration.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                             numberOfComponents);
  mdspan3d_mut spanSurfaceDot(surfaceConcentrationDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                              numberOfComponents);
  mdspan3d_const spanPore(poreConcentration.data(), maxChemisorptionSites, numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanPoreDot(poreConcentrationDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                           numberOfComponents);

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double solidLoadingDensity = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      solidLoadingDensity += fraction * geometries[ads].solidToFluidVolumeRatio *
                             geometries[ads].voidFraction * particleDensities[ads];
    }
    const double loadingPrefactor = solidLoadingDensity / std::max(1.0e-10, totalVoidFraction[grid]);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      for (size_t site = 0; site < maxChemisorptionSites; ++site)
      {
        double film = 0.0;
        double poreDiffusion = 0.0;
        for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
        {
          const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
          const MultiSiteChemisorption& multisite = physisorptionMixtures[ads].components[comp].chemisorption;
          if (fraction <= 0.0 || site >= multisite.numberOfSites) continue;
          const Chemisorption& kinetics = multisite.sites[site];
          if (!kinetics.usesSurfacePoreTransport()) continue;

          const double solidContactArea = geometries[ads].contactAreas.solidFluidPerSolidVolume;
          const double poreDiffusionLength =
              std::max(1.0e-30, geometries[ads].dimensions.poreDiffusionLength);
          const double poreLdf = 15.0 * std::max(0.0, kinetics.poreDiffusivity) /
                                 std::max(1.0e-30, poreDiffusionLength * poreDiffusionLength);
          film += fraction * solidContactArea * std::max(0.0, kinetics.filmMassTransferCoefficient) *
                  (spanConcentration[grid, comp] - spanSurface[site, grid, comp]);
          poreDiffusion += fraction * poreLdf * (spanSurface[site, grid, comp] - spanPore[site, grid, comp]);
        }
        spanSurfaceDot[site, grid, comp] = film - poreDiffusion;
        spanPoreDot[site, grid, comp] = poreDiffusion - loadingPrefactor * spanChemisorptionDot[site, grid, comp];
      }
    }
  }
}

void computeBulkSpeciesSink(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                            std::span<const Geometry> geometries,
                            std::span<const double> adsorbentVoidFractions, std::span<const double> particleDensities,
                            std::span<const double> particleDiameters, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> totalVoidFraction, std::span<const double> concentration,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration, std::span<double> bulkSpeciesSink,
                            std::span<const double> reactionPhysisorptionSource,
                            std::span<const double> reactionChemisorptionSource)
{
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                      numberOfComponents);
  mdspan3d_const spanSurface(surfaceConcentration.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                             numberOfComponents);

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double solidLoadingDensity = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      solidLoadingDensity += fraction * geometries[ads].solidToFluidVolumeRatio *
                             geometries[ads].voidFraction * particleDensities[ads];
    }

    const double loadingPrefactor = solidLoadingDensity / std::max(1e-10, totalVoidFraction[grid]);
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      const size_t index = grid * numberOfComponents + comp;
      const double reactionPhysisorption =
          reactionPhysisorptionSource.empty() ? 0.0 : reactionPhysisorptionSource[index];
      double sink = loadingPrefactor * (physisorptionDot[index] - reactionPhysisorption);

      for (size_t site = 0; site < maxChemisorptionSites; ++site)
      {
        double directFraction = 0.0;
        double filmSink = 0.0;
        for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
        {
          const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
          const MultiSiteChemisorption& multisite = physisorptionMixtures[ads].components[comp].chemisorption;
          if (fraction <= 0.0 || site >= multisite.numberOfSites) continue;
          const Chemisorption& kinetics = multisite.sites[site];
          if (kinetics.usesSurfacePoreTransport())
          {
            const double fluidContactArea = geometries[ads].contactAreas.fluidSolidPerFluidVolume;
            filmSink += fraction * fluidContactArea * std::max(0.0, kinetics.filmMassTransferCoefficient) *
                        (spanConcentration[grid, comp] - spanSurface[site, grid, comp]);
          }
          else if (kinetics.enabled())
          {
            directFraction += fraction;
          }
        }

        const size_t siteIndex = site * (numberOfGridPoints + 1) * numberOfComponents + index;
        const double reactionChemisorption =
            reactionChemisorptionSource.empty() ? 0.0 : reactionChemisorptionSource[siteIndex];
        sink += loadingPrefactor * directFraction * (spanChemisorptionDot[site, grid, comp] - reactionChemisorption) +
                filmSink;
      }
      bulkSpeciesSink[index] = sink;
    }
  }
}

void computeMassDerivatives(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents,
                            std::span<const double> columnDistances, std::span<const double> fractionOfAdsorbent,
                            std::span<const double> interstitialGasVelocity, std::span<const double> concentration,
                            std::span<double> concentrationDot, std::span<const double> bulkSpeciesSink,
                            MultibedColumn::FluidPhase fluidPhase, std::span<const double> totalVoidFraction)
{
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanConcentrationDot(concentrationDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  auto axialDispersion = [&](size_t grid, size_t comp)
  {
    double coefficient = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      coefficient += fractionOfAdsorbent[grid * numberOfAdsorbents + ads] *
                     physisorptionMixtures[ads].components[comp].axialDispersionCoefficient;
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

  if (fluidPhase == MultibedColumn::FluidPhase::Liquid)
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
                              std::span<const double> totalVoidFraction, std::span<const Geometry> geometries,
                              std::span<const double> particleDensities,
                              std::span<const double> particleDiameters, std::span<const double> fractionOfAdsorbent,
                              double internalDiameter, double outerDiameter, double wallDensity,
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

  // commented out the parts for weno, seems to be unstable
  // std::vector<double> gasTemperatureFlux(numberOfGridPoints + 1);
  // computeWENO(gasTemperature, gasTemperatureFlux);

  // is this extra 1/eps necessary? It's in python not in eqs
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

  // prefactor 4 in python, 2 in eqs
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

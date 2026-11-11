#include "sorption.h"

#include <algorithm>
#include <cmath>
#include <mdspan>
#include <vector>

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

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

void computePhysisorptionEquilibriumSensitivities(
    std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, std::span<const double> fractionOfAdsorbent, const std::vector<bool>& hasAdsorbentOfType,
    size_t maxIsothermTerms, std::span<double> idealGasMolFractions, std::span<double> mixtureSensitivity,
    std::span<const double> totalPressure, std::span<const double> cachedPressure, std::span<const double> moleFraction,
    std::span<const double> gasTemperature, MixturePrediction::DrivingForceInput input,
    std::span<const double> concentration, std::span<const double> pH,
    std::span<double> equilibriumPhysisorptionSensitivity)
{
  const size_t blockSize = numberOfComponents * numberOfComponents;
  const bool gas = input == MixturePrediction::DrivingForceInput::MoleFraction;
  std::fill(equilibriumPhysisorptionSensitivity.begin(), equilibriumPhysisorptionSensitivity.end(), 0.0);

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    // the driving forces of computePhysisorptionEquilibriumLoadings
    double sum = 0.0;
    double totalConcentration = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(
          0.0, gas ? moleFraction[grid * numberOfComponents + comp] : concentration[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
      totalConcentration += std::max(0.0, concentration[grid * numberOfComponents + comp]);
    }
    for (size_t comp = 0; gas && comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] =
          sum > 0.0 ? idealGasMolFractions[comp] / sum : 1.0 / static_cast<double>(numberOfComponents);
    }
    // dp_j/dc_l = (P / C)(delta_jl - y_j) for gas; zero when the node holds no gas.
    const double pressurePerConcentration =
        gas && sum > 0.0 && totalConcentration > 0.0 ? totalPressure[grid] / totalConcentration : 0.0;

    std::span<double> block = equilibriumPhysisorptionSensitivity.subspan(grid * blockSize, blockSize);
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      if (!hasAdsorbentOfType[grid * numberOfAdsorbents + ads]) continue;

      physisorptionMixtures[ads].predictMixtureSensitivities(
          idealGasMolFractions, gas ? totalPressure[grid] : 1.0,
          cachedPressure.subspan((grid * numberOfAdsorbents + ads) * numberOfComponents * maxIsothermTerms,
                                 numberOfComponents * maxIsothermTerms),
          gasTemperature[grid], pH[grid], mixtureSensitivity);

      const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      for (size_t i = 0; i < numberOfComponents; ++i)
      {
        double compositionSlope = 0.0;
        for (size_t j = 0; j < numberOfComponents; ++j)
        {
          compositionSlope += mixtureSensitivity[i + j * numberOfComponents] * idealGasMolFractions[j];
        }
        for (size_t l = 0; l < numberOfComponents; ++l)
        {
          const double slope = mixtureSensitivity[i + l * numberOfComponents];
          block[i + l * numberOfComponents] +=
              fraction * (gas ? pressurePerConcentration * (slope - compositionSlope) : slope);
        }
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

void computePhysisorption(size_t numberOfGridPoints, size_t numberOfComponents,
                          std::span<const double> massTransferCoefficient,
                          std::span<const double> equilibriumPhysisorption, std::span<const double> physisorption,
                          std::span<double> physisorptionDot)
{
  for (size_t index = 0; index < (numberOfGridPoints + 1) * numberOfComponents; ++index)
  {
    physisorptionDot[index] = massTransferCoefficient[index] * (equilibriumPhysisorption[index] - physisorption[index]);
  }
}

std::vector<size_t> chemisorptionSitesPerComponent(const std::vector<MixturePrediction>& mixtures,
                                                   size_t numberOfComponents)
{
  std::vector<size_t> sites(numberOfComponents, 0);
  for (const MixturePrediction& mixture : mixtures)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      sites[comp] = std::max(sites[comp], mixture.components[comp].chemisorption.numberOfSites);
    }
  }
  return sites;
}

void computeChemisorption(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                          size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                          double externalTemperature, std::span<const double> fractionOfAdsorbent,
                          std::span<const double> particleDensities,
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
  const std::vector<size_t> sitesOfComponent = chemisorptionSitesPerComponent(physisorptionMixtures, numberOfComponents);

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    const double temperature = solidTemperature.empty() ? externalTemperature : solidTemperature[grid];
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      for (size_t site = 0; site < sitesOfComponent[comp]; ++site)
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
          double rate =
              kinetics.rate(equilibriumLoading, spanChemisorption[site, grid, comp], drivingConcentration, temperature);
          if (kinetics.usesSurfacePoreTransport())
          {
            // Pore-reservoir concentration removed per unit chemisorbed loading: the reservoir holds the solid
            // volume, so this is the particle density (see computeChemisorptionTransportDerivatives).
            const double gamma = particleDensities[ads];
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

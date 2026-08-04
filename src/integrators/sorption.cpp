#include <algorithm>
#include <mdspan>

#include "sorption.h"
#include "transport.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

void computePhysisorptionEquilibriumLoadings(
    MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxIsothermTerms, std::pair<size_t, size_t>& iastPerformance,
    std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
    std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumPhysisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
    std::span<double> gasTemperature)
{
  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(0.0, moleFraction[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
    }
    if (sum > 0.0)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        idealGasMolFractions[comp] /= sum;
      }
    }
    else
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        idealGasMolFractions[comp] = 1.0 / static_cast<double>(numberOfComponents);
      }
    }

    std::span<double> spanCachedPressure =
        cachedPressure.subspan(grid * numberOfComponents * maxIsothermTerms,
                               numberOfComponents * maxIsothermTerms);
    std::span<double> spanCachedGrandPotential =
        cachedGrandPotential.subspan(grid * maxIsothermTerms, maxIsothermTerms);
    iastPerformance +=
        mixture.predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions,
                               numberOfMolecules, spanCachedPressure, spanCachedGrandPotential,
                               gasTemperature[grid]);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      equilibriumPhysisorption[grid * numberOfComponents + comp] = numberOfMolecules[comp];
    }
  }
}

void computeChemisorptionEquilibriumLoadings(
    MixturePrediction& mixture, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, std::pair<size_t, size_t>& iastPerformance,
    std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
    std::span<double> numberOfMolecules, std::span<const double> totalPressure,
    std::span<double> equilibriumChemisorption, std::span<double> cachedPressure,
    std::span<double> cachedGrandPotential, std::span<const double> moleFraction,
    std::span<double> gasTemperature)
{
  std::fill(equilibriumChemisorption.begin(), equilibriumChemisorption.end(), 0.0);
  const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = std::max(0.0, moleFraction[grid * numberOfComponents + comp]);
      sum += idealGasMolFractions[comp];
    }
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      idealGasMolFractions[comp] = sum > 0.0
                                       ? idealGasMolFractions[comp] / sum
                                       : 1.0 / static_cast<double>(numberOfComponents);
    }

    std::span<double> spanCachedPressure =
        cachedPressure.subspan(grid * numberOfComponents * maxChemisorptionSites,
                               numberOfComponents * maxChemisorptionSites);
    std::span<double> spanCachedGrandPotential =
        cachedGrandPotential.subspan(grid * maxChemisorptionSites, maxChemisorptionSites);
    iastPerformance += mixture.predictMixture(
        idealGasMolFractions, totalPressure[grid], adsorbedMolFractions, numberOfMolecules,
        spanCachedPressure, spanCachedGrandPotential, gasTemperature[grid]);

    for (size_t site = 0; site < mixture.maxIsothermTerms; ++site)
    {
      for (size_t comp = 0; comp < numberOfComponents; ++comp)
      {
        equilibriumChemisorption[site * componentBlockSize + grid * numberOfComponents + comp] =
            mixture.equilibriumSiteLoadings[site * numberOfComponents + comp];
      }
    }
  }
}

void computeSorptionDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, const Geometry& geometry,
    double particleDensity,
    std::span<const double> equilibriumPhysisorption,
    std::span<const double> equilibriumChemisorption, std::span<const double> concentration,
    std::span<const double> physisorption, std::span<double> physisorptionDot,
    std::span<const double> chemisorption, std::span<double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot,
    std::span<const double> solidTemperature, std::span<double> bulkSpeciesSink)
{
  std::fill(physisorptionDot.begin(), physisorptionDot.end(), 0.0);
  std::fill(chemisorptionDot.begin(), chemisorptionDot.end(), 0.0);
  std::fill(surfaceConcentrationDot.begin(), surfaceConcentrationDot.end(), 0.0);
  std::fill(poreConcentrationDot.begin(), poreConcentrationDot.end(), 0.0);

  computePhysisorption(components, numberOfGridPoints, numberOfComponents, equilibriumPhysisorption, physisorption,
                       physisorptionDot);
  computeChemisorption(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                       externalTemperature, geometry, particleDensity, equilibriumChemisorption,
                       concentration, chemisorption, chemisorptionDot, poreConcentration, solidTemperature);
  computeChemisorptionTransportDerivatives(
      components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites, geometry, particleDensity,
      concentration, chemisorptionDot, surfaceConcentration, surfaceConcentrationDot,
      poreConcentration, poreConcentrationDot);
  computeBulkSpeciesSink(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                         geometry, particleDensity, concentration, physisorptionDot,
                         chemisorptionDot, surfaceConcentration, bulkSpeciesSink);
}

void computePhysisorption(const std::vector<Component>& components, size_t numberOfGridPoints,
                          size_t numberOfComponents, std::span<const double> equilibriumPhysisorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot)
{
  mdspan2d_const spanEquilibriumAdsorption(equilibriumPhysisorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanPhysisorption(physisorption.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanPhysisorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  for (size_t grid = 0; grid < numberOfGridPoints + 1; grid++)
  {
    for (size_t comp = 0; comp < numberOfComponents; comp++)
    {
      const double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanPhysisorption[grid, comp];
      spanPhysisorptionDot[grid, comp] = components[comp].massTransferCoefficient * diffAdsorption;
    }
  }
}

void computeChemisorption(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, const Geometry& geometry, double particleDensity,
    std::span<const double> equilibriumChemisorption, std::span<const double> concentration,
    std::span<const double> chemisorption, std::span<double> chemisorptionDot,
    std::span<const double> poreConcentration, std::span<const double> solidTemperature)
{
  mdspan3d_const spanEquilibriumChemisorption(equilibriumChemisorption.data(), maxChemisorptionSites,
                                               numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorption(chemisorption.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                   numberOfComponents);
  mdspan3d_mut spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                    numberOfComponents);
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanPoreConcentration(poreConcentration.data(), maxChemisorptionSites,
                                       numberOfGridPoints + 1, numberOfComponents);

  constexpr double poreInventoryLimitTime = 1.0e-4;

  for (size_t comp = 0; comp < numberOfComponents; comp++)
  {
    const MultiSiteChemisorption& kinetics = components[comp].chemisorption;
    for (size_t site = 0; site < kinetics.numberOfSites; ++site)
    {
      const Chemisorption& siteKinetics = kinetics.sites[site];
      if (!siteKinetics.enabled()) continue;

      for (size_t grid = 0; grid < numberOfGridPoints + 1; grid++)
      {
        const double temperature = solidTemperature.empty() ? externalTemperature : solidTemperature[grid];
        const double drivingConcentration = siteKinetics.usesSurfacePoreTransport()
                                                ? spanPoreConcentration[site, grid, comp]
                                                : spanConcentration[grid, comp];
        double siteEquilibriumLoading = std::max(0.0, spanEquilibriumChemisorption[site, grid, comp]);
        if (siteKinetics.maximumLoading > 0.0)
        {
          siteEquilibriumLoading = std::min(siteEquilibriumLoading, siteKinetics.maximumLoading);
        }
        double rate = siteKinetics.rate(siteEquilibriumLoading, spanChemisorption[site, grid, comp],
                                        drivingConcentration, temperature);
        if (siteKinetics.usesSurfacePoreTransport())
        {
          const double gamma = geometry.loadingPrefactor(particleDensity);
          if (rate > 0.0 && gamma > 0.0)
          {
            rate = std::min(rate, drivingConcentration / (gamma * poreInventoryLimitTime));
          }
        }

        if (spanChemisorption[site, grid, comp] >= siteEquilibriumLoading)
        {
          spanChemisorptionDot[site, grid, comp] = std::min(0.0, rate);
        }
        else if (spanChemisorption[site, grid, comp] <= 0.0)
        {
          spanChemisorptionDot[site, grid, comp] = std::max(0.0, rate);
        }
        else
        {
          spanChemisorptionDot[site, grid, comp] = rate;
        }
      }
    }
  }
}

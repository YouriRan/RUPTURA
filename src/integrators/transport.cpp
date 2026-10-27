#include "transport.h"

#include "sorption.h"

#include <algorithm>
#include <cmath>
#include <mdspan>
#include <vector>

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

void computeChemisorptionTransportDerivatives(
    const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t numberOfAdsorbents, size_t maxChemisorptionSites, std::span<const double> fractionOfAdsorbent,
    std::span<const Geometry> geometries, std::span<const double> particleDensities,
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
  const std::vector<size_t> sitesOfComponent = chemisorptionSitesPerComponent(physisorptionMixtures, numberOfComponents);

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    // The surface and pore reservoirs hold the solid volume (the film flux enters them per solid volume), so the
    // pore reservoir loses the solid mass per solid volume times dq/dt: the solid-volume-weighted particle density,
    // which is rho_p in a pure section.
    double solidLoadingDensity = 0.0;
    double solidVolumeDensity = 0.0;
    for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
    {
      const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
      const double solidVolume = fraction * geometries[ads].solidToFluidVolumeRatio * geometries[ads].voidFraction;
      solidLoadingDensity += solidVolume * particleDensities[ads];
      solidVolumeDensity += solidVolume;
    }
    const double poreLoadingConversion = solidLoadingDensity / std::max(1.0e-30, solidVolumeDensity);

    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      for (size_t site = 0; site < sitesOfComponent[comp]; ++site)
      {
        double film = 0.0;
        double poreDiffusion = 0.0;
        bool poreTransport = false;
        for (size_t ads = 0; ads < numberOfAdsorbents; ++ads)
        {
          const double fraction = fractionOfAdsorbent[grid * numberOfAdsorbents + ads];
          const MultiSiteChemisorption& multisite = physisorptionMixtures[ads].components[comp].chemisorption;
          if (fraction <= 0.0 || site >= multisite.numberOfSites) continue;
          const Chemisorption& kinetics = multisite.sites[site];
          if (!kinetics.usesSurfacePoreTransport()) continue;
          poreTransport = true;

          const double solidContactArea = geometries[ads].contactAreas.solidFluidPerSolidVolume;
          const double poreDiffusionLength =
              std::max(1.0e-30, geometries[ads].dimensions.poreDiffusionLength);
          const double poreLdf = 15.0 * std::max(0.0, kinetics.poreDiffusivity) /
                                 std::max(1.0e-30, poreDiffusionLength * poreDiffusionLength);
          film += fraction * solidContactArea * std::max(0.0, kinetics.filmMassTransferCoefficient) *
                  (spanConcentration[grid, comp] - spanSurface[site, grid, comp]);
          poreDiffusion += fraction * poreLdf * (spanSurface[site, grid, comp] - spanPore[site, grid, comp]);
        }
        // A site that takes up directly from the bulk leaves its surface and pore reservoirs unchanged.
        if (!poreTransport) continue;
        spanSurfaceDot[site, grid, comp] = film - poreDiffusion;
        spanPoreDot[site, grid, comp] =
            poreDiffusion - poreLoadingConversion * spanChemisorptionDot[site, grid, comp];
      }
    }
  }
}

void computeBulkSpeciesSink(const std::vector<MixturePrediction>& physisorptionMixtures, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t numberOfAdsorbents, size_t maxChemisorptionSites,
                            std::span<const Geometry> geometries, std::span<const double> particleDensities,
                            std::span<const double> fractionOfAdsorbent, std::span<const double> totalVoidFraction,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot, std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration, std::span<double> bulkSpeciesSink,
                            std::span<const double> reactionPhysisorptionSource,
                            std::span<const double> reactionChemisorptionSource,
                            std::span<const double> reactionBulkSource)
{
  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                                      numberOfComponents);
  mdspan3d_const spanSurface(surfaceConcentration.data(), maxChemisorptionSites, numberOfGridPoints + 1,
                             numberOfComponents);
  const std::vector<size_t> sitesOfComponent = chemisorptionSitesPerComponent(physisorptionMixtures, numberOfComponents);

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

      for (size_t site = 0; site < sitesOfComponent[comp]; ++site)
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
      if (!reactionBulkSource.empty())
      {
        sink -= reactionBulkSource[index];
      }
      bulkSpeciesSink[index] = sink;
    }
  }
}

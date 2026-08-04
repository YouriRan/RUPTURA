#include "transport.h"

#include <algorithm>
#include <mdspan>

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, const Geometry& geometry, double particleDensity,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot)
{
  if (surfaceConcentration.empty() || poreConcentration.empty()) return;

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites,
                                      numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanSurfaceConcentration(surfaceConcentration.data(), maxChemisorptionSites,
                                          numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanSurfaceConcentrationDot(surfaceConcentrationDot.data(), maxChemisorptionSites,
                                           numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanPoreConcentration(poreConcentration.data(), maxChemisorptionSites,
                                       numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_mut spanPoreConcentrationDot(poreConcentrationDot.data(), maxChemisorptionSites,
                                        numberOfGridPoints + 1, numberOfComponents);

  const double solidContactArea = geometry.contactAreas.solidFluidPerSolidVolume;
  const double poreDiffusionLength = std::max(geometry.dimensions.poreDiffusionLength, 1.0e-30);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const MultiSiteChemisorption& multisite = components[comp].chemisorption;
    for (size_t site = 0; site < multisite.numberOfSites; ++site)
    {
      const Chemisorption& kinetics = multisite.sites[site];
      if (!kinetics.usesSurfacePoreTransport()) continue;

      const double kFilm = std::max(0.0, kinetics.filmMassTransferCoefficient);
      const double poreLdf =
          15.0 * std::max(0.0, kinetics.poreDiffusivity) /
          std::max(poreDiffusionLength * poreDiffusionLength, 1.0e-30);
      const double gamma = geometry.loadingPrefactor(particleDensity);

      for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
      {
        const double film = solidContactArea * kFilm *
                            (spanConcentration[grid, comp] - spanSurfaceConcentration[site, grid, comp]);
        const double poreDiffusion =
            poreLdf * (spanSurfaceConcentration[site, grid, comp] - spanPoreConcentration[site, grid, comp]);
        spanSurfaceConcentrationDot[site, grid, comp] = film - poreDiffusion;
        spanPoreConcentrationDot[site, grid, comp] =
            poreDiffusion - gamma * spanChemisorptionDot[site, grid, comp];
      }
    }
  }
}

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double voidFraction, double particleDensity, double particleDiameter,
    std::span<const double> concentration, std::span<const double> chemisorptionDot,
    std::span<const double> surfaceConcentration, std::span<double> surfaceConcentrationDot,
    std::span<const double> poreConcentration, std::span<double> poreConcentrationDot)
{
  const Geometry geometry = makeGeometry(PackedBedTubeSpec{.voidFraction = voidFraction,
                                                            .particleDiameter = particleDiameter});
  computeChemisorptionTransportDerivatives(components, numberOfGridPoints, numberOfComponents,
                                           maxChemisorptionSites, geometry, particleDensity,
                                           concentration, chemisorptionDot, surfaceConcentration,
                                           surfaceConcentrationDot, poreConcentration, poreConcentrationDot);
}

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            const Geometry& geometry, double particleDensity,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink,
                            std::span<const double> reactionPhysisorptionSource,
                            std::span<const double> reactionChemisorptionSource)
{
  std::fill(bulkSpeciesSink.begin(), bulkSpeciesSink.end(), 0.0);

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanPhysisorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites,
                                      numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  const double loadingPrefactor = geometry.loadingPrefactor(particleDensity);
  const double filmPrefactor = geometry.contactAreas.fluidSolidPerFluidVolume;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      double sink = loadingPrefactor * spanPhysisorptionDot[grid, comp];
      if (!reactionPhysisorptionSource.empty())
      {
        sink -= loadingPrefactor * reactionPhysisorptionSource[grid * numberOfComponents + comp];
      }
      const MultiSiteChemisorption& multisite = components[comp].chemisorption;
      for (size_t site = 0; site < multisite.numberOfSites; ++site)
      {
        const Chemisorption& kinetics = multisite.sites[site];
        if (kinetics.usesSurfacePoreTransport())
        {
          mdspan3d_const spanSurfaceConcentration(surfaceConcentration.data(), maxChemisorptionSites,
                                                  numberOfGridPoints + 1, numberOfComponents);
          sink += filmPrefactor * std::max(0.0, kinetics.filmMassTransferCoefficient) *
                  (spanConcentration[grid, comp] - spanSurfaceConcentration[site, grid, comp]);
        }
        else
        {
          double chemisorptionSource = spanChemisorptionDot[site, grid, comp];
          if (!reactionChemisorptionSource.empty())
          {
            const size_t componentBlockSize = (numberOfGridPoints + 1) * numberOfComponents;
            chemisorptionSource -=
                reactionChemisorptionSource[site * componentBlockSize + grid * numberOfComponents + comp];
          }
          sink += loadingPrefactor * chemisorptionSource;
        }
      }
      spanBulkSpeciesSink[grid, comp] = sink;
    }
  }
}

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            double voidFraction, double particleDensity, double particleDiameter,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink)
{
  const Geometry geometry = makeGeometry(PackedBedTubeSpec{.voidFraction = voidFraction,
                                                            .particleDiameter = particleDiameter});
  computeBulkSpeciesSink(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                         geometry, particleDensity, concentration, physisorptionDot,
                         chemisorptionDot, surfaceConcentration, bulkSpeciesSink, std::span<const double>{},
                         std::span<const double>{});
}

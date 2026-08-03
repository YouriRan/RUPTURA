#include "transport.h"

#include <algorithm>
#include <mdspan>

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

void computeChemisorptionTransportDerivatives(Column& column)
{
  computeChemisorptionTransportDerivatives(
      column.components, column.numberOfGridPoints, column.numberOfComponents, column.maxChemisorptionSites,
      column.geometry.shapeParameters(), column.particleDensity, column.concentration,
      column.chemisorptionDot, column.surfaceConcentration, column.surfaceConcentrationDot,
      column.poreConcentration, column.poreConcentrationDot);
}

void computeChemisorptionTransportDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, const ShapeParameters& geometry, double particleDensity,
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

  const double solidContactArea = geometry.solidFluidContactAreaPerSolidVolume;
  const double poreDiffusionLength = std::max(geometry.poreDiffusionLength, 1.0e-30);

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
  const Geometry geometry{HollowTube{voidFraction, particleDiameter}};
  computeChemisorptionTransportDerivatives(components, numberOfGridPoints, numberOfComponents,
                                           maxChemisorptionSites, geometry.shapeParameters(), particleDensity,
                                           concentration, chemisorptionDot, surfaceConcentration,
                                           surfaceConcentrationDot, poreConcentration, poreConcentrationDot);
}

void computeBulkSpeciesSink(Column& column)
{
  computeBulkSpeciesSink(column.components, column.numberOfGridPoints, column.numberOfComponents,
                         column.maxChemisorptionSites, column.geometry.shapeParameters(), column.particleDensity,
                         column.concentration, column.physisorptionDot, column.chemisorptionDot,
                         column.surfaceConcentration, column.bulkSpeciesSink);
}

void computeBulkSpeciesSink(const std::vector<Component>& components, size_t numberOfGridPoints,
                            size_t numberOfComponents, size_t maxChemisorptionSites,
                            const ShapeParameters& geometry, double particleDensity,
                            std::span<const double> concentration,
                            std::span<const double> physisorptionDot,
                            std::span<const double> chemisorptionDot,
                            std::span<const double> surfaceConcentration,
                            std::span<double> bulkSpeciesSink)
{
  std::fill(bulkSpeciesSink.begin(), bulkSpeciesSink.end(), 0.0);

  mdspan2d_const spanConcentration(concentration.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_const spanPhysisorptionDot(physisorptionDot.data(), numberOfGridPoints + 1, numberOfComponents);
  mdspan3d_const spanChemisorptionDot(chemisorptionDot.data(), maxChemisorptionSites,
                                      numberOfGridPoints + 1, numberOfComponents);
  mdspan2d_mut spanBulkSpeciesSink(bulkSpeciesSink.data(), numberOfGridPoints + 1, numberOfComponents);

  const double loadingPrefactor = geometry.loadingPrefactor(particleDensity);
  const double filmPrefactor = geometry.fluidSolidContactAreaPerFluidVolume;

  for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
  {
    for (size_t comp = 0; comp < numberOfComponents; ++comp)
    {
      double sink = loadingPrefactor * spanPhysisorptionDot[grid, comp];
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
          sink += loadingPrefactor * spanChemisorptionDot[site, grid, comp];
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
  const Geometry geometry{HollowTube{voidFraction, particleDiameter}};
  computeBulkSpeciesSink(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                         geometry.shapeParameters(), particleDensity, concentration, physisorptionDot,
                         chemisorptionDot, surfaceConcentration, bulkSpeciesSink);
}

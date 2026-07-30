#include <algorithm>
#include <mdspan>

#include "sorption.h"
#include "transport.h"

using mdspan2d_const = std::mdspan<const double, std::dextents<size_t, 2>>;
using mdspan2d_mut = std::mdspan<double, std::dextents<size_t, 2>>;
using mdspan3d_const = std::mdspan<const double, std::dextents<size_t, 3>>;
using mdspan3d_mut = std::mdspan<double, std::dextents<size_t, 3>>;

void computeSorptionDerivatives(Column& column)
{
  computeSorptionDerivatives(
      column.components, column.numberOfGridPoints, column.numberOfComponents, column.maxChemisorptionSites,
      column.externalTemperature, column.voidFraction, column.particleDensity, column.particleDiameter,
      column.equilibriumAdsorption,
      column.concentration, column.physisorption, column.physisorptionDot, column.chemisorption,
      column.chemisorptionDot, column.surfaceConcentration, column.surfaceConcentrationDot,
      column.poreConcentration, column.poreConcentrationDot, column.solidTemperature, column.bulkSpeciesSink);
}

void computeSorptionDerivatives(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, double voidFraction, double particleDensity,
    double particleDiameter,
    std::span<const double> equilibriumAdsorption, std::span<const double> concentration,
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

  computePhysisorption(components, numberOfGridPoints, numberOfComponents, equilibriumAdsorption, physisorption,
                       physisorptionDot);
  computeChemisorption(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                       externalTemperature, voidFraction, particleDensity, equilibriumAdsorption, concentration,
                       chemisorption, chemisorptionDot, poreConcentration, solidTemperature);
  computeChemisorptionTransportDerivatives(
      components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites, voidFraction, particleDensity,
      particleDiameter, concentration, chemisorptionDot, surfaceConcentration, surfaceConcentrationDot,
      poreConcentration, poreConcentrationDot);
  computeBulkSpeciesSink(components, numberOfGridPoints, numberOfComponents, maxChemisorptionSites,
                         voidFraction, particleDensity, particleDiameter, concentration, physisorptionDot,
                         chemisorptionDot, surfaceConcentration, bulkSpeciesSink);
}

void computePhysisorption(Column& column)
{
  computePhysisorption(column.components, column.numberOfGridPoints, column.numberOfComponents,
                       column.equilibriumAdsorption, column.physisorption, column.physisorptionDot);
}

void computePhysisorption(const std::vector<Component>& components, size_t numberOfGridPoints,
                          size_t numberOfComponents, std::span<const double> equilibriumAdsorption,
                          std::span<const double> physisorption, std::span<double> physisorptionDot)
{
  mdspan2d_const spanEquilibriumAdsorption(equilibriumAdsorption.data(), numberOfGridPoints + 1, numberOfComponents);
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

void computeChemisorption(Column& column)
{
  computeChemisorption(column.components, column.numberOfGridPoints, column.numberOfComponents,
                       column.maxChemisorptionSites, column.externalTemperature, column.voidFraction,
                       column.particleDensity, column.equilibriumAdsorption, column.concentration,
                       column.chemisorption, column.chemisorptionDot, column.poreConcentration,
                       column.solidTemperature);
}

void computeChemisorption(
    const std::vector<Component>& components, size_t numberOfGridPoints, size_t numberOfComponents,
    size_t maxChemisorptionSites, double externalTemperature, double voidFraction, double particleDensity,
    std::span<const double> equilibriumAdsorption, std::span<const double> concentration,
    std::span<const double> chemisorption, std::span<double> chemisorptionDot,
    std::span<const double> poreConcentration, std::span<const double> solidTemperature)
{
  mdspan2d_const spanEquilibriumAdsorption(equilibriumAdsorption.data(), numberOfGridPoints + 1, numberOfComponents);
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
        const double siteEquilibriumLoading =
            kinetics.equilibriumLoading(site, spanEquilibriumAdsorption[grid, comp]);
        double rate = siteKinetics.rate(siteEquilibriumLoading, spanChemisorption[site, grid, comp],
                                        drivingConcentration, temperature);
        if (siteKinetics.usesSurfacePoreTransport())
        {
          const double epsp = std::clamp(voidFraction, 1.0e-12, 1.0 - 1.0e-12);
          const double gamma = particleDensity * (1.0 - epsp) / epsp;
          if (rate > 0.0 && gamma > 0.0)
          {
            rate = std::min(rate, drivingConcentration / (gamma * poreInventoryLimitTime));
          }
        }

        if (siteKinetics.maximumLoading > 0.0 &&
            spanChemisorption[site, grid, comp] >= siteKinetics.maximumLoading)
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

#include <span>
#include <vector>

#include "column.h"
#include "mixture_prediction.h"

void computeEquilibriumLoadings(Column& column);
void computeEquilibriumLoadings(MixturePrediction& mixture, size_t Ngrid, size_t Ncomp, size_t maxIsothermTerms,
                                std::pair<size_t, size_t>& iastPerformance, double pressureGradient,
                                double columnLength, std::span<double> totalPressure, std::span<double> partialPressure,
                                std::span<double> idealGasMolFractions, std::span<double> adsorbedMolFractions,
                                std::span<double> numberOfMolecules, std::span<double> equilibriumAdsorption,
                                std::span<double> cachedPressure, std::span<double> cachedGrandPotential);

void computeVelocity(Column& column);
void computeVelocityFixedGradient(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp,
                                  double pressureGradient, double columnEntranceVelocity, double resolution,
                                  std::span<const double> prefactorMassTransfer,
                                  std::span<double> interstitialGasVelocity, std::span<const double> totalPressure,
                                  std::span<const double> adsorption, std::span<const double> equilibriumAdsorption,
                                  std::span<const double> partialPressure);
void computeVelocityErgun(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp,
                          double externalTemperature, double voidFraction, double columnEntranceVelocity,
                          double dynamicViscosity, double particleDiameter, double resolution,
                          std::span<double> interstitialGasVelocity, std::span<const double> totalPressure,
                          std::span<const double> partialPressure);

void computeFirstDerivatives(Column& column);
void computeFirstDerivatives(const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double resolution,
                             std::span<const double> prefactorMassTransfer,
                             std::span<const double> interstitialGasVelocity, std::span<const double> partialPressure,
                             std::span<double> partialPressureDot, std::span<const double> adsorption,
                             std::span<double> adsorptionDot, std::span<const double> equilibriumAdsorption,
                             std::span<double> totalPressureDot);

void computeFirstDerivativesEnergyBalance(
    const std::vector<Component>& components, size_t Ngrid, size_t Ncomp, double externalTemperature,
    double voidFraction, double particleDensity, double particleDiameter, double internalDiameter, double outerDiameter,
    double wallDensity, double gasThermalConductivity, double wallThermalConductivity, double heatTransferGasSolid,
    double heatTransferGasWall, double heatTransferWallExternal, double heatCapacityGas, double heatCapacitySolid,
    double heatCapacityWall, double resolution, std::span<const double> prefactorMassTransfer,
    std::span<const double> interstitialGasVelocity, std::span<const double> totalPressure,
    std::span<double> gasTemperature, std::span<double> gasTemperatureDot, std::span<double> solidTemperature,
    std::span<double> solidTemperatureDot, std::span<double> wallTemperature, std::span<double> wallTemperatureDot,
    std::span<double> partialPressure, std::span<double> partialPressureDot, std::span<double> adsorption,
    std::span<double> adsorptionDot, std::span<double> equilibriumAdsorption, std::span<double> moleFraction,
    std::span<double> gasDensity, std::span<double> coeffGasGas, std::span<double> coeffGasSolid,
    std::span<double> coeffGasWall, std::span<double> coeffDiffusion, std::span<double> facePressures,
    std::span<double> energyFlow);

void computeFirstDerivativesWENO(Column& column);

std::vector<double> computeWENO(std::vector<double>& vec, bool clamp);
std::vector<double> computeTVD(std::vector<double>& vec, bool clamp);
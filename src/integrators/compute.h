#include <span>
#include <vector>

#include "column.h"
#include "mixture_prediction.h"

void computeEquilibriumLoadings(Column& column);
void computeEquilibriumLoadings(size_t Ncomp, size_t Ngrid, std::span<double> totalPressure,
                                std::span<double> partialPressure, std::span<double> idealGasMolFractions,
                                std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules,
                                std::span<double> cachedPressure, std::span<double> cachedGrandPotential,
                                std::span<double> equilibriumAdsorption, MixturePrediction mixture,
                                std::pair<size_t, size_t>& iastPerformance, double pressureGradient,
                                double columnLength, size_t maxIsothermTerms);

void computeVelocity(Column& column);
void computeVelocity(size_t Ncomp, size_t Ngrid, double resolution, std::span<double> interstitialGasVelocity,
                     double columnEntranceVelocity, double pressureGradient, std::span<const double> totalPressure,
                     std::span<const double> prefactorMassTransfer, std::span<const double> equilibriumAdsorption,
                     std::span<const double> adsorption, const std::vector<Component>& components,
                     std::span<const double> partialPressure);

void computeFirstDerivatives(Column& column);
void computeFirstDerivatives(size_t Ncomp, size_t Ngrid, double resolution, std::span<const double> partialPressure,
                             std::span<const double> equilibriumAdsorption, std::span<const double> adsorption,
                             std::span<double> adsorptionDot, std::span<double> partialPressureDot,
                             std::span<const double> interstitialGasVelocity,
                             std::span<const double> prefactorMassTransfer, const std::vector<Component>& components);

void computeFirstDerivativesWENO(Column& column);

std::vector<double> computeWENO(std::vector<double>& vec, bool clamp);
std::vector<double> computeTVD(std::vector<double>& vec, bool clamp);
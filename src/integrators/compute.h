#include "breakthrough_state.h"
#include "mixture_prediction.h"

void computeEquilibriumLoadings(BreakthroughState& state);
void computeEquilibriumLoadings(size_t Ncomp, size_t Ngrid, std::vector<double>& totalPressure,
                                std::vector<double>& partialPressure, std::vector<double>& idealGasMolFractions,
                                std::vector<double>& adsorbedMolFractions, std::vector<double>& numberOfMolecules,
                                std::vector<double>& cachedPressure, std::vector<double>& cachedGrandPotential,
                                std::vector<double>& equilibriumAdsorption, MixturePrediction mixture,
                                std::pair<size_t, size_t>& iastPerformance, double pressureGradient,
                                double columnLength, size_t maxIsothermTerms);

void computeVelocity(BreakthroughState& state);
void computeVelocity(size_t Ncomp, size_t Ngrid, double resolution, std::vector<double>& interstitialGasVelocity,
                     double columnEntranceVelocity, double pressureGradient, const std::vector<double>& totalPressure,
                     const std::vector<double>& prefactorMassTransfer, const std::vector<double>& equilibriumAdsorption,
                     const std::vector<double>& adsorption, const std::vector<Component>& components,
                     const std::vector<double>& partialPressure);

void computeFirstDerivatives(BreakthroughState& state);
void computeFirstDerivatives(size_t Ncomp, size_t Ngrid, double resolution, const std::vector<double>& partialPressure,
                             const std::vector<double>& equilibriumAdsorption, const std::vector<double>& adsorption,
                             std::vector<double>& adsorptionDot, std::vector<double>& pressureDot,
                             const std::vector<double>& interstitialGasVelocity,
                             const std::vector<double>& prefactorMassTransfer,
                             const std::vector<Component>& components);

void computeFirstDerivativesWENO(BreakthroughState& state);

std::vector<double> computeWENO(std::vector<double>& vec, bool clamp);
std::vector<double> computeTVD(std::vector<double>& vec, bool clamp);
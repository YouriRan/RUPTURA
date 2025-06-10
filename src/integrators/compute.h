#include "breakthrough_state.h"
#include "mixture_prediction.h"

void computeEquilibriumLoadings(BreakthroughState& state, MixturePrediction& mixture);
void computeVelocity(BreakthroughState& state, std::vector<Component>& components);
void computeFirstDerivatives(BreakthroughState& state, std::vector<Component>& components);
void computeFirstDerivativesWENO(BreakthroughState& state, std::vector<Component>& components);
std::vector<double> computeWENO(std::vector<double>& vec, bool clamp);
std::vector<double> computeTVD(std::vector<double>& vec, bool clamp);
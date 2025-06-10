#include "compute.h"

#include <mdspan>

#include "utils.h"

void computeEquilibriumLoadings(BreakthroughState& state, MixturePrediction& mixture)
{
  size_t Ncomp = state.Ncomp;
  size_t Ngrid = state.Ngrid;

  // calculate new equilibrium loadings Qeqnew corresponding to the new timestep
  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    // estimation of total pressure Pt at each grid point from partial pressures
    state.totalPressure[grid] = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      state.totalPressure[grid] += std::max(0.0, state.partialPressure[grid * Ncomp + comp]);
    }

    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      state.idealGasMolFractions[comp] = std::max(state.partialPressure[grid * Ncomp + comp], 0.0);
      sum += state.idealGasMolFractions[comp];
    }
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      state.idealGasMolFractions[comp] /= sum;
    }

    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    state.iastPerformance +=
        mixture.predictMixture(state.idealGasMolFractions, state.totalPressure[grid], state.adsorbedMolFractions,
                               state.numberOfMolecules, &state.cachedPressure[grid * Ncomp * state.maxIsothermTerms],
                               &state.cachedGrandPotential[grid * state.maxIsothermTerms]);

    for (size_t comp = 0; comp < Ncomp; ++comp)

    {
      state.equilibriumAdsorption[grid * Ncomp + comp] = state.numberOfMolecules[comp];
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (state.totalPressure[0] + state.pressureGradient * state.columnLength < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeVelocity(BreakthroughState& state, std::vector<Component>& components)
{
  double idx2 = 1.0 / (state.resolution * state.resolution);
  size_t Ncomp = state.Ncomp;
  size_t Ngrid = state.Ngrid;

  // first grid point
  state.interstitialGasVelocity[0] = state.columnEntranceVelocity;

  // middle gridpoints
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    // sum = derivative at the actual gridpoint i
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      sum = sum -
            state.prefactorMassTransfer[comp] *
                (state.equilibriumAdsorption[grid * Ncomp + comp] - state.adsorption[grid * Ncomp + comp]) +
            components[comp].D *
                (state.partialPressure[(grid - 1) * Ncomp + comp] - 2.0 * state.partialPressure[grid * Ncomp + comp] +
                 state.partialPressure[(grid + 1) * Ncomp + comp]) *
                idx2;
    }

    // explicit version
    state.interstitialGasVelocity[grid] = state.interstitialGasVelocity[grid - 1] +
                                          state.resolution *
                                              (sum - state.interstitialGasVelocity[grid - 1] * state.pressureGradient) /
                                              state.totalPressure[grid];
  }

  // last grid point
  double sum = 0.0;
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    sum = sum -
          state.prefactorMassTransfer[comp] *
              (state.equilibriumAdsorption[Ngrid * Ncomp + comp] - state.adsorption[Ngrid * Ncomp + comp]) +
          components[comp].D *
              (state.partialPressure[(Ngrid - 1) * Ncomp + comp] - state.partialPressure[Ngrid * Ncomp + comp]) * idx2;
  }

  // explicit version
  state.interstitialGasVelocity[Ngrid] = state.interstitialGasVelocity[Ngrid - 1] +
                                         state.resolution *
                                             (sum - state.interstitialGasVelocity[Ngrid - 1] * state.pressureGradient) /
                                             state.totalPressure[Ngrid];
}

void computeFirstDerivatives(BreakthroughState& state, std::vector<Component>& components)
{
  size_t Ncomp = state.Ncomp;
  size_t Ngrid = state.Ngrid;

  double idx = 1.0 / state.resolution;
  double idx2 = idx * idx;

  std::mdspan<const double, std::dextents<size_t, 2>> spanPartialPressure(state.partialPressure.data(), Ngrid + 1,
                                                                          Ncomp);
  std::mdspan<const double, std::dextents<size_t, 2>> spanEquilibriumAdsorption(state.equilibriumAdsorption.data(),
                                                                                Ngrid + 1, Ncomp);
  std::mdspan<const double, std::dextents<size_t, 2>> spanAdsorption(state.adsorption.data(), Ngrid + 1, Ncomp);

  // first gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    state.adsorptionDot[comp] = components[comp].Kl * diffAdsorption;
    state.pressureDot[comp] = 0.0;
  }

  // middle gridpoints
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      state.adsorptionDot[grid * Ncomp + comp] = components[comp].Kl * diffAdsorption;
      state.pressureDot[grid * Ncomp + comp] =
          (state.interstitialGasVelocity[grid - 1] * spanPartialPressure[grid - 1, comp] -
           state.interstitialGasVelocity[grid] * spanPartialPressure[grid, comp]) *
              idx +
          components[comp].D *
              (spanPartialPressure[grid + 1, comp] - 2.0 * spanPartialPressure[grid, comp] +
               spanPartialPressure[grid - 1, comp]) *
              idx2 -
          state.prefactorMassTransfer[comp] * diffAdsorption;
    }
  }

  // last gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[Ngrid, comp] - spanAdsorption[Ngrid, comp];

    state.adsorptionDot[Ngrid * Ncomp + comp] = components[comp].Kl * diffAdsorption;
    state.pressureDot[Ngrid * Ncomp + comp] =
        (state.interstitialGasVelocity[Ngrid - 1] * spanPartialPressure[Ngrid - 1, comp] -
         state.interstitialGasVelocity[Ngrid] * spanPartialPressure[Ngrid, comp]) *
            idx +
        components[comp].D * (spanPartialPressure[Ngrid - 1, comp] - spanPartialPressure[Ngrid, comp]) * idx2 -
        state.prefactorMassTransfer[comp] * diffAdsorption;
  }
}

std::vector<double> computeWENO(std::vector<double>& vec, bool clamp)
{
  if (vec.empty()) return vec;

  double tol = 1e-10;
  double alpha_0, alpha_1, first_term, second_term;

  size_t size = vec.size();
  std::vector<double> out(vec);

  if (clamp && vec[size - 1] >= 1.0) out[size - 1] = 1.0;

  // For right wall of 1st Node, alpha_1 and second term values are calculated using half-cell approximation
  alpha_0 = (2.0 / 3.0) / std::pow((vec[2] - vec[1]) + tol, 4.0);
  alpha_1 = (1.0 / 3.0) / std::pow(2.0 * (vec[1] - vec[0]) + tol, 4.0);

  first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (vec[2] + vec[1]);
  second_term = (alpha_1 / (alpha_0 + alpha_1)) * (2.0 * vec[1] - vec[0]);
  out[1] = first_term + second_term;

  // For right walls of 2nd to Ngrid-1 node
  for (size_t i = 2; i < size - 1; ++i)
  {
    alpha_0 = (2.0 / 3.0) / std::pow(vec[i + 1] - vec[i] + tol, 4.0);
    alpha_1 = (1.0 / 3.0) / std::pow(vec[i] - vec[i - 1] + tol, 4.0);

    first_term = 0.5 * (alpha_0 / (alpha_0 + alpha_1)) * (vec[i + 1] + vec[i]);
    second_term = (alpha_1 / (alpha_0 + alpha_1)) * ((3.0 / 2.0) * vec[i] - (1.0 / 2.0) * vec[i - 1]);
    out[i] = first_term + second_term;
  }
  return out;
}

std::vector<double> computeTVD(std::vector<double>& vec, bool clamp)
{
  if (vec.empty()) return vec;

  double tol = 1e-10;
  double r_value, flux_limiter;

  size_t size = vec.size();
  std::vector<double> out(vec);

  if (clamp && vec[size - 1] >= 1.0) out[size - 1] = 1.0;

  // For right wall of 1st Node, r_value is calculated using half-cell approximation
  r_value = (2.0 * (vec[1] - vec[0]) + tol) / (vec[2] - vec[1] + tol);
  flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
  out[1] += 0.5 * flux_limiter * (vec[2] - vec[1]);

  // For right walls of 2nd to Ngrid-1 node
  for (size_t i = 2; i < size - 1; ++i)
  {
    r_value = ((vec[i] - vec[i - 1]) + tol) / ((vec[i + 1] - vec[i]) + tol);
    flux_limiter = (r_value + std::abs(r_value)) / (1.0 + std::abs(r_value));
    out[i] = vec[i] + 0.5 * flux_limiter * (vec[i + 1] - vec[i]);
  }
  return out;
}

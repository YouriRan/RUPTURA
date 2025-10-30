#include "compute.h"

#include <mdspan>

#include "utils.h"

void computeEquilibriumLoadings(Column& state)
{
  computeEquilibriumLoadings(state.Ncomp, state.Ngrid, state.totalPressure, state.partialPressure,
                             state.idealGasMolFractions, state.adsorbedMolFractions, state.numberOfMolecules,
                             state.cachedPressure, state.cachedGrandPotential, state.equilibriumAdsorption,
                             state.mixture, state.iastPerformance, state.pressureGradient, state.columnLength,
                             state.maxIsothermTerms);
}

void computeEquilibriumLoadings(size_t Ncomp, size_t Ngrid, std::span<double> totalPressure,
                                std::span<double> partialPressure, std::span<double> idealGasMolFractions,
                                std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules,
                                std::span<double> cachedPressure, std::span<double> cachedGrandPotential,
                                std::span<double> equilibriumAdsorption, MixturePrediction mixture,
                                std::pair<size_t, size_t>& iastPerformance, double pressureGradient,
                                double columnLength, size_t maxIsothermTerms)
{
  // calculate new equilibrium loadings Qeqnew corresponding to the new timestep
  for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  {
    // estimation of total pressure Pt at each grid point from partial pressures
    totalPressure[grid] = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      totalPressure[grid] += std::max(0.0, partialPressure[grid * Ncomp + comp]);
    }

    // compute gas-phase mol-fractions
    // force the gas-phase mol-fractions to be positive and normalized
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      idealGasMolFractions[comp] = std::max(partialPressure[grid * Ncomp + comp], 0.0);
      sum += idealGasMolFractions[comp];
    }
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      idealGasMolFractions[comp] /= sum;
    }

    // use Yi and Pt[i] to compute the loadings in the adsorption mixture via mixture prediction
    iastPerformance += mixture.predictMixture(idealGasMolFractions, totalPressure[grid], adsorbedMolFractions,
                                              numberOfMolecules, &cachedPressure[grid * Ncomp * maxIsothermTerms],
                                              &cachedGrandPotential[grid * maxIsothermTerms]);

    for (size_t comp = 0; comp < Ncomp; ++comp)

    {
      equilibriumAdsorption[grid * Ncomp + comp] = numberOfMolecules[comp];
    }
  }

  // check the total pressure at the outlet, it should not be negative
  if (totalPressure[0] + pressureGradient * columnLength < 0.0)
  {
    throw std::runtime_error("Error: pressure gradient is too large (negative outlet pressure)\n");
  }
}

void computeVelocity(Column& state)
{
  computeVelocity(state.Ncomp, state.Ngrid, state.resolution, state.interstitialGasVelocity,
                  state.columnEntranceVelocity, state.pressureGradient, state.totalPressure,
                  state.prefactorMassTransfer, state.equilibriumAdsorption, state.adsorption, state.components,
                  state.partialPressure);
}

void computeVelocity(size_t Ncomp, size_t Ngrid, double resolution, std::span<double> interstitialGasVelocity,
                     double columnEntranceVelocity, double pressureGradient, std::span<const double> totalPressure,
                     std::span<const double> prefactorMassTransfer, std::span<const double> equilibriumAdsorption,
                     std::span<const double> adsorption, const std::vector<Component>& components,
                     std::span<const double> partialPressure)
{
  double idx2 = 1.0 / (resolution * resolution);

  // first grid point
  interstitialGasVelocity[0] = columnEntranceVelocity;

  // middle grid points
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      // mass transfer term
      sum -=
          prefactorMassTransfer[comp] * (equilibriumAdsorption[grid * Ncomp + comp] - adsorption[grid * Ncomp + comp]);

      // diffusion term
      sum += components[comp].D *
             (partialPressure[(grid - 1) * Ncomp + comp] - 2.0 * partialPressure[grid * Ncomp + comp] +
              partialPressure[(grid + 1) * Ncomp + comp]) *
             idx2;
    }

    // explicit update
    interstitialGasVelocity[grid] =
        interstitialGasVelocity[grid - 1] +
        resolution * (sum - interstitialGasVelocity[grid - 1] * pressureGradient) / totalPressure[grid];
  }

  // last grid point
  {
    double sum = 0.0;
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      sum -= prefactorMassTransfer[comp] *
             (equilibriumAdsorption[Ngrid * Ncomp + comp] - adsorption[Ngrid * Ncomp + comp]);

      sum += components[comp].D *
             (partialPressure[(Ngrid - 1) * Ncomp + comp] - partialPressure[Ngrid * Ncomp + comp]) * idx2;
    }

    interstitialGasVelocity[Ngrid] =
        interstitialGasVelocity[Ngrid - 1] +
        resolution * (sum - interstitialGasVelocity[Ngrid - 1] * pressureGradient) / totalPressure[Ngrid];
  }
}

void computeFirstDerivatives(Column& state)
{
  computeFirstDerivatives(state.Ncomp, state.Ngrid, state.resolution, state.partialPressure,
                          state.equilibriumAdsorption, state.adsorption, state.adsorptionDot, state.partialPressureDot,
                          state.interstitialGasVelocity, state.prefactorMassTransfer, state.components);
}

void computeFirstDerivatives(size_t Ncomp, size_t Ngrid, double resolution, std::span<const double> partialPressure,
                             std::span<const double> equilibriumAdsorption, std::span<const double> adsorption,
                             std::span<double> adsorptionDot, std::span<double> partialPressureDot,
                             std::span<const double> interstitialGasVelocity,
                             std::span<const double> prefactorMassTransfer, const std::vector<Component>& components)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  std::mdspan<const double, std::dextents<size_t, 2>> spanPartialPressure(partialPressure.data(), Ngrid + 1, Ncomp);
  std::mdspan<const double, std::dextents<size_t, 2>> spanEquilibriumAdsorption(equilibriumAdsorption.data(), Ngrid + 1,
                                                                                Ncomp);
  std::mdspan<const double, std::dextents<size_t, 2>> spanAdsorption(adsorption.data(), Ngrid + 1, Ncomp);

  // first gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[0, comp] - spanAdsorption[0, comp];
    adsorptionDot[comp] = components[comp].Kl * diffAdsorption;
    partialPressureDot[comp] = 0.0;
  }

  // middle gridpoints
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      double diffAdsorption = spanEquilibriumAdsorption[grid, comp] - spanAdsorption[grid, comp];

      adsorptionDot[grid * Ncomp + comp] = components[comp].Kl * diffAdsorption;

      partialPressureDot[grid * Ncomp + comp] =
          (interstitialGasVelocity[grid - 1] * spanPartialPressure[grid - 1, comp] -
           interstitialGasVelocity[grid] * spanPartialPressure[grid, comp]) *
              idx +
          components[comp].D *
              (spanPartialPressure[grid + 1, comp] - 2.0 * spanPartialPressure[grid, comp] +
               spanPartialPressure[grid - 1, comp]) *
              idx2 -
          prefactorMassTransfer[comp] * diffAdsorption;
    }
  }

  // last gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double diffAdsorption = spanEquilibriumAdsorption[Ngrid, comp] - spanAdsorption[Ngrid, comp];

    adsorptionDot[Ngrid * Ncomp + comp] = components[comp].Kl * diffAdsorption;

    partialPressureDot[Ngrid * Ncomp + comp] =
        (interstitialGasVelocity[Ngrid - 1] * spanPartialPressure[Ngrid - 1, comp] -
         interstitialGasVelocity[Ngrid] * spanPartialPressure[Ngrid, comp]) *
            idx +
        components[comp].D * (spanPartialPressure[Ngrid - 1, comp] - spanPartialPressure[Ngrid, comp]) * idx2 -
        prefactorMassTransfer[comp] * diffAdsorption;
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

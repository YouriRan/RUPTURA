#include "ruptura2/particle_transport.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace ruptura2 {
namespace {
constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();

void finite_nonnegative(double value, const char* label) {
  if (!std::isfinite(value) || value < 0.0)
    throw std::invalid_argument(std::string(label) + " must be finite and nonnegative");
}

bool dense_solve(std::vector<double> matrix, std::vector<double> rhs,
                 std::vector<double>& result) {
  const std::size_t n = rhs.size();
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot = k;
    for (std::size_t row = k + 1; row < n; ++row)
      if (std::abs(matrix[row * n + k]) > std::abs(matrix[pivot * n + k])) pivot = row;
    if (!std::isfinite(matrix[pivot * n + k]) || matrix[pivot * n + k] == 0.0)
      return false;
    if (pivot != k) {
      for (std::size_t col = k; col < n; ++col)
        std::swap(matrix[k * n + col], matrix[pivot * n + col]);
      std::swap(rhs[k], rhs[pivot]);
    }
    for (std::size_t row = k + 1; row < n; ++row) {
      const double ratio = matrix[row * n + k] / matrix[k * n + k];
      for (std::size_t col = k + 1; col < n; ++col)
        matrix[row * n + col] -= ratio * matrix[k * n + col];
      rhs[row] -= ratio * rhs[k];
    }
  }
  result.assign(n, 0.0);
  for (std::size_t row = n; row-- > 0;) {
    double value = rhs[row];
    for (std::size_t col = row + 1; col < n; ++col)
      value -= matrix[row * n + col] * result[col];
    result[row] = value / matrix[row * n + row];
    if (!std::isfinite(result[row])) return false;
  }
  return true;
}
}  // namespace

RadialParticleMesh::RadialParticleMesh(std::size_t cells, int radial_power)
    : power(radial_power), faces(cells + 1), centers(cells), weights(cells) {
  if (cells == 0 || (power != 1 && power != 2))
    throw std::invalid_argument("Particle mesh requires at least one cell and power 1 or 2");
  for (std::size_t face = 0; face <= cells; ++face)
    faces[face] = static_cast<double>(face) / static_cast<double>(cells);
  for (std::size_t cell = 0; cell < cells; ++cell) {
    centers[cell] = 0.5 * (faces[cell] + faces[cell + 1]);
    weights[cell] = std::pow(faces[cell + 1], power + 1) -
                    std::pow(faces[cell], power + 1);
  }
}

double particle_mean(const RadialParticleMesh& mesh, std::span<const double> loading) {
  if (loading.size() != mesh.weights.size())
    throw std::invalid_argument("Particle loading size does not match radial mesh");
  return std::inner_product(mesh.weights.begin(), mesh.weights.end(), loading.begin(), 0.0);
}

void particle_diffusion_rhs(const RadialParticleMesh& mesh,
                            std::span<const double> loading,
                            double q_surface, double diffusivity_over_radius_squared,
                            std::span<double> derivative) {
  const std::size_t n = mesh.weights.size();
  if (loading.size() != n || derivative.size() != n)
    throw std::invalid_argument("Particle RHS sizes do not match radial mesh");
  finite_nonnegative(diffusivity_over_radius_squared, "D_s/R^2");
  if (!std::isfinite(q_surface)) throw std::invalid_argument("Particle surface loading must be finite");
  const double multiplier = static_cast<double>(mesh.power + 1) *
                            diffusivity_over_radius_squared;
  double inner_flux = 0.0;
  for (std::size_t cell = 0; cell < n; ++cell) {
    const double outer_flux = cell + 1 < n
        ? std::pow(mesh.faces[cell + 1], mesh.power) *
          (loading[cell + 1] - loading[cell]) /
          (mesh.centers[cell + 1] - mesh.centers[cell])
        : (q_surface - loading[cell]) / (1.0 - mesh.centers[cell]);
    derivative[cell] = multiplier * (outer_flux - inner_flux) / mesh.weights[cell];
    inner_flux = outer_flux;
  }
}

double arrhenius_from_reference(double reference_coefficient, double activation_energy,
                               double temperature, double reference_temperature,
                               double gas_constant) {
  finite_nonnegative(reference_coefficient, "Reference transport coefficient");
  if (!std::isfinite(activation_energy) || !std::isfinite(temperature) ||
      !std::isfinite(reference_temperature) || !std::isfinite(gas_constant) ||
      temperature <= 0.0 || reference_temperature <= 0.0 || gas_constant <= 0.0)
    throw std::invalid_argument("Arrhenius temperature/reference temperature/gas constant must be positive and finite");
  if (reference_coefficient == 0.0) return 0.0;
  const double exponent = -activation_energy / gas_constant *
                          (1.0 / temperature - 1.0 / reference_temperature);
  // Work in log space so an overflowing exp(exponent) does not discard a
  // perfectly representable result when the reference coefficient is small.
  const double result = std::exp(std::log(reference_coefficient) + exponent);
  if (!std::isfinite(result)) throw std::overflow_error("Temperature-dependent particle transport coefficient overflow");
  return result;
}

ParticleSurfaceResult solve_particle_surface(
    std::span<const double> cb, std::span<const double> q,
    std::span<const double> a, std::span<const double> kf,
    const ParticleEquilibrium& equilibrium, ParticleSurfaceOptions options) {
  const std::size_t n = cb.size();
  if (q.size() != n || a.size() != n || kf.size() != n || !equilibrium)
    throw std::invalid_argument("Particle film closure requires matching vectors and an equilibrium callback");
  if (!(options.relative_tolerance > 0.0) || !std::isfinite(options.relative_tolerance) ||
      !(options.absolute_flux_tolerance > 0.0) || !std::isfinite(options.absolute_flux_tolerance) ||
      options.newton_iterations < 0 || options.coordinate_sweeps < 0)
    throw std::invalid_argument("Invalid particle surface solver tolerances/iteration limits");
  std::vector<double> upper(n), cs(cb.begin(), cb.end());
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < n; ++i) {
    finite_nonnegative(cb[i], "Bulk concentration");
    finite_nonnegative(q[i], "Particle loading");
    finite_nonnegative(a[i], "Loading transfer coefficient");
    finite_nonnegative(kf[i], "Film coefficient");
    if (a[i] > 0.0 && kf[i] == 0.0)
      throw std::invalid_argument("An active finite-film species requires kf>0; kf=0 requires an impermeable zero-flux particle boundary");
    const long double bound = static_cast<long double>(cb[i]) + (a[i] == 0.0 ? 0.0L :
        static_cast<long double>(a[i]) * static_cast<long double>(q[i]) / kf[i]);
    if (bound > std::numeric_limits<double>::max())
      throw std::overflow_error("Particle surface concentration bound overflow");
    upper[i] = static_cast<double>(bound);
    if (a[i] > 0.0 && upper[i] > 0.0) active.push_back(i);
  }
  auto eq_at = [&](const std::vector<double>& c) {
    auto qe = equilibrium(c);
    if (qe.size() != n) throw std::runtime_error("Particle equilibrium callback returned wrong component count");
    for (double value : qe)
      if (!std::isfinite(value) || value < 0.0)
        throw std::runtime_error("Particle equilibrium callback must return finite nonnegative loadings");
    return qe;
  };
  auto residual = [&](std::size_t i, double c, double qe) {
    return a[i] * (qe - q[i]) - kf[i] * (cb[i] - c);
  };
  auto tolerance = [&](std::size_t i, double c, double qe) {
    const double particle_flux = a[i] * (qe - q[i]);
    const double film_flux = kf[i] * (cb[i] - c);
    const double scale = std::max({std::abs(particle_flux), std::abs(film_flux),
                                  kf[i] * cb[i]});
    // Include only cancellation-level rounding allowance, not an arbitrary
    // equilibrium-loading scale that could hide an unbalanced film flux.
    return options.absolute_flux_tolerance + options.relative_tolerance * scale +
           32.0 * machine_epsilon * a[i] * std::max(qe, q[i]);
  };
  auto norm = [&](const std::vector<double>& c, const std::vector<double>& qe) {
    double value = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      const double ri = residual(i, c[i], qe[i]);
      if (!std::isfinite(ri)) throw std::overflow_error("Particle film residual overflow");
      value = std::max(value, std::abs(ri) / tolerance(i, c[i], qe[i]));
    }
    return value;
  };
  auto qe = eq_at(cs);
  double error = norm(cs, qe);
  int iterations = 0;

  auto newton_step = [&]() {
    const std::size_t count = active.size();
    if (count == 0) return false;
    std::vector<double> jacobian(count * count), rhs(count), row_scale(count);
    for (std::size_t row = 0; row < count; ++row) {
      const std::size_t i = active[row];
      row_scale[row] = tolerance(i, cs[i], qe[i]);
      rhs[row] = -residual(i, cs[i], qe[i]) / row_scale[row];
    }
    for (std::size_t column = 0; column < count; ++column) {
      const std::size_t j = active[column];
      auto shifted = cs;
      const double step = 1.0e-6 * std::max({cs[j], cb[j], 1.0e-12});
      shifted[j] += step;
      const auto shifted_q = eq_at(shifted);
      for (std::size_t row = 0; row < count; ++row) {
        const std::size_t i = active[row];
        jacobian[row * count + column] =
            (a[i] * ((shifted_q[i] - qe[i]) / step) + (i == j ? kf[i] : 0.0)) /
            row_scale[row];
      }
    }
    std::vector<double> direction;
    if (!dense_solve(std::move(jacobian), std::move(rhs), direction)) return false;
    double step_length = 1.0;
    for (int search = 0; search < 32; ++search, step_length *= 0.5) {
      auto trial = cs;
      for (std::size_t j = 0; j < count; ++j) {
        const std::size_t i = active[j];
        trial[i] = std::clamp(cs[i] + step_length * direction[j], 0.0, upper[i]);
      }
      auto trial_q = eq_at(trial);
      const double trial_error = norm(trial, trial_q);
      if (trial_error <= 1.0 || trial_error < error * (1.0 - 1.0e-4 * step_length)) {
        cs.swap(trial);
        qe.swap(trial_q);
        error = trial_error;
        return true;
      }
    }
    return false;
  };

  for (int iteration = 0; error > 1.0 && iteration < options.newton_iterations; ++iteration) {
    ++iterations;
    if (!newton_step()) break;
  }

  // Coordinate brackets guarantee nonnegative concentrations for continuous
  // own-species monotone isotherms. Competitive coupling is re-evaluated at
  // every trial, and all residuals are checked after each complete sweep.
  for (int sweep = 0; error > 1.0 && sweep < options.coordinate_sweeps; ++sweep) {
    ++iterations;
    for (std::size_t i : active) {
      double lower = 0.0, higher = upper[i];
      auto trial = cs;
      trial[i] = lower;
      auto lower_q = eq_at(trial);
      double fl = residual(i, lower, lower_q[i]);
      trial[i] = higher;
      auto higher_q = eq_at(trial);
      double fh = residual(i, higher, higher_q[i]);
      if (fl > tolerance(i, lower, lower_q[i]) || fh < -tolerance(i, higher, higher_q[i]))
        throw std::runtime_error("No nonnegative particle film root for component " + std::to_string(i) +
                                 "; check particle loading and equilibrium domain");
      if (std::abs(fl) <= tolerance(i, lower, lower_q[i])) {
        cs[i] = lower;
      } else if (std::abs(fh) <= tolerance(i, higher, higher_q[i])) {
        cs[i] = higher;
      } else {
        for (int bisect = 0; bisect < 100; ++bisect) {
          const double midpoint = lower + 0.5 * (higher - lower);
          trial[i] = midpoint;
          const auto midpoint_q = eq_at(trial);
          const double fm = residual(i, midpoint, midpoint_q[i]);
          cs[i] = midpoint;
          if (std::abs(fm) <= 0.1 * tolerance(i, midpoint, midpoint_q[i]) ||
              midpoint == lower || midpoint == higher) break;
          if (fm > 0.0) higher = midpoint;
          else lower = midpoint;
        }
      }
    }
    qe = eq_at(cs);
    error = norm(cs, qe);
    if (error > 1.0) (void)newton_step();
  }
  if (error > 1.0)
    throw std::runtime_error("Coupled particle film closure did not converge after " +
        std::to_string(iterations) + " iterations; scaled flux residual=" + std::to_string(error));
  return {std::move(cs), std::move(qe), error, iterations};
}

}  // namespace ruptura2

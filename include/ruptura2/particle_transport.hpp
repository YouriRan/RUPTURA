#pragma once

#include <cstddef>
#include <functional>
#include <span>
#include <vector>

namespace ruptura2 {

// Uniform cell-centred mesh in x=r/R. power=2 is a sphere; power=1 is a
// cylinder with radial diffusion only (impermeable ends / sufficiently long).
// weights are exact fractions of particle volume and sum to one.
struct RadialParticleMesh {
  explicit RadialParticleMesh(std::size_t cells, int radial_power);
  int power;
  std::vector<double> faces;
  std::vector<double> centers;
  std::vector<double> weights;
};

double particle_mean(const RadialParticleMesh& mesh,
                     std::span<const double> loading);

// q_surface is the equilibrium loading at the physical exterior x=1.
// D_s/R^2 has units 1/s. Face fluxes telescope exactly in the volume-weighted
// mean. The origin has zero flux and the exterior uses a half-cell distance.
void particle_diffusion_rhs(const RadialParticleMesh& mesh,
                            std::span<const double> loading,
                            double q_surface, double diffusivity_over_radius_squared,
                            std::span<double> derivative);

// The supplied coefficient is its value at Tref, not an Arrhenius prefactor.
double arrhenius_from_reference(double reference_coefficient,
                               double activation_energy,
                               double temperature,
                               double reference_temperature,
                               double gas_constant = 8.31446261815324);

using ParticleEquilibrium =
    std::function<std::vector<double>(std::span<const double>)>;

struct ParticleSurfaceOptions {
  double relative_tolerance{1.0e-10};
  double absolute_flux_tolerance{1.0e-14};
  int newton_iterations{40};
  int coordinate_sweeps{160};
};

struct ParticleSurfaceResult {
  std::vector<double> concentration;
  std::vector<double> equilibrium_loading;
  double max_scaled_residual{0.0};
  int iterations{0};
};

// Solve all surface concentrations together using the complete competitive
// equilibrium callback. Each species obeys
//   a_i [qeq_i(cs)-q_i] = kf_i [cb_i-cs_i].
// a_i = rho_p D_s/(dr/2) for HSDM, or rho_p R k_LDF/(power+1) for LDF.
// The entire concentration vector, including inactive species, is passed to
// equilibrium. Desorption cs>cb is permitted. The returned equilibrium and
// surface concentrations satisfy the same, unclipped flux balance.
// a=0 sets cs=cb. kf=0 with a>0 is rejected: an impermeable particle requires
// a zero-flux boundary and cannot be represented by this equilibrium closure.
// Throws on invalid data or failure to find a nonnegative converged solution.
ParticleSurfaceResult solve_particle_surface(
    std::span<const double> bulk_concentration,
    std::span<const double> particle_loading,
    std::span<const double> loading_transfer_coefficient,
    std::span<const double> film_coefficient,
    const ParticleEquilibrium& equilibrium,
    ParticleSurfaceOptions options = {});

}  // namespace ruptura2

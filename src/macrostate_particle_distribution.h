#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <vector>

/**
 * \brief Inclusive particle-number range for one MPD component.
 *
 * Bounds are ordered like the non-carrier entries in Components. The last
 * bound is the fastest-varying dimension in the C-ordered distribution file.
 */
struct MPDComponentBounds
{
  std::string component;  ///< Optional component name used to verify dimension ordering.
  std::size_t nMin{0};    ///< Inclusive minimum particle number.
  std::size_t nMax{0};    ///< Inclusive maximum particle number.
  std::size_t deltaN{1};  ///< Particle-number increment.
};

/**
 * \brief Input needed to reweight a macrostate particle distribution.
 */
struct MPDSettings
{
  std::string fileName;                             ///< One/two-column C-ordered distribution file.
  double referenceTemperature{0.0};                 ///< Reference temperature in K.
  double referenceFugacity{0.0};                    ///< Common reference fugacity in Pa.
  double referenceFrameworkMass{0.0};               ///< Framework mass represented by the MPD in kg.
  std::vector<MPDComponentBounds> componentBounds;  ///< Bounds for each non-carrier component.
};

/**
 * \brief Loads and reweights a rank-N macrostate particle distribution.
 *
 * The input contains one row for every macrostate in C order. Column one is
 * Pi_ref(n); optional column two is the conditional mean Hamiltonian in J.
 */
class MacrostateParticleDistribution
{
 public:
  explicit MacrostateParticleDistribution(MPDSettings settings);

  /**
   * \brief Reweight the reference distribution and return mean particle counts.
   *
   * gasMoleFractions follows the MPD dimension order, fugacity is in Pa, and
   * temperature is in K. Normalization is a log-sum-exp in two passes: the
   * log weights and their maximum, then exp() of the states within 50 of it.
   * A non-empty covariance (rank x rank) receives <n_a n_b> - <n_a><n_b>,
   * which is d<n_a>/d ln f_b.
   */
  std::vector<double> meanParticleNumbers(std::span<const double> gasMoleFractions, double fugacity,
                                          double temperature, std::span<double> covariance = {}) const;

  [[nodiscard]] std::size_t rank() const noexcept { return settings_.componentBounds.size(); }
  [[nodiscard]] std::size_t numberOfMacrostates() const noexcept { return numberOfMacrostates_; }
  [[nodiscard]] bool hasMeanEnergies() const noexcept { return !states_.meanEnergies.empty(); }
  [[nodiscard]] double referenceFrameworkMass() const noexcept { return settings_.referenceFrameworkMass; }
  [[nodiscard]] const MPDSettings& settings() const noexcept { return settings_; }

 private:
  /// Populated macrostates (Pi_ref > 0) in file order, one contiguous array per quantity.
  struct States
  {
    std::vector<long double> logProbabilities;         ///< log Pi_ref(n).
    std::vector<std::vector<double>> particleNumbers;  ///< n_j of every state, one array per component.
    std::vector<double> meanEnergies;                  ///< <H>_n,ref in J; empty without an energy column.
  };

  MPDSettings settings_;
  std::vector<std::size_t> extents_;
  std::vector<std::size_t> strides_;
  std::size_t numberOfMacrostates_{0};  ///< Dense macrostate count, including Pi_ref = 0 entries.

  States states_;
  std::vector<States> emptyComponentStates_;  ///< Per component j, the states with n_j = 0.

  [[nodiscard]] std::size_t particleNumber(std::size_t linearIndex, std::size_t component) const noexcept;
};

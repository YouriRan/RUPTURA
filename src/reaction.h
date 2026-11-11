#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include "temperature_law.h"

/**
 * \brief Homogeneous or adsorbed-phase reaction coupled to column component balances.
 */
struct Reaction
{
  /**
   * \brief State variable on which the reaction takes place.
   */
  enum class Phase
  {
    Physisorbed = 0,
    Chemisorbed = 1,
    PoreConcentration = 2,
    BulkConcentration = 3  ///< Bulk-fluid concentrations; the rate is per fluid volume (catalyst amount folded in).
  };

  /**
   * \brief Kinetic expression used to turn local activities into a net reaction rate.
   */
  enum class Style
  {
    GeneralPowerLaw = 0,
    LangmuirHinshelwood = 1,             ///< k+ prod theta_r^alpha (irreversible), on the competitive coverages theta.
    LangmuirHinshelwoodHougenWatson = 2  ///< k+ prod theta_r^alpha - k- prod theta_p^beta.
  };

  Phase phase{Phase::PoreConcentration};
  Style style{Style::GeneralPowerLaw};
  size_t site{0};
  std::vector<size_t> reactants{};
  std::vector<size_t> products{};
  std::vector<double> reactantStoichiometry{};
  std::vector<double> productStoichiometry{};
  std::vector<double> forwardOrders{};
  std::vector<double> backwardOrders{};
  /// Forward rate coefficient k+, usually an Arrhenius law; its unit depends on the reaction orders.
  TemperatureLaw forwardRateCoefficient{0.0};
  /// Equilibrium constant K, usually a VantHoff law; the backward coefficient is k- = k+/K.
  TemperatureLaw equilibriumConstant{1.0};
  /// Heat released per mole of reaction in J/mol (positive for exothermic), Q_r = -Delta H_r.
  double heatOfReaction{0.0};
  /// LH/LHHW adsorption constants K_j(T), one per component (zero for a component that does not adsorb), in the
  /// inverse unit of the activity (m^3/mol for concentrations). Coverages are theta_j = K_j c_j / (1 + sum_i K_i c_i).
  std::vector<TemperatureLaw> adsorptionConstants{};
  double rateLimitTime{1.0e-4};
  /// Multibed adsorbent indices the reaction is restricted to; empty means everywhere. The rate at a node is
  /// weighted by the summed local fraction of these adsorbents.
  std::vector<size_t> adsorbents{};

  [[nodiscard]] bool enabled() const noexcept { return forwardRateCoefficient.value != 0.0; }
  [[nodiscard]] double equilibriumConstantAt(double temperature) const;
  [[nodiscard]] double forwardRateConstant(double temperature) const;
  /// Net rate at the activities and temperature; coverage is a work array of at least activity.size() values for the
  /// LH and LHHW coverages.
  [[nodiscard]] double rate(std::span<const double> activity, double temperature, std::span<double> coverage) const;
  [[nodiscard]] std::string repr() const;
};

[[nodiscard]] bool reactionsRequirePoreConcentration(const std::vector<Reaction>& reactions) noexcept;

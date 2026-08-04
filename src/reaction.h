#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <vector>

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
    PoreConcentration = 2
  };

  /**
   * \brief Kinetic expression used to turn local activities into a net reaction rate.
   */
  enum class Style
  {
    GeneralPowerLaw = 0,
    LangmuirHinshelwood = 1,
    LangmuirHinshelwoodHougenWatson = 2
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
  double forwardRateCoefficient{0.0};
  double forwardActivationEnergy{0.0};
  double equilibriumConstant{1.0};
  double gibbsFreeEnergy{0.0};
  double rateLimitTime{1.0e-4};

  [[nodiscard]] bool enabled() const noexcept { return forwardRateCoefficient != 0.0; }
  [[nodiscard]] double stoichiometricCoefficient(size_t component) const noexcept;
  [[nodiscard]] double equilibriumConstantAt(double temperature) const;
  [[nodiscard]] double forwardRateConstant(double temperature) const;
  [[nodiscard]] double rate(std::span<const double> activity, double temperature) const;
  [[nodiscard]] std::string repr() const;
};

[[nodiscard]] bool reactionsRequirePoreConcentration(const std::vector<Reaction>& reactions) noexcept;

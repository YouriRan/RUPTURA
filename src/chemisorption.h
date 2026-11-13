#pragma once

#include <cstddef>
#include <array>
#include <optional>
#include <span>
#include <string>

#include "isotherm.h"

struct Chemisorption
{
  enum class Type
  {
    None = 0,
    FirstOrder = 1,
    PseudoNth = 2,
    Avrami = 3,
    General = 4,
    Elovich = 5,
    DualNth = 6
  };

  struct Branch
  {
    double fraction{0.5};
    double adsorptionRate{0.0}, desorptionRate{0.0};
    double adsorptionActivationEnergy{0.0}, desorptionActivationEnergy{0.0};
    double adsorptionOrder{1.0}, desorptionOrder{1.0};
    std::optional<double> heat;
  };
  struct Rates
  {
    std::array<double, 2> loading{0.0, 0.0};
    double adsorptionAge{0.0}, desorptionAge{0.0};
  };

  Type type{Type::None};
  double rateCoefficient{0.0};
  double order{1.0};
  double maximumLoading{0.0};
  double heatOfChemisorption{0.0};
  double adsorptionRateCoefficient{0.0};
  double adsorptionActivationEnergy{0.0};
  double desorptionRateCoefficient{0.0};
  double desorptionActivationEnergy{0.0};
  size_t poreConcentrationOrder{1};
  size_t capacityOrder{1};
  size_t desorptionOrder{1};
  double elovichAlpha{0.0};
  double elovichBeta{0.0};
  // Optional coverage-dependent desorption for FirstOrder (Stadler et al., Front. Chem. Eng. 4 (2022) 1000064,
  // Eqs. 7-10): when q > q*, the coefficient is k_des = k1 exp(beta (q/q_ref)/(R T)) instead of rateCoefficient.
  double desorptionMassTransferCoefficient{0.0};  ///< k1 in 1/s; zero keeps the symmetric FirstOrder law.
  double desorptionCoverageEnergy{0.0};           ///< beta in J/mol.
  double coverageReferenceLoading{0.0};           ///< q_ref in mol/kg; Stadler et al. use q*(20 bar).
  double filmMassTransferCoefficient{0.0};
  double poreDiffusivity{0.0};
  bool usePoreSurfaceTransport{false};
  std::optional<Isotherm> isotherm{std::nullopt};
  // All loading-based laws share this kernel for gas and liquid. Equilibrium
  // is supplied by the column in the appropriate pressure/concentration units.
  bool directionalRates{false};
  double referenceTemperature{298.15};
  double desorptionOrderValue{1.0};
  bool fractionalRates{false};
  double referenceLoading{1.0};
  double drivingRegularizer{1.0e-12}, ageRegularizer{1.0e-6}, historyDrivingRegularizer{0.0};
  bool continuousHistory{false};
  double adsorptionMemoryErasure{-1.0}, desorptionMemoryErasure{-1.0};
  std::array<Branch, 2> branches{};
  // Optional chemical-equilibrium method; empty inherits the column method.
  std::string mixtureMethod;
  bool physicalEquilibrium{false};
  double equilibriumFraction{1.0}, equilibriumHeat{0.0};
  std::optional<double> equilibriumReferenceTemperature;
  std::optional<double> adsorptionReferenceTemperature;

  [[nodiscard]] std::size_t branchCount() const noexcept { return type == Type::DualNth ? 2 : 1; }
  [[nodiscard]] double branchFraction(std::size_t b) const noexcept
  { return type == Type::DualNth ? branches[b].fraction : (b == 0 ? 1.0 : 0.0); }
  [[nodiscard]] double branchHeat(std::size_t b) const noexcept
  { return type == Type::DualNth ? branches[b].heat.value_or(heatOfChemisorption) : heatOfChemisorption; }
  [[nodiscard]] Rates evaluate(double equilibriumLoading, std::span<const double> loadings,
                               double concentration, double temperature, double elapsedTime = 0.0,
                               double adsorptionAge = 0.0, double desorptionAge = 0.0) const;

  [[nodiscard]] double rate(double equilibriumLoading, double loading, double concentration = 0.0,
                            double temperature = 298.15, double elapsedTime = 0.0) const;
  [[nodiscard]] bool enabled() const noexcept { return type != Type::None; }
  [[nodiscard]] bool usesEquilibriumLoading() const noexcept
  {
    return type == Type::FirstOrder || type == Type::PseudoNth || type == Type::Avrami || type == Type::DualNth;
  }
  [[nodiscard]] bool usesSurfacePoreTransport() const noexcept
  {
    return (type == Type::General || type == Type::Elovich) && usePoreSurfaceTransport;
  }
  [[nodiscard]] std::string repr() const;
};

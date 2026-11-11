#pragma once

#include <cstddef>
#include <optional>
#include <string>

#include "isotherm.h"
#include "temperature_law.h"

struct Chemisorption
{
  enum class Type
  {
    None = 0,
    FirstOrder = 1,
    PseudoNth = 2,
    Avrami = 3,
    General = 4,
    Elovich = 5
  };

  Type type{Type::None};
  // The rate coefficients may depend on temperature (Arrhenius laws); rate() evaluates them at its temperature.
  TemperatureLaw rateCoefficient{0.0};  ///< FirstOrder, PseudoNth, and Avrami rate coefficient k.
  double order{1.0};
  double maximumLoading{0.0};
  double heatOfChemisorption{0.0};  ///< Heat released per mole chemisorbed on this site, in J/mol.
  TemperatureLaw adsorptionRateCoefficient{0.0};  ///< General adsorption coefficient k+.
  TemperatureLaw desorptionRateCoefficient{0.0};  ///< General desorption coefficient k-.
  size_t poreConcentrationOrder{1};
  size_t capacityOrder{1};
  size_t desorptionOrder{1};
  TemperatureLaw elovichAlpha{0.0};  ///< Elovich alpha in m^3/(kg s).
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

  [[nodiscard]] double rate(double equilibriumLoading, double loading, double concentration = 0.0,
                            double temperature = 298.15) const;
  [[nodiscard]] bool enabled() const noexcept { return type != Type::None; }
  [[nodiscard]] bool usesEquilibriumLoading() const noexcept
  {
    return type == Type::FirstOrder || type == Type::PseudoNth || type == Type::Avrami;
  }
  [[nodiscard]] bool usesSurfacePoreTransport() const noexcept
  {
    return (type == Type::General || type == Type::Elovich) && usePoreSurfaceTransport;
  }
  [[nodiscard]] std::string repr() const;
};

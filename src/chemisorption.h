#pragma once

#include <cstddef>
#include <string>

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
  double rateCoefficient{0.0};
  size_t order{1};
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
  double filmMassTransferCoefficient{0.0};
  double poreDiffusivity{0.0};
  bool usePoreSurfaceTransport{false};

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

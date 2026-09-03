
#include "chemisorption.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>

#include "utils.h"

double Chemisorption::rate(double equilibriumLoading, double loading, double concentration, double temperature,
                           double elapsedTime) const
{
  if (type == Type::None) return 0.0;

  const double qeq = maximumLoading > 0.0 ? std::min(equilibriumLoading, maximumLoading) : equilibriumLoading;
  const double driving = qeq - loading;

  switch (type)
  {
    case Type::None:
      return 0.0;
    case Type::FirstOrder:
      return rateCoefficient * driving;
    case Type::PseudoNth:
      return rateCoefficient * std::pow(driving, order);
    case Type::Avrami:
      if (elapsedTime <= 0.0)
      {
        return order == 1.0 ? rateCoefficient * driving : 0.0;
      }
      return order * std::pow(rateCoefficient, order) * std::pow(elapsedTime, order - 1.0) * driving;
    case Type::General:
    {
      const double c = std::max(0.0, concentration);
      const double q = std::max(0.0, loading);
      const double freeCapacity = std::max(0.0, qeq - q);
      const double ka =
          arrhenius(adsorptionRateCoefficient, adsorptionActivationEnergy, temperature);
      const double kd =
          arrhenius(desorptionRateCoefficient, desorptionActivationEnergy, temperature);
      const double adsorption =
          ka * std::pow(c, poreConcentrationOrder) * std::pow(freeCapacity, capacityOrder);
      const double desorption = kd * std::pow(q, desorptionOrder);
      return adsorption - desorption;
    }
    case Type::Elovich:
      return elovichAlpha * std::max(0.0, concentration) *
             std::exp(-elovichBeta * std::max(0.0, loading));
  }

  return 0.0;
}

std::string Chemisorption::repr() const
{
  std::string typeName = "None";
  switch (type)
  {
    case Type::None:
      typeName = "None";
      break;
    case Type::FirstOrder:
      typeName = "FirstOrder";
      break;
    case Type::PseudoNth:
      typeName = "PseudoNth";
      break;
    case Type::Avrami:
      typeName = "Avrami";
      break;
    case Type::General:
      typeName = "General";
      break;
    case Type::Elovich:
      typeName = "Elovich";
      break;
  }

  std::string text = std::format(
      "    chemisorption:\n"
      "        type:                  {}\n"
      "        rate coefficient:      {:8.4e} [1/s]\n"
      "        order:                 {} [-]\n"
      "        maximum loading:       {:8.4e} [mol/kg]\n"
      "        heat of chemisorption: {:8.4e} [J/mol]\n",
      typeName, rateCoefficient, order, maximumLoading, heatOfChemisorption);
  if (type == Type::General)
  {
    text += std::format(
        "        adsorption k0:         {:8.4e}\n"
        "        adsorption Ea:         {:8.4e} [J/mol]\n"
        "        desorption k0:         {:8.4e}\n"
        "        desorption Ea:         {:8.4e} [J/mol]\n"
        "        pore conc. order:      {} [-]\n"
        "        capacity order:        {} [-]\n"
        "        desorption order:      {} [-]\n"
        "        film coefficient:      {:8.4e} [m/s]\n"
        "        pore diffusivity:      {:8.4e} [m^2/s]\n",
        adsorptionRateCoefficient, adsorptionActivationEnergy, desorptionRateCoefficient, desorptionActivationEnergy,
        poreConcentrationOrder, capacityOrder, desorptionOrder, filmMassTransferCoefficient, poreDiffusivity);
    text += std::format("        pore/surface transport: {}\n", usePoreSurfaceTransport ? "enabled" : "disabled");
  }
  else if (type == Type::Elovich)
  {
    text += std::format(
        "        alpha:                 {:8.4e}\n"
        "        beta:                  {:8.4e} [kg/mol]\n"
        "        film coefficient:      {:8.4e} [m/s]\n"
        "        pore diffusivity:      {:8.4e} [m^2/s]\n"
        "        pore/surface transport: {}\n",
        elovichAlpha, elovichBeta, filmMassTransferCoefficient, poreDiffusivity,
        usePoreSurfaceTransport ? "enabled" : "disabled");
  }
  if (isotherm.has_value())
  {
    text += "        equilibrium model:\n";
    text += isotherm->repr();
  }
  return text;
}

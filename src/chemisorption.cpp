
#include "chemisorption.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>

#include "utils.h"

double Chemisorption::rate(double equilibriumLoading, double loading, double concentration, double temperature) const
{
  if (type == Type::None) return 0.0;

  const double qeq = maximumLoading > 0.0 ? std::min(equilibriumLoading, maximumLoading) : equilibriumLoading;
  const double driving = qeq - loading;

  switch (type)
  {
    case Type::None:
      return 0.0;
    case Type::FirstOrder:
      if (driving < 0.0 && desorptionMassTransferCoefficient > 0.0)
      {
        const double coverage = std::max(0.0, loading) / coverageReferenceLoading;
        const double exponent = desorptionCoverageEnergy * coverage / (R * std::max(temperature, 1.0e-10));
        return desorptionMassTransferCoefficient * std::exp(std::min(exponent, 700.0)) * driving;
      }
      return rateCoefficient(temperature) * driving;
    case Type::PseudoNth:
      // Sign-preserving k |Δq|^n sign(Δq): desorbs for q > q*, and stays finite at Δq = 0 for n < 1.
      return rateCoefficient(temperature) * std::copysign(std::pow(std::abs(driving), order), driving);
    case Type::Avrami:
    {
      // Autonomous (Avrami-Erofeev) form of q/q* = 1 - exp(-(k t)^n), without a clock (B-PHY-12): the loading sets
      // the equivalent Avrami time s = k t, s^n = -ln(1 - q/q*) = ln(q*/(q* - q)), and dq/dt = n k s^(n-1) (q* - q).
      // At q = 0 the law has a zero (n > 1) or an infinite (n < 1) rate; the floor on s starts the uptake and moves
      // the curve in time by (1/n - 1) s_min / k. Release (q > q*) is linear driving force.
      const double k = rateCoefficient(temperature);
      if (driving <= 0.0) return k * driving;
      constexpr double minimumAvramiTime = 1.0e-3;
      const double logarithm = qeq > 0.0 ? std::log(qeq / std::min(driving, qeq)) : 0.0;
      return order * k * std::pow(std::max(logarithm, std::pow(minimumAvramiTime, order)), (order - 1.0) / order) *
             driving;
    }
    case Type::General:
    {
      const double c = std::max(0.0, concentration);
      const double q = std::max(0.0, loading);
      const double freeCapacity = std::max(0.0, qeq - q);
      const double ka = adsorptionRateCoefficient(temperature);
      const double kd = desorptionRateCoefficient(temperature);
      const double adsorption =
          ka * std::pow(c, poreConcentrationOrder) * std::pow(freeCapacity, capacityOrder);
      const double desorption = kd * std::pow(q, desorptionOrder);
      return adsorption - desorption;
    }
    case Type::Elovich:
      return elovichAlpha(temperature) * std::max(0.0, concentration) *
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
      "        rate coefficient:      {:8.4e} [1/s]  {}\n"
      "        order:                 {} [-]\n"
      "        maximum loading:       {:8.4e} [mol/kg]\n"
      "        heat of chemisorption: {:8.4e} [J/mol]\n",
      typeName, rateCoefficient.value, rateCoefficient.repr(), order, maximumLoading, heatOfChemisorption);
  if (type == Type::FirstOrder && desorptionMassTransferCoefficient > 0.0)
  {
    text += std::format(
        "        desorption k1:         {:8.4e} [1/s]\n"
        "        desorption coverage energy: {:8.4e} [J/mol]\n"
        "        coverage reference loading: {:8.4e} [mol/kg]\n",
        desorptionMassTransferCoefficient, desorptionCoverageEnergy, coverageReferenceLoading);
  }
  if (type == Type::General)
  {
    text += std::format(
        "        adsorption k:          {:8.4e}  {}\n"
        "        desorption k:          {:8.4e}  {}\n"
        "        pore conc. order:      {} [-]\n"
        "        capacity order:        {} [-]\n"
        "        desorption order:      {} [-]\n"
        "        film coefficient:      {:8.4e} [m/s]\n"
        "        pore diffusivity:      {:8.4e} [m^2/s]\n",
        adsorptionRateCoefficient.value, adsorptionRateCoefficient.repr(), desorptionRateCoefficient.value,
        desorptionRateCoefficient.repr(), poreConcentrationOrder, capacityOrder, desorptionOrder, filmMassTransferCoefficient, poreDiffusivity);
    text += std::format("        pore/surface transport: {}\n", usePoreSurfaceTransport ? "enabled" : "disabled");
  }
  else if (type == Type::Elovich)
  {
    text += std::format(
        "        alpha:                 {:8.4e}  {}\n"
        "        beta:                  {:8.4e} [kg/mol]\n"
        "        film coefficient:      {:8.4e} [m/s]\n"
        "        pore diffusivity:      {:8.4e} [m^2/s]\n"
        "        pore/surface transport: {}\n",
        elovichAlpha.value, elovichAlpha.repr(), elovichBeta, filmMassTransferCoefficient, poreDiffusivity,
        usePoreSurfaceTransport ? "enabled" : "disabled");
  }
  if (isotherm.has_value())
  {
    text += "        equilibrium model:\n";
    text += isotherm->repr();
  }
  return text;
}


#include "chemisorption.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>
#include <stdexcept>

#include "utils.h"

namespace {
double signed_power(double force, double order, double regularizer)
{
  const double magnitude = std::abs(force);
  const double value = order < 1.0
      ? std::pow(regularizer, order) * std::expm1(order * std::log1p(magnitude / regularizer))
      : std::pow(magnitude, order);
  return ((force > 0.0) - (force < 0.0)) * value;
}
}

Chemisorption::Rates Chemisorption::evaluate(double equilibriumLoading, std::span<const double> loadings,
    double concentration, double temperature, double elapsedTime, double adsorptionAge, double desorptionAge) const
{
  if (loadings.size() != branchCount() || !(temperature > 0.0) || !std::isfinite(temperature))
    throw std::invalid_argument("Chemisorption kinetic state or temperature is invalid");
  Rates result;
  if (type == Type::None) return result;
  const double qeq = maximumLoading > 0.0 ? std::min(equilibriumLoading, maximumLoading) : equilibriumLoading;
  const double q = loadings[0];
  if (type == Type::General)
  {
    const double c = std::max(0.0, concentration), freeCapacity = std::max(0.0, qeq - std::max(q, 0.0));
    result.loading[0] = arrhenius(adsorptionRateCoefficient, adsorptionActivationEnergy, temperature) *
        std::pow(c, poreConcentrationOrder) * std::pow(freeCapacity, capacityOrder) -
        arrhenius(desorptionRateCoefficient, desorptionActivationEnergy, temperature) *
        std::pow(std::max(q, 0.0), desorptionOrder);
  }
  else if (type == Type::Elovich)
    result.loading[0] = elovichAlpha * std::max(0.0, concentration) * std::exp(-elovichBeta * std::max(0.0, q));
  else
  {
    for (std::size_t b = 0; b < branchCount(); ++b)
    {
      const double fraction = branchFraction(b), drive = fraction * qeq - loadings[b];
      if (fraction == 0.0) { result.loading[b] = 0.0; continue; }
      const bool adsorption = drive > 0.0, desorption = drive < 0.0;
      const double factor = (1.0 / referenceTemperature - 1.0 / temperature) / R;
      double ka = directionalRates ? adsorptionRateCoefficient * std::exp(adsorptionActivationEnergy * factor) : rateCoefficient;
      double kd = directionalRates ? desorptionRateCoefficient * std::exp(desorptionActivationEnergy * factor) : rateCoefficient;
      double na = type == Type::PseudoNth || type == Type::Avrami ? order : 1.0;
      double nd = directionalRates ? desorptionOrderValue : na;
      if (type == Type::DualNth)
      {
        const auto& branch = branches[b];
        ka = branch.adsorptionRate * std::exp(branch.adsorptionActivationEnergy * factor);
        kd = branch.desorptionRate * std::exp(branch.desorptionActivationEnergy * factor);
        na = branch.adsorptionOrder; nd = branch.desorptionOrder;
      }
      const double k = adsorption ? ka : kd, n = adsorption ? na : nd;
      if (type == Type::FirstOrder && desorption && desorptionMassTransferCoefficient > 0.0)
      {
        const double exponent = desorptionCoverageEnergy * std::max(0.0, q) / coverageReferenceLoading / (R * temperature);
        result.loading[b] = desorptionMassTransferCoefficient * std::exp(std::min(exponent, 700.0)) * drive;
      }
      else if (type == Type::Avrami && continuousHistory)
      {
        const double age = std::max(adsorption ? adsorptionAge : desorptionAge, 0.0);
        result.loading[b] = n * k * std::pow(age + ageRegularizer, n - 1.0) * drive;
        const double eraseAds = adsorptionMemoryErasure >= 0.0 ? adsorptionMemoryErasure : kd;
        const double eraseDes = desorptionMemoryErasure >= 0.0 ? desorptionMemoryErasure : ka;
        double activity = 1.0;
        if (historyDrivingRegularizer > 0.0) {
          const double ratio = std::abs(drive) / std::hypot(drive, historyDrivingRegularizer);
          activity = ratio * ratio;
        }
        result.adsorptionAge = activity * (adsorption ? ka : (desorption ? -eraseAds * adsorptionAge : 0.0));
        result.desorptionAge = activity * (desorption ? kd : (adsorption ? -eraseDes * desorptionAge : 0.0));
      }
      else if (type == Type::Avrami)
      {
        result.loading[b] = elapsedTime <= 0.0 ? (n == 1.0 ? k * drive : 0.0)
            : n * std::pow(k, n) * std::pow(elapsedTime, n - 1.0) * drive;
        result.adsorptionAge = 1.0;
      }
      else
      {
        const double ref = fractionalRates ? fraction * referenceLoading : 1.0;
        result.loading[b] = ref * k * signed_power(drive / ref, n, drivingRegularizer);
      }
    }
  }
  for (double value : result.loading)
    if (!std::isfinite(value)) throw std::runtime_error("Chemisorption rate overflow; check parameters and units");
  if (!std::isfinite(result.adsorptionAge) || !std::isfinite(result.desorptionAge))
    throw std::runtime_error("Chemisorption history rate overflow");
  return result;
}

double Chemisorption::rate(double equilibriumLoading, double loading, double concentration, double temperature,
                          double elapsedTime) const
{
  if (type == Type::DualNth || continuousHistory)
    throw std::invalid_argument("Stateful chemisorption requires the shared FDM/FVM column");
  return evaluate(equilibriumLoading, std::span<const double>(&loading, 1), concentration, temperature, elapsedTime).loading[0];
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
    case Type::DualNth:
      typeName = "DualNth";
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

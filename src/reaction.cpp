#include "reaction.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <stdexcept>
#include <string>
#include <vector>

#include "utils.h"

namespace
{
constexpr double tiny = 1.0e-300;

double finiteNonNegative(double value) noexcept
{
  if (!std::isfinite(value)) return 0.0;
  return std::max(0.0, value);
}

double productPower(std::span<const double> activity, const std::vector<size_t>& components,
                    const std::vector<double>& orders)
{
  double product = 1.0;
  for (size_t i = 0; i < components.size(); ++i)
  {
    const double order = i < orders.size() ? orders[i] : 1.0;
    if (order == 0.0) continue;
    product *= std::pow(finiteNonNegative(activity[components[i]]), order);
  }
  return std::isfinite(product) ? product : 0.0;
}

}  // namespace

double Reaction::equilibriumConstantAt(double temperature) const
{
  return std::max(tiny, equilibriumConstant(temperature));
}

double Reaction::forwardRateConstant(double temperature) const { return forwardRateCoefficient(temperature); }

double Reaction::rate(std::span<const double> activity, double temperature, std::span<double> coverage) const
{
  if (!enabled()) return 0.0;

  const double kf = forwardRateConstant(temperature);
  const double kb = kf / equilibriumConstantAt(temperature);

  switch (style)
  {
    case Style::GeneralPowerLaw:
      return kf * productPower(activity, reactants, forwardOrders) -
             kb * productPower(activity, products, backwardOrders);
    case Style::LangmuirHinshelwood:
    case Style::LangmuirHinshelwoodHougenWatson:
    {
      // Competitive coverages theta_j = K_j c_j / (1 + sum_i K_i c_i), one adsorption constant per component, in the
      // caller's work array.
      if (coverage.size() < activity.size())
      {
        throw std::runtime_error("Error: the reaction coverage work array is smaller than the activity vector");
      }
      const std::span<double> theta = coverage.first(activity.size());
      std::ranges::fill(theta, 0.0);
      double denominator = 1.0;
      for (size_t comp = 0; comp < std::min(adsorptionConstants.size(), activity.size()); ++comp)
      {
        theta[comp] = adsorptionConstants[comp](temperature) * finiteNonNegative(activity[comp]);
        denominator += theta[comp];
      }
      for (double& value : theta) value /= denominator;

      // LH is irreversible; LHHW subtracts the backward rate on the product coverages.
      const double forward = kf * productPower(theta, reactants, forwardOrders);
      if (style == Style::LangmuirHinshelwood) return forward;
      return forward - kb * productPower(theta, products, backwardOrders);
    }
  }
  return 0.0;
}

std::string Reaction::repr() const
{
  auto phaseName = [&]
  {
    switch (phase)
    {
      case Phase::Physisorbed:
        return "Physisorbed";
      case Phase::Chemisorbed:
        return "Chemisorbed";
      case Phase::PoreConcentration:
        return "PoreConcentration";
      case Phase::BulkConcentration:
        return "BulkConcentration";
    }
    return "Unknown";
  }();

  auto styleName = [&]
  {
    switch (style)
    {
      case Style::GeneralPowerLaw:
        return "GeneralPowerLaw";
      case Style::LangmuirHinshelwood:
        return "LangmuirHinshelwood";
      case Style::LangmuirHinshelwoodHougenWatson:
        return "LangmuirHinshelwoodHougenWatson";
    }
    return "Unknown";
  }();

  std::string text = std::format(
      "    reaction:\n"
      "        phase:                 {}\n"
      "        style:                 {}\n"
      "        forward k:             {:8.4e}  {}\n"
      "        equilibrium constant:  {:8.4e}  {}\n"
      "        heat of reaction:      {:8.4e} [J/mol]\n",
      phaseName, styleName, forwardRateCoefficient.value, forwardRateCoefficient.repr(), equilibriumConstant.value,
      equilibriumConstant.repr(), heatOfReaction);
  for (size_t comp = 0; comp < adsorptionConstants.size(); ++comp)
  {
    if (adsorptionConstants[comp].value == 0.0) continue;
    text += std::format("        {:<23}{:8.4e}  {}\n", std::format("adsorption K[{}]:", comp),
                        adsorptionConstants[comp].value, adsorptionConstants[comp].repr());
  }
  return text;
}

bool reactionsRequirePoreConcentration(const std::vector<Reaction>& reactions) noexcept
{
  return std::any_of(reactions.begin(), reactions.end(),
                     [](const Reaction& reaction)
                     { return reaction.phase == Reaction::Phase::PoreConcentration; });
}

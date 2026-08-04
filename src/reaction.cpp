#include "reaction.h"

#include <algorithm>
#include <cmath>
#include <format>
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

double boundedExp(double exponent)
{
  return std::exp(std::clamp(exponent, -700.0, 700.0));
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

std::vector<double> langmuirActivities(const Reaction& reaction, std::span<const double> activity,
                                       double temperature)
{
  const double K = reaction.equilibriumConstantAt(temperature);
  std::vector<double> theta(activity.size(), 0.0);
  for (size_t comp = 0; comp < activity.size(); ++comp)
  {
    const double Kc = K * finiteNonNegative(activity[comp]);
    theta[comp] = Kc / std::max(1.0 + Kc, tiny);
  }
  return theta;
}

std::vector<double> lhhwActivities(const Reaction& reaction, std::span<const double> activity,
                                   double temperature)
{
  const double K = reaction.equilibriumConstantAt(temperature);
  double denominator = 1.0;
  for (double value : activity)
  {
    denominator += K * finiteNonNegative(value);
  }
  denominator = std::max(denominator, tiny);

  std::vector<double> theta(activity.size(), 0.0);
  for (size_t comp = 0; comp < activity.size(); ++comp)
  {
    theta[comp] = K * finiteNonNegative(activity[comp]) / denominator;
  }
  return theta;
}

}  // namespace

double Reaction::stoichiometricCoefficient(size_t component) const noexcept
{
  double coefficient = 0.0;
  for (size_t i = 0; i < reactants.size(); ++i)
  {
    if (reactants[i] == component) coefficient -= reactantStoichiometry[i];
  }
  for (size_t i = 0; i < products.size(); ++i)
  {
    if (products[i] == component) coefficient += productStoichiometry[i];
  }
  return coefficient;
}

double Reaction::equilibriumConstantAt(double temperature) const
{
  const double exponent = -gibbsFreeEnergy / (R * std::max(temperature, 1.0e-10));
  return std::max(tiny, equilibriumConstant * boundedExp(exponent));
}

double Reaction::forwardRateConstant(double temperature) const
{
  return arrhenius(forwardRateCoefficient, forwardActivationEnergy, temperature);
}

double Reaction::rate(std::span<const double> activity, double temperature) const
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
    {
      const std::vector<double> theta = langmuirActivities(*this, activity, temperature);
      return kf * productPower(theta, reactants, forwardOrders) -
             kb * productPower(theta, products, backwardOrders);
    }
    case Style::LangmuirHinshelwoodHougenWatson:
    {
      const std::vector<double> theta = lhhwActivities(*this, activity, temperature);
      return kf * productPower(theta, reactants, forwardOrders) -
             kb * productPower(theta, products, backwardOrders);
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

  return std::format(
      "    reaction:\n"
      "        phase:                 {}\n"
      "        style:                 {}\n"
      "        forward k0:            {:8.4e}\n"
      "        forward Ea:            {:8.4e} [J/mol]\n"
      "        equilibrium constant:  {:8.4e}\n"
      "        delta Gibbs energy:    {:8.4e} [J/mol]\n",
      phaseName, styleName, forwardRateCoefficient, forwardActivationEnergy,
      equilibriumConstant, gibbsFreeEnergy);
}

bool reactionsRequirePoreConcentration(const std::vector<Reaction>& reactions) noexcept
{
  return std::any_of(reactions.begin(), reactions.end(),
                     [](const Reaction& reaction)
                     { return reaction.phase == Reaction::Phase::PoreConcentration; });
}

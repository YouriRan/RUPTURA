#include "temperature_law.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <format>
#include <utility>

#include "utils.h"

namespace
{
constexpr std::array<std::pair<TemperatureLaw::Type, std::string_view>, 6> lawNames{{
    {TemperatureLaw::Type::Constant, "Constant"},
    {TemperatureLaw::Type::VantHoff, "VantHoff"},
    {TemperatureLaw::Type::Arrhenius, "Arrhenius"},
    {TemperatureLaw::Type::Exponential, "Exponential"},
    {TemperatureLaw::Type::Linear, "Linear"},
    {TemperatureLaw::Type::Inverse, "Inverse"},
}};

double boundedExp(double exponent) { return std::exp(std::clamp(exponent, -700.0, 700.0)); }
}  // namespace

double TemperatureLaw::operator()(double temperature) const
{
  switch (type)
  {
    case Type::Constant:
      return value;
    case Type::VantHoff:
    case Type::Arrhenius:
    {
      // 1/(RT) - 1/(R T_ref); without T_ref the second term is zero (the uncentred form).
      const double inverseTemperature = 1.0 / (R * std::max(temperature, 1.0e-10));
      const double inverseReferenceTemperature =
          referenceTemperature.has_value() ? 1.0 / (R * std::max(*referenceTemperature, 1.0e-10)) : 0.0;
      const double energy = type == Type::VantHoff ? coefficient : -coefficient;
      return value * boundedExp(energy * (inverseTemperature - inverseReferenceTemperature));
    }
    case Type::Exponential:
      return value * boundedExp(coefficient * (1.0 - temperature / referenceTemperature.value_or(temperature)));
    case Type::Linear:
      return value + coefficient * (temperature - referenceTemperature.value_or(temperature));
    case Type::Inverse:
      return value +
             coefficient * (1.0 - referenceTemperature.value_or(temperature) / std::max(temperature, 1.0e-10));
  }
  return value;
}

std::string TemperatureLaw::repr() const
{
  if (type == Type::Constant) return "constant";

  std::string s(name(type));
  switch (type)
  {
    case Type::VantHoff:
      s += std::format(", Q = {:8.4e} J/mol", coefficient);
      break;
    case Type::Arrhenius:
      s += std::format(", E = {:8.4e} J/mol", coefficient);
      break;
    case Type::Exponential:
      s += std::format(", chi = {:8.4e}", coefficient);
      break;
    case Type::Linear:
      s += std::format(", slope = {:8.4e} /K", coefficient);
      break;
    case Type::Inverse:
      s += std::format(", alpha = {:8.4e}", coefficient);
      break;
    case Type::Constant:
      break;
  }
  if (energyInherited) s += " (heat of adsorption)";
  s += referenceTemperature.has_value() ? std::format(", T_ref = {:8.4e} K", *referenceTemperature)
                                        : std::string(", uncentred");
  return s;
}

std::string_view TemperatureLaw::name(Type type) noexcept
{
  for (const auto& [lawType, lawName] : lawNames)
  {
    if (lawType == type) return lawName;
  }
  return "Constant";
}

std::string_view TemperatureLaw::coefficientKey(Type type) noexcept
{
  switch (type)
  {
    case Type::Constant:
      return "";
    case Type::VantHoff:
      return "Energy";
    case Type::Arrhenius:
      return "ActivationEnergy";
    case Type::Exponential:
      return "Chi";
    case Type::Linear:
      return "Slope";
    case Type::Inverse:
      return "Alpha";
  }
  return "";
}

std::optional<TemperatureLaw::Type> TemperatureLaw::fromName(std::string_view lawName) noexcept
{
  for (const auto& [lawType, candidate] : lawNames)
  {
    if (candidate.size() == lawName.size() &&
        std::equal(candidate.begin(), candidate.end(), lawName.begin(),
                   [](unsigned char a, unsigned char b) { return std::tolower(a) == std::tolower(b); }))
    {
      return lawType;
    }
  }
  return std::nullopt;
}

std::string TemperatureLaw::describeMask(unsigned allowedLaws)
{
  std::string s;
  for (const auto& [lawType, lawName] : lawNames)
  {
    if (lawType == Type::Constant || (allowedLaws & mask(lawType)) == 0) continue;
    if (!s.empty()) s += " or ";
    s += lawName;
  }
  return s;
}

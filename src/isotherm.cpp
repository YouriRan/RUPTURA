#include "isotherm.h"

#include <array>
#include <cstdlib>
#include <print>
#include <stdexcept>
#include <utility>

Isotherm::Isotherm(Isotherm::Type t, const std::vector<double>& parameters,
                   std::vector<std::pair<size_t, TemperatureLaw>> temperatureLaws)
    : type(t),
      parameters(std::move(parameters)),
      numberOfParameters(parameters.size()),
      temperatureLaws(std::move(temperatureLaws))
{
  for (const auto& [index, law] : this->temperatureLaws)
  {
    if (index >= this->parameters.size())
    {
      throw std::out_of_range("Error: temperature law for isotherm parameter " + std::to_string(index) +
                              ", which does not exist");
    }
    this->parameters[index] = law.value;
  }
}
Isotherm::Isotherm(size_t t, const std::vector<double>& parameters)
    : type(Isotherm::Type(t)), parameters(std::move(parameters)), numberOfParameters(parameters.size())
{
}

const TemperatureLaw* Isotherm::temperatureLaw(size_t index) const noexcept
{
  for (const auto& [lawIndex, law] : temperatureLaws)
  {
    if (lawIndex == index) return &law;
  }
  return nullptr;
}

std::span<const TemperatureParameterSpec> Isotherm::temperatureParameters(Type type) noexcept
{
  using Law = TemperatureLaw::Type;
  constexpr unsigned vantHoff = TemperatureLaw::mask(Law::VantHoff);
  constexpr unsigned capacity = TemperatureLaw::mask(Law::Linear) | TemperatureLaw::mask(Law::Exponential);
  constexpr unsigned exponent = TemperatureLaw::mask(Law::Linear) | TemperatureLaw::mask(Law::Inverse);

  // One entry per Parameters entry. The heat-linked affinity inherits the heat of adsorption in a VantHoff law
  // without "Energy"; Redlich-Peterson and GAB have two affinity-like parameters and need explicit energies.
  static constexpr std::array<TemperatureParameterSpec, 1> henry{{{"Henry coefficient", vantHoff, true}}};
  static constexpr std::array<TemperatureParameterSpec, 2> freundlich{
      {{"Freundlich coefficient", vantHoff, true}, {"exponent", exponent, false}}};
  static constexpr std::array<TemperatureParameterSpec, 2> langmuir{
      {{"saturation loading", capacity, false}, {"affinity", vantHoff, true}}};
  static constexpr std::array<TemperatureParameterSpec, 3> langmuirPH{
      {{"saturation loading", capacity, false}, {"affinity", vantHoff, true}, {"pH_0", 0, false}}};
  static constexpr std::array<TemperatureParameterSpec, 3> heterogeneous{
      {{"saturation loading", capacity, false}, {"affinity", vantHoff, true}, {"exponent", exponent, false}}};
  static constexpr std::array<TemperatureParameterSpec, 3> redlichPeterson{
      {{"linear constant", vantHoff, false}, {"affinity", vantHoff, false}, {"exponent", exponent, false}}};
  static constexpr std::array<TemperatureParameterSpec, 3> gab{
      {{"monolayer loading", capacity, false}, {"C", vantHoff, false}, {"K", vantHoff, false}}};
  // The humidity coefficients are constant.
  static constexpr std::array<TemperatureParameterSpec, 4> humidityLangmuir{{{"saturation loading", capacity, false},
                                                                             {"affinity", vantHoff, true},
                                                                             {"humidity enhancement", 0, false},
                                                                             {"water affinity", 0, false}}};
  static constexpr std::array<TemperatureParameterSpec, 3> moistureSwingLangmuir{
      {{"saturation loading", capacity, false}, {"affinity", vantHoff, true}, {"humidity exponent", 0, false}}};
  static constexpr std::array<TemperatureParameterSpec, 4> gabLangmuir{{{"saturation loading", capacity, false},
                                                                        {"affinity", vantHoff, true},
                                                                        {"water enhancement", 0, false},
                                                                        {"water blocking", 0, false}}};
  static constexpr std::array<TemperatureParameterSpec, 5> gabHeterogeneous{{{"saturation loading", capacity, false},
                                                                             {"affinity", vantHoff, true},
                                                                             {"exponent", exponent, false},
                                                                             {"water enhancement", 0, false},
                                                                             {"water blocking", 0, false}}};
  static constexpr std::array<TemperatureParameterSpec, 5> gabModifiedToth{{{"saturation loading", capacity, false},
                                                                            {"affinity", vantHoff, true},
                                                                            {"exponent", exponent, false},
                                                                            {"water capacity coefficient", 0, false},
                                                                            {"water affinity coefficient", 0, false}}};
  static constexpr std::array<TemperatureParameterSpec, 4> exponentialHumidityToth{
      {{"saturation loading", capacity, false},
       {"affinity", vantHoff, true},
       {"exponent", exponent, false},
       {"humidity capacity coefficient", 0, false}}};

  switch (type)
  {
    case Type::Henry:
      return henry;
    case Type::Freundlich:
      return freundlich;
    case Type::Langmuir:
      return langmuir;
    case Type::Langmuir_pH:
      return langmuirPH;
    case Type::Sips:
    case Type::Langmuir_Freundlich:
    case Type::Toth:
      return heterogeneous;
    case Type::Redlich_Peterson:
      return redlichPeterson;
    case Type::GAB:
      return gab;
    case Type::Langmuir_humidity:
      return humidityLangmuir;
    case Type::Langmuir_moisture_swing:
      return moistureSwingLangmuir;
    case Type::GAB_Langmuir:
      return gabLangmuir;
    case Type::GAB_Sips:
    case Type::GAB_Toth:
      return gabHeterogeneous;
    case Type::GAB_modified_Toth:
      return gabModifiedToth;
    case Type::Toth_exponential_humidity:
      return exponentialHumidityToth;
    default:
      return {};
  }
}

void Isotherm::print() const { std::print("{}", repr()); }

std::string Isotherm::repr() const
{
  std::string s;
  auto out = std::back_inserter(s);

  // A temperature-dependent parameter prints its value at the reference temperature and its law.
  auto appendParameter = [&](std::string_view name, size_t index)
  {
    const TemperatureLaw* law = temperatureLaw(index);
    if (law == nullptr)
    {
      std::format_to(out, "        {:<7}{:8.4e}\n", name, parameters[index]);
    }
    else
    {
      std::format_to(out, "        {:<7}{:8.4e}  {}\n", name, law->value, law->repr());
    }
  };

  switch (type)
  {
    case Isotherm::Type::Langmuir:
    {
      s += "    Langmuir isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      break;
    }
    case Isotherm::Type::Langmuir_pH:
    {
      s += "    pH-dependent Langmuir isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("k:", 1);
      appendParameter("pH_0:", 2);
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      s += "    Anti-Langmuir isotherm\n";
      appendParameter("a:", 0);
      appendParameter("b:", 1);
      break;
    }
    case Isotherm::Type::BET:
    {
      s += "    BET isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("c:", 2);
      break;
    }
    case Isotherm::Type::Henry:
    {
      s += "    Henry isotherm\n";
      appendParameter("a:", 0);
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      s += "    Freundlich isotherm\n";
      appendParameter("a:", 0);
      appendParameter("nu:", 1);
      break;
    }
    case Isotherm::Type::Sips:
    {
      s += "    Sips isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      s += "    Langmuir-Freundlich isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      s += "    Redlich-Peterson isotherm\n";
      appendParameter("a:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      break;
    }
    case Isotherm::Type::Toth:
    {
      s += "    Toth isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      break;
    }
    case Isotherm::Type::Unilan:
    {
      s += "    Unilan isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("eta:", 2);
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      s += "    O'Brian & Myers isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("sigma:", 2);
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      s += "    Quadratic isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("c:", 2);
      break;
    }
    case Isotherm::Type::Temkin:
    {
      s += "    Temkin isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("c:", 2);
      break;
    }
    case Isotherm::Type::BingelWalton:
    {
      s += "    Bingel&Walton isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("a:", 1);
      appendParameter("b:", 2);
      break;
    }
    case Isotherm::Type::GAB:
    {
      s += "    Guggenheim-Anderson-de Boer (GAB) isotherm\n";
      appendParameter("q_mono:", 0);
      appendParameter("C:", 1);
      appendParameter("K:", 2);
      break;
    }
    case Isotherm::Type::Langmuir_humidity:
    {
      s += "    Humidity-enhanced Langmuir isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("alpha:", 2);
      appendParameter("k_w:", 3);
      break;
    }
    case Isotherm::Type::Langmuir_moisture_swing:
    {
      s += "    Moisture-swing Langmuir isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("n_w:", 2);
      break;
    }
    case Isotherm::Type::GAB_Langmuir:
    {
      s += "    GAB-coadsorption Langmuir isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("alpha:", 2);
      appendParameter("omega:", 3);
      break;
    }
    case Isotherm::Type::GAB_Sips:
    case Isotherm::Type::GAB_Toth:
    {
      s += type == Isotherm::Type::GAB_Sips ? "    GAB-coadsorption Sips isotherm\n"
                                             : "    GAB-coadsorption Toth isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      appendParameter("alpha:", 3);
      appendParameter("omega:", 4);
      break;
    }
    case Isotherm::Type::GAB_modified_Toth:
    {
      s += "    GAB-modified Toth isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      appendParameter("gamma:", 3);
      appendParameter("beta:", 4);
      break;
    }
    case Isotherm::Type::Toth_exponential_humidity:
    {
      s += "    Exponential-humidity Toth isotherm\n";
      appendParameter("q_sat:", 0);
      appendParameter("b:", 1);
      appendParameter("nu:", 2);
      appendParameter("a:", 3);
      break;
    }
    default:
      break;
  }

  return s;
}

bool Isotherm::isUnphysical() const
{
  switch (type)
  {
    case Isotherm::Type::Langmuir:
    {
      if (parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10) return true;
      return false;
    }
    case Isotherm::Type::Langmuir_pH:
    {
      if (parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ||
          !std::isfinite(parameters[2])) return true;
      return false;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      if (parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10) return true;
      return false;
    }
    case Isotherm::Type::BET:
    {
      return false;
    }
    case Isotherm::Type::Henry:
    {
      if (parameters[0] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Freundlich:
    {
      if (parameters[0] < 0.0 || parameters[1] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Sips:
    {
      if (parameters[0] < 0 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ||
          parameters[2] < 0.0 || parameters[2] > 100.0)
        return true;
      return false;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      if (parameters[0] < 1.0e-20 || parameters[0] > 1.0e20 || parameters[1] < 0.0 || parameters[1] > 1.0e10 ||
          parameters[2] < 0.0 || parameters[2] > 100.0)
        return true;
      return false;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      if (parameters[0] < 0.0 || parameters[1] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Toth:
    {
      if (parameters[0] < 0 || parameters[1] < 0.0 || parameters[2] < 0.0 || parameters[2] > 100.0) return true;
      return false;
    }
    case Isotherm::Type::Unilan:
    {
      if (parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      if (parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Quadratic:
    {
      if (parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::Temkin:
    {
      if (parameters[0] <= 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
      return false;
    }
    case Isotherm::Type::BingelWalton:
    {
      if (parameters[0] <= 0.0 || (parameters[1] + parameters[2]) < 1e-3) return true;
      return false;
    }
    case Isotherm::Type::GAB:
    {
      return parameters.size() != 3 || parameters[0] <= 0.0 || parameters[1] <= 0.0 || parameters[2] <= 0.0;
    }
    case Isotherm::Type::Langmuir_humidity:
    case Isotherm::Type::Langmuir_moisture_swing:
    case Isotherm::Type::GAB_Langmuir:
    case Isotherm::Type::GAB_Sips:
    case Isotherm::Type::GAB_Toth:
    case Isotherm::Type::GAB_modified_Toth:
    case Isotherm::Type::Toth_exponential_humidity:
    {
      return parameters[0] < 0.0 || std::any_of(parameters.begin(), parameters.end(),
                                                 [](double parameter) { return !std::isfinite(parameter); });
    }
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

void Isotherm::randomize(double maximumLoading)
{
  // affinity-like parameters are drawn as U^(20 (1 - U')), which lies in (0, 1] with half of the draws between 1e-7
  // and 1e-1 (b in 1/Pa)
  switch (type)
  {
    case Isotherm::Type::Langmuir:
    case Isotherm::Type::Langmuir_pH:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      break;
    }
    case Isotherm::Type::BET:
    {
      parameters[0] = 10.0 * RandomNumber::Uniform();
      parameters[1] = 10.0 * RandomNumber::Uniform();
      parameters[2] = 10.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Henry:
    {
      parameters[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      parameters[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[1] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Sips:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Toth:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Unilan:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      parameters[0] = 2.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = 10.0 * RandomNumber::Uniform();
      parameters[2] = 10.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Temkin:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (1.0 - RandomNumber::Uniform()));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::BingelWalton:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = 0.1 + 2.0 * RandomNumber::Uniform();
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::GAB:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = 0.01 + 2.0 * RandomNumber::Uniform();
      parameters[2] = 0.5 + 0.49 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Langmuir_humidity:
    case Isotherm::Type::Langmuir_moisture_swing:
    case Isotherm::Type::GAB_Langmuir:
    case Isotherm::Type::GAB_Sips:
    case Isotherm::Type::GAB_Toth:
    case Isotherm::Type::GAB_modified_Toth:
    case Isotherm::Type::Toth_exponential_humidity:
      throw std::runtime_error("Error: humidity-dependent isotherms cannot be fitted to pure-component data");
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

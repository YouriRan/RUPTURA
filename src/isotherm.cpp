#include "isotherm.h"

#include <cstdlib>

Isotherm::Isotherm(Isotherm::Type t, const std::vector<double>& values, size_t numberOfValues, bool nonIsothermal)
    : type(t), parameters(values), numberOfParameters(numberOfValues), nonIsothermal(nonIsothermal)
{
}
Isotherm::Isotherm(size_t t, const std::vector<double>& values, size_t numberOfValues, bool nonIsothermal)
    : type(Isotherm::Type(t)), parameters(values), numberOfParameters(numberOfValues), nonIsothermal(nonIsothermal)
{
}

void Isotherm::print() const { std::cout << repr(); }

std::string Isotherm::repr() const
{
  std::string s;
  auto out = std::back_inserter(s);

  auto appendParameter = [&](std::string_view name, double value)
  { std::format_to(out, "        {:<7}{:8.4e}\n", name, value); };

  switch (type)
  {
    case Isotherm::Type::Langmuir:
    {
      s += nonIsothermal ? "    Langmuir isotherm (non-isothermal)\n" : "    Langmuir isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter(nonIsothermal ? "b_0:" : "b:", parameters[1]);
      if (nonIsothermal) appendParameter("ΔH:", parameters[2]);
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      s += "    Anti-Langmuir isotherm\n";
      appendParameter("a:", parameters[0]);
      appendParameter("b:", parameters[1]);
      break;
    }
    case Isotherm::Type::BET:
    {
      s += "    BET isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("c:", parameters[2]);
      break;
    }
    case Isotherm::Type::Henry:
    {
      s += "    Henry isotherm\n";
      appendParameter("a:", parameters[0]);
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      s += "    Freundlich isotherm\n";
      appendParameter("a:", parameters[0]);
      appendParameter("nu:", parameters[1]);
      break;
    }
    case Isotherm::Type::Sips:
    {
      s += nonIsothermal ? "    Sips isotherm (non-isothermal)\n" : "    Sips isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("nu:", parameters[2]);
      if (nonIsothermal) appendParameter("ΔH:", parameters[3]);
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      s += nonIsothermal ? "    Langmuir-Freundlich isotherm (non-isothermal)\n" : "    Langmuir-Freundlich isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter(nonIsothermal ? "b_0:" : "b:", parameters[1]);
      appendParameter("nu:", parameters[2]);
      if (nonIsothermal) appendParameter("ΔH:", parameters[3]);
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      s += "    Redlich-Peterson isotherm\n";
      appendParameter("a:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("nu:", parameters[2]);
      break;
    }
    case Isotherm::Type::Toth:
    {
      s += "    Toth isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("nu:", parameters[2]);
      break;
    }
    case Isotherm::Type::Unilan:
    {
      s += "    Unilan isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("eta:", parameters[2]);
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      s += "    O'Brian & Myers isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("sigma:", parameters[2]);
      break;
    }
    case Isotherm::Type::Quadratic:
    {
      s += "    Quadratic isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("c:", parameters[2]);
      break;
    }
    case Isotherm::Type::Temkin:
    {
      s += "    Temkin isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("b:", parameters[1]);
      appendParameter("c:", parameters[2]);
      break;
    }
    case Isotherm::Type::BingelWalton:
    {
      s += "    Bingel&Walton isotherm\n";
      appendParameter("q_sat:", parameters[0]);
      appendParameter("a:", parameters[1]);
      appendParameter("b:", parameters[2]);
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
      if (parameters[0] < 0.0 || parameters[1] < 0.0 || parameters[2] < 0.0) return true;
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
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

void Isotherm::randomize(double maximumLoading)
{
  switch (type)
  {
    case Isotherm::Type::Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Anti_Langmuir:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
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
      parameters[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      break;
    }
    case Isotherm::Type::Freundlich:
    {
      parameters[0] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[1] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Sips:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Langmuir_Freundlich:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Redlich_Peterson:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Toth:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::Unilan:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
      parameters[2] = 0.1 + 2.0 * RandomNumber::Uniform();
      break;
    }
    case Isotherm::Type::OBrien_Myers:
    {
      parameters[0] = 1.1 * maximumLoading * RandomNumber::Uniform();
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
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
      parameters[1] = std::pow(RandomNumber::Uniform(), 10.0 * 2.0 * (RandomNumber::Uniform() - 1.0));
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
    default:
      throw std::runtime_error("Error: unkown isotherm type");
  }
}

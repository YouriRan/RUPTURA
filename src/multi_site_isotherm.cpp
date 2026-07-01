#include "multi_site_isotherm.h"

#include <cmath>
#include <print>
#include <sstream>

#include "special_functions.h"

MultiSiteIsotherm::MultiSiteIsotherm(std::vector<Isotherm> sites)
{
  numberOfSites = sites.size();
  numberOfParameters = 0;

  for (const Isotherm& site : sites)
  {
    add(site);
  }
}

void MultiSiteIsotherm::print() const { std::print("{}", repr()); }

std::string MultiSiteIsotherm::repr() const
{
  std::string s;
  s += "    number of isotherm sites:  " + std::to_string(numberOfSites) + "\n";
  for (size_t i = 0; i < numberOfSites; ++i)
  {
    s += sites[i].repr();
  }
  return s;
}

void MultiSiteIsotherm::add(const Isotherm& isotherm)
{
  siteParameterIndex.push_back(numberOfParameters);
  sites.push_back(isotherm);

  numberOfParameters += isotherm.numberOfParameters;
  for (size_t i = 0; i < isotherm.numberOfParameters; ++i)
  {
    parameterIndices.emplace_back(sites.size() - 1, i);
  }
}

void MultiSiteIsotherm::setParameters(std::vector<double> params)
{
  for (size_t i = 0; i < params.size(); ++i)
  {
    std::pair<size_t, size_t> index = parameterIndices[i];
    sites[index.first].parameters[index.second] = params[i];
  }
}

std::vector<double> MultiSiteIsotherm::getParameters()
{
  std::vector<double> params;
  for (size_t i = 0; i < numberOfParameters; ++i)
  {
    std::pair<size_t, size_t> index = parameterIndices[i];
    params.push_back(sites[index.first].parameters[index.second]);
  }
  return params;
}

// returns the inverse-pressure (1/P) that corresponds to the given reduced_grand_potential reducedGrandPotential
// advantage: for isotherms with zero equilibrium constant the result would be infinite, but the inverse is zero
double MultiSiteIsotherm::inversePressureForPsi(double reduced_grand_potential, double& cachedP0, double scale) const
{
  const double tiny = 1.0e-15;

  double left_bracket;
  double right_bracket;

  // For a single Langmuir or Langmuir-Freundlich site, the inverse can be handled analytically
  if (numberOfSites == 1)
  {
    return sites[0].inversePressureForPsi(reduced_grand_potential, cachedP0, scale);
  }

  // from here on, work with pressure, and return 1.0 / pressure at the end of the routine
  double p_start;
  if (cachedP0 <= 0.0)
  {
    p_start = 5.0;
  }
  else
  {
    // use the last value of Pi0
    p_start = cachedP0;
  }

  // use bisection algorithm
  double s = psiForPressure(p_start, scale);

  size_t nr_steps = 0;
  left_bracket = p_start;
  right_bracket = p_start;

  if (s < reduced_grand_potential)
  {
    // find the bracket on the right
    do
    {
      right_bracket *= 2.0;
      s = psiForPressure(right_bracket, scale);

      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("reduced_grand_potential: {}\n", reduced_grand_potential);
        std::print("reducedGrandPotential: {}\n", s);
        std::print("p_start: {}\n", p_start);
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (s < reduced_grand_potential);
  }
  else
  {
    // find the bracket on the left
    do
    {
      left_bracket *= 0.5;
      s = psiForPressure(left_bracket, scale);

      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("reduced_grand_potential: {}\n", reduced_grand_potential);
        std::print("reducedGrandPotential: {}\n", s);
        std::print("p_start: {}\n", p_start);
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (s > reduced_grand_potential);
  }

  do
  {
    double middle = 0.5 * (left_bracket + right_bracket);
    s = psiForPressure(middle, scale);

    if (s > reduced_grand_potential)
      right_bracket = middle;
    else
      left_bracket = middle;

    ++nr_steps;
    if (nr_steps > 100000)
    {
      std::print("Left bracket: {}\n", left_bracket);
      std::print("Right bracket: {}\n", right_bracket);
      throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);

  double middle = 0.5 * (left_bracket + right_bracket);

  //  Store the last value of Pi0
  cachedP0 = middle;

  return 1.0 / middle;
}

double MultiSiteIsotherm::fitness() const
{
  const double penaltyCost = 50.0;
  for (size_t i = 0; i < numberOfSites; ++i)
  {
    if (sites[i].isUnphysical()) return penaltyCost;
  }
  return 0.0;
}

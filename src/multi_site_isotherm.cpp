#include "multi_site_isotherm.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <print>
#include <sstream>

#include "special_functions.h"

MultiSiteIsotherm::MultiSiteIsotherm(std::vector<Isotherm> sites)
{
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
  s += "    number of isotherm sites:  " + std::to_string(sites.size()) + "\n";
  for (const Isotherm& site : sites)
  {
    s += site.repr();
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

std::vector<double> MultiSiteIsotherm::getParameters() const
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
double MultiSiteIsotherm::inversePressureForPsi(double reduced_grand_potential, double& cachedP0) const
{
  const auto firstActiveSite =
      std::find_if(sites.begin(), sites.end(), [](const Isotherm& site) { return site.enabled(); });
  if (firstActiveSite == sites.end()) return 0.0;

  // A multisite isotherm can contain retained, inactive sites for SIAST index alignment. If only one site is active,
  // its inverse can still be handled directly.
  const auto secondActiveSite =
      std::find_if(std::next(firstActiveSite), sites.end(), [](const Isotherm& site) { return site.enabled(); });
  if (secondActiveSite == sites.end())
  {
    return firstActiveSite->inversePressureForPsi(reduced_grand_potential, cachedP0);
  }

  return inversePressureByBisection([this](double pressure) { return psiForPressure(pressure); },
                                    reduced_grand_potential, cachedP0);
}

double MultiSiteIsotherm::fitness() const
{
  const double penaltyCost = 50.0;
  for (const Isotherm& site : sites)
  {
    if (site.isUnphysical()) return penaltyCost;
  }
  return 0.0;
}

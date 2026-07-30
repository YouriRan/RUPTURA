#include "multi_site_chemisorption.h"

#include <algorithm>
#include <format>
#include <utility>

MultiSiteChemisorption::MultiSiteChemisorption(std::vector<Chemisorption> sites) : sites(std::move(sites))
{
  numberOfSites = this->sites.size();
}

void MultiSiteChemisorption::add(const Chemisorption& chemisorption)
{
  sites.push_back(chemisorption);
  numberOfSites = sites.size();
}

bool MultiSiteChemisorption::enabled() const noexcept
{
  return std::any_of(sites.begin(), sites.end(), [](const Chemisorption& site) { return site.enabled(); });
}

bool MultiSiteChemisorption::usesSurfacePoreTransport() const noexcept
{
  return std::any_of(sites.begin(), sites.end(),
                     [](const Chemisorption& site) { return site.usesSurfacePoreTransport(); });
}

double MultiSiteChemisorption::maximumLoading() const noexcept
{
  double maximum = 0.0;
  for (const Chemisorption& site : sites)
  {
    if (site.enabled()) maximum += std::max(0.0, site.maximumLoading);
  }
  return maximum;
}

double MultiSiteChemisorption::equilibriumLoading(size_t site, double totalEquilibriumLoading) const noexcept
{
  if (site >= sites.size() || !sites[site].usesEquilibriumLoading()) return 0.0;

  double totalCapacity = 0.0;
  size_t equilibriumSites = 0;
  bool allSitesHaveCapacity = true;
  for (const Chemisorption& candidate : sites)
  {
    if (!candidate.usesEquilibriumLoading()) continue;
    ++equilibriumSites;
    if (candidate.maximumLoading > 0.0)
    {
      totalCapacity += candidate.maximumLoading;
    }
    else
    {
      allSitesHaveCapacity = false;
    }
  }

  if (allSitesHaveCapacity && totalCapacity > 0.0)
  {
    return totalEquilibriumLoading * sites[site].maximumLoading / totalCapacity;
  }
  return equilibriumSites == 0 ? 0.0 : totalEquilibriumLoading / static_cast<double>(equilibriumSites);
}

std::string MultiSiteChemisorption::repr() const
{
  std::string text = std::format("    number of chemisorption sites: {}\n", numberOfSites);
  for (size_t site = 0; site < sites.size(); ++site)
  {
    text += std::format("    chemisorption site {}:\n", site);
    text += sites[site].repr();
  }
  return text;
}

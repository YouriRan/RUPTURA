#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "chemisorption.h"

/**
 * \brief Collection of independent chemisorption kinetic sites.
 */
struct MultiSiteChemisorption
{
  MultiSiteChemisorption() noexcept = default;
  MultiSiteChemisorption(std::vector<Chemisorption> sites);

  size_t numberOfSites{0};
  std::vector<Chemisorption> sites{};

  void add(const Chemisorption& chemisorption);

  [[nodiscard]] bool enabled() const noexcept;
  [[nodiscard]] bool usesSurfacePoreTransport() const noexcept;
  [[nodiscard]] double maximumLoading() const noexcept;

  /**
   * \brief Partitions a component equilibrium loading over kinetic sites.
   *
   * Only equilibrium-driven kinetic models participate. Positive capacities
   * are used as weights when every participating site defines one; otherwise
   * the equilibrium loading is divided equally.
   */
  [[nodiscard]] double equilibriumLoading(size_t site, double totalEquilibriumLoading) const noexcept;

  [[nodiscard]] std::string repr() const;
};

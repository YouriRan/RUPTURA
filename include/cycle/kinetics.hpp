#pragma once
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "cycle/types.hpp"
namespace cycle
{
// Arrays use branch/component/node-major ordering, with node varying fastest.
class ChemistryKinetics
{
 public:
  ChemistryKinetics(const Json& config, std::size_t n_components);
  std::size_t n_components() const noexcept { return nc_; }
  std::size_t n_branches() const noexcept { return nb_; }
  std::size_t n_loading_blocks() const noexcept { return nb_; }
  std::size_t n_history_blocks() const noexcept { return nh_; }
  const std::string& model() const noexcept { return model_; }
  const Vec& fractions() const noexcept { return fractions_; }
  const std::optional<Vec>& heat_release() const noexcept { return heat_; }
  // Import-only compatibility adapter. Rates delegate to ::Chemisorption.
  Json canonical_site(std::size_t component) const;
  std::optional<Vec> equilibrium(const Vec& y, const Vec& pressure, const Vec& temperature) const;
  Vec split_loading(const Vec& total, std::size_t nz) const;
  Vec initial_history(std::size_t nz) const;
  std::pair<Vec, Vec> evaluate(const Vec& qeq, const Vec& loads, const Vec& history, const Vec& temperature) const;
  std::optional<Vec> heat_rate(const Vec& rates, std::size_t nz,
                               const std::optional<Vec>& fallback = std::nullopt) const;

 private:
  std::size_t nc_, nb_ = 1, nh_ = 0;
  std::string model_, basis_, eq_model_;
  double tref_ = 298.15, drive_eps_ = 1e-12, age_eps_ = 1e-6, eq_tref_ = 298.15;
  Vec ka_, kd_, eaa_, ead_, fractions_, oa_, od_, na_, nd_, eq_heat_, eq_k_, eq_qsat_, eq_b_;
  std::optional<Vec> qref_, heat_, erase_a_, erase_d_;
  bool has_eq_ = false;
};
using ChemistryModel = ChemistryKinetics;
}  // namespace cycle

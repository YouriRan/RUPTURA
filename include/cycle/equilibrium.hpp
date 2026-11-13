#pragma once
#include <memory>
#include <string>
#include <vector>

#include "cycle/types.hpp"
namespace cycle
{
class MaterialEquilibrium
{
 public:
  MaterialEquilibrium(const Json&, double R, std::size_t ncomp, std::size_t nz, const std::string& phase = "phy");
  ~MaterialEquilibrium();
  MaterialEquilibrium(const MaterialEquilibrium&);
  MaterialEquilibrium& operator=(const MaterialEquilibrium&);
  MaterialEquilibrium(MaterialEquilibrium&&) noexcept;
  MaterialEquilibrium& operator=(MaterialEquilibrium&&) noexcept;
  // Component-major partial pressures and returned mol/kg; Ts is node-major.
  Vec loading(const Vec& partial_pressures, const Vec& Ts) const;
  const Vec& qmax_effective() const;
  const std::vector<bool>& bounded_mask() const;
  Json cache_diagnostics() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
double water_saturation_pressure(double temperature, const Json& config);
double relative_humidity(double water_partial_pressure, double temperature, const Json& config);
}  // namespace cycle

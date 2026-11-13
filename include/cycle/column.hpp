#pragma once
#include <memory>
#include <string>
#include <vector>

#include "cycle/types.hpp"
namespace cycle
{
struct Derived
{
  Vec c, q, Tg, Ts, Tw, Tj, ctot, y, P, rho;
};
struct StateLayout
{
  std::size_t n_components = 0, n_axial = 0, n_loading_groups = 0, n_history_groups = 0;
  bool transient_htf = false;
  std::size_t n_blocks() const { return n_components * (1 + n_loading_groups + n_history_groups) + 3 + transient_htf; }
  std::size_t size() const { return n_blocks() * n_axial; }
};
class Column
{
 public:
  explicit Column(const Json& cfg);
  ~Column();
  Column(Column&&) noexcept;
  Column& operator=(Column&&) noexcept;
  Column(const Column&);
  Column& operator=(const Column&);
  Json cfg;
  std::vector<std::string> comp_names, material_names, loading_group_names, history_group_names;
  std::size_t n_comp = 0, Nz = 0;
  double R = 8.31446261815324, L = 0, area = 0, dint = 0, dout = 0, dz = 0, Cpg = 1000;
  int thermal_mode = 1;
  bool gas_solid_thermal_equilibrium = false;
  bool liquid = false;
  double liquid_density = 1000., liquid_heat_capacity = 4180.;
  Vec M, z, widths, epsilon_node, storage_epsilon_node, epsilon_face;
  StateLayout layout;
  std::size_t state_size() const { return layout.size(); }
  void set_step(const Json& step, const Json& context = Json::object());
  Vec initial_state(const Json& definition = Json());
  Derived derived(const Vec& state) const;
  Json endpoint(const std::string& side, const Vec& state) const;
  Json source_states(double t = 0) const;
  Json port_state(const std::string& side, double t = 0) const;
  std::pair<double, double> liquid_resistance() const;
  Vec face_velocities(const Vec& state, double t = 0, const Json& overrides = Json::object()) const;
  Vec rhs(double t, const Vec& state, const Json& overrides = Json::object()) const;
  Vec explicit_transport_rhs(double t, const Vec& state, const Json& overrides = Json::object()) const;
  bool has_reactions() const;
  Vec reaction_component_rates() const;
  double reaction_energy_rate() const;
  Json boundary_component_rates(double t, const Vec& state, const Json& overrides = Json::object()) const;
  Json boundary_enthalpy_rates(double t, const Vec& state, const Json& overrides = Json::object()) const;
  Vec gas_component_heat_capacity(double temperature) const;
  Vec gas_component_enthalpy(double temperature) const;
  Vec gas_component_internal_energy(double temperature) const;
  Json thermal_control(double t = 0) const;
  Json jacket_record(double t, const Vec& state) const;
  Json thermal_utility_rates(double t, const Vec& state, const Json& overrides = Json::object(),
                             const Vec* state_derivative = nullptr) const;
  Json inventories(const Vec& state) const;
  std::pair<Vec, Vec> loading_profiles(const Vec&) const;
  Vec equilibrium_loading(const Vec&) const;
  double stored_energy(const Vec& state) const;
  double energy_inventory(const Vec& state) const { return stored_energy(state); }
  Json state_record(const Vec& state) const;
  Json state_records(const Vec& state) const { return state_record(state); }
  Vec state_atol(double default_value = 1e-8) const;
  double metric(const std::string& name, const Vec& state, double t = 0, const Json& options = Json::object()) const;
  double observable(const std::string& name, const Vec& state, double t = 0, const Json& options = Json::object()) const
  {
    return metric(name, state, t, options);
  }

 private:
  struct Impl;
  std::unique_ptr<Impl> p;
  std::pair<Vec, Vec> liquid_hydraulics(double, const Json&) const;
  std::pair<Vec, Vec> reservoir_concentrations(double t, const Vec& state, const Json& overrides) const;
  std::pair<Vec, Json> external_wall_heat_density(const Vec& Tw, double t, const Vec& Tj = {}) const;
};
using PSAFixedBedModel = Column;
double resolve_time_profile(const Json& value, double t);
}  // namespace cycle

#pragma once
#include <map>
#include <optional>

#include "cycle/types.hpp"
namespace cycle
{
Json configured_sources(const Json& cfg, std::size_t n_comp, double t = 0.0, const Json& context = Json::object());
double valve_molar_flow(double pa, double pb, double temperature, double molar_mass,
                        const Json& options = Json::object(), double gas_constant = 8.31446261815324);
std::pair<double, double> compressor_outlet(double pin, double pout, double temperature, double cp,
                                            double efficiency = .75, double gas_constant = 8.31446261815324);
struct BufferDefinition
{
  std::string name, thermal_mode;
  double volume, initial_pressure, initial_temperature, ua, ambient_temperature;
  Vec initial_composition, initial_concentration, cp_molar;
  double maximum_volume = 1e100;
  std::vector<std::string> aliases;
  Json pressure_control;
};
class DynamicBufferNetwork
{
 public:
  Json cfg;
  std::vector<std::string> comp_names;
  std::map<std::string, std::size_t> comp_index, name_to_index;
  std::map<std::string, std::string> alias_to_name;
  std::vector<BufferDefinition> definitions;
  Vec M;
  double R;
  std::size_t n_comp, block_size, size;
  bool enabled, liquid = false;
  double liquid_density = 1000., liquid_heat_capacity = 4180.;
  DynamicBufferNetwork(const Json&, const std::vector<std::string>&, const Vec&, double);
  Vec initial_state() const;
  Json state_record(const std::string&, const Vec&) const;
  Json records(const Vec&) const;
  Json source_states(const Vec&, double = 0, const Json& context = Json::object()) const;
  std::optional<std::string> buffer_name_for_port(const Json&) const;
  std::vector<std::string> controller_tags() const;
  std::pair<Vec, Json> conservative_derivatives(const Vec&, const Json& component_rates, const Json& enthalpy_rates,
                                                const Json& sources, const Json& volume_rates = Json::object()) const;
  double metric(const std::string&, const Vec&, const std::string& buffer, const Json& component = nullptr) const;
};
}  // namespace cycle

#pragma once
#include <map>
#include <memory>

#include "cycle/column.hpp"
#include "cycle/process.hpp"
namespace cycle
{
class ColumnNetwork
{
 public:
  Json cfg, phase = Json::object(), context = Json::object(), actions = Json::object(), connections = Json::array();
  std::vector<Column> models;
  std::vector<std::string> column_names, comp_names, stream_tags;
  std::map<std::string, std::size_t> column_indices, column_starts, comp_index;
  std::unique_ptr<DynamicBufferNetwork> buffers;
  std::size_t buffer_start = 0, size = 0, n_comp = 0;
  Vec M, cp_molar;
  double R = 8.31446261815324;
  static const std::vector<std::string> utility_names;
  explicit ColumnNetwork(const Json&);
  Column& model(const std::string&);
  const Column& model(const std::string&) const;
  Vec column_state(const std::string&, const Vec&) const;
  Vec initial_state();
  void set_phase(const Json&, const Json& context = Json::object());
  Json endpoint_pressures(const Vec&) const;
  std::pair<Vec, Json> evaluate(double, const Vec&);
  Vec rhs(double t, const Vec& s) { return evaluate(t, s).first; }
  Vec explicit_transport_rhs(double t, const Vec& state);
  Json inventories(const Vec&) const;
  Json state_records(const Vec&) const;
  Vec state_atol(double default_value = 1e-8) const;
  double metric(const std::string&, const Vec&, double t = 0, const Json& options = Json::object());

 private:
  bool legacy_two = false;
  Vec transport_cache_state_;
  Json transport_boundary_overrides_;
  double transport_cache_time_ = std::numeric_limits<double>::quiet_NaN();
  Json solve_liquid_hydraulics(double, const Vec&, const Json&) const;
  Json endpoint(const Json&) const;
  Json refresh_context(double, const Vec&);
  Json record(const Json&, const Vec&, const std::vector<Derived>&, const Json&) const;
};
}  // namespace cycle

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include "cycle/diagnostics.hpp"
#include "cycle/npz.hpp"
#include "cycle/main_input.hpp"
#include "cycle/plotting.hpp"
#include "cycle/runner.hpp"

namespace cycle
{
namespace
{
// ---------------------------------------------------------------------------
// Optional Matplotlib plotting backend.
//
// The bundled tools/plot_cycle_results.py only reads the saved NPZ and CSV results.
// Python remains absent from the numerical path: if it is missing, or the
// script fails for any reason, the caller falls back to the built-in SVG
// overview and says so.
// ---------------------------------------------------------------------------
std::string shell_quote(const std::string& text)
{
  std::string quoted = "'";
  for (char c : text)
  {
    if (c == '\'')
      quoted += "'\\''";
    else
      quoted += c;
  }
  return quoted + "'";
}
std::string environment(const char* key)
{
  const char* value = std::getenv(key);
  return value ? std::string(value) : std::string();
}
// Walk upwards from a starting directory looking for tools/plot_cycle_results.py.
std::filesystem::path search_upwards(std::filesystem::path directory)
{
  std::error_code code;
  directory = std::filesystem::absolute(directory, code);
  if (code) return {};
  for (int depth = 0; depth < 32 && !directory.empty(); ++depth)
  {
    auto candidate = directory / "tools" / "plot_cycle_results.py";
    if (std::filesystem::is_regular_file(candidate)) return candidate;
    auto parent = directory.parent_path();
    if (parent == directory) break;
    directory = parent;
  }
  return {};
}
std::filesystem::path find_plot_script(const std::filesystem::path& output_dir)
{
  auto explicit_path = environment("RUPTURA_PLOT_SCRIPT");
  if (!explicit_path.empty())
  {
    if (!std::filesystem::is_regular_file(explicit_path))
      throw std::runtime_error("RUPTURA_PLOT_SCRIPT does not name a file: " + explicit_path);
    return explicit_path;
  }
  auto root = environment("RUPTURA_CYCLE_CPP_ROOT");
  if (!root.empty())
  {
    auto candidate = std::filesystem::path(root) / "tools" / "plot_cycle_results.py";
    if (std::filesystem::is_regular_file(candidate)) return candidate;
  }
  std::error_code code;
  auto executable = std::filesystem::read_symlink("/proc/self/exe", code);
  if (!code)
  {
    auto installed = executable.parent_path().parent_path() / "share" / "ruptura" / "tools" / "plot_cycle_results.py";
    if (std::filesystem::is_regular_file(installed)) return installed;
    auto found = search_upwards(executable.parent_path());
    if (!found.empty()) return found;
  }
  auto found = search_upwards(output_dir);
  if (!found.empty()) return found;
  return search_upwards(std::filesystem::current_path(code));
}
std::string python_executable()
{
  auto chosen = environment("RUPTURA_PYTHON");
  return chosen.empty() ? std::string("python3") : chosen;
}
double sum(const Vec& x) { return std::accumulate(x.begin(), x.end(), 0.); }
double number(const Json& value, const char* key, double fallback = 0.)
{
  return value.contains(key) ? value.at(key).get<double>() : fallback;
}
Vec slice(const Vec& v, std::size_t offset, std::size_t size)
{
  if (offset + size > v.size()) throw std::invalid_argument("State slice exceeds array");
  return Vec(v.begin() + offset, v.begin() + offset + size);
}
Json clean(Json value)
{
  if (value.is_number_float() && !std::isfinite(value.get<double>())) return nullptr;
  if (value.is_object() || value.is_array())
    for (auto& x : value) x = clean(x);
  return value;
}
void atomic_text(const std::filesystem::path& path, const std::string& text)
{
  static std::atomic<unsigned long> sequence{0};
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  auto temp = path;
  temp += ".tmp-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + "-" +
          std::to_string(sequence++);
  try
  {
    {
      std::ofstream out(temp, std::ios::binary);
      if (!out) throw std::runtime_error("Cannot write " + temp.string());
      out << text;
      out.close();
      if (!out) throw std::runtime_error("Failed writing " + temp.string());
    }
    std::filesystem::rename(temp, path);
  }
  catch (...)
  {
    std::error_code ignored;
    std::filesystem::remove(temp, ignored);
    throw;
  }
}
std::string csv_field(const Json& value)
{
  if (value.is_null()) return {};
  std::string raw;
  if (value.is_object() || value.is_array())
  {
    nlohmann::json sorted = clean(value);
    raw = sorted.dump();
  }
  else if (value.is_string())
    raw = value.get<std::string>();
  else if (value.is_boolean())
    raw = value.get<bool>() ? "True" : "False";
  else if (value.is_number_float() && !std::isfinite(value.get<double>()))
    raw = std::isnan(value.get<double>()) ? "nan" : value.get<double>() > 0 ? "inf" : "-inf";
  else
    raw = value.dump();
  if (raw.find_first_of(",\"\r\n") == std::string::npos) return raw;
  std::string escaped = "\"";
  for (char c : raw)
  {
    escaped += c;
    if (c == '"') escaped += '"';
  }
  return escaped + '"';
}
Json matrix(const Vec& values, std::size_t rows, std::size_t cols)
{
  if (values.size() != rows * cols) throw std::invalid_argument("Output matrix dimensions differ");
  Json result = Json::array();
  for (std::size_t row = 0; row < rows; ++row) result.push_back(slice(values, row * cols, cols));
  return result;
}
std::string slug(const std::string& name) { return std::regex_replace(name, std::regex("[^A-Za-z0-9_.-]"), "_"); }
std::filesystem::path phase_path(const std::filesystem::path& base, int cycle, int phase, const std::string& name)
{
  std::ostringstream out;
  out << "cycle_" << std::setw(4) << std::setfill('0') << cycle << "_phase_" << std::setw(2) << phase << '_'
      << slug(name) << ".npz";
  return base / "phases" / out.str();
}
void add_phase_time(Json& arrays, double start)
{
  Vec relative = arrays.at("time").get<Vec>();
  for (double& t : relative) t -= start;
  arrays["phase_time"] = relative;
}
}  // namespace

Json CycleRunner::stream_record(const Vec& delta) const
{
  Json record = Json::object();
  for (const auto& tag : tags)
  {
    auto values = ledger_value(delta, {"tag", tag});
    if (std::any_of(values.begin(), values.end(), [](double x) { return std::abs(x) > 1e-30; })) record[tag] = values;
  }
  return record;
}

Json CycleRunner::cycle_metrics(const Vec& start, const Vec& end, const Vec& delta, double duration)
{
  const auto feed_tag = cfg.value("FEED_TAG", std::string("feed")),
             light_tag = cfg.value("LIGHT_PRODUCT_TAG", std::string("light_product")),
             heavy_tag = cfg.value("HEAVY_PRODUCT_TAG", std::string("heavy_product"));
  Vec feed(network.n_comp, 0.);
  Json feed_by_stage = Json::object();
  const auto feed_tags = cfg.value("FRESH_FEED_TAGS", Json::array({feed_tag}));
  for (const auto& tag_value : feed_tags)
  {
    const auto tag = tag_value.get<std::string>();
    auto amount = ledger_slices.count({"external_tag_in", tag}) ? ledger_value(delta, {"external_tag_in", tag})
                                                                : Vec(network.n_comp, 0.);
    feed_by_stage[tag] = sum(amount);
    for (size_t i = 0; i < feed.size(); ++i) feed[i] += amount[i];
  }
  const auto gross = ledger_value(delta, {"external_tag_out", light_tag}),
             gross_heavy = ledger_value(delta, {"external_tag_out", heavy_tag});
  const auto recycled_light = ledger_value(delta, {"recycle", light_tag}),
             recycled_heavy = ledger_value(delta, {"recycle", heavy_tag});
  const auto ncomp = network.n_comp;
  Vec net(ncomp), heavy(ncomp), fresh(ncomp, 0.);
  const double inventory_tolerance = 10. * number(cfg, "LEDGER_ATOL", 1e-10);
  bool light_available = true, heavy_available = true;
  double light_deficit = 0., heavy_deficit = 0.;
  for (std::size_t i = 0; i < ncomp; ++i)
  {
    net[i] = gross[i] - recycled_light[i];
    heavy[i] = gross_heavy[i] - recycled_heavy[i];
    light_available = light_available && net[i] >= -inventory_tolerance;
    heavy_available = heavy_available && heavy[i] >= -inventory_tolerance;
    light_deficit += std::max(-net[i], 0.);
    heavy_deficit += std::max(-heavy[i], 0.);
  }
  for (const auto& tag : tags)
  {
    auto inlet = ledger_value(delta, {"external_tag_in", tag}), recycle = ledger_value(delta, {"recycle", tag});
    for (std::size_t i = 0; i < ncomp; ++i) fresh[i] += inlet[i] - recycle[i];
  }
  const auto inventories0 = network.inventories(start), inventories1 = network.inventories(end);
  const auto inv0 = inventories0.at("total").get<Vec>(), inv1 = inventories1.at("total").get<Vec>();
  const auto incoming = ledger_value(delta, {"external_in"}), outgoing = ledger_value(delta, {"external_out"});
  const auto reaction=ledger_slices.count({"reaction"})?ledger_value(delta,{"reaction"}):Vec(ncomp,0.);
  Vec balance(ncomp);
  double mass_error = 0., throughput_error = 0.;
  for (std::size_t i = 0; i < ncomp; ++i)
  {
    balance[i] = inv0[i] + incoming[i] - outgoing[i] + reaction[i] - inv1[i];
    mass_error =
        std::max(mass_error, std::abs(balance[i]) / std::max(std::abs(inv0[i]) + incoming[i] + outgoing[i], 1e-12));
    throughput_error = std::max(throughput_error, std::abs(balance[i]) / std::max(incoming[i] + outgoing[i], 1e-12));
  }
  const auto light_idx = network.comp_index.at(cfg.value("TARGET_LIGHT_COMPONENT", network.comp_names.back())),
             heavy_idx = network.comp_index.at(cfg.value("TARGET_HEAVY_COMPONENT", network.comp_names.front()));
  auto fraction = [](double numerator, double denominator) -> Json
  { return denominator > 1e-30 ? Json(numerator / denominator) : Json(); };
  Json utilities = Json::object();
  for (const auto& name : utility_names)
  {
    std::string key = name;
    if (key.size() >= 2 && key.substr(key.size() - 2) == "_W") key.resize(key.size() - 2);
    utilities[key + "_J"] = ledger_value(delta, {"utility", name}).at(0);
  }
  auto utility = [&](const char* key) { return number(utilities, key); };
  double mass = 0., volume = 0.;
  Json geometry = Json::object();
  for (const auto& name : network.column_names)
  {
    const auto& model = network.model(name);
    const auto& inventory = inventories0.at("columns").at(name);
    const double column_mass = sum(inventory.at("material_mass_kg").get<Vec>());
    mass += column_mass;
    volume += model.area * model.L;
    geometry[name] = {{"length_m", model.L},
                      {"area_m2", model.area},
                      {"volume_m3", model.area * model.L},
                      {"adsorbent_mass_kg", column_mass}};
  }
  auto tagged_consumption = [&](const char* attribute, std::vector<std::string> aliases)
  {
    if (cfg.contains(attribute)) aliases = {cfg.at(attribute).get<std::string>()};
    double total = 0.;
    for (const auto& tag : aliases)
      if (ledger_slices.count({"tag", tag})) total += sum(ledger_value(delta, {"tag", tag}));
    return total;
  };
  Json result = {{"cycle_duration_s", duration},
                 {"feed_mol", sum(feed)},
                 {"feed_component_mol", feed},
                 {"fresh_inlet_component_mol", fresh},
                 {"recovery_basis", "all_external_fresh_inlets"},
                 {"light_product_available", light_available},
                 {"heavy_product_available", heavy_available},
                 {"process_feasible", light_available && heavy_available},
                 {"gross_light_product_mol", sum(gross)},
                 {"light_product_mol", sum(net)},
                 {"heavy_product_mol", sum(heavy)},
                 {"net_light_product_component_mol", net},
                 {"heavy_product_component_mol", heavy},
                 {"gross_light_product_component_mol", gross},
                 {"light_product_purity", light_available ? fraction(net[light_idx], sum(net)) : Json()},
                 {"gross_light_product_purity", fraction(gross[light_idx], sum(gross))},
                 {"heavy_product_purity", heavy_available ? fraction(heavy[heavy_idx], sum(heavy)) : Json()},
                 {"light_component_recovery", light_available ? fraction(net[light_idx], fresh[light_idx]) : Json()},
                 {"gross_light_component_recovery", fraction(gross[light_idx], fresh[light_idx])},
                 {"heavy_component_recovery", heavy_available ? fraction(heavy[heavy_idx], fresh[heavy_idx]) : Json()},
                 {"product_buffer_consumption_mol", sum(recycled_light)},
                 {"product_buffer_deficit_mol", light_deficit},
                 {"heavy_product_deficit_mol", heavy_deficit},
                 {"per_column_geometry", geometry},
                 {"purge_consumption_mol", tagged_consumption("PURGE_TAG", {"purge", "regeneration_purge"})},
                 {"repressurization_consumption_mol",
                  tagged_consumption("REPRESSURIZATION_GAS_TAG", {"repressurization_gas", "pressurization_gas"})},
                 {"mass_balance_error", mass_error},
                 {"mass_balance_throughput_error", throughput_error},
                 {"component_balance_residual_mol", balance},
                 {"reaction_component_mol",reaction},
                 {"external_in_component_mol", incoming},
                 {"external_out_component_mol", outgoing},
                 {"total_adsorbent_mass_kg", mass},
                 {"total_bed_volume_m3", volume},
                 {"number_of_columns", network.models.size()},
                 {"productivity_mol_kg_s", light_available ? fraction(net[light_idx], mass * duration) : Json()},
                 {"thermal_net_heat_J_cycle", utility("net_heat_J")},
                 {"thermal_heat_input_J_cycle", utility("heating_J") + utility("buffer_heating_J")},
                 {"cooling_duty_J_cycle", utility("cooling_J") + utility("buffer_cooling_J")},
                 {"electric_heater_work_J_cycle", utility("electric_heater_J")},
                 {"purge_and_feed_gas_heating_J_cycle", utility("gas_heating_J")},
                 {"purge_and_feed_gas_cooling_J_cycle", utility("gas_cooling_J")},
                 {"ambient_heat_loss_J_cycle", utility("heat_loss_J")},
                 {"utilities_J", utilities}};
  if (network.models.front().liquid)
  {
    double vin = ledger_value(delta, {"volume_external_in"})[0], vout = ledger_value(delta, {"volume_external_out"})[0];
    double feed_volume = ledger_value(delta, {"volume_tag_in", feed_tag})[0];
    double product_volume =
        ledger_value(delta, {"volume_tag_out", light_tag})[0] - ledger_value(delta, {"volume_tag_in", light_tag})[0];
    Vec concentrations = net;
    for (double& value : concentrations) value = product_volume > 0. ? value / product_volume : 0.;
    result.update({{"external_liquid_in_m3", vin},
                   {"external_liquid_out_m3", vout},
                   {"fresh_feed_volume_m3", feed_volume},
                   {"light_product_volume_m3", product_volume},
                   {"light_product_concentration_mol_m3", concentrations},
                   {"purity_basis", "tracked solutes; implicit solvent excluded"},
                   {"pump_work_J_cycle", utility("pump_J")},
                   {"thermal_energy_basis",
                    "sensible and adsorption energy; liquid pressure work and viscous heating excluded"}});
  }
  if (cfg.contains("FRESH_FEED_TAGS"))
    result.update({{"fresh_feed_mol", sum(feed)},
                   {"fresh_feed_by_stage_mol", feed_by_stage},
                   {"adsorption_feed_mol", sum(ledger_value(delta, {"external_tag_in", feed_tag}))},
                   {"average_fresh_feed_mol_s", duration > 0 ? sum(feed) / duration : 0.}});
  if (inventories0.contains("energy_J") && inventories1.contains("energy_J"))
  {
    double energy0 = number(inventories0, "energy_J"), energy1 = number(inventories1, "energy_J"),
           enthalpy_in = ledger_value(delta, {"energy", "external_enthalpy_in_W"}).at(0),
           enthalpy_out = ledger_value(delta, {"energy", "external_enthalpy_out_W"}).at(0);
    const double reaction_energy=ledger_slices.count({"energy","reaction_energy_W"})?ledger_value(delta,{"energy","reaction_energy_W"})[0]:0.;
    double external_heat = utility("net_heat_J") + utility("buffer_heating_J") - utility("buffer_cooling_J"),
           shaft_work = utility("compressor_J") + utility("vacuum_J");
    double residual = energy0 + reaction_energy + enthalpy_in - enthalpy_out + external_heat + shaft_work - energy1,
           throughput = std::abs(enthalpy_in) + std::abs(enthalpy_out) + std::abs(external_heat) + std::abs(shaft_work),
           scale = std::max({std::abs(energy0), std::abs(energy1), throughput, 1e-12});
    result.update(Json{{"stored_energy_start_J", energy0},
                       {"stored_energy_end_J", energy1},
                       {"external_enthalpy_in_J", enthalpy_in},
                       {"external_enthalpy_out_J", enthalpy_out},
                       {"reaction_energy_J",reaction_energy},
                       {"energy_balance_residual_J", residual},
                       {"energy_balance_error", std::abs(residual) / scale},
                       {"energy_balance_throughput_error", std::abs(residual) / std::max(throughput, 1e-12)}});
  }
  if (utility("compressor_J") > 0) result["network_compression_work_J_cycle"] = utility("compressor_J");
  if (utility("vacuum_J") > 0) result["vacuum_work_J_cycle"] = utility("vacuum_J");
  return result;
}

std::pair<Json, bool> CycleRunner::css(const Vec& start, const Vec& end, const Json& metrics,
                                       const Json& previous_metrics, const Json& compositions,
                                       const Json& previous_compositions)
{
  const double rtol = number(cfg, "CSS_STATE_RTOL", 1e-5);
  const auto absolute = network.state_atol(number(cfg, "CSS_STATE_ATOL", 1e-8));
  if (start.size() != physical_size || end.size() != physical_size)
    throw std::invalid_argument("CSS state dimensions differ");
  Vec scaled(physical_size);
  double state_error = 0.;
  for (std::size_t i = 0; i < physical_size; ++i)
  {
    scaled[i] = std::abs(end[i] - start[i]) /
                std::max(absolute[i] + rtol * std::max(std::abs(start[i]), std::abs(end[i])), 1e-300);
    state_error = std::max(state_error, scaled[i]);
  }
  const auto inv0 = network.inventories(start).at("total").get<Vec>(),
             inv1 = network.inventories(end).at("total").get<Vec>();
  double inventory_error = 0.;
  for (std::size_t i = 0; i < inv0.size(); ++i)
    inventory_error = std::max(inventory_error,
                               std::abs(inv1[i] - inv0[i]) / std::max({std::abs(inv0[i]), std::abs(inv1[i]), 1e-12}));
  const double infinity = std::numeric_limits<double>::infinity();
  double metric_change = previous_metrics.is_null() ? infinity : 0.;
  for (const char* name :
       {"light_product_purity", "light_component_recovery", "heavy_product_purity", "heavy_component_recovery"})
  {
    Json old = previous_metrics.is_null() ? Json() : previous_metrics.value(name, Json()),
         current = metrics.value(name, Json());
    if (old.is_null() && current.is_null()) continue;
    metric_change =
        std::max(metric_change,
                 !old.is_null() && !current.is_null() ? std::abs(current.get<double>() - old.get<double>()) : infinity);
  }
  std::set<std::string> used_sources;
  std::function<void(const Json&)> collect = [&](const Json& v)
  {
    if (v.is_object())
    {
      for (auto it = v.begin(); it != v.end(); ++it)
        if (it.key() == "composition_source" && it.value().is_string())
          used_sources.insert(it.value().get<std::string>());
        else
          collect(it.value());
    }
    else if (v.is_array())
      for (const auto& x : v) collect(x);
  };
  collect(phases);
  double source_error = 0.;
  for (const auto& name : used_sources)
  {
    bool now = compositions.contains(name), before = previous_compositions.contains(name);
    if (!now && !before) continue;
    if (!now || !before)
    {
      source_error = infinity;
      continue;
    }
    const auto a = compositions.at(name).get<Vec>(), b = previous_compositions.at(name).get<Vec>();
    if (a.size() != b.size()) throw std::invalid_argument("CSS source-composition dimensions differ");
    for (std::size_t i = 0; i < a.size(); ++i) source_error = std::max(source_error, std::abs(a[i] - b[i]));
  }
  Json column_errors = Json::object(), field_errors = Json::object();
  for (const auto& name : network.column_names)
  {
    const auto offset = network.column_starts.at(name), size = network.model(name).state_size();
    column_errors[name] = *std::max_element(scaled.begin() + offset, scaled.begin() + offset + size);
    const auto& model = network.model(name);
    const auto nz = model.Nz, gas = model.n_comp * nz;
    auto maximum = [&](std::size_t begin, std::size_t end)
    { return begin < end ? *std::max_element(scaled.begin() + offset + begin, scaled.begin() + offset + end) : 0.; };
    field_errors[name] = {{"gas_concentration", maximum(0, gas)},
                          {"loadings", maximum(gas, gas * (1 + model.layout.n_loading_groups))},
                          {"kinetic_history", maximum(gas * (1 + model.layout.n_loading_groups),
                                                      size - (3 + model.layout.transient_htf) * nz)},
                          {"gas_temperature", maximum(size - 3 * nz, size - 2 * nz)},
                          {"solid_temperature", maximum(size - 2 * nz, size - nz)},
                          {"wall_temperature", maximum(size - nz, size)}};
    if (model.layout.transient_htf) field_errors[name]["jacket_temperature"] = maximum(size - 4 * nz, size - 3 * nz);
  }
  Json result = {{"css_state_error", state_error},
                 {"css_inventory_error", inventory_error},
                 {"css_metric_change", std::isfinite(metric_change) ? Json(metric_change) : Json()},
                 {"css_source_composition_change", std::isfinite(source_error) ? Json(source_error) : Json()},
                 {"css_column_errors", column_errors},
                 {"css_field_errors", field_errors}};
  const double tolerance = number(cfg, "CSS_METRIC_RTOL", 1e-4);
  // Unlike purity fractions, duration, productivity and duties need relative
  // checks. This matters for criterion-terminated phases with varying times.
  double performance_change = previous_metrics.is_null() ? infinity : 0.;
  double performance_error = performance_change;
  for (const char* key : {"cycle_duration_s", "productivity_mol_kg_s", "thermal_heat_input_J_cycle",
                          "cooling_duty_J_cycle", "network_compression_work_J_cycle", "vacuum_work_J_cycle"})
  {
    const std::string name = key;
    const bool duty = name.ends_with("_J_cycle");
    Json a = metrics.value(key, duty ? Json(0.) : Json()),
         b = previous_metrics.is_null() ? Json() : previous_metrics.value(key, duty ? Json(0.) : Json());
    if (a.is_null() && b.is_null()) continue;
    double change = infinity;
    double error = infinity;
    if (a.is_number() && b.is_number())
    {
      const double scale = std::max(std::abs(a.get<double>()), std::abs(b.get<double>()));
      const double diff = std::abs(a.get<double>() - b.get<double>());
      const double absolute_tolerance = duty                         ? number(cfg, "CSS_DUTY_ATOL_J", 1e-6)
                                        : name == "cycle_duration_s" ? number(cfg, "CSS_DURATION_ATOL_S", 1e-6)
                                                                     : number(cfg, "CSS_PRODUCTIVITY_ATOL", 1e-12);
      change = diff / std::max(scale, 1e-12);
      error = diff / std::max(absolute_tolerance + tolerance * scale, 1e-300);
    }
    performance_change = std::max(performance_change, change);
    performance_error = std::max(performance_error, error);
  }
  result["css_performance_relative_change"] = std::isfinite(performance_change) ? Json(performance_change) : Json();
  result["css_performance_error"] = std::isfinite(performance_error) ? Json(performance_error) : Json();
  const double balance = number(metrics, "mass_balance_error", 0.);
  bool satisfied = state_error <= 1. && inventory_error <= number(cfg, "CSS_INVENTORY_RTOL", 1e-5) &&
                   metric_change <= tolerance && source_error <= tolerance && performance_error <= 1. &&
                   std::isfinite(balance) && balance <= number(cfg, "CSS_MASS_BALANCE_TOLERANCE", 1e-3);
  return {result, satisfied};
}

void CycleRunner::update_convergence_report()
{
  auto report = assess_convergence(progress.at("cycle_rows"), cfg);
  auto& summary = progress["summary"];
  summary["convergence_report"] = report;
  summary["plot_cycle"] = report["representative_cycle"];
  summary["convergence_status"] = report["status"];
  summary["reporting_metrics"] = summary["final_metrics"];
  summary["reporting_metrics_basis"] = "last completed cycle";
  const auto& records = progress.value("cycle_records", Json::array());
  if (!records.empty() && report.contains("window_start_cycle"))
  {
    std::size_t begin = 0;
    while (begin < records.size() && records[begin]["cycle"].get<int>() < report["window_start_cycle"].get<int>())
      ++begin;
    // Periodic candidates: aggregate an integer number of complete periods.
    if (report["detected_period_cycles"].is_number())
    {
      const auto period = report["detected_period_cycles"].get<std::size_t>();
      begin += (records.size() - begin) % period;
    }
    if (begin < records.size())
    {
      Vec total(ledger_size, 0.);
      double duration = 0.;
      for (std::size_t i = begin; i < records.size(); ++i)
      {
        auto ledger = records[i]["ledger"].get<Vec>();
        for (std::size_t k = 0; k < ledger_size; ++k) total[k] += ledger[k];
        duration += records[i]["metrics"]["cycle_duration_s"].get<double>();
      }
      auto metrics =
          cycle_metrics(records[begin]["start"].get<Vec>(), records.back()["end"].get<Vec>(), total, duration);
      summary["window_aggregate_metrics"] = metrics;
      summary["window_aggregate_cycles"] = records.size() - begin;
      summary["window_aggregate_start_cycle"] = records[begin]["cycle"];
      summary["window_aggregate_end_cycle"] = records.back()["cycle"];
      if (!summary.value("converged", false) && run_mode == "css")
      {
        summary["reporting_metrics"] = metrics;
        summary["reporting_metrics_basis"] =
            "tail-window summed amounts, duties and duration; ratios of sums (diagnostic only)";
      }
    }
  }
  if (!summary.value("converged", false) && run_mode == "css" &&
      cfg.value("SIMULATION_TYPE", std::string()) != "Breakthrough")
    std::cout << "WARNING: " << report["warning"].get<std::string>()
              << " Plot selection: " << report["selection_reason"] << "; cycle " << summary["plot_cycle"] << '\n';
  // Stable key order also when restoring legacy checkpoints lacking metadata.
  nlohmann::json canonical_summary = summary;
  summary = canonical_summary;
}

void CycleRunner::write_csv(const std::string& filename, const Json& rows) const
{
  std::set<std::string> all;
  for (const auto& row : rows)
    for (auto it = row.begin(); it != row.end(); ++it) all.insert(it.key());
  std::vector<std::string> fields;
  for (const char* key : {"cycle", "phase_index", "phase_name", "stream", "component"})
    if (all.erase(key)) fields.emplace_back(key);
  fields.insert(fields.end(), all.begin(), all.end());
  std::ostringstream text;
  for (std::size_t i = 0; i < fields.size(); ++i)
  {
    if (i) text << ',';
    text << csv_field(Json(fields[i]));
  }
  text << "\r\n";
  for (const auto& row : rows)
  {
    for (std::size_t i = 0; i < fields.size(); ++i)
    {
      if (i) text << ',';
      text << csv_field(row.value(fields[i], Json()));
    }
    text << "\r\n";
  }
  atomic_text(output_dir / filename, text.str());
}

Json CycleRunner::output_arrays(int requested_cycle, int requested_phase) const
{
  Json times = Json::array(), cycles = Json::array(), phase_indices = Json::array(), states = Json::array();
  std::vector<Vec> selected;
  std::vector<Json> selected_pressures;
  Vec local_times;
  std::map<std::pair<int, int>, double> phase_starts;
  // Scalar phase records preserve true local-time origins even if an older
  // profile archive was removed before restarting in the middle of a ramp.
  for (const auto& row : progress.value("phase_rows", Json::array()))
    phase_starts[{row["cycle"].get<int>(), row["phase_index"].get<int>() - 1}] = row["start_time_s"].get<double>();
  if (progress.contains("cycle") && progress.contains("phase") && progress.contains("phase_start_time"))
    phase_starts.try_emplace(std::make_pair(progress["cycle"].get<int>(), progress["phase"].get<int>()),
                             progress["phase_start_time"].get<double>());
  const auto& history = progress.at("history_state");
  for (std::size_t i = 0; i < history.size(); ++i)
  {
    int cycle = progress.at("history_cycle").at(i).get<int>(), phase = progress.at("history_phase").at(i).get<int>();
    const double time = progress.at("history_t").at(i);
    phase_starts.try_emplace({cycle, phase}, time);
    if ((requested_cycle >= 0 && cycle != requested_cycle) || (requested_phase >= 0 && phase != requested_phase))
      continue;
    times.push_back(progress.at("history_t").at(i));
    cycles.push_back(cycle);
    phase_indices.push_back(phase);
    selected.push_back(history[i].get<Vec>());
    selected_pressures.push_back(progress.contains("history_pressure") ? progress["history_pressure"][i] : Json());
    local_times.push_back(time - phase_starts.at({cycle, phase}));
    states.push_back(selected.back());
  }
  Json transpose = Json::array();
  for (std::size_t i = 0; !compact_output && i < physical_size; ++i)
  {
    Vec row;
    for (const auto& state : selected) row.push_back(state.at(i));
    transpose.push_back(row);
  }
  Json result = {{"time", times},
                 {"t", times},
                 {"state", states},
                 {"y", transpose},
                 {"cycle", cycles},
                 {"phase", phase_indices},
                 {"column_names", network.column_names},
                 {"comp_names", network.comp_names}};
  if (compact_output)
  {
    result.erase("t");
    result.erase("y");
  }
  // Optional literature diagnostics: reconstruct rates with the same native
  // network equations, not endpoint composition or finite differences of inventory.
  // These reference cases use prescribed sources (no history-dependent recycle).
  if (cfg.value("EXPORT_VALIDATION_FLUXES", false))
  {
    std::function<void(const Json&)> check_sources = [&](const Json& value)
    {
      if (value.is_object())
      {
        if (value.contains("composition_source"))
          throw std::invalid_argument(
              "Validation flux export does not support history-dependent composition_source controls");
        for (const auto& child : value) check_sources(child);
      }
      else if (value.is_array())
      {
        for (const auto& child : value) check_sources(child);
      }
    };
    check_sources(phases);
    ColumnNetwork recorder(cfg);
    if (recorder.buffers->size != 0)
      throw std::invalid_argument(
          "Validation flux export currently requires prescribed sources and direct column links");
    for (const auto& name : network.column_names)
      for (const auto side : {"left", "right"})
        for (const auto direction : {"in", "out"})
          result[name + "_" + side + "_" + direction + "_mol_s"] = Json::array();
    std::set<std::string> links;
    for (const auto& phase : phases)
      for (const auto& link : phase.value("connections", Json::array()))
        if (link.contains("name")) links.insert(link["name"].get<std::string>());
    for (const auto& link : links) result["connection_" + link + "_mol_s"] = Json::array();
    int last_cycle = -1, last_phase = -1;
    for (std::size_t sample = 0; sample < selected.size(); ++sample)
    {
      const int ci = cycles[sample], pi = phase_indices[sample];
      const auto physical = slice(selected[sample], 0, physical_size);
      if (ci != last_cycle || pi != last_phase)
      {
        Json ctx = {{"global_time_offset", times[sample].get<double>() - local_times[sample]},
                    {"phase_start_pressures", recorder.endpoint_pressures(physical)}};
        recorder.set_phase(phases[pi], ctx);
        last_cycle = ci;
        last_phase = pi;
      }
      const auto diag = recorder.evaluate(local_times[sample], physical).second;
      for (const auto& name : network.column_names)
        for (const auto side : {"left", "right"})
          for (const auto direction : {"in", "out"})
          {
            const std::string key = std::string(side) + "_" + direction;
            result[name + "_" + key + "_mol_s"].push_back(diag["port_rates"][name][key]);
          }
      for (const auto& link : links)
      {
        double flow = 0.;
        for (const auto& transfer : diag["connections"])
          if (transfer["name"] == link) flow = transfer["signed_molar_flow"];
        result["connection_" + link + "_mol_s"].push_back(flow);
      }
    }
  }
  for (const auto& name : network.column_names)
  {
    const auto& model = network.model(name);
    Column recorder(model);
    bool jacket_output = model.layout.transient_htf;
    for (const auto& phase : phases)
    {
      auto thermal = phase.at("columns").at(name).value("thermal", Json::object());
      auto mode = thermal.is_object() ? thermal.value("mode", std::string()) : std::string();
      jacket_output |= mode == "htf" || mode == "heat_transfer_fluid" || mode == "cooling_fluid" || mode == "jacket" ||
                       mode == "coolant" || mode == "heating_jacket" || mode == "cooling_jacket";
    }
    result[name + "_z"] = model.z;
    for (const char* field : {"c", "q", "Tg", "Ts", "Tw", "P", "state"}) result[name + "_" + field] = Json::array();
    if (jacket_output)
      for (const char* field : {"Tj", "htf_inlet_temperature", "htf_outlet_temperature", "htf_active",
                                "htf_wall_heat_W", "htf_boundary_heat_W", "htf_storage_rate_W"})
        result[name + "_" + field] = Json::array();
    for (std::size_t sample = 0; sample < selected.size(); ++sample)
    {
      const auto& state = selected[sample];
      auto local = network.column_state(name, state);
      auto derived = model.derived(local);
      if (model.liquid)
      {
        if (selected_pressures[sample].is_null() || !selected_pressures[sample].contains(name))
          throw std::runtime_error("Recorded liquid pressure profile is unavailable");
        derived.P = selected_pressures[sample][name].get<Vec>();
      }
      result[name + "_c"].push_back(matrix(derived.c, model.n_comp, model.Nz));
      result[name + "_q"].push_back(matrix(derived.q, model.n_comp, model.Nz));
      result[name + "_Tg"].push_back(derived.Tg);
      result[name + "_Ts"].push_back(derived.Ts);
      result[name + "_Tw"].push_back(derived.Tw);
      if (jacket_output)
      {
        // Only thermal controls are needed for this reconstruction. Closed
        // dummy gas ports avoid resolving historical gas-network sources.
        recorder.set_step(
            {{"left_port", {{"mode", "closed"}}},
             {"right_port", {{"mode", "closed"}}},
             {"thermal", phases[phase_indices[sample].get<int>()]["columns"][name].value("thermal", Json::object())}});
        auto jacket = recorder.jacket_record(local_times[sample], local);
        result[name + "_Tj"].push_back(
            jacket.value("temperature_profile", Json(Vec(model.Nz, jacket["temperature"].get<double>()))));
        result[name + "_htf_inlet_temperature"].push_back(jacket["temperature"]);
        result[name + "_htf_outlet_temperature"].push_back(jacket.value("outlet_temperature", jacket["temperature"]));
        result[name + "_htf_active"].push_back(jacket["mode"] == "htf" ? 1. : 0.);
        result[name + "_htf_wall_heat_W"].push_back(jacket.value("wall_heat_W", 0.));
        result[name + "_htf_boundary_heat_W"].push_back(jacket.value("boundary_heat_W", 0.));
        result[name + "_htf_storage_rate_W"].push_back(jacket.value("storage_rate_W", 0.));
      }
      result[name + "_P"].push_back(derived.P);
      if (!compact_output) result[name + "_state"].push_back(local);
    }
    if (compact_output) result.erase(name + "_state");
  }
  return result;
}

void CycleRunner::write_outputs(bool final)
{
  write_csv("cycle_metrics.csv", progress.at("cycle_rows"));
  if (!compact_output) write_csv("cycle_summary.csv", progress.at("cycle_rows"));
  write_csv("phase_metrics.csv", progress.at("phase_rows"));
  if (!compact_output) write_csv("step_summary.csv", progress.at("phase_rows"));
  write_csv("streams.csv", progress.at("stream_rows"));
  Json balances = Json::array();
  for (const auto& row : progress.at("cycle_rows"))
    for (std::size_t i = 0; i < network.n_comp; ++i)
      balances.push_back({{"cycle", row.at("cycle")},
                          {"component", network.comp_names[i]},
                          {"external_in_mol", row.at("external_in_component_mol").at(i)},
                          {"external_out_mol", row.at("external_out_component_mol").at(i)},
                          {"balance_residual_mol", row.at("component_balance_residual_mol").at(i)},
                          {"reaction_mol",row.value("reaction_component_mol",Json(Vec(network.n_comp,0.))).at(i)},
                          {"fresh_inlet_mol", row.at("fresh_inlet_component_mol").at(i)},
                          {"net_light_product_mol", row.at("net_light_product_component_mol").at(i)},
                          {"net_heavy_product_mol", row.at("heavy_product_component_mol").at(i)}});
  write_csv("component_balances.csv", balances);
  if (!final) return;
  if (phase_archives)
    for (const auto& row : progress.at("phase_rows"))
    {
      int cycle = row.at("cycle").get<int>(), phase = row.at("phase_index").get<int>();
      auto path = phase_path(output_dir, cycle, phase, row.at("phase_name").get<std::string>());
      if (!std::filesystem::exists(path))
      {
        auto arrays = output_arrays(cycle, phase - 1);
        if (arrays["time"].empty()) continue;
        arrays["termination_reason"] = row.at("termination_reason");
        add_phase_time(arrays, row.at("start_time_s").get<double>());
        write_npz(path, arrays, compress_output);
      }
    }
  write_npz(output_dir / "states.npz", output_arrays(), compress_output);
  if (progress.at("summary").contains("convergence_report"))
  {
    const auto& report = progress.at("summary").at("convergence_report");
    atomic_text(output_dir / "convergence_report.json", clean(report).dump(2));
    Json statistics = Json::array();
    for (const char* kind : {"error_statistics", "kpi_statistics"})
    {
      const auto group = report.value(kind, Json::object());
      for (const auto& [key, value] : group.items())
      {
        Json row = value;
        row["quantity"] = key;
        row["kind"] = kind;
        statistics.push_back(row);
      }
    }
    write_csv("window_statistics.csv", statistics);
    if (selected_archive && report["representative_cycle"].is_number())
    {
      auto arrays = output_arrays(report["representative_cycle"].get<int>());
      arrays["selection_reason"] = report["selection_reason"];
      arrays["convergence_status"] = report["status"];
      if (!arrays["time"].empty()) write_npz(output_dir / "selected_cycle.npz", arrays, compress_output);
    }
  }
  if (progress.contains("summary") && !progress.at("summary").is_null())
    atomic_text(output_dir / "summary.json", clean(progress.at("summary")).dump(2));
  Json column_slices = Json::object(), phase_names = Json::array();
  for (const auto& name : network.column_names)
  {
    const auto start = network.column_starts.at(name);
    column_slices[name] = Json::array({start, start + network.model(name).state_size()});
  }
  for (const auto& phase : phases) phase_names.push_back(phase.at("name"));
  auto current = progress.at("state").get<Vec>();
  Json metadata = {
      {"version", "2.0-integrated-fvm"},
      {"solver", solver_method},
      {"solver_version", solver_method == "CVODE" ? cvode_version() : "native BDF/NDF"},
      {"component_names", network.comp_names},
      {"column_names", network.column_names},
      {"column_slices", column_slices},
      {"phases", phase_names},
      {"fingerprint", fingerprint},
      {"input_document", cfg.value("INPUT_DOCUMENT", Json())},
      {"velocity_basis", "interstitial"},
      {"time_units", "s"},
      {"pressure_units", "Pa"},
      {"concentration_units", "mol/m3_fluid"},
      {"loading_units", "mol/kg_adsorbent"},
      {"output_recording",
       {{"profile", compact_output ? "compact" : "full"},
        {"compressed_npz", compress_output},
        {"profile_cycles_retained", compact_output ? history_cycles : 0},
        {"sample_interval_s", compact_output ? output_interval : 0},
        {"samples_per_phase_target", output_samples},
        {"first_profile_cycle", progress["history_cycle"].front()},
        {"last_profile_cycle", progress["history_cycle"].back()},
        {"profile_history_incomplete", progress.value("profile_history_incomplete", false)},
        {"note",
         "Sampled profiles do not determine integrated cycle KPIs or CSS. Phase/segment endpoints are retained."}}},
      // Fixed reference state for the normalised plot panels. It is a plotting
      // convenience only; no equation uses it.
      {"plot_reference",
       {{"pressure_Pa", cfg.value("PLOT_REFERENCE_PRESSURE", cfg.value("P_FEED", cfg.value("P_init", 1e5)))},
        {"temperature_K", cfg.value("PLOT_REFERENCE_TEMPERATURE", cfg.value("T_feed", cfg.value("T_init", 298.15)))},
        {"feed_composition", cfg.value("y_feed", Json::array())},
        {"R", cfg.value("R", 8.31446261815324)}}},
      {"column_records", network.state_records(slice(current, 0, physical_size))}};
  atomic_text(output_dir / "metadata.json", clean(metadata).dump(2));
}

void CycleRunner::commit_phase(const Json& phase, const PhaseSolution& solution)
{
  if (solution.times.empty()) throw std::invalid_argument("Cannot commit phase without final time");
  double duration = solution.times.back();
  const auto streams = stream_record(solution.ledger);
  const int cycle = progress.at("cycle").get<int>(), phase_index = progress.at("phase").get<int>();
  const double phase_start = progress.at("phase_start_time").get<double>();
  for (auto it = streams.begin(); it != streams.end(); ++it)
  {
    auto amount = it.value().get<Vec>();
    double total = sum(amount);
    if (total > 1e-30)
    {
      for (double& value : amount) value /= total;
      progress["stream_compositions"][it.key()] = amount;
    }
  }
  Json row = {{"cycle", cycle},
              {"phase_index", phase_index + 1},
              {"phase_name", phase.at("name")},
              {"step_name", phase.at("name")},
              {"start_time_s", phase_start},
              {"duration_s", duration},
              {"termination_reason", solution.reason},
              {"action_results", solution.action_results},
              {"stream_moles", streams}};
  if (cfg.contains("FRESH_FEED_TAGS"))
  {
    Json columns = Json::object();
    for (const auto& name : network.column_names)
      for (const auto& utility : utility_names)
      {
        std::string key = utility;
        if (key.ends_with("_W")) key.resize(key.size() - 2);
        columns[name][key + "_J"] = ledger_value(solution.ledger, {"column_utility", name, utility})[0];
      }
    row["column_utilities_J"] = columns;
  }
  progress["phase_rows"].push_back(row);
  for (const auto& tag : tags)
  {
    auto incoming = ledger_value(solution.ledger, {"external_tag_in", tag}),
         outgoing = ledger_value(solution.ledger, {"external_tag_out", tag}),
         recycle = ledger_value(solution.ledger, {"recycle", tag});
    auto any = [](const Vec& v) { return std::any_of(v.begin(), v.end(), [](double x) { return x != 0.; }); };
    if (!any(incoming) && !any(outgoing) && !any(recycle)) continue;
    for (std::size_t i = 0; i < network.n_comp; ++i)
      progress["stream_rows"].push_back({{"cycle", cycle},
                                         {"phase_index", phase_index + 1},
                                         {"phase_name", phase.at("name")},
                                         {"stream", tag},
                                         {"component", network.comp_names[i]},
                                         {"external_in_mol", incoming[i]},
                                         {"external_out_mol", outgoing[i]},
                                         {"product_recycle_withdrawal_mol", recycle[i]}});
  }
  if (phase_archives)
  {
    auto arrays = output_arrays(cycle, phase_index);
    arrays["termination_reason"] = solution.reason;
    add_phase_time(arrays, phase_start);
    write_npz(phase_path(output_dir, cycle, phase_index + 1, phase.at("name").get<std::string>()), arrays,
              compress_output);
  }
  progress["global_time"] = phase_start + duration;
  progress["phase"] = phase_index + 1;
  progress["elapsed"] = 0.;
  progress["done"] = Json::object();
  progress["phase_start_time"] = progress["global_time"];
  progress["phase_ledger_start"] = slice(progress.at("state").get<Vec>(), physical_size, ledger_size);
  save_checkpoint();
  write_outputs();
}

void write_legacy_breakthrough_plots(const std::filesystem::path& input_path)
{
  const auto document = read_input_document(input_path);
  auto output = std::filesystem::path(document.value("Output", Json::object()).value("Directory", "results"));
  if (output.is_relative()) output = input_path.parent_path() / output;
  const auto script = find_plot_script(output);
  if (script.empty())
    throw std::runtime_error("Legacy plotting needs tools/plot_cycle_results.py; set RUPTURA_PLOT_SCRIPT");
  std::filesystem::create_directories(output / "plots");
  const auto log = output / "plots" / "plot_log.txt";
  const std::string command = shell_quote(python_executable()) + ' ' + shell_quote(script.string()) + ' ' +
      shell_quote(output.string()) + " --legacy-input " + shell_quote(input_path.string()) +
      " > " + shell_quote(log.string()) + " 2>&1";
  if (std::system(command.c_str()) != 0)
    throw std::runtime_error("Legacy plotting failed; see " + log.string() +
                             ". NumPy and Matplotlib must be available to " + python_executable());
  if (!std::filesystem::is_regular_file(output / "plots" / "plot_manifest.json"))
    throw std::runtime_error("Legacy plotting produced no plot manifest: " + log.string());
}

void CycleRunner::write_plots()
{
  if (plot_backend == "svg")
  {
    write_svg_plots();
    return;
  }
  std::string message;
  if (write_matplotlib_plots(&message)) return;
  if (plot_backend == "matplotlib")
    throw std::runtime_error("Matplotlib plotting was requested but failed: " + message +
                             "\nInstall NumPy and Matplotlib for " + python_executable() +
                             ", or select --plot-backend svg.");
  std::cerr << "Matplotlib plots unavailable (" << message << "); writing the built-in SVG overview instead."
            << std::endl;
  write_svg_plots();
}

bool CycleRunner::write_matplotlib_plots(std::string* message)
{
  auto report = [&](const std::string& text)
  {
    if (message) *message = text;
    return false;
  };
  if (!std::filesystem::is_regular_file(output_dir / "states.npz") &&
      !std::filesystem::is_regular_file(output_dir / "cycle_metrics.csv"))
    return report("no states.npz or cycle_metrics.csv has been written yet");
  std::filesystem::path script;
  try
  {
    script = find_plot_script(output_dir);
  }
  catch (const std::exception& e)
  {
    return report(e.what());
  }
  if (script.empty())
    return report(
        "tools/plot_cycle_results.py was not found; set RUPTURA_PLOT_SCRIPT or "
        "RUPTURA_CYCLE_CPP_ROOT");
  std::filesystem::create_directories(output_dir / "plots");
  const auto log = output_dir / "plots" / "plot_log.txt";
  const std::string command = shell_quote(python_executable()) + ' ' + shell_quote(script.string()) + ' ' +
                              shell_quote(output_dir.string()) + " > " + shell_quote(log.string()) + " 2>&1";
  const int status = std::system(command.c_str());
  if (status != 0)
  {
    std::string detail = "exit status " + std::to_string(status);
    std::ifstream stream(log);
    std::string line, last;
    while (std::getline(stream, line))
      if (!line.empty()) last = line;
    if (!last.empty()) detail += ": " + last;
    return report(detail);
  }
  const auto manifest = output_dir / "plots" / "plot_manifest.json";
  if (!std::filesystem::is_regular_file(manifest)) return report("the plot script wrote no manifest");
  return true;
}

void CycleRunner::write_svg_plots()
{
  if (!progress.is_object() || !progress.contains("history_state"))
    throw std::runtime_error("No recorded states available for plotting; run or load a checkpoint first");
  const auto times = progress.at("history_t").get<Vec>();
  if (times.empty()) throw std::runtime_error("Cannot plot an empty state history");
  auto xml = [](const std::string& text)
  {
    std::string result;
    for (char c : text)
    {
      if (c == '&')
        result += "&amp;";
      else if (c == '<')
        result += "&lt;";
      else if (c == '>')
        result += "&gt;";
      else if (c == '\"')
        result += "&quot;";
      else
        result += c;
    }
    return result;
  };
  const double width = 1000., panel_height = 240., height = 100. + network.column_names.size() * 2. * panel_height;
  std::ostringstream svg;
  svg << std::setprecision(12) << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1000\" height=\"" << height
      << "\" viewBox=\"0 0 1000 " << height
      << "\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/><g font-family=\"sans-serif\" "
         "font-size=\"13\" fill=\"#202832\"><text x=\"55\" y=\"28\" font-size=\"20\">RUPTURA "
         "Cycles: recorded column histories</text>";
  const auto summary = progress.value("summary", Json::object());
  if (summary.is_object() && summary.contains("convergence_report"))
  {
    const auto& report = summary["convergence_report"];
    svg << "<text x=\"55\" y=\"46\" font-size=\"11\" fill=\"#8a341c\">"
        << xml(report.value("status", std::string("unknown"))) << "; "
        << xml(report.value("selection_reason", std::string())) << "; selected cycle "
        << report.value("representative_cycle", Json()).dump() << "</text>";
    const auto stats = report.value("error_statistics", Json::object()).value("css_state_error", Json::object());
    if (stats.contains("min"))
      svg << "<text x=\"55\" y=\"64\" font-size=\"11\">Tail CSS error min/mean/max: " << stats["min"] << " / "
          << stats["mean"] << " / " << stats["max"] << "; pass &lt;= 1</text>";
  }
  const double xmin = times.front(), xmax = times.back() > times.front() ? times.back() : times.front() + 1.;
  std::size_t panel = 0;
  for (const auto& name : network.column_names)
  {
    Vec pressure, temperature;
    for (const auto& state : progress.at("history_state"))
    {
      auto d = network.model(name).derived(network.column_state(name, state.get<Vec>()));
      pressure.push_back(sum(d.P) / d.P.size() / 1e5);
      temperature.push_back(sum(d.Tg) / d.Tg.size());
    }
    for (int field = 0; field < 2; ++field)
    {
      const auto& values = field == 0 ? pressure : temperature;
      double lo = *std::min_element(values.begin(), values.end()), hi = *std::max_element(values.begin(), values.end());
      double pad = hi > lo ? .06 * (hi - lo) : std::max(std::abs(lo) * .02, 1.);
      lo -= pad;
      hi += pad;
      double top = 100. + panel * panel_height, left = 85., right = width - 35., bottom = top + 165.;
      auto sx = [&](double x) { return left + (right - left) * (x - xmin) / (xmax - xmin); };
      auto sy = [&](double y) { return bottom - (bottom - top) * (y - lo) / (hi - lo); };
      svg << "<text x=\"" << left << "\" y=\"" << top - 12. << "\" font-weight=\"bold\">" << xml(name) << ": "
          << (field == 0 ? "Mean bed pressure [bar]" : "Mean gas temperature [K]") << "</text>";
      for (int tick = 0; tick <= 4; ++tick)
      {
        double xx = xmin + (xmax - xmin) * tick / 4., yy = lo + (hi - lo) * tick / 4.;
        svg << "<path d=\"M " << left << ' ' << sy(yy) << " H " << right
            << "\" stroke=\"#e4e8ec\" fill=\"none\"/><text x=\"" << left - 8. << "\" y=\"" << sy(yy) + 4.
            << "\" text-anchor=\"end\">" << std::setprecision(5) << yy << std::setprecision(12) << "</text><text x=\""
            << sx(xx) << "\" y=\"" << bottom + 19. << "\" text-anchor=\"middle\">" << std::setprecision(5) << xx
            << std::setprecision(12) << "</text>";
      }
      svg << "<path d=\"M " << left << ' ' << top << " V " << bottom << " H " << right
          << "\" fill=\"none\" stroke=\"#536271\"/><polyline fill=\"none\" stroke=\""
          << (field == 0 ? "#1769aa" : "#b34528") << "\" stroke-width=\"1.7\" points=\"";
      for (std::size_t i = 0; i < times.size(); ++i) svg << sx(times[i]) << ',' << sy(values[i]) << ' ';
      svg << "\"/><text x=\"" << (left + right) / 2. << "\" y=\"" << bottom + 40.
          << "\" text-anchor=\"middle\">Elapsed time [s]</text>";
      ++panel;
    }
  }
  svg << "</g></svg>\n";
  const auto file = output_dir / "plots" / "cycle_overview.svg";
  atomic_text(file, svg.str());
  Json manifest = {{"files", Json::array({file.string()})},
                   {"backend", "native SVG"},
                   {"notes", Json::array({"Overview uses arithmetic means of the recorded axial node pressures "
                                          "and gas temperatures.",
                                          "Complete component, loading, temperature, pressure and phase "
                                          "histories are retained in states.npz and phases/*.npz."})}};
  atomic_text(output_dir / "plots" / "plot_manifest.json", manifest.dump(2));
}

}  // namespace cycle

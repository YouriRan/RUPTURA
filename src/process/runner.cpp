#include "cycle/runner.hpp"

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include "cycle/diagnostics.hpp"
#include "cycle/npz.hpp"
namespace cycle
{
namespace
{
const std::vector<std::string> directions = {"left_in", "left_out", "right_in", "right_out"};
const std::set<std::string> controls = {"duration",
                                        "step_time",
                                        "stop",
                                        "minimum_duration",
                                        "maximum_duration",
                                        "timeout_policy",
                                        "failure_maximum_simulated_time"};
constexpr double eps = std::numeric_limits<double>::epsilon();
void require(bool good, const std::string& message)
{
  if (!good) throw std::invalid_argument(message);
}
double num(const Json& j, const std::string& k, double def = 0)
{
  return j.contains(k) && !j[k].is_null() ? j[k].get<double>() : def;
}
std::string str(const Json& j, const std::string& k, const std::string& def = "")
{
  return j.contains(k) && j[k].is_string() ? j[k].get<std::string>() : def;
}
std::string lower(std::string s)
{
  for (auto& c : s) c = std::tolower(static_cast<unsigned char>(c));
  return s;
}
Vec vec(const Json& j) { return j.get<Vec>(); }
Vec slice(const Vec& v, size_t first, size_t count)
{
  require(first + count <= v.size(), "Invalid state slice");
  return Vec(v.begin() + first, v.begin() + first + count);
}
Vec subtract(const Vec& a, const Vec& b)
{
  require(a.size() == b.size(), "Vector dimensions differ");
  Vec r(a.size());
  for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] - b[i];
  return r;
}
double sum(const Vec& v) { return std::accumulate(v.begin(), v.end(), 0.0); }
void plus(Vec& a, const Vec& b, double factor = 1)
{
  require(a.size() == b.size(), "Vector dimensions differ");
  for (size_t i = 0; i < a.size(); ++i) a[i] += factor * b[i];
}
bool any(const Vec& v, double tol = 0)
{
  return std::any_of(v.begin(), v.end(), [tol](double x) { return std::abs(x) > tol; });
}
std::string joule_key(const std::string& key)
{
  return key.ends_with("_W") ? key.substr(0, key.size() - 2) + "_J" : key + "_J";
}
size_t component_index(const Json& c, const std::vector<std::string>& names)
{
  if (c.is_number_integer())
  {
    auto i = c.get<int>();
    require(i >= 0 && size_t(i) < names.size(), "Invalid component index");
    return size_t(i);
  }
  auto name = c.get<std::string>();
  auto it = std::find(names.begin(), names.end(), name);
  require(it != names.end(), "Unknown component: " + name);
  return size_t(it - names.begin());
}
void atomic_text(const std::filesystem::path& path, const std::string& text)
{
  std::filesystem::create_directories(path.parent_path());
  auto tmp = path;
  tmp += ".tmp." + std::to_string(getpid());
  try
  {
    {
      std::ofstream f(tmp, std::ios::binary);
      if (!f) throw std::runtime_error("Cannot write " + tmp.string());
      f << text;
      f.flush();
      if (!f) throw std::runtime_error("Cannot complete " + tmp.string());
    }
    int fd = ::open(tmp.c_str(), O_RDONLY);
    if (fd >= 0)
    {
      ::fsync(fd);
      ::close(fd);
    }
    std::filesystem::rename(tmp, path);
  }
  catch (...)
  {
    std::filesystem::remove(tmp);
    throw;
  }
}
Json read_json(const std::filesystem::path& p)
{
  std::ifstream f(p);
  if (!f) throw std::runtime_error("Cannot read " + p.string());
  Json j;
  f >> j;
  return j;
}
Json sorted_json(const Json& j)
{
  if (j.is_object())
  {
    Json r = Json::object();
    std::set<std::string> keys;
    for (const auto& kv : j.items()) keys.insert(kv.key());
    for (const auto& k : keys) r[k] = sorted_json(j[k]);
    return r;
  }
  if (j.is_array())
  {
    Json r = Json::array();
    for (const auto& v : j) r.push_back(sorted_json(v));
    return r;
  }
  return j;
}
std::string fingerprint_of(const Json& j)
{
  uint64_t h = 14695981039346656037ull;
  for (unsigned char c : sorted_json(j).dump())
  {
    h ^= c;
    h *= 1099511628211ull;
  }
  std::ostringstream s;
  s << "native-fnv1a-" << std::hex << h;
  return s.str();
}
void collect_strings(const Json& value, const std::function<void(const std::string&, const Json&)>& fn)
{
  if (value.is_object())
  {
    for (const auto& [key, item] : value.items())
    {
      fn(key, item);
      collect_strings(item, fn);
    }
  }
  else if (value.is_array())
    for (const auto& item : value) collect_strings(item, fn);
}
}  // namespace
StageRule StageRule::from_mapping(const Json& data, double default_failure)
{
  StageRule r;
  if (data.contains("duration") && !data["duration"].is_null())
    r.duration = num(data, "duration");
  else if (data.contains("step_time") && !data["step_time"].is_null())
    r.duration = num(data, "step_time");
  r.stop = data.value("stop", Json());
  r.minimum = num(data, "minimum_duration");
  if (data.contains("maximum_duration") && !data["maximum_duration"].is_null())
    r.maximum = num(data, "maximum_duration");
  r.failure_time = num(data, "failure_maximum_simulated_time", default_failure);
  r.timeout_policy = lower(str(data, "timeout_policy", "error"));
  require(!r.duration || (std::isfinite(*r.duration) && *r.duration > 0), "duration must be finite and positive");
  require(std::isfinite(r.minimum) && r.minimum >= 0, "minimum_duration must be finite and nonnegative");
  require(!r.maximum || (std::isfinite(*r.maximum) && *r.maximum > 0 && *r.maximum >= r.minimum),
          "maximum_duration must be positive and at least minimum_duration");
  require(!r.duration || (r.stop.is_null() && r.minimum == 0 && !r.maximum),
          "Use duration for a timed stage or stop/minimum_duration/maximum_duration for a hybrid "
          "stage");
  require(r.timeout_policy == "error" || r.timeout_policy == "advance", "timeout_policy must be error or advance");
  require(std::isfinite(r.failure_time) && r.failure_time > 0,
          "failure_maximum_simulated_time must be finite and positive");
  require(!r.maximum || !r.stop.is_null(), "maximum_duration requires stop; use duration for timed stages");
  return r;
}
double StageRule::next_boundary(double elapsed) const
{
  double value = failure_time;
  if (duration) value = std::min(value, *duration);
  if (minimum > elapsed + 16 * eps * std::max({std::abs(elapsed), std::abs(minimum), 1e-300}))
    value = std::min(value, minimum);
  if (maximum) value = std::min(value, *maximum);
  return value;
}
CycleRunner::CycleRunner(const Json& config) : cfg(config), network(config)
{
  physical_size = network.size;
  require(network.models.size() >= 1 && network.models.size() <= 4, "RUPTURA cycles supports one to four columns");
  canonicalize_phases();
  validate_phases();
  require(!cfg.contains("SOLVER_METHOD") || cfg["SOLVER_METHOD"].is_string(),
          "BreakthroughIntegrator must be a solver name, not null or another type.");
  auto method = lower(str(cfg, "SOLVER_METHOD", "BDF"));
  const std::map<std::string,std::string> methods={{"bdf","BDF"},{"bdf1","BDF1"},{"cvode","CVODE"},
    {"ssprk3","SSPRK3"},{"semiimplicitssprk3","SemiImplicitSSPRK3"}};
  require(methods.contains(method),"Unsupported shared time integrator: "+method);
  solver_method = methods.at(method);
  if (solver_method == "CVODE" && !cvode_available())
    throw std::runtime_error(
        "CVODE was requested, but this executable was built without "
        "SUNDIALS CVODE support. Build with CVODE enabled or select BDF "
        "explicitly; no fallback was performed.");
  // The effective method participates in checkpoint compatibility. A resumed
  // run must retain its original solver; a fresh run may select either method.
  cfg["SOLVER_METHOD"] = solver_method;
  output_dir = std::filesystem::absolute(str(cfg, "OUTPUT_DIR", "results"));
  checkpoint_path = output_dir / "checkpoint_latest.json";
  // Plot backend selection never changes the simulation and is deliberately
  // kept out of the YAML, so one input file serves both package editions.
  const char* backend_env = std::getenv("RUPTURA_PLOT_BACKEND");
  plot_backend = lower(str(cfg, "PLOT_BACKEND", backend_env ? backend_env : "auto"));
  if (plot_backend.empty()) plot_backend = "auto";
  require(plot_backend == "auto" || plot_backend == "matplotlib" || plot_backend == "svg",
          "Plot backend must be auto, matplotlib or svg");
  const char* env = std::getenv("RUPTURA_RESTART");
  std::string environ = env ? lower(env) : "";
  restart = cfg.value("RESTART", false) || cfg.value("RESTART_FROM_CHECKPOINT", false) ||
            !str(cfg, "RESTART_FROM").empty() || environ == "1" || environ == "true" || environ == "yes";
  restart_path =
      str(cfg, "RESTART_FROM").empty() ? checkpoint_path : std::filesystem::absolute(str(cfg, "RESTART_FROM"));
  run_mode = lower(str(cfg, "RUN_MODE", "css"));
  require(run_mode == "css" || run_mode == "fixed_cycles", "RUN_MODE must be css or fixed_cycles");
  cycle_limit =
      int(num(cfg, run_mode == "fixed_cycles" ? "N_CYCLES" : "MAX_CYCLES", run_mode == "fixed_cycles" ? 1 : 200));
  required_streak = int(num(cfg, "CSS_CONSECUTIVE_CYCLES", 3));
  require(cycle_limit > 0 && required_streak > 0, "Cycle counts must be positive");
  assess_convergence(Json::array(), cfg);  // validates reporting options for both run modes
  for (const char* key :
       {"CSS_STATE_RTOL", "CSS_STATE_ATOL", "CSS_INVENTORY_RTOL", "CSS_METRIC_RTOL", "CSS_MASS_BALANCE_TOLERANCE",
        "CSS_DUTY_ATOL_J", "CSS_DURATION_ATOL_S", "CSS_PRODUCTIVITY_ATOL"})
    if (cfg.contains(key))
      require(std::isfinite(num(cfg, key)) && num(cfg, key) > 0, std::string(key) + " must be finite positive");
  segment_time = num(cfg, "CHECKPOINT_INTERVAL_S", 5);
  require(std::isfinite(segment_time) && segment_time > 0, "CHECKPOINT_INTERVAL_S must be finite and positive");
  const auto output_mode = lower(str(cfg, "OUTPUT_PROFILE", "full"));
  require(output_mode == "compact" || output_mode == "full", "OUTPUT_PROFILE must be compact or full");
  compact_output = output_mode == "compact";
  checkpoint_enabled = cfg.value("CHECKPOINT_ENABLED", true);
  phase_archives = cfg.value("WRITE_PHASE_ARCHIVES", !compact_output);
  selected_archive = cfg.value("WRITE_SELECTED_CYCLE", !compact_output);
  compress_output = cfg.value("COMPRESS_NPZ", compact_output);
  for (const auto* key : {"OUTPUT_SAMPLES_PER_PHASE", "OUTPUT_HISTORY_CYCLES"})
    if (cfg.contains(key))
    {
      const double value = num(cfg, key);
      require(std::isfinite(value) && value >= 0 && value < std::numeric_limits<int>::max() - 2 &&
                  value == std::floor(value),
              std::string(key) + " must be a nonnegative integer within the supported range");
    }
  output_interval = num(cfg, "OUTPUT_SAMPLE_INTERVAL_S", 5);
  output_samples = cfg.value("OUTPUT_SAMPLES_PER_PHASE", 200);
  require(std::isfinite(output_interval) && output_interval > 0 && output_samples >= 2,
          "Output sample interval must be positive and samples_per_phase >= 2");
  history_cycles = std::max({cfg.value("OUTPUT_HISTORY_CYCLES", 0), cfg.value("CSS_DIAGNOSTIC_WINDOW", 20),
                             cfg.value("CSS_MAXIMUM_PERIOD", 4)}) +
                   2;
  init_ledger();
  atol = network.state_atol(num(cfg, "SOLVER_ATOL", 1e-9));
  atol.insert(atol.end(), ledger_size, num(cfg, "LEDGER_ATOL", 1e-10));
  rtol = num(cfg, "SOLVER_RTOL", 1e-6);
  const std::set<std::string> excluded = {"OUTPUT_DIR",
                                          "RESTART",
                                          "RESTART_FROM",
                                          "RESTART_FROM_CHECKPOINT",
                                          "ALLOW_EXISTING_OUTPUT_DIR",
                                          "N_CYCLES",
                                          "MAX_CYCLES",
                                          "RUN_MODE",
                                          "CHECKPOINT_INTERVAL_S",
                                          "PLOT",
                                          "AUTO_PLOT",
                                          "PLOT_BACKEND",
                                          "INPUT_DOCUMENT",
                                          "CASE_DIR",
                                          "INPUT_FILE",
                                          "MAKE_PLOTS",
                                          "NOTES",
                                          "CASE_NAME",
                                          "DESCRIPTION"};
  restart_config_ = Json::object();
  for (const auto& [key, value] : cfg.items())
    if (!key.starts_with("_") && !excluded.count(key) && key != "OUTPUT_PROFILE" && key != "OUTPUT_SAMPLE_INTERVAL_S" &&
        key != "OUTPUT_SAMPLES_PER_PHASE" && key != "OUTPUT_HISTORY_CYCLES" && key != "CHECKPOINT_ENABLED" &&
        key != "WRITE_PHASE_ARCHIVES" && key != "WRITE_SELECTED_CYCLE" && key != "COMPRESS_NPZ")
      restart_config_[key] = value;
  restart_config_ = sorted_json(restart_config_);
  fingerprint = fingerprint_of(restart_config_);
}
void CycleRunner::canonicalize_phases()
{
  if (cfg.contains("CYCLE_PHASES") && !cfg["CYCLE_PHASES"].empty())
    phases = cfg["CYCLE_PHASES"];
  else if (cfg.contains("TWO_COLUMN_PHASES") && !cfg["TWO_COLUMN_PHASES"].empty())
    phases = cfg["TWO_COLUMN_PHASES"];
  else
  {
    phases = Json::array();
    size_t i = 0;
    for (const auto& step : cfg.value("CYCLE_STEPS", Json::array()))
    {
      Json phase = Json::object(), action = Json::object();
      for (const auto& [key, value] : step.items())
        if (controls.count(key))
          phase[key] = value;
        else
          action[key] = value;
      phase["name"] = step.value("name", Json("step_" + std::to_string(++i)));
      phase["columns"] = {{network.column_names.front(), action}};
      phases.push_back(phase);
    }
  }
  size_t i = 0;
  for (auto& phase : phases)
  {
    ++i;
    if (!phase.contains("name")) phase["name"] = "phase_" + std::to_string(i);
    if (!phase.contains("columns")) phase["columns"] = Json::object();
    if (phase.contains("connection") && !phase.contains("connections"))
      phase["connections"] = Json::array({phase["connection"]});
    phase.erase("connection");
    if (!phase.contains("connections")) phase["connections"] = Json::array();
    for (const auto& name : network.column_names)
      if (!phase["columns"].contains(name))
        phase["columns"][name] = {
            {"name", "hold"}, {"left_port", {{"mode", "closed"}}}, {"right_port", {{"mode", "closed"}}}};
  }
}
void CycleRunner::validate_criterion(const Json& node) const
{
  require(node.is_object(), "Stop definition must be a mapping");
  if (node.contains("logic"))
  {
    auto logic = lower(str(node, "logic"));
    require(logic == "all" || logic == "any", "Stop logic must be all or any");
    require(node.contains("criteria") && node["criteria"].is_array() && !node["criteria"].empty(),
            "Logical stop requires nonempty criteria");
    for (const auto& c : node["criteria"]) validate_criterion(c);
    return;
  }
  require(lower(str(node, "type")) != "callable", "Python callable criteria are not available in native inputs");
  require(node.contains("metric") && node.contains("value"), "Stop criterion requires metric and value");
  auto op = str(node, "operator", ">=");
  require(op == ">=" || op == ">" || op == "<=" || op == "<", "Unsupported criterion operator");
  require(std::isfinite(num(node, "value")), "Criterion target must be finite");
  if (node.contains("scale"))
    require(std::isfinite(num(node, "scale")) && std::abs(num(node, "scale")) > 0,
            "Criterion scale must be finite and positive");
}
void CycleRunner::validate_phases()
{
  require(phases.is_array() && !phases.empty(), "Provide at least one cycle phase or step");
  std::set<std::string> names;
  double failure = num(cfg, "STEP_FAILURE_MAXIMUM_SIMULATED_TIME", 2e5);
  for (auto& phase : phases)
  {
    require(names.insert(str(phase, "name")).second, "Cycle phase names must be unique");
    require(phase["columns"].size() == network.models.size(), "Phase contains unknown column names");
    auto rule = StageRule::from_mapping(phase, failure);
    bool action_enabled = false;
    for (auto& [name, action] : phase["columns"].items())
    {
      require(network.column_indices.count(name), "Unknown column " + name);
      auto r = StageRule::from_mapping(action, failure);
      action_enabled |= r.enabled();
      if (!action.contains("left_port")) action["left_port"] = {{"mode", "closed"}};
      if (!action.contains("right_port")) action["right_port"] = {{"mode", "closed"}};
      if (!r.stop.is_null()) validate_criterion(r.stop);
    }
    require(rule.enabled() || action_enabled, "Phase requires duration, stop, or action stopping rules");
    if (!phase["connections"].empty())
    {
      require(rule.enabled(), "Connected phases require shared phase duration or stop");
      require(!action_enabled, "Connected columns require phase termination, not independent closure");
    }
    if (!rule.stop.is_null()) validate_criterion(rule.stop);
  }
}
void CycleRunner::init_ledger()
{
  std::set<std::string> tagset(network.stream_tags.begin(), network.stream_tags.end());
  for (const auto& t : {"feed", "light_product", "heavy_product", "product", "purge", "repressurization_gas", "waste"})
    tagset.insert(t);
  for (const auto& key :
       {"FEED_TAG", "LIGHT_PRODUCT_TAG", "HEAVY_PRODUCT_TAG", "PURGE_TAG", "REPRESSURIZATION_GAS_TAG"})
    if (cfg.contains(key)) tagset.insert(str(cfg, key));
  auto collect = [&](const std::string& key, const Json& value)
  {
    if (key == "stream_tags" && value.is_object())
      for (const auto& v : value)
        if (v.is_string()) tagset.insert(v.get<std::string>());
    if ((key.ends_with("tag") || key == "composition_source") && value.is_string())
      tagset.insert(value.get<std::string>());
  };
  collect_strings(phases, collect);
  collect_strings(cfg.value("SOURCES", Json::object()), collect);
  collect_strings(cfg.value("DYNAMIC_BUFFERS", Json::object()), collect);
  for (const auto& name : network.column_names)
    for (const auto& d : directions) tagset.insert(name + "." + d);
  tags.assign(tagset.begin(), tagset.end());
  std::set<std::string> uset(ColumnNetwork::utility_names.begin(), ColumnNetwork::utility_names.end());
  for (const auto& u : {"net_heat_W", "heating_W", "cooling_W", "electric_heater_W", "gas_heating_W", "gas_cooling_W",
                        "heat_loss_W", "compressor_W", "vacuum_W", "buffer_heating_W", "buffer_cooling_W"})
    uset.insert(u);
  utility_names.assign(uset.begin(), uset.end());
  auto add = [&](std::vector<std::string> key, size_t count)
  {
    ledger_slices[key] = {ledger_size, count};
    ledger_size += count;
  };
  for (const auto& name : network.column_names)
    for (const auto& d : directions) add({"port", name, d}, network.n_comp);
  for (const auto& kind : {"tag", "external_tag_in", "external_tag_out", "recycle"})
    for (const auto& tag : tags) add({kind, tag}, network.n_comp);
  add({"external_in"}, network.n_comp);
  add({"external_out"}, network.n_comp);
  if(std::any_of(network.models.begin(),network.models.end(),[](const auto& m){return m.has_reactions();})){add({"reaction"},network.n_comp);add({"energy","reaction_energy_W"},1);}
  if (network.models.front().liquid)
  {
    for (const auto direction : {"in", "out"})
    {
      add({std::string("volume_external_") + direction}, 1);
      for (const auto& tag : tags) add({std::string("volume_tag_") + direction, tag}, 1);
    }
  }
  for (const auto& u : utility_names) add({"utility", u}, 1);
  for (const auto& col : network.column_names)
    for (const auto& u : utility_names) add({"column_utility", col, u}, 1);
  for (const auto& k : {"external_enthalpy_in_W", "external_enthalpy_out_W"}) add({"energy", k}, 1);
}
Vec CycleRunner::ledger_value(const Vec& delta, const std::vector<std::string>& key) const
{
  auto it = ledger_slices.find(key);
  require(it != ledger_slices.end(), "Unknown ledger key");
  return slice(delta, it->second.first, it->second.second);
}
Vec CycleRunner::ledger_rates(const Json& diagnostics) const
{
  Vec result(ledger_size, 0), zero(network.n_comp, 0);
  auto add = [&](std::vector<std::string> key, const Json& value)
  {
    Vec values = value.is_array() ? vec(value) : Vec{value.get<double>()};
    auto it = ledger_slices.find(key);
    if (it == ledger_slices.end())
    {
      require(!any(values), "Undeclared stream or utility; include its tag in configuration");
      return;
    }
    auto [offset, count] = it->second;
    require(values.size() == count, "Ledger component count differs");
    for (size_t i = 0; i < count; ++i) result[offset + i] += values[i];
  };
  const Json ports = diagnostics.value("port_rates", Json::object());
  for (const auto& [name, rates] : ports.items())
    for (const auto& direction : directions) add({"port", name, direction}, rates.value(direction, Json(zero)));
  Json provenance = diagnostics.value("provenance", diagnostics.value("source_provenance", Json::array()));
  const Json streams = diagnostics.value("external_stream_rates", Json::object());
  if (provenance.empty())
    for (const auto& [tag, rates] : streams.items()) add({"tag", tag}, rates);
  add({"external_in"}, diagnostics.value("external_in", Json(zero)));
  add({"external_out"}, diagnostics.value("external_out", Json(zero)));
  add({"reaction"},diagnostics.value("reaction_component_rates",Json(zero)));
  add({"energy","reaction_energy_W"},diagnostics.value("reaction_energy_W",Json(0.)));
  for (const auto& item : provenance)
  {
    auto origin = str(item, "origin", "external"), direction = str(item, "direction"),
         tag = str(item, "tag", direction == "in" ? "feed" : "waste");
    auto values = item.value("component_rates", Json(zero));
    add({"tag", tag}, values);
    if ((origin == "external" || origin == "product_recycle") && (direction == "in" || direction == "out"))
      add({"external_tag_" + direction, tag}, values);
    if (origin == "product_recycle" && direction == "in") add({"recycle", str(item, "product_tag", tag)}, values);
  }
  if (network.models.front().liquid)
    for (auto& entry : diagnostics.value("volume_provenance", Json::array()))
    {
      auto direction = entry.at("direction").get<std::string>();
      add({"volume_external_" + direction}, entry.at("rate_m3_s"));
      add({"volume_tag_" + direction, entry.at("tag").get<std::string>()}, entry.at("rate_m3_s"));
    }
  const Json utilities = diagnostics.value("utility_rates_W", Json::object());
  for (const auto& [name, rate] : utilities.items()) add({"utility", name}, rate);
  const Json column_utilities = diagnostics.value("column_utility_rates_W", Json::object());
  for (const auto& [column, rates] : column_utilities.items())
    for (const auto& [name, rate] : rates.items()) add({"column_utility", column, name}, rate);
  for (const auto& name : {"external_enthalpy_in_W", "external_enthalpy_out_W"})
    add({"energy", name}, diagnostics.value(name, Json(0.0)));
  return result;
}
Vec CycleRunner::rhs(double t, const Vec& augmented)
{
  auto [physical, diagnostics] = network.evaluate(t, slice(augmented, 0, physical_size));
  auto ledger = ledger_rates(diagnostics);
  physical.insert(physical.end(), ledger.begin(), ledger.end());
  for (double value : physical)
    if (!std::isfinite(value)) throw std::runtime_error("Nonfinite physical or ledger derivative");
  return physical;
}
Json CycleRunner::criterion_context(double t, const Vec& augmented, const std::string& column) const
{
  auto delta = subtract(slice(augmented, physical_size, ledger_size), phase_ledger_start_);
  Json context = active_context_, columns = Json::object();
  auto selected = column.empty() ? network.column_names.front() : column;
  for (const auto& name : network.column_names)
  {
    Json ds = Json::object(), ts = Json::object();
    const auto& action = current_phase_["columns"][name];
    for (const auto& d : directions)
    {
      auto amounts = ledger_value(delta, {"port", name, d});
      ds[d] = amounts;
      auto tag = str(action.value("stream_tags", Json::object()), d, name + "." + d);
      auto existing = ts.contains(tag) ? vec(ts[tag]) : Vec(network.n_comp, 0);
      plus(existing, amounts);
      ts[tag] = existing;
    }
    columns[name] = {{"directions", ds}, {"tags", ts}};
  }
  context["column_accumulators"] = columns;
  if (!column.empty())
    context["step_accumulators"] = columns[column];
  else
  {
    Json alltags = Json::object();
    for (const auto& tag : tags) alltags[tag] = ledger_value(delta, {"tag", tag});
    context["step_accumulators"] = {{"directions", columns[selected]["directions"]}, {"tags", alltags}};
  }
  context["column_utilities"] = Json::object();
  for (const auto& name : network.column_names)
    for (const auto& u : utility_names)
      context["column_utilities"][name][joule_key(u)] = ledger_value(delta, {"column_utility", name, u})[0];
  if (!column.empty())
    context["step_utilities"] = context["column_utilities"][column];
  else
    for (const auto& u : utility_names)
      context["step_utilities"][joule_key(u)] = ledger_value(delta, {"utility", u})[0];
  context["elapsed_time"] = t;
  return context;
}
double CycleRunner::criterion_actual(const Json& node, double t, const Vec& augmented, const std::string& column)
{
  auto key = lower(str(node, "metric"));
  auto context = criterion_context(t, augmented, column);
  auto physical = slice(augmented, 0, physical_size);
  const std::set<std::string> cumulative = {"cumulative_stream_moles", "cumulative_component_moles",
                                            "cumulative_port_moles", "bed_volumes_processed"};
  const std::map<std::string, std::string> utilities = {{"cumulative_net_heat", "net_heat_J"},
                                                        {"cumulative_heating_energy", "heating_J"},
                                                        {"cumulative_cooling_energy", "cooling_J"},
                                                        {"cumulative_electric_heater_energy", "electric_heater_J"},
                                                        {"cumulative_gas_heating_energy", "gas_heating_J"},
                                                        {"cumulative_gas_cooling_energy", "gas_cooling_J"},
                                                        {"cumulative_heat_loss", "heat_loss_J"}};
  if (cumulative.count(key))
  {
    auto explicit_column = str(node, "column");
    Json accum =
        explicit_column.empty() ? context["step_accumulators"] : context["column_accumulators"][explicit_column];
    Vec values(network.n_comp, 0);
    if (key != "cumulative_port_moles")
    {
      auto tag = str(node, "stream", str(node, "tag"));
      require(!tag.empty(), key + " requires stream or tag");
      if (accum["tags"].contains(tag)) values = vec(accum["tags"][tag]);
    }
    else
    {
      auto flow_key = str(node, "flow_key");
      if (flow_key.empty())
      {
        auto port = lower(str(node, "port")), direction = lower(str(node, "direction"));
        require((port == "left" || port == "right") && (direction == "in" || direction == "out"),
                "cumulative_port_moles requires flow_key or port and direction");
        flow_key = port + "_" + direction;
      }
      if (accum["directions"].contains(flow_key)) values = vec(accum["directions"][flow_key]);
    }
    double amount = node.contains("component") && !node["component"].is_null()
                        ? values[component_index(node["component"], network.comp_names)]
                        : sum(values);
    if (key == "bed_volumes_processed")
    {
      const auto& model = network.model(!explicit_column.empty() ? explicit_column
                                        : column.empty()         ? network.column_names.front()
                                                                 : column);
      double pressure = num(node, "reference_pressure", num(model.cfg, "P_FEED", num(model.cfg, "P_HIGH", 1e5))),
             temperature =
                 num(node, "reference_temperature", num(model.cfg, "T_feed", num(model.cfg, "T_init", 298.15)));
      auto basis = lower(str(node, "volume_basis", "void"));
      double volume = 0;
      if (basis == "void" || basis == "bed_void" || basis == "void_volume")
      {
        for (size_t z = 0; z < model.Nz; ++z) volume += model.area * model.widths[z] * model.epsilon_node[z];
      }
      else if (basis == "total" || basis == "bed" || basis == "column")
        volume = model.area * model.L;
      else
        throw std::invalid_argument("volume_basis must be void or total");
      return amount * model.R * temperature / std::max(pressure, 1e-300) / std::max(volume, 1e-300);
    }
    return amount;
  }
  if (utilities.count(key))
  {
    auto explicit_column = str(node, "column");
    return explicit_column.empty() ? num(context["step_utilities"], utilities.at(key))
                                   : num(context["column_utilities"][explicit_column], utilities.at(key));
  }
  // Buffer criteria belong to the shared finite-buffer network even when a
  // local column action owns the stopping rule.
  if (key.starts_with("buffer_"))
  {
    require(network.buffers && network.buffers->enabled, key + " requires configured finite buffers");
    return network.buffers->metric(key, slice(physical, network.buffer_start, network.buffers->size),
                                   str(node, "buffer"), node.value("component", Json()));
  }
  if (key == "elapsed_time" || key == "stage_time" || key == "time") return t;
  Json opts = node;
  for (const auto& k : {"metric", "operator", "value", "scale", "label"}) opts.erase(k);
  if (!opts.contains("column") && !column.empty()) opts["column"] = column;
  return network.metric(str(node, "metric"), physical, t, opts);
}
double CycleRunner::criterion_residual(const Json& node, double t, const Vec& augmented, const std::string& column)
{
  if (node.contains("logic"))
  {
    bool all = lower(str(node, "logic")) == "all";
    double value = all ? std::numeric_limits<double>::infinity() : -std::numeric_limits<double>::infinity();
    for (const auto& child : node["criteria"])
    {
      double r = criterion_residual(child, t, augmented, column);
      value = all ? std::min(value, r) : std::max(value, r);
    }
    return value;
  }
  double target = num(node, "value"),
         scale = node.contains("scale") ? std::abs(num(node, "scale")) : std::max(std::abs(target), 1.0),
         actual = criterion_actual(node, t, augmented, column);
  require(std::isfinite(actual), "Criterion metric is nonfinite");
  auto op = str(node, "operator", ">=");
  return op == ">=" || op == ">" ? (actual - target) / scale : (target - actual) / scale;
}
bool CycleRunner::criterion_satisfied_roundoff(const Json& node, double t, const Vec& augmented,
                                               const std::string& column)
{
  if (node.contains("logic"))
  {
    bool all = lower(str(node, "logic")) == "all";
    for (const auto& c : node["criteria"])
    {
      bool good = criterion_satisfied_roundoff(c, t, augmented, column);
      if (all && !good) return false;
      if (!all && good) return true;
    }
    return all;
  }
  double target = num(node, "value"),
         scale = node.contains("scale") ? std::abs(num(node, "scale")) : std::max(std::abs(target), 1.0);
  return criterion_residual(node, t, augmented, column) >= -32 * eps * std::abs(target) / scale;
}
std::string CycleRunner::rule_status(const StageRule& rule, double elapsed, const Vec& augmented,
                                     const std::string& column)
{
  auto tol = [&](double target) { return 16 * eps * std::max({std::abs(elapsed), std::abs(target), 1e-300}); };
  if (rule.duration && elapsed >= *rule.duration - tol(*rule.duration)) return "duration";
  if (!rule.stop.is_null() && elapsed >= rule.minimum - tol(rule.minimum) &&
      criterion_satisfied_roundoff(rule.stop, elapsed, augmented, column))
    return "criterion";
  if (rule.maximum && elapsed >= *rule.maximum - tol(*rule.maximum))
  {
    if (rule.timeout_policy == "advance") return "maximum_duration";
    throw std::runtime_error("Physical criterion was not reached by maximum_duration (timeout_policy=error)");
  }
  if (elapsed >= rule.failure_time - tol(rule.failure_time))
    throw std::runtime_error("Physical criterion was not reached before failure_maximum_simulated_time");
  return "";
}
void CycleRunner::validate_state(const Vec& physical) const
{
  require(physical.size() == physical_size, "Accepted physical state has incompatible dimension");
  for (double v : physical)
    if (!std::isfinite(v)) throw std::runtime_error("Accepted physical state contains nonfinite values");
  for (size_t c = 0; c < network.models.size(); ++c)
  {
    const auto& model = network.models[c];
    size_t offset = network.column_starts.at(network.column_names[c]),
           nonnegative = model.state_size() - (3 + model.layout.transient_htf) * model.Nz;
    for (size_t i = 0; i < nonnegative; ++i)
      if (physical[offset + i] < -20 * atol[offset + i])
        throw std::runtime_error("Accepted state in column " + network.column_names[c] +
                                 " has negative concentration, loading or kinetic history at row " + std::to_string(i) +
                                 "; value=" + Json(physical[offset + i]).dump() +
                                 "; atol=" + Json(atol[offset + i]).dump() + "; no state was clipped or committed");
    for (size_t i = nonnegative; i < model.state_size(); ++i)
      if (physical[offset + i] <= 0) throw std::runtime_error("Accepted column state has nonpositive temperature");
  }
  if (network.buffers && network.buffers->enabled)
  {
    size_t block = network.buffers->block_size;
    for (size_t offset = network.buffer_start; offset < physical_size; offset += block)
    {
      for (size_t c = 0; c < network.n_comp; ++c)
        if (physical[offset + c] < -20 * atol[offset + c])
          throw std::runtime_error("Finite buffer exhausted; reduce withdrawal, add supply or terminate on inventory");
      if (physical[offset + network.n_comp] <= 0)
        throw std::runtime_error("Accepted buffer has nonpositive temperature");
    }
  }
}
void CycleRunner::activate_phase(const Json& phase, const Json& context, const Json& done)
{
  current_phase_ = phase;
  for (const auto& [name, value] : done.items())
    if (network.column_indices.count(name))
    {
      current_phase_["columns"][name].erase("flow_control");
      current_phase_["columns"][name]["left_port"] = {{"mode", "closed"}};
      current_phase_["columns"][name]["right_port"] = {{"mode", "closed"}};
      if (phase["columns"][name].contains("wait_action"))
        current_phase_["columns"][name].update(phase["columns"][name]["wait_action"]);
    }
  active_context_ = context;
  network.set_phase(current_phase_, context);
}
bool CycleRunner::needs_pressure_anchor(const Json& phase)
{
  for (const auto& action : phase["columns"])
    for (const auto& side : {"left_port", "right_port"})
    {
      auto port = action.value(side, Json::object());
      auto pressure = port.value("pressure", Json());
      if (pressure.is_object() && str(pressure, "initial_value") == "phase_start_pressure") return true;
    }
  return false;
}
PhaseSolution CycleRunner::integrate_phase(const Vec& physical_state, const Json& phase, Json context, Vec augmented,
                                           double elapsed, Json done, Vec phase_ledger_start, const PhaseCommit& commit)
{
  if (augmented.empty())
  {
    augmented = physical_state;
    augmented.resize(physical_size + ledger_size, 0);
  }
  require(augmented.size() == physical_size + ledger_size, "Augmented state dimension mismatch");
  validate_state(slice(augmented, 0, physical_size));
  if (needs_pressure_anchor(phase) && !context.contains("phase_start_pressures"))
  {
    require(elapsed == 0, "Resuming anchored pressure ramp requires original phase_start_pressures in context");
    context["phase_start_pressures"] = network.endpoint_pressures(slice(augmented, 0, physical_size));
  }
  if (phase_ledger_start.empty()) phase_ledger_start = slice(augmented, physical_size, ledger_size);
  phase_ledger_start_ = phase_ledger_start;
  double failure = num(cfg, "STEP_FAILURE_MAXIMUM_SIMULATED_TIME", 2e5);
  auto phase_rule = StageRule::from_mapping(phase, failure);
  std::map<std::string, StageRule> action_rules;
  for (const auto& [name, action] : phase["columns"].items())
    action_rules[name] = StageRule::from_mapping(action, failure);
  activate_phase(phase, context, done);
  PhaseSolution result;
  result.times.push_back(elapsed);
  result.augmented.push_back(augmented);
  while (result.reason.empty())
  {
    if (done.contains(":phase")) result.reason = str(done[":phase"], "reason");
    if (result.reason.empty() && phase_rule.enabled()) result.reason = rule_status(phase_rule, elapsed, augmented);
    bool changed = false;
    for (const auto& name : network.column_names)
    {
      const auto& rule = action_rules.at(name);
      if (done.contains(name) || !rule.enabled()) continue;
      auto status = rule_status(rule, elapsed, augmented, name);
      if (!status.empty())
      {
        done[name] = {{"reason", status}, {"elapsed_s", elapsed}};
        changed = true;
      }
    }
    if (!phase_rule.enabled())
    {
      bool all = true;
      for (const auto& [name, rule] : action_rules) all &= done.contains(name) || !rule.enabled();
      if (all) result.reason = "all_actions_completed";
    }
    if (!result.reason.empty()) break;
    if (changed) activate_phase(phase, context, done);
    double next_time = std::min(elapsed + segment_time, phase_rule.next_boundary(elapsed));
    std::vector<OdeEvent> events;
    std::vector<std::string> event_columns;
    auto add_rule = [&](const StageRule& rule, const std::string& column)
    {
      next_time = std::min(next_time, rule.next_boundary(elapsed));
      if (!rule.stop.is_null() &&
          elapsed >= rule.minimum - 16 * eps * std::max({std::abs(elapsed), std::abs(rule.minimum), 1e-300}))
      {
        auto definition = rule.stop;
        events.push_back({column.empty() ? ":phase" : column, [this, definition, column](double t, const Vec& state)
                          { return criterion_residual(definition, t, state, column); }, 1, true});
        event_columns.push_back(column);
      }
    };
    if (phase_rule.enabled()) add_rule(phase_rule, "");
    for (const auto& name : network.column_names)
      if (!done.contains(name) && action_rules.at(name).enabled()) add_rule(action_rules.at(name), name);
    if (next_time <= elapsed) throw std::runtime_error("Phase exhausted its simulated-time failure guard");
    BdfOptions opts;
    opts.relative_tolerance = rtol;
    opts.absolute_tolerance = atol;
    opts.maximum_step = num(cfg, "SOLVER_MAX_STEP", std::numeric_limits<double>::infinity());
    opts.physical_dimension = physical_size;
    opts.first_step = num(cfg, "SOLVER_FIRST_STEP", 0.);
    opts.apply_nonnegative_constraints=cfg.value("SOLVER_NONNEGATIVE_CONSTRAINTS",true);
    for(std::size_t i=0;i<physical_size;++i)opts.nonnegative_indices.push_back(i);
    if (opts.apply_nonnegative_constraints)
    {
      opts.constraint_offsets.assign(physical_size + ledger_size, 0.);
      for (std::size_t i = 0; i < physical_size; ++i)
      {
        opts.constraint_offsets[i] = 10.*atol[i];
      }
    }
    opts.cvode_spgmr = cfg.value("CVODE_LINEAR_SOLVER", std::string("Dense")) == "SPGMR";
    for (const auto& m : network.models)
      opts.cvode_preconditioner_bandwidth = std::max(opts.cvode_preconditioner_bandwidth, m.state_size() - 1);
    if (opts.cvode_spgmr)
    {
      opts.state_permutation.resize(physical_size + ledger_size);
      std::iota(opts.state_permutation.begin(), opts.state_permutation.end(), 0);
      opts.cvode_preconditioner_bandwidth = 0;
      for (const auto& name : network.column_names)
      {
        const auto& model = network.model(name);
        const auto start = network.column_starts.at(name), blocks = model.layout.n_blocks();
        for (std::size_t node = 0; node < model.Nz; ++node)
          for (std::size_t block = 0; block < blocks; ++block)
            opts.state_permutation[start + node * blocks + block] = start + block * model.Nz + node;
        opts.cvode_preconditioner_bandwidth = std::max(opts.cvode_preconditioner_bandwidth, 4 * blocks);
      }
    }
    opts.cvode_krylov_dimension = cfg.value("CVODE_KRYLOV_DIMENSION", 30);
    opts.cvode_maximum_order = cfg.value("CVODE_MAXIMUM_ORDER", 5);
    opts.adaptive_rk=cfg.value("RK_ADAPTIVE",true);
    opts.initialization_steps=cfg.value("RK_INITIALIZATION_STEPS",std::size_t(0));
    opts.event_time_tolerance=cfg.value("EVENT_TIME_TOL",1e-6);
    auto jacobian_mode = cfg.value("SOLVER_JACOBIAN_MODE", std::string("structured_fd"));
    require(jacobian_mode == "dense" || jacobian_mode == "structured_fd",
            "JacobianMode must be dense or structured_fd");
    if ((solver_method == "BDF" || solver_method == "BDF1") && jacobian_mode == "structured_fd" && network.connections.empty() &&
        network.buffers->size == 0)
    {
      opts.jacobian_rows.resize(physical_size);
      const auto radius = cfg.value("ADVECTION_SCHEME", std::string("TVD")) == "WENO" ? 3 : 2;
      for (const auto& name : network.column_names)
      {
        const auto& model = network.model(name);
        auto offset = network.column_starts.at(name);
        auto blocks = model.layout.n_blocks();
        auto local_radius = (model.thermal_mode != 1 && !model.layout.transient_htf) ? model.Nz : std::size_t(radius);
        for (std::size_t block = 0; block < blocks; ++block)
          for (std::size_t node = 0; node < model.Nz; ++node)
          {
            auto& rows = opts.jacobian_rows[offset + block * model.Nz + node];
            for (std::size_t rb = 0; rb < blocks; ++rb)
              for (std::size_t rn = node > local_radius ? node - local_radius : 0;
                   rn < std::min(model.Nz, node + local_radius + 1); ++rn)
                rows.push_back(offset + rb * model.Nz + rn);
            for (std::size_t row = physical_size; row < physical_size + ledger_size; ++row) rows.push_back(row);
          }
      }
    }
    opts.store_accepted_steps = !compact_output;
    if (compact_output)
    {
      double horizon = phase_rule.duration.value_or(phase_rule.maximum.value_or(failure));
      if (!phase_rule.enabled())
      {
        horizon = 0;
        for (const auto& [name, rule] : action_rules)
          if (rule.enabled()) horizon = std::max(horizon, rule.duration.value_or(rule.maximum.value_or(failure)));
      }
      const double spacing = std::max(output_interval, horizon / output_samples);
      for (double k = std::floor(elapsed / spacing) + 1; k * spacing < next_time; ++k)
        opts.sample_times.push_back(k * spacing);
    }
    opts.maximum_steps = size_t(num(cfg, "SOLVER_MAXIMUM_STEPS", 1000000));
    Vec cache_state, cache_derivative;
    double cache_time = 0.;
    bool cache_valid = false;
    const bool use_cache = cfg.value("SOLVER_CACHE_RHS", true);
    const OdeRhs augmented_rhs = [this, &cache_state, &cache_derivative, &cache_time, &cache_valid, use_cache](
                                     double t, const Vec& state, Vec& deriv)
    {
      if (use_cache && cache_valid && t == cache_time &&
          std::equal(state.begin(), state.begin() + physical_size, cache_state.begin()))
      {
        deriv = cache_derivative;
        return;
      }
      deriv = rhs(t, state);
      if (use_cache)
      {
        cache_time = t;
        cache_state.assign(state.begin(), state.begin() + physical_size);
        cache_derivative = deriv;
        cache_valid = true;
      }
    };
    const OdeAcceptedCallback check_accepted = [this](double, const Vec& state, const BdfStatistics&)
    {
      try
      {
        validate_state(slice(state, 0, physical_size));
      }
      catch (const std::exception& e)
      {
        throw std::runtime_error(std::string("Accepted solver step: ") + e.what());
      }
      return true;
    };
    BdfResult sol;
    if(solver_method=="CVODE")sol=integrate_cvode(augmented_rhs,elapsed,augmented,next_time,opts,check_accepted,events);
    else if(solver_method=="BDF1")sol=integrate_bdf1(augmented_rhs,elapsed,augmented,next_time,opts,check_accepted,events);
    else if(solver_method=="SSPRK3")sol=integrate_ssprk3(augmented_rhs,elapsed,augmented,next_time,opts,check_accepted,events);
    else if(solver_method=="SemiImplicitSSPRK3") {
      const OdeRhs explicit_rhs=[&](double t,const Vec& state,Vec& dy){
        augmented_rhs(t,state,dy);
        const auto physical=network.explicit_transport_rhs(t,slice(state,0,physical_size));
        std::copy(physical.begin(),physical.end(),dy.begin());
      };
      sol=integrate_semi_implicit_ssprk3(augmented_rhs,explicit_rhs,elapsed,augmented,next_time,opts,check_accepted,events);
    } else sol=integrate_bdf(augmented_rhs,elapsed,augmented,next_time,opts,check_accepted,events);
    if (!sol.success) throw std::runtime_error(solver_method + " failed: " + sol.message);
    for (const auto& state : sol.states) try
      {
        validate_state(slice(state, 0, physical_size));
      }
      catch (const std::exception& e)
      {
        throw std::runtime_error(std::string("Recorded sample: ") + e.what());
      }
    if (sol.time <= elapsed && sol.events.empty()) throw std::runtime_error("Stage event failed to advance time");
    if (cfg.value("WRITE_MAIN_FORMAT", false)) append_main_frames(sol.times, sol.states);
    elapsed = sol.time;
    augmented = sol.state;
    for (const auto& hit : sol.events)
      done[event_columns.at(hit.index).empty() ? ":phase" : event_columns.at(hit.index)] = {{"reason", "criterion"},
                                                                                            {"elapsed_s", elapsed}};
    Vec times;
    std::vector<Vec> values;
    for (size_t i = 1; i < sol.times.size(); ++i)
    {
      times.push_back(sol.times[i]);
      values.push_back(sol.states[i]);
      if (!compact_output || !commit)
      {
        result.times.push_back(sol.times[i]);
        result.augmented.push_back(sol.states[i]);
      }
    }
    if (compact_output && commit)
    {
      result.times = {elapsed};
      result.augmented.clear();
    }
    result.statistics.accepted_steps += sol.statistics.accepted_steps;
    result.statistics.rejected_steps += sol.statistics.rejected_steps;
    result.statistics.rhs_evaluations += sol.statistics.rhs_evaluations;
    result.statistics.jacobian_evaluations += sol.statistics.jacobian_evaluations;
    result.statistics.lu_factorizations += sol.statistics.lu_factorizations;
    result.statistics.largest_order = std::max(result.statistics.largest_order, sol.statistics.largest_order);
    result.statistics.linear_solver = sol.statistics.linear_solver;
    if (commit) commit(elapsed, augmented, done, times, values);
    bool close = false;
    for (const auto& hit : sol.events) close |= !event_columns.at(hit.index).empty();
    if (close) activate_phase(phase, context, done);
  }
  result.final_augmented = augmented;
  result.ledger = subtract(slice(augmented, physical_size, ledger_size), phase_ledger_start);
  for (const auto& name : network.column_names)
    if (done.contains(name)) result.action_results[name] = done[name];
  return result;
}
Json CycleRunner::fresh_progress()
{
  auto state = network.initial_state();
  state.resize(physical_size + ledger_size, 0);
  Json initial_context = {{"stream_compositions", Json::object()}};
  if (needs_pressure_anchor(phases[0])) initial_context["phase_start_pressures"] = network.endpoint_pressures(state);
  network.set_phase(phases[0], initial_context);
  auto pressure = capture_liquid_pressures(0., state);
  return {{"cycle", 1},
          {"phase", 0},
          {"elapsed", 0.0},
          {"global_time", 0.0},
          {"cycle_start_time", 0.0},
          {"phase_start_time", 0.0},
          {"state", state},
          {"cycle_start", state},
          {"phase_ledger_start", slice(state, physical_size, ledger_size)},
          {"done", Json::object()},
          {"stream_compositions", Json::object()},
          {"previous_compositions", Json::object()},
          {"previous_metrics", Json()},
          {"css_streak", 0},
          {"cycle_rows", Json::array()},
          {"cycle_records", Json::array()},
          {"phase_rows", Json::array()},
          {"stream_rows", Json::array()},
          {"history_pressure", Json::array({pressure})},
          {"history_t", Json::array({0.0})},
          {"history_state", Json::array({slice(state, 0, physical_size)})},
          {"history_cycle", Json::array({1})},
          {"history_phase", Json::array({0})},
          {"completed", false},
          {"summary", Json()}};
}
void CycleRunner::prepare_output()
{
  if (restart)
  {
    require(std::filesystem::is_regular_file(restart_path), "Checkpoint does not exist: " + restart_path.string());
    return;
  }
  if (!cfg.value("ALLOW_EXISTING_OUTPUT_DIR", true) && std::filesystem::exists(output_dir) &&
      !std::filesystem::is_empty(output_dir))
    throw std::runtime_error("Output directory is not empty and ALLOW_EXISTING_OUTPUT_DIR is false: " +
                             output_dir.string());
  std::filesystem::create_directories(output_dir);
  for (const auto& name : {"checkpoint_latest.json",
                           "summary.json",
                           "cycle_metrics.csv",
                           "cycle_summary.csv",
                           "phase_metrics.csv",
                           "step_summary.csv",
                           "streams.csv",
                           "states.json",
                           "states.npz",
                           "metadata.json",
                           "component_balances.csv",
                           "manifest.json",
                           "tea_summary.json",
                           "tea_summary.csv",
                           "convergence_report.json",
                           "selected_cycle.npz",
                           "window_statistics.csv",
                           "tea_handoff.json",
                           "tea_handoff.csv",
                           "history_trends.csv"})
    std::filesystem::remove(output_dir / name);
  auto path = output_dir / "phases";
  if (std::filesystem::is_directory(path))
    for (const auto& entry : std::filesystem::directory_iterator(path))
    {
      auto name = entry.path().filename().string();
      if (entry.is_regular_file() && name.starts_with("cycle_") && name.find("_phase_") != std::string::npos &&
          (entry.path().extension() == ".json" || entry.path().extension() == ".npz"))
        std::filesystem::remove(entry.path());
    }
  if (phase_archives)
    std::filesystem::create_directories(path);
  else if (std::filesystem::is_directory(path) && std::filesystem::is_empty(path))
    std::filesystem::remove(path);
  const auto manifest = output_dir / "plots" / "plot_manifest.json";
  if (std::filesystem::is_regular_file(manifest))
  {
    const auto record = read_json(manifest);
    for (const auto& item : record.value("files", Json::array()))
    {
      const auto file = output_dir / "plots" / std::filesystem::path(item.get<std::string>()).filename();
      if (file.extension() == ".png" || file.extension() == ".pdf" || file.extension() == ".svg")
        std::filesystem::remove(file);
    }
    std::filesystem::remove(manifest);
  }
  std::filesystem::remove(output_dir / "plots" / "plot_log.txt");
}
void CycleRunner::save_checkpoint()
{
  if (compact_output) trim_history();
  if (!checkpoint_enabled) return;
  Json checkpoint = Json::object();
  for (const auto& [key, value] : progress.items())
    if (!compact_output || !key.starts_with("history_")) checkpoint[key] = value;
  if (compact_output)
  {
    // One canonical history archive; checkpoints contain no profile trajectories.
    // Replacing it before the checkpoint is safe: restart discards any uncommitted suffix.
    const auto stamp = std::to_string(progress["history_t"].size()) + ":" + progress["history_t"].back().dump() + ":" +
                       progress["history_cycle"].front().dump() + ":" + progress["history_phase"].back().dump();
    if (stamp != archive_stamp)
    {
      Json arrays = {{"time", progress["history_t"]},
                     {"state", progress["history_state"]},
                     {"cycle", progress["history_cycle"]},
                     {"phase", progress["history_phase"]}};
      if (network.models.front().liquid)
        for (const auto& name : network.column_names)
        {
          arrays[name + "_P"] = Json::array();
          for (const auto& record : progress["history_pressure"]) arrays[name + "_P"].push_back(record.at(name));
        }
      write_npz(output_dir / "states.npz", arrays, compress_output);
      archive_stamp = stamp;
    }
    checkpoint["history_last"] = {progress["history_t"].back(), progress["history_cycle"].back(),
                                  progress["history_phase"].back()};
    const auto trend_path = output_dir / "history_trends.csv";
    checkpoint["trend_bytes"] = std::filesystem::exists(trend_path) ? std::filesystem::file_size(trend_path) : 0;
  }
  if (cfg.value("WRITE_MAIN_FORMAT", false))
  {
    Json bytes = Json::object();
    auto directory = output_dir / "main";
    if (std::filesystem::exists(directory))
      for (const auto& entry : std::filesystem::recursive_directory_iterator(directory))
        if (entry.is_regular_file())
          bytes[entry.path().lexically_relative(directory).generic_string()] = entry.file_size();
    checkpoint["main_file_bytes"] = bytes;
  }
  checkpoint["checkpoint_schema"] = compact_output ? "ruptura_native_cycle_2" : "ruptura_native_cycle_1";
  checkpoint["fingerprint"] = fingerprint;
  checkpoint["solver"] = solver_method;
  checkpoint["solver_version"] = solver_method == "CVODE" ? cvode_version() : "native "+solver_method;
  checkpoint["restart_config"] = restart_config_;
  checkpoint["physical_size"] = physical_size;
  checkpoint["ledger_size"] = ledger_size;
  checkpoint["column_names"] = network.column_names;
  atomic_text(checkpoint_path, checkpoint.dump());
}
Json CycleRunner::load_checkpoint()
{
  auto p = read_json(restart_path);
  auto saved_config = p.value("restart_config", Json());
  if (saved_config.is_object())
  {
    for (const auto* key :
         {"OUTPUT_PROFILE", "OUTPUT_SAMPLE_INTERVAL_S", "OUTPUT_SAMPLES_PER_PHASE", "OUTPUT_HISTORY_CYCLES",
          "CHECKPOINT_ENABLED", "WRITE_PHASE_ARCHIVES", "WRITE_SELECTED_CYCLE", "COMPRESS_NPZ"})
      saved_config.erase(key);
    // Older native checkpoints preserved input casing and could omit the
    // default BDF method. Normalize that representation, not the algorithm.
    require(!saved_config.contains("SOLVER_METHOD") || saved_config["SOLVER_METHOD"].is_string(),
            "Checkpoint SOLVER_METHOD must be a string");
    const auto saved_method = lower(str(saved_config, "SOLVER_METHOD", "BDF"));
    {
      const std::map<std::string,std::string> methods={{"bdf","BDF"},{"bdf1","BDF1"},{"cvode","CVODE"},{"ssprk3","SSPRK3"},{"semiimplicitssprk3","SemiImplicitSSPRK3"}};
      if(methods.contains(saved_method))saved_config["SOLVER_METHOD"]=methods.at(saved_method);
    }
    saved_config = sorted_json(saved_config);
  }
  if (p.contains("solver"))
    require(p["solver"].is_string() && lower(p["solver"].get<std::string>()) == lower(solver_method),
            "Checkpoint solver metadata disagrees with the selected solver");
  if (p.contains("summary") && p["summary"].is_object() && p["summary"].contains("solver"))
    require(
        p["summary"]["solver"].is_string() && lower(p["summary"]["solver"].get<std::string>()) == lower(solver_method),
        "Checkpoint summary solver disagrees with the selected solver");
  require((str(p, "checkpoint_schema") == "ruptura_native_cycle_1" ||
           str(p, "checkpoint_schema") == "ruptura_native_cycle_2") &&
              saved_config == restart_config_ && p.value("physical_size", size_t(0)) == physical_size &&
              p.value("ledger_size", size_t(0)) == ledger_size &&
              p.value("column_names", std::vector<std::string>()) == network.column_names,
          "Checkpoint incompatible with materials, geometry, schedule, solver or components. "
          "Native JSON checkpoints are distinct from Python NPZ checkpoints.");
  for (const auto& k : {"state", "cycle_start"})
  {
    auto state = vec(p.at(k));
    require(state.size() == physical_size + ledger_size, "Checkpoint state dimension incompatible");
    validate_state(slice(state, 0, physical_size));
    for (double x : state) require(std::isfinite(x), "Checkpoint contains nonfinite ledger");
  }
  require(p["phase_ledger_start"].size() == ledger_size, "Checkpoint ledger start dimension incompatible");
  if (!p.contains("history_state"))
  {
    restore_history(p, restart_path.parent_path() / "states.npz", true);
    auto trend_source = restart_path.parent_path() / "history_trends.csv";
    const auto bytes = p.value("trend_bytes", std::uintmax_t(0));
    if (std::filesystem::exists(trend_source))
    {
      require(std::filesystem::file_size(trend_source) >= bytes, "Truncated history_trends.csv");
      std::filesystem::create_directories(output_dir);
      auto target = output_dir / "history_trends.csv";
      if (std::filesystem::absolute(trend_source).lexically_normal() != target.lexically_normal())
        std::filesystem::copy_file(trend_source, target, std::filesystem::copy_options::overwrite_existing);
      std::filesystem::resize_file(target, bytes);
    }
  }
  if (cfg.value("WRITE_MAIN_FORMAT", false) && p.contains("main_file_bytes"))
    for (auto& [relative, size] : p["main_file_bytes"].items())
    {
      auto raw = std::filesystem::path(relative);
      require(!raw.is_absolute() && relative.find("..") == std::string::npos, "Invalid checkpoint output path");
      auto file = output_dir / "main" / raw;
      require(std::filesystem::exists(file) && std::filesystem::file_size(file) >= size.get<std::uintmax_t>(),
              "Missing or truncated main-format output during restart");
      std::filesystem::resize_file(file, size.get<std::uintmax_t>());
    }
  for (const auto& state : p["history_state"]) validate_state(vec(state));
  require(p["history_t"].size() == p["history_state"].size() &&
              p["history_phase"].size() == p["history_state"].size() &&
              p["history_cycle"].size() == p["history_state"].size(),
          "Checkpoint history dimensions incompatible");
  require(num(p, "cycle") >= 1 && num(p, "phase") >= 0 && num(p, "phase") <= phases.size() && num(p, "elapsed") >= 0,
          "Checkpoint schedule position invalid");
  return p;
}
Json CycleRunner::validate()
{
  auto physical = network.initial_state();
  validate_state(physical);
  for (const auto& phase : phases)
  {
    Json context = {{"stream_compositions", Json::object()}};
    if (needs_pressure_anchor(phase)) context["phase_start_pressures"] = network.endpoint_pressures(physical);
    network.set_phase(phase, context);
    auto [derivative, diagnostics] = network.evaluate(0, physical);
    require(derivative.size() == physical_size, "Invalid initial derivative dimension");
    for (double d : derivative)
      require(std::isfinite(d), "Nonfinite initial derivative in phase " + str(phase, "name"));
    ledger_rates(diagnostics);
  }
  Json names = Json::array();
  for (const auto& phase : phases) names.push_back(phase["name"]);
  return {{"valid", true},
          {"columns", network.column_names},
          {"phases", names},
          {"solver", solver_method},
          {"solver_version", solver_method == "CVODE" ? cvode_version() : "native "+solver_method},
          {"physical_size", physical_size},
          {"ledger_size", ledger_size}};
}
Json CycleRunner::run()
{
  validate();
  prepare_output();
  if (!restart && cfg.value("WRITE_MAIN_FORMAT", false)) std::filesystem::remove_all(output_dir / "main");
  progress = restart ? load_checkpoint() : fresh_progress();
  auto& p = progress;
  if (compact_output && !restart) append_trends({0.}, {vec(p["state"])}, 1, 0);
  if (!p.contains("cycle_records")) p["cycle_records"] = Json::array();
  if (p.value("completed", false))
  {
    bool extend =
        cycle_limit > p["cycle"].get<int>() && (run_mode == "fixed_cycles" || !p["summary"].value("converged", false));
    if (!extend)
    {
      // Older native BDF checkpoints predate explicit solver metadata.
      if (!p["summary"].contains("solver")) p["summary"]["solver"] = solver_method;
      if (!p["summary"].contains("solver_version"))
        p["summary"]["solver_version"] = solver_method == "CVODE" ? cvode_version() : "native "+solver_method;
      update_convergence_report();
      save_checkpoint();
      write_outputs(true);
      return p["summary"];
    }
    p["completed"] = false;
    p["summary"] = Json();
    p["cycle"] = p["cycle"].get<int>() + 1;
    p["phase"] = 0;
    p["elapsed"] = 0.0;
    p["done"] = Json::object();
    p["cycle_start"] = p["state"];
    p["cycle_start_time"] = p["global_time"];
    p["phase_start_time"] = p["global_time"];
    p["phase_ledger_start"] = slice(vec(p["state"]), physical_size, ledger_size);
  }
  save_checkpoint();
  while (p["cycle"].get<int>() <= cycle_limit)
  {
    while (p["phase"].get<size_t>() < phases.size())
    {
      int cycle = p["cycle"].get<int>(), index = p["phase"].get<int>();
      Json phase = phases[index];
      if (phase.contains("duration_fraction_of_previous"))
      {
        require(!p["phase_rows"].empty(), "Post-saturation observation requires a completed stopping stage");
        phase["duration"] = num(p["phase_rows"].back(), "duration_s") * num(phase, "duration_fraction_of_previous");
        phases[index]["duration"] = phase["duration"];
      }
      Json context = {{"cycle_index", cycle},
                      {"phase_index", index},
                      {"stream_compositions", p["stream_compositions"]},
                      {"global_time_offset", p["phase_start_time"]}};
      bool switching_record = false;
      if (num(p, "elapsed") == 0 &&
          (p["history_phase"].back() != p["phase"] || p["history_cycle"].back() != p["cycle"]))
      {
        switching_record = true;
        p["history_t"].push_back(p["global_time"]);
        p["history_state"].push_back(slice(vec(p["state"]), 0, physical_size));
        p["history_cycle"].push_back(cycle);
        p["history_phase"].push_back(index);
      }
      if (needs_pressure_anchor(phase))
      {
        bool found = p.contains("phase_start_pressures") && num(p, "elapsed") > 0;
        if (found) context["phase_start_pressures"] = p["phase_start_pressures"];
        for (size_t i = 0; i < p["history_t"].size(); ++i)
          if (!found && p["history_cycle"][i] == cycle && p["history_phase"][i] == index)
          {
            context["phase_start_pressures"] = network.endpoint_pressures(vec(p["history_state"][i]));
            found = true;
            break;
          }
        require(found, "Checkpoint lacks original switching state for anchored ramp");
        p["phase_start_pressures"] = context["phase_start_pressures"];
      }
      if (switching_record)
      {
        activate_phase(phase, context, p["done"]);
        p["history_pressure"].push_back(capture_liquid_pressures(0., vec(p["state"])));
        if (compact_output) append_trends({0.}, {vec(p["state"])}, cycle, index);
      }
      auto commit =
          [&](double elapsed, const Vec& augmented, const Json& done, const Vec& times, const std::vector<Vec>& states)
      {
        p["elapsed"] = elapsed;
        p["state"] = augmented;
        p["done"] = done;
        p["global_time"] = num(p, "phase_start_time") + elapsed;
        if (compact_output) append_trends(times, states, cycle, index);
        for (size_t i = 0; i < times.size(); ++i)
        {
          p["history_t"].push_back(num(p, "phase_start_time") + times[i]);
          p["history_state"].push_back(slice(states[i], 0, physical_size));
          p["history_pressure"].push_back(capture_liquid_pressures(times[i], states[i]));
          p["history_cycle"].push_back(p["cycle"]);
          p["history_phase"].push_back(p["phase"]);
        }
        save_checkpoint();
      };
      auto state = vec(p["state"]);
      auto solution = integrate_phase(slice(state, 0, physical_size), phase, context, state, num(p, "elapsed"),
                                      p["done"], vec(p["phase_ledger_start"]), commit);
      p["state"] = solution.final_augmented;
      commit_phase(phase, solution);
    }
    auto start = vec(p["cycle_start"]), end = vec(p["state"]),
         delta = subtract(slice(end, physical_size, ledger_size), slice(start, physical_size, ledger_size));
    auto metrics = cycle_metrics(slice(start, 0, physical_size), slice(end, 0, physical_size), delta,
                                 num(p, "global_time") - num(p, "cycle_start_time"));
    auto [css_record, reached] = css(slice(start, 0, physical_size), slice(end, 0, physical_size), metrics,
                                     p["previous_metrics"], p["stream_compositions"], p["previous_compositions"]);
    p["css_streak"] = reached ? p["css_streak"].get<int>() + 1 : 0;
    bool converged = p["css_streak"].get<int>() >= required_streak;
    Json row = {{"cycle", p["cycle"]}};
    row.update(metrics);
    row.update(css_record);
    row["css_period_errors"] = Json::object();
    for (int lag = 2; lag <= cfg.value("CSS_MAXIMUM_PERIOD", 4); ++lag)
    {
      const auto& records = p["cycle_records"];
      if (records.size() < std::size_t(lag)) continue;
      const auto& earlier = records[records.size() - lag];
      auto [record, passed] = css(earlier["end"].get<Vec>(), slice(end, 0, physical_size), metrics, earlier["metrics"],
                                  p["stream_compositions"], earlier["compositions"]);
      record["criteria_satisfied"] = passed;
      row["css_period_errors"][std::to_string(lag)] = record;
    }
    row["css_streak"] = p["css_streak"];
    row["converged"] = converged;
    p["cycle_rows"].push_back(row);
    p["cycle_records"].push_back({{"cycle", p["cycle"]},
                                  {"start", slice(start, 0, physical_size)},
                                  {"end", slice(end, 0, physical_size)},
                                  {"ledger", delta},
                                  {"metrics", metrics},
                                  {"compositions", p["stream_compositions"]}});
    const auto keep =
        std::size_t(std::max(cfg.value("CSS_DIAGNOSTIC_WINDOW", 20), cfg.value("CSS_MAXIMUM_PERIOD", 4)) + 1);
    while (p["cycle_records"].size() > keep) p["cycle_records"].erase(p["cycle_records"].begin());
    p["previous_metrics"] = metrics;
    p["previous_compositions"] = p["stream_compositions"];
    if (cfg.value("SIMULATION_TYPE", std::string()) != "Breakthrough")
      std::cout << "Cycle " << p["cycle"] << ": CSS state error " << css_record["css_state_error"] << "; streak "
                << p["css_streak"] << "/" << required_streak << std::endl;
    bool finished = (run_mode == "css" && converged) || p["cycle"].get<int>() >= cycle_limit;
    if (finished)
    {
      std::string reason = run_mode == "css" && converged ? "css_reached"
                           : run_mode == "fixed_cycles"   ? "fixed_cycles_completed"
                                                          : "maximum_cycles_reached";
      p["summary"] = {{"completed", run_mode == "fixed_cycles" || converged},
                      {"converged", converged},
                      {"cycles", p["cycle"]},
                      {"run_mode", run_mode},
                      {"termination_reason", reason},
                      {"output_dir", output_dir.string()},
                      {"final_metrics", metrics},
                      {"metrics", metrics},
                      {"css", css_record},
                      {"stream_moles", stream_record(delta)},
                      {"stream_compositions", p["stream_compositions"]},
                      {"number_of_columns", network.models.size()},
                      {"column_names", network.column_names}};
      if (cfg.value("SIMULATION_TYPE", std::string()) == "Breakthrough")
      {
        p["summary"]["termination_reason"] =
            phases.size() > 1 ? "stopping_and_observation_completed" : "breakthrough_stage_completed";
        p["summary"]["css_applicable"] = false;
      }
      p["summary"]["solver"] = solver_method;
      p["summary"]["solver_version"] = solver_method == "CVODE" ? cvode_version() : "native "+solver_method;
      p["completed"] = true;
      update_convergence_report();
      save_checkpoint();
      write_outputs(true);
      return p["summary"];
    }
    p["cycle"] = p["cycle"].get<int>() + 1;
    p["phase"] = 0;
    p["cycle_start"] = p["state"];
    p["cycle_start_time"] = p["global_time"];
    p["phase_start_time"] = p["global_time"];
    save_checkpoint();
    write_outputs();
  }
  throw std::runtime_error("Invalid checkpoint cycle position");
}
}  // namespace cycle

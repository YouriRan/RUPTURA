#include "cycle/input.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <regex>
#include <set>
#include <stdexcept>

#include "yaml.h"

namespace cycle {
namespace {
constexpr double R = 8.31446261815324, PI = 3.141592653589793238462643383279502884;
void require(bool value, const std::string& message) {
  if (!value) throw std::invalid_argument(message);
}
void object(const Json& value, const std::string& path) {
  require(value.is_object(), path + " must be a mapping");
}
void allowed(const Json& value, const std::set<std::string>& keys, const std::string& path) {
  object(value, path);
  for (auto i = value.begin(); i != value.end(); ++i)
    require(keys.count(i.key()), "Unknown " + path + " option: " + i.key());
}
std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char x) { return std::tolower(x); });
  return s;
}
std::string upper(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char x) { return std::toupper(x); });
  return s;
}
double number(const Json& v) {
  if (v.is_number()) return v.get<double>();
  if (v.is_string()) {
    auto s = v.get<std::string>();
    std::size_t n = 0;
    double x = std::stod(s, &n);
    require(n == s.size(), "Expected numeric value");
    return x;
  }
  throw std::invalid_argument("Expected numeric value: " + v.dump());
}
void update(Json& target, const Json& source) {
  object(source, "mapping");
  for (auto i = source.begin(); i != source.end(); ++i) target[i.key()] = i.value();
}
const Json& maps() {
  static const Json x = Json::parse(
      R"MAPS({"simulation": {"name": "CASE_NAME", "description": "DESCRIPTION", "run_mode": "RUN_MODE", "cycles": "N_CYCLES", "maximum_cycles": "MAX_CYCLES"}, "gas": {"viscosity_Pa_s": "mu", "heat_capacity_J_kg_K": "Cpg", "thermal_conductivity_W_m_K": "lambda_g"}, "bed": {"type": "CONTACTOR", "length_m": "L", "nodes": "Nz", "inner_diameter_m": "dint", "outer_diameter_m": "dout", "area_m2": "COLUMN_AREA", "void_fraction": "eps_b", "gas_storage_void_fraction": "GAS_STORAGE_VOID_FRACTION", "adsorbent_density_kg_m3": "rho_s", "particle_diameter_m": "dp", "solid_heat_capacity_J_kg_K": "Cps", "gas_solid_heat_transfer_W_m2_K": "hgs"}, "thermal": {"wall_density_kg_m3": "rho_w", "wall_heat_capacity_J_kg_K": "Cpw", "wall_conductivity_W_m_K": "lambda_w", "gas_wall_heat_transfer_W_m2_K": "hgw", "ambient_heat_transfer_W_m2_K": "hw_inf", "ambient_temperature_K": "T_inf"}, "kinetics": {"ldf_rate_s_inv": "km_phy", "adsorption_heat_J_mol": "dH_phy", "ldf_activation_energy_J_mol": "LDF_ACTIVATION_ENERGY_J_MOL", "reference_temperature_K": "LDF_REFERENCE_TEMPERATURE", "physisorption_enabled": "PHYSISORPTION_ENABLED", "chemisorption": "CHEMISORPTION"}, "solver": {"method": "SOLVER_METHOD", "rtol": "SOLVER_RTOL", "atol": "SOLVER_ATOL", "maximum_step_s": "SOLVER_MAX_STEP", "reconstruction": "ADVECTION_SCHEME", "jacobian_mode": "SOLVER_JACOBIAN_MODE", "cache_rhs": "SOLVER_CACHE_RHS", "event_time_tolerance_s": "EVENT_TIME_TOL", "checkpoint_interval_s": "CHECKPOINT_INTERVAL_S", "dispersion_closure": "DISPERSION_CLOSURE"}, "output": {"profile": "OUTPUT_PROFILE", "sample_interval_s": "OUTPUT_SAMPLE_INTERVAL_S", "samples_per_phase": "OUTPUT_SAMPLES_PER_PHASE", "history_cycles": "OUTPUT_HISTORY_CYCLES", "checkpoint": "CHECKPOINT_ENABLED", "phase_archives": "WRITE_PHASE_ARCHIVES", "selected_cycle_archive": "WRITE_SELECTED_CYCLE", "compress_npz": "COMPRESS_NPZ", "directory": "OUTPUT_DIR", "plots": "MAKE_PLOTS", "reference_pressure_Pa": "PLOT_REFERENCE_PRESSURE", "reference_temperature_K": "PLOT_REFERENCE_TEMPERATURE"}, "css": {"state_rtol": "CSS_STATE_RTOL", "state_atol": "CSS_STATE_ATOL", "inventory_rtol": "CSS_INVENTORY_RTOL", "metric_rtol": "CSS_METRIC_RTOL", "consecutive_cycles": "CSS_CONSECUTIVE_CYCLES"}, "performance": {"feed_tag": "FEED_TAG", "light_product_tag": "LIGHT_PRODUCT_TAG", "heavy_product_tag": "HEAVY_PRODUCT_TAG", "light_component": "TARGET_LIGHT_COMPONENT", "heavy_component": "TARGET_HEAVY_COMPONENT"}, "utilities": {"ambient_pressure_Pa": "AMBIENT_PRESSURE", "feed_supply_pressure_Pa": "FEED_SUPPLY_PRESSURE", "product_delivery_pressure_Pa": "PRODUCT_DELIVERY_PRESSURE", "product_buffer_pressure_Pa": "PRODUCT_BUFFER_PRESSURE", "plant_fresh_feed_target_Nm3_h": "PLANT_FRESH_FEED_TARGET_NM3_H", "charge_feed_compression": "CHARGE_FEED_COMPRESSION", "charge_recycle_compression": "CHARGE_RECYCLE_COMPRESSION", "compressor_efficiency": "COMPRESSOR_ISENTROPIC_EFFICIENCY", "vacuum_efficiency": "VACUUM_PUMP_EFFICIENCY", "compression_model": "COMPRESSION_WORK_MODEL", "heat_capacity_ratio": "COMPRESSION_HEAT_CAPACITY_RATIO", "heating_efficiency": "HEATING_UTILITY_EFFICIENCY", "cooling_COP": "COOLING_COP", "heat_recovery_fraction": "HEAT_RECOVERY_FRACTION", "electricity_cost_per_kWh": "ELECTRICITY_COST_PER_KWH", "heat_cost_per_kWh": "HEAT_UTILITY_COST_PER_KWH_TH", "cooling_cost_per_kWh": "COOLING_UTILITY_COST_PER_KWH_TH", "annual_operating_hours": "ANNUAL_OPERATING_HOURS", "availability": "PROCESS_AVAILABILITY", "standard_temperature_K": "STANDARD_TEMPERATURE_K", "standard_pressure_Pa": "STANDARD_PRESSURE_PA", "cooling_cost_basis": "COOLING_COST_BASIS"}})MAPS");
  return x;
}
const Json& port_aliases() {
  static const Json x = Json::parse(
      R"ALIASES({"pressure_Pa": "pressure", "temperature_K": "temperature", "velocity_m_s": "velocity", "duration_s": "duration", "minimum_duration_s": "minimum_duration", "maximum_duration_s": "maximum_duration", "failure_maximum_simulated_time_s": "failure_maximum_simulated_time", "initial_pressure_Pa": "initial_pressure", "initial_temperature_K": "initial_temperature", "outlet_pressure_Pa": "outlet_pressure", "conductance_mol_s_Pa": "conductance"})ALIASES");
  return x;
}
Json aliases(const Json& v) {
  if (v.is_object()) {
    Json out = Json::object();
    for (auto i = v.begin(); i != v.end(); ++i) {
      auto key = port_aliases().value(i.key(), i.key());
      require(!out.contains(key), "Duplicate aliases for " + key);
      out[key] = aliases(i.value());
    }
    return out;
  }
  if (v.is_array()) {
    Json out = Json::array();
    for (const auto& x : v) out.push_back(aliases(x));
    return out;
  }
  return v;
}
void mapped(const Json& data, const std::string& section, Json& target) {
  object(data, section);
  const Json diagnostics = {{"diagnostic_window", "CSS_DIAGNOSTIC_WINDOW"},
      {"diagnostic_minimum_cycles", "CSS_DIAGNOSTIC_MINIMUM_CYCLES"},
      {"maximum_period", "CSS_MAXIMUM_PERIOD"},
      {"plot_cycle_selection", "CSS_PLOT_CYCLE_SELECTION"},
      {"mass_balance_tolerance", "CSS_MASS_BALANCE_TOLERANCE"}};
  for (auto i = data.begin(); i != data.end(); ++i) {
    if (section == "css" && diagnostics.contains(i.key())) {
      target[diagnostics[i.key()].get<std::string>()] = i.value();
      continue;
    }
    require(maps().at(section).contains(i.key()), "Unknown " + section + " option: " + i.key());
    target[maps().at(section).at(i.key()).get<std::string>()] = i.value();
  }
}
Json composition(const Json& v, const std::vector<std::string>& names, const std::string& path) {
  Json out = Json::array();
  if (v.is_object()) {
    require(v.size() == names.size(), path + " requires every component");
    for (const auto& n : names) {
      require(v.contains(n), path + " is missing " + n);
      out.push_back(number(v.at(n)));
    }
  } else {
    require(v.is_array() && v.size() == names.size(),
            path + " requires one fraction per component");
    for (const auto& x : v) out.push_back(number(x));
  }
  double sum = 0;
  for (const auto& x : out) {
    double a = x.get<double>();
    require(std::isfinite(a) && a >= 0, path + " fractions must be finite and nonnegative");
    sum += a;
  }
  require(std::abs(sum - 1.) <= 1e-8, path + " fractions must sum to 1");
  return out;
}
Json compositions(const Json& v, const std::vector<std::string>& names,
                  const std::string& path = "input") {
  if (v.is_object()) {
    Json out = Json::object();
    const std::set<std::string> ck = {"composition", "initial_composition", "fallback_composition",
                                      "y_init", "y_feed"};
    for (auto i = v.begin(); i != v.end(); ++i)
      out[i.key()] = ck.count(i.key()) ? composition(i.value(), names, path + "." + i.key())
                                       : compositions(i.value(), names, path + "." + i.key());
    return out;
  }
  if (v.is_array()) {
    Json out = Json::array();
    for (const auto& x : v) out.push_back(compositions(x, names, path));
    return out;
  }
  return v;
}
Json capacity_temperature(const Json& data) {
  if (data.is_null()) return nullptr;
  allowed(data, {"model", "reference_temperature_K", "coefficient"}, "capacity_temperature");
  auto model = lower(data.value("model", std::string("exponential")));
  require(model == "exponential", "capacity_temperature.model must be exponential");
  double tref = number(data.at("reference_temperature_K")), a = number(data.at("coefficient"));
  require(std::isfinite(tref) && tref > 0 && std::isfinite(a),
          "Invalid capacity temperature parameters");
  return Json{{"model", model}, {"reference_temperature_K", tref}, {"coefficient", a}};
}
void equilibrium(const Json& data, const std::vector<std::string>& names, Json& target,
                 const std::string& label) {
  allowed(data, {"method", "reference_temperature_K", "components", "name"}, "equilibrium");
  auto method = upper(data.value("method", std::string("IAST")));
  require(std::set<std::string>{"SPI", "SCI", "SEI", "IAST", "SIAST"}.count(method),
          "Unsupported equilibrium.method");
  double tref = number(data.value("reference_temperature_K", Json(298.15)));
  require(std::isfinite(tref) && tref > 0, "equilibrium.reference_temperature_K must be positive");
  const auto& comps = data.at("components");
  object(comps, "equilibrium.components");
  require(comps.size() == names.size(), "equilibrium.components must specify every component");
  Json parsed = Json::object(), capacities = Json::array(), prefactors = Json::array(),
       heats = Json::array();
  std::set<std::size_t> counts;
  std::set<std::string> models;
  std::size_t maxsites = 0;
  for (const auto& name : names) {
    require(comps.contains(name), "equilibrium missing " + name);
    const auto& comp = comps.at(name);
    allowed(comp, {"model", "sites"}, "equilibrium.components." + name);
    Json sites = Json::array();
    double capacity = 0;
    for (const auto& site : comp.at("sites")) {
      object(site, name + ".sites");
      Json raw = site;
      require(raw.contains("capacity_mol_kg") || raw.contains("qsat"),
              name + ": site capacity required");
      double q = number(raw.value("capacity_mol_kg", raw.value("qsat", Json())));
      raw.erase("capacity_mol_kg");
      raw.erase("qsat");
      double heat = number(raw.value("heat_J_mol", raw.value("dH", Json(0.))));
      raw.erase("heat_J_mol");
      raw.erase("dH");
      double k0 = 0.;
      if (raw.contains("affinity_Pa_inv")) {
        require(raw.value("affinity_basis", std::string("pressure")) == "pressure",
                name + ": concentration affinity requires k0 in m3/mol, not affinity_Pa_inv");
        require(!raw.contains("k0"), name + ": use affinity_Pa_inv or k0, not both");
        k0 = number(raw.at("affinity_Pa_inv")) * std::exp(-heat / (R * tref));
        raw.erase("affinity_Pa_inv");
      } else {
        require(raw.contains("k0"), name + ": affinity is required");
        k0 = number(raw.at("k0"));
        raw.erase("k0");
      }
      require(std::isfinite(q) && std::isfinite(k0) && std::isfinite(heat) && q >= 0 && k0 >= 0,
              "Invalid site parameters for " + name);
      if (raw.contains("capacity_temperature"))
        raw["capacity_temperature"] = capacity_temperature(raw.at("capacity_temperature"));
      Json s = {{"qsat", q}, {"k0", k0}, {"dH", heat}};
      update(s, raw);
      sites.push_back(s);
      capacity += q;
    }
    require(!sites.empty(), name + ": at least one site is required");
    auto model = comp.value("model", std::string("langmuir"));
    parsed[name] = {{"model", model}, {"sites", sites}};
    capacities.push_back(capacity);
    prefactors.push_back(sites[0]["k0"]);
    heats.push_back(sites[0]["dH"]);
    counts.insert(sites.size());
    models.insert(model);
    maxsites = std::max(maxsites, sites.size());
  }
  Json record = {{"adsorbent", data.value("name", label)}, {"mixture_mode", method}};
  if (method == "SCI") {
    require(counts.size() == 1 && models.size() == 1, "SCI requires matching model and site count");
    record["model"] = *models.begin();
    record["n_sites"] = *counts.begin();
    record["sites"] = Json::array();
    for (std::size_t k = 0; k < *counts.begin(); ++k) {
      Json site = Json::object();
      for (const auto& n : names) site[n] = parsed[n]["sites"][k];
      record["sites"].push_back(site);
    }
  } else {
    record["n_sites"] = maxsites;
    record["components"] = parsed;
  }
  update(target, Json{{"MIXTURE_MODE", method},
                      {"PHY_MIXTURE_MODE", method},
                      {"PHY_ISOTHERM_INPUT", Json::array({record})},
                      {"qmax_phy", Json::array({capacities})},
                      {"kj0_phy", Json::array({prefactors})},
                      {"dH_phy", Json::array({heats})},
                      {"ISOTHERM_REFERENCE_TEMPERATURE", tref}});
}
Json column_sections(const Json& data, const std::vector<std::string>& names,
                     const std::string& label = "column") {
  allowed(data,
          {"bed", "thermal", "kinetics", "equilibrium", "initial", "monolith", "layers", "material",
           "advanced"},
          label);
  Json out = Json::object();
  if (data.contains("equilibrium")) equilibrium(data.at("equilibrium"), names, out, label);
  for (const std::string s : {"bed", "kinetics"})
    if (data.contains(s)) mapped(data.at(s), s, out);
  if (data.contains("thermal")) {
    Json t = data.at("thermal");
    if (t.contains("htf")) {
      out["HTF"] = t.at("htf");
      t.erase("htf");
    }
    if (t.contains("mode")) {
      Json modes = {{"isothermal", 1}, {"adiabatic", 2}, {"wall_exchange", 3}};
      auto m = t.at("mode").get<std::string>();
      require(modes.contains(m), "Invalid thermal.mode: " + m);
      out["THERMAL_MODE"] = modes.at(m);
      t.erase("mode");
    }
    mapped(t, "thermal", out);
  }
  if (data.contains("initial")) {
    Json ini = aliases(data.at("initial")), x = Json::object();
    Json rename = {{"gas_temperature_K", "gas_temperature"},
                   {"solid_temperature_K", "solid_temperature"},
                   {"wall_temperature_K", "wall_temperature"},
                   {"jacket_temperature_K", "jacket_temperature"}};
    for (auto i = ini.begin(); i != ini.end(); ++i) x[rename.value(i.key(), i.key())] = i.value();
    out["INITIAL_STATE"] = x;
    for (const auto& p :
         std::vector<std::pair<std::string, std::string>>{{"pressure", "P_init"},
                                                          {"temperature", "T_init"},
                                                          {"gas_temperature", "T_init"},
                                                          {"composition", "y_init"},
                                                          {"loading", "INITIAL_LOADING_MODE"}})
      if (x.contains(p.first)) out[p.second] = x.at(p.first);
  }
  for (const auto& p : std::vector<std::pair<std::string, std::string>>{
           {"monolith", "MONOLITH"}, {"layers", "LAYERS"}, {"material", "MATERIAL"}})
    if (data.contains(p.first)) out[p.second] = data.at(p.first);
  if (data.contains("advanced")) update(out, data.at("advanced"));
  if (out.contains("dint") && !out.contains("COLUMN_AREA")) {
    double d = number(out.at("dint"));
    out["COLUMN_AREA"] = PI * d * d / 4.;
  }
  return out;
}
Json yaml_scalar(const yaml_node_t& n) {
  std::string s(reinterpret_cast<char*>(n.data.scalar.value), n.data.scalar.length);
  if (n.data.scalar.style != YAML_PLAIN_SCALAR_STYLE) return s;
  auto l = lower(s);
  if (l == "null" || l == "~" || s.empty()) return nullptr;
  if (l == "true" || l == "yes" || l == "on") return true;
  if (l == "false" || l == "no" || l == "off") return false;
  std::string a = s;
  a.erase(std::remove(a.begin(), a.end(), '_'), a.end());
  static const std::regex integer(R"([+-]?(?:0|[1-9][0-9]*|0[0-7]+|0x[0-9a-fA-F]+))"),
      floating(R"([+-]?(?:[0-9]+\.[0-9]*|\.[0-9]+)(?:[eE][+-][0-9]+)?)");
  try {
    if (std::regex_match(a, integer)) return std::stoll(a, nullptr, 0);
    if (std::regex_match(a, floating)) return std::stod(a);
  } catch (const std::exception&) {
    throw std::invalid_argument("YAML number out of range: " + s);
  }
  if (l == ".inf" || l == "+.inf" || l == "-.inf" || l == ".nan")
    throw std::invalid_argument("Nonfinite YAML numbers are not accepted");
  return s;
}
Json yaml_node(yaml_document_t& doc, int id, std::set<int>& visiting) {
  require(visiting.insert(id).second, "Recursive YAML alias is not supported");
  auto* n = yaml_document_get_node(&doc, id);
  require(n, "Malformed YAML node");
  Json out;
  if (n->type == YAML_SCALAR_NODE)
    out = yaml_scalar(*n);
  else if (n->type == YAML_SEQUENCE_NODE) {
    out = Json::array();
    for (auto p = n->data.sequence.items.start; p < n->data.sequence.items.top; ++p)
      out.push_back(yaml_node(doc, *p, visiting));
  } else if (n->type == YAML_MAPPING_NODE) {
    out = Json::object();
    for (auto p = n->data.mapping.pairs.start; p < n->data.mapping.pairs.top; ++p) {
      auto* key = yaml_document_get_node(&doc, p->key);
      require(key && key->type == YAML_SCALAR_NODE, "YAML mapping keys must be scalar strings");
      std::string name(reinterpret_cast<char*>(key->data.scalar.value), key->data.scalar.length);
      require(!out.contains(name), "Duplicate YAML key: " + name);
      require(name != "<<", "YAML merge keys are not supported; use explicit fields");
      out[name] = yaml_node(doc, p->value, visiting);
    }
  } else
    throw std::invalid_argument("Unsupported YAML node");
  visiting.erase(id);
  return out;
}
}  // namespace
Json read_input_document(const std::filesystem::path& path) {
  std::ifstream f(path, std::ios::binary);
  require(bool(f), "Cannot read input: " + path.string());
  std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  if (lower(path.extension().string()) == ".json") {
    std::vector<std::set<std::string>> keys;
    return Json::parse(text,[&](int,Json::parse_event_t event,Json& value){
      if(event==Json::parse_event_t::object_start)keys.emplace_back();
      else if(event==Json::parse_event_t::key)require(keys.back().insert(value.get<std::string>()).second,"Duplicate JSON key: "+value.get<std::string>());
      else if(event==Json::parse_event_t::object_end)keys.pop_back();
      return true;
    });
  }
  require(lower(path.extension().string()) == ".yaml" || lower(path.extension().string()) == ".yml",
          "Native input must be YAML or JSON; Python cfg.py is not executable input");
  yaml_parser_t parser;
  require(yaml_parser_initialize(&parser), "Cannot initialize YAML parser");
  yaml_parser_set_input_string(&parser, reinterpret_cast<const unsigned char*>(text.data()),
                               text.size());
  yaml_document_t doc;
  if (!yaml_parser_load(&parser, &doc)) {
    auto message = "YAML parse error at line " + std::to_string(parser.problem_mark.line + 1) +
                   ": " + (parser.problem ? parser.problem : "unknown error");
    yaml_parser_delete(&parser);
    throw std::invalid_argument(message);
  }
  Json result;
  try {
    auto* root = yaml_document_get_root_node(&doc);
    require(root, "Empty YAML document");
    std::set<int> visiting;
    result = yaml_node(doc, int(root - doc.nodes.start + 1), visiting);
    yaml_document_delete(&doc);
  } catch (...) {
    yaml_document_delete(&doc);
    yaml_parser_delete(&parser);
    throw;
  }
  yaml_document_t extra;
  if (!yaml_parser_load(&parser, &extra)) {
    yaml_parser_delete(&parser);
    throw std::invalid_argument("Invalid trailing YAML document");
  }
  bool multiple = yaml_document_get_root_node(&extra) != nullptr;
  yaml_document_delete(&extra);
  yaml_parser_delete(&parser);
  require(!multiple, "Only one YAML document is allowed");
  return result;
}
Json normalize_config(const Json& data, const std::filesystem::path& base_dir) {
  std::set<std::string> sections = {
      "schema_version", "components", "feed",    "initial", "equilibrium", "monolith", "materials",
      "layers",         "columns",    "sources", "buffers", "cycle",       "advanced", "notes"};
  for (auto i = maps().begin(); i != maps().end(); ++i) sections.insert(i.key());
  allowed(data, sections, "input section");
  require(number(data.value("schema_version", Json(1))) == 1,
          "Only schema_version: 1 is supported");
  const auto& components = data.at("components");
  require(components.is_array() && !components.empty(),
          "components must be a nonempty ordered list");
  std::vector<std::string> names;
  Json masses = Json::array(), dispersion = Json::array();
  bool anycp = false, allcp = true, anyads = false;
  for (const auto& comp : components) {
    allowed(comp,
            {"name", "molar_mass_kg_mol", "dispersion_m2_s", "gas_heat_capacity_J_mol_K",
             "adsorbed_heat_capacity_J_mol_K"},
            "component");
    auto name = comp.at("name").get<std::string>();
    require(std::find(names.begin(), names.end(), name) == names.end(),
            "Duplicate component name: " + name);
    names.push_back(name);
    masses.push_back(number(comp.at("molar_mass_kg_mol")));
    dispersion.push_back(number(comp.value("dispersion_m2_s", Json(0.))));
    anycp |= comp.contains("gas_heat_capacity_J_mol_K");
    allcp &= comp.contains("gas_heat_capacity_J_mol_K");
    anyads |= comp.contains("adsorbed_heat_capacity_J_mol_K");
  }
  Json c = Json::parse(
      R"DEFAULTS({"mu": 1.8e-05, "Cpg": 1000.0, "lambda_g": 0.025, "L": 0.1, "Nz": 10, "dint": 0.01, "dout": 0.012, "eps_b": 0.4, "rho_s": 1000.0, "dp": 0.002, "Cps": 900.0, "hgs": 20.0, "rho_w": 7800.0, "Cpw": 450.0, "lambda_w": 0.0, "hgw": 10.0, "hw_inf": 10.0, "T_inf": 298.15, "THERMAL_MODE": 1, "P_init": 100000.0, "T_init": 298.15, "T_feed": 298.15, "P_FEED": 100000.0, "INITIAL_LOADING_MODE": "equilibrium", "VELOCITY_INPUT_BASIS": "interstitial", "PHYSISORPTION_ENABLED": true, "IAST_METHOD": "AUTO", "IAST_FAST_FALLBACK": true, "SOLVER_METHOD": "BDF", "SOLVER_RTOL": 1e-06, "SOLVER_ATOL": 1e-08, "SOLVER_MAX_STEP": 0.1, "SOLVER_JACOBIAN_MODE": "structured_fd", "ADVECTION_SCHEME": "UPWIND", "RUN_MODE": "fixed_cycles", "N_CYCLES": 2, "MAX_CYCLES": 100, "CSS_STATE_RTOL": 1e-05, "CSS_STATE_ATOL": 1e-08, "CSS_CONSECUTIVE_CYCLES": 3, "SOURCES": {}, "DYNAMIC_BUFFERS": {}, "OUTPUT_DIR": "results", "MAKE_PLOTS": true, "RESTART": false, "RESTART_FROM_CHECKPOINT": false, "ALLOW_EXISTING_OUTPUT_DIR": true, "FEED_TAG": "feed", "LIGHT_PRODUCT_TAG": "light_product", "HEAVY_PRODUCT_TAG": "heavy_product", "CASE_NAME": "RUPTURA cycle", "DESCRIPTION": ""})DEFAULTS");
  c["R"] = R;
  c["comp_names"] = names;
  c["M"] = masses;
  c["Dj"] = dispersion;
  c["y_feed"] = Vec(names.size(), 1. / double(names.size()));
  c["y_init"] = c["y_feed"];
  c["km_phy"] = Json::array({Vec(names.size(), 0.01)});
  c["TARGET_HEAVY_COMPONENT"] = names.front();
  c["TARGET_LIGHT_COMPONENT"] = names.back();
  c["NOTES"] = data.value("notes", Json::array());
  require(!anycp || allcp, "Supply gas heat capacity for every component or omit all");
  if (anycp) {
    c["GAS_COMPONENT_CP_MOLAR"] = Json::array();
    for (const auto& x : components)
      c["GAS_COMPONENT_CP_MOLAR"].push_back(number(x.at("gas_heat_capacity_J_mol_K")));
  }
  if (anyads) {
    c["ADSORBED_CP_MOLAR"] = Json::array();
    for (const auto& x : components)
      c["ADSORBED_CP_MOLAR"].push_back(number(x.value("adsorbed_heat_capacity_J_mol_K", Json(0.))));
  }
  for (auto i = maps().begin(); i != maps().end(); ++i)
    if (data.contains(i.key()) && i.key() != "bed" && i.key() != "thermal" && i.key() != "kinetics")
      mapped(data.at(i.key()), i.key(), c);
  Json common = Json::object();
  for (const std::string k :
       {"bed", "thermal", "kinetics", "equilibrium", "initial", "monolith", "layers"})
    if (data.contains(k)) common[k] = data.at(k);
  update(c, column_sections(common, names));
  require(data.contains("equilibrium"), "Explicit common equilibrium is required");
  if (data.contains("feed")) {
    Json feed = aliases(data.at("feed"));
    allowed(feed, {"pressure", "temperature", "composition", "velocity_basis"}, "feed");
    for (const auto& p : std::vector<std::pair<std::string, std::string>>{
             {"pressure", "P_FEED"},
             {"temperature", "T_feed"},
             {"composition", "y_feed"},
             {"velocity_basis", "VELOCITY_INPUT_BASIS"}})
      if (feed.contains(p.first)) c[p.second] = feed.at(p.first);
  }
  auto columns = data.value("columns", Json{{"A", Json::object()}});
  object(columns, "columns");
  require(columns.size() >= 1 && columns.size() <= 4, "Cycle supports one to four columns");
  c["COLUMNS"] = Json::object();
  for (auto i = columns.begin(); i != columns.end(); ++i)
    c["COLUMNS"][i.key()] = column_sections(i.value(), names, i.key());
  if (data.contains("materials")) {
    object(data.at("materials"), "materials");
    c["MATERIALS"] = Json::object();
    for (auto i = data.at("materials").begin(); i != data.at("materials").end(); ++i)
      c["MATERIALS"][i.key()] = column_sections(i.value(), names, i.key());
  }
  c["SOURCES"] = aliases(data.value("sources", Json::object()));
  c["DYNAMIC_BUFFERS"] = aliases(data.value("buffers", Json::object()));
  c["BUFFERS"] = c["DYNAMIC_BUFFERS"];
  auto cycle = data.value("cycle", Json::object());
  allowed(cycle, {"phases", "steps", "failure_maximum_simulated_time_s"}, "cycle");
  require(!(cycle.contains("phases") && cycle.contains("steps")),
          "Specify either cycle.phases or cycle.steps");
  if (cycle.contains("phases"))
    c["CYCLE_PHASES"] = aliases(cycle.at("phases"));
  else if (cycle.contains("steps"))
    c["CYCLE_STEPS"] = aliases(cycle.at("steps"));
  else
    throw std::invalid_argument("cycle.phases or cycle.steps is required");
  if (cycle.contains("failure_maximum_simulated_time_s"))
    c["STEP_FAILURE_MAXIMUM_SIMULATED_TIME"] = cycle.at("failure_maximum_simulated_time_s");
  if (data.contains("advanced")) update(c, data.at("advanced"));
  c = compositions(c, names);
  c["y_feed"] = composition(c.at("y_feed"), names, "feed.composition");
  c["y_init"] = composition(c.at("y_init"), names, "initial.composition");
  for (const auto& x : c.at("M"))
    require(std::isfinite(number(x)) && number(x) > 0, "Molar masses must be positive");
  require(c["VELOCITY_INPUT_BASIS"] == "interstitial" || c["VELOCITY_INPUT_BASIS"] == "superficial",
          "Invalid feed.velocity_basis");
  require(c["RUN_MODE"] == "fixed_cycles" || c["RUN_MODE"] == "css", "Invalid simulation.run_mode");
  for (const std::string k : {"N_CYCLES", "MAX_CYCLES", "Nz", "CSS_CONSECUTIVE_CYCLES"}) {
    double x = number(c.at(k));
    require(!c.at(k).is_boolean() && std::isfinite(x) && x >= 1 && std::floor(x) == x,
            k + " must be a positive integer");
    c[k] = static_cast<long long>(x);
  }
  if (!c.contains("COLUMN_AREA")) {
    double d = number(c.at("dint"));
    c["COLUMN_AREA"] = PI * d * d / 4.;
  }
  auto base = std::filesystem::absolute(base_dir).lexically_normal();
  auto output = std::filesystem::path(c.at("OUTPUT_DIR").get<std::string>());
  c["OUTPUT_DIR"] = std::filesystem::absolute(output.is_absolute() ? output : base / output)
                        .lexically_normal()
                        .string();
  c["CASE_DIR"] = base.string();
  c["INPUT_DOCUMENT"] = data;
  return c;
}
Json load_config(const std::filesystem::path& supplied) {
  auto path = std::filesystem::absolute(supplied).lexically_normal();
  if (std::filesystem::is_directory(path)) {
    bool found = false;
    for (const std::string n : {"simulation.yaml", "simulation.yml", "simulation.json"})
      if (std::filesystem::is_regular_file(path / n)) {
        path /= n;
        found = true;
        break;
      }
    require(found, "No simulation.yaml or simulation.json in " + path.string());
  }
  auto c = normalize_config(read_input_document(path), path.parent_path());
  c["INPUT_FILE"] = path.string();
  return c;
}
}  // namespace cycle

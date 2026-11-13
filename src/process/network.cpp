#include "cycle/network.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <numeric>
#include <set>
#include <stdexcept>
namespace cycle {
const std::vector<std::string> ColumnNetwork::utility_names = {
    "heating_W",        "cooling_W",        "compressor_W",      "vacuum_W",
    "buffer_heating_W", "buffer_cooling_W", "net_heat_W",        "gas_heating_W",
    "gas_cooling_W",    "pump_W",          "heat_loss_W",      "electric_heater_W", "isothermal_heat_W"};
namespace {
const std::vector<std::string> sides = {"left", "right"},
                               directions = {"left_in", "left_out", "right_in", "right_out"};
std::string lo(std::string s) {
  for (char& c : s) c = std::tolower(static_cast<unsigned char>(c));
  return s;
}
double nget(const Json& j, const std::string& k, double v = 0) {
  return j.contains(k) && !j[k].is_null() ? j[k].get<double>() : v;
}
std::string sget(const Json& j, const std::string& k, const std::string& v = "") {
  return j.contains(k) && !j[k].is_null() ? j[k].get<std::string>() : v;
}
Json obj(const Json& j, const std::string& k) {
  return j.contains(k) && !j[k].is_null() ? j[k] : Json::object();
}
Vec slice(const Vec& v, size_t start, size_t count) {
  if (start + count > v.size()) throw std::invalid_argument("State size mismatch");
  return Vec(v.begin() + start, v.begin() + start + count);
}
Vec flat(Json j, size_t n) {
  if (j.is_number()) return Vec(n, j.get<double>());
  if (j.is_array() && j.size() == 1 && j[0].is_array()) j = j[0];
  Vec a = j.get<Vec>();
  if (a.size() != n) throw std::invalid_argument("Component array size mismatch");
  return a;
}
double sum(const Vec& v) { return std::accumulate(v.begin(), v.end(), 0.); }
double dot(const Vec& a, const Vec& b) {
  double r = 0;
  for (size_t i = 0; i < a.size(); ++i) r += a[i] * b[i];
  return r;
}
Vec mul(Vec v, double a) {
  for (auto& x : v) x *= a;
  return v;
}
Vec diff(Vec a, const Vec& b) {
  for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
  return a;
}
void plus(Vec& a, const Vec& b) {
  for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
}
void addvec(Json& j, const std::string& k, const Vec& v) {
  Vec a = j.contains(k) ? j[k].get<Vec>() : Vec(v.size(), 0);
  plus(a, v);
  j[k] = a;
}
void addnum(Json& j, const std::string& k, double v) { j[k] = nget(j, k) + v; }
bool allclose(const Vec& a, const Vec& b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i)
    if (!std::isfinite(a[i]) || !std::isfinite(b[i]) ||
        std::abs(a[i] - b[i]) > 1e-8 + 1e-5 * std::abs(b[i]))
      return false;
  return true;
}
Json zero_utils() {
  Json j = Json::object();
  for (auto& k : ColumnNetwork::utility_names) j[k] = 0.;
  return j;
}
}  // namespace
ColumnNetwork::ColumnNetwork(const Json& c) : cfg(c) {
  Json cols = cfg.value("COLUMNS", Json());
  legacy_two =
      cols.is_null() && cfg.contains("TWO_COLUMN_PHASES") && !cfg["TWO_COLUMN_PHASES"].empty();
  if (cols.is_null())
    cols = legacy_two ? Json{{"A", Json::object()}, {"B", Json::object()}}
                      : Json{{"A", Json::object()}};
  if (!cols.is_object() || cols.empty() || cols.size() > 4)
    throw std::invalid_argument("COLUMNS must map one to four names");
  for (auto& [name, raw] : cols.items()) {
    if (name.empty() || name.find_first_of(".:") != std::string::npos)
      throw std::invalid_argument("Invalid column name");
    Json local = cfg;
    if (!raw.is_null()) {
      if (!raw.is_object()) throw std::invalid_argument("Column override must be mapping");
      local.update(raw);
    }
    column_indices[name] = models.size();
    column_starts[name] = size;
    column_names.push_back(name);
    models.emplace_back(local);
    size += models.back().layout.size();
  }
  auto& first = models.front();
  n_comp = first.n_comp;
  comp_names = first.comp_names;
  M = first.M;
  R = first.R;
  for (size_t i = 0; i < n_comp; ++i) comp_index[comp_names[i]] = i;
  cp_molar = flat(first.cfg.value("GAS_COMPONENT_CP_MOLAR",
                                  first.cfg.value("CPG_MOLAR", Json(mul(M, first.Cpg)))),
                  n_comp);
  for (auto& m : models) {
    if (cfg.value("SIMULATION_TYPE",std::string("SwingAdsorption"))!="Breakthrough" && m.has_reactions())
      throw std::invalid_argument("General reaction networks are excluded from cycles");
    if(m.liquid!=first.liquid)throw std::invalid_argument("Connected columns must use the same fluid phase");
    if (m.comp_names != comp_names || !allclose(m.M, M))
      throw std::invalid_argument("Connected columns require same components and M");
    auto cp =
        flat(m.cfg.value("GAS_COMPONENT_CP_MOLAR", m.cfg.value("CPG_MOLAR", Json(mul(m.M, m.Cpg)))),
             n_comp);
    if (!allclose(cp, cp_molar))
      throw std::invalid_argument("Connected columns require consistent gas Cp");
    if (std::abs(m.R - R) > 1e-12 * std::abs(R))
      throw std::invalid_argument("Gas constants must match");
    if (m.cfg.contains("GAS_HEAT_CAPACITY_FUNCTION") &&
        !m.cfg["GAS_HEAT_CAPACITY_FUNCTION"].is_null())
      throw std::invalid_argument("Python callable heat capacity is not supported by native input");
  }
  buffers = std::make_unique<DynamicBufferNetwork>(cfg, comp_names, M, R);
  for (auto& d : buffers->definitions)
    if (!allclose(d.cp_molar, cp_molar))
      throw std::invalid_argument("Buffers and columns need consistent Cp");
  buffer_start = size;
  size += buffers->size;
  std::set<std::string> tags = {"feed",  "light_product",        "heavy_product",
                                "purge", "repressurization_gas", "waste"};
  for (auto& name : column_names)
    for (auto& dir : directions) tags.insert(name + "." + dir);
  std::function<void(const Json&)> collect = [&](const Json& j) {
    if (j.is_object()) {
      auto streams = obj(j, "stream_tags");
      if (streams.is_object())
        for (auto& [k, v] : streams.items())
          if (v.is_string()) tags.insert(v.get<std::string>());
      for (auto& [k, v] : j.items()) {
        if ((k == "tag" || (k.size() >= 4 && k.substr(k.size() - 4) == "_tag")) && v.is_string())
          tags.insert(v.get<std::string>());
        collect(v);
      }
    } else if (j.is_array())
      for (auto& v : j) collect(v);
  };
  for (auto k : {"CYCLE_PHASES", "TWO_COLUMN_PHASES", "CYCLE_STEPS", "SOURCES", "BUFFERS",
                 "DYNAMIC_BUFFERS"})
    if (cfg.contains(k)) collect(cfg[k]);
  auto sources = obj(cfg, "SOURCES");
  for (auto& [name, v] : sources.items()) tags.insert(name);
  stream_tags.assign(tags.begin(), tags.end());
}
Column& ColumnNetwork::model(const std::string& n) { return models.at(column_indices.at(n)); }
const Column& ColumnNetwork::model(const std::string& n) const {
  return models.at(column_indices.at(n));
}
Vec ColumnNetwork::column_state(const std::string& n, const Vec& s) const {
  return slice(s, column_starts.at(n), model(n).layout.size());
}
Vec ColumnNetwork::initial_state() {
  if (cfg.contains("TWO_COLUMN_WARM_START_FILE") && !cfg["TWO_COLUMN_WARM_START_FILE"].is_null())
    throw std::invalid_argument("Use native --restart or explicit initial profiles");
  Json defs =
      cfg.value("INITIAL_COLUMN_STATES", cfg.value("TWO_COLUMN_INITIAL_STATES", Json::object()));
  if (defs.is_null()) defs = Json::object();
  if (defs.is_array()) {
    if (defs.size() != models.size()) throw std::invalid_argument("Initial state count mismatch");
    Json d = Json::object();
    for (size_t i = 0; i < models.size(); ++i) d[column_names[i]] = defs[i];
    defs = d;
  }
  if (!defs.is_object()) throw std::invalid_argument("Initial columns must be mapping");
  for (auto& [n, v] : defs.items())
    if (!column_indices.count(n)) throw std::invalid_argument("Unknown column initial definition");
  Json pressures = cfg.value("TWO_COLUMN_INITIAL_PRESSURES", Json()),
       compositions = cfg.value("TWO_COLUMN_INITIAL_COMPOSITIONS", Json());
  if (legacy_two) {
    if (pressures.is_null())
      pressures = Json::array(
          {cfg.value("P_LOW", cfg.at("P_init")), cfg.value("P_HIGH", cfg.at("P_init"))});
    if (compositions.is_null()) compositions = Json::array({cfg.at("y_init"), cfg.at("y_feed")});
  }
  for (auto& v : {pressures, compositions})
    if (!v.is_null() && v.size() != models.size())
      throw std::invalid_argument("Legacy column initial count mismatch");
  Vec s;
  Json explicit_cols = obj(cfg, "COLUMNS");
  for (size_t i = 0; i < models.size(); ++i) {
    auto& name = column_names[i];
    Json d;
    if (defs.contains(name))
      d = defs[name];
    else {
      d = Json::object();
      auto local = obj(explicit_cols, name);
      if (!pressures.is_null() && !local.contains("P_init")) d["pressure"] = pressures[i];
      if (!compositions.is_null() && !local.contains("y_init")) d["composition"] = compositions[i];
      if (d.empty()) d = nullptr;
    }
    Vec v = models[i].initial_state(d);
    s.insert(s.end(), v.begin(), v.end());
  }
  auto b = buffers->initial_state();
  s.insert(s.end(), b.begin(), b.end());
  return s;
}
Json ColumnNetwork::endpoint(const Json& raw) const {
  std::string kind, name, side;
  if (raw.is_object()) {
    if (raw.contains("column")) {
      kind = "column";
      name = raw["column"].get<std::string>();
      side = lo(sget(raw, "side", "right"));
    } else if (raw.contains("buffer")) {
      kind = "buffer";
      name = raw["buffer"].get<std::string>();
    } else if (raw.contains("source")) {
      kind = "source";
      name = raw["source"].get<std::string>();
    } else
      throw std::invalid_argument("Invalid endpoint mapping");
  } else {
    auto text = raw.get<std::string>();
    auto dotpos = text.rfind('.'), colon = text.find(':');
    if (dotpos != std::string::npos) {
      kind = "column";
      name = text.substr(0, dotpos);
      side = lo(text.substr(dotpos + 1));
    } else if (colon != std::string::npos) {
      kind = lo(text.substr(0, colon));
      name = text.substr(colon + 1);
    } else if (buffers->alias_to_name.count(text)) {
      kind = "buffer";
      name = buffers->alias_to_name.at(text);
    } else {
      kind = "source";
      name = text;
    }
  }
  if (kind == "source" && buffers->alias_to_name.count(name)) {
    kind = "buffer";
    name = buffers->alias_to_name.at(name);
  }
  if (kind == "column" && (!column_indices.count(name) || (side != "left" && side != "right")))
    throw std::invalid_argument("Invalid column endpoint");
  if (kind == "buffer") {
    if (buffers->alias_to_name.count(name)) name = buffers->alias_to_name.at(name);
    if (!buffers->name_to_index.count(name)) throw std::invalid_argument("Unknown buffer endpoint");
  }
  if (kind != "column" && kind != "buffer" && kind != "source")
    throw std::invalid_argument("Unknown endpoint kind");
  return Json::array({kind, name, kind == "column" ? Json(side) : Json()});
}
void ColumnNetwork::set_phase(const Json& p, const Json& ctx) {
  transport_cache_time_=std::numeric_limits<double>::quiet_NaN();transport_cache_state_.clear();
  phase = p;
  context = ctx;
  Json rawactions = obj(p, "columns");
  if (rawactions.empty() && models.size() == 1 && p.contains("left_port"))
    rawactions = {{column_names.front(), p}};
  for (auto& [name, a] : rawactions.items())
    if (!column_indices.count(name)) throw std::invalid_argument("Unknown phase column");
  actions = Json::object();
  connections = Json::array();
  Json raws = p.value("connections", Json::array());
  if (raws.is_null()) raws = Json::array();
  if (p.contains("connection") && !p["connection"].is_null()) raws.push_back(p["connection"]);
  for (size_t i = 0; i < raws.size(); ++i) {
    Json d = raws[i];
    if (!d.contains("from")) {
      d["from"] = sget(d, "donor") + "." + sget(d, "donor_side", "right");
      d["to"] = sget(d, "receiver") + "." + sget(d, "receiver_side", "right");
      if (!d.contains("check_valve")) d["check_valve"] = true;
    }
    d["_a"] = endpoint(d.at("from"));
    d["_b"] = endpoint(d.at("to"));
    if (d["_a"] == d["_b"]) throw std::invalid_argument("Self connection");
    if (!d.contains("name")) d["name"] = "connection_" + std::to_string(i + 1);
    connections.push_back(d);
  }
  for (auto& name : column_names) {
    Json action = obj(rawactions, name);
    if (!action.contains("name")) action["name"] = sget(p, "name", "phase");
    for (auto& side : sides) {
      auto key = side + "_port";
      if (!action.contains(key)) action[key] = {{"mode", "closed"}};
      Json port = action[key];
      auto mode = lo(sget(port, "mode", "closed"));
      if (mode == "valve" || mode == "compressor") {
        auto source = port.value("source", port.value("buffer", port.value("sink", Json())));
        if (source.is_null())
          throw std::invalid_argument("Valve/compressor must name source or buffer");
        Json d = obj(port, "valve");
        for (auto& [k, v] : port.items())
          if (k != "mode" && k != "source" && k != "buffer" && k != "sink" && k != "valve")
            d[k] = v;
        d["_a"] = endpoint(source);
        d["_b"] = Json::array({"column", name, side});
        d["name"] = name + "." + side;
        if (mode == "compressor" && (!d.contains("compressor") || d["compressor"].empty()))
          throw std::invalid_argument("Compressor mapping required");
        if (!d.contains("tag")) d["tag"] = obj(action, "stream_tags").value(side + "_in", source);
        connections.push_back(d);
        action[key] = {{"mode", "closed"}};
      }
    }
    actions[name] = action;
  }
  auto s = initial_state();
  refresh_context(0, slice(s, buffer_start, buffers->size));
}
Json ColumnNetwork::refresh_context(double t, const Vec& b) {
  Json sources = configured_sources(cfg, n_comp, t, context);
  sources.update(obj(context, "source_states"));
  for (auto& d : buffers->definitions) {
    auto r = buffers->state_record(d.name, b);
    sources[d.name] = r;
    for (auto& a : d.aliases) sources[a] = r;
  }
  Json ctx = context;
  ctx["source_states"] = sources;
  for (auto& name : column_names) {
    Json local = ctx;
    local["phase_start_pressure"] = obj(obj(ctx, "phase_start_pressures"), name);
    model(name).set_step(actions[name], local);
  }
  return sources;
}
Json ColumnNetwork::endpoint_pressures(const Vec& s) const {
  Json j = Json::object();
  for (auto& name : column_names) {
    auto d = model(name).derived(column_state(name, s));
    j[name] = {{"left", d.P.front()}, {"right", d.P.back()}};
  }
  return j;
}
Json ColumnNetwork::record(const Json& ep, const Vec& s, const std::vector<Derived>& ds,
                           const Json& sources) const {
  auto kind = ep[0].get<std::string>(), name = ep[1].get<std::string>();
  if (kind == "column") {
    auto& d = ds[column_indices.at(name)];
    auto& m = model(name);
    size_t z = ep[2] == "left" ? 0 : m.Nz - 1;
    Vec y(n_comp),c(n_comp);
    for (size_t i = 0; i < n_comp; ++i){ y[i] = d.y[i * m.Nz + z];c[i]=d.c[i*m.Nz+z];}
    return {{"pressure", d.P[z]},
            {"temperature", d.Tg[z]},
            {"composition", y},{"concentration",c},
            {"kind", kind},
            {"name", name}};
  }
  if (kind == "buffer") return buffers->state_record(name, slice(s, buffer_start, buffers->size));
  Json r = sources.at(name);
  r["name"] = name;
  return r;
}
std::pair<Vec, Json> ColumnNetwork::evaluate(double t, const Vec& s) {
  if (s.size() != size) throw std::invalid_argument("Network state size mismatch");
  if (actions.empty()) throw std::runtime_error("set_phase must precede evaluate");
  Json sources = refresh_context(t, slice(s, buffer_start, buffers->size));
  std::vector<Derived> ds;
  for (auto& name : column_names) ds.push_back(model(name).derived(column_state(name, s)));
  const bool liquid=models.front().liquid;
  if (!liquid)
    for (const auto& link : connections)
      for (const auto key : {"_a", "_b"})
        if (link[key][0] == "column")
        {
          const auto name = link[key][1].get<std::string>();
          if (obj(actions.at(name), "flow_control").value("model", std::string{}) == "fixed_velocity")
            throw std::invalid_argument("Gas FixedVelocity uses prescribed ports; use pressure/velocity ports for connected columns");
        }
  const Json hydraulics=liquid?solve_liquid_hydraulics(t,s,sources):Json::object();
  Json volumes=Json::object(),buffer_volumes=Json::object();
  Json boundary = Json::object(), power = Json::object(), gross = Json::object(),
       bdn = Json::object(), bh = Json::object(), cu = Json::object();
  for (auto& name : column_names) {
    boundary[name] = Json::object();volumes[name]={{"left",0.},{"right",0.}};
    power[name] = Json::object();
    gross[name] = Json::object();
    cu[name] = zero_utils();
    for (auto& side : sides) {
      boundary[name][side] = Vec(n_comp, 0);
      power[name][side] = 0.;
    }
    for (auto& dir : directions) gross[name][dir] = Vec(n_comp, 0);
  }
  for (auto& d : buffers->definitions) {
    bdn[d.name] = Vec(n_comp, 0);
    bh[d.name] = 0.;
  }
  Json diag = {{"port_rates", Json::object()},   {"external_in", Vec(n_comp, 0)},
               {"external_out", Vec(n_comp, 0)}, {"external_stream_rates", Json::object()},
               {"provenance", Json::array()},    {"utility_rates_W", zero_utils()},
               {"column_utility_rates_W", cu},   {"internal_transfers", Json::array()},
               {"connections", Json::array()}, {"volume_provenance",Json::array()},   {"external_enthalpy_in_W", 0.},
               {"external_enthalpy_out_W", 0.}};
  auto H = [&](double T) { return models.front().gas_component_enthalpy(T); };
  auto ledger = [&](const std::string& tag, const std::string& direction, const Vec& amounts,
                    const std::string& source, const std::string& origin, const Json& product) {
    if (std::all_of(amounts.begin(), amounts.end(), [](double x) { return x == 0; })) return;
    Json entry = {{"tag", tag},
                  {"direction", direction},
                  {"component_rates", amounts},
                  {"source", source},
                  {"origin", origin}};
    if (!product.is_null()) entry["product_tag"] = product;
    diag["provenance"].push_back(entry);
    if (origin == "external" || origin == "product_recycle") {
      addvec(diag["external_stream_rates"], tag, amounts);
      addvec(diag, "external_" + direction, amounts);
    }
  };
  auto add = [&](const Json& ep, const Vec& amounts, double enthalpy, const Json& def) {
    auto kind = ep[0].get<std::string>(), name = ep[1].get<std::string>();
    if (kind == "column") {
      auto side = ep[2].get<std::string>();
      addvec(boundary[name], side, amounts);
      addnum(power[name], side, enthalpy);
      auto dir = sum(amounts) >= 0 ? "in" : "out";
      Vec a = amounts;
      for (auto& v : a) v = std::abs(v);
      addvec(gross[name], side + "_" + dir, a);
    } else if (kind == "buffer") {
      addvec(bdn, name, amounts);
      addnum(bh, name, enthalpy);
    } else {
      bool positive = sum(amounts) >= 0;
      std::string dir = positive ? "out" : "in";
      auto tag = sget(def, dir + "_tag", sget(def, "tag", sget(sources[name], dir + "_tag", name)));
      ledger(tag, dir, positive ? amounts : mul(amounts, -1), name, "external", nullptr);
      addnum(diag, "external_enthalpy_" + dir + "_W", std::abs(enthalpy));
    }
  };
  for (size_t ci = 0; ci < models.size(); ++ci) {
    auto& name = column_names[ci];
    auto& m = models[ci];
    auto v = column_state(name, s);
    const auto hydraulic=liquid?hydraulics["columns"][name]:Json::object();
    auto rates = m.boundary_component_rates(t, v,hydraulic);
    auto enthalpies=m.boundary_enthalpy_rates(t,v,hydraulic);
    for (auto& side : sides) {
      if(liquid){bool connected=false;for(auto& link:connections)for(auto endpoint_key:{"_a","_b"}){const auto& ep=link[endpoint_key];connected|=ep[0]=="column"&&ep[1]==name&&ep[2]==side;}if(connected)continue;double U=hydraulic[side]["volume_flux"];volumes[name][side]=(side=="left"?1.:-1.)*U*m.area;}
      auto incoming = rates.at(side + "_in").get<Vec>(),
           outgoing = rates.at(side + "_out").get<Vec>();
      addvec(gross[name], side + "_in", incoming);
      addvec(gross[name], side + "_out", outgoing);
      auto port = actions[name][side + "_port"], resolved = m.port_state(side, t);
      size_t z = side == "left" ? 0 : m.Nz - 1;
      double Tin = nget(resolved, "temperature", ds[ci].Tg[z]), hin = dot(incoming, H(Tin)),
             hout = dot(outgoing, H(ds[ci].Tg[z]));
      if(liquid){hin=enthalpies[side+"_in"];hout=enthalpies[side+"_out"];}
      addvec(boundary[name], side, diff(incoming, outgoing));
      addnum(power[name], side, hin - hout);
      auto b = buffers->buffer_name_for_port(port);
      if (b) {
        if (sget(resolved, "mode") == "velocity")
          throw std::invalid_argument("Finite buffer velocity port needs valve or compressor");
        addvec(bdn, *b, diff(outgoing, incoming));
        addnum(bh, *b, hout - hin);
        if(liquid)addnum(buffer_volumes,*b,-nget(volumes[name],side));
        for (auto dir : {"in", "out"})
          ledger(sget(obj(actions[name], "stream_tags"), side + "_" + dir,
                      name + "." + side + "_" + dir),
                 dir, std::string(dir) == "in" ? incoming : outgoing, *b, "buffer", nullptr);
      } else {
        if(liquid){double q=nget(volumes[name],side);if(q!=0.)diag["volume_provenance"].push_back({{"direction",q>0.?"in":"out"},{"tag",sget(obj(actions[name],"stream_tags"),side+(q>0.?"_in":"_out"),name+"."+side+(q>0.?"_in":"_out"))},{"rate_m3_s",std::abs(q)}});}
        addnum(diag, "external_enthalpy_in_W", hin);
        addnum(diag, "external_enthalpy_out_W", hout);
        for (auto dir : {"in", "out"}) {
          auto tag = sget(obj(actions[name], "stream_tags"), side + "_" + dir,
                          name + "." + side + "_" + dir);
          Json history =
              std::string(dir) == "in" ? port.value("composition_source", Json()) : Json();
          std::string origin = (history == "light_product" || history == "heavy_product")
                                   ? "product_recycle"
                                   : "external";
          ledger(tag, dir, std::string(dir) == "in" ? incoming : outgoing,
                 sget(port, "source", name + "." + side), origin, history);
        }
        if(liquid && nget(volumes[name],side)>0.) {
        double eta=nget(cfg,"LIQUID_PUMP_EFFICIENCY",.75);
        if(!(eta>0.&&eta<=1.))throw std::invalid_argument("PumpEfficiency must lie in (0,1]");
        double pressure=hydraulics["columns"][name][side]["pressure"];
        double power=nget(volumes[name],side)*std::max(pressure-nget(resolved,"pressure"),0.)/eta;
        addnum(diag["utility_rates_W"],"pump_W",power);addnum(diag["column_utility_rates_W"][name],"pump_W",power);
      }
      if (sum(incoming) > 0 && !liquid) {
          auto src = obj(sources, sget(port, "source"));
          double reference = nget(src, "temperature",
                                  nget(actions[name], "gas_utility_reference_temperature",
                                       nget(cfg, "T_feed", 298.15))),
                 conditioning = dot(incoming, diff(H(Tin), H(reference)));
          addnum(diag["utility_rates_W"], "gas_heating_W", std::max(conditioning, 0.));
          addnum(diag["utility_rates_W"], "gas_cooling_W", std::max(-conditioning, 0.));
          addnum(diag["column_utility_rates_W"][name], "gas_heating_W", std::max(conditioning, 0.));
          addnum(diag["column_utility_rates_W"][name], "gas_cooling_W",
                 std::max(-conditioning, 0.));
        }
      }
    }
  }
  for (auto& def : connections) {
    Json a = def["_a"], b = def["_b"], ra = record(a, s, ds, sources),
         rb = record(b, s, ds, sources);
    if(liquid) {
      double q=hydraulics["connections"][def["name"].get<std::string>()];
      Json donor=q>=0? a:b,receiver=q>=0?b:a,dr=q>=0?ra:rb;
      Vec rate=mul(dr.at("concentration").get<Vec>(),std::abs(q));
      double Hrate=std::abs(q)*models.front().liquid_density*models.front().liquid_heat_capacity*nget(dr,"temperature");
      auto endpoint_pressure=[&](const Json& ep,const Json& resource){return ep[0]=="column"?hydraulics["columns"][ep[1].get<std::string>()][ep[2].get<std::string>()]["pressure"].get<double>():nget(resource,"pressure");};
      double pa=endpoint_pressure(a,ra),pb=endpoint_pressure(b,rb),eta=nget(cfg,"LIQUID_PUMP_EFFICIENCY",.75);
      if(!(eta>0.&&eta<=1.))throw std::invalid_argument("PumpEfficiency must lie in (0,1]");
      double pump=std::abs(q)*std::max(q>=0?pb-pa:pa-pb,0.)/eta;
      addnum(diag["utility_rates_W"],"pump_W",pump);if(receiver[0]=="column")addnum(diag["column_utility_rates_W"][receiver[1].get<std::string>()],"pump_W",pump);
      add(donor,mul(rate,-1.),-Hrate,def);add(receiver,rate,Hrate,def);
      auto volume=[&](const Json& ep,double value){auto kind=ep[0].get<std::string>(),name=ep[1].get<std::string>();if(kind=="column")addnum(volumes[name],ep[2].get<std::string>(),value);else if(kind=="buffer")addnum(buffer_volumes,name,value);};
      volume(donor,-std::abs(q));volume(receiver,std::abs(q));
      if(donor[0]=="source")diag["volume_provenance"].push_back({{"direction","in"},{"tag",sget(def,"in_tag",sget(def,"tag",donor[1]))},{"rate_m3_s",std::abs(q)}});
      if(receiver[0]=="source")diag["volume_provenance"].push_back({{"direction","out"},{"tag",sget(def,"out_tag",sget(def,"tag",receiver[1]))},{"rate_m3_s",std::abs(q)}});
      diag["connections"].push_back({{"name",def["name"]},{"donor",donor},{"receiver",receiver},{"component_rates",rate},{"volume_flow_m3_s",std::abs(q)},{"signed_molar_flow",q*(q>=0?sum(ra.at("concentration").get<Vec>()):sum(rb.at("concentration").get<Vec>()))},{"enthalpy_out_W",Hrate},{"enthalpy_in_W",Hrate}});
      continue;
    }
    Json up = nget(ra, "pressure") >= nget(rb, "pressure") ? ra : rb;
    double cp = dot(up["composition"].get<Vec>(), cp_molar);
    Json kw = Json::object();
    for (auto k : {"flow_model", "conductance", "area_m2", "Cd", "gamma", "opening", "check_valve",
                   "maximum_molar_flow"})
      if (def.contains(k)) kw[k] = def[k];
    if (!kw.contains("gamma")) kw["gamma"] = cp / (cp - R);
    if (def.contains("Cv") || def.contains("cv"))
      throw std::invalid_argument("Use physical area_m2 and Cd; Cv unsupported");
    double flow = 0, work = 0;
    Json donor, receiver, dr, rr;
    Json comp = obj(def, "compressor");
    if (!comp.empty()) {
      double discharge = nget(comp, "outlet_pressure",
                              std::max(nget(rb, "pressure"), nget(ra, "pressure"))),
             icp = dot(ra["composition"].get<Vec>(), cp_molar);
      if (!def.contains("gamma")) kw["gamma"] = icp / (icp - R);
      if (nget(ra, "pressure") > 0) {
        auto [tout, wpm] =
            compressor_outlet(nget(ra, "pressure"), discharge, nget(ra, "temperature"), icp,
                              nget(comp, "efficiency", .75), R);
        if (discharge < nget(rb, "pressure"))
          flow = 0;
        else if (comp.contains("flow_mol_s")) {
          flow = nget(comp, "flow_mol_s");
          if (flow < 0 || !std::isfinite(flow))
            throw std::invalid_argument("Invalid compressor flow");
          double opening = nget(kw, "opening", 1), lim = nget(kw, "maximum_molar_flow", INFINITY);
          if (opening < 0 || opening > 1 || !std::isfinite(opening) || lim <= 0 || std::isnan(lim))
            throw std::invalid_argument("Invalid compressor opening or limit");
          flow = std::min(flow * opening, lim);
        } else {
          kw["check_valve"] = true;
          flow = valve_molar_flow(discharge, nget(rb, "pressure"), tout,
                                  dot(ra["composition"].get<Vec>(), M), kw, R);
        }
        work = flow * wpm;
      }
      donor = a;
      receiver = b;
      dr = ra;
      rr = rb;
    } else {
      flow = valve_molar_flow(nget(ra, "pressure"), nget(rb, "pressure"), nget(up, "temperature"),
                              dot(up["composition"].get<Vec>(), M), kw, R);
      donor = flow >= 0 ? a : b;
      receiver = flow >= 0 ? b : a;
      dr = flow >= 0 ? ra : rb;
      rr = flow >= 0 ? rb : ra;
    }
    Vec rate = mul(dr["composition"].get<Vec>(), std::abs(flow));
    double hd = dot(rate, H(nget(dr, "temperature"))), hr = hd + work;
    add(donor, mul(rate, -1), -hd, def);
    add(receiver, rate, hr, def);
    addnum(diag["utility_rates_W"], "compressor_W", work);
    Json transfer = {{"name", def["name"]},
                     {"donor", donor},
                     {"receiver", receiver},
                     {"component_rates", rate},
                     {"molar_flow", std::abs(flow)},
                     {"signed_molar_flow", flow},
                     {"pressure_difference", nget(ra, "pressure") - nget(rb, "pressure")},
                     {"enthalpy_out_W", hd},
                     {"enthalpy_in_W", hr},
                     {"compressor_W", work}};
    diag["connections"].push_back(transfer);
    if (donor[0] != "source" && receiver[0] != "source") {
      diag["internal_transfers"].push_back(transfer);
      for (int ii = 0; ii < 2; ++ii) {
        auto ep = ii == 0 ? donor : receiver, other = ii == 0 ? receiver : donor;
        std::string dir = ii == 0 ? "out" : "in";
        if (ep[0] != "column") continue;
        auto name = ep[1].get<std::string>(), side = ep[2].get<std::string>();
        auto dt = sget(def, dir + "_tag", dir == "in" ? sget(def, "tag") : "");
        auto tag = sget(obj(actions[name], "stream_tags"), side + "_" + dir,
                        dt.empty() ? name + "." + side + "_" + dir : dt);
        ledger(tag, dir, rate, other[1].get<std::string>(),
               other[0] == "buffer" ? "buffer" : "column", nullptr);
      }
    }
  }
  Vec derivative(size);
  Json overridesall = Json::object();
  for (size_t ci = 0; ci < models.size(); ++ci) {
    auto& name = column_names[ci];
    auto& m = models[ci];
    auto& d = ds[ci];
    Json overrides = Json::object();
    for (auto& side : sides) {
      size_t z = side == "left" ? 0 : m.Nz - 1;
      double sign = side == "left" ? 1 : -1;
      Vec flux = mul(boundary[name][side].get<Vec>(), sign / m.area), c(n_comp);
      for (size_t i = 0; i < n_comp; ++i) c[i] = d.c[i * m.Nz + z];
      overrides[side] = {{"component_flux", flux},
                         {"enthalpy_flux", sign * nget(power[name], side) / m.area},
                         {"velocity", sum(flux) / std::max(m.epsilon_node[z] * sum(c), 1e-300)},
                         {"concentration", c},
                         {"temperature", d.Tg[z]}};
    }
    if(liquid)for(auto& side:sides){double sign=side=="left"?1.:-1.;auto z=side=="left"?0:m.Nz-1;double U=sign*nget(volumes[name],side)/m.area;overrides[side]["velocity"]=U/m.epsilon_node[z];overrides[side]["volume_flux"]=U;overrides[side]["pressure"]=hydraulics["columns"][name][side]["pressure"];}
    overridesall[name] = overrides;
    auto v = column_state(name, s), dv = m.rhs(t, v, overrides);
    std::copy(dv.begin(), dv.end(), derivative.begin() + column_starts[name]);
    addvec(diag,"reaction_component_rates",m.reaction_component_rates());
    addnum(diag,"reaction_energy_W",m.reaction_energy_rate());
    auto rates = gross[name];
    for (auto& side : sides)
      rates[side] = mul(boundary[name][side].get<Vec>(), side == "left" ? 1 : -1);
    diag["port_rates"][name] = rates;
    auto u = m.thermal_utility_rates(t, v, overrides, &dv);
    for (auto& [k, val] : u.items())
      if (k != "gas_heating_W" && k != "gas_cooling_W") {
        addnum(diag["utility_rates_W"], k, val.get<double>());
        addnum(diag["column_utility_rates_W"][name], k, val.get<double>());
      }
  }
  auto [bd, bdiag] =
      buffers->conservative_derivatives(slice(s, buffer_start, buffers->size), bdn, bh, sources,buffer_volumes);
  std::copy(bd.begin(), bd.end(), derivative.begin() + buffer_start);
  for (auto k : {"external_in", "external_out"}) addvec(diag, k, bdiag[k].get<Vec>());
  for (auto k : {"external_enthalpy_in_W", "external_enthalpy_out_W"})
    addnum(diag, k, bdiag[k].get<double>());
  for (auto& p : bdiag["provenance"]) diag["provenance"].push_back(p);
  for (auto& [k, v] : bdiag["external_stream_rates"].items())
    addvec(diag["external_stream_rates"], k, v.get<Vec>());
  for (auto& [k, v] : bdiag["utility_rates_W"].items())
    addnum(diag["utility_rates_W"], k, v.get<double>());
  diag["boundary_overrides"] = overridesall;
  transport_cache_state_=s;transport_cache_time_=t;transport_boundary_overrides_=overridesall;
  return {derivative, diag};
}
Vec ColumnNetwork::explicit_transport_rhs(double t, const Vec& state) {
  if(transport_cache_time_!=t||transport_cache_state_!=state)evaluate(t,state);
  Vec result(size,0.);
  for (const auto& name : column_names) {
    const auto derivative = model(name).explicit_transport_rhs(t,column_state(name,state),
        transport_boundary_overrides_.at(name));
    std::copy(derivative.begin(),derivative.end(),result.begin()+column_starts.at(name));
  }
  return result;
}
Json ColumnNetwork::inventories(const Vec& s) const {
  if (s.size() != size) throw std::invalid_argument("Network state size mismatch");
  Json cols = Json::object();
  Vec gas(n_comp, 0), ads(n_comp, 0), pore(n_comp, 0), bm(n_comp, 0);
  double e = 0;
  for (auto& n : column_names) {
    auto r = model(n).inventories(column_state(n, s));
    cols[n] = r;
    plus(gas, r["gas"].get<Vec>());
    plus(ads, r["adsorbed"].get<Vec>());
    plus(pore, r["surface_and_pore"].get<Vec>());
    e += nget(r, "energy_J");
  }
  auto records = buffers->records(slice(s, buffer_start, buffers->size));
  for (auto& [n, r] : records.items()) {
    plus(bm, r["component_moles"].get<Vec>());
    e += nget(r, "energy_J");
  }
  plus(gas, bm);
  Vec total = gas;
  plus(total, ads);
  plus(total, pore);
  return {{"gas", gas},      {"adsorbed", ads},    {"total", total},    {"energy_J", e},
          {"surface_and_pore", pore}, {"columns", cols}, {"buffers", records}, {"buffer_moles", bm}};
}
Json ColumnNetwork::state_records(const Vec& s) const {
  if (s.size() != size) throw std::invalid_argument("Network state size mismatch");
  Json cols = Json::object();
  for (auto& n : column_names) cols[n] = model(n).state_record(column_state(n, s));
  return {{"columns", cols}, {"buffers", buffers->records(slice(s, buffer_start, buffers->size))}};
}
Vec ColumnNetwork::state_atol(double def) const {
  Vec v(size, def);
  for (auto& n : column_names) {
    auto a = model(n).state_atol(def);
    std::copy(a.begin(), a.end(), v.begin() + column_starts.at(n));
  }
  for (size_t i = 0; i < buffers->definitions.size(); ++i)
    v[buffer_start + i * buffers->block_size + n_comp] = std::max(def, 1e-5);
  if(buffers->liquid)for(std::size_t i=0;i<buffers->definitions.size();++i)v[buffer_start+i*buffers->block_size+n_comp+1]=cfg.value("ATOL_LIQUID_VOLUME",1e-12);
  return v;
}
double ColumnNetwork::metric(const std::string& metric, const Vec& s, double t,
                             const Json& options) {
  auto k = lo(metric);
  Json o = options;
  Json component = o.value("component", Json());
  if (o.contains("column") && !o["column"].is_null()) {
    auto name = o["column"].get<std::string>();
    bool fm = k == "port_velocity" || k == "absolute_port_velocity" ||
              k == "port_superficial_velocity" || k == "absolute_port_superficial_velocity" ||
              k == "port_total_molar_flow",
         om = k == "outlet_mole_fraction" || k == "outlet_concentration_ratio";
    Json selected = o.value("port", o.value("side", Json()));
    if (fm || (om && selected.is_null())) {
      auto diag = evaluate(t, s).second, rates = diag["port_rates"][name];
      if (om) {
        double l = sum(rates["left_out"].get<Vec>()), r = sum(rates["right_out"].get<Vec>());
        if (std::max(l, r) <= 0) return 0;
        selected = l > r ? "left" : "right";
        o["port"] = selected;
      } else {
        if (selected != "left" && selected != "right")
          throw std::invalid_argument("Port flow metric needs left/right port");
        auto side = selected.get<std::string>();
        if (k == "port_total_molar_flow")
          return sum(diff(rates[side + "_in"].get<Vec>(), rates[side + "_out"].get<Vec>()));
        double v = diag["boundary_overrides"][name][side]["velocity"];
        if (k.find("superficial") != std::string::npos)
          v *= side == "left" ? model(name).epsilon_node.front() : model(name).epsilon_node.back();
        return k.rfind("absolute", 0) == 0 ? std::abs(v) : v;
      }
    }
    return model(name).metric(metric, column_state(name, s), t, o);
  }
  if (o.contains("buffer") && !o["buffer"].is_null())
    return buffers->metric(metric, slice(s, buffer_start, buffers->size),
                           o["buffer"].get<std::string>(), component);
  if (k == "combined_total_inventory" || k == "total_inventory") {
    auto v = inventories(s)["total"].get<Vec>();
    if (component.is_null()) return sum(v);
    return v.at(component.is_string() ? comp_index.at(component.get<std::string>())
                                      : component.get<size_t>());
  }
  if (k.find("pressure_difference") != std::string::npos ||
      k.find("equalization_molar_flow") != std::string::npos) {
    double v;
    if (k.find("column_average") != std::string::npos) {
      if (models.size() != 2)
        throw std::invalid_argument("Column-average difference requires two columns");
      v = models[0].metric("bed_average_pressure", column_state(column_names[0], s)) -
          models[1].metric("bed_average_pressure", column_state(column_names[1], s));
    } else {
      auto tr = evaluate(t, s).second["connections"];
      Json selected = Json::array();
      for (auto& x : tr)
        if (!o.contains("connection") || o["connection"].is_null() || x["name"] == o["connection"])
          selected.push_back(x);
      if (selected.size() != 1)
        throw std::invalid_argument("Connection metric requires one selected connection");
      v = selected[0][k.find("pressure_difference") != std::string::npos ? "pressure_difference"
                                                                         : "signed_molar_flow"];
    }
    return k.rfind("absolute", 0) == 0 ? std::abs(v) : v;
  }
  if (models.size() == 1) {
    o["column"] = column_names.front();
    return this->metric(metric, s, t, o);
  }
  throw std::invalid_argument("Network metric needs column, buffer or connection");
}
}  // namespace cycle

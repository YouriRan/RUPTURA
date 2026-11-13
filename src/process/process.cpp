#include "cycle/process.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
namespace cycle {
namespace {
double num(const Json& j, const char* k, double d = 0) {
  return j.contains(k) && !j[k].is_null() ? j[k].get<double>() : d;
}
std::string str(const Json& j, const char* k, const std::string& d = "") {
  return j.contains(k) && !j[k].is_null() ? j[k].get<std::string>() : d;
}
std::string low(std::string s) {
  for (char& c : s) c = std::tolower(static_cast<unsigned char>(c));
  return s;
}
Vec arr(Json j, size_t n) {
  if (j.is_number()) return Vec(n, j.get<double>());
  if (j.is_array() && j.size() == 1 && j[0].is_array()) j = j[0];
  auto v = j.get<Vec>();
  if (v.size() != n) throw std::invalid_argument("Component vector has incorrect size");
  return v;
}
Vec norm(const Json& j, size_t n) {
  Vec a = arr(j, n);
  double sum = 0;
  for (double v : a) {
    if (!std::isfinite(v) || v < 0)
      throw std::invalid_argument("Composition must be finite and nonnegative");
    sum += v;
  }
  if (sum <= 0) throw std::invalid_argument("Composition sum must be positive");
  for (double& v : a) v /= sum;
  return a;
}
double dot(const Vec& a, const Vec& b) {
  double x = 0;
  for (size_t i = 0; i < a.size(); ++i) x += a[i] * b[i];
  return x;
}
Vec scale(Vec a, double x) {
  for (auto& v : a) v *= x;
  return a;
}
void add(Vec& a, const Vec& b, double x = 1) {
  for (size_t i = 0; i < a.size(); ++i) a[i] += x * b[i];
}
void addj(Json& j, const std::string& k, const Vec& a) {
  Vec v = j.contains(k) ? j[k].get<Vec>() : Vec(a.size(), 0);
  add(v, a);
  j[k] = v;
}
}  // namespace
Json configured_sources(const Json& cfg, size_t nc, double, const Json&) {
  Json s = Json::object();
  if (cfg.contains("y_feed"))
    s["feed"] = {{"composition", cfg["y_feed"]},
                 {"pressure", num(cfg, "P_FEED", num(cfg, "P_HIGH", num(cfg, "P_init", 1e5)))},
                 {"temperature", num(cfg, "T_feed", num(cfg, "T_init", 298.15))},
                 {"kind", "reservoir"}};
  if (cfg.contains("y_init"))
    s["initial_gas"] = {{"composition", cfg["y_init"]},
                        {"pressure", num(cfg, "P_init", 1e5)},
                        {"temperature", num(cfg, "T_init", 298.15)},
                        {"kind", "reservoir"}};
  Json defs = cfg.value("SOURCES", Json::object());
  if (defs.is_null()) defs = Json::object();
  if (!defs.is_object()) throw std::invalid_argument("SOURCES must be a mapping");
  for (auto& [name, d] : defs.items()) {
    if (!d.is_object()) throw std::invalid_argument("Source must be a mapping");
    if (low(str(d, "kind", "reservoir")) == "buffer") continue;
    const bool liquid = low(str(cfg, "FLUID_PHASE", "gas")) == "liquid";
    Vec y = norm(liquid ? d.value("composition", cfg.value("y_feed", Json(Vec(nc, 1.))))
                        : d.at("composition"), nc);
    double p = d.at("pressure").get<double>(), t = num(d, "temperature", 298.15);
    if (!std::isfinite(p) || p < 0 || !std::isfinite(t) || t <= 0)
      throw std::invalid_argument("Invalid source P or T");
    s[name] = {{"composition", y}, {"pressure", p}, {"temperature", t}, {"kind", "reservoir"}};
  }
  if (low(str(cfg,"FLUID_PHASE","gas"))=="liquid") {
    for (auto& [name,state] : s.items()) {
      Json c = name=="feed" ? cfg.value("LIQUID_FEED_CONCENTRATION",Json(Vec(nc,0.))) : cfg.value("LIQUID_INITIAL_CONCENTRATION",Json(Vec(nc,0.)));
      if(defs.contains(name)) c=defs[name].value("concentration",c);
      state["concentration"]=arr(c,nc);
    }
  }
  return s;
}
double valve_molar_flow(double pa, double pb, double temp, double mw, const Json& o, double R) {
  if (!std::isfinite(pa) || !std::isfinite(pb) || !std::isfinite(temp) || !std::isfinite(mw) ||
      !std::isfinite(R) || std::min(pa, pb) < 0 || std::min({temp, mw, R}) <= 0)
    throw std::invalid_argument("Invalid valve P,T,M or R");
  double opening = num(o, "opening", 1), maximum = num(o, "maximum_molar_flow", INFINITY);
  if (!std::isfinite(opening) || opening < 0 || opening > 1 || std::isnan(maximum) || maximum <= 0)
    throw std::invalid_argument("Invalid valve opening or maximum flow");
  double delta = pa - pb;
  if (delta == 0 || opening == 0 || (o.value("check_valve", false) && delta < 0)) return 0;
  double sign = delta > 0 ? 1 : -1, up = std::max(pa, pb), down = std::min(pa, pb), flow = 0;
  auto kind = low(str(o, "flow_model", "linear"));
  if (kind == "linear" || kind == "sqrt" || kind == "square_root") {
    double k = num(o, "conductance", 1e-9);
    if (!std::isfinite(k) || k < 0) throw std::invalid_argument("Invalid valve conductance");
    flow = k * (kind == "linear" ? std::abs(delta) : std::sqrt(std::abs(delta)));
  } else if (kind == "ideal_gas" || kind == "orifice" || kind == "compressible") {
    if (!o.contains("area_m2")) throw std::invalid_argument("Ideal gas valve requires area_m2");
    double area = num(o, "area_m2"), cd = num(o, "Cd", .8), g = num(o, "gamma", 1.4);
    if (!std::isfinite(area) || !std::isfinite(cd) || !std::isfinite(g) || area <= 0 || cd <= 0 ||
        cd > 1 || g <= 1)
      throw std::invalid_argument("Invalid orifice area, Cd or gamma");
    double ratio = down / up, critical = std::pow(2 / (g + 1), g / (g - 1)), factor;
    if (ratio <= critical)
      factor = std::sqrt(g) * std::pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
    else {
      double lr = std::log(ratio), term = std::exp(2 * lr / g) * (-std::expm1((g - 1) * lr / g));
      factor = std::sqrt(std::max(2 * g / (g - 1) * term, 0.));
    }
    flow = cd * area * up / std::sqrt(mw * R * temp) * factor;
  } else
    throw std::invalid_argument("Unknown valve flow model");
  return sign * std::min(opening * flow, maximum);
}
std::pair<double, double> compressor_outlet(double pin, double pout, double tin, double cp,
                                            double eta, double R) {
  if (!std::isfinite(pin) || !std::isfinite(pout) || !std::isfinite(tin) || !std::isfinite(cp) ||
      !std::isfinite(eta) || !std::isfinite(R) || std::min({pin, pout, tin, R}) <= 0 || cp <= R ||
      eta <= 0 || eta > 1)
    throw std::invalid_argument("Invalid compressor P,T,Cp or efficiency");
  if (pout <= pin) return {tin, 0};
  double rise = tin * std::expm1((R / cp) * std::log(pout / pin)) / eta;
  return {tin + rise, cp * rise};
}
DynamicBufferNetwork::DynamicBufferNetwork(const Json& c, const std::vector<std::string>& names,
                                           const Vec& m, double gas)
    : cfg(c),
      comp_names(names),
      M(m),
      R(gas),
      n_comp(names.size()),
      block_size(n_comp + 1),
      size(0),
      enabled(false) {
  liquid=low(str(cfg,"FLUID_PHASE","gas"))=="liquid";
  liquid_density=num(cfg,"LIQUID_DENSITY",1000.);liquid_heat_capacity=num(cfg,"LIQUID_HEAT_CAPACITY",4180.);

  for (size_t i = 0; i < n_comp; ++i) comp_index[names[i]] = i;
  Json defs = cfg.value("DYNAMIC_BUFFERS", cfg.value("BUFFERS", Json::object()));
  if (defs.is_null()) defs = Json::object();
  if (!defs.is_object()) throw std::invalid_argument("DYNAMIC_BUFFERS must be a mapping");
  Vec defaultcp = scale(M, num(cfg, "Cpg", 1000));
  defaultcp =
      arr(cfg.value("BUFFER_COMPONENT_CP_MOLAR",
                    cfg.value("GAS_COMPONENT_CP_MOLAR", cfg.value("CPG_MOLAR", Json(defaultcp)))),
          n_comp);
  for (auto cp : defaultcp)
    if (!std::isfinite(cp) || cp <= (liquid?0.:R))
      throw std::invalid_argument("Buffer component Cp must exceed R");
  for (auto& [name, d] : defs.items()) {
    if (!d.is_object()) throw std::invalid_argument("Buffer definition must be a mapping");
    BufferDefinition b;
    b.name = name;
    b.volume = num(d, "volume_m3", num(d, "volume"));
    b.initial_pressure = num(d, "initial_pressure", num(cfg, "P_init", 1e5));
    b.initial_temperature = num(d, "initial_temperature", num(cfg, "T_init", 298.15));
    b.initial_composition = norm(d.value("initial_composition", cfg.at("y_init")), n_comp);
    b.initial_concentration=arr(d.value("initial_concentration",Json(Vec(n_comp,0.))),n_comp);
    b.maximum_volume=num(d,"maximum_volume",1e100);
    b.cp_molar = arr(d.value("cp_molar", Json(defaultcp)), n_comp);
    b.thermal_mode = low(str(d, "thermal_mode", "isothermal"));
    b.ua = num(d, "UA_W_K", num(d, "ua"));
    b.ambient_temperature = num(d, "ambient_temperature", b.initial_temperature);
    b.aliases = d.value("aliases", d.value("source_names", std::vector<std::string>{}));
    b.pressure_control = d.value("pressure_control", Json::object());
    if (b.pressure_control.is_null()) b.pressure_control = Json::object();
    if (!std::isfinite(b.volume) || b.volume <= 0 || !std::isfinite(b.initial_pressure) ||
        b.initial_pressure < 0 || !std::isfinite(b.initial_temperature) ||
        b.initial_temperature <= 0 || !std::isfinite(b.ua) || b.ua < 0 ||
        !std::isfinite(b.ambient_temperature) || b.ambient_temperature <= 0)
      throw std::invalid_argument("Invalid buffer volume, P,T or UA");
    for (auto cp : b.cp_molar)
      if (!std::isfinite(cp) || cp <= (liquid?0.:R)) throw std::invalid_argument("Buffer Cp must exceed R");
    if (b.thermal_mode != "isothermal" && b.thermal_mode != "adiabatic" &&
        b.thermal_mode != "nonadiabatic" && b.thermal_mode != "non_adiabatic")
      throw std::invalid_argument("Unknown buffer thermal mode");
    if (b.thermal_mode != "isothermal" && b.initial_pressure == 0)
      throw std::invalid_argument("Nonisothermal buffer needs positive initial P");
    auto ctl = b.pressure_control;
    if(liquid&&!ctl.empty())throw std::invalid_argument("Liquid buffers use volume controls, not ideal-gas pressure controllers");
    if (!ctl.empty()) {
      if (!ctl.contains("setpoint_pressure") && !ctl.contains("pressure"))
        throw std::invalid_argument("Controller setpoint required");
      double p = num(ctl, "setpoint_pressure", num(ctl, "pressure")),
             k = num(ctl, "conductance_mol_s_pa"), dead = num(ctl, "deadband_pa"),
             lim = num(ctl, "maximum_molar_flow", INFINITY);
      if (!std::isfinite(p) || p <= 0 || !std::isfinite(k) || k < 0 || !std::isfinite(dead) ||
          dead < 0 || std::isnan(lim) || lim <= 0)
        throw std::invalid_argument("Invalid buffer controller");
    }
    auto aliases = b.aliases;
    aliases.insert(aliases.begin(), name);
    for (auto& a : aliases) {
      if (alias_to_name.count(a)) throw std::invalid_argument("Duplicate buffer alias");
      alias_to_name[a] = name;
    }
    name_to_index[name] = definitions.size();
    definitions.push_back(b);
  }
  if(liquid)block_size=n_comp+2;
  size = definitions.size() * block_size;
  enabled = !definitions.empty();
}
Vec DynamicBufferNetwork::initial_state() const {
  Vec s;
  for (auto& d : definitions) {
    auto n = scale(d.initial_composition,
                   d.initial_pressure * d.volume / std::max(R * d.initial_temperature, 1e-300));
    if(liquid)n=scale(d.initial_concentration,d.volume);
    s.insert(s.end(), n.begin(), n.end());
    s.push_back(d.initial_temperature);
    if(liquid)s.push_back(d.volume);
  }
  return s;
}
Json DynamicBufferNetwork::state_record(const std::string& name, const Vec& s) const {
  const auto& d = definitions.at(name_to_index.at(name));
  size_t start = name_to_index.at(name) * block_size;
  if (s.size() != size) throw std::invalid_argument("Buffer state size mismatch");
  Vec n(n_comp);
  for (size_t i = 0; i < n_comp; ++i) n[i] = std::max(s[start + i], 0.);
  double T = std::max(s[start + n_comp], 1.), total = std::accumulate(n.begin(), n.end(), 0.);
  if(liquid) {
    double V=s[start+n_comp+1];
    if(!std::isfinite(V)||V<=0.||V>d.maximum_volume)throw std::runtime_error("Liquid buffer volume is outside its configured bounds");
    return {{"name",name},{"kind","buffer"},{"pressure",d.initial_pressure},{"temperature",T},{"volume",V},{"component_moles",n},{"total_moles",total},{"concentration",scale(n,1./V)},{"composition",total>1e-30?scale(n,1./total):d.initial_composition},{"energy_J",V*liquid_density*liquid_heat_capacity*T}};
  }
  Vec y = total > 1e-30 ? scale(n, 1 / total) : d.initial_composition, cv = d.cp_molar;
  for (auto& x : cv) x -= R;
  return {{"kind", "buffer"},
          {"buffer", name},
          {"component_moles", n},
          {"total_moles", total},
          {"composition", y},
          {"temperature", T},
          {"pressure", total * R * T / d.volume},
          {"volume_m3", d.volume},
          {"energy_J", dot(n, cv) * T},
          {"cp_molar", d.cp_molar}};
}
Json DynamicBufferNetwork::records(const Vec& s) const {
  Json j = Json::object();
  for (auto& d : definitions) j[d.name] = state_record(d.name, s);
  return j;
}
Json DynamicBufferNetwork::source_states(const Vec& s, double t, const Json& ctx) const {
  Json j = configured_sources(cfg, n_comp, t, ctx);
  for (auto& d : definitions) {
    auto rec = state_record(d.name, s);
    j[d.name] = rec;
    for (auto& a : d.aliases) j[a] = rec;
  }
  return j;
}
std::optional<std::string> DynamicBufferNetwork::buffer_name_for_port(const Json& p) const {
  for (auto k : {"buffer", "buffer_name"})
    if (p.contains(k) && !p[k].is_null() && p[k] != "") {
      auto n = p[k].get<std::string>();
      if (!name_to_index.count(n)) throw std::invalid_argument("Unknown dynamic buffer " + n);
      return n;
    }
  for (auto k : {"source", "sink"})
    if (p.contains(k) && p[k].is_string()) {
      auto it = alias_to_name.find(p[k].get<std::string>());
      if (it != alias_to_name.end()) return it->second;
    }
  return {};
}
std::vector<std::string> DynamicBufferNetwork::controller_tags() const {
  std::vector<std::string> t;
  for (auto& d : definitions)
    for (auto k : {"export_tag", "makeup_tag"})
      if (d.pressure_control.contains(k) && !d.pressure_control[k].is_null()) {
        auto s = d.pressure_control[k].get<std::string>();
        if (std::find(t.begin(), t.end(), s) == t.end()) t.push_back(s);
      }
  return t;
}
std::pair<Vec, Json> DynamicBufferNetwork::conservative_derivatives(const Vec& s, const Json& rates,
                                                                    const Json& hrates,
                                                                    const Json& sources, const Json& volume_rates) const {
  Vec out, ein(n_comp, 0), eout(n_comp, 0);
  double hin = 0, hout = 0, heat = 0, cool = 0;
  Json streams = Json::object(), prov = Json::array();
  for (auto& d : definitions) {
    auto rec = state_record(d.name, s);
    Vec dn = arr(rates.value(d.name, Json(Vec(n_comp, 0))), n_comp);
    double power = hrates.value(d.name, 0.);
    const auto& ctl = d.pressure_control;
    double T = rec["temperature"].get<double>(), p = rec["pressure"].get<double>();
    if(liquid&&!ctl.empty())throw std::invalid_argument("Liquid buffers use volume controls, not ideal-gas pressure controllers");
    if (!ctl.empty()) {
      double target = num(ctl, "setpoint_pressure", num(ctl, "pressure")), diff = p - target,
             dead = num(ctl, "deadband_pa"), gain = num(ctl, "conductance_mol_s_pa"),
             lim = num(ctl, "maximum_molar_flow", INFINITY);
      if (diff > dead) {
        double flow = std::min(
            gain * std::min(diff - dead, std::max(p - num(ctl, "export_pressure"), 0.)), lim);
        auto amounts = scale(rec["composition"].get<Vec>(), flow);
        double H = dot(amounts, d.cp_molar) * T;
        add(dn, amounts, -1);
        power -= H;
        hout += H;
        auto tag = str(ctl, "export_tag", "waste");
        addj(streams, tag, amounts);
        add(eout, amounts);
        prov.push_back({{"tag", tag},
                        {"direction", "out"},
                        {"source", d.name},
                        {"component_rates", amounts},
                        {"origin", "external"}});
      } else if (diff < -dead && ctl.contains("makeup_source") && !ctl["makeup_source"].is_null() &&
                 ctl["makeup_source"] != "") {
        auto source_name = ctl["makeup_source"].get<std::string>();
        const auto& source = sources.at(source_name);
        if (source.value("kind", "") == "buffer")
          throw std::invalid_argument(
              "Connect finite buffers explicitly instead of controller makeup");
        double avail = std::max(source["pressure"].get<double>() - p, 0.),
               flow = std::min(gain * std::min(-diff - dead, avail), lim);
        auto amounts = scale(source["composition"].get<Vec>(), flow);
        double H = dot(amounts, d.cp_molar) * source["temperature"].get<double>();
        add(dn, amounts);
        power += H;
        hin += H;
        auto tag = str(ctl, "makeup_tag", source_name);
        addj(streams, tag, amounts);
        add(ein, amounts);
        prov.push_back({{"tag", tag},
                        {"direction", "in"},
                        {"source", source_name},
                        {"component_rates", amounts},
                        {"origin", "external"}});
      }
    }
    double q = (d.thermal_mode == "nonadiabatic" || d.thermal_mode == "non_adiabatic")
                   ? d.ua * (d.ambient_temperature - T)
                   : 0;
    if(liquid) {
      double dv=num(volume_rates,d.name.c_str());
      double V=rec.at("volume"), capacity=liquid_density*liquid_heat_capacity;
      double dt=(power+q-capacity*T*dv)/(capacity*V);
      if(d.thermal_mode=="isothermal"){q=capacity*T*dv-power;dt=0.;}
      out.insert(out.end(),dn.begin(),dn.end());out.push_back(dt);out.push_back(dv);
      heat+=std::max(q,0.);cool+=std::max(-q,0.);continue;
    }
    Vec cv = d.cp_molar;
    for (auto& v : cv) v -= R;
    double dt;
    if (d.thermal_mode == "isothermal") {
      q = T * dot(cv, dn) - power;
      dt = 0;
    } else
      dt = (power + q - T * dot(cv, dn)) /
           std::max(dot(rec["component_moles"].get<Vec>(), cv), 1e-12);
    heat += std::max(q, 0.);
    cool += std::max(-q, 0.);
    out.insert(out.end(), dn.begin(), dn.end());
    out.push_back(dt);
  }
  return {out,
          {{"external_in", ein},
           {"external_out", eout},
           {"external_stream_rates", streams},
           {"provenance", prov},
           {"utility_rates_W", {{"buffer_heating_W", heat}, {"buffer_cooling_W", cool}}},
           {"external_enthalpy_in_W", hin},
           {"external_enthalpy_out_W", hout}}};
}
double DynamicBufferNetwork::metric(const std::string& metric, const Vec& s,
                                    const std::string& name, const Json& component) const {
  auto r = state_record(name, s);
  auto k = low(metric);
  if(k=="buffer_volume")return liquid?r["volume"]:r["volume_m3"];
  if (k == "buffer_pressure") return r["pressure"];
  if (k == "buffer_temperature") return r["temperature"];
  if (k == "buffer_total_moles") return r["total_moles"];
  if (component.is_null()) throw std::invalid_argument("Buffer metric requires component");
  size_t i =
      component.is_string() ? comp_index.at(component.get<std::string>()) : component.get<size_t>();
  if(k=="buffer_concentration")return liquid?r["concentration"].at(i).get<double>():r["component_moles"].at(i).get<double>()/r["volume_m3"].get<double>();
  if (k == "buffer_component_moles") return r["component_moles"].at(i);
  if (k == "buffer_mole_fraction") return r["composition"].at(i);
  throw std::invalid_argument("Unsupported buffer metric");
}
}  // namespace cycle

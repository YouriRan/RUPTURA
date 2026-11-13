#include "cycle/equilibrium.hpp"
#include "cycle/main_models.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
namespace cycle {
namespace {
constexpr double tiny = 1e-300, ceiling = 1e300;
double clip(double x, double lo = 0., double hi = ceiling) {
  return std::isnan(x) ? lo : std::clamp(x, lo, hi);
}
double val(const Json& j, const std::string& k, double d) {
  if (!j.is_object() || !j.contains(k) || !j[k].is_number()) return d;
  double v = j[k].get<double>();
  return std::isfinite(v) ? v : d;
}
std::string lower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
  auto a = s.find_first_not_of(" \r\n\t"), b = s.find_last_not_of(" \r\n\t");
  return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}
std::string upper(std::string s) {
  s = lower(s);
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::toupper(c); });
  return s;
}
bool flag(const Json& j, const std::string& k, bool d) {
  if (!j.contains(k)) return d;
  if (j[k].is_boolean()) return j[k].get<bool>();
  if (j[k].is_number()) return j[k].get<double>() != 0.;
  if (j[k].is_string()) {
    auto s = lower(j[k].get<std::string>());
    return s == "true" || s == "yes" || s == "on" || s == "1";
  }
  return d;
}
double exp_safe(double x, double cap = 700.) { return std::exp(clip(x, -cap, cap)); }
double multiply(double a, double b, double cap = ceiling) {
  a = clip(a, 0., cap);
  b = clip(b, 0., cap);
  return a > 0. && b > 0. ? std::exp(std::min(std::log(a) + std::log(b), std::log(cap))) : 0.;
}
double power(double b, double e, double r = 0., bool local_c1=false) {
  b = clip(b);
  r = clip(r);
  e = clip(e, 0., 1e6);
  if (b <= 0.) return 0.;
  if(local_c1 && r>0. && e<1.) {
    if(b>=r)return std::exp(clip(e*std::log(b),-700.,std::log(ceiling)));
    return std::exp(clip(std::log(b)+(e-1.)*std::log(r),-700.,std::log(ceiling)))*((2.-e)+(e-1.)*b/r);
  }
  if (r <= 0.) return std::exp(clip(e * std::log(b), -700., std::log(ceiling)));
  if (e <= 0.) return 0.;
  double delta = e * std::log1p(b / r), f = delta > 50. ? delta : std::log(std::expm1(delta));
  return std::exp(clip(e * std::log(r) + f, -700., std::log(ceiling)));
}
double bound_ratio(double q, double n, double d) {
  return clip(q * n / std::max(d, tiny), 0., std::max(q, 0.));
}
bool powered(const std::string& m) {
  return m == "freundlich" || m == "langmuir_freundlich" || m == "sips" || m == "toth" ||
         m == "modified_redlich";
}
bool rh_model(const std::string& m) {
  return m == "humidity_enhanced_langmuir" || m == "moisture_swing_langmuir" ||
         m == "gab_coadsorption_langmuir" || m == "gab_coadsorption_toth" ||
         m == "gab_coadsorption_sips" || m == "gab_modified_toth_liu" || m == "gab_water" || m == "gab_main" ||
         m == "humidity_toth_climeworks";
}
std::string model_name(std::string m) {
  m = lower(m);
  std::replace(m.begin(), m.end(), '-', '_');
  static const std::map<std::string, std::string> aliases = {
      {"none", "zero"},
      {"nonadsorbing", "zero"},
      {"carrier", "zero"},
      {"lf", "langmuir_freundlich"},
      {"antilangmuir", "anti_langmuir"},
      {"modified redlich", "modified_redlich"},
      {"redlich", "modified_redlich"},
      {"wang_lackner_moisture_swing", "moisture_swing_langmuir"},
      {"gab", "gab_water"},
      {"water_gab", "gab_water"},
      {"liu_modified_toth", "gab_modified_toth_liu"},
      {"rh_dependent_toth", "humidity_toth_climeworks"},
      {"ph_dependent_langmuir", "ph_langmuir"}};
  auto it = aliases.find(m);
  if (it != aliases.end()) m = it->second;
  if (m != "zero" && m != "henry" && m != "freundlich" && m != "langmuir" &&
      m != "langmuir_freundlich" && m != "sips" && m != "anti_langmuir" && m != "toth" &&
      m != "modified_redlich" && m != "ph_langmuir" && !rh_model(m))
    throw std::invalid_argument("Unknown isotherm model: " + m);
  return m;
}
std::string iast_method(std::string s) {
  s = lower(s);
  s.erase(
      std::remove_if(s.begin(), s.end(), [](char c) { return c == '-' || c == '_' || c == ' '; }),
      s.end());
  if (s == "auto") return "AUTO";
  if (s == "fast" || s == "fastias" || s == "fastiast" || s == "newton") return "FASTIAST";
  if (s == "bisection" || s == "nestedloopbisection" || s == "nestedbisection") return "BISECTION";
  throw std::invalid_argument("IAST_METHOD must be AUTO, FASTIAST or BISECTION");
}
double array_component(const Json& c, const std::string& key, std::size_t j, double d = 0.) {
  if (!c.contains(key)) return d;
  auto a = c[key];
  if (a.is_number()) return a.get<double>();
  if (a.is_array() && !a.empty() && a[0].is_array()) a = a[0];
  if (!a.is_array() || j >= a.size()) throw std::invalid_argument(key + " component shape");
  return a[j].get<double>();
}
Vec dense_solve(Vec a, Vec b) {
  std::size_t n = b.size();
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t p = k;
    for (std::size_t i = k + 1; i < n; ++i)
      if (std::abs(a[i * n + k]) > std::abs(a[p * n + k])) p = i;
    if (std::abs(a[p * n + k]) <= tiny) throw std::runtime_error("FastIAST singular Jacobian");
    if (p != k) {
      for (std::size_t j = k; j < n; ++j) std::swap(a[p * n + j], a[k * n + j]);
      std::swap(b[p], b[k]);
    }
    for (std::size_t i = k + 1; i < n; ++i) {
      double f = a[i * n + k] / a[k * n + k];
      a[i * n + k] = 0.;
      for (std::size_t j = k + 1; j < n; ++j) a[i * n + j] -= f * a[k * n + j];
      b[i] -= f * b[k];
    }
  }
  for (int i = int(n) - 1; i >= 0; --i) {
    for (std::size_t j = i + 1; j < n; ++j) b[i] -= a[i * n + j] * b[j];
    b[i] /= a[i * n + i];
    if (!std::isfinite(b[i])) throw std::runtime_error("FastIAST nonfinite correction");
  }
  return b;
}
struct Site {
  std::string model = "langmuir", numode = "none";
  double q = 0., k = 0., heat = 0., nu = 1., anu = 0., cc = 0., ctref = 1.;
  bool inverse = false, concentration_affinity = false;
  Json extra = Json::object();
};
struct LocalSite {
  const Site* ref;
  double q, k, nu;
};
const std::array<std::pair<double, double>, 48>& quadrature() {
  static const auto a = [] {
    std::array<std::pair<double, double>, 48> out{};
    constexpr int n = 48;
    for (int i = 0; i < n / 2; ++i) {
      double x = std::cos(3.14159265358979323846 * (i + .75) / (n + .5)), pp = 0.;
      for (int k = 0; k < 30; ++k) {
        double p1 = 1., p2 = 0.;
        for (int j = 1; j <= n; ++j) {
          double p3 = p2;
          p2 = p1;
          p1 = ((2 * j - 1) * x * p2 - (j - 1) * p3) / j;
        }
        pp = n * (x * p1 - p2) / (x * x - 1.);
        double next = x - p1 / pp;
        if (std::abs(next - x) < 2e-16) {
          x = next;
          break;
        }
        x = next;
      }
      double w = 2. / ((1. - x * x) * pp * pp);
      out[i] = {-x, w};
      out[n - 1 - i] = {x, w};
    }
    return out;
  }();
  return a;
}
double basic_loading(const LocalSite& a, double p, double reg) {
  auto& m = a.ref->model;
  double q = clip(a.q), k = clip(a.k), n = std::max(clip(a.nu, 0., 1e6), tiny);
  p = clip(p);
  if (m == "zero" || q <= 0. || p <= 0.) return 0.;
  double kp = multiply(k, p), x = 0.;
  if (m == "henry") return clip(q * kp);
  if (m == "freundlich") return clip(q * multiply(k, power(p, n, reg,a.ref->extra.value("local_c1",false))));
  if (m == "langmuir"||m=="ph_langmuir") return bound_ratio(q, kp, 1. + kp);
  if (m == "langmuir_freundlich") {
    x = multiply(k, power(p, n, reg));
    return bound_ratio(q, x, 1. + x);
  }
  if (m == "sips") {
    x = power(kp, n, multiply(k, reg));
    return bound_ratio(q, x, 1. + x);
  }
  if (m == "anti_langmuir") return bound_ratio(q, kp, std::max(1. - kp, 1e-12));
  if (m == "toth")
    return bound_ratio(q, kp, std::max(power(1. + power(kp, n, multiply(k, reg)), 1. / n), 1e-12));
  if (m == "modified_redlich") return bound_ratio(q, kp, 1. + power(kp, n, multiply(k, reg)));
  throw std::invalid_argument("Not a standard pure isotherm: " + m);
}
double numerical_psi(const LocalSite& a, double p, double reg) {
  if (p <= 0.) return 0.;
  double result = 0.;
  if (reg > 0. && powered(a.ref->model)) {
    double ratio = p / reg,
           umax = std::isfinite(ratio) ? std::log1p(ratio) : std::log(p) - std::log(reg);
    for (const auto& [x, w] : quadrature()) {
      double u = .5 * umax * (x + 1.), d = reg * std::expm1(u);
      result += .5 * umax * w * basic_loading(a, d, reg) / (-std::expm1(-u));
    }
  } else
    for (const auto& [x, w] : quadrature()) {
      double d = .5 * p * (x + 1.);
      result += .5 * p * w * basic_loading(a, d, 0.) / std::max(d, tiny);
    }
  return clip(result);
}
double spreading(const LocalSite& a, double p, double reg) {
  const auto& m = a.ref->model;
  double q = clip(a.q), k = clip(a.k), n = std::max(clip(a.nu, 0., 1e6), tiny);
  p = clip(p);
  if (m == "zero" || q <= 0. || k <= 0. || p <= 0.) return 0.;
  if (reg > 0. && powered(m)) return numerical_psi(a, p, reg);
  double kp = clip(k * p);
  if (m == "henry") return clip(q * kp);
  if (m == "freundlich") return clip(q * k / n * power(p, n));
  if (m == "langmuir") return clip(q * std::log1p(kp));
  if (m == "langmuir_freundlich") return clip(q / n * std::log1p(k * power(p, n)));
  if (m == "sips") return clip(q / n * std::log1p(power(kp, n)));
  if (m == "toth") {
    double theta = clip(kp / std::max(power(1. + power(kp, n), 1. / n), tiny), 0., 1. - 1e-15),
           tn = power(theta, n),
           result = q * (theta - theta / n * std::log(std::max(1. - tn, tiny))), term = q * theta,
           acc = 0.;
    for (int i = 0; i < 100; ++i) {
      term *= tn;
      acc += n;
      result -= term / std::max(acc * (acc + 1.), tiny);
    }
    return clip(result);
  }
  return numerical_psi(a, p, 0.);
}
double inverse_bisection(const std::function<double(double)>& psi, double target, double cached,
                         double cap, double tol, int iters) {
  double l = 0., r = clip(cached, 0., cap);
  if (r <= 0.) r = std::min(1., cap);
  int steps = 0;
  while (psi(r) < target) {
    if (r >= cap || ++steps > iters) throw std::runtime_error("IAST inverse bracket failed");
    r = std::min(2 * r, cap);
  }
  for (int i = 0; i < iters; ++i) {
    double mid = .5 * (l + r);
    if (psi(mid) > target)
      r = mid;
    else
      l = mid;
    if ((r - l) / std::max(r + l, tiny) <= tol) return .5 * (l + r);
  }
  throw std::runtime_error("IAST inverse failed to converge");
}
double inverse_site(const LocalSite& a, double psi, double cached, double cap, double tol,
                    int iters, double reg) {
  if (psi <= 0.) return 0.;
  const auto& m = a.ref->model;
  double q = a.q, k = a.k, n = a.nu;
  if (q <= 0. || k <= 0.) throw std::runtime_error("IAST inverse on zero site");
  if (!(reg > 0. && powered(m))) {
    if (m == "langmuir") return std::min(std::expm1(std::min(psi / q, 700.)) / k, cap);
    if (m == "henry") return std::min(psi / std::max(q * k, tiny), cap);
    if (m == "freundlich") return std::min(power(psi * n / std::max(q * k, tiny), 1. / n), cap);
    if (m == "langmuir_freundlich")
      return std::min(power(std::expm1(std::min(psi * n / q, 700.)) / k, 1. / n), cap);
    if (m == "sips")
      return std::min(power(std::expm1(std::min(psi * n / q, 700.)), 1. / n) / k, cap);
  }
  return inverse_bisection([&](double d) { return spreading(a, d, reg); }, psi, cached, cap, tol,
                           iters);
}
}  // namespace
double water_saturation_pressure(double T, const Json& c) {
  std::string m = lower(c.value("WATER_PSAT_MODEL", std::string("iapws")));
  if (m == "constant" || m == "manual") return c.at("WATER_PSAT_CONSTANT").get<double>();
  if (m == "antoine") {
    double t = clip(T, 293., 473.), a = t <= 343. ? 6.20963 : (t <= 373. ? 5.40221 : 5.20389),
           b = t <= 343. ? 2354.731 : (t <= 373. ? 1838.675 : 1733.926),
           cc = t <= 343. ? 7.559 : (t <= 373. ? -31.737 : -39.485);
    return std::pow(10., a - b / (t + cc)) * 1e5;
  }
  if (m != "iapws" && m != "wagner" && m != "iapws_wagner")
    throw std::invalid_argument("Unsupported water saturation model");
  double t = clip(T, 273.16, 647.096 - 1e-9), th = 1. - t / 647.096;
  return 22.064e6 * exp_safe(647.096 / t *
                             (-7.85951783 * th + 1.84408259 * std::pow(th, 1.5) -
                              11.7866497 * std::pow(th, 3) + 22.6807411 * std::pow(th, 3.5) -
                              15.9618719 * std::pow(th, 4) + 1.80122502 * std::pow(th, 7.5)));
}
double relative_humidity(double pw, double T, const Json& c) {
  auto m = lower(c.value("RH_MODE", c.value("WATER_RH_MODE", std::string("from_PH2O"))));
  double rh = 0.;
  if (m == "direct" || m == "fixed" || m == "constant")
    rh = val(c, "RH_VALUE", 0.);
  else if (m == "from_ph2o" || m == "from_partial_pressure" || m == "partial_pressure")
    rh = pw / std::max(water_saturation_pressure(T, c), tiny);
  else
    throw std::invalid_argument("Unsupported RH mode");
  if (std::isnan(rh)) rh = 0.;
  if (std::isinf(rh)) rh = rh > 0. ? val(c, "RH_MAX", 1.) : 0.;
  if (!flag(c, "RH_CLIP", true)) return rh;
  return rh > 0. ? clip(rh, val(c, "RH_MIN", 1e-6), val(c, "RH_MAX", 1.)) : 0.;
}

struct MaterialEquilibrium::Impl {
  Json cfg;
  double R, reg, driving_cap, qcap, tol, active_tol;
  std::size_t nc, nz, ns;
  int iters;
  std::string mode, method;
  std::vector<std::string> names;
  std::vector<Site> sites;
  Vec capacities;
  std::vector<bool> bounded;
  bool needs_rh = false, cache_enabled = true;
  int water = -1;
  mutable Vec cache_p, cache_t, cache_q;
  mutable std::size_t hits = 0, misses = 0;
  std::optional<MainMaterialModels> main_model;
  Impl(const Json& c, double gas, std::size_t n, std::size_t z, std::string phase)
      : cfg(c), R(gas), nc(n), nz(z) {
    if (nc < 1 || nz < 1 || !std::isfinite(R) || R <= 0.)
      throw std::invalid_argument("Positive counts and gas constant required");
    phase = lower(phase);
    if (phase != "phy" && phase != "chem") throw std::invalid_argument("phase must be phy or chem");
    names = cfg.at("comp_names").get<std::vector<std::string>>();
    if (names.size() != nc) throw std::invalid_argument("Component count mismatch");
    if(phase=="phy" && cfg.contains("MAIN_EQUILIBRIUM")) {
      main_model.emplace(cfg.at("MAIN_EQUILIBRIUM"));capacities.assign(nc,0.);bounded.assign(nc,false);
      for(std::size_t j=0;j<nc;++j)for(const auto& site:main_model->components[j].isotherm.sites)
        if(!site.parameters.empty())capacities[j]+=std::max(site.parameters[0],0.);
      if (auto maximum = main_model->mpd_maximum_loadings()) {
        capacities = std::move(*maximum);
        bounded.assign(nc, true);
      }
      return;
    }
    reg = std::max(val(cfg, "ISOTHERM_POWER_REG", val(cfg, "P_ISOTHERM_POWER_REG", 1e-8)), 0.);
    driving_cap = val(cfg, "ISOTHERM_DRIVING_CEILING",
                      val(cfg, "P_ISOTHERM_CEILING", val(cfg, "P_CEILING", ceiling)));
    qcap = val(cfg, "ISOTHERM_QEQ_CEILING", val(cfg, "Q_STATE_CEILING", ceiling));
    tol = std::max(val(cfg, "IAST_REL_TOL", 1e-12), 1e-15);
    active_tol = val(cfg, "MIXTURE_ACTIVE_TOL", 1e-15);
    iters = std::max(int(val(cfg, "IAST_MAX_ITER", 250)), 10);
    cache_enabled = flag(cfg, "ENABLE_EQUILIBRIUM_CACHE", true);
    Json entry;
    auto key = upper(phase) + "_ISOTHERM_INPUT";
    if (cfg.contains(key) && !cfg[key].is_null()) {
      if (cfg[key].size() != 1)
        throw std::invalid_argument("One material per equilibrium instance required");
      entry = cfg[key][0];
      mode = upper(entry.value("mixture_mode", std::string("")));
      ns = entry.value("n_sites", 0);
    } else {
      mode = upper(cfg.value(phase + "_MIXTURE_MODE", std::string("SCI")));
      ns = std::size_t(array_component(cfg, "N_sites_" + phase, 0, 1.));
    }
    if (mode != "SPI" && mode != "SCI" && mode != "SEI" && mode != "IAST" && mode != "SIAST")
      throw std::invalid_argument("Unsupported mixture method");
    if (ns < 1) throw std::invalid_argument("At least one isotherm site required");
    method = iast_method(cfg.value("IAST_METHOD", std::string("AUTO")));
    sites.resize(ns * nc);
    capacities.assign(nc, 0.);
    bounded.assign(nc, true);
    auto multi = [&](const std::string& k, std::size_t s, std::size_t j, double d) {
      if (!cfg.contains(k)) return d;
      const auto& a = cfg.at(k).at(0);
      if (a.size() == ns && a.at(0).size() == nc) return a.at(s).at(j).get<double>();
      if (a.size() == nc && a.at(0).size() == ns) return a.at(j).at(s).get<double>();
      throw std::invalid_argument("Invalid multisite tensor: " + k);
    };
    for (std::size_t s = 0; s < ns; ++s)
      for (std::size_t j = 0; j < nc; ++j) {
        auto& a = sites[s * nc + j];
        Json raw = Json::object();
        std::string nm = "none";
        double qb = array_component(cfg, "qmax_" + phase, j),
               kb = array_component(cfg, "kj0_" + phase, j),
               hb = array_component(cfg, phase == "chem" ? "dH_chem_eq" : "dH_phy", j);
        if (!entry.is_null()) {
          if (mode == "SCI") {
            a.model = model_name(entry.value("model", std::string("langmuir")));
            raw = entry.at("sites").at(s).at(names[j]);
            nm = entry.at("sites").at(s).value("nu_mode", std::string("none"));
          } else {
            const auto& component = entry.at("components").at(names[j]);
            a.model = model_name(component.value("model", std::string("langmuir")));
            raw = component.at("sites").at(s);
            nm = raw.value("nu_mode", std::string("none"));
          }
          a.q = val(raw, "q_scale", val(raw, "qsat", qb / double(ns)));
          a.k = val(raw, "k0", kb);
          a.heat = val(raw, "dH", hb);
        } else {
          if (cfg.contains(phase + "_ISOTHERM_MODEL")) {
            const auto& m = cfg.at(phase + "_ISOTHERM_MODEL");
            a.model = model_name(mode == "SCI" ? m.at(0).get<std::string>()
                                               : m.at(0).at(j).get<std::string>());
          }
          a.q = multi("qmax_" + phase + "_multi", s, j, qb);
          a.k = multi("kj0_" + phase + "_multi", s, j, kb);
          a.heat =
              multi((phase == "chem" ? "dH_chem_eq" : "dH_phy") + std::string("_multi"), s, j, hb);
          raw["nu0"] = multi("nu0_" + phase + "_multi", s, j, 1.);
          raw["A_nu"] = multi("A_nu_" + phase + "_multi", s, j, 0.);
          if (cfg.contains("nu_temperature_mode_" + phase))
            nm = cfg.at("nu_temperature_mode_" + phase).at(0).at(s).get<std::string>();
        }
        if (raw.contains("nu0") && raw.contains("freundlich_n"))
          throw std::invalid_argument("Use nu0 or freundlich_n, not both");
        a.nu = val(raw, "freundlich_n", val(raw, "nu0", 1.));
        a.anu = val(raw, "A_nu", 0.);
        a.numode = lower(nm);
        a.extra = raw;
        const std::string affinity_basis = lower(raw.value("affinity_basis", std::string("pressure")));
        if (affinity_basis != "pressure" && affinity_basis != "concentration")
          throw std::invalid_argument("affinity_basis must be pressure or concentration");
        a.concentration_affinity = affinity_basis == "concentration";
        if (a.concentration_affinity && a.model != "langmuir" && a.model != "zero")
          throw std::invalid_argument("Concentration affinity currently requires Langmuir/zero");
        std::string cv =
            lower(raw.value("exponent_convention",
                            raw.value("nu_convention", cfg.value("ISOTHERM_EXPONENT_CONVENTION",
                                                                 std::string("direct")))));
        if (raw.contains("freundlich_n") && !raw.contains("exponent_convention") &&
            !raw.contains("nu_convention"))
          cv = "inverse_n";
        a.inverse = cv == "inverse_n" || cv == "inverse-n" || cv == "reciprocal" ||
                    cv == "one_over_n" || cv == "1/n" || cv == "pysorb";
        if (!a.inverse && cv != "direct" && cv != "actual" && cv != "exponent" && cv != "nu")
          throw std::invalid_argument("Invalid exponent convention");
        if (raw.contains("capacity_temperature") && !raw["capacity_temperature"].is_null()) {
          auto& law = raw["capacity_temperature"];
          if (lower(law.value("model", std::string("exponential"))) != "exponential")
            throw std::invalid_argument("Capacity temperature must be exponential");
          a.ctref = law.at("reference_temperature_K");
          a.cc = law.at("coefficient");
          if (!std::isfinite(a.ctref) || a.ctref <= 0. || !std::isfinite(a.cc))
            throw std::invalid_argument("Invalid capacity temperature parameters");
        }
        if (a.model == "ph_langmuir" && lower(cfg.value("FLUID_PHASE",std::string("gas")))!="liquid")
          throw std::invalid_argument("pH-dependent isotherms require liquid operation");
        if (rh_model(a.model)) {
          if (mode != "SPI") throw std::invalid_argument("RH-dependent models require SPI");
          needs_rh = true;
        }
        if (mode == "SEI" && a.model != "langmuir" && a.model != "zero")
          throw std::invalid_argument("SEI requires Langmuir/zero");
        if (mode == "SCI" && (a.model == "zero" || a.model == "henry" || a.model == "freundlich"))
          throw std::invalid_argument("Unsupported SCI model");
        if (a.model.rfind("gab_", 0) == 0) {
          double qm = val(raw, "qm_w", val(raw, "q_m", 0.)),
                 qblock = val(raw, "qH2O_block_max", val(raw, "q_h2o_block_max", std::max(qm, 0.)));
          if (qm < 0. || qblock < 0. || val(raw, "Kw", val(raw, "K_w", 0.)) < 0.)
            throw std::invalid_argument("GAB capacities and Kw must be nonnegative");
          if (a.model == "gab_modified_toth_liu" &&
              (val(raw, "gamma", 0.) < 0. || val(raw, "beta", 0.) < 0. ||
               val(raw, "qsat_cap", a.q) < a.q))
            throw std::invalid_argument("Invalid Liu parameters");
        }
        if (a.model == "humidity_toth_climeworks") {
          double qi = val(raw, "qinf0", val(raw, "q_inf_0", a.q));
          if (qi < 0. || val(raw, "b0", 0.) < 0. || val(raw, "A", 0.) < 0.)
            throw std::invalid_argument("Invalid Climeworks parameters");
          if (a.q > 0. && std::abs(qi - a.q) > val(cfg, "ISOTHERM_SPECIAL_PARAM_RTOL", 1e-10) *
                                                   std::max({std::abs(a.q), std::abs(qi), 1.}))
            throw std::invalid_argument("Climeworks qsat must equal qinf0");
        }
        bool finite = a.model != "henry" && a.model != "freundlich" && a.model != "anti_langmuir" &&
                      a.model != "modified_redlich";
        bounded[j] = finite;
        double capacity = a.q;
        if (a.model == "gab_modified_toth_liu")
          capacity = std::max(capacity, val(raw, "qsat_cap", capacity));
        if (a.model == "humidity_toth_climeworks")
          capacity = val(raw, "qinf0", val(raw, "q_inf_0", capacity));
        if (finite)
          capacities[j] += std::max(capacity, 0.);
        else
          capacities[j] = std::numeric_limits<double>::infinity();
      }
    if (mode != "SCI")
      for (std::size_t s = 0; s < ns; ++s) {
        std::string nm = "none";
        for (std::size_t j = 0; j < nc; ++j) {
          auto candidate = sites[s * nc + j].numode;
          if (nm == "none")
            nm = candidate;
          else if (candidate != "none" && candidate != nm)
            throw std::invalid_argument("Inconsistent nu_mode across components");
        }
        for (std::size_t j = 0; j < nc; ++j) sites[s * nc + j].numode = nm;
      }
    if (needs_rh) {
      auto wm = cfg.value("WATER_COMPONENT", std::string("H2O"));
      auto it = std::find(names.begin(), names.end(), wm);
      if (it != names.end()) water = int(it - names.begin());
      auto rm = lower(cfg.value("RH_MODE", cfg.value("WATER_RH_MODE", std::string("from_PH2O"))));
      if (water < 0 && rm != "direct" && rm != "fixed" && rm != "constant")
        throw std::invalid_argument("RH model requires water component or direct RH");
    }
  }
  std::vector<LocalSite> local(double temperature) const {
    double t = clip(temperature,
                    val(cfg, "ISOTHERM_T_FLOOR", std::max(val(cfg, "T_FLOOR", 150.), 1e-12)),
                    val(cfg, "ISOTHERM_T_CEILING", val(cfg, "T_CEILING", 2000.))),
           ec = val(cfg, "ISOTHERM_EXP_ARG_CLIP", 700.);
    std::vector<LocalSite> out;
    out.reserve(sites.size());
    for (const auto& a : sites) {
      double k = clip(a.k * exp_safe(a.heat / (R * t), ec), 0.,
                      val(cfg, "ISOTHERM_K_CEILING", ceiling)),
             n = a.nu;
      // Concentration-basis k0 has units m3/mol; the solver passes partial
      // pressures. Retain the full 1/(R*T) conversion at every temperature.
      if (a.concentration_affinity) k /= R * t;
      if(a.model=="ph_langmuir")k/=1.+std::pow(10.,val(a.extra,"pH0",7.)-val(cfg,"PH_VALUE",7.));
      if (a.numode == "a_over_t")
        n += a.anu / t;
      else if (a.numode == "a_times_t")
        n += a.anu * t;
      else if (a.numode != "none" && a.numode != "" && a.numode != "constant")
        throw std::invalid_argument("Unsupported nu temperature mode");
      if (a.inverse) {
        if (!std::isfinite(n) || n <= 0.)
          throw std::invalid_argument("inverse_n requires positive parameter");
        n = 1. / n;
      }
      n = clip(n, std::max(val(cfg, "ISOTHERM_NU_FLOOR", 1e-12), tiny),
               val(cfg, "ISOTHERM_NU_CEILING", 1e6));
      out.push_back({&a, clip(a.q * exp_safe(a.cc * (t / a.ctref - 1.), ec), 0., qcap), k, n});
    }
    return out;
  }
  double gab(double rh, double qm, double cw, double kw) const {
    kw = clip(kw, val(cfg, "GAB_KW_MIN", 0.), val(cfg, "GAB_KW_MAX", .999999));
    cw = clip(cw, val(cfg, "GAB_CW_MIN", 0.), val(cfg, "GAB_CW_MAX", 1e6));
    double krh = kw * rh, f = val(cfg, "GAB_DEN_FLOOR", 1e-12);
    return clip(
        std::max(qm, 0.) * cw * krh / (std::max(1. - krh, f) * std::max(1. - krh + cw * krh, f)),
        0., val(cfg, "GAB_QH2O_CEILING", 1e100));
  }
  double pure(const LocalSite& a, double p, double t, double rh) const;
  Vec sci(const std::vector<LocalSite>& a, const Vec& p) const;
  Vec sei(const std::vector<LocalSite>& a, const Vec& p) const;
  Vec fast_iast(const std::vector<LocalSite>& a, const Vec& p, int site, const Vec& warm,
                Vec& pstar) const;
  Vec iast(const std::vector<LocalSite>& a, const Vec& p, int site, const Vec& warm,
           Vec& pstar) const;
  Vec loading(const Vec& p, const Vec& t) const;
};
double MaterialEquilibrium::Impl::pure(const LocalSite& a, double p, double t, double rh) const {
  const auto& m = a.ref->model;
  const auto& e = a.ref->extra;
  if (!rh_model(m)) return basic_loading(a, p, reg);
  rh = rh > 0. ? clip(rh, val(cfg, "RH_MIN", 1e-6), val(cfg, "RH_MAX", 1.)) : 0.;
  double q = a.q, k = a.k, n = a.nu;
  auto langmuir = [&](double capacity, double affinity) {
    double x = multiply(affinity, p);
    return bound_ratio(capacity, x, 1. + x);
  };
  if (m == "humidity_enhanced_langmuir") {
    double kw = std::max(val(e, "kw", val(e, "K_w", 0.)), 0.),
           alpha = std::max(val(e, "alpha", 0.), 0.);
    return langmuir(q, clip(k * (1. + alpha * kw * rh / std::max(1. + kw * rh, tiny)), 0.,
                            val(cfg, "ISOTHERM_K_CEILING", ceiling)));
  }
  if (m == "moisture_swing_langmuir") {
    double nms = clip(val(e, "n", n), val(cfg, "MOISTURE_SWING_N_MIN", 0.),
                      val(cfg, "MOISTURE_SWING_N_MAX", 20.));
    return langmuir(
        q, clip(k * exp_safe(-nms * std::log(rh), val(cfg, "ISOTHERM_EXP_ARG_CLIP", 700.)), 0.,
                val(cfg, "ISOTHERM_K_CEILING", ceiling)));
  }
  if(m=="gab_main") {
    double C=val(e,"C0",1.)*exp_safe(val(e,"dHC",0.)/(R*t)),K=val(e,"K0",.5)*exp_safe(val(e,"dHK",0.)/(R*t)),x=K*rh;
    if(x>=1.)throw std::domain_error("GAB requires K(T)*relative humidity < 1");
    return q*C*x/((1.-x)*(1.-x+C*x));
  }
  if (m == "gab_water")
    return clip(gab(rh, val(e, "qm_w", val(e, "q_m", q)), val(e, "Cw", val(e, "C_w", 1.)),
                    val(e, "Kw", val(e, "K_w", 0.))),
                0., q);
  if (m == "humidity_toth_climeworks") {
    double T = clip(t, val(cfg, "ISOTHERM_T_FLOOR", 150.), val(cfg, "ISOTHERM_T_CEILING", 2000.)),
           T0 = val(e, "T0", 298.15), ec = val(cfg, "ISOTHERM_EXP_ARG_CLIP", 700.);
    double qm = clip(val(e, "qinf0", val(e, "q_inf_0", q)) *
                         exp_safe(val(e, "chi", 0.) * (1. - T / std::max(T0, tiny)), ec),
                     0., qcap),
           b = clip((val(e, "A", 0.) * rh * exp_safe(val(e, "B", 0.) * rh, ec) + val(e, "b0", 0.)) *
                        exp_safe(val(e, "dH0", val(e, "DeltaH0", 0.)) / (R * T), ec),
                    0., val(cfg, "ISOTHERM_K_CEILING", ceiling));
    double tau = (val(e, "C", 0.) * rh * rh + val(e, "D", 0.) * rh + val(e, "F", 0.)) * rh +
                 val(e, "tau0", 1.) +
                 ((val(e, "G", 0.) * rh * rh + val(e, "H", 0.) * rh + val(e, "J", 0.)) * rh +
                  val(e, "alpha_het", 0.)) *
                     (1. - T0 / T);
    tau = clip(tau, val(cfg, "ISOTHERM_TAU_FLOOR", val(cfg, "ISOTHERM_NU_FLOOR", 1e-12)),
               val(cfg, "ISOTHERM_TAU_CEILING", val(cfg, "ISOTHERM_NU_CEILING", 1e6)));
    double x = clip(b * p, 0., val(cfg, "ISOTHERM_X_CEILING", ceiling)),
           den = std::max(power(1. + power(x, tau, multiply(b, reg)), 1. / tau), 1e-12);
    return clip(qm * x / den, 0., qcap);
  }
  double qm = val(e, "qm_w", val(e, "q_m", 0.)),
         qw = gab(rh, qm, val(e, "Cw", val(e, "C_w", 1.)), val(e, "Kw", val(e, "K_w", 0.)));
  if (m == "gab_modified_toth_liu") {
    double ns = q / std::max(1. - std::max(val(e, "gamma", 0.), 0.) * qw,
                             val(cfg, "LIU_MODIFIED_TOTH_DEN_FLOOR", 1e-8)),
           b = clip(k * (1. + std::max(val(e, "beta", 0.), 0.) * qw), 0.,
                    val(cfg, "ISOTHERM_K_CEILING", ceiling)),
           x = clip(b * p),
           den = std::max(power(1. + power(x, n, multiply(b, reg)), 1. / n), 1e-12);
    return clip(ns * x / den, 0., std::min(std::max(val(e, "qsat_cap", q), q), qcap));
  }
  auto gm = lower(e.value("mode", cfg.value("GAB_COADS_MODE", std::string("both"))));
  bool enhance =
           gm == "enhancement" || gm == "enhance" || gm == "both" || gm == "enhancement_blocking",
       block = gm == "blocking" || gm == "block" || gm == "both" || gm == "enhancement_blocking";
  double keff = clip(enhance ? k * (1. + std::max(val(e, "alpha", 0.), 0.) * qw) : k, 0.,
                     val(cfg, "ISOTHERM_K_CEILING", ceiling)),
         theta = block ? clip(qw / std::max(val(e, "qH2O_block_max",
                                                val(e, "q_h2o_block_max", std::max(qm, tiny))),
                                            tiny),
                              0., 1.)
                       : 0.,
         qeff = q * (1. - theta), x = clip(keff * p);
  if (m == "gab_coadsorption_langmuir") return clip(qeff * x / std::max(1. + x, tiny), 0., q);
  if (m == "gab_coadsorption_sips") {
    x = power(x, n, multiply(keff, reg));
    return clip(qeff * x / std::max(1. + x, tiny), 0., q);
  }
  return clip(qeff * x / std::max(power(1. + power(x, n, multiply(keff, reg)), 1. / n), 1e-12), 0.,
              q);
}
Vec MaterialEquilibrium::Impl::sci(const std::vector<LocalSite>& a, const Vec& p) const {
  Vec q(nc, 0.), kp(nc), terms(nc);
  for (std::size_t s = 0; s < ns; ++s) {
    const auto& m = a[s * nc].ref->model;
    double sum = 0.;
    for (std::size_t j = 0; j < nc; ++j) {
      const auto& x = a[s * nc + j];
      kp[j] = multiply(x.k, p[j]);
      if (m == "langmuir" || m == "anti_langmuir")
        terms[j] = kp[j];
      else if (m == "langmuir_freundlich")
        terms[j] = multiply(x.k, power(p[j], x.nu, reg));
      else
        terms[j] = power(kp[j], x.nu, multiply(x.k, reg));
      sum += terms[j];
    }
    for (std::size_t j = 0; j < nc; ++j) {
      const auto& x = a[s * nc + j];
      double den = m == "anti_langmuir" ? std::max(1. - sum, 1e-12) : 1. + sum;
      if (m == "toth") den = std::max(power(1. + sum, 1. / x.nu), 1e-12);
      q[j] += bound_ratio(x.q, m == "toth" ? kp[j] : terms[j], den);
    }
  }
  return q;
}
Vec MaterialEquilibrium::Impl::sei(const std::vector<LocalSite>& a, const Vec& p) const {
  Vec q(nc, 0.);
  for (std::size_t s = 0; s < ns; ++s) {
    std::vector<std::size_t> active;
    for (std::size_t j = 0; j < nc; ++j) {
      auto& x = a[s * nc + j];
      if (x.ref->model != "zero" && p[j] > active_tol && x.q > 0. && x.k > 0.) active.push_back(j);
    }
    std::stable_sort(active.begin(), active.end(),
                     [&](auto i, auto j) { return a[s * nc + i].q < a[s * nc + j].q; });
    std::size_t n = active.size();
    if (n == 0) continue;
    if (n == 1) {
      auto j = active[0];
      q[j] += basic_loading(a[s * nc + j], p[j], reg);
      continue;
    }
    Vec ratio(n, 1.), a1(n, 0.), a2(n, 0.), ap(n, 1.);
    auto qs = [&](std::size_t j) { return a[s * nc + active[j]].q; };
    auto kd = [&](std::size_t i) {
      auto j = active[i];
      return clip(a[s * nc + j].k * p[j]);
    };
    for (std::size_t i = 1; i < n; ++i) ratio[i] = qs(i) / std::max(qs(i - 1), tiny);
    a2[n - 1] = 1. + kd(n - 1);
    a1[n - 1] = power(a2[n - 1], ratio[n - 1]);
    for (int i = int(n) - 2; i > 0; --i) {
      a2[i] = a1[i + 1] + kd(i);
      a1[i] = power(a2[i], ratio[i]);
    }
    a2[0] = a1[1] + kd(0);
    a1[0] = a2[0];
    double beta = std::max(a2[0], tiny);
    for (std::size_t i = 1; i < n; ++i) ap[i] = (a1[i] / std::max(a2[i], tiny)) * ap[i - 1];
    for (std::size_t i = 0; i < n; ++i)
      q[active[i]] += clip(qs(i) * kd(i) * ap[i] / beta, 0., qs(i));
  }
  return q;
}
Vec MaterialEquilibrium::Impl::fast_iast(const std::vector<LocalSite>& a, const Vec& p, int site,
                                         const Vec& warm, Vec& pstar) const {
  Vec result(nc, 0.);
  pstar.assign(nc, 0.);
  double total = std::accumulate(p.begin(), p.end(), 0.);
  if (total <= active_tol) return result;
  std::size_t begin = site < 0 ? 0 : std::size_t(site), end = site < 0 ? ns : begin + 1;
  auto psi = [&](std::size_t j, double d) {
    double x = 0.;
    for (std::size_t s = begin; s < end; ++s) x += spreading(a[s * nc + j], d, reg);
    return clip(x);
  };
  auto load = [&](std::size_t j, double d) {
    double x = 0.;
    for (std::size_t s = begin; s < end; ++s) x += basic_loading(a[s * nc + j], d, reg);
    return clip(x);
  };
  double ftol = std::max(val(cfg, "IAST_FAST_REL_TOL", tol), 1e-15),
         xtol = std::max(val(cfg, "IAST_FAST_X_TOL", 1e-10), 1e-15),
         maxstep = std::max(val(cfg, "IAST_FAST_MAX_LOG_STEP", 8.), 1e-6);
  int maxit = std::max(int(val(cfg, "IAST_FAST_MAX_ITER", 50)), 1),
      lines = std::max(int(val(cfg, "IAST_FAST_LINESEARCH_STEPS", 16)), 1);
  auto inverse = [&](std::size_t j, double target, double cached) {
    if (end - begin == 1)
      return inverse_site(a[begin * nc + j], target, cached, driving_cap, ftol, iters, reg);
    return inverse_bisection([&](double d) { return psi(j, d); }, target, cached, driving_cap, ftol,
                             iters);
  };
  std::vector<std::size_t> active;
  double initial = 0.;
  for (std::size_t j = 0; j < nc; ++j) {
    if (p[j] <= active_tol) continue;
    double x = psi(j, total);
    if (x <= active_tol) continue;
    active.push_back(j);
    initial += p[j] / total * x;
  }
  if (active.empty() || initial <= active_tol) return result;
  if (active.size() == 1) {
    auto j = active[0];
    result[j] = load(j, p[j]);
    pstar[j] = std::max(p[j], tiny);
    return result;
  }
  std::size_t n = active.size();
  Vec z(n);
  bool goodwarm = warm.size() == nc;
  for (auto j : active)
    if (!goodwarm || !std::isfinite(warm[j]) || warm[j] <= 0.) {
      goodwarm = false;
      break;
    }
  double cached = 0.;
  for (std::size_t i = 0; i < n; ++i) {
    double x =
        goodwarm ? clip(warm[active[i]], tiny, driving_cap) : inverse(active[i], initial, cached);
    if (!std::isfinite(x) || x <= 0.)
      throw std::runtime_error("Invalid FastIAST starting pressure");
    cached = x;
    z[i] = std::log(clip(x, tiny, driving_cap));
  }
  struct Eval {
    Vec residual, jacobian, pressure, qpure, x;
    double pe = 0., xe = 0., merit = 0.;
  };
  auto evaluate = [&](const Vec& logp) {
    Eval e;
    e.residual.resize(n);
    e.jacobian.assign(n * n, 0.);
    e.pressure.resize(n);
    e.qpure.resize(n);
    e.x.resize(n);
    Vec spread(n);
    double mean = 0., absmean = 0., sx = 0.;
    for (std::size_t i = 0; i < n; ++i) {
      e.pressure[i] =
          std::exp(clip(logp[i], std::log(tiny), std::log(std::max(driving_cap, tiny))));
      spread[i] = psi(active[i], e.pressure[i]);
      e.qpure[i] = load(active[i], e.pressure[i]);
      e.x[i] = p[active[i]] / std::max(e.pressure[i], tiny);
      mean += spread[i];
      absmean += std::abs(spread[i]);
      sx += e.x[i];
      if (!std::isfinite(spread[i]) || !std::isfinite(e.qpure[i]) || !std::isfinite(e.x[i]))
        throw std::runtime_error("Nonfinite FastIAST residual");
    }
    std::size_t ref = std::max_element(e.qpure.begin(), e.qpure.end()) - e.qpure.begin(), row = 0;
    for (std::size_t i = 0; i < n; ++i)
      if (i != ref) {
        double scale = std::max({std::abs(spread[i]), std::abs(spread[ref]), active_tol});
        e.residual[row] = (spread[i] - spread[ref]) / scale;
        e.jacobian[row * n + i] = e.qpure[i] / scale;
        e.jacobian[row * n + ref] = -e.qpure[ref] / scale;
        ++row;
      }
    e.residual[n - 1] = sx - 1.;
    for (std::size_t i = 0; i < n; ++i) e.jacobian[(n - 1) * n + i] = -e.x[i];
    mean /= n;
    absmean = std::max(absmean / n, active_tol);
    for (std::size_t i = 0; i < n; ++i) {
      e.pe = std::max(e.pe, std::abs(spread[i] - mean) / absmean);
      e.merit = std::max(e.merit, std::abs(e.residual[i]));
    }
    e.xe = std::abs(sx - 1.);
    e.merit = std::max({e.merit, e.pe, e.xe});
    return e;
  };
  for (int iteration = 0; iteration <= maxit; ++iteration) {
    auto e = evaluate(z);
    if (e.pe <= ftol && e.xe <= xtol) {
      double sx = std::accumulate(e.x.begin(), e.x.end(), 0.), iq = 0.;
      if (!std::isfinite(sx) || sx <= active_tol)
        throw std::runtime_error("Invalid FastIAST composition");
      for (std::size_t i = 0; i < n; ++i)
        if (e.x[i] / sx > active_tol && e.qpure[i] > active_tol) iq += e.x[i] / sx / e.qpure[i];
      if (!std::isfinite(iq) || iq <= active_tol)
        throw std::runtime_error("Invalid FastIAST loading");
      for (std::size_t i = 0; i < n; ++i) {
        result[active[i]] = clip(e.x[i] / sx / iq);
        pstar[active[i]] = e.pressure[i];
      }
      return result;
    }
    if (iteration >= maxit) break;
    Vec rhs = e.residual;
    for (auto& v : rhs) v = -v;
    auto delta = dense_solve(e.jacobian, rhs);
    double largest = 0.;
    for (auto v : delta) largest = std::max(largest, std::abs(v));
    if (largest > maxstep)
      for (auto& v : delta) v *= maxstep / largest;
    bool accepted = false;
    double alpha = 1.;
    for (int l = 0; l < lines; ++l) {
      auto trial = z;
      for (std::size_t i = 0; i < n; ++i) trial[i] += alpha * delta[i];
      double merit = std::numeric_limits<double>::infinity();
      try {
        merit = evaluate(trial).merit;
      } catch (const std::runtime_error&) {
      }
      if (std::isfinite(merit) && merit < e.merit * (1. - 1e-4 * alpha)) {
        z = std::move(trial);
        accepted = true;
        break;
      }
      alpha *= .5;
    }
    if (!accepted) throw std::runtime_error("FastIAST line search failed");
  }
  throw std::runtime_error("FastIAST iteration limit reached");
}
Vec MaterialEquilibrium::Impl::iast(const std::vector<LocalSite>& a, const Vec& p, int site,
                                    const Vec& warm, Vec& pstar) const {
  pstar.clear();
  bool fast = site < 0 && ns > 1;
  for (std::size_t j = 0; j < nc; ++j) {
    auto& m = a[(site < 0 ? 0 : std::size_t(site)) * nc + j].ref->model;
    if (m == "zero") continue;
    if ((reg > 0. && powered(m)) || (m != "langmuir" && m != "henry" && m != "freundlich" &&
                                     m != "langmuir_freundlich" && m != "sips"))
      fast = true;
  }
  if (method == "FASTIAST" || (method == "AUTO" && fast)) {
    try {
      return fast_iast(a, p, site, warm, pstar);
    } catch (const std::runtime_error&) {
      pstar.clear();
      if (method == "FASTIAST" && !flag(cfg, "IAST_FAST_FALLBACK", true)) throw;
    }
  }
  Vec result(nc, 0.);
  double total = std::accumulate(p.begin(), p.end(), 0.);
  if (total <= active_tol) return result;
  std::size_t begin = site < 0 ? 0 : std::size_t(site), end = site < 0 ? ns : begin + 1;
  auto psi = [&](std::size_t j, double d) {
    double x = 0.;
    for (std::size_t s = begin; s < end; ++s) x += spreading(a[s * nc + j], d, reg);
    return clip(x);
  };
  auto load = [&](std::size_t j, double d) {
    double x = 0.;
    for (std::size_t s = begin; s < end; ++s) x += basic_loading(a[s * nc + j], d, reg);
    return clip(x);
  };
  auto inverse = [&](std::size_t j, double target, double cached) {
    if (end - begin == 1)
      return inverse_site(a[begin * nc + j], target, cached, driving_cap, tol, iters, reg);
    return inverse_bisection([&](double d) { return psi(j, d); }, target, cached, driving_cap, tol,
                             iters);
  };
  double initial = 0.;
  std::vector<std::size_t> active;
  for (std::size_t j = 0; j < nc; ++j) {
    if (p[j] <= active_tol) continue;
    double x = psi(j, total);
    if (x <= active_tol) continue;
    active.push_back(j);
    initial += p[j] / total * x;
  }
  if (active.empty() || initial <= active_tol) return result;
  if (active.size() == 1) {
    auto j = active[0];
    result[j] = load(j, p[j]);
    return result;
  }
  Vec cache(nc, 0.);
  auto sumx = [&](double target) {
    double x = 0.;
    for (auto j : active) {
      cache[j] = inverse(j, target, cache[j]);
      x += p[j] / std::max(cache[j], tiny);
    }
    return x;
  };
  double l = initial, r = initial, sx = sumx(initial);
  int steps = 0;
  if (sx > 1.)
    while (sx > 1.) {
      r *= 2.;
      if (!std::isfinite(r) || r > ceiling || ++steps > iters)
        throw std::runtime_error("IAST right bracket failed");
      sx = sumx(r);
    }
  else if (sx < 1.)
    while (sx < 1.) {
      l *= .5;
      if (l <= tiny) {
        l = 0.;
        break;
      }
      if (++steps > iters) throw std::runtime_error("IAST left bracket failed");
      sx = sumx(l);
    }
  double target = 0.;
  if (l > 0.) {
    bool converged = false;
    for (int i = 0; i < iters; ++i) {
      double mid = .5 * (l + r);
      sx = sumx(mid);
      if (sx > 1.)
        l = mid;
      else
        r = mid;
      if (std::abs(l - r) / std::max(std::abs(l + r), tiny) <= tol) {
        converged = true;
        break;
      }
    }
    if (!converged) throw std::runtime_error("IAST did not converge");
    target = .5 * (l + r);
  }
  sx = sumx(target);
  double iq = 0.;
  for (auto j : active) {
    cache[j] = inverse(j, target, cache[j]);
    result[j] = p[j] / std::max(cache[j], tiny) / std::max(sx, tiny);
    double pureq = load(j, cache[j]);
    if (result[j] > active_tol && pureq > active_tol) iq += result[j] / pureq;
  }
  if (iq <= active_tol) return Vec(nc, 0.);
  for (auto j : active) result[j] = clip(result[j] / iq);
  return result;
}
Vec MaterialEquilibrium::Impl::loading(const Vec& pp, const Vec& t) const {
  if(main_model)return main_model->equilibrium(pp,t,nz);
  if (pp.size() != nc * nz || (t.size() != nz && t.size() != 1))
    throw std::invalid_argument("Equilibrium array dimensions");
  for (double x : pp)
    if (!std::isfinite(x)) throw std::invalid_argument("Equilibrium pressure must be finite");
  for (double x : t)
    if (!std::isfinite(x) || x <= 0.)
      throw std::invalid_argument("Equilibrium temperature must be positive");
  Vec pressure(pp.size());
  std::transform(pp.begin(), pp.end(), pressure.begin(),
                 [&](double x) { return clip(x, 0., driving_cap); });
  if (cache_enabled && pressure == cache_p && t == cache_t) {
    ++hits;
    return cache_q;
  }
  if (cache_enabled) ++misses;
  Vec out(nc * nz, 0.), p(nc);
  std::map<int, std::pair<Vec, Vec>> warmcache;
  for (std::size_t z = 0; z < nz; ++z) {
    double T = t[t.size() == 1 ? 0 : z];
    auto a = local(T);
    for (std::size_t j = 0; j < nc; ++j) p[j] = pressure[j * nz + z];
    Vec q(nc, 0.);
    if (mode == "SCI")
      q = sci(a, p);
    else if (mode == "SEI")
      q = sei(a, p);
    else if (mode == "SPI") {
      double rh = needs_rh ? relative_humidity(water < 0 ? 0. : p[water], T, cfg) : 0.;
      for (std::size_t s = 0; s < ns; ++s)
        for (std::size_t j = 0; j < nc; ++j) q[j] += pure(a[s * nc + j], p[j], T, rh);
    } else {
      int first = mode == "IAST" ? -1 : 0, last = mode == "IAST" ? 0 : int(ns);
      for (int s = first; s < last; ++s) {
        Vec warm, pstar;
        auto cached = warmcache.find(s);
        if (flag(cfg, "IAST_WARM_START_GRID", true) && cached != warmcache.end()) {
          bool valid = false;
          const auto& old = cached->second;
          double ratio = std::max(val(cfg, "IAST_WARM_START_MAX_RATIO", 100.), 1.);
          for (std::size_t j = 0; j < nc; ++j)
            if (old.first[j] > 0.) {
              valid = true;
              if (p[j] <= active_tol || old.second[j] <= active_tol ||
                  std::max(p[j] / old.second[j], old.second[j] / p[j]) > ratio) {
                valid = false;
                break;
              }
            }
          if (valid) warm = old.first;
        }
        auto qs = iast(a, p, s, warm, pstar);
        for (std::size_t j = 0; j < nc; ++j) q[j] += qs[j];
        if (flag(cfg, "IAST_WARM_START_GRID", true) && !pstar.empty() &&
            std::any_of(pstar.begin(), pstar.end(), [](double v) { return v > 0.; }))
          warmcache[s] = {pstar, p};
        else
          warmcache.erase(s);
      }
    }
    for (std::size_t j = 0; j < nc; ++j) out[j * nz + z] = clip(q[j], 0., qcap);
  }
  if (cache_enabled) {
    cache_p = std::move(pressure);
    cache_t = t;
    cache_q = out;
  }
  return out;
}
MaterialEquilibrium::MaterialEquilibrium(const Json& c, double r, std::size_t n, std::size_t z,
                                         const std::string& phase)
    : impl_(std::make_unique<Impl>(c, r, n, z, phase)) {}
MaterialEquilibrium::~MaterialEquilibrium() = default;
MaterialEquilibrium::MaterialEquilibrium(MaterialEquilibrium&&) noexcept = default;
MaterialEquilibrium& MaterialEquilibrium::operator=(MaterialEquilibrium&&) noexcept = default;
MaterialEquilibrium::MaterialEquilibrium(const MaterialEquilibrium& x)
    : impl_(std::make_unique<Impl>(*x.impl_)) {}
MaterialEquilibrium& MaterialEquilibrium::operator=(const MaterialEquilibrium& x) {
  if (this != &x) impl_ = std::make_unique<Impl>(*x.impl_);
  return *this;
}
Vec MaterialEquilibrium::loading(const Vec& p, const Vec& t) const { return impl_->loading(p, t); }
const Vec& MaterialEquilibrium::qmax_effective() const { return impl_->capacities; }
const std::vector<bool>& MaterialEquilibrium::bounded_mask() const { return impl_->bounded; }
Json MaterialEquilibrium::cache_diagnostics() const {
  return {{"hits", impl_->hits}, {"misses", impl_->misses}};
}
}  // namespace cycle

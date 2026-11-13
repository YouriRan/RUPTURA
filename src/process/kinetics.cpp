#include "cycle/kinetics.hpp"
#include "chemisorption.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>
namespace cycle {
namespace {
constexpr double R = 8.31446261815324;
[[noreturn]] void fail(const std::string& key, const std::string& reason) {
  throw std::invalid_argument("CHEMISORPTION." + key + " " + reason);
}
std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
double num(const Json& v, const std::string& name) {
  if (!v.is_number()) fail(name, "must contain numbers");
  double x = v.get<double>();
  if (!std::isfinite(x)) fail(name, "must be finite");
  return x;
}
double scalar(const Json& cfg, const char* key, double d, bool positive = false) {
  double x = num(cfg.value(key, Json(d)), key);
  if (positive && x <= 0) fail(key, "must be positive");
  return x;
}
Vec broadcast(const Json& v, std::size_t nb, std::size_t nc, const std::string& name,
              bool nonnegative = false, bool positive = false) {
  Vec out(nb * nc);
  auto put = [&](std::size_t i, const Json& item) {
    double x = num(item, name);
    if (positive && x <= 0) fail(name, "must be positive");
    if (nonnegative && x < 0) fail(name, "must be nonnegative");
    out[i] = x;
  };
  if (v.is_number()) {
    for (std::size_t i = 0; i < out.size(); ++i) put(i, v);
    return out;
  }
  if (!v.is_array() || v.empty()) fail(name, "has incompatible broadcast shape");
  if (!v.front().is_array()) {
    const bool branch = v.size() == nb && v.size() != nc;
    if (!branch && v.size() != 1 && v.size() != nc) fail(name, "has incompatible broadcast shape");
    for (std::size_t b = 0; b < nb; ++b)
      for (std::size_t c = 0; c < nc; ++c) put(b * nc + c, v[branch ? b : (v.size() == 1 ? 0 : c)]);
  } else {
    if (v.size() != 1 && v.size() != nb) fail(name, "has incompatible broadcast shape");
    const auto width = v.front().size();
    if (width != 1 && width != nc) fail(name, "has incompatible broadcast shape");
    for (const auto& row : v)
      if (!row.is_array() || row.size() != width) fail(name, "must be a rectangular matrix");
    for (std::size_t b = 0; b < nb; ++b)
      for (std::size_t c = 0; c < nc; ++c)
        put(b * nc + c, v[v.size() == 1 ? 0 : b][width == 1 ? 0 : c]);
  }
  return out;
}
Vec component(const Json& v, std::size_t nc, const std::string& name, bool nonnegative = false,
              bool positive = false) {
  if (v.is_array() && !v.empty() && v.front().is_array())
    fail(name, "has incompatible broadcast shape");
  return broadcast(v, 1, nc, name, nonnegative, positive);
}
Vec nodes(const Vec& v, std::size_t nz, const std::string& name, bool positive) {
  if (v.size() != 1 && v.size() != nz) fail(name, "has incompatible broadcast shape");
  Vec out(nz);
  for (std::size_t z = 0; z < nz; ++z) {
    double x = v[v.size() == 1 ? 0 : z];
    if (!std::isfinite(x)) fail(name, "must be finite");
    if (positive ? x <= 0 : x < 0)
      fail(name, positive ? "must be positive" : "must be nonnegative");
    out[z] = x;
  }
  return out;
}
double checked(double x, const char* kind) {
  if (!std::isfinite(x))
    throw std::invalid_argument(std::string("Chemical ") + kind +
                                " overflow; check parameters and temperature units");
  return x;
}
}  // namespace
ChemistryKinetics::ChemistryKinetics(const Json& cfg, std::size_t nc) : nc_(nc) {
  if (!cfg.is_object()) throw std::invalid_argument("CHEMISORPTION must be a mapping");
  if (!nc_) throw std::invalid_argument("Chemistry requires at least one component");
  model_ = lower(cfg.value("model", std::string("pfo")));
  std::replace(model_.begin(), model_.end(), '-', '_');
  if (model_ == "pseudo_first_order") model_ = "pfo";
  if (model_ == "pseudo_second_order") model_ = "pso";
  if (model_ != "pfo" && model_ != "pso" && model_ != "avrami" && model_ != "dual_nth")
    fail("model", "must be pfo, pso, avrami, or dual_nth");
  nb_ = model_ == "dual_nth" ? 2 : 1;
  nh_ = model_ == "avrami" ? 2 : 0;
  auto param = [&](const Json& v, const std::string& name, bool nonnegative = false,
                   bool positive = false) {
    return broadcast(v, nb_, nc_, name, nonnegative, positive);
  };
  if (!cfg.contains("k_ads")) fail("k_ads", "is required; see model-specific rate units");
  ka_ = param(cfg["k_ads"], "k_ads", true);
  kd_ = param(cfg.value("k_des", cfg["k_ads"]), "k_des", true);
  eaa_ = param(cfg.value("Ea_ads", Json(0.)), "Ea_ads");
  ead_ = param(cfg.value("Ea_des", Json(0.)), "Ea_des");
  tref_ = scalar(cfg, "T_ref", 298.15, true);
  if (model_ == "dual_nth") {
    if (!cfg.contains("fractions")) fail("fractions", "is required for dual_nth");
    Json f = cfg["fractions"];
    if (f.is_array() && f.size() == 2 && !f.front().is_array())
      f = Json::array({Json::array({f[0]}), Json::array({f[1]})});
    fractions_ = param(f, "fractions", true);
    for (std::size_t c = 0; c < nc_; ++c)
      if (std::abs(fractions_[c] + fractions_[nc_ + c] - 1.) > 1e-12)
        fail("fractions", "must sum to one for every component");
    oa_ = param(cfg.value("order_ads", cfg.value("order", Json(1.))), "order_ads", false, true);
    od_ = param(cfg.value("order_des", cfg.value("order", Json(1.))), "order_des", false, true);
  } else {
    fractions_.assign(nb_ * nc_, 1.);
    oa_.assign(nb_ * nc_, model_ == "pso" ? 2. : 1.);
    od_ = oa_;
  }
  basis_ = lower(cfg.value("rate_basis", std::string("loading")));
  if (basis_ != "loading" && basis_ != "fractional")
    fail("rate_basis", "must be loading or fractional");
  if (model_ == "avrami" && basis_ != "loading") fail("rate_basis", "must be loading for Avrami");
  const Json eq = cfg.value("equilibrium", Json());
  if (!eq.is_null() && !eq.is_object()) fail("equilibrium", "must be a mapping");
  if (basis_ == "fractional") {
    auto ref = cfg.value("q_reference", eq.is_object() ? eq.value("qsat", Json()) : Json());
    if (ref.is_null()) fail("q_reference", "or equilibrium.qsat is required for fractional rates");
    auto r = component(ref, nc_, "q_reference", false, true);
    qref_ = fractions_;
    for (std::size_t b = 0; b < nb_; ++b)
      for (std::size_t c = 0; c < nc_; ++c) (*qref_)[b * nc_ + c] *= r[c];
  }
  drive_eps_ = scalar(cfg, "driving_regularizer", 1e-12, true);
  auto heat = cfg.value("heat_release_J_mol",
                        eq.is_object() ? eq.value("heat_release_J_mol", Json()) : Json());
  if (!heat.is_null()) heat_ = param(heat, "heat_release_J_mol");
  if (model_ == "avrami") {
    na_ = param(cfg.value("n_ads", cfg.value("n", Json(1.))), "n_ads", false, true);
    nd_ = param(cfg.value("n_des", cfg.value("n", Json(1.))), "n_des", false, true);
    age_eps_ = scalar(cfg, "age_regularizer", 1e-6, true);
    if (cfg.value("history_policy", std::string("continuous_erasure")) != "continuous_erasure")
      fail("history_policy", "must be continuous_erasure; stages do not reset history");
    if (cfg.contains("memory_erasure_rate_ads"))
      erase_a_ = param(cfg["memory_erasure_rate_ads"], "memory_erasure_rate_ads", false, true);
    if (cfg.contains("memory_erasure_rate_des"))
      erase_d_ = param(cfg["memory_erasure_rate_des"], "memory_erasure_rate_des", false, true);
  }
  has_eq_ = !eq.is_null();
  if (!has_eq_) return;
  eq_model_ = lower(eq.value("model", std::string("competitive_langmuir")));
  if (eq_model_ != "langmuir" && eq_model_ != "competitive_langmuir" &&
      eq_model_ != "extended_langmuir" && eq_model_ != "henry")
    fail("equilibrium.model",
         "must be langmuir, competitive_langmuir, extended_langmuir, or henry");
  eq_tref_ = scalar(eq, "T_ref", tref_, true);
  eq_heat_ = component(eq.value("heat_release_J_mol", eq.value("heat_affinity_J_mol", Json(0.))),
                       nc_, "equilibrium.heat_release_J_mol");
  if (eq_model_ == "henry") {
    if (!eq.contains("K")) fail("equilibrium.K", "is required for Henry equilibrium");
    eq_k_ = component(eq["K"], nc_, "equilibrium.K", true);
  } else {
    if (!eq.contains("qsat") || !eq.contains("b"))
      fail("equilibrium", "requires qsat and b for Langmuir equilibrium");
    eq_qsat_ = component(eq["qsat"], nc_, "equilibrium.qsat", true);
    eq_b_ = component(eq["b"], nc_, "equilibrium.b", true);
  }
}
std::optional<Vec> ChemistryKinetics::equilibrium(const Vec& y, const Vec& pressure,
                                                  const Vec& temperature) const {
  if (!has_eq_) return std::nullopt;
  if (y.empty() || y.size() % nc_) fail("equilibrium y", "must have shape (components,nodes)");
  const auto nz = y.size() / nc_;
  auto p = nodes(pressure, nz, "equilibrium pressure", false),
       ts = nodes(temperature, nz, "equilibrium temperature", true);
  for (double yi : y)
    if (!std::isfinite(yi) || yi < 0) fail("equilibrium y", "must be finite and nonnegative");
  Vec out(y.size()), bp(nc_);
  for (std::size_t z = 0; z < nz; ++z) {
    double sum = 0.;
    for (std::size_t c = 0; c < nc_; ++c) sum += y[c * nz + z];
    if (std::abs(sum - 1.) > 2e-8) fail("equilibrium y", "must sum to one");
    double denominator = 0.;
    for (std::size_t c = 0; c < nc_; ++c) {
      double f = checked(std::exp(eq_heat_[c] / R * (1. / ts[z] - 1. / eq_tref_)), "equilibrium"),
             partial = checked(y[c * nz + z] * p[z], "equilibrium");
      if (eq_model_ == "henry")
        out[c * nz + z] = checked(checked(eq_k_[c] * f, "equilibrium") * partial, "equilibrium");
      else {
        bp[c] = checked(checked(eq_b_[c] * f, "equilibrium") * partial, "equilibrium");
        if (eq_model_ != "langmuir") denominator = checked(denominator + bp[c], "equilibrium");
      }
    }
    // Match NumPy: first sum bp, then add unity (critical at zero driving force).
    if (eq_model_ != "henry")
      for (std::size_t c = 0; c < nc_; ++c)
        out[c * nz + z] = checked(eq_qsat_[c] * bp[c], "equilibrium") /
                          (1. + (eq_model_ == "langmuir" ? bp[c] : denominator));
  }
  return out;
}
Vec ChemistryKinetics::split_loading(const Vec& q, std::size_t nz) const {
  if (!nz || q.size() != nc_ * nz) fail("total loading", "must have shape (components,nodes)");
  Vec out(nb_ * q.size());
  for (std::size_t b = 0; b < nb_; ++b)
    for (std::size_t c = 0; c < nc_; ++c)
      for (std::size_t z = 0; z < nz; ++z)
        out[(b * nc_ + c) * nz + z] = fractions_[b * nc_ + c] * q[c * nz + z];
  return out;
}
Vec ChemistryKinetics::initial_history(std::size_t nz) const {
  if (!nz) throw std::invalid_argument("n_nodes must be positive");
  return Vec(nh_ * nc_ * nz, 0.);
}
std::pair<Vec, Vec> ChemistryKinetics::evaluate(const Vec& qeq, const Vec& loads,
                                                const Vec& history, const Vec& temperature) const {
  if (qeq.empty() || qeq.size() % nc_) fail("qeq_total", "must have shape (components,nodes)");
  const auto nz = qeq.size() / nc_;
  if (loads.size() != nb_ * qeq.size())
    fail("branch loading", "shape is incompatible with kinetics");
  if (history.size() != nh_ * qeq.size()) fail("history", "shape is incompatible with kinetics");
  for (double q : qeq)
    if (!std::isfinite(q) || q < 0) fail("targets", "must be finite and nonnegative");
  for (const auto* v : {&loads, &history})
    for (double x : *v)
      if (!std::isfinite(x)) fail("states", "must be finite");
  auto ts = nodes(temperature, nz, "solid temperature", true);
  Vec rate(loads.size()), dh(history.size(), 0.);
  for (std::size_t c = 0; c < nc_; ++c) {
    Chemisorption law;
    law.type = model_ == "dual_nth" ? Chemisorption::Type::DualNth : model_ == "avrami" ? Chemisorption::Type::Avrami
        : model_ == "pso" ? Chemisorption::Type::PseudoNth : Chemisorption::Type::FirstOrder;
    law.directionalRates = true; law.referenceTemperature = tref_; law.drivingRegularizer = drive_eps_;
    law.fractionalRates = basis_ == "fractional";
    for (std::size_t b = 0; b < nb_; ++b) {
      const auto bc = b * nc_ + c;
      auto& branch = law.branches[b];
      branch.fraction = fractions_[bc]; branch.adsorptionRate = ka_[bc]; branch.desorptionRate = kd_[bc];
      branch.adsorptionActivationEnergy = eaa_[bc]; branch.desorptionActivationEnergy = ead_[bc];
      branch.adsorptionOrder = oa_[bc]; branch.desorptionOrder = od_[bc];
      if (qref_ && fractions_[bc] > 0.0) law.referenceLoading = (*qref_)[bc] / fractions_[bc];
    }
    law.adsorptionRateCoefficient = ka_[c]; law.desorptionRateCoefficient = kd_[c];
    law.adsorptionActivationEnergy = eaa_[c]; law.desorptionActivationEnergy = ead_[c];
    law.order = model_ == "avrami" ? na_[c] : oa_[c];
    law.desorptionOrderValue = model_ == "avrami" ? nd_[c] : od_[c];
    law.continuousHistory = model_ == "avrami"; law.ageRegularizer = age_eps_;
    if (erase_a_) law.adsorptionMemoryErasure = (*erase_a_)[c];
    if (erase_d_) law.desorptionMemoryErasure = (*erase_d_)[c];
    for (std::size_t z = 0; z < nz; ++z) {
      std::array<double,2> q{};
      for (std::size_t b = 0; b < nb_; ++b) q[b] = loads[(b * nc_ + c) * nz + z];
      const auto rates = law.evaluate(qeq[c*nz+z], std::span<const double>(q.data(),nb_), 0.0, ts[z], 0.0,
          nh_ ? history[c*nz+z] : 0.0, nh_ ? history[qeq.size()+c*nz+z] : 0.0);
      for (std::size_t b = 0; b < nb_; ++b) rate[(b * nc_ + c)*nz+z] = rates.loading[b];
      if (nh_) { dh[c*nz+z] = rates.adsorptionAge; dh[qeq.size()+c*nz+z] = rates.desorptionAge; }
    }
  }
  return {std::move(rate), std::move(dh)};
}
std::optional<Vec> ChemistryKinetics::heat_rate(const Vec& rates, std::size_t nz,
                                                const std::optional<Vec>& fallback) const {
  auto heat = heat_;
  if (!heat && fallback) {
    if (nb_ > 1 && fallback->size() == nb_ * nc_) {
      heat = *fallback;
      for (double h : *heat)
        if (!std::isfinite(h)) fail("heat fallback", "must be finite");
    } else
      heat = broadcast(Json(*fallback), nb_, nc_, "heat fallback");
  }
  if (!heat) return std::nullopt;
  if (rates.size() != nb_ * nc_ * nz) fail("heat rates", "shape is incompatible with kinetics");
  Vec out(nz, 0.);
  for (std::size_t bc = 0; bc < nb_ * nc_; ++bc)
    for (std::size_t z = 0; z < nz; ++z) out[z] += rates[bc * nz + z] * (*heat)[bc];
  return out;
}
}  // namespace cycle

namespace cycle {
Json ChemistryKinetics::canonical_site(std::size_t c) const {
  if (c >= nc_) throw std::invalid_argument("Chemisorption component index");
  Json params = {{"ReferenceTemperature",tref_},{"RateBasis",basis_},{"DrivingRegularizer",drive_eps_}};
  Json site = {{"RateEquation",model_=="dual_nth"?"DualNth":model_=="avrami"?"Avrami":model_=="pso"?"PseudoNth":"FirstOrder"},
    {"MaximumLoading",1e100}};
  if (qref_) for (std::size_t b=0;b<nb_;++b) if (fractions_[b*nc_+c]>0.0)
    params["ReferenceLoading"] = (*qref_)[b*nc_+c]/fractions_[b*nc_+c];
  if (model_=="dual_nth") {
    params["Branches"]=Json::array();
    for(std::size_t b=0;b<nb_;++b) {const auto bc=b*nc_+c;
      Json branch={{"Fraction",fractions_[bc]},{"AdsorptionRateCoefficient",ka_[bc]},{"DesorptionRateCoefficient",kd_[bc]},
        {"AdsorptionActivationEnergy",eaa_[bc]},{"DesorptionActivationEnergy",ead_[bc]},
        {"AdsorptionOrder",oa_[bc]},{"DesorptionOrder",od_[bc]}};
      if(heat_)branch["HeatOfAdsorption"]=(*heat_)[bc];params["Branches"].push_back(branch);
    }
  } else {
    params.update({{"AdsorptionRateCoefficient",ka_[c]},{"DesorptionRateCoefficient",kd_[c]},
      {"AdsorptionActivationEnergy",eaa_[c]},{"DesorptionActivationEnergy",ead_[c]},
      {"Order",model_=="avrami"?na_[c]:oa_[c]},{"DesorptionOrder",model_=="avrami"?nd_[c]:od_[c]}});
    if(model_=="avrami") {
      params.update({{"HistoryPolicy","ContinuousErasure"},{"AgeRegularizer",age_eps_}});
      if(erase_a_)params["MemoryErasureRateAds"]=(*erase_a_)[c];
      if(erase_d_)params["MemoryErasureRateDes"]=(*erase_d_)[c];
    }
    if(heat_)site["HeatOfAdsorption"]=(*heat_)[c];
  }
  if(has_eq_) {
    site["ReferenceTemperature"]=eq_tref_;site["EquilibriumHeat"]=eq_heat_[c];
    site["MixturePredictionMethod"]=eq_model_=="competitive_langmuir"||eq_model_=="extended_langmuir"?"SCI":"SPI";
    if(eq_model_=="henry") {site["Type"]="Henry";site["Parameters"]=Json::array({eq_k_[c]});}
    else {site["Type"]="Langmuir";site["Parameters"]=Json::array({eq_qsat_[c],eq_b_[c]});site["MaximumLoading"]=std::max(eq_qsat_[c],1e-300);}
  } else {site["Type"]="Henry";site["Parameters"]=Json::array({0.0});site["EquilibriumSource"]="Physisorption";}
  site["AdsorptionReferenceTemperature"]=has_eq_?eq_tref_:tref_;
  site["RateEquationParameters"]=params;
  return site;
}
}

#include <algorithm>
#include <cmath>
#include <numeric>
#include <optional>
#include <stdexcept>

#include "cycle/network.hpp"
namespace cycle
{
namespace
{
double root(double dp, double a, double b)
{
  if (dp == 0.) return 0.;
  if (a <= 0. && b <= 0.) throw std::invalid_argument("Liquid hydraulic resistance must be positive");
  return std::copysign(2. * std::abs(dp) / (a + std::hypot(a, 2. * std::sqrt(b * std::abs(dp)))), dp);
}
struct Edge
{
  int a = -1, b = -1;
  double pa = 0., pb = 0., linear = 0., quadratic = 0., opening = 1.;
  std::optional<double> imposed;
  bool check = false;
  std::string name;
};
}  // namespace
Json ColumnNetwork::solve_liquid_hydraulics(double t, const Vec& state, const Json& sources) const
{
  const std::size_t n = 2 * models.size();
  Vec P(n), external(n, 0.);
  std::vector<bool> fixed(n, false), referenced(n, false);
  std::vector<Edge> edges;
  std::vector<std::vector<std::size_t>> adjacent(n);
  auto index = [&](const Json& ep) -> int
  {
    if (ep[0] != "column") return -1;
    return static_cast<int>(2 * column_indices.at(ep[1].get<std::string>()) + (ep[2] == "right" ? 1 : 0));
  };
  auto reservoir_pressure = [&](const Json& ep)
  {
    if (ep[0] == "source") return sources.at(ep[1].get<std::string>()).at("pressure").get<double>();
    return buffers
        ->state_record(ep[1].get<std::string>(),
                       Vec(state.begin() + buffer_start, state.begin() + buffer_start + buffers->size))
        .at("pressure")
        .get<double>();
  };
  for (std::size_t ci = 0; ci < models.size(); ++ci)
  {
    auto [a, b] = models[ci].liquid_resistance();
    edges.push_back(
        {int(2 * ci), int(2 * ci + 1), 0., 0., a, b, 1., std::nullopt, false, "column_" + column_names[ci]});
    for (std::size_t k = 0; k < 2; ++k)
    {
      std::size_t i = 2 * ci + k;
      P[i] = models[ci].cfg.value("P_init", 1e5);
      auto port = models[ci].port_state(k == 0 ? "left" : "right", t);
      if (port["mode"] == "pressure")
      {
        fixed[i] = referenced[i] = true;
        P[i] = port["pressure"];
      }
      else if (port["mode"] == "velocity")
      {
        P[i] = port["pressure"];  // Absolute reference for a floating, flow-controlled column.
        double eps = k == 0 ? models[ci].epsilon_node.front() : models[ci].epsilon_node.back();
        double U = port["velocity"].get<double>() * eps;
        if (models[ci].cfg.value("VELOCITY_INPUT_BASIS", std::string("interstitial")) == "superficial")
          U = port["velocity"];
        external[i] = (k == 0 ? 1. : -1.) * models[ci].area * U;
      }
    }
  }
  for (auto& d : connections)
  {
    Edge e;
    e.name = d.at("name");
    e.a = index(d.at("_a"));
    e.b = index(d.at("_b"));
    if (e.a < 0) e.pa = reservoir_pressure(d.at("_a"));
    if (e.b < 0) e.pb = reservoir_pressure(d.at("_b"));
    Json options = d;
    if (d.contains("valve")) options.update(d["valve"]);
    if (d.contains("compressor")) throw std::invalid_argument("Liquid connections use pumps or incompressible valves");
    e.opening = options.value("opening", 1.);
    e.check = options.value("check_valve", false);
    if (options.contains("volume_flow_m3_s"))
      e.imposed = resolve_time_profile(options["volume_flow_m3_s"], t);
    else
    {
      auto model = options.value("flow_model", std::string("linear"));
      if (model != "linear" && model != "incompressible")
        throw std::invalid_argument("Liquid valves require linear/incompressible volumetric conductance");
      double conductance = options.value("volume_conductance", options.value("conductance", 0.));
      if (conductance <= 0. || !std::isfinite(conductance))
        throw std::invalid_argument("Liquid valve conductance must be positive in m3/s/Pa");
      e.linear = 1. / conductance;
    }
    edges.push_back(e);
  }
  for (const auto& e : edges)
  {
    if (e.imposed) continue;
    if (e.a >= 0 && e.b >= 0)
    {
      adjacent[e.a].push_back(e.b);
      adjacent[e.b].push_back(e.a);
    }
    else if (e.a >= 0)
      referenced[e.a] = true;
    else if (e.b >= 0)
      referenced[e.b] = true;
  }
  std::vector<bool> seen(n, false);
  for (std::size_t i = 0; i < n; ++i)
    if (!seen[i])
    {
      std::vector<std::size_t> component = {i};
      seen[i] = true;
      bool reference = false;
      for (std::size_t k = 0; k < component.size(); ++k)
      {
        auto j = component[k];
        reference |= referenced[j];
        for (auto q : adjacent[j])
          if (!seen[q])
          {
            seen[q] = true;
            component.push_back(q);
          }
      }
      if (!reference) fixed[component.front()] = true;
    }
  std::vector<std::size_t> unknown;
  for (std::size_t i = 0; i < n; ++i)
    if (!fixed[i]) unknown.push_back(i);
  auto edge_flow = [&](const Edge& e, const Vec& p)
  {
    double dp = (e.a < 0 ? e.pa : p[e.a]) - (e.b < 0 ? e.pb : p[e.b]);
    if (e.check && dp < 0.) return 0.;
    return e.opening * (e.imposed ? *e.imposed : root(dp, e.linear, e.quadratic));
  };
  auto residual = [&](const Vec& p)
  {
    Vec r = external;
    for (auto& e : edges)
    {
      double q = edge_flow(e, p);
      if (e.a >= 0) r[e.a] -= q;
      if (e.b >= 0) r[e.b] += q;
    }
    return r;
  };
  bool converged = false;
  for (int iteration = 0; iteration < 60; ++iteration)
  {
    Vec full = residual(P), r(unknown.size());
    double norm = 0.;
    for (std::size_t i = 0; i < unknown.size(); ++i)
    {
      r[i] = full[unknown[i]];
      norm = std::max(norm, std::abs(r[i]));
    }
    if (norm < 1e-13)
    {
      converged = true;
      break;
    }
    std::size_t m = unknown.size();
    std::vector<Vec> J(m, Vec(m));
    for (std::size_t j = 0; j < m; ++j)
    {
      Vec trial = P;
      double h = std::max(1e-3, std::abs(P[unknown[j]]) * 1e-6);
      trial[unknown[j]] += h;
      auto rr = residual(trial);
      for (std::size_t i = 0; i < m; ++i) J[i][j] = (rr[unknown[i]] - r[i]) / h;
    }
    for (std::size_t j = 0; j < m; ++j)
    {
      std::size_t pivot = j;
      for (std::size_t k = j + 1; k < m; ++k)
        if (std::abs(J[k][j]) > std::abs(J[pivot][j])) pivot = k;
      if (std::abs(J[pivot][j]) < 1e-30)
        throw std::invalid_argument("Liquid network is hydraulically under-specified or inconsistent");
      std::swap(J[j], J[pivot]);
      std::swap(r[j], r[pivot]);
      for (std::size_t k = j + 1; k < m; ++k)
      {
        double f = J[k][j] / J[j][j];
        for (std::size_t l = j; l < m; ++l) J[k][l] -= f * J[j][l];
        r[k] -= f * r[j];
      }
    }
    Vec delta(m);
    for (std::size_t j = m; j-- > 0;)
    {
      double x = -r[j];
      for (std::size_t k = j + 1; k < m; ++k) x -= J[j][k] * delta[k];
      delta[j] = x / J[j][j];
    }
    double factor = 1.;
    bool accepted = false;
    for (int k = 0; k < 25; ++k)
    {
      Vec trial = P;
      for (std::size_t j = 0; j < m; ++j) trial[unknown[j]] += factor * delta[j];
      auto rr = residual(trial);
      double next = 0.;
      for (auto j : unknown) next = std::max(next, std::abs(rr[j]));
      if (next < norm)
      {
        P = trial;
        accepted = true;
        break;
      }
      factor *= .5;
    }
    if (!accepted) break;
  }
  if (!converged) throw std::runtime_error("Incompressible hydraulic network failed to converge");
  for (double pressure : P)
    if (!std::isfinite(pressure) || pressure < 0.)
      throw std::invalid_argument("Liquid boundary settings require a negative absolute pressure");
  // Check all unanchored flow constraints, including the equation omitted for a floating pressure reference.
  auto all = residual(P);
  for (std::size_t i = 0; i < n; ++i)
    if (!referenced[i] && std::abs(all[i]) > 1e-11)
      throw std::invalid_argument("Incompressible volume balance is inconsistent with imposed flows");
  Json out = {{"columns", Json::object()}, {"connections", Json::object()}};
  for (std::size_t ci = 0; ci < models.size(); ++ci)
  {
    double q = edge_flow(edges[ci], P);
    out["columns"][column_names[ci]] = {{"left", {{"volume_flux", q / models[ci].area}, {"pressure", P[2 * ci]}}},
                                        {"right", {{"volume_flux", q / models[ci].area}, {"pressure", P[2 * ci + 1]}}}};
  }
  for (std::size_t i = models.size(); i < edges.size(); ++i) out["connections"][edges[i].name] = edge_flow(edges[i], P);
  return out;
}
}  // namespace cycle

#include "cycle/diagnostics.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace cycle {
Json sample_statistics(const Vec& values) {
  Vec v;
  for (double x : values) if (std::isfinite(x)) v.push_back(x);
  if (v.empty()) return {{"count", 0}, {"min", nullptr}, {"max", nullptr},
                         {"mean", nullptr}, {"median", nullptr}, {"std", nullptr}};
  double mean = std::accumulate(v.begin(), v.end(), 0.) / v.size(), variance = 0.;
  for (double x : v) variance += (x - mean) * (x - mean);
  std::sort(v.begin(), v.end());
  return {{"count", v.size()}, {"min", v.front()}, {"max", v.back()}, {"mean", mean},
          {"median", (v[(v.size()-1)/2] + v[v.size()/2]) / 2.},
          {"std", std::sqrt(variance / v.size())}, {"peak_to_peak", v.back() - v.front()}};
}

Json assess_convergence(const Json& rows, const Json& cfg) {
  if (!rows.is_array()) throw std::invalid_argument("Cycle diagnostics require an array");
  for (const char* key : {"CSS_DIAGNOSTIC_WINDOW", "CSS_DIAGNOSTIC_MINIMUM_CYCLES", "CSS_MAXIMUM_PERIOD"})
    if (cfg.contains(key) && (!cfg[key].is_number() || !std::isfinite(cfg[key].get<double>()) ||
        cfg[key].get<double>() != std::floor(cfg[key].get<double>())))
      throw std::invalid_argument(std::string(key) + " must be an integer");
  const int window = cfg.value("CSS_DIAGNOSTIC_WINDOW", 20);
  const int minimum = cfg.value("CSS_DIAGNOSTIC_MINIMUM_CYCLES", 6);
  const int max_period = cfg.value("CSS_MAXIMUM_PERIOD", 4);
  const std::string selection = cfg.value("CSS_PLOT_CYCLE_SELECTION", std::string("minimum_error"));
  if (window < 2 || minimum < 2 || minimum > window || max_period < 1 || max_period > 20 ||
      (selection != "minimum_error" && selection != "last"))
    throw std::invalid_argument("Invalid CSS diagnostic window, minimum cycles, maximum period or plot selection");
  Json out = {{"status", "insufficient_history"}, {"strict_css", false},
              {"requested_window_cycles", window}, {"minimum_diagnostic_cycles", minimum},
              {"state_error_pass_threshold", 1.0},
              {"state_error_definition", "max(abs(end-start)/(state_atol + state_rtol*max(abs(start),abs(end))))"},
              {"oscillation_detected", false}, {"detected_period_cycles", nullptr},
              {"physical_periodicity_verified", false}, {"representative_cycle", nullptr},
              {"primary_ranking_eligible", false},
              {"statistics_note", "Final-window descriptive spread, not statistical confidence intervals or numerical error bounds."}};
  if (rows.empty()) return out;
  const std::size_t first = rows.size() > std::size_t(window) ? rows.size() - window : 0;
  const std::size_t count = rows.size() - first;
  out["window_start_cycle"] = rows[first].at("cycle");
  out["window_end_cycle"] = rows.back().at("cycle");
  out["cycles_evaluated"] = count;
  out["window_complete"] = count == std::size_t(window);
  out["history_sufficient"] = count >= std::size_t(minimum);
  bool strict = rows.back().value("converged", false), balances = true;
  std::set<std::string> keys;
  for (std::size_t i = first; i < rows.size(); ++i)
    for (const auto& [key, value] : rows[i].items())
      if (value.is_number() && key != "cycle" && key != "css_streak") keys.insert(key);
  Json errors = Json::object(), kpis = Json::object();
  Vec e;
  std::size_t chosen = rows.size()-1;
  double best = std::numeric_limits<double>::infinity();
  for (std::size_t i = first; i < rows.size(); ++i) {
    const auto v = rows[i].value("css_state_error", Json());
    if (v.is_number() && std::isfinite(v.get<double>())) {
      e.push_back(v.get<double>());
      if (v.get<double>() <= best) { best = v; chosen = i; } // latest tie
    }
    auto b = rows[i].value("mass_balance_error", Json());
    balances &= b.is_number() && std::isfinite(b.get<double>()) &&
                b.get<double>() <= cfg.value("CSS_MASS_BALANCE_TOLERANCE", 1e-3);
  }
  for (const auto& key : keys) {
    Vec values;
    for (std::size_t i = first; i < rows.size(); ++i) {
      const auto v = rows[i].value(key, Json());
      if (v.is_number()) values.push_back(v.get<double>());
    }
    auto stats = sample_statistics(values);
    stats["missing_cycles"] = count - stats["count"].get<std::size_t>();
    (key.starts_with("css_") ? errors : kpis)[key] = stats;
  }
  out["error_statistics"] = errors;
  out["kpi_statistics"] = kpis;
  out["window_mass_balance_checks_passed"] = balances;
  int turns = 0, sign = 0;
  if (e.size() == count && count >= std::size_t(minimum)) {
    const auto [lo, hi] = std::minmax_element(e.begin(), e.end());
    double floor = std::max(1e-12, 1e-3 * (*hi - *lo));
    for (std::size_t i = 1; i < e.size(); ++i) {
      double diff = e[i] - e[i-1];
      int next = diff > floor ? 1 : diff < -floor ? -1 : 0;
      if (next) { if (sign && sign != next) ++turns; sign = next; }
    }
    const auto half = e.size()/2;
    double a = std::accumulate(e.begin(), e.begin()+half, 0.) / half;
    double b = std::accumulate(e.begin()+half, e.end(), 0.) / (e.size()-half);
    out["tail_error_mean_relative_shift"] = (b-a) / std::max({std::abs(a), std::abs(b), 1e-12});
    out["oscillation_detected"] = turns >= 2;
  }
  out["error_direction_changes"] = turns;
  int period = 0;
  Json lags = Json::object();
  for (int p = 2; p <= max_period; ++p) {
    const auto key = std::to_string(p);
    int streak = 0;
    Vec values;
    for (std::size_t i = first; i < rows.size(); ++i) {
      Json l = rows[i].value("css_period_errors", Json::object()).value(key, Json::object());
      streak = l.value("criteria_satisfied", false) ? streak+1 : 0;
      if (l.contains("css_state_error") && l["css_state_error"].is_number())
        values.push_back(l["css_state_error"]);
    }
    lags[key] = {{"state_error", sample_statistics(values)}, {"final_passing_streak", streak}};
    if (!period && count >= std::size_t(minimum) &&
        streak >= std::max(p, cfg.value("CSS_CONSECUTIVE_CYCLES", 3))) period = p;
  }
  out["period_diagnostics"] = lags;
  std::string status = strict ? "strict_css" : count < std::size_t(minimum) ? "insufficient_history"
      : period ? "periodic_candidate" : turns >= 2 ? "oscillatory_nonconverged" : "nonconverged";
  if (period && !strict) out["detected_period_cycles"] = period;
  if (strict || selection == "last") chosen = rows.size()-1;
  // A repeated multi-cycle orbit is represented by a complete final period,
  // never by the lowest one-cycle residual alone.
  if (period && !strict) chosen = rows.size()-1;
  Json plot_cycles = Json::array();
  if (period && !strict) {
    for (std::size_t i = rows.size()-period; i < rows.size(); ++i) plot_cycles.push_back(rows[i]["cycle"]);
  } else plot_cycles.push_back(rows[chosen]["cycle"]);
  out["status"] = status;
  out["strict_css"] = strict;
  out["representative_cycle"] = rows[chosen]["cycle"];
  out["plot_cycles"] = plot_cycles;
  out["selection_reason"] = strict ? "last strict-CSS cycle" : period ? "complete candidate period"
      : selection == "last" ? "last completed cycle (not CSS)" : "minimum-residual cycle in final window (not CSS)";
  out["primary_ranking_eligible"] = strict && balances && rows.back().value("process_feasible", false);
  out["warning"] = strict ? "" : "Strict CSS was not reached. Selected plots and tail-window metrics are diagnostic; numerical refinement is needed before physical or ranking claims.";
  return out;
}
}

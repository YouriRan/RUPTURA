#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "cycle/npz.hpp"
#include "cycle/runner.hpp"
namespace cycle
{
Json CycleRunner::capture_liquid_pressures(double t, const Vec& augmented)
{
  if (!network.models.front().liquid) return Json();
  auto physical = Vec(augmented.begin(), augmented.begin() + physical_size);
  network.evaluate(t, physical);
  Json record = Json::object();
  for (const auto& name : network.column_names)
    record[name] = network.model(name).derived(network.column_state(name, physical)).P;
  return record;
}
void CycleRunner::trim_history()
{
  const int first = progress.at("cycle").get<int>() - history_cycles + 1;
  std::size_t n = 0;
  while (n < progress["history_cycle"].size() && progress["history_cycle"][n].get<int>() < first) ++n;
  for (const auto* key : {"history_t", "history_cycle", "history_phase", "history_state", "history_pressure"})
    if (n) progress[key].erase(progress[key].begin(), progress[key].begin() + n);
}
void CycleRunner::restore_history(Json& p, const std::filesystem::path& archive, bool committed)
{
  for (const auto* key : {"history_t", "history_cycle", "history_phase", "history_state", "history_pressure"})
    p[key] = Json::array();
  if (std::filesystem::exists(archive))
  {
    Json pressure = Json::object();
    if (network.models.front().liquid)
      for (const auto& name : network.column_names)
      {
        try
        {
          auto a = read_npz_array(archive, name + "_P");
          pressure[name] = Json::array();
          for (std::size_t k = 0; k < a.shape.at(0); ++k)
            pressure[name].push_back(
                Vec(a.data.begin() + k * network.model(name).Nz, a.data.begin() + (k + 1) * network.model(name).Nz));
        }
        catch (const std::exception&)
        {
          throw std::runtime_error(
              "Liquid restart requires recorded pressure profiles; use a checkpoint from this integrated build");
        }
      }
    auto state = read_npz_array(archive, "state"), time = read_npz_array(archive, "time"),
         cycles = read_npz_array(archive, "cycle"), indices = read_npz_array(archive, "phase");
    const auto n = time.data.size();
    if (state.shape != std::vector<std::size_t>{n, physical_size} || cycles.data.size() != n ||
        indices.data.size() != n)
      throw std::runtime_error("Profile archive dimensions incompatible");
    const auto last = p.value("history_last", Json::array());
    for (std::size_t i = 0; i < n; ++i)
    {
      const double t = time.data[i];
      const int c = int(cycles.data[i]), ph = int(indices.data[i]);
      if (!std::isfinite(t) || c < 1 || ph < 0 || ph >= int(phases.size()))
        throw std::runtime_error("Invalid profile coordinates");
      if (committed && !last.empty() &&
          (c > last[1].get<int>() || (c == last[1].get<int>() && ph > last[2].get<int>()) || t > last[0].get<double>()))
        continue;
      p["history_t"].push_back(t);
      p["history_cycle"].push_back(c);
      p["history_phase"].push_back(ph);
      p["history_state"].push_back(
          Vec(state.data.begin() + i * physical_size, state.data.begin() + (i + 1) * physical_size));
      Json record = Json::object();
      for (auto& [name, rows] : pressure.items()) record[name] = rows[i];
      p["history_pressure"].push_back(record);
    }
  }
  if (p["history_t"].empty())
  {
    if (!p.contains("state")) throw std::runtime_error("No profiles available for plotting");
    std::cerr << "Profile history is unavailable; numerical restart is retained, but older profile plots cannot be "
                 "recovered.\n";
    p["profile_history_incomplete"] = true;
    p["history_t"].push_back(p["global_time"]);
    p["history_cycle"].push_back(p["cycle"]);
    p["history_phase"].push_back(std::min(p["phase"].get<int>(), int(phases.size()) - 1));
    auto state = p["state"].get<Vec>();
    state.resize(physical_size);
    p["history_state"].push_back(state);
    p["history_pressure"].push_back(Json::object());
  }
}
void CycleRunner::append_trends(const Vec& times, const std::vector<Vec>& states, int cycle, int phase)
{
  auto path = output_dir / "history_trends.csv";
  const bool header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
  std::ofstream f(path, std::ios::app);
  f.exceptions(std::ios::badbit | std::ios::failbit);
  f << std::setprecision(17);
  if (header)
  {
    f << "time,cycle,phase";
    for (const auto& name : network.column_names)
    {
      for (const auto* field :
           {"P", "Tg", "Ts", "Tw", "Tj", "htf_inlet_temperature", "htf_outlet_temperature", "htf_active"})
        f << ',' << name << '_' << field;
      for (const auto* end : {"left", "right"})
        for (std::size_t j = 0; j < network.n_comp; ++j) f << ',' << name << "_c_" << end << '_' << j;
    }
    f << '\n';
  }
  auto mean = [](const Vec& v) { return std::accumulate(v.begin(), v.end(), 0.) / v.size(); };
  for (std::size_t i = 0; i < times.size(); ++i)
  {
    f << progress["phase_start_time"].get<double>() + times[i] << ',' << cycle << ',' << phase;
    if (network.models.front().liquid)
    {
      auto physical = Vec(states[i].begin(), states[i].begin() + physical_size);
      network.evaluate(times[i], physical);
    }
    for (const auto& name : network.column_names)
    {
      auto local = network.column_state(name, states[i]);
      const auto& model = network.model(name);
      auto d = model.derived(local);
      Column recorder(model);
      recorder.set_step({{"left_port", {{"mode", "closed"}}},
                         {"right_port", {{"mode", "closed"}}},
                         {"thermal", phases[phase]["columns"][name].value("thermal", Json::object())}});
      auto jacket = recorder.jacket_record(times[i], local);
      auto tj =
          jacket.value("temperature_profile", Json(Vec(model.Nz, jacket["temperature"].get<double>()))).get<Vec>();
      f << ',' << mean(d.P) << ',' << mean(d.Tg) << ',' << mean(d.Ts) << ',' << mean(d.Tw) << ',' << mean(tj) << ','
        << jacket["temperature"].get<double>() << ','
        << jacket.value("outlet_temperature", jacket["temperature"]).get<double>() << ','
        << (jacket["mode"] == "htf" ? 1 : 0);
      for (auto end : {std::size_t(0), model.Nz - 1})
        for (std::size_t j = 0; j < network.n_comp; ++j) f << ',' << d.c[j * model.Nz + end];
    }
    f << '\n';
  }
  f.flush();
}
void CycleRunner::load_plot_progress()
{
  const auto path = output_dir / "checkpoint_latest.json";
  if (std::filesystem::exists(path))
  {
    std::ifstream f(path);
    f >> progress;
    if (!progress.contains("history_state")) restore_history(progress, output_dir / "states.npz", true);
    if (progress.value("summary", Json()).is_null())
    {
      const auto& rows = progress.at("cycle_rows");
      progress["summary"] = {{"completed", false},
                             {"converged", !rows.empty() && rows.back().value("converged", false)},
                             {"cycles", rows.size()},
                             {"termination_reason", "plot_snapshot_of_incomplete_run"},
                             {"output_dir", output_dir.string()},
                             {"final_metrics", rows.empty() ? Json::object() : rows.back()}};
      // An extended/interrupted run must not reuse an older final report whose
      // representative cycle may already have left the rolling profile window.
      if (!rows.empty())
        update_convergence_report();
      else
        std::filesystem::remove(output_dir / "convergence_report.json");
    }
    // Also makes an interrupted run's raw recording readable by the regular plotter.
    write_outputs(true);
  }
  else
  {
    progress = Json::object();
    restore_history(progress, output_dir / "states.npz", false);
    std::ifstream f(output_dir / "summary.json");
    if (f) f >> progress["summary"];
    progress["cycle_rows"] = Json::array();
  }
}
}  // namespace cycle

namespace cycle
{
namespace
{
std::string safe_file_label(std::string label)
{
  for (char& c : label)
    if (!std::isalnum(static_cast<unsigned char>(c)) && c != '-' && c != '_') c = '_';
  return label;
}
void output_header(std::ofstream& f, const std::vector<std::string>& labels)
{
  for (std::size_t i = 0; i < labels.size(); ++i) f << "# column " << i + 1 << ": " << labels[i] << '\n';
}
}  // namespace
void CycleRunner::append_main_frames(const Vec& times, const std::vector<Vec>& states)
{
  const std::vector<std::string> ch = {"Dimensionless time, tau = tv/L [-]",
                                       "Time, t [min]",
                                       "Column position, z [m]",
                                       "Interstitial fluid velocity, v [m/s]",
                                       "Total pressure, p_t [Pa]",
                                       "Fluid temperature, T_g [K]",
                                       "Fluid temperature time derivative, dT_g/dt [K/s]",
                                       "Solid temperature, T_s [K]",
                                       "Solid temperature time derivative, dT_s/dt [K/s]",
                                       "Wall temperature, T_w [K]",
                                       "Wall temperature time derivative, dT_w/dt [K/s]",
                                       "Fluid density [kg/m^3]"};
  const std::vector<std::string> bh = {"Dimensionless time, tau = tv/L [-]",
                                       "Time, t [min]",
                                       "Column position, z [m]",
                                       "Concentration, c_i [mol/m^3]",
                                       "Concentration time derivative, dc_i/dt [mol/m^3/s]",
                                       "Mole fraction, y_i [-]",
                                       "Physisorption, q_phy_i [mol/kg]",
                                       "Physisorption time derivative, dq_phy_i/dt [mol/kg/s]",
                                       "Chemisorption, q_chem_i [mol/kg]",
                                       "Chemisorption time derivative, dq_chem_i/dt [mol/kg/s]",
                                       "Partial pressure or liquid concentration driving force",
                                       "Equilibrium physisorption, q_phy_i^* [mol/kg]",
                                       "Normalized outlet composition or concentration [-]",
                                       "Equilibrium chemisorption (not recorded)"};
  for (const auto& name : network.column_names)
  {
    auto directory = output_dir / "main" / safe_file_label(name);
    std::filesystem::create_directories(directory);
    auto column_path = directory / "column.data";
    bool first = !std::filesystem::exists(column_path);
    std::ofstream column(column_path, std::ios::app);
    column.exceptions(std::ios::badbit | std::ios::failbit);
    column << std::setprecision(17);
    if (first) output_header(column, ch);
    std::vector<std::ofstream> component;
    for (std::size_t j = 0; j < network.n_comp; ++j)
    {
      auto path =
          directory / ("component_" + std::to_string(j) + "_" + safe_file_label(network.comp_names[j]) + ".data");
      bool fresh = !std::filesystem::exists(path);
      component.emplace_back(path, std::ios::app);
      component.back().exceptions(std::ios::badbit | std::ios::failbit);
      component.back() << std::setprecision(17);
      if (fresh)
      {
        component.back() << "# component index: " << j << "\n# component name: " << network.comp_names[j] << '\n';
        output_header(component.back(), bh);
      }
    }
    auto& m = network.model(name);
    const auto start = network.column_starts.at(name);
    double reference_v = m.cfg.value("u_ads", m.cfg.value("u_in", .01));
    Vec feed = m.cfg.at("y_feed").get<Vec>();
    if (m.liquid) feed = m.cfg.at("LIQUID_FEED_CONCENTRATION").get<Vec>();
    for (std::size_t sample = 0; sample < times.size(); ++sample)
    {
      Vec physical(states[sample].begin(), states[sample].begin() + physical_size);
      auto [all_deriv, diag] = network.evaluate(times[sample], physical);
      auto local = network.column_state(name, physical), deriv = network.column_state(name, all_deriv);
      auto vf = m.face_velocities(local, times[sample], diag["boundary_overrides"][name]);
      auto d = m.derived(local);
      auto [qphy, qchem] = m.loading_profiles(local);
      auto [dqphy, dqchem] = m.loading_profiles(deriv);
      auto eq = m.equilibrium_loading(local);
      double time = progress.value("phase_start_time", 0.) + times[sample], tau = time * reference_v / m.L;
      for (std::size_t i = 0; i < m.Nz; ++i)
      {
        double v = .5 * (vf[i] * m.epsilon_face[i] + vf[i + 1] * m.epsilon_face[i + 1]) / m.epsilon_node[i];
        column << tau << ' ' << time / 60. << ' ' << m.z[i] << ' ' << v << ' ' << d.P[i] << ' ' << d.Tg[i] << ' '
               << deriv[local.size() - 3 * m.Nz + i] << ' ' << d.Ts[i] << ' ' << deriv[local.size() - 2 * m.Nz + i]
               << ' ' << d.Tw[i] << ' ' << deriv[local.size() - m.Nz + i] << ' ' << d.rho[i] << '\n';
        for (std::size_t j = 0; j < m.n_comp; ++j)
        {
          auto ix = j * m.Nz + i;
          double ratio = feed[j] > 0 ? (m.liquid ? d.c[ix] : d.y[ix]) / feed[j] : 0.;
          component[j] << tau << ' ' << time / 60. << ' ' << m.z[i] << ' ' << d.c[ix] << ' ' << deriv[ix] << ' '
                       << d.y[ix] << ' ' << qphy[ix] << ' ' << dqphy[ix] << ' ' << qchem[ix] << ' ' << dqchem[ix] << ' '
                       << (m.liquid ? d.c[ix] : d.y[ix] * d.P[i]) << ' ' << eq[ix] << ' ' << ratio
                       << (m.layout.n_loading_groups == m.material_names.size() ? " 0\n" : " nan\n");
        }
      }
      column << '\n';
      for (auto& f : component) f << '\n';
    }
  }
}
}  // namespace cycle

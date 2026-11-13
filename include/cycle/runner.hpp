#pragma once
#include <filesystem>
#include <functional>
#include <map>
#include <optional>

#include "cycle/integrator.hpp"
#include "cycle/network.hpp"
namespace cycle
{
struct StageRule
{
  std::optional<double> duration, maximum;
  Json stop;
  double minimum = 0, failure_time = 2e5;
  std::string timeout_policy = "error";
  static StageRule from_mapping(const Json&, double default_failure = 2e5);
  bool enabled() const { return duration.has_value() || !stop.is_null(); }
  double next_boundary(double elapsed) const;
};
struct PhaseSolution
{
  Vec times, final_augmented, ledger;
  std::vector<Vec> augmented;
  std::string reason;
  Json action_results = Json::object();
  BdfStatistics statistics;
};
using PhaseCommit = std::function<void(double, const Vec&, const Json&, const Vec&, const std::vector<Vec>&)>;
class CycleRunner
{
 public:
  explicit CycleRunner(const Json& config);
  Json cfg, phases, progress;
  std::string solver_method = "BDF";
  // auto: Matplotlib when it is usable, otherwise the built-in SVG overview.
  std::string plot_backend = "auto";
  ColumnNetwork network;
  size_t physical_size = 0, ledger_size = 0;
  std::vector<std::string> tags, utility_names;
  std::filesystem::path output_dir, checkpoint_path;
  Json validate();
  Json run();
  PhaseSolution integrate_phase(const Vec& physical_state, const Json& phase, Json context = Json::object(),
                                Vec augmented = {}, double elapsed = 0, Json done = Json::object(),
                                Vec phase_ledger_start = {}, const PhaseCommit& commit = {});
  Json cycle_metrics(const Vec& start, const Vec& end, const Vec& delta, double duration);
  std::pair<Json, bool> css(const Vec& start, const Vec& end, const Json& metrics, const Json& previous_metrics,
                            const Json& compositions, const Json& previous_compositions);
  Vec ledger_rates(const Json&) const;
  Vec rhs(double, const Vec&);
  double criterion_residual(const Json&, double, const Vec&, const std::string& column = "");
  Json stream_record(const Vec& delta) const;
  Vec ledger_value(const Vec&, const std::vector<std::string>& key) const;
  void validate_state(const Vec&) const;
  void write_plots();
  // Individual backends. Matplotlib runs the bundled tools/plot_results.py on
  // the saved NPZ/CSV results; it never integrates anything.
  bool write_matplotlib_plots(std::string* message = nullptr);
  void write_svg_plots();
  void update_convergence_report();
  void load_plot_progress();

 private:
  std::map<std::vector<std::string>, std::pair<size_t, size_t>> ledger_slices;
  Vec atol, phase_ledger_start_;
  Json current_phase_, active_context_, restart_config_;
  double rtol = 1e-6, segment_time = 5;
  int cycle_limit = 1, required_streak = 3;
  bool restart = false;
  bool compact_output = false, checkpoint_enabled = true, phase_archives = true, selected_archive = true,
       compress_output = false;
  double output_interval = 5;
  int output_samples = 200, history_cycles = 22;
  std::string archive_stamp;
  Json capture_liquid_pressures(double, const Vec&);
  void trim_history();
  void restore_history(Json&, const std::filesystem::path&, bool trim_to_checkpoint);
  void append_main_frames(const Vec&, const std::vector<Vec>&);
  void append_trends(const Vec&, const std::vector<Vec>&, int, int);
  std::string run_mode = "css", fingerprint;
  std::filesystem::path restart_path;
  void canonicalize_phases();
  void init_ledger();
  void validate_phases();
  void validate_criterion(const Json&) const;
  void activate_phase(const Json&, const Json&, const Json&);
  static bool needs_pressure_anchor(const Json&);
  Json criterion_context(double, const Vec&, const std::string& column = "") const;
  double criterion_actual(const Json&, double, const Vec&, const std::string&);
  bool criterion_satisfied_roundoff(const Json&, double, const Vec&, const std::string&);
  std::string rule_status(const StageRule&, double, const Vec&, const std::string& column = "");
  Json fresh_progress();
  void save_checkpoint();
  Json load_checkpoint();
  void prepare_output();
  void commit_phase(const Json&, const PhaseSolution&);
  void write_outputs(bool final = false);
  Json output_arrays(int cycle = -1, int phase = -1) const;
  void write_csv(const std::string&, const Json&) const;
};
using MultiColumnCycleRunner = CycleRunner;
}  // namespace cycle

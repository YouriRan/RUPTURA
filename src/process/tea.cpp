#include "cycle/tea.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
namespace cycle {
namespace {
constexpr double nan = std::numeric_limits<double>::quiet_NaN();
double number(const Json& o, const char* k, double d = 0.) {
  return o.contains(k) ? o.at(k).get<double>() : d;
}
bool truth(const Json& v) {
  if (v.is_null()) return false;
  if (v.is_boolean()) return v.get<bool>();
  if (v.is_number()) return v.get<double>() != 0.;
  if (v.is_string()) return !v.get<std::string>().empty();
  return !v.empty();
}
bool flag(const Json& o, const char* k, bool d = false) {
  return o.contains(k) ? truth(o.at(k)) : d;
}
std::string lower(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
std::string upper(std::string s) {
  for (auto& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}
Json clean(Json v) {
  if (v.is_number_float() && !std::isfinite(v.get<double>())) return nullptr;
  if (v.is_object() || v.is_array())
    for (auto& x : v) x = clean(x);
  return v;
}
std::string csv_text(const Json& v) {
  if (v.is_null()) return {};
  std::string raw = v.is_string()    ? v.get<std::string>()
                    : v.is_boolean() ? (v.get<bool>() ? "True" : "False")
                                     : v.dump();
  if (raw.find_first_of(",\"\r\n") == std::string::npos) return raw;
  std::string out = "\"";
  for (char c : raw) {
    out += c;
    if (c == '"') out += '"';
  }
  return out + '"';
}
}  // namespace
double compression_work(double moles, double temperature, double inlet_pressure,
                        double outlet_pressure, double efficiency, double gas_constant,
                        const std::string& model, double gamma) {
  for (double x : {moles, temperature, inlet_pressure, outlet_pressure, efficiency, gas_constant})
    if (!std::isfinite(x)) throw std::invalid_argument("Compression inputs must be finite");
  if (moles < 0 || std::min({temperature, inlet_pressure, outlet_pressure, gas_constant}) <= 0)
    throw std::invalid_argument("Compression requires nonnegative moles and positive T/P/R");
  if (!(efficiency > 0 && efficiency <= 1))
    throw std::invalid_argument("Compression efficiency must lie in (0, 1]");
  if (moles == 0 || outlet_pressure <= inlet_pressure) return 0.;
  double ratio = outlet_pressure / inlet_pressure;
  std::string selected = lower(model);
  auto first = selected.find_first_not_of(" \t\r\n");
  selected = first == std::string::npos
                 ? std::string()
                 : selected.substr(first, selected.find_last_not_of(" \t\r\n") - first + 1);
  std::replace(selected.begin(), selected.end(), '-', '_');
  double reversible;
  if (selected == "isothermal" || selected == "reversible_isothermal")
    reversible = moles * gas_constant * temperature * std::log(ratio);
  else if (selected == "isentropic" || selected == "adiabatic") {
    if (gamma <= 1) throw std::invalid_argument("COMPRESSION_HEAT_CAPACITY_RATIO must exceed one");
    reversible = moles * gas_constant * temperature * gamma / (gamma - 1.) *
                 (std::pow(ratio, (gamma - 1.) / gamma) - 1.);
  } else
    throw std::invalid_argument("COMPRESSION_WORK_MODEL must be 'isothermal' or 'isentropic'");
  return reversible / std::max(efficiency, 1e-12);
}
Json build_tea_record(const Json& cfg, const Json& summary) {
  const Json metrics = summary.value("reporting_metrics", summary.value("final_metrics", summary.value("metrics", Json::object())));
  if (!metrics.is_object()) throw std::invalid_argument("Cycle metrics must be a mapping");
  double r = number(cfg, "R", 8.31446261815324), temperature = number(cfg, "T_feed", 298.15),
         feed_moles = number(metrics, "feed_mol");
  double heavy_moles = std::max(
      number(metrics, "gross_heavy_product_mol", number(metrics, "heavy_product_mol")), 0.);
  double high_pressure =
      number(cfg, "P_HIGH", number(metrics, "high_pressure_Pa", number(cfg, "P_init", 1e5)));
  double low_pressure =
      number(cfg, "P_LOW", number(metrics, "low_pressure_Pa", number(cfg, "P_init", 1e5)));
  double supply_pressure = number(cfg, "FEED_SUPPLY_PRESSURE", high_pressure),
         ambient_pressure = number(cfg, "AMBIENT_PRESSURE", 1e5);
  double compressor_efficiency = number(cfg, "COMPRESSOR_ISENTROPIC_EFFICIENCY", .75),
         vacuum_efficiency = number(cfg, "VACUUM_PUMP_EFFICIENCY", .60);
  if (!(compressor_efficiency > 0 && compressor_efficiency <= 1 && vacuum_efficiency > 0 &&
        vacuum_efficiency <= 1))
    throw std::invalid_argument("Compression and vacuum efficiencies must lie in (0, 1]");
  std::string work_model = cfg.value("COMPRESSION_WORK_MODEL", std::string("isothermal"));
  double gamma = number(cfg, "COMPRESSION_HEAT_CAPACITY_RATIO", 1.4);
  auto work = [&](double mol, double pin, double pout, double eff) {
    return compression_work(mol, temperature, pin, pout, eff, r, work_model, gamma);
  };
  double feed_work = flag(cfg, "CHARGE_FEED_COMPRESSION")
                         ? work(feed_moles, supply_pressure, high_pressure, compressor_efficiency)
                         : 0.;
  double ambient_feed_work =
             work(feed_moles, ambient_pressure, high_pressure, compressor_efficiency),
         vacuum_work = work(heavy_moles, low_pressure, ambient_pressure, vacuum_efficiency);
  feed_work = number(metrics, "feed_compression_work_J_cycle", feed_work);
  vacuum_work = number(metrics, "vacuum_work_J_cycle", vacuum_work);
  double thermal_heat = number(metrics, "thermal_heat_input_J_cycle"),
         cooling_duty = number(metrics, "cooling_duty_J_cycle"),
         electric_heater = number(metrics, "electric_heater_work_J_cycle");
  double gas_heating = number(metrics, "purge_and_feed_gas_heating_J_cycle"),
         gas_cooling = number(metrics, "purge_and_feed_gas_cooling_J_cycle"),
         ambient_heat_loss = number(metrics, "ambient_heat_loss_J_cycle");
  double recovery_fraction = number(cfg, "HEAT_RECOVERY_FRACTION");
  if (!(recovery_fraction >= 0 && recovery_fraction <= 1))
    throw std::invalid_argument("HEAT_RECOVERY_FRACTION must lie in [0, 1]");
  double non_electric_heat = std::max(thermal_heat - electric_heater, 0.),
         recoverable_heat = non_electric_heat + gas_heating,
         total_cooling = cooling_duty + gas_cooling;
  double recovered_heat = recovery_fraction * std::min(recoverable_heat, total_cooling),
         heating_efficiency = number(cfg, "HEATING_UTILITY_EFFICIENCY", .90),
         cooling_cop = number(cfg, "COOLING_COP", 3.);
  if (heating_efficiency <= 0 || cooling_cop <= 0)
    throw std::invalid_argument("HEATING_UTILITY_EFFICIENCY and COOLING_COP must be positive");
  double external_cooling = std::max(total_cooling - recovered_heat, 0.);
  std::string cooling_basis = lower(cfg.value("COOLING_COST_BASIS", std::string("electricity")));
  if (cooling_basis != "electricity" && cooling_basis != "purchased_cooling")
    throw std::invalid_argument("COOLING_COST_BASIS must be 'electricity' or 'purchased_cooling'");
  double cooling_electricity = cooling_basis == "electricity" ? external_cooling / cooling_cop : 0.;
  double light_total = std::max(number(metrics, "light_product_mol"), 0.),
         delivery_pressure = number(cfg, "PRODUCT_DELIVERY_PRESSURE", high_pressure),
         buffer_pressure = number(cfg, "PRODUCT_BUFFER_PRESSURE", high_pressure);
  double product_work =
             work(light_total, buffer_pressure, delivery_pressure, compressor_efficiency),
         recycle_moles = number(metrics, "product_buffer_consumption_mol");
  double recycle_work =
      flag(cfg, "CHARGE_RECYCLE_COMPRESSION") &&
              !metrics.contains("network_compression_work_J_cycle")
          ? work(recycle_moles, number(cfg, "RECYCLE_SUCTION_PRESSURE", low_pressure),
                 number(cfg, "RECYCLE_DISCHARGE_PRESSURE", high_pressure), compressor_efficiency)
          : 0.;
  product_work = number(metrics, "product_compression_work_J_cycle", product_work);
  recycle_work = number(metrics, "recycle_compression_work_J_cycle", recycle_work);
  double network_work = number(metrics, "network_compression_work_J_cycle");
  for (double value : {feed_work, vacuum_work, product_work, recycle_work, network_work})
    if (!std::isfinite(value) || value < 0)
      throw std::invalid_argument("Measured machine work integrals must be finite and nonnegative");
  double compression_and_vacuum_work =
             feed_work + vacuum_work + product_work + recycle_work + network_work,
         total_electric = compression_and_vacuum_work + electric_heater + cooling_electricity,
         total_work = feed_work + vacuum_work;
  auto names = cfg.at("comp_names").get<std::vector<std::string>>();
  if (names.empty()) throw std::invalid_argument("comp_names cannot be empty");
  std::string default_target = std::find(names.begin(), names.end(), "CH4") != names.end()
                                   ? "CH4"
                                   : names.back(),
              target = cfg.value("TARGET_LIGHT_COMPONENT", default_target);
  auto pos = std::find(names.begin(), names.end(), target);
  if (pos == names.end())
    throw std::invalid_argument("TARGET_LIGHT_COMPONENT is absent from comp_names");
  auto index = static_cast<std::size_t>(pos - names.begin());
  bool product_available = flag(metrics, "light_product_available", true);
  double target_moles;
  std::string amount_basis;
  if (metrics.contains("net_light_product_component_mol") &&
      !metrics.at("net_light_product_component_mol").is_null()) {
    auto amounts = metrics.at("net_light_product_component_mol").get<Vec>();
    if (amounts.size() != names.size())
      throw std::invalid_argument("net_light_product_component_mol has invalid shape");
    for (double v : amounts) {
      if (!std::isfinite(v))
        throw std::invalid_argument(
            "net_light_product_component_mol must contain finite component amounts");
      if (v < -1e-9) product_available = false;
    }
    target_moles = product_available ? std::max(amounts[index], 0.) : nan;
    amount_basis = "measured net external component export";
  } else if (metrics.contains("light_product_mol") && metrics.contains("light_product_purity")) {
    target_moles =
        product_available && !metrics.at("light_product_purity").is_null()
            ? number(metrics, "light_product_mol") * number(metrics, "light_product_purity")
            : nan;
    amount_basis = "net product total times its measured target purity";
  } else {
    double feed_target =
        metrics.contains("feed_component_mol") && !metrics.at("feed_component_mol").is_null()
            ? metrics.at("feed_component_mol").at(index).get<double>()
            : feed_moles * cfg.at("y_feed").at(index).get<double>();
    auto rec = metrics.value("light_component_recovery", Json(0.));
    target_moles = !rec.is_null() && product_available ? rec.get<double>() * feed_target : nan;
    amount_basis = "legacy recovery-based reconstruction";
  }
  double target_kg = target_moles * cfg.at("M").at(index).get<double>(),
         standard_t = number(cfg, "STANDARD_TEMPERATURE_K", 273.15),
         standard_p = number(cfg, "STANDARD_PRESSURE_PA", 101325.);
  if (std::min(standard_t, standard_p) <= 0)
    throw std::invalid_argument("Standard gas reference T and P must be positive");
  double target_nm3 = target_moles * r * standard_t / standard_p, work_kwh = total_work / 3.6e6,
         electric_kwh = total_electric / 3.6e6;
  double purchased_heat =
             std::max(recoverable_heat - recovered_heat, 0.) / std::max(heating_efficiency, 1e-12),
         heat_kwh = purchased_heat / 3.6e6, cooling_kwh = external_cooling / 3.6e6;
  double electricity_price = number(cfg, "ELECTRICITY_COST_PER_KWH"),
         heating_price = number(cfg, "HEAT_UTILITY_COST_PER_KWH_TH"),
         cooling_price = number(cfg, "COOLING_UTILITY_COST_PER_KWH_TH");
  if (std::min({electricity_price, heating_price, cooling_price}) < 0)
    throw std::invalid_argument("Utility prices cannot be negative");
  double utility_cost = electric_kwh * electricity_price + heat_kwh * heating_price +
                        (cooling_basis == "purchased_cooling" ? cooling_kwh * cooling_price : 0.);
  double operating_hours = number(cfg, "ANNUAL_OPERATING_HOURS", 8000.),
         availability = number(cfg, "PROCESS_AVAILABILITY", 1.);
  if (operating_hours < 0 || !(availability >= 0 && availability <= 1))
    throw std::invalid_argument(
        "ANNUAL_OPERATING_HOURS must be non-negative and PROCESS_AVAILABILITY must lie in [0, 1]");
  double duration = number(metrics, "cycle_duration_s"),
         annual_cycles = duration > 0 ? operating_hours * 3600. * availability / duration : 0.;
  double bed_volume = number(cfg, "COLUMN_AREA") * number(cfg, "L"),
         adsorbent_mass = (1. - number(cfg, "eps_b")) * number(cfg, "rho_s") * bed_volume;
  auto definitions = cfg.value("COLUMNS", Json::object());
  int columns =
      static_cast<int>(number(metrics, "number_of_columns",
                              !definitions.empty() ? static_cast<double>(definitions.size())
                                                   : number(cfg, "NUMBER_OF_COLUMNS", 1)));
  double total_volume = number(metrics, "total_bed_volume_m3", columns * bed_volume),
         total_mass = number(metrics, "total_adsorbent_mass_kg", columns * adsorbent_mass),
         purity_target = number(cfg, "PRODUCT_PURITY_SPECIFICATION", nan);
  double purity = metrics.contains("light_product_purity") &&
                          !metrics.at("light_product_purity").is_null() && product_available
                      ? number(metrics, "light_product_purity")
                      : nan;
  double recovery = metrics.contains("light_component_recovery") &&
                            !metrics.at("light_component_recovery").is_null() && product_available
                        ? number(metrics, "light_component_recovery")
                        : nan;
  std::string level = cfg.value("CYCLE_LEVEL", std::string("unspecified"));
  bool feasible = flag(metrics, "process_feasible", product_available),
       dynamic_buffers = flag(metrics, "dynamic_buffers_enabled",
                              flag(cfg, "DYNAMIC_BUFFERS", flag(cfg, "BUFFERS"))),
       converged = flag(summary, "converged");
  bool recommended = converged && feasible && product_available && std::isfinite(target_moles) &&
                     target_moles > 0 && level == "optimized";
  Json record = {
      {"material", cfg.value("MATERIAL", std::string("unspecified"))},
      {"scenario", cfg.value("SCENARIO", std::string("unspecified"))},
      {"scenario_group", cfg.value("SCENARIO_GROUP", std::string("unspecified"))},
      {"cycle_level", level},
      {"recommended_for_primary_tea", recommended},
      {"process_feasible", feasible},
      {"number_of_columns", columns},
      {"product_buffer_included",
       dynamic_buffers || flag(cfg, "PRODUCT_BUFFER_INCLUDED") || recycle_moles > 0},
      {"cycle_duration_s", duration},
      {"bed_volume_each_m3", metrics.contains("total_bed_volume_m3") ? Json() : Json(bed_volume)},
      {"total_bed_volume_m3", total_volume},
      {"adsorbent_mass_each_kg",
       metrics.contains("total_adsorbent_mass_kg") ? Json() : Json(adsorbent_mass)},
      {"total_adsorbent_mass_kg", total_mass},
      {"per_column_geometry", metrics.value("per_column_geometry", Json::object())},
      {"target_product_available", product_available},
      {"target_light_component", target},
      {"product_amount_basis", amount_basis},
      {"net_target_product_mol_cycle", target_moles},
      {"net_target_product_kg_cycle", target_kg},
      {"net_target_product_standard_m3_cycle", target_nm3},
      {"standard_temperature_K", standard_t},
      {"standard_pressure_Pa", standard_p},
      {"cycle_converged", converged},
      {"cooling_cost_basis", cooling_basis},
      {"feed_pressure_bar", high_pressure / 1e5},
      {"regeneration_pressure_bar", low_pressure / 1e5},
      {"feed_supply_pressure_bar", supply_pressure / 1e5},
      {"feed_compression_charged", flag(cfg, "CHARGE_FEED_COMPRESSION")},
      {"net_methane_product_mol_cycle", target_moles},
      {"net_methane_product_kg_cycle", target_kg},
      {"net_methane_product_Nm3_cycle", target_nm3},
      {"methane_purity", purity},
      {"methane_recovery", recovery},
      {"product_purity_target", purity_target},
      {"product_purity_specification_met", std::isfinite(purity_target) && purity >= purity_target},
      {"gross_light_product_mol_cycle", number(metrics, "gross_light_product_mol")},
      {"purge_consumption_mol_cycle", number(metrics, "purge_consumption_mol")},
      {"repressurization_consumption_mol_cycle",
       number(metrics, "repressurization_consumption_mol")},
      {"product_buffer_consumption_mol_cycle", recycle_moles},
      {"net_light_product_mol_cycle", number(metrics, "light_product_mol")},
      {"product_buffer_deficit_mol_cycle", number(metrics, "product_buffer_deficit_mol")},
      {"feed_compression_work_kWh_cycle", feed_work / 3.6e6},
      {"feed_compression_if_from_ambient_kWh_cycle", ambient_feed_work / 3.6e6},
      {"vacuum_work_kWh_cycle", vacuum_work / 3.6e6},
      {"product_buffer_compression_work_kWh_cycle", product_work / 3.6e6},
      {"recycle_compression_work_kWh_cycle", recycle_work / 3.6e6},
      {"network_compression_work_kWh_cycle", network_work / 3.6e6},
      {"compression_and_vacuum_work_kWh_cycle", compression_and_vacuum_work / 3.6e6},
      {"screening_total_work_kWh_cycle", work_kwh},
      {"thermal_heat_delivered_kWh_th_cycle", thermal_heat / 3.6e6},
      {"hot_purge_and_feed_gas_heating_kWh_th_cycle", gas_heating / 3.6e6},
      {"cold_purge_and_feed_gas_cooling_kWh_th_cycle", gas_cooling / 3.6e6},
      {"recoverable_heat_kWh_th_cycle", recovered_heat / 3.6e6},
      {"ambient_heat_loss_kWh_th_cycle", ambient_heat_loss / 3.6e6},
      {"thermal_heat_purchased_kWh_th_cycle", heat_kwh},
      {"wall_cooling_duty_kWh_th_cycle", cooling_duty / 3.6e6},
      {"cooling_duty_kWh_th_cycle", cooling_kwh},
      {"gross_cooling_duty_kWh_th_cycle", total_cooling / 3.6e6},
      {"recovered_heat_kWh_th_cycle", recovered_heat / 3.6e6},
      {"cooling_electricity_kWh_cycle", cooling_electricity / 3.6e6},
      {"electric_heater_work_kWh_cycle", electric_heater / 3.6e6},
      {"total_purchased_electricity_kWh_cycle", electric_kwh},
      {"utility_cost_per_cycle", utility_cost},
      {"annual_cycles", annual_cycles},
      {"annual_purchased_electricity_kWh", electric_kwh * annual_cycles},
      {"annual_purchased_heat_kWh_th", heat_kwh * annual_cycles},
      {"annual_cooling_duty_kWh_th", cooling_kwh * annual_cycles},
      {"annual_utility_cost", utility_cost * annual_cycles},
      {"annual_net_methane_product_kg", target_kg * annual_cycles},
      {"utility_cost_per_kg_CH4", target_kg > 0 ? utility_cost / target_kg : nan},
      {"screening_total_work_if_feed_from_ambient_kWh_cycle",
       (ambient_feed_work + vacuum_work) / 3.6e6},
      {"screening_specific_energy_kWh_kg_CH4", target_kg > 0 ? work_kwh / target_kg : nan},
      {"screening_specific_energy_kWh_Nm3_CH4", target_nm3 > 0 ? work_kwh / target_nm3 : nan},
      {"purchased_electricity_kWh_kg_CH4", target_kg > 0 ? electric_kwh / target_kg : nan},
      {"purchased_heat_kWh_th_kg_CH4", target_kg > 0 ? heat_kwh / target_kg : nan},
      {"mass_balance_error", number(metrics, "mass_balance_error", nan)},
      {"energy_model", "efficiency-adjusted ideal-gas " + work_model +
                           " compression plus explicit cycle thermal duties"},
      {"dynamic_buffer_model",
       dynamic_buffers ? "finite well-mixed ideal-gas buffers included in cycle balances"
                       : "legacy cycle-integrated product-buffer accounting"},
      {"tea_status", recommended
                         ? "primary cycle-to-TEA handoff; not a complete equipment-cost model"
                         : "screening diagnostic only; primary TEA requires feasible product "
                           "export, CSS and an appropriate optimized case"},
      {"machine_work_basis",
       {{"feed", metrics.contains("feed_compression_work_J_cycle")
                     ? "measured"
                     : "pressure-level screening estimate if enabled"},
        {"vacuum", metrics.contains("vacuum_work_J_cycle") ? "measured"
                                                           : "pressure-level screening estimate"},
        {"product_delivery", metrics.contains("product_compression_work_J_cycle")
                                 ? "measured"
                                 : "pressure-level screening estimate"},
        {"recycle", metrics.contains("recycle_compression_work_J_cycle")
                        ? "measured"
                        : "screening estimate only without explicit network compression"},
        {"network", metrics.contains("network_compression_work_J_cycle")
                        ? "measured explicit device transfers"
                        : "no measured work recorded"}}}};
  record["target_product_purity"] = purity;
  if (cfg.contains("FRESH_FEED_TAGS")) record.update({
      {"fresh_feed_mol_reporting_interval", feed_moles},
      {"adsorption_feed_mol_reporting_interval", number(metrics, "adsorption_feed_mol", feed_moles)},
      {"fresh_feed_by_stage_mol", metrics.value("fresh_feed_by_stage_mol", Json::object())},
      {"average_fresh_feed_mol_s", duration > 0 ? feed_moles / duration : 0.},
      {"average_fresh_feed_Nm3_h", duration > 0 ? feed_moles / duration * r * standard_t / standard_p * 3600. : 0.},
      {"plant_fresh_feed_target_Nm3_h", number(cfg, "PLANT_FRESH_FEED_TARGET_NM3_H")},
      {"plant_capacity_basis", "cycle-average tagged fresh feed including adsorption, repressurization, cooling and waiting time; replication estimates capacity, not continuity"},
      {"product_delivery_pressure_Pa", delivery_pressure},
      {"product_buffer_pressure_Pa", buffer_pressure}});
  if (summary.contains("convergence_report")) {
    record["convergence_report"] = summary["convergence_report"];
    record["metrics_basis"] = summary.value("reporting_metrics_basis", std::string("last cycle"));
    record["reporting_interval_cycles"] = summary.value("converged", false) ? Json(1) : summary.value("window_aggregate_cycles", Json(1));
    record["interval_note"] = "For window reporting, *_cycle quantities refer to the full reporting interval, not an individual physical cycle. Annualization uses the same interval duration.";
  }
  record["target_product_recovery"] = recovery;
  record["annual_net_target_product_kg"] = target_kg * annual_cycles;
  record["utility_cost_per_kg_target"] = target_kg > 0 ? Json(utility_cost / target_kg) : Json();
  record["purchased_electricity_kWh_kg_target"] =
      target_kg > 0 ? Json(electric_kwh / target_kg) : Json();
  record["purchased_heat_kWh_th_kg_target"] = target_kg > 0 ? Json(heat_kwh / target_kg) : Json();
  record["column_geometry_basis"] = metrics.contains("total_adsorbent_mass_kg")
                                        ? "measured cycle model inventories"
                                        : "legacy uniform-bed approximation";
  if (upper(target) != "CH4") {
    std::vector<std::string> remove;
    for (auto it = record.begin(); it != record.end(); ++it)
      if (it.key().find("methane") != std::string::npos ||
          it.key().find("CH4") != std::string::npos)
        remove.push_back(it.key());
    for (const auto& key : remove) record.erase(key);
  }
  return clean(std::move(record));
}
Json write_tea_record(const Json& cfg, const Json& summary,
                      const std::filesystem::path& output_path) {
  const auto record = build_tea_record(cfg, summary);
  auto path = output_path.empty() ? std::filesystem::path(cfg.at("OUTPUT_DIR").get<std::string>()) /
                                        "tea_handoff.json"
                                  : output_path;
  if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Cannot write TEA handoff: " + path.string());
  out << std::setw(2) << record << '\n';
  if (!out) throw std::runtime_error("Failed writing TEA handoff: " + path.string());
  auto csv_path = path;
  csv_path.replace_extension(".csv");
  std::ofstream csv(csv_path);
  if (!csv) throw std::runtime_error("Cannot write TEA CSV: " + csv_path.string());
  bool first = true;
  for (auto it = record.begin(); it != record.end(); ++it) {
    if (!first) csv << ',';
    first = false;
    csv << csv_text(Json(it.key()));
  }
  csv << "\r\n";
  first = true;
  for (const auto& v : record) {
    if (!first) csv << ',';
    first = false;
    csv << csv_text(v);
  }
  csv << "\r\n";
  if (!csv) throw std::runtime_error("Failed writing TEA CSV: " + csv_path.string());
  return record;
}
}  // namespace cycle

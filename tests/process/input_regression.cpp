#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "cycle/input.hpp"
using cycle::Json;
void check(bool value, const char* message) {
  if (!value) throw std::runtime_error(message);
}
template <class F>
void rejects(F f) {
  bool bad = false;
  try {
    f();
  } catch (const std::exception&) {
    bad = true;
  }
  check(bad, "Invalid input was accepted");
}
int main() {
  try {
    Json d = Json::parse(R"JSON({
      "schema_version":1,
      "components":[{"name":"B","molar_mass_kg_mol":0.028},{"name":"A","molar_mass_kg_mol":0.044}],
      "equilibrium":{"reference_temperature_K":300.0,"components":{
        "A":{"sites":[{"capacity_mol_kg":2.0,"affinity_Pa_inv":0.00001,"heat_J_mol":10000.0}]},
        "B":{"sites":[{"capacity_mol_kg":1.0,"affinity_Pa_inv":0.00002,"heat_J_mol":20000.0}]}}},
      "initial":{"composition":{"A":0.1,"B":0.9}},
      "columns":{"Z":{},"A":{"initial":{"composition":{"A":0.2,"B":0.8}}}},
      "cycle":{"phases":[{"name":"idle","duration_s":1.0}]}
    })JSON");
    auto c = cycle::normalize_config(d);
    check(c["comp_names"][0] == "B", "Component ordering changed");
    check(c["COLUMNS"].begin().key() == "Z", "Column ordering changed");
    check(c["y_init"][0] == .9 && c["COLUMNS"]["A"]["y_init"][0] == .8,
          "Named composition mapping failed");
    check(c["CYCLE_PHASES"][0]["duration"] == 1., "Duration alias missing");
    double k0 = c["PHY_ISOTHERM_INPUT"][0]["components"]["A"]["sites"][0]["k0"];
    check(std::abs(k0 * std::exp(10000. / (8.31446261815324 * 300.)) / 1e-5 - 1.) < 1e-14,
          "Affinity conversion failed");
    auto bad = d;
    bad["bed"] = {{"lenght_m", 1.}};
    rejects([&] { cycle::normalize_config(bad); });
    bad = d;
    bad["initial"]["composition"]["A"] = .4;
    rejects([&] { cycle::normalize_config(bad); });
    bad = d;
    bad["cycle"]["steps"] = Json::array();
    rejects([&] { cycle::normalize_config(bad); });
    bad = d;
    bad["columns"]["A"]["initial"]["pressure_Pa"] = 100.;
    bad["columns"]["A"]["initial"]["pressure"] = 200.;
    rejects([&] { cycle::normalize_config(bad); });
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    auto path = std::filesystem::temp_directory_path() /
                ("ruptura_yaml_" + std::to_string(stamp) + ".yaml");
    {
      std::ofstream f(path);
      f << "description: 'quoted: text'\nnumbers: [1, 2.0, 1.0e-06]\nshared: &x [3, 4]\ncopy: *x\n";
    }
    auto y = cycle::read_input_document(path);
    check(y["description"] == "quoted: text" && y["numbers"][2] == 1e-6 && y["shared"] == y["copy"],
          "Native YAML parser failed");
    {
      std::ofstream f(path);
      f << "key: 1\nkey: 2\n";
    }
    rejects([&] { cycle::read_input_document(path); });
    {
      std::ofstream f(path);
      f << "key: 1\n---\nkey: 2\n";
    }
    rejects([&] { cycle::read_input_document(path); });
    std::filesystem::remove(path);
    std::cout << "Native input normalization regression passed\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

#include <iostream>

#include "cycle/input.hpp"
#include "cycle/network.hpp"
int main() {
  try {
    cycle::Json request;
    std::cin >> request;
    cycle::ColumnNetwork n(cycle::load_config(request.at("input").get<std::string>()));
    auto initial = n.initial_state();
    cycle::Json out = {{"initial", initial},
                       {"inventories", n.inventories(initial)},
                       {"atol", n.state_atol()},
                       {"columns", n.column_names},
                       {"checks", cycle::Json::array()}};
    for (auto& c : request.at("commands")) {
      auto s = c.at("state").get<cycle::Vec>();
      n.set_phase(c.at("phase"), c.value("context", cycle::Json::object()));
      auto [d, diag] = n.evaluate(c.value("time", 0.), s);
      out["checks"].push_back(
          {{"rhs", d}, {"diagnostics", diag}, {"inventories", n.inventories(s)}});
    }
    std::cout << out.dump() << '\n';
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

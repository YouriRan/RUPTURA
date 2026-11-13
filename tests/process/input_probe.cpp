#include <iostream>

#include "cycle/input.hpp"
int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "Usage: input_probe simulation.yaml\n";
      return 2;
    }
    std::cout << cycle::load_config(argv[1]).dump(2) << '\n';
    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}

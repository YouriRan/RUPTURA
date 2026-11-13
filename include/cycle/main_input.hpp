#pragma once
#include "cycle/input.hpp"
namespace cycle
{
// Public input follows RUPTURA main conventions. Legacy package mappings are private.
Json canonicalize_chemisorption_input(const Json&);
bool is_process_input(const Json&);
Json normalize_main_input(const Json&, const std::filesystem::path& base_dir = ".");
Json load_main_input(const std::filesystem::path&);
Json process_capabilities();
void write_main_outputs(const Json& config, const std::filesystem::path& destination);
}  // namespace cycle

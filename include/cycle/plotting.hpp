#pragma once
#include <filesystem>

namespace cycle
{
// Read legacy column/component text outputs and write headless PNG/PDF figures.
void write_legacy_breakthrough_plots(const std::filesystem::path& input_path);
}

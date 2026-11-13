#pragma once
#include <filesystem>

#include "cycle/types.hpp"
namespace cycle
{
Json read_input_document(const std::filesystem::path& path);
Json normalize_config(const Json& document, const std::filesystem::path& base_dir = ".");
Json load_config(const std::filesystem::path& path);
}  // namespace cycle

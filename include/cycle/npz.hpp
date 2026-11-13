#pragma once
#include <filesystem>
#include <map>

#include "cycle/types.hpp"
namespace cycle
{
struct NpyArray
{
  std::vector<std::size_t> shape;
  std::vector<double> data;
};
// Native ZIP/NPY writer with optional lossless DEFLATE; little-endian C-order arrays. JSON
// supports rectangular numeric/bool/UTF-8-text arrays and scalar values.
// The explicit overload supports empty shapes and IEEE NaN/Inf doubles.
// ZIP32 only: archives exceeding 4 GiB are rejected explicitly.
void write_npz(const std::filesystem::path&, const std::map<std::string, NpyArray>&, bool compress = false);
void write_npz(const std::filesystem::path&, const Json&, bool compress = false);
// Read a single native numeric array; validates NPY layout, dimensions and ZIP CRC.
NpyArray read_npz_array(const std::filesystem::path&, const std::string& key);
}  // namespace cycle

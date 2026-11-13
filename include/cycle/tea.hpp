#pragma once
#include <filesystem>
#include <string>

#include "cycle/types.hpp"
namespace cycle
{
double compression_work(double moles, double temperature, double inlet_pressure, double outlet_pressure,
                        double efficiency, double gas_constant, const std::string& model = "isothermal",
                        double heat_capacity_ratio = 1.4);
// Config is the normalized legacy mapping. Summary contains final_metrics or metrics.
// Unavailable/nonfinite diagnostics are represented by JSON null.
Json build_tea_record(const Json& cfg, const Json& summary);
Json write_tea_record(const Json& cfg, const Json& summary, const std::filesystem::path& output_path = {});
}  // namespace cycle

#pragma once
#include <string>
#include <vector>

#include "cycle/types.hpp"
namespace cycle
{
Vec component_array(const Json& value, std::size_t n, const std::string& name);
Json channel_geometry(const Json& options, double total_area);
// Material-by-node arrays use material-major indexing material*nz+node.
struct ColumnGeometry
{
  Vec z, widths, edges;
  double area = 0, L = 0;
  std::vector<std::string> names, material_names;
  std::vector<Json> materials, monoliths;
  Vec interfaces, fractions;
  Vec epsilon_material, storage_epsilon_material, mass_material, heat_material;
  Vec a_material, linear_resistance, quadratic_resistance;
  Vec epsilon, storage_epsilon, face_epsilon, epsilon_face, mass_density;
  Vec heat_capacity, hgs_as, resistance_linear, resistance_quadratic, material_mass_kg;
  ColumnGeometry(const Json& cfg, const Vec& coordinates, const Vec& cell_widths, double column_area);
  ColumnGeometry(const Json& cfg, std::size_t nz, double length, double column_area);
};
using Geometry = ColumnGeometry;
}  // namespace cycle

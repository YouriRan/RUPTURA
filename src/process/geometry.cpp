#include "cycle/geometry.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "cycle/spatial.hpp"
namespace cycle
{
namespace
{
double number(const Json& j, const char* key, double fallback)
{
  return j.contains(key) ? j.at(key).get<double>() : fallback;
}
double primitive(double x, double b, double w)
{
  if (w == 0) return std::max(x - b, 0.0);
  double lo = b - w / 2, hi = b + w / 2;
  return x <= lo ? 0 : x < hi ? (x - lo) * (x - lo) / (2 * w) : x - b;
}
Vec coordinates(std::size_t n, double L)
{
  if (n < 2 || !std::isfinite(L) || L <= 0)
    throw std::invalid_argument("Column requires two or more nodes and positive length");
  Vec z(n);
  for (std::size_t i = 0; i < n; ++i) z[i] = L * i / (n - 1);
  return z;
}
Vec widths_for(std::size_t n, double L)
{
  auto e = control_volume_faces(coordinates(n, L));
  Vec w(n);
  for (std::size_t i = 0; i < n; ++i) w[i] = e[i + 1] - e[i];
  return w;
}
}  // namespace
Vec component_array(const Json& value, std::size_t n, const std::string& name)
{
  Vec result;
  const Json* a = &value;
  while (a->is_array() && a->size() == 1) a = &a->at(0);
  if (a->is_number())
    result.assign(n, a->get<double>());
  else if (a->is_array())
    for (const auto& item : *a)
    {
      const Json* p = &item;
      while (p->is_array() && p->size() == 1) p = &p->at(0);
      if (!p->is_number()) throw std::invalid_argument(name + " must have one finite value per component");
      result.push_back(p->get<double>());
    }
  if (result.size() != n || std::any_of(result.begin(), result.end(), [](double v) { return !std::isfinite(v); }))
    throw std::invalid_argument(name + " must have one finite value per component");
  return result;
}
Json channel_geometry(const Json& options, double total_area)
{
  constexpr double pi = 3.141592653589793238462643383279502884;
  std::string shape = options.value("shape", std::string("square"));
  std::transform(shape.begin(), shape.end(), shape.begin(), [](unsigned char c) { return std::tolower(c); });
  if (shape == "circle") shape = "circular";
  if (shape == "triangle") shape = "triangular";
  if (shape == "hexagon") shape = "hexagonal";
  double size = number(options, "channel_size_m", 1e-3), coat = number(options, "coating_thickness_m", 0);
  if (!(size > 0) || coat < 0 || !std::isfinite(size) || !std::isfinite(coat) || !(total_area > 0))
    throw std::invalid_argument("Monolith channel size and total area must be positive; coating thickness nonnegative");
  double a, p, r, alpha;
  if (shape == "circular")
  {
    a = pi * size * size / 4;
    p = pi * size;
    r = size / 2;
    alpha = 32;
  }
  else if (shape == "square")
  {
    a = size * size;
    p = 4 * size;
    r = size / 2;
    alpha = 28.45415377;
  }
  else if (shape == "triangular")
  {
    a = std::sqrt(3.0) * size * size / 4;
    p = 3 * size;
    r = size / (2 * std::sqrt(3.0));
    alpha = 80.0 / 3;
  }
  else if (shape == "hexagonal")
  {
    a = 3 * std::sqrt(3.0) * size * size / 2;
    p = 6 * size;
    r = std::sqrt(3.0) * size / 2;
    alpha = 30.108;
  }
  else
    throw std::invalid_argument("Monolith shape must be circular, square, triangular, or hexagonal");
  double count, eps;
  if (options.contains("channel_count"))
  {
    count = options.at("channel_count").get<double>();
    if (!std::isfinite(count) || count < 1 || count != std::floor(count))
      throw std::invalid_argument("channel_count must be a positive integer");
    eps = count * a / total_area;
    if (options.contains("open_fraction") && std::abs(eps - options.at("open_fraction").get<double>()) >
                                                 1e-8 + 1e-5 * std::abs(options.at("open_fraction").get<double>()))
      throw std::invalid_argument("channel_count and open_fraction describe different channel areas");
  }
  else
  {
    eps = number(options, "open_fraction", .65);
    count = eps * total_area / a;
  }
  double coating = eps * (std::pow(1 + coat / r, 2) - 1), substrate = 1 - eps - coating,
         solid = number(options, "coating_solid_fraction", 1);
  if (!(eps > 0 && eps < 1) || substrate < -1e-12)
    throw std::invalid_argument("Monolith free channel and coating areas exceed total column area");
  if (!(solid > 0 && solid <= 1)) throw std::invalid_argument("coating_solid_fraction must lie in (0,1]");
  return Json{{"shape", shape},
              {"epsilon", eps},
              {"coating_fraction", coating},
              {"active_fraction", coating * solid},
              {"substrate_fraction", std::max(substrate, 0.0)},
              {"substrate_density", number(options, "substrate_density_kg_m3", 2500)},
              {"substrate_cp", number(options, "substrate_heat_capacity_J_kg_K", 800)},
              {"interfacial_area", count * p / total_area},
              {"hydraulic_diameter", 4 * a / p},
              {"alpha", alpha},
              {"channel_count", count},
              {"channel_area", a},
              {"channel_perimeter", p}};
}
ColumnGeometry::ColumnGeometry(const Json& cfg, std::size_t nz, double length, double column_area)
    : ColumnGeometry(cfg, coordinates(nz, length), widths_for(nz, length), column_area)
{
}
ColumnGeometry::ColumnGeometry(const Json& cfg, const Vec& coordinates_, const Vec& cell_widths, double column_area)
    : z(coordinates_), widths(cell_widths), area(column_area), L(cfg.at("L").get<double>())
{
  const auto nz = z.size();
  if (nz < 2 || widths.size() != nz || !(area > 0) || !(L > 0)) throw std::invalid_argument("Invalid column geometry");
  edges = control_volume_faces(z);
  edges.front() = 0;
  edges.back() = L;
  for (double w : widths)
    if (!(w > 0) || !std::isfinite(w)) throw std::invalid_argument("Cell widths must be positive and finite");
  Json defs = cfg.value("MATERIALS", Json{{"default", Json::object()}});
  if (defs.is_null()) defs = Json{{"default", Json::object()}};
  if (!defs.is_object() || defs.empty()) throw std::invalid_argument("MATERIALS must be a nonempty mapping");
  if (cfg.contains("MATERIAL_ORDER"))
    names = cfg.at("MATERIAL_ORDER").get<std::vector<std::string>>();
  else
    for (auto it = defs.begin(); it != defs.end(); ++it) names.push_back(it.key());
  if (names.size() != defs.size()) throw std::invalid_argument("MATERIAL_ORDER must list every material exactly once");
  auto unique = names;
  std::sort(unique.begin(), unique.end());
  if (std::adjacent_find(unique.begin(), unique.end()) != unique.end())
    throw std::invalid_argument("MATERIAL_ORDER contains duplicates");
  material_names = names;
  for (const auto& name : names)
  {
    Json merged = cfg;
    for (auto it = defs.at(name).begin(); it != defs.at(name).end(); ++it) merged[it.key()] = it.value();
    materials.push_back(std::move(merged));
  }
  Json layers = cfg.value("LAYERS", Json());
  if (layers.is_null() || layers.empty())
  {
    if (names.size() != 1 && !cfg.contains("MIX_FRACTIONS"))
      throw std::invalid_argument("Multiple materials require ColumnSections or MixFraction");
    layers = Json::array();
    for (const auto& name : names) layers.push_back({{"material", name}, {"length_m", L / names.size()}});
  }
  Vec ramps;
  double sum = 0;
  for (std::size_t i = 0; i < layers.size(); ++i)
  {
    double length = layers[i].at("length_m").get<double>();
    if (!(length > 0) || !std::isfinite(length)) throw std::invalid_argument("Positive LAYERS lengths must sum to L");
    sum += length;
    if (i + 1 < layers.size())
    {
      interfaces.push_back(sum);
      double w = number(layers[i], "interface_after_m", 0);
      if (w < 0 || !std::isfinite(w)) throw std::invalid_argument("interface_after_m cannot be negative");
      ramps.push_back(w);
    }
  }
  if (std::abs(sum - L) > 1e-12 + 1e-10 * std::abs(L))
    throw std::invalid_argument("Positive LAYERS lengths must sum to L");
  for (std::size_t i = 0; i < ramps.size(); ++i)
    if (interfaces[i] - ramps[i] / 2 < 0 || interfaces[i] + ramps[i] / 2 > L ||
        (i + 1 < ramps.size() && interfaces[i] + ramps[i] / 2 > interfaces[i + 1] - ramps[i + 1] / 2))
      throw std::invalid_argument("Finite transition regions must lie inside bed and cannot overlap");
  Vec cumulative(interfaces.size() * nz);
  for (std::size_t b = 0; b < interfaces.size(); ++b)
    for (std::size_t j = 0; j < nz; ++j)
      cumulative[b * nz + j] =
          (primitive(edges[j + 1], interfaces[b], ramps[b]) - primitive(edges[j], interfaces[b], ramps[b])) / widths[j];
  fractions.assign(names.size() * nz, 0);
  for (std::size_t i = 0; i < layers.size(); ++i)
  {
    std::string name = layers[i].at("material").get<std::string>();
    auto it = std::find(names.begin(), names.end(), name);
    if (it == names.end()) throw std::invalid_argument("Unknown layer material " + name);
    std::size_t m = it - names.begin();
    for (std::size_t j = 0; j < nz; ++j)
      fractions[m * nz + j] +=
          (i == 0 ? 1 : cumulative[(i - 1) * nz + j]) - (i + 1 == layers.size() ? 0 : cumulative[i * nz + j]);
  }
  if (cfg.contains("MIX_FRACTIONS"))
  {
    Vec mix = cfg.at("MIX_FRACTIONS").get<Vec>();
    if (mix.size() != names.size()) throw std::invalid_argument("MixFraction requires one value per adsorbent");
    double total = 0.;
    for (double v : mix)
    {
      if (!std::isfinite(v) || v < 0.) throw std::invalid_argument("MixFraction must be finite and nonnegative");
      total += v;
    }
    if (std::abs(total - 1.) > 1e-10) throw std::invalid_argument("MixFraction values must sum to one");
    for (std::size_t m = 0; m < mix.size(); ++m)
      for (std::size_t i = 0; i < nz; ++i) fractions[m * nz + i] = mix[m];
    interfaces.clear();
  }
  for (double& f : fractions) f = std::max(f, 0.0);
  for (std::size_t j = 0; j < nz; ++j)
  {
    double total = 0;
    for (std::size_t m = 0; m < names.size(); ++m) total += fractions[m * nz + j];
    if (std::abs(total - 1) > 1e-5 + 1e-12)
      throw std::invalid_argument("Layer material fractions do not partition column");
  }
  for (const auto& mat : materials)
  {
    std::string kind = mat.value("CONTACTOR", std::string("fixed_bed"));
    std::transform(kind.begin(), kind.end(), kind.begin(),
                   [](unsigned char c) { return c == '-' ? '_' : std::tolower(c); });
    double rho = mat.at("rho_s").get<double>(), cp = mat.at("Cps").get<double>();
    if (!(rho > 0 && cp > 0)) throw std::invalid_argument("Material rho_s and Cps must be positive");
    double eps, mass, heat, a, linear, quadratic;
    Json mono;
    if (kind == "fixedbed" || kind == "fixed_bed" || kind == "packed_bed")
    {
      eps = mat.at("eps_b").get<double>();
      double dp = mat.at("dp").get<double>();
      double phi = number(mat, "PARTICLE_SPHERICITY", 1.0);
      if (!std::isfinite(phi) || phi <= 0.0 || phi > 1.0)
        throw std::invalid_argument("ParticleSphericity must lie in (0, 1]");
      const double ergun_diameter = phi * dp;
      if (!(eps > 0 && eps < 1 && dp > 0))
        throw std::invalid_argument("Fixed-bed eps_b must lie in (0,1), and dp must be positive");
      mass = (1 - eps) * rho;
      heat = mass * cp;
      a = 6 * (1 - eps) / dp;
      linear = 150 * std::pow(1 - eps, 2) / (std::pow(eps, 3) * ergun_diameter * ergun_diameter);
      quadratic = 1.75 * (1 - eps) / (std::pow(eps, 3) * ergun_diameter);
    }
    else if (kind == "monolith")
    {
      mono = channel_geometry(mat.value("MONOLITH", Json::object()), area);
      eps = mono.at("epsilon");
      mass = mono.at("active_fraction").get<double>() * rho;
      if (!(mass > 0)) throw std::invalid_argument("An adsorbing monolith requires positive coating_thickness_m");
      heat = mass * cp + mono.at("substrate_fraction").get<double>() * mono.at("substrate_density").get<double>() *
                             mono.at("substrate_cp").get<double>();
      a = mono.at("interfacial_area");
      linear = mono.at("alpha").get<double>() / (std::pow(mono.at("hydraulic_diameter").get<double>(), 2) * eps);
      quadratic = 0;
    }
    else
      throw std::invalid_argument("CONTACTOR must be fixed_bed or monolith");
    double storage = number(mat, "GAS_STORAGE_VOID_FRACTION", eps);
    if (!std::isfinite(storage) || storage < eps || storage >= 1)
      throw std::invalid_argument("GAS_STORAGE_VOID_FRACTION must be finite, >= flow void fraction and < 1");
    storage_epsilon_material.push_back(storage);
    epsilon_material.push_back(eps);
    mass_material.push_back(mass);
    heat_material.push_back(heat);
    a_material.push_back(a);
    linear_resistance.push_back(linear);
    quadratic_resistance.push_back(quadratic);
    monoliths.push_back(mono);
  }
  epsilon.assign(nz, 0);
  storage_epsilon.assign(nz, 0);
  heat_capacity.assign(nz, 0);
  hgs_as.assign(nz, 0);
  resistance_linear.assign(nz, 0);
  resistance_quadratic.assign(nz, 0);
  mass_density.resize(names.size() * nz);
  material_mass_kg.assign(names.size(), 0);
  for (std::size_t m = 0; m < names.size(); ++m)
    for (std::size_t j = 0; j < nz; ++j)
    {
      double f = fractions[m * nz + j];
      epsilon[j] += epsilon_material[m] * f;
      storage_epsilon[j] += storage_epsilon_material[m] * f;
      mass_density[m * nz + j] = mass_material[m] * f;
      heat_capacity[j] += heat_material[m] * f;
      hgs_as[j] += a_material[m] * materials[m].at("hgs").get<double>() * f;
      resistance_linear[j] += linear_resistance[m] * f;
      resistance_quadratic[j] += quadratic_resistance[m] * f;
      material_mass_kg[m] += area * mass_density[m * nz + j] * widths[j];
    }
  face_epsilon.resize(nz + 1);
  face_epsilon.front() = epsilon.front();
  face_epsilon.back() = epsilon.back();
  for (std::size_t i = 1; i < nz; ++i) face_epsilon[i] = .5 * (epsilon[i - 1] + epsilon[i]);
  epsilon_face = face_epsilon;
}
}  // namespace cycle

#include "cycle/column.hpp"
#include "cycle/main_models.hpp"
#include "compute.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>

#include "cycle/equilibrium.hpp"
#include "cycle/geometry.hpp"
#include "cycle/kinetics.hpp"
#include "cycle/process.hpp"
#include "cycle/spatial.hpp"
#include "ruptura2/particle_transport.hpp"
namespace cycle
{
namespace
{
constexpr double pi = 3.1415926535897932384626433832795;
std::string lower(std::string s)
{
  for (char& c : s)
  {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (c == '-') c = '_';
  }
  return s;
}
double sum(const Vec& a) { return std::accumulate(a.begin(), a.end(), 0.0); }
double dot(const Vec& a, const Vec& b)
{
  double r = 0;
  for (std::size_t i = 0; i < a.size(); ++i) r += a[i] * b[i];
  return r;
}
Vec slice(const Vec& a, std::size_t start, std::size_t size)
{
  return Vec(a.begin() + start, a.begin() + start + size);
}
void insert(Vec& out, const Vec& v, std::size_t pos) { std::copy(v.begin(), v.end(), out.begin() + pos); }
Vec arr(const Json& value, std::size_t n, const std::string& name) { return component_array(value, n, name); }
Vec normal(const Json& value, std::size_t n, const std::string& name)
{
  Vec a = arr(value, n, name);
  double total = sum(a);
  if (total <= 0) throw std::invalid_argument(name + " must have positive sum");
  for (double& v : a)
  {
    if (v < 0) throw std::invalid_argument(name + " cannot be negative");
    v /= total;
  }
  return a;
}
Vec profile(const Json& value, std::size_t n, const std::string& name) { return arr(value, n, name); }
double safe(double v, double lo, double hi, double nan)
{
  if (std::isnan(v)) return nan;
  return std::clamp(v, lo, hi);
}
Vec safe_array(Vec a, double lo, double hi, double nan)
{
  for (double& v : a) v = safe(v, lo, hi, nan);
  return a;
}
Vec components_at(const Vec& a, std::size_t nc, std::size_t nz, std::size_t i)
{
  Vec v(nc);
  for (std::size_t j = 0; j < nc; ++j) v[j] = a[j * nz + i];
  return v;
}
Vec composition_profile(const Json& value, std::size_t nc, std::size_t nz, const std::string& name, bool normalize)
{
  if (!value.is_array() || value.empty())
    throw std::invalid_argument(name + " requires a component vector or component-by-node array");
  Vec out(nc * nz);
  if (value[0].is_array())
  {
    bool transpose = value.size() == nz && value[0].size() == nc;
    if (!transpose && (value.size() != nc || value[0].size() != nz))
      throw std::invalid_argument(name + " invalid shape");
    for (std::size_t j = 0; j < nc; ++j)
      for (std::size_t i = 0; i < nz; ++i)
        out[j * nz + i] = transpose ? value[i][j].get<double>() : value[j][i].get<double>();
  }
  else
  {
    Vec a = arr(value, nc, name);
    for (std::size_t j = 0; j < nc; ++j) std::fill(out.begin() + j * nz, out.begin() + (j + 1) * nz, a[j]);
  }
  for (std::size_t i = 0; i < nz; ++i)
  {
    double total = 0;
    for (std::size_t j = 0; j < nc; ++j)
    {
      double v = out[j * nz + i];
      if (!std::isfinite(v) || v < 0) throw std::invalid_argument(name + " must be finite nonnegative");
      total += v;
    }
    if (normalize)
    {
      if (total <= 0) throw std::invalid_argument(name + " must have positive sum");
      for (std::size_t j = 0; j < nc; ++j) out[j * nz + i] /= total;
    }
  }
  return out;
}
Json matrix_json(const Vec& a, std::size_t rows, std::size_t cols)
{
  Json result = Json::array();
  for (std::size_t r = 0; r < rows; ++r) result.push_back(slice(a, r * cols, cols));
  return result;
}
double hydraulic_root(double gradient, double a, double b)
{
  double m = std::abs(gradient);
  return std::copysign(2 * m / std::max(a + std::hypot(a, 2 * std::sqrt(std::max(b, 0.0) * m)), 1e-300), gradient);
}
Vec conductive_divergence(const Vec& T, const Vec& k, const Vec& coordinates, const Vec& widths)
{
  std::size_t n = T.size();
  Vec flux(n + 1), out(n);
  for (std::size_t i = 1; i < n; ++i)
  {
    double kf = k[i - 1] + k[i] > 0 ? 2 * k[i - 1] * k[i] / (k[i - 1] + k[i]) : 0;
    flux[i] = kf * (T[i] - T[i - 1]) / (coordinates[i] - coordinates[i - 1]);
  }
  for (std::size_t i = 0; i < n; ++i) out[i] = (flux[i + 1] - flux[i]) / widths[i];
  return out;
}
}  // namespace
double resolve_time_profile(const Json& value, double t)
{
  if (value.is_number()) return value.get<double>();
  if (!value.is_object()) throw std::invalid_argument("Time-dependent value must be numeric or a profile mapping");
  std::string model = value.value("model", std::string{});
  double elapsed = t + value.value("time_offset_s", 0.0);
  if (model == "linear" || model == "exponential")
  {
    std::set<std::string> allowed = {"model", "initial_value", "final_value", "time_offset_s",
                                     model == "linear" ? "ramp_time_s" : "decay_rate_s_inv"};
    for (auto it = value.begin(); it != value.end(); ++it)
      if (!allowed.count(it.key())) throw std::invalid_argument("Unknown time-profile key " + it.key());
    double a = value.at("initial_value"), b = value.at("final_value"),
           k = value.at(model == "linear" ? "ramp_time_s" : "decay_rate_s_inv");
    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(k) || !std::isfinite(elapsed) || k <= 0)
      throw std::invalid_argument("Invalid time-profile parameters");
    return model == "linear" ? a + (b - a) * std::clamp(elapsed / k, 0.0, 1.0)
                             : b + (a - b) * std::exp(-k * std::max(elapsed, 0.0));
  }
  std::set<std::string> allowed = {"times_s", "values", "time_offset_s"};
  for (auto it = value.begin(); it != value.end(); ++it)
    if (!allowed.count(it.key())) throw std::invalid_argument("Unknown tabulated-profile key " + it.key());
  Vec times = value.at("times_s").get<Vec>(), values = value.at("values").get<Vec>();
  if (times.size() < 2 || times.size() != values.size() || times.front() < 0 || !std::isfinite(elapsed))
    throw std::invalid_argument("Invalid tabulated time profile");
  for (std::size_t i = 0; i < times.size(); ++i)
    if (!std::isfinite(times[i]) || !std::isfinite(values[i]) || (i && times[i] <= times[i - 1]))
      throw std::invalid_argument("Invalid tabulated time profile");
  if (elapsed <= times.front()) return values.front();
  if (elapsed >= times.back()) return values.back();
  std::size_t i = std::upper_bound(times.begin(), times.end(), elapsed) - times.begin();
  double f = (elapsed - times[i - 1]) / (times[i] - times[i - 1]);
  return values[i - 1] + f * (values[i] - values[i - 1]);
}
struct Column::Impl
{
  ColumnGeometry geometry;
  std::vector<MaterialEquilibrium> equilibria, particle_equilibria;
  std::vector<std::optional<ChemistryKinetics>> chemistries;
  std::vector<std::optional<MainMaterialModels>> main_models;
  mutable Vec reaction_rates;
  mutable double reaction_heat=0.;
  std::vector<std::vector<std::size_t>> loading_groups;
  std::vector<std::vector<std::vector<std::size_t>>> chemical_groups;
  std::vector<std::size_t> history_start, history_count;
  Vec phy_factors, chem_factors, cp, ads_cp, ads_ref_u, ads_ref_T;
  Json step = Json::object(), context = Json::object();
  std::vector<Json> particles;
  std::vector<std::optional<ruptura2::RadialParticleMesh>> meshes;
  std::vector<std::size_t> radial_start;
  bool step_set = false;
  double mu, Cps, lambda_g, lambda_w, rho_w, Cpw, hgw, hw_inf, T_inf, aw_int, aw_ext, wall_area, wall_volume;
  double T_floor, T_ceiling, P_ceiling, C_ceiling, V_max, boundary_length, Q_ceiling;
  double jacket_area = 0, jacket_rho = 0, jacket_cp = 0;
  Json htf = Json::object();
  bool darcy_only = false, fixed_wall_ends = false;
  Vec wall_end_temperature, cached_liquid_pressure;
  std::string advection_scheme, spatial_method;
  Impl(const Json& cfg, const Vec& z, const Vec& widths, double area) : geometry(cfg, z, widths, area) {}
};
Column::Column(const Json& config) : cfg(config)
{
  R = cfg.value("R", R);
  liquid = lower(cfg.value("FLUID_PHASE", std::string("gas"))) == "liquid";
  liquid_density = cfg.value("LIQUID_DENSITY", 1000.);
  liquid_heat_capacity = cfg.value("LIQUID_HEAT_CAPACITY", 4180.);
  if (!(liquid_density > 0) || !(liquid_heat_capacity > 0))
    throw std::invalid_argument("Liquid density and heat capacity must be positive");
  comp_names = cfg.at("comp_names").get<std::vector<std::string>>();
  n_comp = comp_names.size();
  L = cfg.at("L");
  Nz = cfg.value("Nz", cfg.value("NGRID", 60));
  dint = cfg.at("dint");
  dout = cfg.value("dout", dint);
  area = cfg.value("COLUMN_AREA", pi * dint * dint / 4);
  if (n_comp < 1 || L <= 0 || Nz < 3 || area <= 0 || dint <= 0)
    throw std::invalid_argument("Column requires components, positive geometry and Nz >= 3");
  dz = L / (Nz - 1);
  z.resize(Nz);
  widths.assign(Nz, dz);
  widths.front() = widths.back() = .5 * dz;
  for (std::size_t i = 0; i < Nz; ++i) z[i] = dz * i;
  if (cfg.contains("COLUMN_DISTANCES"))
  {
    z = cfg.at("COLUMN_DISTANCES").get<Vec>();
    if (z.size() != Nz || std::abs(z.front()) > 1e-12 || std::abs(z.back() - L) > 1e-12)
      throw std::invalid_argument("ColumnDistances must match the interval count and span 0 to ColumnLength");
    for (std::size_t i = 1; i < Nz; ++i)
      if (!std::isfinite(z[i]) || z[i] <= z[i - 1])
        throw std::invalid_argument("ColumnDistances must be finite and strictly increasing");
    auto faces = control_volume_faces(z);
    for (std::size_t i = 0; i < Nz; ++i) widths[i] = faces[i + 1] - faces[i];
  }
  M = arr(cfg.at("M"), n_comp, "M");
  for (double v : M)
    if (v <= 0) throw std::invalid_argument("M must be positive");
  p = std::make_unique<Impl>(cfg, z, widths, area);
  auto& geom = p->geometry;
  material_names = geom.names;
  epsilon_node = geom.epsilon;
  storage_epsilon_node = geom.storage_epsilon;
  epsilon_face = geom.face_epsilon;
  thermal_mode = cfg.value("THERMAL_MODE", 1);
  gas_solid_thermal_equilibrium = cfg.value("GAS_SOLID_THERMAL_EQUILIBRIUM", false);
  if (thermal_mode < 1 || thermal_mode > 3) throw std::invalid_argument("THERMAL_MODE must be 1, 2 or 3");
  p->mu = cfg.at("mu");
  p->Cps = cfg.at("Cps");
  Cpg = cfg.value("Cpg", 1000.0);
  p->lambda_g = cfg.value("lambda_g", 0.0);
  p->lambda_w = cfg.value("lambda_w", 0.0);
  p->rho_w = cfg.value("rho_w", 1.0);
  p->Cpw = cfg.value("Cpw", 1.0);
  p->hgw = cfg.value("hgw", 0.0);
  p->hw_inf = cfg.value("hw_inf", 0.0);
  p->T_inf = cfg.value("T_inf", cfg.value("T_init", 298.15));
  const auto pressure_model = lower(cfg.value("PRESSURE_DROP_MODEL", std::string("ergun")));
  if (pressure_model != "ergun" && pressure_model != "darcy")
    throw std::invalid_argument("PRESSURE_DROP_MODEL must be ergun or darcy");
  p->darcy_only = pressure_model == "darcy";
  const auto wall_ends = lower(cfg.value("WALL_END_BOUNDARY", std::string("insulated")));
  if (wall_ends != "insulated" && wall_ends != "fixed_temperature")
    throw std::invalid_argument("WALL_END_BOUNDARY must be insulated or fixed_temperature");
  p->fixed_wall_ends = wall_ends == "fixed_temperature";
  p->wall_end_temperature = arr(cfg.value("WALL_END_TEMPERATURE_K", Json(p->T_inf)), 2, "WALL_END_TEMPERATURE_K");
  for (double v : p->wall_end_temperature)
    if (!std::isfinite(v) || v <= 0) throw std::invalid_argument("Invalid wall-end temperature");
  p->aw_int = cfg.value("aw_int", 4 * dint / std::max(dout * dout - dint * dint, 1e-30));
  p->aw_ext = cfg.value("aw_ext", 4 * dout / std::max(dout * dout - dint * dint, 1e-30));
  p->wall_area = std::max(pi * (dout * dout - dint * dint) / 4, 1e-30);
  p->wall_volume = p->wall_area * L;
  p->T_floor = cfg.value("T_FLOOR", 100.0);
  p->T_ceiling = cfg.value("T_CEILING", 2500.0);
  p->P_ceiling = cfg.value("P_CEILING", 1e9);
  p->C_ceiling = cfg.value("C_GAS_CEILING", 1e7);
  p->V_max = cfg.value("V_MAX_LIMIT", 100.0);
  p->boundary_length = cfg.value("BOUNDARY_HYDRAULIC_LENGTH", .5 * dz);
  p->Q_ceiling = cfg.value("Q_STATE_CEILING", 1e8);
  if (p->mu <= 0 || p->boundary_length <= 0 || p->Cps <= 0)
    throw std::invalid_argument("Column material/hydraulic parameters must be positive");
  Vec cp_default = M;
  for (double& v : cp_default) v *= Cpg;
  p->cp = arr(cfg.value("GAS_COMPONENT_CP_MOLAR", cfg.value("CPG_MOLAR", Json(cp_default))), n_comp,
              "GAS_COMPONENT_CP_MOLAR");
  for (double v : p->cp)
    if (!std::isfinite(v) || v <= (liquid ? 0. : R))
      throw std::invalid_argument(liquid ? "Component heat capacity must be positive"
                                         : "Ideal gas component Cp must exceed R");
  std::tie(p->spatial_method, p->advection_scheme) = spatial_settings(cfg);
  for (auto key : {"COMPRESSIBILITY_FACTOR", "FUGACITY_COEFFICIENT"})
    if (cfg.contains(key))
    {
      Vec a = arr(cfg[key], n_comp, key);
      for (double v : a)
        if (v != 1) throw std::invalid_argument(std::string(key) + " must equal one");
    }
  for (auto key : {"GAS_HEAT_CAPACITY_FUNCTION", "VISCOSITY_FUNCTION", "GAS_THERMAL_CONDUCTIVITY_FUNCTION",
                   "WALL_THERMAL_CONDUCTIVITY_FUNCTION", "WALL_HEAT_CAPACITY_FUNCTION", "DISPERSION_FUNCTION",
                   "COMPRESSIBILITY_FACTOR_FUNCTION", "FUGACITY_COEFFICIENT_FUNCTION"})
    if (cfg.contains(key) && !cfg[key].is_null())
      throw std::invalid_argument(std::string(key) + " Python callable is not a declarative native input");
  layout = {n_comp, Nz, 0, 0};
  p->htf = cfg.value("HTF", Json::object());
  if (!p->htf.is_object()) throw std::invalid_argument("thermal.htf must be a mapping");
  const std::set<std::string> htf_keys = {"model",         "cross_section_area_m2", "volume_m3",
                                          "density_kg_m3", "heat_capacity_J_kg_K",  "initial_temperature_K"};
  for (const auto& [key, value] : p->htf.items())
    if (!htf_keys.count(key)) throw std::invalid_argument("Unknown thermal.htf key: " + key);
  const std::string htf_model = lower(p->htf.value("model", std::string("quasi_steady")));
  if (htf_model != "transient" && htf_model != "quasi_steady")
    throw std::invalid_argument("thermal.htf.model must be transient or quasi_steady");
  layout.transient_htf = htf_model == "transient";
  if (layout.transient_htf)
  {
    if (thermal_mode == 1 || dout <= dint)
      throw std::invalid_argument("Transient HTF requires a non-isothermal column and finite wall");
    if (!std::isfinite(p->rho_w) || !std::isfinite(p->Cpw) || p->rho_w <= 0 || p->Cpw <= 0)
      throw std::invalid_argument("Transient HTF requires finite positive wall density and heat capacity");
    if (p->htf.contains("cross_section_area_m2") == p->htf.contains("volume_m3"))
      throw std::invalid_argument("Transient HTF requires exactly one of cross_section_area_m2 or volume_m3");
    p->jacket_area = p->htf.contains("volume_m3") ? p->htf.at("volume_m3").get<double>() / L
                                                  : p->htf.at("cross_section_area_m2").get<double>();
    p->jacket_rho = p->htf.at("density_kg_m3");
    p->jacket_cp = p->htf.at("heat_capacity_J_kg_K");
    for (double v : {p->jacket_area, p->jacket_rho, p->jacket_cp})
      if (!std::isfinite(v) || v <= 0) throw std::invalid_argument("HTF area, density and Cp must be finite positive");
  }
  for (std::size_t m = 0; m < geom.materials.size(); ++m)
  {
    const Json& mat = geom.materials[m];
    p->main_models.emplace_back(mat.contains("MAIN_RAW_INPUT") ? std::optional<MainMaterialModels>(MainMaterialModels(mat.at("MAIN_RAW_INPUT"),material_names[m])) : std::nullopt);
    const auto& main=p->main_models.back();
    const bool main_chemical=main && main->chemical_sites>0;
    if(main && !main->reactions.empty() && cfg.value("SIMULATION_TYPE",std::string("SwingAdsorption"))!="Breakthrough")
      throw std::invalid_argument("General reaction networks are excluded from cycles");
    p->equilibria.emplace_back(mat, R, n_comp, Nz);
    p->particle_equilibria.emplace_back(mat, R, n_comp, 1);
    const std::string ldf_model = lower(mat.value("LDF_MODEL", std::string("constant")));
    if (ldf_model != "constant" && ldf_model != "pore_diffusion")
      throw std::invalid_argument("LDF_MODEL must be constant or pore_diffusion");
    if (ldf_model == "pore_diffusion")
    {
      double factor = mat.value("LDF_PORE_RATE_S_INV", 0.);
      if (!std::isfinite(factor) || factor <= 0 ||
          (mat.contains("LDF_ACTIVATION_ENERGY_J_MOL") && !mat["LDF_ACTIVATION_ENERGY_J_MOL"].is_null()))
        throw std::invalid_argument("Pore LDF requires positive LDF_PORE_RATE_S_INV and no Arrhenius LDF");
    }
    bool haschem = mat.contains("CHEMISORPTION") && !mat["CHEMISORPTION"].is_null() && !mat["CHEMISORPTION"].empty();
    if(haschem&&main_chemical)throw std::invalid_argument("Use ChemisorptionSites or ChemicalKinetics, not both");
    p->chemistries.emplace_back(
        haschem ? std::optional<ChemistryKinetics>(ChemistryKinetics(mat["CHEMISORPTION"], n_comp)) : std::nullopt);
    double phy = 1, chem = 0;
    Json spec = haschem ? mat["CHEMISORPTION"] : Json::object();
    if (haschem)
    {
      if (spec.contains("capacity_fraction"))
      {
        chem = spec["capacity_fraction"];
        phy = 1 - chem;
        if (chem < 0 || chem > 1 || (spec.contains("equilibrium") && !spec["equilibrium"].empty()))
          throw std::invalid_argument("Invalid chemistry capacity partition");
      }
      else
      {
        chem = 1;
        phy = mat.value("PHYSISORPTION_ENABLED", false) ? 1 : 0;
        if (phy && !spec.value("independent_sites", false))
          throw std::invalid_argument(
              "Combined chemical and physical capacities require independent_sites or "
              "capacity_fraction");
      }
    }
    else if (!mat.value("PHYSISORPTION_ENABLED", true))
      phy = 0;
    const double physical_fraction=mat.value("PHYSISORPTION_FRACTION",cfg.value("PHYSISORPTION_FRACTION",1.0));
    if(!std::isfinite(physical_fraction)||physical_fraction<0.0||physical_fraction>1.0)throw std::invalid_argument("PhysisorptionFraction must lie in [0,1]");
    p->phy_factors.push_back(phy*physical_fraction);
    p->chem_factors.push_back(main_chemical?1.:chem);
    std::vector<std::size_t> groups{layout.n_loading_groups++};
    loading_group_names.push_back(material_names[m] + "/physisorption");
    p->history_start.push_back(layout.n_history_groups);
    p->history_count.push_back(haschem ? p->chemistries[m]->n_history_blocks() : 0);
    if (haschem)
    {
      auto& chem_model = *p->chemistries[m];
      for (std::size_t b = 0; b < chem_model.n_branches(); ++b)
      {
        groups.push_back(layout.n_loading_groups++);
        loading_group_names.push_back(material_names[m] + "/" + chem_model.model() + "/branch_" +
                                      std::to_string(b + 1));
      }
      if (chem_model.n_history_blocks())
      {
        history_group_names.push_back(material_names[m] + "/avrami_adsorption_age");
        history_group_names.push_back(material_names[m] + "/avrami_desorption_age");
      }
      layout.n_history_groups += chem_model.n_history_blocks();
    }
    p->chemical_groups.emplace_back();
    if(main_chemical){
      for(std::size_t site=0;site<main->chemical_sites;++site){
        std::size_t branches=1;
        for(const auto& component:main->components)if(site<component.chemisorption.sites.size())branches=std::max(branches,component.chemisorption.sites[site].branchCount());
        std::vector<std::size_t> site_groups;
        for(std::size_t branch=0;branch<branches;++branch){
          site_groups.push_back(layout.n_loading_groups);groups.push_back(layout.n_loading_groups++);
          loading_group_names.push_back(material_names[m]+"/chemisorption/site_"+std::to_string(site+1)+(branches>1?"/branch_"+std::to_string(branch+1):""));
        }
        p->chemical_groups[m].push_back(site_groups);
      }
      for(std::size_t site=0;site<main->chemical_sites;++site)for(const auto* kind:{"surface","pore","adsorption_age","desorption_age"})
        history_group_names.push_back(material_names[m]+"/"+kind+"/site_"+std::to_string(site+1));
      layout.n_history_groups+=4*main->chemical_sites;
      p->history_count.back()=4*main->chemical_sites;
    }
    Json particle = mat.value("PARTICLE_TRANSPORT", Json::object());
    p->radial_start.push_back(layout.n_history_groups);
    p->particles.push_back(particle);
    p->meshes.push_back(std::nullopt);
    if (!particle.empty())
    {
      if (haschem || main_chemical || (main&&!main->reactions.empty()) || lower(mat.value("CONTACTOR", cfg.value("CONTACTOR", std::string("fixed_bed")))) != "fixed_bed")
        throw std::invalid_argument("Particle transport supports non-reactive physisorption in fixed beds");
      auto model = lower(particle.value("model", std::string("ldf")));
      auto shape = lower(particle.value("geometry", std::string("sphere")));
      if (model != "ldf" && model != "hsdm") throw std::invalid_argument("Particle Model must be LDF or HSDM");
      if (shape != "sphere" && shape != "cylinder_radial" && shape != "immobilized_powder")
        throw std::invalid_argument("Unknown particle geometry");
      double radius = particle.value("radius", .5 * mat.value("dp", cfg.value("dp", .002)));
      if (!std::isfinite(radius) || radius <= 0) throw std::invalid_argument("Particle Radius must be positive");
      p->particles.back()["radius"] = radius;
      if (model == "hsdm")
      {
        if (shape == "immobilized_powder") throw std::invalid_argument("Immobilized powder uses an effective LDF rate");
        std::size_t cells = particle.value("radial_cells", 12);
        if (cells < 2 || cells > 1000) throw std::invalid_argument("RadialCells must be between 2 and 1000");
        p->meshes.back().emplace(cells, shape == "sphere" ? 2 : 1);
        layout.n_history_groups += cells;
        for (std::size_t k = 0; k < cells; ++k)
          history_group_names.push_back(material_names[m] + "/HSDM/radial_" + std::to_string(k));
      }
    }
    p->loading_groups.push_back(groups);
    Vec ads_cp = arr(mat.value("ADSORBED_CP_MOLAR", Json(0.0)), n_comp, "ADSORBED_CP_MOLAR"),
        dh = arr(mat.value("dH_phy", Json(0.0)), n_comp, "dH_phy");
    for (double v : ads_cp)
      if (v < 0) throw std::invalid_argument("ADSORBED_CP_MOLAR cannot be negative");
    p->ads_cp.insert(p->ads_cp.end(), ads_cp.begin(), ads_cp.end());
    Json chemical_eq = spec.value("equilibrium", Json::object());
    if (chemical_eq.is_null()) chemical_eq = Json::object();
    for (std::size_t b = 0; b < groups.size(); ++b)
    {
      double ref = b ? chemical_eq.value("T_ref", spec.value("T_ref", 298.15))
                     : mat.value("ADSORPTION_REFERENCE_TEMPERATURE", 298.15);
      if (ref <= 0) throw std::invalid_argument("Adsorption reference temperature must be positive");

      for (std::size_t j = 0; j < n_comp; ++j)
      {
        double heat = dh[j], reference = ref;
        if (b && p->chemistries[m] && p->chemistries[m]->heat_release()) heat = (*p->chemistries[m]->heat_release())[(b - 1) * n_comp + j];
        if(b&&main_chemical)for(std::size_t site=0;site<p->chemical_groups[m].size();++site)
          for(std::size_t branch=0;branch<p->chemical_groups[m][site].size();++branch)
            if(p->chemical_groups[m][site][branch]==groups[b]&&site<main->components[j].chemisorption.sites.size()) {
              const auto& law=main->components[j].chemisorption.sites[site];
              heat=branch<law.branchCount()?law.branchHeat(branch):law.heatOfChemisorption;
              reference=law.adsorptionReferenceTemperature.value_or(ref);
            }
        p->ads_ref_T.push_back(reference);
        p->ads_ref_u.push_back((liquid ? 0. : p->cp[j] * reference) - heat);
      }
    }
    for (auto key : {"HEAT_OF_ADSORPTION_FUNCTION", "SOLID_HEAT_CAPACITY_FUNCTION", "LDF_RATE_FUNCTION"})
      if (mat.contains(key) && !mat[key].is_null())
        throw std::invalid_argument(std::string(key) + " callable is not a declarative native input");
  }
}
Column::~Column() = default;
Column::Column(Column&&) noexcept = default;
Column& Column::operator=(Column&&) noexcept = default;
Column::Column(const Column& other) : Column(other.cfg)
{
  p->step = other.p->step;
  p->context = other.p->context;
  p->step_set = other.p->step_set;
}
Column& Column::operator=(const Column& other)
{
  if (this != &other)
  {
    Column copy(other);
    *this = std::move(copy);
  }
  return *this;
}
Vec Column::gas_component_heat_capacity(double) const { return p->cp; }
Vec Column::gas_component_enthalpy(double T) const
{
  Vec a = p->cp;
  for (double& v : a) v *= T;
  return a;
}
Vec Column::gas_component_internal_energy(double T) const
{
  Vec a = p->cp;
  for (double& v : a) v = (v - R) * T;
  return a;
}
Vec Column::initial_state(const Json& state_definition)
{
  Json def = state_definition.is_null() ? cfg.value("INITIAL_STATE", Json::object()) : state_definition;
  std::string mode = lower(def.value("mode", std::string("uniform")));
  if (mode == "saved" || mode == "saved_profile" || mode == "warm_start" || mode == "npz")
    throw std::invalid_argument(
        "Saved profile initialization requires native checkpoint state; use --restart or explicit "
        "initialization profiles");
  if (mode != "uniform" && mode != "gas_filled" && mode != "pressurized" && mode != "pressurised" &&
      mode != "evacuated" && mode != "empty" && mode != "custom")
    throw std::invalid_argument("Unsupported INITIAL_STATE mode");
  double defaultP = cfg.value("P_init", 1e5);
  if (mode == "empty")
    defaultP = 0;
  else if (mode == "evacuated")
    defaultP = cfg.value("INITIAL_VACUUM_PRESSURE", defaultP);
  Vec P = profile(def.value("pressure", Json(defaultP)), Nz, "initial pressure");
  double baseT = cfg.value("T_init", 298.15);
  Vec Tg = profile(def.value("gas_temperature", Json(baseT)), Nz, "initial gas temperature"),
      Ts = profile(def.value("solid_temperature", Json(baseT)), Nz, "initial solid temperature"),
      Tw = profile(def.value("wall_temperature", cfg.value("Tw_init", Json(baseT))), Nz, "initial wall temperature");
  for (std::size_t i = 0; i < Nz; ++i)
  {
    if (P[i] < 0 || (P[i] == 0 && thermal_mode != 1) || std::min({Tg[i], Ts[i], Tw[i]}) <= 0)
      throw std::invalid_argument("Invalid initial pressure/temperature");
    if (gas_solid_thermal_equilibrium && std::abs(Tg[i] - Ts[i]) > 1e-8)
      throw std::invalid_argument("Gas-solid thermal equilibrium requires equal initial temperatures");
  }
  Vec y = composition_profile(def.value("composition", cfg.at("y_init")), n_comp, Nz, "initial composition", true),
      c(n_comp * Nz), partial(n_comp * Nz), ctot(Nz), equilibriumP(Nz);
  // Preserve operation order at zero driving force for Avrami age branches.
  for (std::size_t i = 0; i < Nz; ++i)
  {
    ctot[i] = P[i] / (R * Tg[i]);
    equilibriumP[i] = ctot[i] * R * Tg[i];
  }
  for (std::size_t j = 0; j < n_comp; ++j)
    for (std::size_t i = 0; i < Nz; ++i)
    {
      partial[j * Nz + i] = y[j * Nz + i] * equilibriumP[i];
      c[j * Nz + i] = y[j * Nz + i] * ctot[i];
    }
  if (liquid)
  {
    c = composition_profile(
        def.value("concentration", cfg.value("LIQUID_INITIAL_CONCENTRATION", Json(Vec(n_comp, 0.)))), n_comp, Nz,
        "initial liquid concentration", false);
    partial = c;
  }
  {
  }
  Json loading = def.value("loading", cfg.value("INITIAL_LOADING_MODE", Json("zero")));
  bool eq = false;
  Vec q0(n_comp * Nz);
  if (loading.is_string())
  {
    std::string name = lower(loading.get<std::string>());
    eq = name == "eq" || name == "equilibrium";
    if (!eq && name != "zero" && name != "empty")
      throw std::invalid_argument("Initial loading must be zero, equilibrium or explicit array");
  }
  else
    q0 = composition_profile(loading, n_comp, Nz, "initial loading", false);
  Vec state(layout.size(), 0);
  insert(state, c, 0);
  for (std::size_t m = 0; m < p->loading_groups.size(); ++m)
  {
    Vec target = eq ? p->equilibria[m].loading(partial, Ts) : q0;
    const auto& groups = p->loading_groups[m];
    for (std::size_t k = 0; k < target.size(); ++k)
      state[(n_comp + groups[0] * n_comp) * Nz + k] = p->phy_factors[m] * target[k];
    if (p->chemistries[m])
    {
      auto cq = eq ? p->chemistries[m]->equilibrium(y, equilibriumP, Ts) : std::nullopt;
      Vec chemical = cq ? *cq : target;
      for (double& v : chemical) v *= p->chem_factors[m];
      insert(state, p->chemistries[m]->split_loading(chemical, Nz), (n_comp + groups[1] * n_comp) * Nz);
    }
    if(p->main_models[m] && p->main_models[m]->chemical_sites){
      const auto& main=*p->main_models[m];
      for(std::size_t i=0;i<Nz;++i){Vec cb(n_comp);for(std::size_t j=0;j<n_comp;++j)cb[j]=c[j*Nz+i];
        const Vec target=eq?main.chemical_equilibrium(cb,Tg[i],Ts[i]):Vec(main.chemical_sites*n_comp,0.);
        for(std::size_t site=0;site<main.chemical_sites;++site)for(std::size_t j=0;j<n_comp;++j){
          if(site>=main.components[j].chemisorption.sites.size())continue;
          const auto& law=main.components[j].chemisorption.sites[site];
          for(std::size_t branch=0;branch<p->chemical_groups[m][site].size();++branch){
            const auto group=p->chemical_groups[m][site][branch];
            const double fraction=branch<law.branchCount()?law.branchFraction(branch):0.0;
            state[(1+group)*n_comp*Nz+j*Nz+i]=fraction*(eq?target[site*n_comp+j]:q0[j*Nz+i]);
          }
          if(law.usesSurfacePoreTransport()){
            state[(1+layout.n_loading_groups+p->history_start[m]+4*site)*n_comp*Nz+j*Nz+i]=cb[j];
            state[(2+layout.n_loading_groups+p->history_start[m]+4*site)*n_comp*Nz+j*Nz+i]=cb[j];
          }
        }
      }
    }
  }
  for (std::size_t m = 0; m < p->meshes.size(); ++m)
    if (p->meshes[m])
    {
      auto g = p->loading_groups[m][0];
      for (std::size_t k = 0; k < p->meshes[m]->weights.size(); ++k)
        for (std::size_t j = 0; j < n_comp; ++j)
          for (std::size_t i = 0; i < Nz; ++i)
            state[(1 + layout.n_loading_groups + p->radial_start[m] + k) * n_comp * Nz + j * Nz + i] =
                state[(1 + g) * n_comp * Nz + j * Nz + i];
    }
  insert(state, Tg, state.size() - 3 * Nz);
  insert(state, Ts, state.size() - 2 * Nz);
  if (p->fixed_wall_ends)
  {
    Tw.front() = p->wall_end_temperature[0];
    Tw.back() = p->wall_end_temperature[1];
  }
  insert(state, Tw, state.size() - Nz);
  if (layout.transient_htf)
  {
    Vec Tj = profile(def.value("jacket_temperature", p->htf.value("initial_temperature_K", Json(baseT))), Nz,
                     "initial jacket temperature");
    for (double v : Tj)
      if (!std::isfinite(v) || v <= 0) throw std::invalid_argument("Invalid initial jacket temperature");
    insert(state, Tj, state.size() - 4 * Nz);
  }
  return state;
}
Derived Column::derived(const Vec& state) const
{
  if (state.size() != layout.size()) throw std::invalid_argument("Column state size mismatch");
  Derived d;
  d.c = safe_array(slice(state, 0, n_comp * Nz), 0, p->C_ceiling, 0);
  d.q.assign(n_comp * Nz, 0);
  for (std::size_t m = 0; m < p->loading_groups.size(); ++m)
    for (std::size_t i = 0; i < Nz; ++i)
    {
      double total = 0;
      for (std::size_t mm = 0; mm < p->loading_groups.size(); ++mm) total += p->geometry.mass_density[mm * Nz + i];
      double weight = p->geometry.mass_density[m * Nz + i] / std::max(total, 1e-300);
      for (auto g : p->loading_groups[m])
        for (std::size_t j = 0; j < n_comp; ++j) d.q[j * Nz + i] += weight * state[(n_comp + g * n_comp + j) * Nz + i];
    }
  d.q = safe_array(std::move(d.q), -p->Q_ceiling, p->Q_ceiling, 0);
  d.Tg = safe_array(slice(state, state.size() - 3 * Nz, Nz), p->T_floor, p->T_ceiling, p->T_floor);
  d.Ts = safe_array(slice(state, state.size() - 2 * Nz, Nz), p->T_floor, p->T_ceiling, p->T_floor);
  d.Tw = safe_array(slice(state, state.size() - Nz, Nz), p->T_floor, p->T_ceiling, p->T_floor);
  if (layout.transient_htf)
    d.Tj = safe_array(slice(state, state.size() - 4 * Nz, Nz), p->T_floor, p->T_ceiling, p->T_floor);
  d.ctot.assign(Nz, 0);
  d.y.resize(n_comp * Nz);
  d.P.resize(Nz);
  d.rho.assign(Nz, 0);
  Vec y0 = normal(cfg.at("y_init"), n_comp, "y_init");
  for (std::size_t i = 0; i < Nz; ++i)
  {
    for (std::size_t j = 0; j < n_comp; ++j)
    {
      d.ctot[i] += d.c[j * Nz + i];
      d.rho[i] += M[j] * d.c[j * Nz + i];
    }
    for (std::size_t j = 0; j < n_comp; ++j) d.y[j * Nz + i] = d.ctot[i] > 1e-300 ? d.c[j * Nz + i] / d.ctot[i] : y0[j];
    d.P[i] = std::clamp(d.ctot[i] * R * d.Tg[i], 0.0, p->P_ceiling);
    if (liquid)
    {
      d.P[i] = cfg.value("P_init", 1e5);
      d.rho[i] = liquid_density;
    }
  }
  if (liquid && p->cached_liquid_pressure.size() == Nz) d.P = p->cached_liquid_pressure;
  return d;
}
Json Column::endpoint(const std::string& side, const Vec& state) const
{
  if (side != "left" && side != "right") throw std::invalid_argument("Endpoint side must be left or right");
  auto d = derived(state);
  std::size_t i = side == "left" ? 0 : Nz - 1;
  return Json{{"pressure", d.P[i]},
              {"temperature", d.Tg[i]},
              {"composition", components_at(d.y, n_comp, Nz, i)},
              {"concentration", components_at(d.c, n_comp, Nz, i)},
              {"total_concentration", d.ctot[i]},
              {"density", d.rho[i]},
              {"epsilon", epsilon_node[i]},
              {"area", area}};
}
Json Column::source_states(double t) const
{
  Json states = configured_sources(cfg, n_comp, t, p->context);
  if (p->context.contains("source_states"))
    for (auto it = p->context["source_states"].begin(); it != p->context["source_states"].end(); ++it)
      states[it.key()] = it.value();
  return states;
}
void Column::set_step(const Json& step, const Json& context)
{
  if (!step.contains("left_port") || !step.contains("right_port"))
    throw std::invalid_argument("Every cycle step requires left_port and right_port");
  p->step = step;
  if (step.contains("flow_control"))
  {
    const auto& control = step.at("flow_control");
    if (!control.is_object()) throw std::invalid_argument("FlowControl must be an object");
    for (auto& [key, value] : control.items())
      if (key != "model" && key != "pressure_gradient")
        throw std::invalid_argument("Unknown FlowControl setting: " + key);
    const auto model = lower(control.value("model", std::string{}));
    if (model != "fixed_velocity" && model != "fixed_pressure_gradient")
      throw std::invalid_argument("FlowControl.Model must be FixedVelocity or FixedPressureGradient");
    if (model == "fixed_pressure_gradient" && !liquid)
      throw std::invalid_argument("FVM FixedPressureGradient is an incompressible liquid flow mode; gas uses pressure ports or the original main solver");
    for (auto side : {"left_port", "right_port"})
      if (lower(step.at(side).value("mode", std::string{})) != (model == "fixed_velocity" ? "velocity" : "pressure"))
        throw std::invalid_argument("FlowControl requires two velocity ports for FixedVelocity or two pressure ports for FixedPressureGradient");
    if (model == "fixed_pressure_gradient" && !control.contains("pressure_gradient"))
      throw std::invalid_argument("FixedPressureGradient requires PressureGradient in Pa/m");
  }
  p->context = context;
  p->step_set = true;
  port_state("left");
  port_state("right");
}
Json Column::port_state(const std::string& side, double t) const
{
  if (!p->step_set) throw std::runtime_error("set_step must be called first");
  const Json& port = p->step.at(side + "_port");
  std::string mode = lower(port.value("mode", std::string("closed")));
  if (mode == "pressure_reservoir" || mode == "reservoir" || mode == "open") mode = "pressure";
  if (mode == "shut") mode = "closed";
  if (mode != "closed" && mode != "pressure" && mode != "velocity")
    throw std::invalid_argument("Column port mode must be closed, pressure or velocity");
  Json out{{"mode", mode}};
  if (mode == "closed") return out;
  Json source = Json(), states = source_states(t);
  if (port.contains("source") && port.contains("source_mix"))
    throw std::invalid_argument("Port cannot define source and source_mix together");
  if (port.contains("source_mix"))
  {
    const Json& mix = port["source_mix"];
    if (!mix.is_object() || mix.empty()) throw std::invalid_argument("source_mix must be nonempty mapping");
    double total = 0;
    for (auto it = mix.begin(); it != mix.end(); ++it)
    {
      double w = it.value();
      if (!std::isfinite(w) || w < 0) throw std::invalid_argument("Invalid source_mix weights");
      total += w;
    }
    if (total <= 0) throw std::invalid_argument("source_mix weights must have positive sum");
    Vec y(n_comp, 0);
    double pressure = std::numeric_limits<double>::infinity(), enthalpy = 0;
    for (auto it = mix.begin(); it != mix.end(); ++it)
    {
      const Json& s = states.at(it.key());
      if (s.value("kind", std::string{}) == "buffer")
        throw std::invalid_argument("Dynamic buffers must connect directly, not through source_mix");
      double w = it.value().get<double>() / total;
      Vec sy = arr(s.at("composition"), n_comp, "source composition");
      for (std::size_t j = 0; j < n_comp; ++j) y[j] += w * sy[j];
      if (w > 0)
      {
        pressure = std::min(pressure, s.at("pressure").get<double>());
        enthalpy += w * dot(sy, gas_component_enthalpy(s.at("temperature")));
      }
    }
    source = {{"pressure", pressure},
              {"temperature", enthalpy / dot(y, p->cp)},
              {"composition", y},
              {"kind", "mixed_reservoir"}};
  }
  if (port.contains("source"))
  {
    std::string name = port["source"];
    if (states.contains(name))
      source = states[name];
    else if (port.contains("fallback_source") && states.contains(port["fallback_source"].get<std::string>()))
      source = states[port["fallback_source"].get<std::string>()];
    if (source.is_null() && !port.contains("fallback_composition"))
      throw std::invalid_argument("Unavailable port source " + name);
    out["source"] = name;
  }
  double defaultP = cfg.value("P_FEED", cfg.value("P_HIGH", cfg.value("P_init", 1e5))),
         sourceP = source.is_null() ? defaultP : source.value("pressure", defaultP);
  if (mode == "pressure" && !port.contains("pressure") && source.is_null())
    throw std::invalid_argument("Pressure port requires pressure or source");
  Json pressure_value = port.value("pressure", Json(sourceP));
  if (pressure_value.is_object() && pressure_value.contains("initial_value") &&
      pressure_value["initial_value"].is_string() && pressure_value["initial_value"] == "phase_start_pressure")
  {
    if (!p->context.contains("phase_start_pressure") || !p->context["phase_start_pressure"].contains(side))
      throw std::invalid_argument("phase_start_pressure requires original phase snapshot");
    pressure_value["initial_value"] = p->context["phase_start_pressure"][side];
  }
  double pressure = resolve_time_profile(pressure_value, t);
  const auto flow = p->step.value("flow_control", Json::object());
  const auto flow_model = lower(flow.value("model", std::string{}));
  if (flow_model == "fixed_pressure_gradient" && side == "right")
  {
    const double gradient = resolve_time_profile(flow.at("pressure_gradient"), t);
    pressure = port_state("left", t).at("pressure").get<double>() + gradient * L;
  }
  if (!std::isfinite(pressure) || pressure < 0) throw std::invalid_argument("Port pressure must be finite nonnegative");
  Json composition = port.value("composition", source.is_null() ? Json() : source.value("composition", Json()));
  if (port.contains("composition_source"))
  {
    std::string name = port["composition_source"];
    composition = Json();
    if (p->context.contains("stream_compositions") && p->context["stream_compositions"].contains(name))
      composition = p->context["stream_compositions"][name];
    if (composition.is_null() && states.contains(name)) composition = states[name].value("composition", Json());
    if (composition.is_null()) composition = port.value("fallback_composition", Json());
    if (composition.is_null()) throw std::invalid_argument("Unavailable composition_source " + name);
  }
  if (composition.is_null()) composition = port.value("fallback_composition", cfg.at("y_feed"));
  double defaultT = cfg.value("T_feed", cfg.value("T_init", 298.15)),
         sourceT = source.is_null() ? defaultT : source.value("temperature", defaultT),
         temperature = resolve_time_profile(port.value("temperature", Json(sourceT)), t);
  if (!std::isfinite(temperature) || temperature <= 0)
    throw std::invalid_argument("Port temperature must be finite positive");
  out["pressure"] = pressure;
  out["temperature"] = temperature;
  out["composition"] = normal(composition, n_comp, "port composition");
  if (liquid)
  {
    out["concentration"] =
        port.value("concentration", source.is_null() ? cfg.value("LIQUID_FEED_CONCENTRATION", Json(Vec(n_comp, 0.)))
                                                     : source.value("concentration", Json(Vec(n_comp, 0.))));
    Vec c = arr(out["concentration"], n_comp, "liquid reservoir concentration");
    for (double v : c)
      if (!std::isfinite(v) || v < 0) throw std::invalid_argument("Invalid liquid reservoir concentration");
    out["concentration"] = c;
  }
  if (mode == "velocity")
  {
    if (port.contains("velocity") && port.contains("molar_flow_mol_s"))
      throw std::invalid_argument("Velocity port requires velocity OR molar_flow_mol_s");
    double velocity;
    if (port.contains("molar_flow_mol_s"))
    {
      if (liquid) throw std::invalid_argument("Use VolumetricFlowRate for an implicit-solvent liquid feed");
      if (pressure <= 0) throw std::invalid_argument("Molar-flow port requires positive pressure");
      double rate = resolve_time_profile(port.at("molar_flow_mol_s"), t);
      double epsilon = side == "left" ? epsilon_node.front() : epsilon_node.back();
      velocity = rate * R * temperature / (pressure * area * epsilon);
      out["molar_flow_mol_s"] = rate;
    }
    else if (port.contains("volume_flow_m3_s"))
    {
      double eps = side == "left" ? epsilon_node.front() : epsilon_node.back();
      velocity = resolve_time_profile(port.at("volume_flow_m3_s"), t) / (area * eps);
    }
    else
      velocity = resolve_time_profile(port.at("velocity"), t);
    if (!std::isfinite(velocity) || std::abs(velocity) > p->V_max)
      throw std::invalid_argument("Port velocity exceeds configured limit");
    out["velocity"] = velocity;
    if (flow_model == "fixed_velocity" && side == "right")
    {
      const double left_velocity = port_state("left", t).at("velocity").get<double>();
      const bool superficial = cfg.value("VELOCITY_INPUT_BASIS", std::string("interstitial")) == "superficial";
      out["velocity"] = left_velocity * (liquid && !superficial ? epsilon_node.front() / epsilon_node.back() : 1.);
    }
  }
  return out;
}
std::pair<double, double> Column::liquid_resistance() const
{
  double a = 0., b = 0.;
  const auto& g = p->geometry;
  for (std::size_t i = 1; i < Nz; ++i)
  {
    a += .5 * p->mu * (g.resistance_linear[i - 1] + g.resistance_linear[i]) * (z[i] - z[i - 1]) / area;
    if (!p->darcy_only)
      b += .5 * liquid_density * (g.resistance_quadratic[i - 1] + g.resistance_quadratic[i]) * (z[i] - z[i - 1]) /
           (area * area);
  }
  return {a, b};
}
std::pair<Vec, Vec> Column::liquid_hydraulics(double t, const Json& overrides) const
{
  const auto& g = p->geometry;
  Json lp = port_state("left", t), rp = port_state("right", t);
  if (overrides.contains("left") && overrides["left"].contains("pressure"))
    lp["pressure"] = overrides["left"]["pressure"];
  if (overrides.contains("right") && overrides["right"].contains("pressure"))
    rp["pressure"] = overrides["right"]["pressure"];
  double a = 0., b = 0.;
  for (std::size_t i = 1; i < Nz; ++i)
  {
    const double dx = z[i] - z[i - 1];
    a += .5 * p->mu * (g.resistance_linear[i - 1] + g.resistance_linear[i]) * dx;
    if (!p->darcy_only) b += .5 * liquid_density * (g.resistance_quadratic[i - 1] + g.resistance_quadratic[i]) * dx;
  }
  std::optional<double> specified;
  for (const std::string side : {"left", "right"})
  {
    const bool left = side == "left";
    const auto& port = left ? lp : rp;
    double eps = left ? epsilon_node.front() : epsilon_node.back();
    std::optional<double> u;
    if (overrides.contains(side) && overrides[side].contains("volume_flux"))
      u = overrides[side]["volume_flux"].get<double>();
    else if (overrides.contains(side) && overrides[side].contains("velocity"))
      u = eps * overrides[side]["velocity"].get<double>();
    else if (port["mode"] == "closed")
      u = 0.;
    else if (port["mode"] == "velocity")
    {
      u = port["velocity"].get<double>() * eps;
      if (cfg.value("VELOCITY_INPUT_BASIS", std::string("interstitial")) == "superficial")
        *u = port["velocity"].get<double>();
    }
    if (u)
    {
      if (specified && std::abs(*specified - *u) > 1e-10 * std::max({1., std::abs(*specified), std::abs(*u)}))
        throw std::invalid_argument(
            "Incompressible liquid ports require equal axial volume flow; a closed port cannot accompany a nonzero "
            "imposed flow");
      specified = u;
    }
  }
  double U = specified.value_or(hydraulic_root(lp.value("pressure", 0.) - rp.value("pressure", 0.), a, b));
  Vec velocity(Nz + 1), pressure(Nz);
  for (std::size_t f = 0; f <= Nz; ++f) velocity[f] = U / epsilon_face[f];
  pressure[0] = lp.value("pressure", rp.value("pressure", cfg.value("P_init", 1e5)));
  for (std::size_t i = 1; i < Nz; ++i)
  {
    double grad = .5 * p->mu * (g.resistance_linear[i - 1] + g.resistance_linear[i]) * U;
    if (!p->darcy_only)
      grad += .5 * liquid_density * (g.resistance_quadratic[i - 1] + g.resistance_quadratic[i]) * U * std::abs(U);
    pressure[i] = pressure[i - 1] - grad * (z[i] - z[i - 1]);
  }
  if (rp["mode"] == "pressure" && lp["mode"] != "pressure")
  {
    double shift = rp["pressure"].get<double>() - pressure.back();
    for (double& v : pressure) v += shift;
  }
  return {velocity, pressure};
}
Vec Column::face_velocities(const Vec& state, double t, const Json& overrides) const
{
  if (liquid)
  {
    auto flow = liquid_hydraulics(t, overrides);
    p->cached_liquid_pressure = flow.second;
    return flow.first;
  }
  auto d = derived(state);
  if (lower(p->step.value("flow_control", Json::object()).value("model", std::string{})) == "fixed_velocity")
  {
    double velocity = port_state("left", t).at("velocity");
    if (cfg.value("VELOCITY_INPUT_BASIS", std::string("interstitial")) != "interstitial")
      throw std::invalid_argument("Gas FixedVelocity uses interstitial velocity");
    return Vec(Nz + 1, velocity);
  }
  Vec vf(Nz + 1, 0);
  auto& geom = p->geometry;
  for (std::size_t f = 1; f < Nz; ++f)
  {
    double grad = -(d.P[f] - d.P[f - 1]) / (z[f] - z[f - 1]),
           a = .5 * p->mu * (geom.resistance_linear[f - 1] + geom.resistance_linear[f]),
           b = .5 * (geom.resistance_quadratic[f - 1] * d.rho[f - 1] + geom.resistance_quadratic[f] * d.rho[f]);
    vf[f] = hydraulic_root(grad, a, p->darcy_only ? 0. : b) / epsilon_face[f];
  }
  for (const std::string side : {"left", "right"})
  {
    bool left = side == "left";
    std::size_t i = left ? 0 : Nz - 1, f = left ? 0 : Nz;
    Json port = port_state(side, t);
    if (overrides.contains(side))
    {
      const Json& o = overrides[side];
      vf[f] = o.value("velocity", 0.0);
      if (o.contains("component_flux") && !o.contains("velocity"))
      {
        double J = sum(o["component_flux"].get<Vec>());
        double donor = (left ? J >= 0 : J <= 0)
                           ? sum(o.value("concentration", Json(components_at(d.c, n_comp, Nz, i))).get<Vec>())
                           : d.ctot[i];
        vf[f] = J / std::max(epsilon_node[i] * donor, 1e-300);
      }
    }
    else if (port["mode"] == "velocity")
    {
      vf[f] = port["velocity"];
      std::string basis = lower(cfg.value("VELOCITY_INPUT_BASIS", std::string("interstitial")));
      if (!port.contains("molar_flow_mol_s") && (basis == "superficial" || basis == "darcy")) vf[f] /= epsilon_node[i];
    }
    else if (port["mode"] == "pressure")
    {
      double P = port["pressure"], T = port["temperature"], rho = P / (R * T) * dot(port["composition"].get<Vec>(), M),
             grad = (left ? P - d.P[i] : d.P[i] - P) / p->boundary_length;
      vf[f] = hydraulic_root(grad, p->mu * geom.resistance_linear[i],
                             p->darcy_only ? 0. : .5 * (rho + d.rho[i]) * geom.resistance_quadratic[i]) /
              epsilon_node[i];
    }
  }
  return vf;
}
std::pair<Vec, Vec> Column::reservoir_concentrations(double t, const Vec& state, const Json& overrides) const
{
  std::vector<Vec> result;
  for (const std::string side : {"left", "right"})
  {
    Vec c(n_comp, 0);
    if (overrides.contains(side) && overrides[side].contains("concentration"))
    {
      c = arr(overrides[side]["concentration"], n_comp, "boundary concentration");
      for (double v : c)
        if (v < 0) throw std::invalid_argument("Boundary concentration cannot be negative");
    }
    else
    {
      Json port = port_state(side, t);
      if (port["mode"] != "closed")
      {
        c = port["composition"].get<Vec>();
        double pressure = port["pressure"].get<double>();
        const auto basis = cfg.value("INLET_DENSITY_BASIS", std::string("reservoir"));
        if (basis != "reservoir" && basis != "local_pressure")
          throw std::invalid_argument("InletDensityBasis must be Reservoir or LocalPressure");
        if (!liquid && basis == "local_pressure" && port["mode"] == "velocity" && !overrides.contains(side))
        {
          if (port.contains("molar_flow_mol_s"))
            throw std::invalid_argument(
                "LocalPressure inlet density requires a prescribed velocity, not MolarFlowRate");
          const bool left = side == "left";
          const auto node = left ? std::size_t(0) : Nz - 1;
          const auto velocity = face_velocities(state, t, overrides);
          const double U = velocity[left ? 0 : Nz] * epsilon_node[node];
          if (left ? U > 0. : U < 0.)
          {
            const auto d = derived(state);
            const auto& g = p->geometry;
            const double sign = left ? 1. : -1.;
            const double inertia = p->darcy_only ? 0. : .5 * g.resistance_quadratic[node] * U * std::abs(U);
            const double density_per_P = dot(c, M) / (R * port["temperature"].get<double>());
            const double denominator = 1. - sign * p->boundary_length * inertia * density_per_P;
            if (denominator <= 0.)
              throw std::invalid_argument("Prescribed inlet velocity exceeds the local Ergun pressure relation");
            pressure = (d.P[node] +
                        sign * p->boundary_length * (p->mu * g.resistance_linear[node] * U + inertia * d.rho[node])) /
                       denominator;
          }
        }
        double total = pressure / (R * port["temperature"].get<double>());
        for (double& v : c) v *= total;
        if (liquid) c = port.at("concentration").get<Vec>();
      }
    }
    result.push_back(c);
  }
  return {result[0], result[1]};
}
Json Column::boundary_component_rates(double t, const Vec& state, const Json& overrides) const
{
  auto d = derived(state);
  Vec vf = face_velocities(state, t, overrides);
  auto reservoirs = reservoir_concentrations(t, state, overrides);
  Json rates;
  for (const std::string side : {"left", "right"})
  {
    bool left = side == "left";
    std::size_t i = left ? 0 : Nz - 1, f = left ? 0 : Nz;
    Vec in(n_comp, 0), out(n_comp, 0);
    bool incoming = left ? vf[f] >= 0 : vf[f] < 0;
    const Vec& reservoir = left ? reservoirs.first : reservoirs.second;
    for (std::size_t j = 0; j < n_comp; ++j)
    {
      double rate = area * epsilon_node[i] * std::abs(vf[f]);
      if (incoming)
        in[j] = rate * reservoir[j];
      else
        out[j] = rate * d.c[j * Nz + i];
    }
    if (overrides.contains(side) && overrides[side].contains("component_flux"))
    {
      Vec imposed = arr(overrides[side]["component_flux"], n_comp, "component_flux");
      for (std::size_t j = 0; j < n_comp; ++j)
      {
        double rate = area * imposed[j] * (left ? 1 : -1);
        in[j] = std::max(rate, 0.0);
        out[j] = std::max(-rate, 0.0);
      }
    }
    rates[side + "_in"] = in;
    rates[side + "_out"] = out;
  }
  return rates;
}
Json Column::boundary_enthalpy_rates(double t, const Vec& state, const Json& overrides) const
{
  Json rates = boundary_component_rates(t, state, overrides), result;
  auto d = derived(state);
  for (const std::string side : {"left", "right"})
  {
    bool left = side == "left";
    std::size_t i = left ? 0 : Nz - 1;
    Json port = port_state(side, t), o = overrides.value(side, Json::object());
    double incomingT = o.value("temperature", port.value("temperature", d.Tg[i]));
    double hin = dot(rates[side + "_in"].get<Vec>(), gas_component_enthalpy(incomingT)),
           hout = dot(rates[side + "_out"].get<Vec>(), gas_component_enthalpy(d.Tg[i]));
    if (o.contains("enthalpy_flux"))
    {
      double signedflow = area * o["enthalpy_flux"].get<double>() * (left ? 1 : -1);
      hin = std::max(signedflow, 0.0);
      hout = std::max(-signedflow, 0.0);
    }
    if (liquid && !o.contains("enthalpy_flux"))
    {
      auto vf = face_velocities(state, t, overrides);
      double q = area * epsilon_node[i] * vf[left ? 0 : Nz] * (left ? 1. : -1.);
      hin = std::max(q, 0.) * liquid_density * liquid_heat_capacity * incomingT;
      hout = std::max(-q, 0.) * liquid_density * liquid_heat_capacity * d.Tg[i];
    }
    result[side + "_in"] = hin;
    result[side + "_out"] = hout;
  }
  return result;
}
Json Column::thermal_control(double t) const
{
  Json raw = p->step.value("thermal", Json::object());
  if (raw.is_null()) raw = Json::object();
  std::string mode = lower(raw.value("mode", thermal_mode == 1   ? std::string("isothermal")
                                             : thermal_mode == 2 ? std::string("adiabatic")
                                                                 : std::string("ambient")));
  if (mode == "none" || mode == "closed")
    mode = "adiabatic";
  else if (mode == "environment")
    mode = "ambient";
  else if (mode == "coolant" || mode == "heating_jacket" || mode == "cooling_jacket")
    mode = "jacket";
  else if (mode == "heater")
    mode = "electric";
  else if (mode == "specified_heat_flux")
    mode = "heat_flux";
  else if (mode == "heat_transfer_fluid" || mode == "cooling_fluid")
    mode = "htf";
  if (mode != "isothermal" && mode != "adiabatic" && mode != "ambient" && mode != "jacket" && mode != "heat_flux" &&
      mode != "electric" && mode != "htf")
    throw std::invalid_argument("Unknown thermal mode");
  if (thermal_mode == 1 && mode != "isothermal" && mode != "adiabatic")
    throw std::invalid_argument("Active thermal controls require non-isothermal mode");
  const Json temperature = raw.value("temperature", Json(p->T_inf));
  Vec temperature_profile(Nz);
  if (temperature.is_object() && temperature.value("model", std::string()) == "axial_linear")
  {
    for (const auto& [key, value] : temperature.items())
      if (key != "model" && key != "left" && key != "right")
        throw std::invalid_argument("Unknown axial_linear temperature key: " + key);
    const double a = resolve_time_profile(temperature.at("left"), t),
                 b = resolve_time_profile(temperature.at("right"), t);
    for (std::size_t i = 0; i < Nz; ++i) temperature_profile[i] = a + (b - a) * z[i] / L;
    if (mode == "htf" || mode == "isothermal")
      throw std::invalid_argument("axial_linear temperature is a prescribed jacket boundary");
  }
  else
  {
    const double uniform = resolve_time_profile(temperature, t);
    std::fill(temperature_profile.begin(), temperature_profile.end(), uniform);
  }
  for (double v : temperature_profile)
    if (!std::isfinite(v) || v <= 0) throw std::invalid_argument("Invalid prescribed jacket temperature");
  double T = dot(temperature_profile, widths) / L;
  Json hraw = raw.value("heat_transfer_coefficient_W_m2_K",
                        raw.value("heat_transfer_coefficient", Json(mode == "ambient" ? p->hw_inf : 0)));
  Vec hprofile(Nz);
  if (hraw.is_array() && hraw.size() != Nz)
    throw std::invalid_argument("Axial wall-to-fluid coefficient requires one entry per bed node");
  for (std::size_t i = 0; i < Nz; ++i)
  {
    hprofile[i] = resolve_time_profile(hraw.is_array() ? hraw[i] : hraw, t);
    if (!std::isfinite(hprofile[i]) || hprofile[i] < 0)
      throw std::invalid_argument("Wall-to-fluid coefficients must be finite nonnegative");
  }
  double h = dot(hprofile, widths) / L;
  double bedh = resolve_time_profile(raw.value("bed_wall_heat_transfer_W_m2_K", Json(p->hgw)), t),
         power = resolve_time_profile(raw.value("power_W", Json(0.0)), t),
         flux = resolve_time_profile(raw.value("heat_flux_W_m2", Json(0.0)), t);
  if (!std::isfinite(T) || T <= 0 || !std::isfinite(h) || h < 0 || !std::isfinite(bedh) || bedh < 0 ||
      !std::isfinite(power) || !std::isfinite(flux))
    throw std::invalid_argument("Invalid thermal control parameters");
  if (mode == "isothermal" || mode == "adiabatic")
  {
    h = power = flux = 0;
    std::fill(hprofile.begin(), hprofile.end(), 0.);
  }
  Json out{{"mode", mode},
           {"temperature", T},
           {"temperature_profile", temperature_profile},
           {"heat_transfer_coefficient", h},
           {"heat_transfer_coefficient_profile", hprofile},
           {"bed_wall_heat_transfer_coefficient", bedh},
           {"power_W", power},
           {"heat_flux_W_m2", flux},
           {"utility", raw.value("utility", mode)}};
  if (mode == "htf")
  {
    double mdot = resolve_time_profile(raw.value("mass_flow_rate_kg_s", raw.value("mass_flow_rate", Json(0.0))), t),
           cp = resolve_time_profile(
               raw.value("heat_capacity_J_kg_K",
                         raw.value("heat_capacity", Json(layout.transient_htf ? p->jacket_cp : 4180.0))),
               t);
    std::string direction = lower(raw.value("flow_direction", std::string("co_current")));
    if (!std::isfinite(mdot) || !std::isfinite(cp) || mdot < 0 || cp <= 0 || (!layout.transient_htf && mdot == 0) ||
        (direction != "co_current" && direction != "cocurrent" && direction != "counter_current" &&
         direction != "countercurrent"))
      throw std::invalid_argument("Invalid HTF flow parameters");
    if (layout.transient_htf && std::abs(cp - p->jacket_cp) > 1e-10 * p->jacket_cp)
      throw std::invalid_argument("Transient HTF uses constant Cp from thermal.htf; step Cp must match");
    out["mass_flow_rate_kg_s"] = mdot;
    out["heat_capacity_J_kg_K"] = cp;
    out["flow_direction"] = direction;
    out["model"] = layout.transient_htf ? "transient" : "quasi_steady";
  }
  return out;
}
std::pair<Vec, Json> Column::external_wall_heat_density(const Vec& Tw, double t, const Vec& Tj) const
{
  Json control = thermal_control(t);
  Vec density(Nz, 0);
  double T = control["temperature"];
  Vec h = control["heat_transfer_coefficient_profile"].get<Vec>();
  if (control["mode"] == "htf")
  {
    Vec temperatures(Nz), inlets(Nz), derivative(Nz, 0.);
    double wall_heat = 0.;
    double mcp = control["mass_flow_rate_kg_s"].get<double>() * control["heat_capacity_J_kg_K"].get<double>(),
           fluid = T;
    std::string direction = control["flow_direction"];
    bool reverse = direction == "counter_current" || direction == "countercurrent";
    for (std::size_t j = 0; j < Nz; ++j)
    {
      std::size_t i = reverse ? Nz - 1 - j : j;
      inlets[i] = fluid;
      double ua = h[i] * pi * dout * widths[i], outlet, heat;
      if (layout.transient_htf)
      {
        if (Tj.size() != Nz) throw std::invalid_argument("Missing transient HTF state");
        outlet = Tj[i];
        heat = ua * (Tj[i] - Tw[i]);
        temperatures[i] = Tj[i];
        derivative[i] = (mcp * (fluid - outlet) - heat) / (p->jacket_rho * p->jacket_cp * p->jacket_area * widths[i]);
      }
      else
      {
        // expm1 avoids cancellation at small NTU/high HTF flow rate.
        heat = -mcp * (fluid - Tw[i]) * std::expm1(-ua / mcp);
        outlet = fluid - heat / mcp;
        temperatures[i] = ua > 0 ? Tw[i] + heat / ua : fluid;
      }
      wall_heat += heat;
      density[i] = heat / std::max(p->wall_area * widths[i], 1e-300);
      fluid = outlet;
    }
    control["outlet_temperature"] = fluid;
    control["temperature_profile"] = temperatures;
    control["inlet_temperature_profile"] = inlets;
    control["temperature_derivative_K_s"] = derivative;
    control["wall_heat_W"] = wall_heat;
    control["boundary_heat_W"] = layout.transient_htf ? mcp * (T - fluid) : wall_heat;
    control["storage_rate_W"] = layout.transient_htf ? mcp * (T - fluid) - wall_heat : 0.;
  }
  else
  {
    const Vec temperatures = control["temperature_profile"].get<Vec>();
    for (std::size_t i = 0; i < Nz; ++i) density[i] = h[i] * p->aw_ext * (temperatures[i] - Tw[i]);
  }
  double direct =
      control["power_W"].get<double>() / p->wall_volume + control["heat_flux_W_m2"].get<double>() * p->aw_ext;
  for (double& v : density) v += direct;
  return {density, control};
}
Json Column::jacket_record(double t, const Vec& state) const
{
  auto d = derived(state);
  auto record = external_wall_heat_density(d.Tw, t, d.Tj).second;
  if (layout.transient_htf && record["mode"] != "htf") record["temperature_profile"] = d.Tj;
  return record;
}
Vec Column::rhs(double t, const Vec& state, const Json& overrides) const
{
  if (!p->step_set) throw std::runtime_error("set_step must precede integration");
  auto d = derived(state);
  const auto& geom = p->geometry;
  std::size_t ncells = n_comp * Nz;
  Vec result(layout.size(), 0), partial(ncells), dloads(layout.n_loading_groups * ncells, 0),
      dhistory(layout.n_history_groups * ncells, 0);
  Vec reaction_dloads(dloads.size(),0.),extra_sink(ncells,0.),pore_capacity(Nz,0.),pore_energy_rate(Nz,0.),reaction_heat(Nz,0.);
  p->reaction_rates.assign(n_comp,0.);p->reaction_heat=0.;
  for (std::size_t j = 0; j < n_comp; ++j)
    for (std::size_t i = 0; i < Nz; ++i)
      partial[j * Nz + i] = liquid ? d.c[j * Nz + i] : d.y[j * Nz + i] * d.ctot[i] * R * d.Tg[i];
  for (std::size_t m = 0; m < geom.materials.size(); ++m)
  {
    const Json& mat = geom.materials[m];
    const auto& groups = p->loading_groups[m];
    Vec qeq = p->equilibria[m].loading(partial, d.Ts), rates = arr(mat.value("km_phy", Json(0.0)), n_comp, "km_phy");
    bool activated = mat.contains("LDF_ACTIVATION_ENERGY_J_MOL") && !mat["LDF_ACTIVATION_ENERGY_J_MOL"].is_null();
    Vec energy = arr(activated ? mat["LDF_ACTIVATION_ENERGY_J_MOL"] : Json(0.0), n_comp, "LDF_ACTIVATION_ENERGY_J_MOL");
    double tref = mat.value("LDF_REFERENCE_TEMPERATURE", 298.15);
    const std::string ldf_model = lower(mat.value("LDF_MODEL", std::string("constant")));
    if (ldf_model != "constant" && ldf_model != "pore_diffusion")
      throw std::invalid_argument("LDF_MODEL must be constant or pore_diffusion");
    Vec ratio_c, ratio_q;
    double pore_rate = 0, particle_density = 1;
    if (ldf_model == "pore_diffusion")
    {
      if (activated) throw std::invalid_argument("Pore LDF cannot also use Arrhenius LDF");
      pore_rate = mat.at("LDF_PORE_RATE_S_INV").get<double>();
      particle_density = mat.at("rho_s").get<double>();
      if (!std::isfinite(pore_rate) || pore_rate <= 0 || !std::isfinite(particle_density) || particle_density <= 0)
        throw std::invalid_argument("Pore LDF prefactor and density must be finite positive");
      // Evaluate C/q* with a dilute probe at zero concentration, while keeping
      // the actual equilibrium target. This avoids a 0/0 rate at pure feeds.
      Vec probe = partial;
      ratio_c.resize(ncells);
      for (std::size_t j = 0; j < n_comp; ++j)
        for (std::size_t i = 0; i < Nz; ++i)
        {
          const auto ix = j * Nz + i;
          probe[ix] = std::max(probe[ix], R * d.Ts[i] * 1e-12);
          ratio_c[ix] = probe[ix] / (R * d.Ts[i]);
        }
      ratio_q = p->equilibria[m].loading(probe, d.Ts);
    }
    for (std::size_t j = 0; j < n_comp; ++j)
      for (std::size_t i = 0; i < Nz; ++i)
      {
        double k = rates[j];
        if (ldf_model == "pore_diffusion")
        {
          const auto ix = j * Nz + i;
          k = ratio_q[ix] > 0 ? pore_rate * ratio_c[ix] / (particle_density * ratio_q[ix]) : 0.;
        }
        if (activated) k *= std::exp(std::clamp(-energy[j] / R * (1 / d.Ts[i] - 1 / tref), -700.0, 700.0));
        if (!std::isfinite(k) || k < 0) throw std::invalid_argument("LDF rates must be finite nonnegative");
        std::size_t index = groups[0] * ncells + j * Nz + i;
        dloads[index] = k * (p->phy_factors[m] * qeq[j * Nz + i] - state[ncells + index]);
      }
    if (!p->particles[m].empty())
    {
      const auto& spec = p->particles[m];
      const auto& mesh = p->meshes[m];
      double radius = spec.at("radius"), tref = spec.value("reference_temperature", 298.15), rho = mat.at("rho_s");
      Vec ds = arr(spec.value("surface_diffusivity", Json(0.)), n_comp, "SurfaceDiffusivity"),
          ea = arr(spec.value("diffusivity_activation_energy", Json(0.)), n_comp, "DiffusivityActivationEnergy"),
          film = arr(spec.value("film_coefficient", Json(0.)), n_comp, "FilmCoefficient"),
          rate_ea = arr(spec.value("rate_activation_energy", Json(0.)), n_comp, "RateActivationEnergy");
      for (std::size_t i = 0; i < Nz; ++i)
      {
        Vec cb(n_comp), edge(n_comp), coeff(n_comp), k(n_comp), D(n_comp);
        for (std::size_t j = 0; j < n_comp; ++j)
        {
          cb[j] = d.c[j * Nz + i];
          D[j] = ruptura2::arrhenius_from_reference(ds[j], ea[j], d.Ts[i], tref, R);
          k[j] = ruptura2::arrhenius_from_reference(rates[j], rate_ea[j], d.Ts[i], tref, R);
          if (spec.value("ldf_coefficient_source", std::string("specified")) == "diffusivity")
            k[j] = (lower(spec.value("geometry", std::string("sphere"))) == "sphere" ? 15. : 8.) * D[j] /
                   (radius * radius);
          auto meanidx = (1 + groups[0]) * ncells + j * Nz + i;
          edge[j] = mesh
                        ? state[(1 + layout.n_loading_groups + p->radial_start[m] + mesh->weights.size() - 1) * ncells +
                                j * Nz + i]
                        : state[meanidx];
          // Extend the nonlinear film closure continuously for Newton trial
          // states. The radial balance below still uses the unmodified shell
          // states; accepted-state validation never projects the trajectory.
          edge[j] = std::max(edge[j], 0.);
          double factor =
              mesh ? D[j] / (radius * (1. - mesh->centers.back()))
                   : radius * k[j] / (lower(spec.value("geometry", std::string("sphere"))) == "sphere" ? 3. : 2.);
          coeff[j] = rho * factor;
        }
        auto equilibrium = [&](std::span<const double> c)
        {
          Vec driving(n_comp);
          for (std::size_t j = 0; j < n_comp; ++j) driving[j] = c[j] * (liquid ? 1. : R * d.Tg[i]);
          auto out = p->particle_equilibria[m].loading(driving, {d.Ts[i]});
          for (auto& value : out) value *= p->phy_factors[m];
          return out;
        };
        Vec surface = spec.value("film_resistance", false)
                          ? ruptura2::solve_particle_surface(cb, edge, coeff, film, equilibrium).equilibrium_loading
                          : equilibrium(cb);
        for (std::size_t j = 0; j < n_comp; ++j)
        {
          auto meanidx = groups[0] * ncells + j * Nz + i;
          if (mesh)
          {
            Vec q(mesh->weights.size()), dq(q.size());
            for (std::size_t r = 0; r < q.size(); ++r)
              q[r] = state[(1 + layout.n_loading_groups + p->radial_start[m] + r) * ncells + j * Nz + i];
            ruptura2::particle_diffusion_rhs(*mesh, q, surface[j], D[j] / (radius * radius), dq);
            dloads[meanidx] = ruptura2::particle_mean(*mesh, dq);
            for (std::size_t r = 0; r < q.size(); ++r)
              dhistory[(p->radial_start[m] + r) * ncells + j * Nz + i] = geom.fractions[m * Nz + i] > 0 ? dq[r] : 0.;
          }
          else
            dloads[meanidx] = k[j] * (surface[j] - state[ncells + meanidx]);
        }
      }
    }
    if (p->chemistries[m])
    {
      auto& chem = *p->chemistries[m];
      Vec P(Nz);
      for (std::size_t i = 0; i < Nz; ++i) P[i] = d.ctot[i] * R * d.Tg[i];
      auto cq = chem.equilibrium(d.y, P, d.Ts);
      Vec target = cq ? *cq : qeq;
      for (double& v : target) v *= p->chem_factors[m];
      auto evaluated = chem.evaluate(
          target, slice(state, ncells + groups[1] * ncells, chem.n_branches() * ncells),
          slice(state, ncells * (1 + layout.n_loading_groups + p->history_start[m]), chem.n_history_blocks() * ncells),
          d.Ts);
      insert(dloads, evaluated.first, groups[1] * ncells);
      insert(dhistory, evaluated.second, p->history_start[m] * ncells);
      for (std::size_t h = 0; h < chem.n_history_blocks(); ++h)
        for (std::size_t j = 0; j < n_comp; ++j)
          for (std::size_t i = 0; i < Nz; ++i)
            if (geom.fractions[m * Nz + i] <= 0) dhistory[(p->history_start[m] + h) * ncells + j * Nz + i] = 0;
    }
    if(p->main_models[m]) {
      const auto& main=*p->main_models[m];const auto ns=main.chemical_sites;
      for(std::size_t i=0;i<Nz;++i) {
        if(geom.fractions[m*Nz+i]<=0.)continue;
        const double mass=geom.mass_density[m*Nz+i],solid_volume=mass/main.particle_density;
        Vec cb(n_comp),phy(n_comp),phy_dot(n_comp),cq(ns*n_comp,0.0),cq_dot(ns*n_comp,0.0),pore(ns*n_comp),pore_dot(ns*n_comp,0.0);
        for(std::size_t j=0;j<n_comp;++j){cb[j]=d.c[j*Nz+i];phy[j]=state[(1+groups[0])*ncells+j*Nz+i];phy_dot[j]=dloads[groups[0]*ncells+j*Nz+i];}
        const Vec eq=main.chemical_equilibrium(cb,d.Tg[i],d.Ts[i]);
        for(std::size_t site=0;site<ns;++site)for(std::size_t j=0;j<n_comp;++j){
          const auto ix=site*n_comp+j,hi=p->history_start[m]+4*site;
          const auto& site_groups=p->chemical_groups[m][site];
          std::array<double,2> loadings{};
          for(std::size_t branch=0;branch<site_groups.size();++branch){
            loadings[branch]=state[(1+site_groups[branch])*ncells+j*Nz+i];cq[ix]+=loadings[branch];
          }
          pore[ix]=state[(2+layout.n_loading_groups+hi)*ncells+j*Nz+i];
          if(site>=main.components[j].chemisorption.sites.size())continue;
          const auto& law=main.components[j].chemisorption.sites[site];
          const double cs=state[(1+layout.n_loading_groups+hi)*ncells+j*Nz+i];
          auto evaluated=law.evaluate(eq[ix],std::span<const double>(loadings.data(),law.branchCount()),
              law.usesSurfacePoreTransport()?pore[ix]:cb[j],d.Ts[i],
              state[(3+layout.n_loading_groups+hi)*ncells+j*Nz+i],
              state[(3+layout.n_loading_groups+hi)*ncells+j*Nz+i],
              state[(4+layout.n_loading_groups+hi)*ncells+j*Nz+i]);
          if(!law.usesEquilibriumLoading()){
            auto& rate=evaluated.loading[0];
            if(law.usesSurfacePoreTransport()&&rate>0.)rate=std::min(rate,std::max(pore[ix],0.)/(main.particle_density*1e-4));
            if(cq[ix]>=eq[ix])rate=std::min(rate,0.);else if(cq[ix]<=0.)rate=std::max(rate,0.);
          }
          for(std::size_t branch=0;branch<site_groups.size();++branch){
            const double rate=branch<law.branchCount()?evaluated.loading[branch]:0.0;
            dloads[site_groups[branch]*ncells+j*Nz+i]=rate;cq_dot[ix]+=rate;
          }
          if(law.type==Chemisorption::Type::Avrami){
            dhistory[(hi+2)*ncells+j*Nz+i]=evaluated.adsorptionAge;
            dhistory[(hi+3)*ncells+j*Nz+i]=evaluated.desorptionAge;
          }
          if(law.usesSurfacePoreTransport()) {
            const double length=std::max(main.geometry.dimensions.poreDiffusionLength,1e-30);
            const double film=main.geometry.contactAreas.solidFluidPerSolidVolume*law.filmMassTransferCoefficient*(cb[j]-cs);
            const double diff=15.*law.poreDiffusivity/(length*length)*(cs-pore[ix]);
            dhistory[hi*ncells+j*Nz+i]=film-diff;pore_dot[ix]=diff-main.particle_density*cq_dot[ix];
            extra_sink[j*Nz+i]+=solid_volume*film-mass*cq_dot[ix];
          }
        }
        Vec rp(n_comp,0.),rc(ns*n_comp,0.),rpo(ns*n_comp,0.),rb(n_comp,0.),rh(1,0.);
        if(!main.reactions.empty()) {
          const Vec T{d.Ts[i]},fraction{1.};
          ::computeReactionDerivatives(main.components,main.reactions,0,n_comp,ns,d.Ts[i],phy,phy_dot,cq,cq_dot,pore,pore_dot,T,rp,rc,rpo,rh,cb,rb,fraction,1);
        }
        for(std::size_t j=0;j<n_comp;++j){
          dloads[groups[0]*ncells+j*Nz+i]=phy_dot[j];reaction_dloads[groups[0]*ncells+j*Nz+i]=rp[j];
          const double bulk=storage_epsilon_node[i]*geom.fractions[m*Nz+i]*rb[j];extra_sink[j*Nz+i]-=bulk;
          double produced=mass*rp[j]+bulk;
          for(std::size_t site=0;site<ns;++site){
            const auto ix=site*n_comp+j,hi=p->history_start[m]+4*site;
            const auto& site_groups=p->chemical_groups[m][site];
            for(std::size_t branch=0;branch<site_groups.size();++branch){
              double fraction=branch==0?1.0:0.0;
              if(rc[ix]<0.0&&cq[ix]>0.0)fraction=std::max(state[(1+site_groups[branch])*ncells+j*Nz+i],0.0)/cq[ix];
              else if(site<main.components[j].chemisorption.sites.size()) {
                const auto& law=main.components[j].chemisorption.sites[site];fraction=branch<law.branchCount()?law.branchFraction(branch):0.0;
              }
              const double reaction=rc[ix]*fraction;
              dloads[site_groups[branch]*ncells+j*Nz+i]+=reaction;reaction_dloads[site_groups[branch]*ncells+j*Nz+i]=reaction;
            }
            dhistory[(hi+1)*ncells+j*Nz+i]=pore_dot[ix];produced+=mass*rc[ix]+solid_volume*rpo[ix];
            if(!liquid){
              const double csum=state[(1+layout.n_loading_groups+hi)*ncells+j*Nz+i]+pore[ix];
              pore_capacity[i]+=solid_volume*(p->cp[j]-R)*csum;
              pore_energy_rate[i]+=solid_volume*(p->cp[j]-R)*d.Ts[i]*(dhistory[hi*ncells+j*Nz+i]+pore_dot[ix]);
            }
          }
          p->reaction_rates[j]+=area*widths[i]*produced;
        }
        reaction_heat[i]+=mass*rh[0];p->reaction_heat+=area*widths[i]*mass*rh[0];
      }
    }
    double regularizer = mat.value("KINETIC_CONCENTRATION_REGULARIZER", 1e-10);
    if (regularizer <= 0) throw std::invalid_argument("KINETIC_CONCENTRATION_REGULARIZER must be positive");
    for (auto group : groups)
      for (std::size_t j = 0; j < n_comp; ++j)
        for (std::size_t i = 0; i < Nz; ++i)
        {
          std::size_t index = group * ncells + j * Nz + i;
          if (p->particles[m].empty() && !p->main_models[m] && dloads[index] > 0)
            dloads[index] *= d.c[j * Nz + i] / (d.c[j * Nz + i] + regularizer);
          if (geom.fractions[m * Nz + i] <= 0) dloads[index] = 0;
        }
  }
  Vec sink(ncells, 0), ads_capacity(Nz, 0), ads_energy_rate(Nz, 0);
  for (std::size_t m = 0; m < geom.materials.size(); ++m)
    for (auto group : p->loading_groups[m])
      for (std::size_t j = 0; j < n_comp; ++j)
        for (std::size_t i = 0; i < Nz; ++i)
        {
          std::size_t index = group * ncells + j * Nz + i;
          double mass = geom.mass_density[m * Nz + i], cp = p->ads_cp[m * n_comp + j];
          sink[j * Nz + i] += mass * (dloads[index]-reaction_dloads[index]);
          ads_capacity[i] += mass * state[ncells + index] * cp;
          ads_energy_rate[i] +=
              mass * (p->ads_ref_u[group * n_comp + j] + cp * (d.Ts[i] - p->ads_ref_T[group * n_comp + j])) * dloads[index];
        }
  for(std::size_t k=0;k<ncells;++k)sink[k]+=extra_sink[k];
  for(std::size_t i=0;i<Nz;++i){ads_capacity[i]+=pore_capacity[i];ads_energy_rate[i]+=pore_energy_rate[i]-reaction_heat[i];}
  Vec vf = face_velocities(state, t, overrides), volume_flux(Nz + 1), dispersion(ncells, 0);
  for (std::size_t f = 0; f <= Nz; ++f) volume_flux[f] = vf[f] * epsilon_face[f];
  auto reservoirs = reservoir_concentrations(t, state, overrides);
  for (std::size_t m = 0; m < geom.materials.size(); ++m)
  {
    Vec coefficient = arr(geom.materials[m].value("Dj", cfg.value("Dj", Json(0.0))), n_comp, "Dj");
    for (std::size_t j = 0; j < n_comp; ++j)
      for (std::size_t i = 0; i < Nz; ++i) dispersion[j * Nz + i] += geom.fractions[m * Nz + i] * coefficient[j];
  }
  Vec left_flux, right_flux;
  if (overrides.contains("left") && overrides["left"].contains("component_flux"))
    left_flux = arr(overrides["left"]["component_flux"], n_comp, "left component_flux");
  if (overrides.contains("right") && overrides["right"].contains("component_flux"))
    right_flux = arr(overrides["right"]["component_flux"], n_comp, "right component_flux");
  auto transport = conservative_species_face_flux_rhs(
      d.c, volume_flux, sink, storage_epsilon_node, dispersion, reservoirs.first, reservoirs.second, z, epsilon_face,
      p->advection_scheme, geom.interfaces, widths, cfg.value("DISPERSION_CLOSURE", std::string("concentration")),
      left_flux, right_flux, p->spatial_method);
  const Vec& dc = transport.first;
  const Vec& component_flux = transport.second;
  insert(result, dc, 0);
  insert(result, dloads, ncells);
  insert(result, dhistory, (1 + layout.n_loading_groups) * ncells);
  if (thermal_mode != 1)
  {
    Json left = port_state("left", t), right = port_state("right", t), ol = overrides.value("left", Json::object()),
         orr = overrides.value("right", Json::object());
    double leftT = ol.value("temperature", left.value("temperature", d.Tg.front())),
           rightT = orr.value("temperature", right.value("temperature", d.Tg.back()));
    Vec faceT = reconstructed_face_values(d.Tg, z, volume_flux, Vec{leftT}, Vec{rightT}, p->advection_scheme,
                                          geom.interfaces, p->spatial_method),
        enthalpy_flux(Nz + 1, 0);
    for (std::size_t j = 0; j < n_comp; ++j)
      for (std::size_t f = 0; f <= Nz; ++f) enthalpy_flux[f] += component_flux[j * (Nz + 1) + f] * p->cp[j] * faceT[f];
    if (liquid)
      for (std::size_t f = 0; f <= Nz; ++f)
        enthalpy_flux[f] = volume_flux[f] * liquid_density * liquid_heat_capacity * faceT[f];
    if (ol.contains("enthalpy_flux")) enthalpy_flux.front() = ol["enthalpy_flux"];
    if (orr.contains("enthalpy_flux")) enthalpy_flux.back() = orr["enthalpy_flux"];
    Vec kg(Nz);
    for (std::size_t i = 0; i < Nz; ++i) kg[i] = epsilon_node[i] * std::max(p->lambda_g, 1e-300);
    Vec conduction = conductive_divergence(d.Tg, kg, z, widths), kw(Nz, std::max(p->lambda_w, 1e-300)),
        wall_conduction = conductive_divergence(d.Tw, kw, z, widths);
    auto external = external_wall_heat_density(d.Tw, t, d.Tj);
    if (layout.transient_htf && external.second["mode"] == "htf")
      insert(result, external.second["temperature_derivative_K_s"].get<Vec>(), result.size() - 4 * Nz);
    double bedh = thermal_control(t)["bed_wall_heat_transfer_coefficient"];
    for (std::size_t i = 0; i < Nz; ++i)
    {
      double exchange_gs = geom.hgs_as[i] * (d.Tg[i] - d.Ts[i]), exchange_gw = 4 * bedh / dint * (d.Tg[i] - d.Tw[i]),
             uptake = 0, capacity = 0, internal_mass = 0;
      for (std::size_t j = 0; j < n_comp; ++j)
      {
        std::size_t index = j * Nz + i;
        uptake += sink[index] * p->cp[j] * d.Tg[i];
        capacity += storage_epsilon_node[i] * d.c[index] * (p->cp[j] - R);
        internal_mass += storage_epsilon_node[i] * (p->cp[j] - R) * d.Tg[i] * dc[index];
      }
      if (liquid)
      {
        capacity = storage_epsilon_node[i] * liquid_density * liquid_heat_capacity;
        internal_mass = 0.;
        uptake = 0.;
      }
      double gasrhs = -(enthalpy_flux[i + 1] - enthalpy_flux[i]) / widths[i] + conduction[i] - exchange_gs -
                      exchange_gw - uptake,
             dTg = (gasrhs - internal_mass) / std::max(capacity, 1e-12),
             dTs =
                 (exchange_gs + uptake - ads_energy_rate[i]) / std::max(geom.heat_capacity[i] + ads_capacity[i], 1e-12);
      if (gas_solid_thermal_equilibrium)
      {
        double value = -(enthalpy_flux[i + 1] - enthalpy_flux[i]) / widths[i] + conduction[i] - exchange_gw -
                       ads_energy_rate[i] - internal_mass;
        dTg = value / std::max(capacity + geom.heat_capacity[i] + ads_capacity[i], 1e-12);
        dTs = dTg;
      }
      result[result.size() - 3 * Nz + i] = dTg;
      result[result.size() - 2 * Nz + i] = dTs;
      result[result.size() - Nz + i] =
          (wall_conduction[i] + exchange_gw * area / p->wall_area + external.first[i]) / (p->rho_w * p->Cpw);
      // Boundary nodes are thermostat-controlled. The compensating heat is
      // accounted for explicitly in thermal_utility_rates, not discarded.
      if (p->fixed_wall_ends && (i == 0 || i == Nz - 1)) result[result.size() - Nz + i] = 0.;
    }
  }
  return result;
}
Vec Column::explicit_transport_rhs(double t, const Vec& state, const Json& overrides) const
{
  // Pressure-driven compressible transport is itself a stiff pressure mode.
  // Keep it implicit; prescribed-velocity gas and incompressible solute
  // convection form the explicit block of the shared additive RK method.
  if(!liquid && lower(p->step.value("flow_control",Json::object()).value("model",std::string{}))!="fixed_velocity")
    return Vec(state_size(),0.);
  const auto d=derived(state);
  const auto& geom=p->geometry;
  const std::size_t ncells=n_comp*Nz;
  Vec result(state_size(),0.),flow=face_velocities(state,t,overrides);
  for(std::size_t f=0;f<flow.size();++f)flow[f]*=epsilon_face[f];
  const auto reservoirs=reservoir_concentrations(t,state,overrides);
  Vec left_flux,right_flux;
  const Json ol=overrides.value("left",Json::object()),orr=overrides.value("right",Json::object());
  if(ol.contains("component_flux"))left_flux=ol.at("component_flux").get<Vec>();
  if(orr.contains("component_flux"))right_flux=orr.at("component_flux").get<Vec>();
  const auto transport=conservative_species_face_flux_rhs(d.c,flow,Vec(ncells,0.),storage_epsilon_node,Vec{0.},
      reservoirs.first,reservoirs.second,z,epsilon_face,p->advection_scheme,geom.interfaces,widths,
      "concentration",left_flux,right_flux,p->spatial_method);
  std::copy(transport.first.begin(),transport.first.end(),result.begin());
  if(thermal_mode==1)return result;
  const auto left=port_state("left",t),right=port_state("right",t);
  const double tl=ol.value("temperature",left.value("temperature",d.Tg.front())),
               tr=orr.value("temperature",right.value("temperature",d.Tg.back()));
  const Vec faceT=reconstructed_face_values(d.Tg,z,flow,Vec{tl},Vec{tr},p->advection_scheme,geom.interfaces,p->spatial_method);
  Vec enthalpy(Nz+1,0.);
  for(std::size_t f=0;f<=Nz;++f) {
    if(liquid)enthalpy[f]=flow[f]*liquid_density*liquid_heat_capacity*faceT[f];
    else for(std::size_t j=0;j<n_comp;++j)enthalpy[f]+=transport.second[j*(Nz+1)+f]*p->cp[j]*faceT[f];
  }
  if(ol.contains("enthalpy_flux"))enthalpy.front()=ol.at("enthalpy_flux");
  if(orr.contains("enthalpy_flux"))enthalpy.back()=orr.at("enthalpy_flux");
  for(std::size_t i=0;i<Nz;++i) {
    double capacity=0.,mass_energy=0.,ads_capacity=0.;
    if(liquid)capacity=storage_epsilon_node[i]*liquid_density*liquid_heat_capacity;
    else for(std::size_t j=0;j<n_comp;++j) {
      capacity+=storage_epsilon_node[i]*d.c[j*Nz+i]*(p->cp[j]-R);
      mass_energy+=storage_epsilon_node[i]*(p->cp[j]-R)*d.Tg[i]*transport.first[j*Nz+i];
    }
    if(gas_solid_thermal_equilibrium) {
      for(std::size_t m=0;m<geom.materials.size();++m)for(auto group:p->loading_groups[m])
        for(std::size_t j=0;j<n_comp;++j)ads_capacity+=geom.mass_density[m*Nz+i]*state[ncells+group*ncells+j*Nz+i]*p->ads_cp[m*n_comp+j];
      capacity+=geom.heat_capacity[i]+ads_capacity;
    }
    const double value=(-(enthalpy[i+1]-enthalpy[i])/widths[i]-mass_energy)/std::max(capacity,1e-12);
    result[result.size()-3*Nz+i]=value;
    if(gas_solid_thermal_equilibrium)result[result.size()-2*Nz+i]=value;
  }
  return result;
}
double Column::stored_energy(const Vec& state) const
{
  if (state.size() != layout.size()) throw std::invalid_argument("Column state size mismatch");
  double total = 0;
  const auto& geom = p->geometry;
  std::size_t ncells = n_comp * Nz;
  for (std::size_t i = 0; i < Nz; ++i)
  {
    double Tg = state[state.size() - 3 * Nz + i], Ts = state[state.size() - 2 * Nz + i],
           Tw = state[state.size() - Nz + i], density = geom.heat_capacity[i] * Ts;
    for (std::size_t j = 0; j < n_comp; ++j)
      density += liquid ? 0. : storage_epsilon_node[i] * state[j * Nz + i] * (p->cp[j] - R) * Tg;
    if (liquid) density += storage_epsilon_node[i] * liquid_density * liquid_heat_capacity * Tg;
    for (std::size_t m = 0; m < geom.materials.size(); ++m)
      for (auto group : p->loading_groups[m])
        for (std::size_t j = 0; j < n_comp; ++j)
          density += geom.mass_density[m * Nz + i] * state[ncells + group * ncells + j * Nz + i] *
                     (p->ads_ref_u[group * n_comp + j] + p->ads_cp[m * n_comp + j] * (Ts - p->ads_ref_T[group * n_comp + j]));
    if(!liquid)for(std::size_t m=0;m<p->main_models.size();++m)if(p->main_models[m])
      for(std::size_t a=0;a<p->main_models[m]->chemical_sites;++a)for(std::size_t j=0;j<n_comp;++j){
        const auto h=p->history_start[m]+4*a;
        density+=geom.mass_density[m*Nz+i]/p->main_models[m]->particle_density*(p->cp[j]-R)*Ts*
          (state[(1+layout.n_loading_groups+h)*ncells+j*Nz+i]+state[(2+layout.n_loading_groups+h)*ncells+j*Nz+i]);
      }
    total += widths[i] * (area * density + p->wall_area * p->rho_w * p->Cpw * Tw);
    if (layout.transient_htf)
      total += widths[i] * p->jacket_area * p->jacket_rho * p->jacket_cp * state[state.size() - 4 * Nz + i];
  }
  return total;
}
std::pair<Vec, Vec> Column::loading_profiles(const Vec& state) const
{
  Vec physical(n_comp * Nz, 0.), chemical(n_comp * Nz, 0.);
  for (std::size_t i = 0; i < Nz; ++i)
  {
    double mass = 0.;
    for (std::size_t m = 0; m < p->loading_groups.size(); ++m) mass += p->geometry.mass_density[m * Nz + i];
    for (std::size_t m = 0; m < p->loading_groups.size(); ++m)
      for (std::size_t g = 0; g < p->loading_groups[m].size(); ++g)
        for (std::size_t j = 0; j < n_comp; ++j)
        {
          auto& out = g == 0 ? physical : chemical;
          out[j * Nz + i] += p->geometry.mass_density[m * Nz + i] / std::max(mass, 1e-300) *
                             state[(1 + p->loading_groups[m][g]) * n_comp * Nz + j * Nz + i];
        }
  }
  return {physical, chemical};
}
Vec Column::equilibrium_loading(const Vec& state) const
{
  auto d = derived(state);
  Vec drive = d.c, out(n_comp * Nz, 0.);
  if (!liquid)
    for (std::size_t j = 0; j < n_comp; ++j)
      for (std::size_t i = 0; i < Nz; ++i) drive[j * Nz + i] *= R * d.Tg[i];
  for (std::size_t m = 0; m < p->equilibria.size(); ++m)
  {
    auto q = p->equilibria[m].loading(drive, d.Ts);
    for (std::size_t i = 0; i < Nz; ++i)
    {
      double mass = 0.;
      for (std::size_t k = 0; k < p->equilibria.size(); ++k) mass += p->geometry.mass_density[k * Nz + i];
      double weight = p->geometry.mass_density[m * Nz + i] / std::max(mass, 1e-300);
      for (std::size_t j = 0; j < n_comp; ++j) out[j * Nz + i] += weight * p->phy_factors[m] * q[j * Nz + i];
    }
  }
  return out;
}
Json Column::inventories(const Vec& state) const
{
  if (state.size() != layout.size()) throw std::invalid_argument("Column state size mismatch");
  Vec gas(n_comp, 0), ads(n_comp, 0), pore(n_comp,0.), total(n_comp), by_material(p->loading_groups.size() * n_comp, 0);
  std::size_t ncells = n_comp * Nz;
  for (std::size_t j = 0; j < n_comp; ++j)
    for (std::size_t i = 0; i < Nz; ++i) gas[j] += area * widths[i] * storage_epsilon_node[i] * state[j * Nz + i];
  for (std::size_t m = 0; m < p->loading_groups.size(); ++m)
    for (auto group : p->loading_groups[m])
      for (std::size_t j = 0; j < n_comp; ++j)
        for (std::size_t i = 0; i < Nz; ++i)
          by_material[m * n_comp + j] +=
              area * widths[i] * p->geometry.mass_density[m * Nz + i] * state[ncells + group * ncells + j * Nz + i];
  for(std::size_t m=0;m<p->main_models.size();++m)if(p->main_models[m])
    for(std::size_t a=0;a<p->main_models[m]->chemical_sites;++a)for(std::size_t j=0;j<n_comp;++j)for(std::size_t i=0;i<Nz;++i){
      const auto h=p->history_start[m]+4*a;
      pore[j]+=area*widths[i]*p->geometry.mass_density[m*Nz+i]/p->main_models[m]->particle_density*
        (state[(1+layout.n_loading_groups+h)*ncells+j*Nz+i]+state[(2+layout.n_loading_groups+h)*ncells+j*Nz+i]);
    }
  for (std::size_t j = 0; j < n_comp; ++j)
  {
    for (std::size_t m = 0; m < p->loading_groups.size(); ++m) ads[j] += by_material[m * n_comp + j];
    total[j] = gas[j] + ads[j] + pore[j];
  }
  return Json{{"gas", gas},
              {"adsorbed", ads},
              {"surface_and_pore",pore},
              {"total", total},
              {"adsorbed_by_material", matrix_json(by_material, p->loading_groups.size(), n_comp)},
              {"material_mass_kg", p->geometry.material_mass_kg},
              {"energy_J", stored_energy(state)}};
}
Json Column::thermal_utility_rates(double t, const Vec& state, const Json& overrides, const Vec* state_derivative) const
{
  if (thermal_mode == 1)
  {
    Vec own_derivative;
    if (!state_derivative) own_derivative = rhs(t, state, overrides);
    const Vec& deriv = state_derivative ? *state_derivative : own_derivative;
    std::size_t ncells = n_comp * Nz;
    double stored_rate = 0;
    for (std::size_t i = 0; i < Nz; ++i)
    {
      double Tg = state[state.size() - 3 * Nz + i], Ts = state[state.size() - 2 * Nz + i], density = 0;
      for (std::size_t j = 0; j < n_comp; ++j)
        density += liquid ? 0. : storage_epsilon_node[i] * (p->cp[j] - R) * Tg * deriv[j * Nz + i];
      for (std::size_t m = 0; m < p->loading_groups.size(); ++m)
        for (auto group : p->loading_groups[m])
          for (std::size_t j = 0; j < n_comp; ++j)
            density += p->geometry.mass_density[m * Nz + i] *
                       (p->ads_ref_u[group * n_comp + j] + p->ads_cp[m * n_comp + j] * (Ts - p->ads_ref_T[group * n_comp + j])) *
                       deriv[ncells + group * ncells + j * Nz + i];
      if(!liquid)for(std::size_t m=0;m<p->main_models.size();++m)if(p->main_models[m])
        for(std::size_t a=0;a<p->main_models[m]->chemical_sites;++a)for(std::size_t j=0;j<n_comp;++j){
          const auto h=p->history_start[m]+4*a;
          density+=p->geometry.mass_density[m*Nz+i]/p->main_models[m]->particle_density*(p->cp[j]-R)*Ts*
            (deriv[(1+layout.n_loading_groups+h)*ncells+j*Nz+i]+deriv[(2+layout.n_loading_groups+h)*ncells+j*Nz+i]);
        }
      stored_rate += area * widths[i] * density;
    }
    Json b = boundary_enthalpy_rates(t, state, overrides);
    double flow = b["left_in"].get<double>() + b["right_in"].get<double>() - b["left_out"].get<double>() -
                  b["right_out"].get<double>(),
           heat = stored_rate - flow - p->reaction_heat;
    return Json{{"net_heat_W", heat},
                {"isothermal_heat_W", heat},
                {"heating_W", std::max(heat, 0.0)},
                {"cooling_W", std::max(-heat, 0.0)},
                {"electric_heater_W", 0.0},
                {"gas_heating_W", 0.0},
                {"gas_cooling_W", 0.0},
                {"heat_loss_W", 0.0}};
  }
  auto d = derived(state);
  auto external = external_wall_heat_density(d.Tw, t, d.Tj);
  const Json& control = external.second;
  std::string mode = control["mode"];
  double electricity = mode == "electric" ? std::max(control["power_W"].get<double>(), 0.0) : 0, net = 0, positive = 0,
         negative = 0;
  for (std::size_t i = 0; i < Nz; ++i)
  {
    double volume = p->wall_area * widths[i], thermal = external.first[i] - electricity / p->wall_volume;
    net += volume * external.first[i];
    positive += volume * std::max(thermal, 0.0);
    negative += volume * std::max(-thermal, 0.0);
  }
  if (layout.transient_htf && mode == "htf")
  {
    // The system energy inventory includes the jacket. Count boundary enthalpy,
    // not internal wall exchange; this also accounts for heating the HTF holdup.
    const double boundary = control["boundary_heat_W"];
    const double direct = control["power_W"].get<double>() + control["heat_flux_W_m2"].get<double>() * pi * dout * L;
    net = boundary + direct;
    positive = std::max(boundary, 0.) + std::max(direct, 0.);
    negative = std::max(-boundary, 0.) + std::max(-direct, 0.);
  }
  double end_heat = 0.;
  if (p->fixed_wall_ends)
  {
    const Vec kw(Nz, std::max(p->lambda_w, 1e-300));
    const Vec conduction = conductive_divergence(d.Tw, kw, z, widths);
    const double bedh = control["bed_wall_heat_transfer_coefficient"];
    for (std::size_t i : {std::size_t(0), Nz - 1})
    {
      const double bed_wall = 4 * bedh / dint * (d.Tg[i] - d.Tw[i]) * area / p->wall_area;
      const double heat = -p->wall_area * widths[i] * (conduction[i] + bed_wall + external.first[i]);
      end_heat += heat;
      positive += std::max(heat, 0.);
      negative += std::max(-heat, 0.);
    }
    net += end_heat;
  }
  double heating = mode == "ambient" ? 0 : positive + electricity,
         cooling = (mode == "jacket" || mode == "htf" || mode == "heat_flux" || mode == "electric") ? negative : 0,
         loss = mode == "ambient" ? negative : 0,
         refT = p->step.value("gas_utility_reference_temperature",
                              cfg.value("GAS_UTILITY_REFERENCE_TEMPERATURE", cfg.value("T_feed", 298.15))),
         gasheat = 0, gascool = 0;
  Json rates = boundary_component_rates(t, state, overrides);
  for (const std::string side : {"left", "right"})
  {
    Json port = port_state(side, t);
    if (port["mode"] != "closed")
    {
      double condition = dot(rates[side + "_in"].get<Vec>(), p->cp) * (port["temperature"].get<double>() - refT);
      gasheat += std::max(condition, 0.0);
      gascool += std::max(-condition, 0.0);
    }
  }
  return Json{{"net_heat_W", net},        {"heating_W", heating},
              {"cooling_W", cooling},     {"electric_heater_W", electricity},
              {"gas_heating_W", gasheat}, {"gas_cooling_W", gascool},
              {"heat_loss_W", loss}};
}
Json Column::state_record(const Vec& state) const
{
  auto d = derived(state);
  std::size_t ncells = n_comp * Nz;
  Json loads = Json::array(), history = Json::array();
  for (std::size_t group = 0; group < layout.n_loading_groups; ++group)
    loads.push_back(matrix_json(slice(state, ncells + group * ncells, ncells), n_comp, Nz));
  for (std::size_t h = 0; h < layout.n_history_groups; ++h)
    history.push_back(matrix_json(slice(state, (1 + layout.n_loading_groups + h) * ncells, ncells), n_comp, Nz));
  Json record{{"z", z},
              {"concentration", matrix_json(slice(state, 0, ncells), n_comp, Nz)},
              {"loading", matrix_json(d.q, n_comp, Nz)},
              {"loading_branches", loads},
              {"kinetic_history", history},
              {"gas_temperature", slice(state, state.size() - 3 * Nz, Nz)},
              {"solid_temperature", slice(state, state.size() - 2 * Nz, Nz)},
              {"wall_temperature", slice(state, state.size() - Nz, Nz)},
              {"pressure", d.P},
              {"mole_fraction", matrix_json(d.y, n_comp, Nz)},
              {"epsilon", epsilon_node},
              {"gas_storage_void_fraction", storage_epsilon_node},
              {"material_fraction", matrix_json(p->geometry.fractions, p->loading_groups.size(), Nz)},
              {"material_mass_kg", p->geometry.material_mass_kg},
              {"material_names", material_names},
              {"component_names", comp_names},
              {"loading_group_names", loading_group_names},
              {"history_group_names", history_group_names}};
  if (layout.transient_htf) record["jacket_temperature"] = d.Tj;
  return record;
}
Vec Column::state_atol(double default_value) const
{
  Vec out(layout.size(), default_value);
  std::fill(out.begin(), out.begin() + n_comp * Nz, cfg.value("ATOL_CONCENTRATION", default_value));
  const auto thermal_blocks = 3 + layout.transient_htf;
  std::fill(out.begin() + n_comp * Nz, out.end() - thermal_blocks * Nz, cfg.value("ATOL_LOADING", default_value));
  std::fill(out.end() - thermal_blocks * Nz, out.end(), cfg.value("ATOL_TEMPERATURE", std::max(default_value, 1e-6)));
  return out;
}
double Column::metric(const std::string& metric_name, const Vec& state, double t, const Json& options) const
{
  auto d = derived(state);
  std::string key = lower(metric_name), side = options.value("port", options.value("side", std::string{}));
  int component = -1;
  if (options.contains("component") && !options["component"].is_null())
  {
    if (options["component"].is_string())
    {
      auto it = std::find(comp_names.begin(), comp_names.end(), options["component"].get<std::string>());
      if (it == comp_names.end()) throw std::invalid_argument("Unknown metric component");
      component = it - comp_names.begin();
    }
    else
      component = options["component"].get<int>();
    if (component < 0 || static_cast<std::size_t>(component) >= n_comp)
      throw std::invalid_argument("Invalid metric component");
  }
  auto average = [&](const Vec& a) { return dot(a, widths) / L; };
  auto minimum = [](const Vec& a) { return *std::min_element(a.begin(), a.end()); };
  auto maximum = [](const Vec& a) { return *std::max_element(a.begin(), a.end()); };
  auto require_component = [&]()
  {
    if (component < 0) throw std::invalid_argument("Metric requires component");
  };
  auto require_port = [&]()
  {
    if (side != "left" && side != "right") throw std::invalid_argument("Metric requires left or right port");
  };
  if (key == "loading_equilibrium_error")
  {
    require_component();
    Vec driving = d.c;
    if (!liquid)
      for (std::size_t j = 0; j < n_comp; ++j)
        for (std::size_t i = 0; i < Nz; ++i) driving[j * Nz + i] *= R * d.Tg[i];
    double error = 0.;
    for (std::size_t m = 0; m < p->equilibria.size(); ++m)
    {
      auto target = p->equilibria[m].loading(driving, d.Ts);
      for (std::size_t i = 0; i < Nz; ++i)
        if (p->geometry.fractions[m * Nz + i] > 0.)
        {
          double q = state[(1 + p->loading_groups[m][0]) * n_comp * Nz + component * Nz + i],
                 qe = target[component * Nz + i] * p->phy_factors[m], scale = std::max(std::abs(qe), 1e-9);
          error = std::max(error, std::abs(q - qe) / scale);
          if (p->meshes[m])
            for (std::size_t k = 0; k < p->meshes[m]->weights.size(); ++k)
              error = std::max(
                  error,
                  std::abs(
                      state[(1 + layout.n_loading_groups + p->radial_start[m] + k) * n_comp * Nz + component * Nz + i] -
                      qe) /
                      scale);
        }
    }
    return error;
  }
  if (key == "bed_average_pressure") return average(d.P);
  if (key == "bed_minimum_pressure") return minimum(d.P);
  if (key == "bed_maximum_pressure") return maximum(d.P);
  if (key == "left_pressure" || key == "left_end_pressure") return d.P.front();
  if (key == "right_pressure" || key == "right_end_pressure") return d.P.back();
  if (key == "pressure_drop") return d.P.front() - d.P.back();
  if (key == "absolute_pressure_drop") return std::abs(d.P.front() - d.P.back());
  if (key == "pressure_nonuniformity") return (maximum(d.P) - minimum(d.P)) / std::max(average(d.P), 1e-30);
  if (key == "bed_average_temperature" || key == "bed_average_gas_temperature") return average(d.Tg);
  if (key == "bed_average_solid_temperature" || key == "average_solid_temperature") return average(d.Ts);
  if (key == "bed_average_wall_temperature" || key == "average_wall_temperature") return average(d.Tw);
  if (key == "bed_minimum_temperature") return std::min({minimum(d.Tg), minimum(d.Ts), minimum(d.Tw)});
  if (key == "bed_maximum_temperature") return std::max({maximum(d.Tg), maximum(d.Ts), maximum(d.Tw)});
  if (key == "bed_minimum_solid_temperature") return minimum(d.Ts);
  if (key == "bed_maximum_solid_temperature") return maximum(d.Ts);
  if (key == "temperature_nonuniformity")
    return (maximum(d.Tg) - minimum(d.Tg)) / std::max(std::abs(average(d.Tg)), 1e-30);
  if (key == "solid_temperature_nonuniformity")
    return (maximum(d.Ts) - minimum(d.Ts)) / std::max(std::abs(average(d.Ts)), 1e-30);
  if (key == "wall_temperature_nonuniformity")
    return (maximum(d.Tw) - minimum(d.Tw)) / std::max(std::abs(average(d.Tw)), 1e-30);
  if (key == "wall_temperature_approach" || key == "maximum_wall_temperature_approach")
  {
    auto external = external_wall_heat_density(d.Tw, t);
    Vec reference = profile(external.second.value("temperature_profile", external.second["temperature"]), Nz,
                            "utility temperature");
    double result = 0;
    for (std::size_t i = 0; i < Nz; ++i) result = std::max(result, std::abs(d.Tw[i] - reference[i]));
    return result;
  }
  if (key == "bed_average_loading" || key == "average_loading")
  {
    require_component();
    return inventories(state)["adsorbed"][component].get<double>() /
           std::max(sum(p->geometry.material_mass_kg), 1e-300);
  }
  if (key == "bed_maximum_loading")
  {
    require_component();
    return maximum(slice(d.q, component * Nz, Nz));
  }
  if (key == "adsorbed_inventory" || key == "gas_inventory" || key == "total_inventory")
  {
    require_component();
    return inventories(state)[key.substr(0, key.find('_'))][component].get<double>();
  }
  if (key == "port_mole_fraction" || key == "outlet_mole_fraction" || key == "port_concentration" ||
      key == "port_concentration_ratio" || key == "outlet_concentration_ratio")
  {
    require_component();
    if ((key == "outlet_mole_fraction" || key == "outlet_concentration_ratio") && side.empty())
    {
      Vec vf = face_velocities(state, t);
      double left = std::max(-vf.front(), 0.0), right = std::max(vf.back(), 0.0);
      if (left == 0 && right == 0) return 0;
      side = left > right ? "left" : "right";
    }
    require_port();
    std::size_t i = side == "left" ? 0 : Nz - 1;
    if (key == "port_mole_fraction" || key == "outlet_mole_fraction") return d.y[component * Nz + i];
    if (key == "port_concentration") return d.c[component * Nz + i];
    double reference;
    if (options.contains("reference_concentration"))
      reference = options["reference_concentration"];
    else
    {
      if (lower(options.value("reference", std::string("feed"))) != "feed")
        throw std::invalid_argument("Concentration ratio requires feed reference or explicit concentration");
      Vec yfeed = normal(cfg.at("y_feed"), n_comp, "y_feed");
      double P = options.value("reference_pressure",
                               cfg.value("P_FEED", cfg.value("P_HIGH", std::max(d.P.front(), d.P.back())))),
             T = options.value("reference_temperature", cfg.value("T_feed", cfg.value("T_init", d.Tg.front())));
      reference = liquid ? arr(cfg.at("LIQUID_FEED_CONCENTRATION"), n_comp, "feed concentration")[component]
                         : yfeed[component] * P / (R * T);
    }
    return d.c[component * Nz + i] / std::max(reference, 1e-300);
  }
  if (key == "port_velocity" || key == "absolute_port_velocity" || key == "port_superficial_velocity" ||
      key == "absolute_port_superficial_velocity")
  {
    require_port();
    Vec vf = face_velocities(state, t);
    double value = side == "left" ? vf.front() : vf.back();
    if (key.find("superficial") != std::string::npos)
      value *= side == "left" ? epsilon_node.front() : epsilon_node.back();
    return key.rfind("absolute", 0) == 0 ? std::abs(value) : value;
  }
  if (key == "port_total_molar_flow")
  {
    require_port();
    Json rates = boundary_component_rates(t, state);
    return sum(rates[side + "_in"].get<Vec>()) - sum(rates[side + "_out"].get<Vec>());
  }
  throw std::invalid_argument("Unsupported stop metric " + metric_name);
}
}  // namespace cycle

namespace cycle {
bool Column::has_reactions() const {for(const auto& m:p->main_models)if(m&&!m->reactions.empty())return true;return false;}
Vec Column::reaction_component_rates() const {return p->reaction_rates.empty()?Vec(n_comp,0.):p->reaction_rates;}
double Column::reaction_energy_rate() const {return p->reaction_heat;}
}

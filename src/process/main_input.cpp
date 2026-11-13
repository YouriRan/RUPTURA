#include "cycle/main_input.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>

#include "cycle/runner.hpp"
#include "cycle/kinetics.hpp"
namespace cycle
{
namespace
{
std::string lower(std::string s)
{
  for (char& c : s) c = std::tolower(static_cast<unsigned char>(c));
  return s;
}
std::string snake(const std::string& s)
{
  std::string out;
  for (std::size_t i = 0; i < s.size(); ++i)
  {
    auto c = static_cast<unsigned char>(s[i]);
    if (std::isupper(c) && i && std::islower(static_cast<unsigned char>(s[i - 1]))) out += '_';
    out += std::tolower(c);
  }
  return out;
}
void require(bool v, const std::string& message)
{
  if (!v) throw std::invalid_argument(message);
}
Json get(const Json& j, const std::string& key, Json def = nullptr)
{
  const Json* found = nullptr;
  for (auto& [k, v] : j.items())
    if (lower(k) == lower(key))
    {
      require(!found, "Duplicate case-insensitive key: " + key);
      found = &v;
    }
  return found ? *found : def;
}
Json component_values(const Json& raw, const std::vector<std::string>& names, Json def = 0.)
{
  if (raw.is_null()) return Vec(names.size(), def.get<double>());
  if (!raw.is_object()) return raw;
  Json out = Json::array();
  for (auto& n : names) out.push_back(raw.value(n, def));
  for (auto& [k, v] : raw.items())
    require(std::find(names.begin(), names.end(), k) != names.end(), "Unknown component: " + k);
  return out;
}
Json normalized_fractions(const Json& raw, const std::vector<std::string>& names)
{
  auto values = component_values(raw, names);
  if (!values.empty() && values[0].is_array()) return values;
  double total = 0.;
  for (auto& value : values)
  {
    double x = value.get<double>();
    require(std::isfinite(x) && x >= 0., "Mole fractions must be finite and nonnegative");
    total += x;
  }
  require(total > 0., "At least one mole fraction must be positive");
  for (auto& value : values) value = value.get<double>() / total;
  return values;
}
void known_options(const Json& object, const std::set<std::string>& names, const std::string& context)
{
  require(object.is_object(), context + " must be an object");
  for (auto& [key, value] : object.items()) require(names.count(lower(key)), "Unknown " + context + " setting: " + key);
}
double physical_value(const Json& component, const std::string& key, double fallback)
{
  double result = get(component, key, fallback);
  std::optional<double> override;
  for (auto& site : get(component, "PhysisorptionSites", Json::array()))
  {
    auto value = key == "MassTransferCoefficient" ? get(get(site, "RateEquationParameters", Json::object()), key)
                                                  : get(site, key);
    if (!value.is_null())
    {
      double x = value.get<double>();
      require(!override || x == *override, "Conflicting " + key + " across PhysisorptionSites");
      override = x;
    }
    require(lower(get(site, "RateEquation", "FirstOrder").get<std::string>()) == "firstorder",
            "PhysisorptionSites support FirstOrder kinetics; use ParticleTransport for HSDM");
  }
  return override.value_or(result);
}
Json first_number(const Json& value, Json fallback)
{
  if (value.is_array()) return value.empty() ? fallback : value.front();
  return value.is_null() ? fallback : value;
}
const std::map<std::string, std::string> aliases = {{"controllerconductance", "conductance_mol_s_pa"},
                                                    {"flowcontrol", "flow_control"},
                                                    {"pressuregradient", "pressure_gradient"},
                                                    {"deadband", "deadband_pa"},
                                                    {"fallbackcomposition", "fallback_composition"},
                                                    {"failuremaximumsimulatedtime", "failure_maximum_simulated_time"},
                                                    {"timeoffset", "time_offset_s"},
                                                    {"decayrate", "decay_rate_s_inv"},
                                                    {"pressurecontrol", "pressure_control"},
                                                    {"setpointpressure", "setpoint_pressure"},
                                                    {"exporttag", "export_tag"},
                                                    {"importtag", "import_tag"},
                                                    {"maximumflow", "maximum_flow"},
                                                    {"supplysource", "supply_source"},
                                                    {"autofill", "autofill"},
                                                    {"initialmoles", "initial_moles"},
                                                    {"thermal", "thermal"},
                                                    {"massflow", "mass_flow_rate_kg_s"},
                                                    {"outletpressure", "outlet_pressure"},
                                                    {"equilibrium", "equilibrium"},
                                                    {"kads", "k_ads"},
                                                    {"kdes", "k_des"},
                                                    {"eaads", "Ea_ads"},
                                                    {"eades", "Ea_des"},
                                                    {"heatofadsorption", "heat_release_J_mol"},
                                                    {"heataffinity", "heat_affinity_J_mol"},
                                                    {"qsat", "qsat"},
                                                    {"affinity", "b"},
                                                    {"henrycoefficient", "K"},
                                                    {"fractions", "fractions"},
                                                    {"capacityfraction", "capacity_fraction"},
                                                    {"independentsites", "independent_sites"},
                                                    {"orderads", "order_ads"},
                                                    {"orderdes", "order_des"},
                                                    {"order", "order"},
                                                    {"ratebasis", "rate_basis"},
                                                    {"qreference", "q_reference"},
                                                    {"drivingregularizer", "driving_regularizer"},
                                                    {"nads", "n_ads"},
                                                    {"ndes", "n_des"},
                                                    {"n", "n"},
                                                    {"ageregularizer", "age_regularizer"},
                                                    {"historypolicy", "history_policy"},
                                                    {"memoryerasurerateads", "memory_erasure_rate_ads"},
                                                    {"memoryerasureratedes", "memory_erasure_rate_des"},
                                                    {"inletpressure", "inlet_pressure"},
                                                    {"name", "name"},
                                                    {"volumetricconductance", "volume_conductance"},
                                                    {"deltapressure", "pressure_difference"},
                                                    {"jackettemperature", "jacket_temperature"},
                                                    {"walltemperature", "wall_temperature"},
                                                    {"solidtemperature", "solid_temperature"},
                                                    {"initialloading", "loading"},
                                                    {"criteria", "criteria"},
                                                    {"logic", "logic"},
                                                    {"mode", "mode"},
                                                    {"source", "source"},
                                                    {"sink", "sink"},
                                                    {"buffer", "buffer"},
                                                    {"column", "column"},
                                                    {"side", "side"},
                                                    {"leftport", "left_port"},
                                                    {"rightport", "right_port"},
                                                    {"columns", "columns"},
                                                    {"connections", "connections"},
                                                    {"from", "from"},
                                                    {"to", "to"},
                                                    {"pressure", "pressure"},
                                                    {"temperature", "temperature"},
                                                    {"velocity", "velocity"},
                                                    {"columnentrancevelocity", "velocity"},
                                                    {"gasphasemolfractions", "composition"},
                                                    {"liquidphaseconcentrations", "concentration"},
                                                    {"composition", "composition"},
                                                    {"concentration", "concentration"},
                                                    {"molarflowrate", "molar_flow_mol_s"},
                                                    {"volumetricflowrate", "volume_flow_m3_s"},
                                                    {"streamtags", "stream_tags"},
                                                    {"leftin", "left_in"},
                                                    {"leftout", "left_out"},
                                                    {"rightin", "right_in"},
                                                    {"rightout", "right_out"},
                                                    {"duration", "duration"},
                                                    {"minimumduration", "minimum_duration"},
                                                    {"maximumduration", "maximum_duration"},
                                                    {"stopcriterion", "stop"},
                                                    {"timeoutpolicy", "timeout_policy"},
                                                    {"metric", "metric"},
                                                    {"operator", "operator"},
                                                    {"value", "value"},
                                                    {"all", "all"},
                                                    {"any", "any"},
                                                    {"logic", "logic"},
                                                    {"criteria", "criteria"},
                                                    {"component", "component"},
                                                    {"port", "port"},
                                                    {"model", "model"},
                                                    {"initialvalue", "initial_value"},
                                                    {"finalvalue", "final_value"},
                                                    {"ramptime", "ramp_time_s"},
                                                    {"times", "times_s"},
                                                    {"values", "values"},
                                                    {"referenceconcentration", "reference_concentration"},
                                                    {"referencetemperature", "reference_temperature"},
                                                    {"referencepressure", "reference_pressure"},
                                                    {"compositionsource", "composition_source"},
                                                    {"thermaltcontrol", "thermal"},
                                                    {"thermalcontrol", "thermal"},
                                                    {"heattransfercoefficient", "heat_transfer_coefficient"},
                                                    {"massflowrate", "mass_flow_rate_kg_s"},
                                                    {"heatcapacity", "heat_capacity_J_kg_K"},
                                                    {"flowdirection", "flow_direction"},
                                                    {"power", "power_W"},
                                                    {"heatflux", "heat_flux_W_m2"},
                                                    {"conductance", "conductance"},
                                                    {"flowmodel", "flow_model"},
                                                    {"opening", "opening"},
                                                    {"checkvalve", "check_valve"},
                                                    {"maximummolarflow", "maximum_molar_flow"},
                                                    {"area", "area_m2"},
                                                    {"dischargecoefficient", "Cd"},
                                                    {"heatcapacityratio", "gamma"},
                                                    {"tag", "tag"},
                                                    {"intag", "in_tag"},
                                                    {"outtag", "out_tag"},
                                                    {"compressor", "compressor"},
                                                    {"efficiency", "efficiency"},
                                                    {"dischargepressure", "discharge_pressure"},
                                                    {"flowrate", "flow_mol_s"},
                                                    {"volume", "volume"},
                                                    {"initialpressure", "initial_pressure"},
                                                    {"initialtemperature", "initial_temperature"},
                                                    {"initialcomposition", "initial_composition"},
                                                    {"initialconcentration", "initial_concentration"},
                                                    {"thermalmode", "thermal_mode"},
                                                    {"ambienttemperature", "ambient_temperature"},
                                                    {"ua", "UA"},
                                                    {"initialvolume", "initial_volume"},
                                                    {"maximumvolume", "maximum_volume"},
                                                    {"density", "density"},
                                                    {"liquidheatcapacity", "liquid_heat_capacity"},
                                                    {"waitaction", "wait_action"}};
Json translated(const Json& x, const std::vector<std::string>& names, bool named = false)
{
  if (x.is_array())
  {
    Json r = Json::array();
    for (auto& v : x) r.push_back(translated(v, names));
    return r;
  }
  if (!x.is_object()) return x;
  Json r = Json::object();
  for (auto& [k, v] : x.items())
  {
    auto key = lower(k);
    auto it = aliases.find(key);
    std::string target = it == aliases.end() ? k : it->second;
    if (named)
      target = k;
    else if (it == aliases.end())
    {
      bool internal = false;
      for (auto& [a, b] : aliases) internal |= b == k;
      require(internal, "Unknown process setting: " + k);
    }
    require(!r.contains(target), "Duplicate process setting: " + target);
    if (target == "composition" || target == "initial_composition" || target == "fallback_composition")
      r[target] = normalized_fractions(v, names);
    else if (target == "concentration" || target == "initial_concentration")
      r[target] = component_values(v, names);
    else if (target == "columns")
      r[target] = translated(v, names, true);
    else
      r[target] = translated(v, names);
    if ((target == "flow_model" || target == "mode" || target == "model" || target == "metric" ||
         target == "flow_direction" || target == "thermal_mode") &&
        r[target].is_string())
      r[target] = snake(r[target].get<std::string>());
  }
  if (r.contains("all") || r.contains("any"))
  {
    require(r.size() == 1, "Use All or Any as the sole logical stop key");
    const auto key = r.contains("all") ? "all" : "any";
    return Json{{"logic", key}, {"criteria", r[key]}};
  }
  return r;
}
Json equilibrium(const Json& doc, const Json& comps, const std::vector<std::string>& names)
{
  Json e = {{"method", get(doc, "MixturePredictionMethod", "IAST")}, {"components", Json::object()}};
  const std::set<std::string> native={"langmuir","langmuir_freundlich","langmuir-freundlich","sips","toth","henry","freundlich","anti_langmuir","anti-langmuir","modified_redlich","ph-langmuir","gab",
    "humidity_enhanced_langmuir","moisture_swing_langmuir","gab_coadsorption_langmuir","gab_coadsorption_toth","gab_coadsorption_sips","gab_modified_toth_liu","gab_water","humidity_toth_climeworks"};
  const auto mixture_method = lower(e["method"].get<std::string>());
  bool main_adapter=lower(get(doc,"pHMode","Fixed").get<std::string>())!="fixed"||mixture_method=="ei"||mixture_method=="mpd";
  for(const auto& c:comps){std::string first;for(const auto& site:get(c,"PhysisorptionSites",Json::array())){
    const auto model=snake(get(site,"Type","Langmuir").get<std::string>());
    main_adapter|=!native.contains(model)||(!first.empty()&&first!=model);if(first.empty())first=model;
  }}
  if(main_adapter){
    e["method"]="SPI";
    for(const auto& name:names)e["components"][name]={{"model","zero"},{"sites",Json::array({Json{{"capacity_mol_kg",0.},{"k0",0.}}})}};
    e["_main_adapter"]=true;return e;
  }
  for (std::size_t i = 0; i < comps.size(); ++i)
  {
    const auto& c = comps[i];
    Json sites = Json::array();
    std::string common;
    for (auto& site : get(c, "PhysisorptionSites", Json::array()))
    {
      std::string model = snake(get(site, "Type", "Langmuir").get<std::string>());
      if (lower(get(site, "Type", "Langmuir").get<std::string>()) == "ph-langmuir") model = "ph_langmuir";
      if (model == "langmuir-freundlich") model = "langmuir_freundlich";
      Json params = get(site, "Parameters", Json::array());
      require(params.is_array() && !params.empty(), "PhysisorptionSites.Parameters must be an array");
      double q = params[0], k = params.size() > 1 ? params[1].get<double>() : 1.;
      double heat =
          get(site, "EquilibriumHeat",
              get(c, "nonIsothermal", false).get<bool>() ? Json(physical_value(c, "HeatOfAdsorption", 0.)) : Json(0.))
              .get<double>();
      Json ref = get(site, "referenceTemperature", get(c, "referenceTemperature"));
      double k0 = k * (ref.is_number() ? std::exp(-heat / (8.31446261815324 * ref.get<double>())) : 1.);
      if (model == "henry")
      {
        q = 1.;
        k0 =
            params[0].get<double>() * (ref.is_number() ? std::exp(-heat / (8.31446261815324 * ref.get<double>())) : 1.);
      }
      if (model == "freundlich")
      {
        q = params[0];
        k0 = 1.;
        heat = 0.;
      }
      Json s = {{"capacity_mol_kg", q}, {"k0", k0}, {"heat_J_mol", heat}};
      if (params.size() > 2) s["nu0"] = model == "sips" ? 1. / params[2].get<double>() : params[2].get<double>();
      if (model == "freundlich") s["nu0"] = params[1];
      if (!get(site, "CapacityTemperature").is_null())
      {
        auto law = get(site, "CapacityTemperature");
        s["capacity_temperature"] = {{"model", lower(get(law, "Model", "exponential").get<std::string>())},
                                     {"reference_temperature_K", get(law, "referenceTemperature", 298.15)},
                                     {"coefficient", get(law, "Coefficient", 0.)}};
      }
      if (model == "ph_langmuir")
      {
        require(params.size() == 3, "pH-Langmuir uses [capacity, affinity, pH0]");
        s["pH0"] = params[2];
      }
      if (model == "gab")
      {
        require(params.size() == 5, "GAB uses [capacity, C0, K0, dHC, dHK]");
        model = "gab_main";
        s["C0"] = params[1];
        s["K0"] = params[2];
        s["dHC"] = params[3];
        s["dHK"] = params[4];
      }
      Json extras = get(site, "ExtraParameters", Json::object());
      for (auto& [n, v] : extras.items()) s[n] = v;
      if (common.empty()) common = model;
      require(common == model, "FVM multisite parameters currently require the same model within each component");
      sites.push_back(s);
    }
    if (sites.empty())
    {
      common = "zero";
      sites = Json::array({Json{{"capacity_mol_kg", 0.}, {"k0", 0.}}});
    }
    e["components"][names[i]] = {{"model", common}, {"sites", sites}};
  }
  std::size_t maxsites = 1;
  for (auto& [name, c] : e["components"].items()) maxsites = std::max(maxsites, c["sites"].size());
  for (auto& [name, c] : e["components"].items())
    while (c["sites"].size() < maxsites)
      c["sites"].push_back({{"capacity_mol_kg", 0.}, {"k0", 0.}, {"heat_J_mol", 0.}});
  if (lower(e["method"].get<std::string>()) == "sci")
  {
    std::string model = "langmuir";
    std::size_t count = 1;
    for (auto& [name, c] : e["components"].items())
      if (c["model"] != "zero")
      {
        model = c["model"];
        count = c["sites"].size();
        break;
      }
    for (auto& [name, c] : e["components"].items())
      if (c["model"] == "zero")
      {
        c["model"] = model;
        c["sites"] = Json::array();
        for (std::size_t k = 0; k < count; ++k) c["sites"].push_back({{"capacity_mol_kg", 0.}, {"k0", 0.}});
      }
  }
  return e;
}
Json column_seed(const Json& d, const Json& components, const std::vector<std::string>& names,
                 const std::filesystem::path& base)
{
  Json geometry = get(d, "Geometry", Json::object());
  std::string kind = lower(get(geometry, "Type", "PackedBed").get<std::string>());
  require(kind == "monolith" || kind == "packedbed" || kind == "fixedbed",
          "Geometry.Type must be PackedBed or Monolith");
  Json b = {{"type", kind == "monolith" ? "monolith" : "fixed_bed"},
            {"length_m", get(d, "ColumnLength", .3)},
            {"nodes", get(d, "NumberOfGridPoints", 40).get<int>() + 1},
            {"inner_diameter_m", get(geometry, "InternalDiameter", .1)},
            {"outer_diameter_m", get(geometry, "OuterDiameter", get(geometry, "InternalDiameter", .1))},
            {"void_fraction", get(geometry, "ColumnVoidFraction", get(d, "ColumnVoidFraction", .4))},
            {"adsorbent_density_kg_m3", get(d, "ParticleDensity", 1000.)},
            {"particle_diameter_m", get(geometry, "ParticleDiameter", get(d, "ParticleDiameter", .001))},
            {"solid_heat_capacity_J_kg_K", get(d, "heatCapacitySolid", 900.)},
            {"gas_solid_heat_transfer_W_m2_K", get(d, "heatTransferGasSolid", get(d, "heatTransferLiquidSolid", 0.))}};
  Json rates = Json::array(), heats = Json::array(), dispersion = Json::array(), ads_cp = Json::array();
  for (auto& c : components)
  {
    rates.push_back(physical_value(c, "MassTransferCoefficient", 0.));
    heats.push_back(physical_value(c, "HeatOfAdsorption", 0.));
    dispersion.push_back(get(c, "AxialDispersionCoefficient", 0.));
    ads_cp.push_back(get(c, "AdsorbedHeatCapacity", 0.));
  }
  bool energy = get(d, "energyBalance", false);
  std::string thermal_mode = energy ? "wall_exchange" : "isothermal";
  if (!get(d, "ThermalModel").is_null())
  {
    auto model = lower(get(d, "ThermalModel").get<std::string>());
    if (model == "isothermal") thermal_mode = "isothermal";
    else if (model == "adiabatic") thermal_mode = "adiabatic";
    else if (model == "nonadiabatic" || model == "non_adiabatic") thermal_mode = "wall_exchange";
    else require(false, "ThermalModel must be Isothermal, Adiabatic or NonAdiabatic");
    if (!get(d, "energyBalance").is_null())
      require(energy == (thermal_mode != "isothermal"), "energyBalance conflicts with ThermalModel");
  }
  Json t = {{"mode", thermal_mode},
            {"wall_density_kg_m3", get(d, "wallDensity", 7800.)},
            {"wall_heat_capacity_J_kg_K", get(d, "heatCapacityWall", 450.)},
            {"wall_conductivity_W_m_K", get(d, "wallThermalConductivity", 0.)},
            {"gas_wall_heat_transfer_W_m2_K", get(d, "heatTransferGasWall", get(d, "heatTransferLiquidWall", 0.))},
            {"ambient_heat_transfer_W_m2_K", get(d, "heatTransferWallExternal", 0.)},
            {"ambient_temperature_K", get(d, "AmbientTemperature", get(d, "Temperature", 298.15))}};
  Json out = {{"bed", b},
              {"thermal", t},
              {"equilibrium", equilibrium(d, components, names)},
              {"kinetics", {{"ldf_rate_s_inv", rates}, {"adsorption_heat_J_mol", heats}}},
              {"advanced", {{"ADSORBED_CP_MOLAR", ads_cp}, {"Dj", dispersion}}}};
  Json original={{"Components",components}};
  for(const auto key:{"Geometry","Temperature","ParticleDensity","MixturePredictionMethod","Reactions","FluidPhase","pHMode","pHValue","pKw","pHComponent","LiquidDensity"})
    if(!get(d,key).is_null())original[key]=get(d,key);
  if (lower(get(d,"MixturePredictionMethod","IAST").get<std::string>()) == "mpd") {
    require(lower(get(d,"FluidPhase","Gas").get<std::string>()) == "gas",
            "MPD mixture prediction is available only for gas simulations");
    auto settings = get(d,"MPDSettings");
    require(settings.is_object(), "MixturePredictionMethod MPD requires MPDSettings");
    const auto file = get(settings,"FileName");
    require(file.is_string() && !file.get<std::string>().empty(), "MPDSettings.FileName is required");
    auto path = std::filesystem::path(file.get<std::string>());
    if (path.is_relative()) path = base / path;
    path = std::filesystem::absolute(path).lexically_normal();
    for (auto it=settings.begin();it!=settings.end();)
      if (lower(it.key())=="filename") it=settings.erase(it); else ++it;
    settings["FileName"] = path.string();
    original["MPDSettings"] = settings;
    // External distribution bytes are part of the physical model. Include a
    // content digest in restart compatibility, not only the file's path.
    std::ifstream distribution(path,std::ios::binary);
    require(static_cast<bool>(distribution), "MPD distribution file could not be opened: " + path.string());
    std::uint64_t hash=14695981039346656037ull;
    std::array<char,65536> buffer{};
    while (distribution.read(buffer.data(),buffer.size()) || distribution.gcount())
      for (std::streamsize i=0;i<distribution.gcount();++i) {
        hash^=static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]);
        hash*=1099511628211ull;
      }
    require(distribution.eof(), "MPD distribution file could not be read: " + path.string());
    std::ostringstream digest;digest<<"fnv1a64:"<<std::hex<<hash;
    out["advanced"]["MPD_SOURCE_DIGEST"] = digest.str();
  }
  bool main_kinetics=!get(d,"Reactions",Json::array()).empty();
  for(const auto& component:components)main_kinetics|=!get(component,"ChemisorptionSites",Json::array()).empty();
  if(main_kinetics)out["advanced"]["MAIN_RAW_INPUT"]=original;
  if(out["equilibrium"].value("_main_adapter",false)){
    out["equilibrium"].erase("_main_adapter");out["advanced"]["MAIN_EQUILIBRIUM"]=original;
  }
  out["advanced"]["PHYSISORPTION_FRACTION"] = get(d,"PhysisorptionFraction",1.0);
  out["advanced"]["PARTICLE_SPHERICITY"] =
      get(geometry, "ParticleSphericity", get(d, "ParticleSphericity", 1.0));
  Json particle = get(d, "ParticleTransport");
  if (!particle.is_null())
  {
    Json p = Json::object();
    const std::map<std::string, std::string> pm = {{"Model", "model"},
                                                   {"Geometry", "geometry"},
                                                   {"Radius", "radius"},
                                                   {"RadialCells", "radial_cells"},
                                                   {"SurfaceDiffusivity", "surface_diffusivity"},
                                                   {"DiffusivityActivationEnergy", "diffusivity_activation_energy"},
                                                   {"referenceTemperature", "reference_temperature"},
                                                   {"FilmResistance", "film_resistance"},
                                                   {"FilmCoefficient", "film_coefficient"},
                                                   {"LDFCoefficientSource", "ldf_coefficient_source"},
                                                   {"RateActivationEnergy", "rate_activation_energy"}};
    for (auto& [k, v] : particle.items())
    {
      auto it = std::find_if(pm.begin(), pm.end(), [&](auto& kv) { return lower(kv.first) == lower(k); });
      require(it != pm.end(), "Unknown ParticleTransport option: " + k);
      p[it->second] = v;
      if (v.is_object()) p[it->second] = component_values(v, names);
    }
    for (const auto key : {"model", "geometry", "ldf_coefficient_source"})
      if (p.contains(key)) p[key] = snake(p[key].get<std::string>());
    out["advanced"]["PARTICLE_TRANSPORT"] = p;
  }
  if (kind == "monolith")
  {
    Json m = {{"shape", get(geometry, "ChannelShape", "square")},
              {"channel_size_m", get(geometry, "ChannelSize", get(geometry, "ChannelDiameter", .001))},
              {"coating_thickness_m", get(geometry, "WashcoatThickness", .00005)},
              {"coating_solid_fraction", get(geometry, "WashcoatSolidFraction", 1.)},
              {"substrate_density_kg_m3", get(geometry, "SubstrateDensity", 2500.)},
              {"substrate_heat_capacity_J_kg_K", get(geometry, "SubstrateHeatCapacity", 800.)}};
    if (!get(geometry, "NumberOfChannels").is_null())
      m["channel_count"] = get(geometry, "NumberOfChannels");
    else
      m["open_fraction"] = get(geometry, "ColumnVoidFraction", .6);
    out["monolith"] = m;
  }
  return out;
}
Json phase_document(const Json& raw, const Json& d, const std::vector<std::string>& names)
{
  Json phase = translated(raw, names);
  if (phase.contains("columns"))
  {
    Json coolant = get(d, "Coolant");
    if (!coolant.is_null() && lower(get(d, "ThermalModel", "").get<std::string>()) != "adiabatic")
      for (auto& action : phase["columns"])
        if (!action.contains("thermal"))
          action["thermal"] = {
              {"mode", "htf"},
              {"temperature", get(coolant, "InletTemperature", 298.15)},
              {"mass_flow_rate_kg_s", get(coolant, "MassFlowRate", 0.)},
              {"heat_capacity_J_kg_K", get(coolant, "HeatCapacity", 4180.)},
              {"heat_transfer_coefficient", get(coolant, "HeatTransferCoefficient", 0.)},
              {"flow_direction", snake(get(coolant, "FlowDirection", "CoCurrent").get<std::string>())}};
    return phase;
  }
  Json action = {{"left_port",
                  {{"mode", "velocity"},
                   {"source", "feed"},
                   {"velocity", get(raw, "ColumnEntranceVelocity", get(d, "ColumnEntranceVelocity", .01))}}},
                 {"right_port", {{"mode", "pressure"}, {"source", "outlet"}}},
                 {"stream_tags", {{"left_in", "feed"}, {"right_out", "light_product"}, {"left_out", "heavy_product"}}}};
  std::string bc = lower(get(raw, "BoundaryCondition", get(d, "BoundaryCondition", "InletVelocityOutletPressure")).get<std::string>());
  if (bc == "constantvelocity") bc = "fixedvelocity";
  if (bc == "constantpressuregradient") bc = "fixedpressuregradient";
  if (bc == "inletpressureoutletpressure" || bc == "fixedpressuregradient")
    action["left_port"] = {{"mode", "pressure"}, {"source", "feed"}};
  else if (bc != "inletvelocityoutletpressure" && bc != "fixedvelocity")
    require(false,
            "For FVM, specify pressure/velocity/closed ports or use "
            "InletVelocityOutletPressure/InletPressureOutletPressure/FixedVelocity/FixedPressureGradient");
  if (bc == "fixedvelocity")
  {
    action["flow_control"] = {{"model", "fixed_velocity"}};
    action["right_port"] = {{"mode", "velocity"}, {"source", "outlet"}, {"velocity", 0.}};
  }
  if (bc == "fixedpressuregradient")
  {
    require(lower(get(d, "FluidPhase", "Gas").get<std::string>()) == "liquid",
            "FVM FixedPressureGradient requires FluidPhase Liquid; gas uses pressure ports or the original main solver");
    auto gradient = get(raw, "PressureGradient", get(d, "PressureGradient"));
    require(!gradient.is_null(), "FixedPressureGradient requires PressureGradient in Pa/m");
    action["flow_control"] = {{"model", "fixed_pressure_gradient"}, {"pressure_gradient", gradient}};
  }
  if (!get(raw, "InletPressure").is_null()) action["left_port"]["pressure"] = get(raw, "InletPressure");
  if (!get(raw, "Temperature").is_null()) action["left_port"]["temperature"] = get(raw, "Temperature");
  if (!get(raw, "GasPhaseMolFractions").is_null())
    action["left_port"]["composition"] = normalized_fractions(get(raw, "GasPhaseMolFractions"), names);
  if (!get(d, "InletMolarFlowRate").is_null())
  {
    action["left_port"].erase("velocity");
    action["left_port"]["molar_flow_mol_s"] = get(d, "InletMolarFlowRate");
  }
  if (!get(d, "InletVolumetricFlowRate").is_null())
  {
    action["left_port"].erase("velocity");
    action["left_port"]["volume_flow_m3_s"] = get(d, "InletVolumetricFlowRate");
  }
  for (const auto key : {"ThermalControl", "LeftPort", "RightPort", "StreamTags", "FlowControl"})
    if (!get(raw, key).is_null()) action[aliases.at(lower(key))] = translated(get(raw, key), names);
  Json coolant = get(d, "Coolant");
  if (!coolant.is_null() && !action.contains("thermal") &&
      lower(get(d, "ThermalModel", "").get<std::string>()) != "adiabatic")
    action["thermal"] = {{"mode", "htf"},
                         {"temperature", get(coolant, "InletTemperature", 298.15)},
                         {"mass_flow_rate_kg_s", get(coolant, "MassFlowRate", 0.)},
                         {"heat_capacity_J_kg_K", get(coolant, "HeatCapacity", 4180.)},
                         {"heat_transfer_coefficient", get(coolant, "HeatTransferCoefficient", 0.)},
                         {"flow_direction", snake(get(coolant, "FlowDirection", "CoCurrent").get<std::string>())}};
  Json clean = Json::object();
  for (const auto key :
       {"name", "duration", "stop", "minimum_duration", "maximum_duration", "timeout_policy", "connections"})
    if (phase.contains(key)) clean[key] = phase[key];
  clean["columns"] = {{"A", action}};
  return clean;
}
}  // namespace
Json canonicalize_chemisorption_input(const Json& input)
{
  Json result = input;
  auto put=[](Json& object,const std::string& name,Json value){
    for(auto it=object.begin();it!=object.end();)if(lower(it.key())==lower(name))it=object.erase(it);else++it;
    object[name]=std::move(value);
  };
  auto chemical = get(result, "ChemicalKinetics");
  if (!chemical.is_null()) {
    auto components = get(result,"Components",Json::array());
    require(!components.empty(), "ChemicalKinetics import requires Components");
    std::vector<std::string> names;for(const auto& component:components)names.push_back(get(component,"Name").get<std::string>());
    auto spec = translated(chemical,names);
    std::function<void(Json&)> fix_reference = [&](Json& value) {
      if(value.is_object()) {
        if(value.contains("reference_temperature")){value["T_ref"]=value["reference_temperature"];value.erase("reference_temperature");}
        for(auto& child:value)fix_reference(child);
      } else if(value.is_array())for(auto& child:value)fix_reference(child);
    };
    fix_reference(spec);
    ChemistryKinetics legacy(spec,names.size());
    const bool partition = spec.contains("capacity_fraction");
    const double fraction = partition ? spec.at("capacity_fraction").get<double>() : 1.0;
    require(std::isfinite(fraction)&&fraction>=0.0&&fraction<=1.0,"ChemicalKinetics CapacityFraction must lie in [0,1]");
    require(!partition || !spec.contains("equilibrium") || spec.at("equilibrium").empty(),"CapacityFraction cannot accompany separate chemical equilibrium");
    const bool physical = get(result,"PhysisorptionEnabled",true);
    require(partition || !physical || spec.value("independent_sites",false),
      "Combined chemical and physical capacities require IndependentSites or CapacityFraction");
    for(std::size_t c=0;c<components.size();++c) {
      require(get(components[c],"ChemisorptionSites",Json::array()).empty(),"Use ChemisorptionSites; do not combine them with legacy ChemicalKinetics");
      auto site=legacy.canonical_site(c);
      if(partition)site["EquilibriumFraction"]=fraction;
      if(!site.contains("HeatOfAdsorption"))site["HeatOfAdsorption"]=get(components[c],"HeatOfAdsorption",0.0);
      put(components[c],"ChemisorptionSites",Json::array({site}));
    }
    put(result,"Components",components);
    if(get(result,"SpatialDiscretization").is_null())result["SpatialDiscretization"]="FDM";
    if(partition){put(result,"PhysisorptionFraction",1.0-fraction);put(result,"PhysisorptionEnabled",true);}
    for(auto it=result.begin();it!=result.end();)if(lower(it.key())=="chemicalkinetics")it=result.erase(it);else++it;
  }
  // Material and column overrides may inherit the common component definitions.
  for(const auto* key:{"Columns","Adsorbents"}) {
    auto children=get(result,key);
    if(children.is_object()||children.is_array()) {
      for(auto& child:children)if(child.is_object()&&!get(child,"ChemicalKinetics").is_null()) {
        if(get(child,"Components").is_null())child["Components"]=get(result,"Components",Json::array());
        child=canonicalize_chemisorption_input(child);
      }
      if(!children.is_null())put(result,key,children);
    }
  }
  return result;
}

bool is_process_input(const Json& d)
{
  const auto spatial = lower(get(d, "SpatialDiscretization", "").get<std::string>());
  if (spatial == "fvm" || spatial == "fdm" || !get(d,"Columns").is_null() || !get(d,"CyclicSteadyState").is_null() || !get(d,"ChemicalKinetics").is_null()) return true;
  for(const auto& component:get(d,"Components",Json::array()))for(const auto& site:get(component,"ChemisorptionSites",Json::array())) {
    const auto equation=lower(get(site,"RateEquation","").get<std::string>());
    if(equation=="dualnth"||equation=="dual_nth"||!get(site,"EquilibriumSource").is_null()||!get(site,"MixturePredictionMethod").is_null()||!get(site,"ReferenceTemperature").is_null()||!get(site,"EquilibriumHeat").is_null())return true;
    const auto parameters=get(site,"RateEquationParameters",Json::object());
    for(const auto* key:{"HistoryPolicy","AdsorptionRateCoefficient","DesorptionRateCoefficient","AdsorptionActivationEnergy","DesorptionActivationEnergy","RateBasis","DesorptionOrder"})if(!get(parameters,key).is_null())return true;
  }
  return false;
}
Json normalize_main_input(const Json& document, const std::filesystem::path& base)
{
  const Json d = canonicalize_chemisorption_input(document);
  const std::set<std::string> known = {"particlesphericity", "thermalmodel", "fluidsolidthermalequilibrium", "inletdensitybasis",
                                       "watersaturationmodel",
                                       "physisorptionenabled", "physisorptionfraction",
                                       "chemicalkinetics",
                                       "notes",
                                       "provenance",
                                       "boundaryhydrauliclength",
                                       "exportvalidationfluxes",
                                       "kineticconcentrationregularizer",
                                       "gasconstant",
                                       "adsorbents",
                                       "adsorptionreferencetemperature",
                                       "ambienttemperature",
                                       "automaticstopping",
                                       "boundarycondition",
                                       "breakthroughintegrator",
                                       "buffers",
                                       "cvodeabsolutetoleranceconcentration",
                                       "cvodeabsolutetoleranceloading",
                                       "cvodeabsolutetolerancetemperature",
                                       "cvodekrylovdimension",
                                       "cvodelinearsolver",
                                       "cvodemaximumorder",
                                       "cvodemaximumtimestep",
                                       "cvoderelativetolerance",
                                       "columndistances",
                                       "columnentrancevelocity",
                                       "columnerror",
                                       "columnlength",
                                       "columnloading",
                                       "columnpressure",
                                       "columnsections",
                                       "columnvoidfraction",
                                       "columns",
                                       "components",
                                       "convectionscheme",
                                       "coolant",
                                       "cyclicsteadystate",
                                       "debugforcemultibed",
                                       "displayname",
                                       "dynamicviscosity",
                                       "fluidphase",
                                       "gassolidthermalequilibrium",
                                       "gasstoragevoidfraction",
                                       "geometry",
                                       "iastmethod",
                                       "influxtemperature",
                                       "initialconditions",
                                       "initialpressure",
                                       "inletmolarflowrate",
                                       "inletpressure",
                                       "inletvolumetricflowrate",
                                       "ldfmodel",
                                       "liquiddensity",
                                       "mpdsettings",
                                       "mixturepredictionmethod",
                                       "momentummodel",
                                       "numberofcycles",
                                       "numberofgridpoints",
                                       "numberofinittimesteps",
                                       "numberofpressurepoints",
                                       "outletpressure",
                                       "output",
                                       "particledensity",
                                       "particlediameter",
                                       "particletransport",
                                       "performance",
                                       "poreldfrate",
                                       "pressureend",
                                       "pressuregradient",
                                       "pressurescale",
                                       "pressurestart",
                                       "reactions",
                                       "readcolumnfile",
                                       "relativehumidity",
                                       "relativehumiditymode",
                                       "restart",
                                       "simulationtime",
                                       "simulationtype",
                                       "solver",
                                       "sources",
                                       "spatialdiscretization",
                                       "stopcriterion",
                                       "swingadsorptionphases",
                                       "temperature",
                                       "timestep",
                                       "utilities",
                                       "wallendboundary",
                                       "wallendtemperatures",
                                       "watercomponent",
                                       "writeinterval",
                                       "energybalance",
                                       "gasthermalconductivity",
                                       "heatcapacitygas",
                                       "heatcapacityliquid",
                                       "heatcapacitysolid",
                                       "heatcapacitywall",
                                       "heattransfergassolid",
                                       "heattransfergaswall",
                                       "heattransferliquidsolid",
                                       "heattransferliquidwall",
                                       "heattransferwallexternal",
                                       "liquidthermalconductivity",
                                       "phcomponent",
                                       "phmode",
                                       "phvalue",
                                       "pkw",
                                       "walldensity",
                                       "wallthermalconductivity"};
  for (auto& [key, value] : d.items()) require(known.count(lower(key)), "Unknown RUPTURA input setting: " + key);
  auto boundary = lower(get(d, "BoundaryCondition", "InletVelocityOutletPressure").get<std::string>());
  require(boundary == "inletvelocityoutletpressure" || boundary == "inletpressureoutletpressure" ||
          boundary == "fixedvelocity" || boundary == "constantvelocity" ||
          boundary == "fixedpressuregradient" || boundary == "constantpressuregradient",
          "BoundaryCondition must be InletVelocityOutletPressure, InletPressureOutletPressure, FixedVelocity or FixedPressureGradient");
  require(get(d, "ReadColumnFile").is_null(),
          "FVM restart uses Restart and its complete checkpoint, including particle and coolant states");
  if(!get(d,"Reactions",Json::array()).empty()) {
    require(get(d,"SimulationType","Breakthrough")=="Breakthrough","General reaction networks are excluded from cycles");
    require(get(d,"ParticleTransport").is_null(),"ParticleTransport is restricted to nonreactive physisorption");
    if(get(d,"energyBalance",false).get<bool>())for(const auto& reaction:get(d,"Reactions"))
      require(lower(get(reaction,"Phase","").get<std::string>())!="bulkconcentration","BulkConcentration reactions do not provide reaction heat; use energyBalance false, as in original main");
  }
  require(get(d, "InletMolarFlowRate").is_null() || get(d, "ColumnEntranceVelocity").is_null(),
          "Specify InletMolarFlowRate or ColumnEntranceVelocity, not both");
  require(get(d, "InletVolumetricFlowRate").is_null() || get(d, "ColumnEntranceVelocity").is_null(),
          "Specify InletVolumetricFlowRate or ColumnEntranceVelocity, not both");
  auto components = get(d, "Components");
  require(components.is_array() && !components.empty(), "Components must be a nonempty ordered array");
  std::vector<std::string> names;
  Json public_components = Json::array(), feed = Json::array(), initial = Json::array(), cin = Json::array(),
       c0 = Json::array();
  std::string fluid = lower(get(d, "FluidPhase", "Gas").get<std::string>());
  require(fluid == "gas" || fluid == "liquid", "FluidPhase must be Gas or Liquid");
  bool liquid = fluid == "liquid";
  for (auto& x : components)
  {
    std::string name = get(x, "Name").get<std::string>();
    require(std::find(names.begin(), names.end(), name) == names.end(), "Duplicate component: " + name);
    names.push_back(name);
    double mw = get(x, "MolecularWeight", .028),
           cp = get(x, "GasHeatCapacity", get(d, "heatCapacityGas", 1000.).get<double>() * mw);
    if (!liquid)
      require(cp > 8.31446261815324,
              "GasHeatCapacity must exceed R for every FVM gas component, including the initial carrier: " + name);
    Json c = {{"name", name},
              {"molar_mass_kg_mol", mw},
              {"dispersion_m2_s", get(x, "AxialDispersionCoefficient", 0.)},
              {"gas_heat_capacity_J_mol_K", cp},
              {"adsorbed_heat_capacity_J_mol_K", get(x, "AdsorbedHeatCapacity", 0.)}};
    public_components.push_back(c);
    feed.push_back(get(x, "GasPhaseMolFraction", liquid ? 1. : 0.));
    initial.push_back(get(x, "CarrierGas", false).get<bool>() ? 1. : 0.);
    cin.push_back(get(x, "LiquidPhaseConcentration", 0.));
    c0.push_back(get(x, "InitialLiquidPhaseConcentration", 0.));
  }
  double s = 0;
  for (auto& x : initial) s += x.get<double>();
  if (s == 0) initial = feed;
  feed = normalized_fractions(feed, names);
  initial = normalized_fractions(initial, names);
  Json seed = column_seed(d, components, names, base);
  seed["components"] = public_components;
  if (!get(d, "ColumnDistances").is_null())
  {
    seed["advanced"]["COLUMN_DISTANCES"] = get(d, "ColumnDistances");
    seed["bed"]["nodes"] = get(d, "NumberOfGridPoints", get(d, "ColumnDistances").size() - 1).get<int>() + 1;
  }
  Json adsorbents = get(d, "Adsorbents", Json::array());
  if (!adsorbents.empty())
  {
    seed["materials"] = Json::object();
    Json order = Json::array(), mix = Json::array();
    for (auto& ads : adsorbents)
    {
      std::string name = get(ads, "Name").get<std::string>();
      require(!seed["materials"].contains(name), "Duplicate adsorbent: " + name);
      Json local = d;
      local.erase("Adsorbents");
      local.update(ads);
      for (const auto key : {"ColumnVoidFraction", "ParticleDiameter", "ParticleSphericity"})
        if (!get(ads, key).is_null()) local["Geometry"][key] = get(ads, key);
      Json comps = components;
      auto overrides = get(ads, "Components", get(ads, "ComponentParameters", Json::array()));
      if (overrides.is_object())
      {
        Json list = Json::array();
        for (auto& [name, entry] : overrides.items())
        {
          Json c = entry;
          c["Name"] = name;
          list.push_back(c);
        }
        overrides = list;
      }
      for (auto& override : overrides)
      {
        std::string cname = get(override, "Name").get<std::string>();
        auto it = std::find(names.begin(), names.end(), cname);
        require(it != names.end(), "Adsorbent names an unknown component: " + cname);
        comps[std::size_t(it - names.begin())].update(override);
      }
      seed["materials"][name] = column_seed(local, comps, names, base);
      order.push_back(name);
      if (!get(ads, "MixFraction").is_null()) mix.push_back(get(ads, "MixFraction"));
    }
    seed["advanced"]["MATERIAL_ORDER"] = order;
    Json sections = get(d, "ColumnSections", Json::array());
    if (!sections.empty())
    {
      require(mix.empty(), "MixFraction cannot be combined with ColumnSections");
      Json beds = Json::array();
      Vec pure, interfaces;
      bool explicit_grid = false;
      for (auto& section : sections) explicit_grid |= !get(section, "NumberOfGridPoints").is_null();
      Vec coordinates = {0.};
      double position = 0.;
      bool want_bed = true;
      for (auto& section : sections)
      {
        bool is_interface = !get(section, "InterfaceLength").is_null();
        require(is_interface != want_bed, "ColumnSections must alternate adsorbent Length and InterfaceLength entries");
        double length = get(section, is_interface ? "InterfaceLength" : "Length", 0.);
        require(std::isfinite(length) && (is_interface ? length >= 0. : length > 0.), "Invalid ColumnSections length");
        if (is_interface)
          interfaces.push_back(length);
        else
        {
          require(beds.size() < order.size() && get(section, "Adsorbent") == order[beds.size()],
                  "ColumnSections adsorbent order must match Adsorbents");
          beds.push_back(section);
          pure.push_back(length);
        }
        if (explicit_grid)
        {
          auto count = get(section, "NumberOfGridPoints", 0);
          require(count.is_number_integer(), "ColumnSections NumberOfGridPoints must be an integer");
          int n = count.get<int>();
          require(n >= 0 && (length == 0. ? n == 0 : n > 0),
                  "Set NumberOfGridPoints for every positive-length section");
          for (int k = 1; k <= n; ++k) coordinates.push_back(position + length * k / n);
        }
        position += length;
        want_bed = !want_bed;
      }
      require(!want_bed && beds.size() == order.size(),
              "ColumnSections must end with an adsorbent and include every material");
      seed["layers"] = Json::array();
      for (std::size_t i = 0; i < beds.size(); ++i)
      {
        double before = i ? interfaces[i - 1] : 0., after = i < interfaces.size() ? interfaces[i] : 0.;
        seed["layers"].push_back(
            {{"material", order[i]}, {"length_m", pure[i] + .5 * (before + after)}, {"interface_after_m", after}});
      }
      seed["bed"]["length_m"] = position;
      if (explicit_grid)
      {
        require(get(d, "ColumnDistances").is_null(), "Use ColumnDistances or per-section NumberOfGridPoints");
        seed["advanced"]["COLUMN_DISTANCES"] = coordinates;
        seed["bed"]["nodes"] = coordinates.size();
      }
    }
    else if (!mix.empty())
    {
      require(mix.size() == adsorbents.size(), "Set MixFraction for every adsorbent");
      seed["advanced"]["MIX_FRACTIONS"] = mix;
    }
    else
    {
      seed["layers"] = Json::array();
      double length = 0.;
      for (auto& ads : adsorbents)
      {
        double l = get(ads, "AdsorbentLength", 0.);
        require(l > 0., "Multiple adsorbents require ColumnSections, AdsorbentLength or MixFraction");
        seed["layers"].push_back({{"material", get(ads, "Name")}, {"length_m", l}});
        length += l;
      }
      seed["bed"]["length_m"] = length;
    }
  }

  seed["gas"] = {
      {"viscosity_Pa_s", get(d, "DynamicViscosity", liquid ? .001 : 1.8e-5)},
      {"heat_capacity_J_kg_K", get(d, "heatCapacityGas", 1000.)},
      {"thermal_conductivity_W_m_K", get(d, "gasThermalConductivity", get(d, "liquidThermalConductivity", .025))}};
  seed["feed"] = {{"pressure_Pa", get(d, "InletPressure", 1e5)},
                  {"temperature_K", get(d, "InfluxTemperature", get(d, "Temperature", 298.15))},
                  {"composition", feed},
                  {"velocity_basis", "interstitial"}};
  Json ini = get(d, "InitialConditions", Json::object());
  known_options(ini,
                {"type", "pressure", "temperature", "gastemperature", "solidtemperature", "walltemperature",
                 "jackettemperature", "gasphasemolfractions", "liquidphaseconcentrations", "loading"},
                "InitialConditions");
  Json explicit_initial = Json::object();
  for (auto& pair : std::map<std::string, std::string>{{"Pressure", "pressure"},
                                                       {"GasTemperature", "gas_temperature"},
                                                       {"SolidTemperature", "solid_temperature"},
                                                       {"WallTemperature", "wall_temperature"},
                                                       {"JacketTemperature", "jacket_temperature"},
                                                       {"Loading", "loading"}})
    if (!get(ini, pair.first).is_null()) explicit_initial[pair.second] = get(ini, pair.first);
  for (auto [public_key, internal_key] : std::map<std::string, std::string>{{"GasTemperature", "gas_temperature"},
          {"SolidTemperature", "solid_temperature"}, {"WallTemperature", "wall_temperature"}})
    if (get(ini, public_key).is_null() && !get(ini, "Temperature").is_null())
      explicit_initial[internal_key] = get(ini, "Temperature");
  if (explicit_initial.contains("loading") && explicit_initial["loading"].is_object())
    explicit_initial["loading"] = component_values(explicit_initial["loading"], names);
  if (!get(ini, "GasPhaseMolFractions").is_null())
    explicit_initial["composition"] = normalized_fractions(get(ini, "GasPhaseMolFractions"), names);
  auto initial_type = lower(get(ini, "Type", "Custom").get<std::string>());
  require(initial_type == "custom" || initial_type == "empty" || initial_type == "cleansolvent" ||
          initial_type == "carriergas" || initial_type == "feedpressurized" || initial_type == "gaspressurized" ||
          initial_type == "evacuated", "Unknown InitialConditions.Type");
  require(!liquid || initial_type == "custom" || initial_type == "empty" || initial_type == "cleansolvent",
          "Liquid initial conditions use Custom, Empty or CleanSolvent; gas presets require FluidPhase Gas");
  require(liquid || initial_type != "cleansolvent", "CleanSolvent requires FluidPhase Liquid");
  if (initial_type == "carriergas")
  {
    bool carrier = false;
    for (const auto& c : components) carrier |= get(c, "CarrierGas", false).get<bool>();
    require(carrier, "CarrierGas initial conditions require a marked carrier component");
    require(get(ini, "GasPhaseMolFractions").is_null(), "CarrierGas preset uses the marked carrier composition");
  }
  if (initial_type == "feedpressurized")
  {
    require(get(ini, "GasPhaseMolFractions").is_null(), "FeedPressurized preset uses the feed composition");
    explicit_initial["composition"] = feed;
  }
  if (initial_type == "gaspressurized")
    require(!get(ini, "GasPhaseMolFractions").is_null(), "GasPressurized requires GasPhaseMolFractions");
  if (initial_type == "empty" || initial_type == "cleansolvent")
  {
    require(get(ini, "Loading", "zero") == "zero", "Empty/CleanSolvent requires Loading zero");
    explicit_initial["loading"] = "zero";
    if (liquid)
      require(get(ini, "LiquidPhaseConcentrations").is_null(), "Empty/CleanSolvent uses zero tracked solute concentrations");
  }
  if (initial_type == "evacuated")
  {
    explicit_initial["pressure"] = get(ini, "Pressure", get(d, "InitialPressure", 1.));
    const auto& pressure = explicit_initial["pressure"];
    if (pressure.is_array())
      for (const auto& value : pressure) require(value.get<double>() > 0., "Evacuated requires positive residual pressure");
    else require(pressure.get<double>() > 0., "Evacuated requires positive residual pressure");
  }

  seed["initial"] = {
      {"pressure_Pa", first_number(get(ini, "Pressure"), get(d, "InitialPressure", get(d, "InletPressure", 1e5)))},
      {"gas_temperature_K",
       first_number(get(ini, "GasTemperature", get(ini, "Temperature")), get(d, "Temperature", 298.15))},
      {"solid_temperature_K", first_number(get(ini, "SolidTemperature"), get(d, "Temperature", 298.15))},
      {"wall_temperature_K", first_number(get(ini, "WallTemperature"), get(d, "Temperature", 298.15))},
      {"composition", initial},
      {"loading", "zero"}};
  seed["sources"] = {{"outlet",
                      {{"pressure_Pa", get(d, "OutletPressure", get(d, "InletPressure", 1e5))},
                       {"temperature_K", get(d, "Temperature", 298.15)},
                       {"composition", feed}}}};
  Json sources = get(d, "Sources", Json::object());
  for (auto& [name, x] : sources.items()) seed["sources"][name] = translated(x, names);
  Json advanced = seed["advanced"];
  advanced["SIMULATION_TYPE"] = get(d, "SimulationType", "Breakthrough");
  if (!get(get(d, "Performance", Json::object()), "FreshFeedTags").is_null())
    advanced["FRESH_FEED_TAGS"] = get(get(d, "Performance", Json::object()), "FreshFeedTags");
  advanced["FLUID_PHASE"] = liquid ? "liquid" : "gas";
  advanced["LIQUID_DENSITY"] = get(d, "LiquidDensity", 1000.);
  advanced["LIQUID_HEAT_CAPACITY"] = get(d, "heatCapacityLiquid", 4180.);
  advanced["LIQUID_FEED_CONCENTRATION"] = cin;
  advanced["LIQUID_INITIAL_CONCENTRATION"] = component_values(get(ini, "LiquidPhaseConcentrations", c0), names);
  if (liquid && (initial_type == "empty" || initial_type == "cleansolvent"))
    advanced["LIQUID_INITIAL_CONCENTRATION"] = Vec(names.size(), 0.);
  if (!get(d, "Coolant").is_null())
  {
    const auto c = get(d, "Coolant");
    seed["thermal"]["htf"] = {{"model", snake(get(c, "Model", "Transient").get<std::string>())},
                              {"volume_m3", get(c, "Volume", .001)},
                              {"density_kg_m3", get(c, "Density", 1000.)},
                              {"heat_capacity_J_kg_K", get(c, "HeatCapacity", 4180.)},
                              {"initial_temperature_K", get(c, "InitialTemperature", get(d, "Temperature", 298.15))}};
  }
  Json css = get(d, "CyclicSteadyState", Json::object());
  const auto spatial = lower(get(d, "SpatialDiscretization", "FVM").get<std::string>());
  require(spatial == "fdm" || spatial == "fvm", "Shared simulations require SpatialDiscretization FDM or FVM");
  advanced["SPATIAL_METHOD"] = spatial == "fdm" ? "FDM" : "FVM";
  advanced["SIMULATION_TYPE"] = get(d, "SimulationType", "Breakthrough");
  bool untilcss = get(css, "Enabled", false);
  seed["simulation"] = {{"name", get(d, "DisplayName", "RUPTURA")},
                        {"run_mode", untilcss ? "css" : "fixed_cycles"},
                        {"cycles", get(d, "NumberOfCycles", 1)},
                        {"maximum_cycles", get(css, "MaximumCycles", 200)}};
  seed["css"] = {{"state_rtol", get(css, "StateRelativeTolerance", 1e-4)},
                 {"state_atol", get(css, "StateAbsoluteTolerance", 1e-7)},
                 {"inventory_rtol", get(css, "InventoryRelativeTolerance", 1e-5)},
                 {"metric_rtol", get(css, "MetricRelativeTolerance", 1e-4)},
                 {"consecutive_cycles", get(css, "ConsecutiveCycles", 3)}};
  std::string method = get(d, "BreakthroughIntegrator", "BDF");
  if (method == "RungeKutta3" || method == "SSP-RK3") method = "SSPRK3";
  if (method == "SemiImplicitSSP-RK3") method = "SemiImplicitSSPRK3";
  require(method == "BDF" || method == "BDF1" || method == "CVODE" || method == "SSPRK3" || method == "SemiImplicitSSPRK3",
          "BreakthroughIntegrator must be BDF, BDF1, CVODE, SSPRK3 or SemiImplicitSSPRK3 for shared simulations");
  std::string scheme = get(d, "ConvectionScheme", "TVD");
  require(scheme == "TVD" || scheme == "Upwind" || scheme == "WENO3" || scheme == "WENO5" || scheme == "WENO",
          "ConvectionScheme must be Upwind, TVD, WENO3 or WENO5");
  Json solver = get(d, "Solver", Json::object());
  seed["solver"] = {{"method", method},
                    {"rtol", get(solver, "RelativeTolerance", get(d, "CVODERelativeTolerance", 1e-6))},
                    {"atol", get(solver, "AbsoluteTolerance", 1e-8)},
                    {"maximum_step_s", get(solver, "MaximumTimeStep", get(d, "TimeStep", 1.))},
                    {"reconstruction", scheme == "Upwind"  ? "UPWIND"
                                       : scheme == "WENO5" ? "WENO"
                                                           : scheme},
                    {"jacobian_mode", get(solver, "JacobianMode", "structured_fd")},
                    {"cache_rhs", get(solver, "CacheRHS", true)},
                    {"event_time_tolerance_s", get(solver, "EventTimeTolerance", 1e-6)}};
  const std::map<std::string, std::string> cssmap = {{"DiagnosticWindow", "CSS_DIAGNOSTIC_WINDOW"},
                                                     {"DiagnosticMinimumCycles", "CSS_DIAGNOSTIC_MINIMUM_CYCLES"},
                                                     {"MaximumPeriod", "CSS_MAXIMUM_PERIOD"},
                                                     {"PlotCycleSelection", "CSS_PLOT_CYCLE_SELECTION"},
                                                     {"MassBalanceTolerance", "CSS_MASS_BALANCE_TOLERANCE"}};
  for (auto& [key, target] : cssmap)
    if (!get(css, key).is_null()) advanced[target] = get(css, key);
  Json output = get(d, "Output", Json::object());
  seed["output"] = {{"directory", get(output, "Directory", "results")},
                    {"profile", get(output, "Profile", "compact")},
                    {"plots", get(output, "Plots", true)},
                    {"sample_interval_s", get(d, "WriteInterval", 1.)},
                    {"checkpoint", get(output, "Checkpoint", true)},
                    {"phase_archives", false},
                    {"selected_cycle_archive", false},
                    {"compress_npz", true}};
  if (!get(solver, "FirstTimeStep").is_null()) advanced["SOLVER_FIRST_STEP"] = get(solver, "FirstTimeStep");
  advanced["RK_INITIALIZATION_STEPS"] = get(d, "NumberOfInitTimeSteps", 0);
  advanced["RK_ADAPTIVE"] = get(solver, "Adaptive", true);
  if (!get(solver, "MaximumSteps").is_null()) advanced["SOLVER_MAXIMUM_STEPS"] = get(solver, "MaximumSteps");
  advanced["SOLVER_NONNEGATIVE_CONSTRAINTS"]=get(solver,"NonnegativeConstraints",true);
  if(!get(solver,"IsothermPowerRegularizer").is_null())advanced["ISOTHERM_POWER_REG"]=get(solver,"IsothermPowerRegularizer");
  advanced["T_FLOOR"] = get(solver, "MinimumTemperature", 1.);
  advanced["T_CEILING"] = get(solver, "MaximumTemperature", 5000.);
  advanced["CVODE_LINEAR_SOLVER"] = get(d, "CVODELinearSolver", "Dense");
  advanced["CVODE_KRYLOV_DIMENSION"] = get(d, "CVODEKrylovDimension", 30);
  advanced["CVODE_MAXIMUM_ORDER"] = get(d, "CVODEMaximumOrder", 5);
  advanced["ATOL_CONCENTRATION"] = get(d, "CVODEAbsoluteToleranceConcentration", 1e-8);
  advanced["ATOL_LOADING"] = get(d, "CVODEAbsoluteToleranceLoading", 1e-9);
  advanced["ATOL_TEMPERATURE"] = get(d, "CVODEAbsoluteToleranceTemperature", 1e-5);
  if (get(d, "CVODEMaximumTimeStep", 0.).get<double>() > 0.)
    seed["solver"]["maximum_step_s"] = get(d, "CVODEMaximumTimeStep");
  advanced["WRITE_MAIN_FORMAT"] = get(output, "MainFormatFiles", true);
  advanced["RESTART"] = get(d, "Restart", false);
  advanced["RESTART_FROM_CHECKPOINT"] = get(d, "Restart", false);
  seed["output"]["samples_per_phase"] = get(output, "SamplesPerPhase", 1000000);
  seed["output"]["history_cycles"] = get(output, "HistoryCycles", 22);
  seed["solver"]["checkpoint_interval_s"] = get(output, "CheckpointInterval", 60.);
  Json phases = Json::array();
  std::string type = get(d, "SimulationType", "Breakthrough");
  if (type == "Breakthrough")
  {
    Json time = get(d, "SimulationTime", 100.), automatic = get(d, "AutomaticStopping", Json::object());
    const bool auto_time = time.is_string() && lower(time.get<std::string>()) == "auto";
    require(time.is_number() || auto_time, "SimulationTime must be a positive number or auto");
    double maximum = auto_time ? get(automatic, "MaximumTime", 2e5).get<double>() : time.get<double>();
    Json p = {{"Name", "Breakthrough"}, {"Duration", maximum}};
    if (auto_time || !get(d, "StopCriterion").is_null())
    {
      p.erase("Duration");
      p["MaximumDuration"] = maximum;
      p["MinimumDuration"] = get(automatic, "MinimumTime", 0.);
      Json stop = get(d, "StopCriterion");
      if (stop.is_null())
      {
        double threshold = get(automatic, "Threshold", .999),
               loading_tol = get(automatic, "LoadingRelativeTolerance", .001);
        require(threshold > 0. && threshold <= 1. && loading_tol > 0.,
                "Automatic stopping tolerances must be positive; Threshold must not exceed one");
        auto selected = get(automatic, "Components", names).get<std::vector<std::string>>();
        Json criteria = Json::array();
        double total = 0.;
        for (auto& v : feed) total += v.get<double>();
        for (const auto& name : selected)
        {
          auto it = std::find(names.begin(), names.end(), name);
          require(it != names.end(), "AutomaticStopping names an unknown component: " + name);
          auto j = std::size_t(it - names.begin());
          double reference = liquid ? cin[j].get<double>() : feed[j].get<double>() / total;
          if (reference <= 0.) continue;
          criteria.push_back({{"metric", liquid ? "outlet_concentration_ratio" : "outlet_mole_fraction"},
                              {"component", name},
                              {"operator", ">="},
                              {"value", liquid ? threshold : reference * threshold},
                              {"column", "A"}});
          criteria.push_back({{"metric", "loading_equilibrium_error"},
                              {"component", name},
                              {"operator", "<="},
                              {"value", loading_tol},
                              {"column", "A"}});
        }
        require(!criteria.empty(), "AutomaticStopping requires at least one positively fed component");
        stop = {{"logic", "all"}, {"criteria", criteria}};
      }
      p["StopCriterion"] = stop;
    }
    phases.push_back(phase_document(p, d, names));
    double tail = auto_time ? get(automatic, "ExtraFraction", .1).get<double>() : 0.;
    require(tail >= 0. && std::isfinite(tail), "ExtraFraction must be finite and nonnegative");
    if (tail > 0.)
    {
      auto tail_phase = phase_document(Json{{"Name", "Post-saturation observation"}, {"Duration", 1e-12}}, d, names);
      tail_phase["duration_fraction_of_previous"] = tail;
      phases.push_back(tail_phase);
    }
  }
  else
  {
    require(type == "SwingAdsorption", "Shared FVM mode supports Breakthrough or SwingAdsorption");
    for (auto& p : get(d, "SwingAdsorptionPhases", Json::array())) phases.push_back(phase_document(p, d, names));
    require(!phases.empty(), "SwingAdsorptionPhases must contain at least one phase");
  }
  Json performance = get(d, "Performance", Json::object());
  const std::map<std::string, std::string> perfmap = {{"FeedTag", "feed_tag"},
                                                      {"LightProductTag", "light_product_tag"},
                                                      {"HeavyProductTag", "heavy_product_tag"},
                                                      {"LightComponent", "light_component"},
                                                      {"HeavyComponent", "heavy_component"}};
  seed["performance"] = Json::object();
  for (auto& [a, b] : perfmap)
    if (!get(performance, a).is_null()) seed["performance"][b] = get(performance, a);
  Json utilities = get(d, "Utilities", Json::object());
  seed["utilities"] = Json::object();
  const std::map<std::string, std::string> umap = {{"ProductDeliveryPressure", "product_delivery_pressure_Pa"},
                                                   {"ProductBufferPressure", "product_buffer_pressure_Pa"},
                                                   {"PlantFreshFeedTarget", "plant_fresh_feed_target_Nm3_h"},
                                                   {"StandardTemperature", "standard_temperature_K"},
                                                   {"StandardPressure", "standard_pressure_Pa"},
                                                   {"CompressionModel", "compression_model"},
                                                   {"HeatCapacityRatio", "heat_capacity_ratio"},
                                                   {"ElectricityCostPerKWh", "electricity_cost_per_kWh"},
                                                   {"HeatCostPerKWh", "heat_cost_per_kWh"},
                                                   {"CoolingCostPerKWh", "cooling_cost_per_kWh"},
                                                   {"CoolingCostBasis", "cooling_cost_basis"},
                                                   {"AnnualOperatingHours", "annual_operating_hours"},
                                                   {"Availability", "availability"},
                                                   {"AmbientPressure", "ambient_pressure_Pa"},
                                                   {"FeedSupplyPressure", "feed_supply_pressure_Pa"},
                                                   {"ChargeFeedCompression", "charge_feed_compression"},
                                                   {"ChargeRecycleCompression", "charge_recycle_compression"},
                                                   {"CompressorEfficiency", "compressor_efficiency"},
                                                   {"VacuumEfficiency", "vacuum_efficiency"},
                                                   {"HeatingEfficiency", "heating_efficiency"},
                                                   {"CoolingCOP", "cooling_COP"},
                                                   {"HeatRecoveryFraction", "heat_recovery_fraction"}};
  for (auto& [a, b] : umap)
    if (!get(utilities, a).is_null()) seed["utilities"][b] = get(utilities, a);
  if (liquid)
  {
    advanced["LIQUID_PUMP_EFFICIENCY"] = get(utilities, "PumpEfficiency", .75);
    seed["utilities"]["charge_feed_compression"] = false;
    seed["utilities"]["charge_recycle_compression"] = false;
  }
  advanced["PHYSISORPTION_ENABLED"] = get(d, "PhysisorptionEnabled", true);
  advanced["PHYSISORPTION_FRACTION"] = get(d,"PhysisorptionFraction",1.0);
  for(const auto& component:get(d,"Components",Json::array()))
    if(!get(component,"ChemisorptionSites",Json::array()).empty())advanced["CHEMISORPTION_SCHEMA_VERSION"] = 2;
  seed["cycle"] = {{"phases", phases}};
  seed["columns"] = {{"A", Json::object()}};
  seed["buffers"] = translated(get(d, "Buffers", Json::object()), names, true);
  seed["advanced"] = advanced;
  const std::map<std::string, std::string> physical_options = {
      {"InletDensityBasis", "INLET_DENSITY_BASIS"},
      {"IASTMethod", "IAST_METHOD"},
      {"WaterSaturationModel", "WATER_PSAT_MODEL"},
      {"GasConstant", "R"},
      {"KineticConcentrationRegularizer", "KINETIC_CONCENTRATION_REGULARIZER"},
      {"ExportValidationFluxes", "EXPORT_VALIDATION_FLUXES"},
      {"BoundaryHydraulicLength", "BOUNDARY_HYDRAULIC_LENGTH"},
      {"pHValue", "PH_VALUE"},
      {"MomentumModel", "PRESSURE_DROP_MODEL"},
      {"GasSolidThermalEquilibrium", "GAS_SOLID_THERMAL_EQUILIBRIUM"},
      {"GasStorageVoidFraction", "GAS_STORAGE_VOID_FRACTION"},
      {"WallEndBoundary", "WALL_END_BOUNDARY"},
      {"WallEndTemperatures", "WALL_END_TEMPERATURE_K"},
      {"LDFModel", "LDF_MODEL"},
      {"PoreLDFRate", "LDF_PORE_RATE_S_INV"},
      {"WaterComponent", "WATER_COMPONENT"},
      {"RelativeHumidityMode", "RH_MODE"},
      {"RelativeHumidity", "RH_VALUE"},
      {"AdsorptionReferenceTemperature", "ADSORPTION_REFERENCE_TEMPERATURE"}};
  for (auto& [key, target] : physical_options)
    if (!get(d, key).is_null())
    {
      auto value = get(d, key);
      if ((target == "INLET_DENSITY_BASIS" || target == "PRESSURE_DROP_MODEL" || target == "LDF_MODEL" ||
           target == "WALL_END_BOUNDARY") &&
          value.is_string())
        value = snake(value.get<std::string>());
      seed["advanced"][target] = value;
    }
  if (!get(d, "FluidSolidThermalEquilibrium").is_null())
  {
    auto common = get(d, "FluidSolidThermalEquilibrium");
    require(get(d, "GasSolidThermalEquilibrium").is_null() || common == get(d, "GasSolidThermalEquilibrium"),
            "FluidSolidThermalEquilibrium conflicts with GasSolidThermalEquilibrium");
    seed["advanced"]["GAS_SOLID_THERMAL_EQUILIBRIUM"] = common;
  }
  if (!get(solver, "DispersionClosure").is_null())
    seed["solver"]["dispersion_closure"] = get(solver, "DispersionClosure");
  Json cfg = normalize_config(seed, base);
  cfg["INITIAL_STATE"].update(explicit_initial);
  cfg["INITIAL_STATE"]["concentration"] = advanced["LIQUID_INITIAL_CONCENTRATION"];
  cfg["NOTES"] = get(d, "Notes", Json::array());
  cfg["INPUT_DOCUMENT"] = d;
  cfg["INPUT_STYLE"] = "RUPTURA-main";
  auto output_path = std::filesystem::weakly_canonical(std::filesystem::path(cfg["OUTPUT_DIR"].get<std::string>())),
       case_path = std::filesystem::weakly_canonical(std::filesystem::absolute(base));
  auto relation = case_path.lexically_relative(output_path);
  require(relation.empty() || *relation.begin() == "..",
          "Output.Directory must be a dedicated results folder, not the input folder or its parent");
  Json explicit_columns = get(d, "Columns");
  if (!explicit_columns.is_null())
  {
    require(explicit_columns.is_object() && !explicit_columns.empty(), "Columns must be a named object");
    Json cols = Json::object();
    for (auto& [name, override] : explicit_columns.items())
    {
      Json local = d;
      local.erase("Columns");
      local["SimulationType"] = "Breakthrough";
      local["SimulationTime"] = 1.;
      local.erase("CyclicSteadyState");
      local.update(override);
      auto child = normalize_main_input(local, base);
      // Physical overrides are parsed as a single operation, then inherit the
      // enclosing simulation's scope. A cycle column must not acquire the
      // breakthrough-only reaction capability through this parsing step.
      child["SIMULATION_TYPE"] = cfg["SIMULATION_TYPE"];
      for (auto it = child.begin(); it != child.end();)
      {
        auto key = it.key();
        if (key == "INPUT_DOCUMENT" || key == "INPUT_STYLE" || key == "INPUT_FILE" || key == "CASE_DIR" ||
            key == "OUTPUT_DIR" || key == "NOTES" || key == "CASE_NAME" || key == "DESCRIPTION" || key == "RUN_MODE" ||
            key == "N_CYCLES" || key == "MAX_CYCLES" || key == "CYCLE_PHASES" || key == "COLUMNS" || key == "RESTART" ||
            key == "RESTART_FROM_CHECKPOINT" || key == "MAKE_PLOTS" || key == "CHECKPOINT_ENABLED" ||
            key.starts_with("OUTPUT_") || key.starts_with("SOLVER_") || key.starts_with("CSS_") ||
            key.starts_with("WRITE_"))
          it = child.erase(it);
        else
          ++it;
      }
      cols[name] = child;
    }
    cfg["COLUMNS"] = cols;
  }
  return cfg;
}
Json load_main_input(const std::filesystem::path& path)
{
  auto c = normalize_main_input(read_input_document(path), std::filesystem::absolute(path).parent_path());
  c["INPUT_FILE"] = std::filesystem::absolute(path).string();
  return c;
}
Json process_capabilities()
{
  return Json{{"fluid_phases", {"Gas", "Liquid"}},
              {"contactors", {"PackedBed", "Monolith"}},
              {"simulation_types", {"Breakthrough", "SwingAdsorption"}},
              {"discretization", {"FDM", "FVM"}},
              {"convection", {"Upwind", "TVD", "WENO3", "WENO5"}},
              {"integrators", {"BDF1", "BDF", "CVODE", "SSPRK3", "SemiImplicitSSPRK3"}},
              {"mixture_methods", {"SPI", "SCI", "EI", "SEI", "IAST", "SIAST", "MPD"}},
              {"mpd", {{"fluid_phase", "Gas"}, {"configuration", "MPDSettings"},
                        {"discretization", {"FDM", "FVM"}}}},
              {"flow_controls", {{"FixedVelocity", {"Gas", "Liquid"}}, {"FixedPressureGradient", {"Liquid"}}}},
              {"thermal_models", {"Isothermal", "Adiabatic", "NonAdiabatic"}},
              {"fluid_solid_temperature", {"Separate", "Common"}},
              {"initial_condition_presets", {"Empty", "CleanSolvent", "CarrierGas", "FeedPressurized", "GasPressurized", "Evacuated", "Custom"}},
              {"particle_transport_fixed_beds", {"LDF", "HSDM"}},
              {"chemisorption", {{"configuration", "ChemisorptionSites"}, {"fluid_phases", {"Gas", "Liquid"}},
                {"rate_equations", {"FirstOrder", "PseudoNth", "Avrami", "General", "Elovich", "DualNth"}},
                {"history_policies", {"ElapsedTime", "ContinuousErasure"}}, {"legacy_import", "ChemicalKinetics"}}}};
}
}  // namespace cycle
namespace cycle
{
void write_main_outputs(const Json& config, const std::filesystem::path& destination)
{
  auto root = std::filesystem::path(config.at("OUTPUT_DIR").get<std::string>()) / "main";
  if (!std::filesystem::is_directory(root)) return;
  std::vector<std::filesystem::path> columns;
  for (auto& entry : std::filesystem::directory_iterator(root))
    if (entry.is_directory()) columns.push_back(entry.path());
  std::sort(columns.begin(), columns.end());
  for (std::size_t i = 0; i < columns.size(); ++i)
  {
    auto target = columns.size() == 1 ? destination : destination / "columns" / columns[i].filename();
    std::filesystem::create_directories(target);
    for (auto& entry : std::filesystem::directory_iterator(columns[i]))
      if (entry.is_regular_file())
        std::filesystem::copy_file(entry.path(), target / entry.path().filename(),
                                   std::filesystem::copy_options::overwrite_existing);
    if (i == 0 && columns.size() > 1)
      for (auto& entry : std::filesystem::directory_iterator(columns[i]))
        if (entry.is_regular_file())
          std::filesystem::copy_file(entry.path(), destination / entry.path().filename(),
                                     std::filesystem::copy_options::overwrite_existing);
  }
}
}  // namespace cycle

#include "inputreader.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string_view>

#include "json.h"

/**
 * \brief Parses a JSON stream and reports the byte offset on parse failure.
 */
nlohmann::json readJson(std::ifstream& fileInput)
{
  try
  {
    return nlohmann::json::parse(fileInput);
  }
  catch (const nlohmann::json::parse_error& ex)
  {
    std::cerr << "parse error at byte " << ex.byte << std::endl;
    throw;
  }
}

/**
 * \brief Case-insensitive string comparison for user-facing input tokens.
 */
bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
  return str1.size() == str2.size() &&
         std::equal(str1.begin(), str1.end(), str2.begin(),
                    [](unsigned char a, unsigned char b) { return std::tolower(a) == std::tolower(b); });
}

/**
 * \brief Finds a JSON object key using exact lookup first, then case-insensitive lookup.
 */
static const nlohmann::json* findKeyCaseInsensitive(const nlohmann::json& object, const std::string& key)
{
  if (!object.is_object())
  {
    return nullptr;
  }

  auto it = object.find(key);
  if (it != object.end())
  {
    return &(*it);
  }

  for (auto it2 = object.begin(); it2 != object.end(); ++it2)
  {
    if (caseInSensStringCompare(it2.key(), key))
    {
      return &(*it2);
    }
  }

  return nullptr;
}

static bool containsKeyCaseInsensitive(const nlohmann::json& object, const std::string& key)
{
  return findKeyCaseInsensitive(object, key) != nullptr;
}

static const nlohmann::json& requireKeyCaseInsensitive(const nlohmann::json& object, const std::string& key,
                                                       const std::string& context)
{
  const nlohmann::json* value = findKeyCaseInsensitive(object, key);
  if (value == nullptr)
  {
    throw std::runtime_error("Error: required key '" + key + "' missing" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }
  return *value;
}

template <typename T>
static T getNumberOrThrow(const nlohmann::json& value, const std::string& key, const std::string& context)
{
  try
  {
    return value.get<T>();
  }
  catch (const nlohmann::json::exception&)
  {
    throw std::runtime_error("Error: key '" + key + "' has invalid value" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }
}

static std::string getStringOrThrow(const nlohmann::json& value, const std::string& key, const std::string& context)
{
  try
  {
    return value.get<std::string>();
  }
  catch (const nlohmann::json::exception&)
  {
    throw std::runtime_error("Error: key '" + key + "' has invalid value" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }
}

static bool getBoolOrThrow(const nlohmann::json& value, const std::string& key, const std::string& context)
{
  if (value.is_boolean())
  {
    return value.get<bool>();
  }

  // Accept legacy-like encodings from older input files.
  if (value.is_string())
  {
    std::string s = value.get<std::string>();
    if (caseInSensStringCompare(s, "yes")) return true;
    if (caseInSensStringCompare(s, "no")) return false;
    if (caseInSensStringCompare(s, "true")) return true;
    if (caseInSensStringCompare(s, "false")) return false;
  }

  throw std::runtime_error("Error: key '" + key + "' has invalid boolean value" +
                           (context.empty() ? "" : (" (" + context + ")")));
}

template <typename T>
static std::vector<T> getNumberListOrThrow(const nlohmann::json& value, const std::string& key,
                                           const std::string& context)
{
  if (!value.is_array())
  {
    throw std::runtime_error("Error: key '" + key + "' must be an array" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  try
  {
    return value.get<std::vector<T>>();
  }
  catch (const nlohmann::json::exception&)
  {
    throw std::runtime_error("Error: key '" + key + "' must be an array of numbers" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }
}

static std::vector<double> requireDoubleParameterCount(std::vector<double> values, std::size_t n,
                                                       const std::string& name, const std::string& context)
{
  if (values.size() < n)
  {
    throw std::runtime_error("Error: " + name + " requires " + std::to_string(n) + " parameters" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  values.resize(n);
  return values;
}

template <typename T>
static void readOptionalNumber(const nlohmann::json& object, const std::string& key, T& target)
{
  if (containsKeyCaseInsensitive(object, key))
  {
    target = getNumberOrThrow<T>(requireKeyCaseInsensitive(object, key, ""), key, "");
  }
}

static void readOptionalString(const nlohmann::json& object, const std::string& key, std::string& target)
{
  if (containsKeyCaseInsensitive(object, key))
  {
    target = getStringOrThrow(requireKeyCaseInsensitive(object, key, ""), key, "");
  }
}

static void readOptionalString(const nlohmann::json& object, const std::string& key, std::optional<std::string>& target)
{
  if (containsKeyCaseInsensitive(object, key))
  {
    target = getStringOrThrow(requireKeyCaseInsensitive(object, key, ""), key, "");
  }
}

static void readOptionalBool(const nlohmann::json& object, const std::string& key, bool& target)
{
  if (containsKeyCaseInsensitive(object, key))
  {
    target = getBoolOrThrow(requireKeyCaseInsensitive(object, key, ""), key, "");
  }
}

template <typename T>
static T parseMappedStringOrThrow(const nlohmann::json& object, const std::string& key,
                                  std::initializer_list<std::pair<std::string_view, T>> mapping)
{
  const nlohmann::json* value = findKeyCaseInsensitive(object, key);
  if (value == nullptr)
  {
    throw std::runtime_error("Error: required key '" + key + "' missing");
  }

  std::string str = getStringOrThrow(*value, key, "");
  for (const auto& [name, result] : mapping)
  {
    if (caseInSensStringCompare(str, std::string{name}))
    {
      return result;
    }
  }

  throw std::runtime_error("Error: invalid " + key + " '" + str + "'");
}

template <typename T>
static void readOptionalMappedString(const nlohmann::json& object, const std::string& key, T& target,
                                     std::initializer_list<std::pair<std::string_view, T>> mapping)
{
  const nlohmann::json* value = findKeyCaseInsensitive(object, key);
  if (value == nullptr)
  {
    return;
  }

  std::string str = getStringOrThrow(*value, key, "");
  for (const auto& [name, result] : mapping)
  {
    if (caseInSensStringCompare(str, std::string{name}))
    {
      target = result;
      return;
    }
  }

  throw std::runtime_error("Error: invalid " + key + " '" + str + "'");
}

struct IsothermSpec
{
  Isotherm::Type type;
  std::size_t parameterCount;
};

static const IsothermSpec* findIsothermSpec(const std::string& typeString)
{
  static const std::vector<std::pair<std::string_view, IsothermSpec>> specs{
      {"Langmuir", {Isotherm::Type::Langmuir, 2}},
      {"Anti-Langmuir", {Isotherm::Type::Anti_Langmuir, 2}},
      {"Anti_Langmuir", {Isotherm::Type::Anti_Langmuir, 2}},
      {"BET", {Isotherm::Type::BET, 3}},
      {"Henry", {Isotherm::Type::Henry, 1}},
      {"Freundlich", {Isotherm::Type::Freundlich, 2}},
      {"Sips", {Isotherm::Type::Sips, 3}},
      {"Langmuir-Freundlich", {Isotherm::Type::Langmuir_Freundlich, 3}},
      {"Langmuir_Freundlich", {Isotherm::Type::Langmuir_Freundlich, 3}},
      {"Redlich-Peterson", {Isotherm::Type::Redlich_Peterson, 3}},
      {"Redlich_Peterson", {Isotherm::Type::Redlich_Peterson, 3}},
      {"Toth", {Isotherm::Type::Toth, 3}},
      {"Unilan", {Isotherm::Type::Unilan, 3}},
      {"O'Brian&Myers", {Isotherm::Type::OBrien_Myers, 3}},
      {"OBrien_Myers", {Isotherm::Type::OBrien_Myers, 3}},
      {"OBrien&Myers", {Isotherm::Type::OBrien_Myers, 3}},
      {"Quadratic", {Isotherm::Type::Quadratic, 3}},
      {"Temkin", {Isotherm::Type::Temkin, 3}},
      {"Bingel&Walton", {Isotherm::Type::BingelWalton, 3}},
      {"BingelWalton", {Isotherm::Type::BingelWalton, 3}},
  };

  for (const auto& [name, spec] : specs)
  {
    if (caseInSensStringCompare(typeString, std::string{name}))
    {
      return &spec;
    }
  }
  return nullptr;
}

static void addIsothermSiteFromJson(Component& component, const std::string& typeString, const nlohmann::json& params,
                                    const std::string& context)
{
  const IsothermSpec* spec = findIsothermSpec(typeString);
  if (spec == nullptr)
  {
    throw std::runtime_error("Error: unknown isotherm type '" + typeString + "'" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  std::vector<double> values = requireDoubleParameterCount(getNumberListOrThrow<double>(params, typeString, context),
                                                           spec->parameterCount, typeString, context);

  component.isotherm.add(Isotherm(spec->type, values, spec->parameterCount));
}

static void requireConfigured(bool condition, const std::string& message)
{
  if (!condition)
  {
    throw std::runtime_error(message);
  }
}

InputReader::InputReader(const std::string fileName) : components()
{
  components.reserve(16);

  if (!std::filesystem::exists(fileName))
  {
    throw std::runtime_error("Required input file '" + fileName + "' does not exist");
  }

  std::ifstream fileInput{fileName};
  if (!fileInput)
  {
    throw std::runtime_error("Required input file '" + fileName + "' could not be opened");
  }

  const nlohmann::json parsed_data = readJson(fileInput);

  // Track presence separately from value, because some numeric defaults are valid runtime values.
  const bool hasInletPressure = containsKeyCaseInsensitive(parsed_data, "InletPressure") ||
                                containsKeyCaseInsensitive(parsed_data, "TotalPressure");
  const bool hasOutletPressure = containsKeyCaseInsensitive(parsed_data, "OutletPressure");
  const bool hasColumnEntranceVelocity = containsKeyCaseInsensitive(parsed_data, "ColumnEntranceVelocity");
  const bool hasPressureGradient = containsKeyCaseInsensitive(parsed_data, "PressureGradient");

  readOptionalMappedString(parsed_data, "SimulationType", simulationType,
                           {{"Breakthrough", SimulationType::Breakthrough},
                            {"MixturePrediction", SimulationType::MixturePrediction},
                            {"Fitting", SimulationType::Fitting},
                            {"SwingAdsorption", SimulationType::SwingAdsorption},
                            {"Test", SimulationType::Test}});

  readOptionalMappedString(parsed_data, "MixturePredictionMethod", mixturePredictionMethod,
                           {{"IAST", 0}, {"SIAST", 1}, {"EI", 2}, {"SEI", 3}});

  readOptionalMappedString(parsed_data, "IASTMethod", IASTMethod, {{"FastIAS", 0}, {"Bisection", 1}});

  readOptionalMappedString(parsed_data, "BreakthroughIntegrator", breakthroughIntegrator,
                           {{"RungeKutta3", 0}, {"CVODE", 1}, {"SIRK3", 3}});

  readOptionalMappedString(parsed_data, "VelocityProfile", velocityProfile,
                           {{"FixedPressureGradient", 0}, {"Ergun", 1}, {"FixedVelocity", 2}});

  readOptionalMappedString(
      parsed_data, "BoundaryCondition", boundaryCondition,
      {{"InletPressureInletVelocity", 0}, {"InletPressureOutletPressure", 1}, {"InletVelocityOutletPressure", 2}});

  readOptionalMappedString(parsed_data, "PressureScale", pressureScale, {{"Log", 0}, {"Linear", 1}, {"Normal", 1}});

  readOptionalString(parsed_data, "ReadColumnFile", readColumnFile);
  readOptionalString(parsed_data, "DisplayName", displayName);

  readOptionalNumber<double>(parsed_data, "Temperature", temperature);
  readOptionalNumber<double>(parsed_data, "ColumnVoidFraction", columnVoidFraction);
  readOptionalNumber<double>(parsed_data, "DynamicViscosity", dynamicViscosity);
  readOptionalNumber<double>(parsed_data, "ParticleDiameter", particleDiameter);
  readOptionalNumber<double>(parsed_data, "ParticleDensity", particleDensity);

  if (containsKeyCaseInsensitive(parsed_data, "InletPressure"))
  {
    inletPressure =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "InletPressure", ""), "InletPressure", "");
  }
  else if (containsKeyCaseInsensitive(parsed_data, "TotalPressure"))
  {
    // Backward-compatible alias for old breakthrough input files.
    inletPressure =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "TotalPressure", ""), "TotalPressure", "");
  }

  readOptionalNumber<double>(parsed_data, "OutletPressure", outletPressure);
  readOptionalNumber<double>(parsed_data, "PressureStart", pressureStart);
  readOptionalNumber<double>(parsed_data, "PressureEnd", pressureEnd);
  readOptionalNumber<size_t>(parsed_data, "NumberOfPressurePoints", numberOfPressurePoints);
  readOptionalNumber<double>(parsed_data, "PressureGradient", pressureGradient);
  readOptionalNumber<double>(parsed_data, "ColumnEntranceVelocity", columnEntranceVelocity);
  readOptionalNumber<size_t>(parsed_data, "NumberOfInitTimeSteps", numberOfInitTimeSteps);
  readOptionalNumber<double>(parsed_data, "TimeStep", timeStep);
  readOptionalNumber<size_t>(parsed_data, "PrintEvery", printEvery);
  readOptionalNumber<size_t>(parsed_data, "WriteEvery", writeEvery);
  readOptionalNumber<double>(parsed_data, "ColumnLength", columnLength);
  readOptionalNumber<size_t>(parsed_data, "NumberOfGridPoints", numberOfGridPoints);
  readOptionalNumber<size_t>(parsed_data, "ColumnPressure", columnPressure);
  readOptionalNumber<size_t>(parsed_data, "ColumnLoading", columnLoading);
  readOptionalNumber<size_t>(parsed_data, "ColumnError", columnError);

  if (containsKeyCaseInsensitive(parsed_data, "NumberOfTimeSteps"))
  {
    const nlohmann::json& v = requireKeyCaseInsensitive(parsed_data, "NumberOfTimeSteps", "");
    if (v.is_string())
    {
      std::string s = v.get<std::string>();
      if (caseInSensStringCompare(s, "auto"))
      {
        autoNumberOfTimeSteps = true;
      }
      else
      {
        throw std::runtime_error("Error: NumberOfTimeSteps must be integer or 'auto'");
      }
    }
    else
    {
      numberOfTimeSteps = getNumberOrThrow<std::size_t>(v, "NumberOfTimeSteps", "");
      autoNumberOfTimeSteps = false;
    }
  }

  // InfluxTemperature defaults to Temperature when not specified.
  if (containsKeyCaseInsensitive(parsed_data, "InfluxTemperature"))
  {
    influxTemperature = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "InfluxTemperature", ""),
                                                 "InfluxTemperature", "");
  }
  else
  {
    influxTemperature = temperature;
  }

  readOptionalNumber<double>(parsed_data, "internalDiameter", internalDiameter);
  readOptionalNumber<double>(parsed_data, "outerDiameter", outerDiameter);
  readOptionalNumber<double>(parsed_data, "wallDensity", wallDensity);
  readOptionalNumber<double>(parsed_data, "gasThermalConductivity", gasThermalConductivity);
  readOptionalNumber<double>(parsed_data, "wallThermalConductivity", wallThermalConductivity);
  readOptionalNumber<double>(parsed_data, "heatTransferGasWall", heatTransferGasWall);
  readOptionalNumber<double>(parsed_data, "heatTransferGasSolid", heatTransferGasSolid);
  readOptionalNumber<double>(parsed_data, "heatTransferWallExternal", heatTransferWallExternal);
  readOptionalNumber<double>(parsed_data, "heatCapacityGas", heatCapacityGas);
  readOptionalNumber<double>(parsed_data, "heatCapacitySolid", heatCapacitySolid);
  readOptionalNumber<double>(parsed_data, "heatCapacityWall", heatCapacityWall);
  readOptionalBool(parsed_data, "energyBalance", energyBalance);

  if (containsKeyCaseInsensitive(parsed_data, "swingTemperatures"))
  {
    swingTemperatures = getNumberListOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "swingTemperatures", ""),
                                                     "swingTemperatures", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "swingPressures"))
  {
    swingPressures = getNumberListOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "swingPressures", ""),
                                                  "swingPressures", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "swingSteps"))
  {
    swingSteps =
        getNumberListOrThrow<size_t>(requireKeyCaseInsensitive(parsed_data, "swingSteps", ""), "swingSteps", "");
  }

  // Components
  if (containsKeyCaseInsensitive(parsed_data, "Components"))
  {
    const nlohmann::json& comps = requireKeyCaseInsensitive(parsed_data, "Components", "");
    auto parseComponentObject = [&](std::size_t componentId, const nlohmann::json& item)
    {
      if (!item.is_object())
      {
        throw std::runtime_error("Error: each component entry must be an object");
      }

      std::string context = "Component " + std::to_string(componentId);

      std::string componentName = getStringOrThrow(requireKeyCaseInsensitive(item, "Name", context), "Name", context);

      components.emplace_back(componentId, componentName);
      Component& comp = components.back();

      readOptionalString(item, "FileName", comp.filename);
      readOptionalBool(item, "CarrierGas", comp.isCarrierGas);
      readOptionalNumber<double>(item, "GasPhaseMolFraction", comp.Yi0);
      readOptionalNumber<double>(item, "MassTransferCoefficient", comp.Kl);
      readOptionalNumber<double>(item, "AxialDispersionCoefficient", comp.D);
      readOptionalNumber<double>(item, "MolecularWeight", comp.molecularWeight);
      readOptionalNumber<double>(item, "HeatOfAdsorption", comp.heatOfAdsorption);
      readOptionalNumber<size_t>(item, "NumberOfIsothermSites", comp.isotherm.numberOfSites);

      // Preferred format: "IsothermSites": [{"Type": "...", "Parameters": [...]}, ...].
      if (containsKeyCaseInsensitive(item, "IsothermSites"))
      {
        const nlohmann::json& sites = requireKeyCaseInsensitive(item, "IsothermSites", context);
        if (!sites.is_array())
        {
          throw std::runtime_error("Error: IsothermSites must be an array (" + context + ")");
        }

        for (std::size_t siteId = 0; siteId < sites.size(); ++siteId)
        {
          const nlohmann::json& site = sites[siteId];
          std::string siteContext = context + ", IsothermSite " + std::to_string(siteId);

          if (site.is_object() && containsKeyCaseInsensitive(site, "Type") &&
              containsKeyCaseInsensitive(site, "Parameters"))
          {
            std::string typeString =
                getStringOrThrow(requireKeyCaseInsensitive(site, "Type", siteContext), "Type", siteContext);
            const nlohmann::json& params = requireKeyCaseInsensitive(site, "Parameters", siteContext);
            addIsothermSiteFromJson(comp, typeString, params, siteContext);
          }
          else if (site.is_object() && site.size() == 1)
          {
            // Compact format: {"Langmuir": [...]}.
            auto it = site.begin();
            addIsothermSiteFromJson(comp, it.key(), it.value(), siteContext);
          }
          else
          {
            throw std::runtime_error("Error: invalid IsothermSites entry (" + siteContext + ")");
          }
        }
      }
      else
      {
        // Legacy format: isotherm type names appear directly as component keys.
        for (auto it = item.begin(); it != item.end(); ++it)
        {
          if (findIsothermSpec(it.key()) != nullptr)
          {
            addIsothermSiteFromJson(comp, it.key(), it.value(), context);
          }
        }
      }
    };

    components.clear();

    if (comps.is_array())
    {
      components.reserve(comps.size());
      for (std::size_t componentId = 0; componentId < comps.size(); ++componentId)
      {
        parseComponentObject(componentId, comps[componentId]);
      }
    }
    else if (comps.is_object())
    {
      // Object form is accepted, but component IDs follow JSON insertion order.
      components.reserve(comps.size());
      std::size_t componentId = 0;
      for (auto it = comps.begin(); it != comps.end(); ++it, ++componentId)
      {
        parseComponentObject(componentId, it.value());
      }
    }
    else
    {
      throw std::runtime_error("Error: 'Components' must be an array or object");
    }
  }

  // Normalize feed gas fractions for non-fitting runs.
  if (simulationType != SimulationType::Fitting)
  {
    double sum = 0.0;
    for (size_t j = 0; j < components.size(); ++j)
    {
      sum += components[j].Yi0;
    }
    if (std::abs(sum - 1.0) > 1e-15)
    {
      std::cout << "Normalizing: Gas-phase molfractions did not sum exactly to unity!\n\n";
      for (size_t j = 0; j < components.size(); ++j)
      {
        components[j].Yi0 /= sum;
      }
    }
  }

  numberOfCarrierGases = 0;
  carrierGasComponent = 0;
  for (size_t j = 0; j < components.size(); ++j)
  {
    if (components[j].isCarrierGas)
    {
      carrierGasComponent = j;
      std::vector<double> values{1.0, 0.0};
      Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, 2);
      components[carrierGasComponent].isotherm.add(isotherm);
      components[carrierGasComponent].isotherm.numberOfSites = 1;

      ++numberOfCarrierGases;
    }
  }

  if ((mixturePredictionMethod == 2) || (mixturePredictionMethod == 3))
  {
    for (size_t i = 0; i < components.size(); ++i)
    {
      for (size_t j = 0; j < components[i].isotherm.numberOfSites; ++j)
      {
        if (components[i].isotherm.sites[j].type != Isotherm::Type::Langmuir)
        {
          throw std::runtime_error("Error: Explicit mixture prediction must use single Langmuir isotherms");
        }
      }
    }
  }

  maxIsothermTerms = 0;
  if (!components.empty())
  {
    std::vector<Component>::iterator maxIsothermTermsIterator =
        std::max_element(components.begin(), components.end(), [](Component& lhs, Component& rhs)
                         { return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites; });
    maxIsothermTerms = maxIsothermTermsIterator->isotherm.numberOfSites;
  }

  if (simulationType == SimulationType::Breakthrough)
  {
    static constexpr size_t fixedPressureGradient = 0;
    static constexpr size_t ergun = 1;
    static constexpr size_t fixedVelocity = 2;

    static constexpr size_t inletPressureInletVelocity = 0;
    static constexpr size_t inletPressureOutletPressure = 1;
    static constexpr size_t inletVelocityOutletPressure = 2;

    // Boundary conditions define which externally supplied quantities are mandatory.
    const bool boundaryNeedsInletPressure =
        boundaryCondition == inletPressureInletVelocity || boundaryCondition == inletPressureOutletPressure;
    const bool boundaryNeedsOutletPressure =
        boundaryCondition == inletPressureOutletPressure || boundaryCondition == inletVelocityOutletPressure;
    const bool boundaryNeedsInletVelocity =
        boundaryCondition == inletPressureInletVelocity || boundaryCondition == inletVelocityOutletPressure;

    requireConfigured(numberOfCarrierGases != 0, "Error: no carrier gas component present");
    requireConfigured(numberOfCarrierGases == 1,
                      "Error: multiple carrier gas component present (there can be only one)");

    requireConfigured(temperature >= 0.0, "Error: temperature not set (Use e.g.: 'Temperature 300')");
    requireConfigured(columnVoidFraction > 0.0,
                      "Error: void-fraction of the column not set or invalid (Use e.g.: 'ColumnVoidFraction 0.4')");
    requireConfigured(particleDensity >= 0.0, "Error: particle density not set (Use e.g.: 'ParticleDensity 1408.2')");
    requireConfigured(numberOfTimeSteps != 0 || autoNumberOfTimeSteps,
                      "Error: number of time steps not set (Use e.g.: 'NumberOfTimeSteps 5000000')");
    requireConfigured(numberOfGridPoints != 0,
                      "Error: number of grid points not set (Use e.g.: 'NumberOfGridPoints 50')");
    requireConfigured(columnLength > 0.0, "Error: column length not set or invalid (Use e.g.: 'ColumnLength 0.3')");

    if (boundaryNeedsInletPressure)
    {
      requireConfigured(hasInletPressure && inletPressure >= 0.0,
                        "Error: inlet pressure not set (Use e.g.: 'InletPressure 1e5')");
    }

    if (boundaryNeedsOutletPressure)
    {
      requireConfigured(hasOutletPressure && outletPressure >= 0.0,
                        "Error: outlet pressure not set (Use e.g.: 'OutletPressure 1e5')");
    }

    if (boundaryNeedsInletVelocity || velocityProfile == fixedVelocity)
    {
      requireConfigured(hasColumnEntranceVelocity && columnEntranceVelocity >= 0.0,
                        "Error: inlet velocity not set (Use e.g.: 'ColumnEntranceVelocity 0.1')");
    }

    if (velocityProfile == fixedPressureGradient)
    {
      requireConfigured(hasPressureGradient,
                        "Error: pressure gradient not set (Use e.g.: 'PressureGradient -1e5')");

      requireConfigured(inletPressure + columnLength * pressureGradient >= 0.0,
                        "Error: pressure profile becomes negative: InletPressure + ColumnLength * PressureGradient < 0");

      if (boundaryNeedsOutletPressure)
      {
        throw std::runtime_error(
            "Error: FixedPressureGradient cannot be used with a boundary condition that sets OutletPressure");
      }
    }

    if (velocityProfile == ergun && hasPressureGradient)
    {
      std::cerr << "Warning: PressureGradient is ignored when VelocityProfile is Ergun" << std::endl;
    }
  }
}
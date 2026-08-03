#include "inputreader.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <print>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

#include "json.h"
#include "utils.h"

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
    std::print(stderr, "parse error at byte {}\n", ex.byte);
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

template <typename T>
static void readOptionalNumber(const nlohmann::json& object, const std::string& key, std::optional<T>& target)
{
  if (containsKeyCaseInsensitive(object, key))
  {
    target = getNumberOrThrow<T>(requireKeyCaseInsensitive(object, key, ""), key, "");
  }
}

static void readOptionalNonNegativeInteger(const nlohmann::json& object, const std::string& key, size_t& target)
{
  const nlohmann::json* value = findKeyCaseInsensitive(object, key);
  if (value == nullptr) return;

  try
  {
    if (value->is_number_unsigned())
    {
      target = value->get<size_t>();
      return;
    }
    if (value->is_number_integer())
    {
      const std::int64_t integer = value->get<std::int64_t>();
      if (integer >= 0)
      {
        target = static_cast<size_t>(integer);
        return;
      }
    }
  }
  catch (const nlohmann::json::exception&)
  {
  }

  throw std::runtime_error("Error: key '" + key + "' must be a non-negative integer");
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

static void requireOnlyKnownKeys(const nlohmann::json& object, std::initializer_list<std::string_view> allowedKeys,
                                 const std::string& context, bool allowIsothermTypes);

static void requireOnlyExactKeys(const nlohmann::json& object,
                                 std::initializer_list<std::string_view> allowedKeys,
                                 const std::string& context)
{
  for (auto it = object.begin(); it != object.end(); ++it)
  {
    const bool known = std::any_of(allowedKeys.begin(), allowedKeys.end(),
                                   [&](std::string_view key) { return it.key() == key; });
    if (!known)
    {
      throw std::runtime_error("Error: unknown key '" + it.key() + "' (" + context + ")");
    }
  }
}

static void requireExactKeys(const nlohmann::json& object,
                             std::initializer_list<std::string_view> requiredKeys,
                             const std::string& context)
{
  requireOnlyExactKeys(object, requiredKeys, context);
  for (const std::string_view key : requiredKeys)
  {
    if (!object.contains(std::string{key}))
    {
      throw std::runtime_error("Error: required key '" + std::string{key} + "' missing (" + context + ")");
    }
  }
}

static Isotherm readIsothermModel(const nlohmann::json& value, bool nonIsothermal,
                                  const std::string& context);

static Chemisorption::Type parseChemisorptionType(const std::string& typeString)
{
  if (caseInSensStringCompare(typeString, "None")) return Chemisorption::Type::None;
  if (caseInSensStringCompare(typeString, "FirstOrder")) return Chemisorption::Type::FirstOrder;
  if (caseInSensStringCompare(typeString, "PseudoNth")) return Chemisorption::Type::PseudoNth;
  if (caseInSensStringCompare(typeString, "Avrami")) return Chemisorption::Type::Avrami;
  if (caseInSensStringCompare(typeString, "General")) return Chemisorption::Type::General;
  if (caseInSensStringCompare(typeString, "Elovich")) return Chemisorption::Type::Elovich;

  throw std::runtime_error("Error: invalid Chemisorption Type '" + typeString + "'");
}

static Chemisorption readChemisorptionSite(const nlohmann::json& value, bool nonIsothermal,
                                           const std::string& context)
{
  if (!value.is_object())
  {
    throw std::runtime_error("Error: Chemisorption site must be an object (" + context + ")");
  }

  requireOnlyExactKeys(value, {"Type", "Parameters"}, context);

  Chemisorption chemisorption;
  chemisorption.type = parseChemisorptionType(
      getStringOrThrow(requireKeyCaseInsensitive(value, "Type", context), "Type", context));

  const nlohmann::json& params = requireKeyCaseInsensitive(value, "Parameters", context);
  if (!params.is_object())
  {
    throw std::runtime_error("Error: Chemisorption Parameters must be an object (" + context + ")");
  }

  const std::string parametersContext = context + " Chemisorption Parameters";
  switch (chemisorption.type)
  {
    case Chemisorption::Type::None:
      throw std::runtime_error("Error: None is not a valid ChemisorptionSites model (" + context + ")");
    case Chemisorption::Type::FirstOrder:
      requireExactKeys(params,
                       {"rateCoefficient", "maximumLoading", "heatOfChemisorption", "Isotherm"},
                       parametersContext);
      break;
    case Chemisorption::Type::PseudoNth:
    case Chemisorption::Type::Avrami:
      requireExactKeys(params,
                       {"rateCoefficient", "order", "maximumLoading", "heatOfChemisorption", "Isotherm"},
                       parametersContext);
      break;
    case Chemisorption::Type::General:
      requireExactKeys(params,
                       {"maximumLoading", "heatOfChemisorption", "adsorptionRateCoefficient",
                        "adsorptionActivationEnergy", "desorptionRateCoefficient",
                        "desorptionActivationEnergy", "poreConcentrationOrder", "capacityOrder",
                        "desorptionOrder", "filmMassTransferCoefficient", "poreDiffusivity",
                        "usePoreSurfaceTransport", "Isotherm"},
                       parametersContext);
      break;
    case Chemisorption::Type::Elovich:
      requireExactKeys(params,
                       {"maximumLoading", "heatOfChemisorption", "alpha", "beta",
                        "filmMassTransferCoefficient", "poreDiffusivity", "usePoreSurfaceTransport",
                        "Isotherm"},
                       parametersContext);
      break;
  }

  readOptionalNumber<double>(params, "rateCoefficient", chemisorption.rateCoefficient);
  readOptionalNonNegativeInteger(params, "order", chemisorption.order);
  readOptionalNumber<double>(params, "maximumLoading", chemisorption.maximumLoading);
  readOptionalNumber<double>(params, "heatOfChemisorption", chemisorption.heatOfChemisorption);
  readOptionalNumber<double>(params, "adsorptionRateCoefficient", chemisorption.adsorptionRateCoefficient);
  readOptionalNumber<double>(params, "adsorptionActivationEnergy", chemisorption.adsorptionActivationEnergy);
  readOptionalNumber<double>(params, "desorptionRateCoefficient", chemisorption.desorptionRateCoefficient);
  readOptionalNumber<double>(params, "desorptionActivationEnergy", chemisorption.desorptionActivationEnergy);
  readOptionalNonNegativeInteger(params, "poreConcentrationOrder", chemisorption.poreConcentrationOrder);
  readOptionalNonNegativeInteger(params, "capacityOrder", chemisorption.capacityOrder);
  readOptionalNonNegativeInteger(params, "desorptionOrder", chemisorption.desorptionOrder);
  readOptionalNumber<double>(params, "alpha", chemisorption.elovichAlpha);
  readOptionalNumber<double>(params, "beta", chemisorption.elovichBeta);
  readOptionalNumber<double>(params, "filmMassTransferCoefficient", chemisorption.filmMassTransferCoefficient);
  readOptionalNumber<double>(params, "poreDiffusivity", chemisorption.poreDiffusivity);
  readOptionalBool(params, "usePoreSurfaceTransport", chemisorption.usePoreSurfaceTransport);

  chemisorption.isotherm = readIsothermModel(
      requireKeyCaseInsensitive(params, "Isotherm", context), nonIsothermal, context + " Isotherm");

  if (chemisorption.type == Chemisorption::Type::General)
  {
    if (chemisorption.adsorptionRateCoefficient < 0.0 || chemisorption.desorptionRateCoefficient < 0.0 ||
        chemisorption.filmMassTransferCoefficient < 0.0 || chemisorption.poreDiffusivity < 0.0)
    {
      throw std::runtime_error("Error: General Chemisorption kinetic and transport coefficients must be non-negative (" +
                               context + ")");
    }
  }
  else if (chemisorption.type == Chemisorption::Type::Elovich)
  {
    if (chemisorption.elovichAlpha < 0.0 || chemisorption.elovichBeta < 0.0 ||
        chemisorption.filmMassTransferCoefficient < 0.0 || chemisorption.poreDiffusivity < 0.0)
    {
      throw std::runtime_error("Error: Elovich Chemisorption kinetic and transport coefficients must be non-negative (" +
                               context + ")");
    }
  }
  else if (chemisorption.rateCoefficient < 0.0)
  {
    throw std::runtime_error("Error: Chemisorption rateCoefficient must be non-negative (" + context + ")");
  }
  if (chemisorption.maximumLoading <= 0.0)
  {
    throw std::runtime_error("Error: Chemisorption maximumLoading must be positive (" + context + ")");
  }
  return chemisorption;
}

static void readChemisorption(MultiSiteChemisorption& chemisorption, const nlohmann::json& object,
                              bool nonIsothermal, const std::string& context)
{
  const nlohmann::json* sites = findKeyCaseInsensitive(object, "ChemisorptionSites");

  if (sites != nullptr)
  {
    if (!sites->is_array())
    {
      throw std::runtime_error("Error: ChemisorptionSites must be an array (" + context + ")");
    }
    chemisorption = MultiSiteChemisorption{};
    for (size_t site = 0; site < sites->size(); ++site)
    {
      chemisorption.add(
          readChemisorptionSite((*sites)[site], nonIsothermal,
                                context + ", ChemisorptionSite " + std::to_string(site)));
    }
  }

  if (containsKeyCaseInsensitive(object, "NumberOfChemisorptionSites"))
  {
    size_t declaredSites = 0;
    readOptionalNonNegativeInteger(object, "NumberOfChemisorptionSites", declaredSites);
    if (declaredSites != chemisorption.numberOfSites)
    {
      throw std::runtime_error("Error: NumberOfChemisorptionSites does not match ChemisorptionSites (" + context + ")");
    }
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
  bool nonIsothermalImplemented;
};

static const IsothermSpec* findIsothermSpec(const std::string& typeString)
{
  static const std::vector<std::pair<std::string_view, IsothermSpec>> specs{
      {"Langmuir", {Isotherm::Type::Langmuir, 2, true}},
      {"Anti-Langmuir", {Isotherm::Type::Anti_Langmuir, 2, false}},
      {"Anti_Langmuir", {Isotherm::Type::Anti_Langmuir, 2, false}},
      {"BET", {Isotherm::Type::BET, 3, false}},
      {"Henry", {Isotherm::Type::Henry, 1, false}},
      {"Freundlich", {Isotherm::Type::Freundlich, 2, false}},
      {"Sips", {Isotherm::Type::Sips, 3, true}},
      {"Langmuir-Freundlich", {Isotherm::Type::Langmuir_Freundlich, 3, true}},
      {"Langmuir_Freundlich", {Isotherm::Type::Langmuir_Freundlich, 3, true}},
      {"Redlich-Peterson", {Isotherm::Type::Redlich_Peterson, 3, false}},
      {"Redlich_Peterson", {Isotherm::Type::Redlich_Peterson, 3, false}},
      {"Toth", {Isotherm::Type::Toth, 3, false}},
      {"Unilan", {Isotherm::Type::Unilan, 3, false}},
      {"O'Brian&Myers", {Isotherm::Type::OBrien_Myers, 3, false}},
      {"OBrien_Myers", {Isotherm::Type::OBrien_Myers, 3, false}},
      {"OBrien&Myers", {Isotherm::Type::OBrien_Myers, 3, false}},
      {"Quadratic", {Isotherm::Type::Quadratic, 3, false}},
      {"Temkin", {Isotherm::Type::Temkin, 3, false}},
      {"Bingel&Walton", {Isotherm::Type::BingelWalton, 3, false}},
      {"BingelWalton", {Isotherm::Type::BingelWalton, 3, false}},
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

static bool matchesAllowedKey(std::string_view key, std::initializer_list<std::string_view> allowedKeys)
{
  for (const std::string_view allowedKey : allowedKeys)
  {
    if (caseInSensStringCompare(std::string{key}, std::string{allowedKey}))
    {
      return true;
    }
  }

  return false;
}

static void requireOnlyKnownKeys(const nlohmann::json& object, std::initializer_list<std::string_view> allowedKeys,
                                 const std::string& context, bool allowIsothermTypeKeys = false)
{
  if (!object.is_object())
  {
    throw std::runtime_error("Error: input settings must be a JSON object" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  for (auto it = object.begin(); it != object.end(); ++it)
  {
    if (matchesAllowedKey(it.key(), allowedKeys))
    {
      continue;
    }

    if (allowIsothermTypeKeys && findIsothermSpec(it.key()) != nullptr)
    {
      continue;
    }

    throw std::runtime_error("Error: unknown input key '" + it.key() + "'" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }
}

static void requireOnlyKnownGeneralSettings(const nlohmann::json& object)
{
  requireOnlyKnownKeys(object,
                       {"SimulationType",
                        "MixturePredictionMethod",
                        "IASTMethod",
                        "BreakthroughIntegrator",
                        "BoundaryCondition",
                        "VelocityProfile",
                        "PressureScale",
                        "ReadColumnFile",
                        "DisplayName",
                        "Temperature",
                        "ColumnVoidFraction",
                        "DynamicViscosity",
                        "ParticleDiameter",
                        "ParticleDensity",
                        "InletPressure",
                        "TotalPressure",
                        "OutletPressure",
                        "PressureStart",
                        "PressureEnd",
                        "NumberOfPressurePoints",
                        "PressureGradient",
                        "ColumnEntranceVelocity",
                        "NumberOfInitTimeSteps",
                        "TimeStep",
                        "PrintEvery",
                        "WriteEvery",
                        "ColumnLength",
                        "ColumnDistances",
                        "NumberOfGridPoints",
                        "ColumnPressure",
                        "ColumnLoading",
                        "ColumnError",
                        "NumberOfTimeSteps",
                        "Geometry",
                        "InfluxTemperature",
                        "internalDiameter",
                        "outerDiameter",
                        "wallDensity",
                        "gasThermalConductivity",
                        "wallThermalConductivity",
                        "heatTransferGasWall",
                        "heatTransferGasSolid",
                        "heatTransferWallExternal",
                        "heatCapacityGas",
                        "heatCapacitySolid",
                        "heatCapacityWall",
                        "energyBalance",
                        "swingTemperatures",
                        "swingPressures",
                        "swingSteps",
                        "Interfaces",
                        "ColumnSections",
                        "Adsorbents",
                        "Components"},
                       "general settings");
}

static void requireOnlyKnownComponentSettings(const nlohmann::json& object, const std::string& context)
{
  requireOnlyKnownKeys(object,
                       {"Name",
                        "FileName",
                        "CarrierGas",
                        "GasPhaseMolFraction",
                        "MassTransferCoefficient",
                        "AxialDispersionCoefficient",
                        "MolecularWeight",
                        "HeatOfAdsorption",
                        "NumberOfPhysisorptionSites",
                        "referenceTemperature",
                        "nonIsothermal",
                        "NumberOfChemisorptionSites",
                        "ChemisorptionSites",
                        "PhysisorptionSites"},
                       context);
}

static void requireOnlyKnownIsothermSiteSettings(const nlohmann::json& object, const std::string& context)
{
  requireOnlyKnownKeys(object, {"Type", "Parameters"}, context);
}

static void requireOnlyKnownAdsorbentSettings(const nlohmann::json& object, const std::string& context)
{
  requireOnlyKnownKeys(object,
                       {"Name",
                        "ParticleDiameter",
                        "ColumnVoidFraction",
                        "ParticleDensity",
                        "AdsorbentLength",
                        "ComponentParameters",
                        "Components"},
                       context);
}

static void requireOnlyKnownColumnSectionSettings(const nlohmann::json& object, const std::string& context)
{
  requireOnlyKnownKeys(object, {"Adsorbent", "Length", "InterfaceLength", "NumberOfGridPoints"}, context);
}

static void requireOnlyKnownAdsorbentComponentSettings(const nlohmann::json& object, const std::string& context)
{
  requireOnlyKnownKeys(object,
                       {"Name",
                        "MassTransferCoefficient",
                        "AxialDispersionCoefficient",
                        "HeatOfAdsorption",
                        "NumberOfPhysisorptionSites",
                        "referenceTemperature",
                        "nonIsothermal",
                        "NumberOfChemisorptionSites",
                        "ChemisorptionSites",
                        "PhysisorptionSites"},
                       context);
}

static Isotherm readIsothermModel(const nlohmann::json& value, bool nonIsothermal,
                                  const std::string& context)
{
  if (!value.is_object())
  {
    throw std::runtime_error("Error: Isotherm must be an object (" + context + ")");
  }
  requireOnlyExactKeys(value, {"Type", "Parameters"}, context);

  const std::string typeString =
      getStringOrThrow(requireKeyCaseInsensitive(value, "Type", context), "Type", context);
  const IsothermSpec* spec = findIsothermSpec(typeString);
  if (spec == nullptr)
  {
    throw std::runtime_error("Error: unknown isotherm type '" + typeString + "' (" + context + ")");
  }
  if (nonIsothermal && !spec->nonIsothermalImplemented)
  {
    throw std::logic_error("Error: nonIsothermal not implemented for " + typeString);
  }

  const nlohmann::json& parameters = requireKeyCaseInsensitive(value, "Parameters", context);
  std::vector<double> values = requireDoubleParameterCount(
      getNumberListOrThrow<double>(parameters, typeString, context), spec->parameterCount,
      typeString, context);
  return Isotherm(spec->type, values, nonIsothermal);
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

  if (component.nonIsothermal)
  {
    if (!spec->nonIsothermalImplemented)
    {
      throw std::logic_error("Error: nonIsothermal not implemented for " + typeString);
    }
  }

  component.isotherm.add(Isotherm(spec->type, values, component.nonIsothermal));
  component.isotherm.numberOfSites = component.isotherm.sites.size();
}

static void readPhysisorptionSites(Component& comp, const nlohmann::json& item,
                                   const std::string& context)
{
  if (!containsKeyCaseInsensitive(item, "PhysisorptionSites")) return;

  const nlohmann::json& sites = requireKeyCaseInsensitive(item, "PhysisorptionSites", context);
  if (!sites.is_array())
  {
    throw std::runtime_error("Error: PhysisorptionSites must be an array (" + context + ")");
  }

  for (std::size_t siteId = 0; siteId < sites.size(); ++siteId)
  {
    const nlohmann::json& site = sites[siteId];
    std::string siteContext = context + ", PhysisorptionSite " + std::to_string(siteId);

    if (site.is_object() && containsKeyCaseInsensitive(site, "Type") &&
        containsKeyCaseInsensitive(site, "Parameters"))
    {
      requireOnlyKnownIsothermSiteSettings(site, siteContext);

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
      throw std::runtime_error("Error: invalid PhysisorptionSites entry (" + siteContext + ")");
    }
  }
}

static Component parseComponentObject(std::size_t componentId, const nlohmann::json& item)
{
  if (!item.is_object())
  {
    throw std::runtime_error("Error: each component entry must be an object");
  }

  std::string context = "Component " + std::to_string(componentId);
  requireOnlyKnownComponentSettings(item, context);

  std::string componentName = getStringOrThrow(requireKeyCaseInsensitive(item, "Name", context), "Name", context);

  Component comp(componentId, componentName);

  readOptionalString(item, "FileName", comp.filename);
  readOptionalBool(item, "CarrierGas", comp.isCarrierGas);
  readOptionalNumber<double>(item, "GasPhaseMolFraction", comp.initialGasMoleFraction);
  readOptionalNumber<double>(item, "MassTransferCoefficient", comp.massTransferCoefficient);
  readOptionalNumber<double>(item, "AxialDispersionCoefficient", comp.axialDispersionCoefficient);
  readOptionalNumber<double>(item, "MolecularWeight", comp.molecularWeight);
  readOptionalNumber<double>(item, "HeatOfAdsorption", comp.heatOfAdsorption);
  readOptionalNumber<double>(item, "referenceTemperature", comp.referenceTemperature);
  readOptionalBool(item, "nonIsothermal", comp.nonIsothermal);
  readChemisorption(comp.chemisorption, item, comp.nonIsothermal, context);

  if (comp.nonIsothermal)
  {
    auto tmp = requireKeyCaseInsensitive(item, "HeatOfAdsorption",
                                         "Heat of Adsorption should be set for non-isothermal adsorption models.");
  }

  readPhysisorptionSites(comp, item, context);
  if (containsKeyCaseInsensitive(item, "NumberOfPhysisorptionSites"))
  {
    size_t declaredSites = 0;
    readOptionalNonNegativeInteger(item, "NumberOfPhysisorptionSites", declaredSites);
    if (declaredSites != comp.isotherm.numberOfSites)
    {
      throw std::runtime_error(
          "Error: NumberOfPhysisorptionSites does not match PhysisorptionSites (" + context + ")");
    }
  }
  return comp;
}

static void applyAdsorbentComponentParameters(Component& comp, const nlohmann::json& params,
                                              const std::string& context)
{
  if (!params.is_object())
  {
    throw std::runtime_error("Error: component parameters must be an object (" + context + ")");
  }

  requireOnlyKnownAdsorbentComponentSettings(params, context);

  readOptionalNumber<double>(params, "MassTransferCoefficient", comp.massTransferCoefficient);
  readOptionalNumber<double>(params, "AxialDispersionCoefficient", comp.axialDispersionCoefficient);
  readOptionalNumber<double>(params, "HeatOfAdsorption", comp.heatOfAdsorption);
  readOptionalNumber<double>(params, "referenceTemperature", comp.referenceTemperature);
  readOptionalBool(params, "nonIsothermal", comp.nonIsothermal);
  readChemisorption(comp.chemisorption, params, comp.nonIsothermal, context);

  const bool hasIsothermOverride = containsKeyCaseInsensitive(params, "PhysisorptionSites");

  if (hasIsothermOverride)
  {
    comp.isotherm = MultiSiteIsotherm{};
  }
  readPhysisorptionSites(comp, params, context);
  if (containsKeyCaseInsensitive(params, "NumberOfPhysisorptionSites"))
  {
    size_t declaredSites = 0;
    readOptionalNonNegativeInteger(params, "NumberOfPhysisorptionSites", declaredSites);
    if (declaredSites != comp.isotherm.numberOfSites)
    {
      throw std::runtime_error(
          "Error: NumberOfPhysisorptionSites does not match PhysisorptionSites (" + context + ")");
    }
  }
}

static void requireConfigured(bool condition, const std::string& message)
{
  if (!condition)
  {
    throw std::runtime_error(message);
  }
}

static const nlohmann::json* findFirstKeyCaseInsensitive(
    const nlohmann::json& object, std::initializer_list<std::string_view> keys)
{
  for (std::string_view key : keys)
  {
    const nlohmann::json* value = findKeyCaseInsensitive(object, std::string{key});
    if (value != nullptr)
    {
      return value;
    }
  }
  return nullptr;
}

static double geometryDoubleOrDefault(const nlohmann::json& object,
                                      std::initializer_list<std::string_view> keys, double defaultValue)
{
  const nlohmann::json* value = findFirstKeyCaseInsensitive(object, keys);
  return value == nullptr ? defaultValue : getNumberOrThrow<double>(*value, std::string{*keys.begin()}, "Geometry");
}

static double geometryRequiredDouble(const nlohmann::json& object,
                                     std::initializer_list<std::string_view> keys)
{
  const nlohmann::json* value = findFirstKeyCaseInsensitive(object, keys);
  if (value == nullptr)
  {
    throw std::runtime_error("Error: required Geometry key '" + std::string{*keys.begin()} + "' missing");
  }
  return getNumberOrThrow<double>(*value, std::string{*keys.begin()}, "Geometry");
}

static std::size_t geometryRequiredSize(const nlohmann::json& object,
                                        std::initializer_list<std::string_view> keys)
{
  const nlohmann::json* value = findFirstKeyCaseInsensitive(object, keys);
  if (value == nullptr)
  {
    throw std::runtime_error("Error: required Geometry key '" + std::string{*keys.begin()} + "' missing");
  }

  try
  {
    if (value->is_number_unsigned())
    {
      return value->get<std::size_t>();
    }
    if (value->is_number_integer())
    {
      const std::int64_t integer = value->get<std::int64_t>();
      if (integer > 0)
      {
        return static_cast<std::size_t>(integer);
      }
    }
  }
  catch (const nlohmann::json::exception&)
  {
  }

  throw std::runtime_error("Error: Geometry key '" + std::string{*keys.begin()} +
                           "' must be a positive integer");
}

static std::string geometryStringOrDefault(const nlohmann::json& object,
                                           std::initializer_list<std::string_view> keys,
                                           const std::string& defaultValue)
{
  const nlohmann::json* value = findFirstKeyCaseInsensitive(object, keys);
  return value == nullptr ? defaultValue : getStringOrThrow(*value, std::string{*keys.begin()}, "Geometry");
}

static Geometry readGeometry(const nlohmann::json& parsedData, double columnVoidFraction,
                             double particleDiameter, double internalDiameter, double outerDiameter)
{
  const nlohmann::json* geometryJson = findKeyCaseInsensitive(parsedData, "Geometry");
  if (geometryJson == nullptr)
  {
    return Geometry{HollowTube{columnVoidFraction, particleDiameter, internalDiameter, outerDiameter}};
  }
  if (!geometryJson->is_object())
  {
    throw std::runtime_error("Error: Geometry must be a JSON object");
  }

  requireOnlyKnownKeys(
      *geometryJson,
      {"Type",
       "GeometryType",
       "Kind",
       "ColumnVoidFraction",
       "VoidFraction",
       "ParticleDiameter",
       "InternalDiameter",
       "internalDiameter",
       "OuterDiameter",
       "outerDiameter",
       "ChannelShape",
       "MONOLITH_CHANNEL_SHAPE",
       "channel_shape",
       "InternalChannelDimension",
       "ChannelDimension",
       "ChannelDiameter",
       "d_int",
       "dint",
       "d_out",
       "dout",
       "NumberOfChannels",
       "NChannels",
       "N_channels",
       "N_MONOLITH_CHANNELS",
       "WashcoatThickness",
       "washcoat_thickness",
       "t_wc",
       "WashcoatVolumePerChannelVolume",
       "WASHCOAT_VOLUME_PER_CHANNEL_VOLUME"},
      "Geometry");

  const std::string geometryType =
      geometryStringOrDefault(*geometryJson, {"Type", "GeometryType", "Kind"}, "HollowTube");

  if (caseInSensStringCompare(geometryType, "HollowTube") ||
      caseInSensStringCompare(geometryType, "PackedBed") ||
      caseInSensStringCompare(geometryType, "PackedColumn"))
  {
    const double eps = geometryDoubleOrDefault(
        *geometryJson, {"ColumnVoidFraction", "VoidFraction"}, columnVoidFraction);
    const double dp = geometryDoubleOrDefault(
        *geometryJson, {"ParticleDiameter"}, particleDiameter);
    const double di = geometryDoubleOrDefault(
        *geometryJson, {"InternalDiameter", "internalDiameter"}, internalDiameter);
    const double douter = geometryDoubleOrDefault(
        *geometryJson, {"OuterDiameter", "outerDiameter"}, outerDiameter);
    return Geometry{HollowTube{eps, dp, di, douter}};
  }

  if (caseInSensStringCompare(geometryType, "Monolith"))
  {
    const std::string channelShape =
        geometryStringOrDefault(*geometryJson, {"ChannelShape", "MONOLITH_CHANNEL_SHAPE", "channel_shape"},
                                "circular");
    const double dInt = geometryRequiredDouble(
        *geometryJson, {"InternalChannelDimension", "ChannelDimension", "ChannelDiameter", "d_int", "dint"});
    const double dOut = geometryRequiredDouble(
        *geometryJson, {"OuterDiameter", "outerDiameter", "d_out", "dout"});
    const std::size_t numberOfChannels = geometryRequiredSize(
        *geometryJson, {"NumberOfChannels", "NChannels", "N_channels", "N_MONOLITH_CHANNELS"});
    const double washcoatThickness = geometryDoubleOrDefault(
        *geometryJson, {"WashcoatThickness", "washcoat_thickness", "t_wc"}, 0.0);
    const double washcoatVolume = geometryDoubleOrDefault(
        *geometryJson,
        {"WashcoatVolumePerChannelVolume", "WASHCOAT_VOLUME_PER_CHANNEL_VOLUME"},
        -1.0);
    return Geometry{Monolith{Monolith::parseChannelShape(channelShape), dInt, dOut,
                             numberOfChannels, washcoatThickness, washcoatVolume}};
  }

  throw std::runtime_error("Error: Geometry Type must be HollowTube or Monolith");
}

static void applyGeometryLegacyFields(const Geometry& geometry, double& columnVoidFraction,
                                      double& particleDiameter, double& internalDiameter,
                                      double& outerDiameter)
{
  const ShapeParameters& shape = geometry.shapeParameters();
  columnVoidFraction = shape.voidFraction;
  if (shape.solidCharacteristicLength > 0.0)
  {
    particleDiameter = shape.solidCharacteristicLength;
  }
  internalDiameter = shape.internalDiameter;
  outerDiameter = shape.outerDiameter;
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
  requireOnlyKnownGeneralSettings(parsed_data);

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

  readOptionalMappedString(parsed_data, "BoundaryCondition", boundaryCondition,
                           {{"InletPressureInletVelocity", 0},
                            {"InletVelocityInletPressure", 0},
                            {"InletPressureOutletPressure", 1},
                            {"OutletPressureInletPressure", 1},
                            {"InletVelocityOutletPressure", 2},
                            {"OutletPressureInletVelocity", 2},
                            {"FixedVelocity", 3},
                            {"FixedPressureInletVelocity", 4}});

  size_t legacyVelocityProfile = 1;
  readOptionalMappedString(parsed_data, "VelocityProfile", legacyVelocityProfile,
                           {{"FixedPressureGradient", 1}, {"Ergun", 1}, {"FixedVelocity", 3}});
  if (legacyVelocityProfile == 3)
  {
    boundaryCondition = legacyVelocityProfile;
  }

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
  if (containsKeyCaseInsensitive(parsed_data, "ColumnDistances"))
  {
    columnDistances =
        getNumberListOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ColumnDistances", ""), "ColumnDistances", "");
    if (columnDistances.size() < 2)
    {
      throw std::runtime_error("Error: ColumnDistances must contain at least two node positions");
    }
  }
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

  geometry = readGeometry(parsed_data, columnVoidFraction, particleDiameter, internalDiameter, outerDiameter);
  applyGeometryLegacyFields(geometry, columnVoidFraction, particleDiameter, internalDiameter, outerDiameter);

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

    components.clear();

    if (comps.is_array())
    {
      components.reserve(comps.size());
      for (std::size_t componentId = 0; componentId < comps.size(); ++componentId)
      {
        components.push_back(parseComponentObject(componentId, comps[componentId]));
      }
    }
    else if (comps.is_object())
    {
      // Object form is accepted, but component IDs follow JSON insertion order.
      components.reserve(comps.size());
      std::size_t componentId = 0;
      for (auto it = comps.begin(); it != comps.end(); ++it, ++componentId)
      {
        components.push_back(parseComponentObject(componentId, it.value()));
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
      sum += components[j].initialGasMoleFraction;
    }
    if (std::abs(sum - 1.0) > 1e-15)
    {
      std::print("Normalizing: Gas-phase molfractions did not sum exactly to unity!\n\n");
      for (size_t j = 0; j < components.size(); ++j)
      {
        components[j].initialGasMoleFraction /= sum;
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
      Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, false);
      components[carrierGasComponent].isotherm.add(isotherm);
      components[carrierGasComponent].isotherm.numberOfSites = 1;

      ++numberOfCarrierGases;
    }
  }

  adsorbentComponents.clear();
  adsorbentLengths.clear();
  adsorbentInterfaceLengths.clear();
  adsorbentGridPoints.clear();
  adsorbentVoidFractions.clear();
  adsorbentParticleDensities.clear();
  adsorbentParticleDiameters.clear();

  if (containsKeyCaseInsensitive(parsed_data, "Adsorbents"))
  {
    const nlohmann::json& adsorbentsJson = requireKeyCaseInsensitive(parsed_data, "Adsorbents", "");
    if (!adsorbentsJson.is_array() || adsorbentsJson.empty())
    {
      throw std::runtime_error("Error: 'Adsorbents' must be a non-empty array");
    }

    std::unordered_map<std::string, std::size_t> componentIndexByName;
    for (std::size_t comp = 0; comp < components.size(); ++comp)
    {
      componentIndexByName.emplace(components[comp].name, comp);
    }

    std::vector<std::string> adsorbentNames;
    adsorbentNames.reserve(adsorbentsJson.size());
    std::unordered_map<std::string, std::size_t> adsorbentIndexByName;

    for (std::size_t ads = 0; ads < adsorbentsJson.size(); ++ads)
    {
      const nlohmann::json& adsorbent = adsorbentsJson[ads];
      const std::string context = "Adsorbent " + std::to_string(ads);
      if (!adsorbent.is_object())
      {
        throw std::runtime_error("Error: each adsorbent entry must be an object (" + context + ")");
      }
      requireOnlyKnownAdsorbentSettings(adsorbent, context);

      std::string adsorbentName = std::to_string(ads);
      readOptionalString(adsorbent, "Name", adsorbentName);
      if (!adsorbentIndexByName.emplace(adsorbentName, ads).second)
      {
        throw std::runtime_error("Error: duplicate adsorbent name '" + adsorbentName + "'");
      }
      adsorbentNames.push_back(adsorbentName);

      double adsorbentParticleDiameter = particleDiameter;
      double adsorbentVoidFraction = columnVoidFraction;
      double adsorbentParticleDensity = particleDensity;
      double adsorbentLength = -1.0;
      readOptionalNumber<double>(adsorbent, "ParticleDiameter", adsorbentParticleDiameter);
      readOptionalNumber<double>(adsorbent, "ColumnVoidFraction", adsorbentVoidFraction);
      readOptionalNumber<double>(adsorbent, "ParticleDensity", adsorbentParticleDensity);
      readOptionalNumber<double>(adsorbent, "AdsorbentLength", adsorbentLength);

      adsorbentParticleDiameters.push_back(adsorbentParticleDiameter);
      adsorbentVoidFractions.push_back(adsorbentVoidFraction);
      adsorbentParticleDensities.push_back(adsorbentParticleDensity);
      adsorbentLengths.push_back(adsorbentLength);

      std::vector<Component> componentsForAdsorbent = components;

      if (containsKeyCaseInsensitive(adsorbent, "ComponentParameters"))
      {
        const nlohmann::json& componentParameters =
            requireKeyCaseInsensitive(adsorbent, "ComponentParameters", context);
        if (!componentParameters.is_object())
        {
          throw std::runtime_error("Error: ComponentParameters must be an object (" + context + ")");
        }

        for (auto it = componentParameters.begin(); it != componentParameters.end(); ++it)
        {
          auto componentIndex = componentIndexByName.find(it.key());
          if (componentIndex == componentIndexByName.end())
          {
            throw std::runtime_error("Error: unknown component '" + it.key() + "' in ComponentParameters (" +
                                     context + ")");
          }

          applyAdsorbentComponentParameters(componentsForAdsorbent[componentIndex->second], it.value(),
                                            context + ", ComponentParameters " + it.key());
        }
      }

      if (containsKeyCaseInsensitive(adsorbent, "Components"))
      {
        const nlohmann::json& adsorbentComponentsJson = requireKeyCaseInsensitive(adsorbent, "Components", context);
        if (!adsorbentComponentsJson.is_array())
        {
          throw std::runtime_error("Error: Adsorbent Components must be an array (" + context + ")");
        }

        for (std::size_t componentId = 0; componentId < adsorbentComponentsJson.size(); ++componentId)
        {
          Component adsorbentComponent = parseComponentObject(componentId, adsorbentComponentsJson[componentId]);
          auto componentIndex = componentIndexByName.find(adsorbentComponent.name);
          if (componentIndex == componentIndexByName.end())
          {
            throw std::runtime_error("Error: unknown component '" + adsorbentComponent.name +
                                     "' in Adsorbent Components (" + context + ")");
          }
          adsorbentComponent.id = componentsForAdsorbent[componentIndex->second].id;
          componentsForAdsorbent[componentIndex->second] = adsorbentComponent;
        }
      }

      adsorbentComponents.push_back(std::move(componentsForAdsorbent));
    }

    if (containsKeyCaseInsensitive(parsed_data, "ColumnSections"))
    {
      const nlohmann::json& sections = requireKeyCaseInsensitive(parsed_data, "ColumnSections", "");
      if (!sections.is_array())
      {
        throw std::runtime_error("Error: ColumnSections must be an array");
      }

      std::vector<double> sectionLengths(adsorbentsJson.size(), -1.0);
      std::vector<size_t> sectionGridPoints(adsorbentsJson.size(), 0);
      std::vector<double> sectionInterfaceLengths;
      std::size_t nextAdsorbentInOrder = 0;
      for (std::size_t sectionId = 0; sectionId < sections.size(); ++sectionId)
      {
        const nlohmann::json& section = sections[sectionId];
        const std::string context = "ColumnSection " + std::to_string(sectionId);
        if (!section.is_object())
        {
          throw std::runtime_error("Error: each ColumnSections entry must be an object (" + context + ")");
        }
        requireOnlyKnownColumnSectionSettings(section, context);

        const bool hasAdsorbent = containsKeyCaseInsensitive(section, "Adsorbent");
        const bool hasInterface = containsKeyCaseInsensitive(section, "InterfaceLength");
        if (hasAdsorbent == hasInterface)
        {
          throw std::runtime_error("Error: ColumnSections entries must define either Adsorbent/Length or "
                                   "InterfaceLength (" +
                                   context + ")");
        }

        if (hasInterface)
        {
          sectionInterfaceLengths.push_back(getNumberOrThrow<double>(
              requireKeyCaseInsensitive(section, "InterfaceLength", context), "InterfaceLength", context));
          continue;
        }

        const std::string adsorbentName =
            getStringOrThrow(requireKeyCaseInsensitive(section, "Adsorbent", context), "Adsorbent", context);
        auto adsorbentIndex = adsorbentIndexByName.find(adsorbentName);
        if (adsorbentIndex == adsorbentIndexByName.end())
        {
          throw std::runtime_error("Error: unknown adsorbent '" + adsorbentName + "' (" + context + ")");
        }
        if (adsorbentIndex->second != nextAdsorbentInOrder)
        {
          throw std::runtime_error("Error: ColumnSections adsorbent order must match Adsorbents order (" + context +
                                   ")");
        }
        sectionLengths[adsorbentIndex->second] =
            getNumberOrThrow<double>(requireKeyCaseInsensitive(section, "Length", context), "Length", context);
        readOptionalNumber<size_t>(section, "NumberOfGridPoints", sectionGridPoints[adsorbentIndex->second]);
        ++nextAdsorbentInOrder;
      }

      if (nextAdsorbentInOrder != adsorbentsJson.size())
      {
        throw std::runtime_error("Error: ColumnSections must define one Length for each adsorbent");
      }
      if (sectionInterfaceLengths.size() != adsorbentsJson.size() - 1)
      {
        throw std::runtime_error("Error: ColumnSections must define one InterfaceLength between neighboring adsorbents");
      }

      adsorbentLengths = std::move(sectionLengths);
      adsorbentGridPoints = std::move(sectionGridPoints);
      adsorbentInterfaceLengths = std::move(sectionInterfaceLengths);
    }
    else if (containsKeyCaseInsensitive(parsed_data, "Interfaces"))
    {
      adsorbentInterfaceLengths =
          getNumberListOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "Interfaces", ""), "Interfaces", "");
    }
    else
    {
      adsorbentInterfaceLengths.assign(adsorbentsJson.size() - 1, 0.0);
    }

    for (std::size_t ads = 0; ads < adsorbentLengths.size(); ++ads)
    {
      if (adsorbentLengths[ads] <= 0.0)
      {
        throw std::runtime_error("Error: AdsorbentLength/ColumnSections Length must be positive for adsorbent '" +
                                 adsorbentNames[ads] + "'");
      }
      if (adsorbentVoidFractions[ads] <= 0.0 || adsorbentVoidFractions[ads] >= 1.0)
      {
        throw std::runtime_error("Error: ColumnVoidFraction must be between 0 and 1 for adsorbent '" +
                                 adsorbentNames[ads] + "'");
      }
      if (adsorbentParticleDensities[ads] < 0.0)
      {
        throw std::runtime_error("Error: ParticleDensity must be non-negative for adsorbent '" + adsorbentNames[ads] +
                                 "'");
      }
      if (adsorbentParticleDiameters[ads] <= 0.0)
      {
        throw std::runtime_error("Error: ParticleDiameter must be positive for adsorbent '" + adsorbentNames[ads] +
                                 "'");
      }
    }

    if (adsorbentInterfaceLengths.size() != adsorbentsJson.size() - 1)
    {
      throw std::runtime_error("Error: Interfaces size must equal Adsorbents size minus one");
    }

    columnLength = std::reduce(adsorbentLengths.begin(), adsorbentLengths.end(), 0.0);

    if (!columnDistances.empty())
    {
      numberOfGridPoints = columnDistances.size() - 1;
    }
    else
    {
      const bool hasSectionGridPoints =
          std::any_of(adsorbentGridPoints.begin(), adsorbentGridPoints.end(), [](size_t n) { return n != 0; });
      if (hasSectionGridPoints)
      {
        if (adsorbentGridPoints.size() != adsorbentLengths.size() ||
            std::any_of(adsorbentGridPoints.begin(), adsorbentGridPoints.end(), [](size_t n) { return n == 0; }))
        {
          throw std::runtime_error("Error: ColumnSections NumberOfGridPoints must be set for every adsorbent section");
        }

        columnDistances.clear();
        columnDistances.push_back(0.0);
        double z = 0.0;
        for (std::size_t ads = 0; ads < adsorbentLengths.size(); ++ads)
        {
          const double dz = adsorbentLengths[ads] / static_cast<double>(adsorbentGridPoints[ads]);
          for (size_t grid = 1; grid <= adsorbentGridPoints[ads]; ++grid)
          {
            columnDistances.push_back(z + static_cast<double>(grid) * dz);
          }
          z = columnDistances.back();
        }
        columnDistances.back() = columnLength;
        numberOfGridPoints = columnDistances.size() - 1;
      }
    }
  }
  else
  {
    if (!columnDistances.empty())
    {
      throw std::runtime_error("Error: ColumnDistances is only supported for multibed columns");
    }
    adsorbentComponents.push_back(components);
    adsorbentLengths.push_back(columnLength);
    adsorbentInterfaceLengths.clear();
    adsorbentGridPoints.push_back(numberOfGridPoints);
    adsorbentVoidFractions.push_back(columnVoidFraction);
    adsorbentParticleDensities.push_back(particleDensity);
    adsorbentParticleDiameters.push_back(particleDiameter);
  }

  if (columnDistances.empty())
  {
    columnDistances = makeUniformColumnDistances(numberOfGridPoints, columnLength);
  }
  validateColumnDistances(columnDistances, numberOfGridPoints, columnLength);

  for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
  {
    for (const Component& comp : componentsForAdsorbent)
    {
      if (!comp.isCarrierGas && comp.isotherm.numberOfSites == 0)
      {
        throw std::runtime_error("Error: non-carrier component '" + comp.name +
                                 "' has no isotherm for at least one adsorbent");
      }
    }
  }

  if ((mixturePredictionMethod == 2) || (mixturePredictionMethod == 3))
  {
    for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
    {
      for (size_t i = 0; i < componentsForAdsorbent.size(); ++i)
      {
        for (size_t j = 0; j < componentsForAdsorbent[i].isotherm.numberOfSites; ++j)
        {
          if (componentsForAdsorbent[i].isotherm.sites[j].type != Isotherm::Type::Langmuir)
          {
            throw std::runtime_error("Error: Explicit mixture prediction must use single Langmuir isotherms");
          }
        }
      }
    }
  }

  bool hasChemisorption = false;
  for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
  {
    for (const Component& component : componentsForAdsorbent)
    {
      hasChemisorption = hasChemisorption || component.chemisorption.numberOfSites > 0;
      if (mixturePredictionMethod == 3)
      {
        for (const Chemisorption& site : component.chemisorption.sites)
        {
          if (site.isotherm.has_value() && site.isotherm->type != Isotherm::Type::Langmuir)
          {
            throw std::runtime_error(
                "Error: SEI chemisorption mixture prediction requires Langmuir isotherms");
          }
        }
      }
    }
  }
  if (hasChemisorption && mixturePredictionMethod == 2)
  {
    throw std::runtime_error(
        "Error: chemisorption mixture prediction supports IAST, SIAST, or SEI; EI is not supported");
  }

  maxIsothermTerms = 0;
  for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
  {
    if (!componentsForAdsorbent.empty())
    {
      std::vector<Component>::const_iterator maxIsothermTermsIterator =
          std::max_element(componentsForAdsorbent.begin(), componentsForAdsorbent.end(),
                           [](const Component& lhs, const Component& rhs)
                           { return lhs.isotherm.numberOfSites < rhs.isotherm.numberOfSites; });
      maxIsothermTerms = std::max(maxIsothermTerms, maxIsothermTermsIterator->isotherm.numberOfSites);
    }
  }

  if (simulationType == SimulationType::Breakthrough)
  {
    static constexpr size_t inletPressureInletVelocity = 0;
    static constexpr size_t inletPressureOutletPressure = 1;
    static constexpr size_t inletVelocityOutletPressure = 2;
    static constexpr size_t fixedVelocity = 3;
    static constexpr size_t fixedPressureInletVelocity = 4;

    // Boundary conditions define which externally supplied quantities are mandatory.
    const bool boundaryNeedsInletPressure =
        boundaryCondition == inletPressureInletVelocity || boundaryCondition == inletPressureOutletPressure ||
        boundaryCondition == fixedPressureInletVelocity;
    const bool boundaryNeedsOutletPressure =
        boundaryCondition == inletPressureOutletPressure || boundaryCondition == inletVelocityOutletPressure;
    const bool boundaryNeedsInletVelocity =
        boundaryCondition == inletPressureInletVelocity || boundaryCondition == inletVelocityOutletPressure ||
        boundaryCondition == fixedVelocity || boundaryCondition == fixedPressureInletVelocity;

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

    if (boundaryNeedsInletVelocity)
    {
      requireConfigured(hasColumnEntranceVelocity && columnEntranceVelocity >= 0.0,
                        "Error: inlet velocity not set (Use e.g.: 'ColumnEntranceVelocity 0.1')");
    }

    if (boundaryCondition == fixedVelocity)
    {
      requireConfigured((hasInletPressure && inletPressure >= 0.0) || (hasOutletPressure && outletPressure >= 0.0),
                        "Error: FixedVelocity requires InletPressure or OutletPressure");
    }

    if (boundaryCondition == fixedPressureInletVelocity)
    {
      requireConfigured(hasPressureGradient, "Error: FixedPressureInletVelocity requires PressureGradient");
      requireConfigured(inletPressure + pressureGradient / columnLength > 0.0,
                        "Error: fixed pressure profile becomes non-positive");
    }
  }
}

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
#include <optional>
#include <print>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "json.h"
#include "reaction.h"
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
    const nlohmann::json& value = requireKeyCaseInsensitive(object, key, "");
    if (value.is_boolean())
    {
      target = value.get<bool>();
      return;
    }

    throw std::runtime_error("Error: key '" + key + "' has invalid boolean value");
  }
}

static void requireOnlyKnownKeys(const nlohmann::json& object, std::initializer_list<std::string_view> allowedKeys,
                                 const std::string& context, bool allowIsothermTypes);

static void requireOnlyExactKeys(const nlohmann::json& object, std::initializer_list<std::string_view> allowedKeys,
                                 const std::string& context)
{
  for (auto it = object.begin(); it != object.end(); ++it)
  {
    const bool known =
        std::any_of(allowedKeys.begin(), allowedKeys.end(), [&](std::string_view key) { return it.key() == key; });
    if (!known)
    {
      throw std::runtime_error("Error: unknown key '" + it.key() + "' (" + context + ")");
    }
  }
}

static void requireExactKeys(const nlohmann::json& object, std::initializer_list<std::string_view> requiredKeys,
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

struct IsothermSpec
{
  Isotherm::Type type;
  std::size_t parameterCount;
  bool nonIsothermalImplemented;
};

static const IsothermSpec* findIsothermSpec(const std::string& typeString);

static void readChemisorption(MultiSiteChemisorption& chemisorption, const nlohmann::json& object, bool nonIsothermal,
                              const std::string& context)
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
      const nlohmann::json& value = (*sites)[site];
      const std::string siteContext = context + ", ChemisorptionSite " + std::to_string(site);
      if (!value.is_object())
      {
        throw std::runtime_error("Error: Chemisorption site must be an object (" + siteContext + ")");
      }

      requireOnlyExactKeys(value, {"Type", "Parameters"}, siteContext);

      const std::string chemisorptionTypeString =
          getStringOrThrow(requireKeyCaseInsensitive(value, "Type", siteContext), "Type", siteContext);
      Chemisorption siteChemisorption;
      if (caseInSensStringCompare(chemisorptionTypeString, "None"))
        siteChemisorption.type = Chemisorption::Type::None;
      else if (caseInSensStringCompare(chemisorptionTypeString, "FirstOrder"))
        siteChemisorption.type = Chemisorption::Type::FirstOrder;
      else if (caseInSensStringCompare(chemisorptionTypeString, "PseudoNth"))
        siteChemisorption.type = Chemisorption::Type::PseudoNth;
      else if (caseInSensStringCompare(chemisorptionTypeString, "Avrami"))
        siteChemisorption.type = Chemisorption::Type::Avrami;
      else if (caseInSensStringCompare(chemisorptionTypeString, "General"))
        siteChemisorption.type = Chemisorption::Type::General;
      else if (caseInSensStringCompare(chemisorptionTypeString, "Elovich"))
        siteChemisorption.type = Chemisorption::Type::Elovich;
      else
      {
        throw std::runtime_error("Error: invalid Chemisorption Type '" + chemisorptionTypeString + "'");
      }

      const nlohmann::json& params = requireKeyCaseInsensitive(value, "Parameters", siteContext);
      if (!params.is_object())
      {
        throw std::runtime_error("Error: Chemisorption Parameters must be an object (" + siteContext + ")");
      }

      const std::string parametersContext = siteContext + " Chemisorption Parameters";
      switch (siteChemisorption.type)
      {
        case Chemisorption::Type::None:
          throw std::runtime_error("Error: None is not a valid ChemisorptionSites model (" + siteContext + ")");
        case Chemisorption::Type::FirstOrder:
          requireExactKeys(params, {"rateCoefficient", "maximumLoading", "heatOfChemisorption", "Isotherm"},
                           parametersContext);
          break;
        case Chemisorption::Type::PseudoNth:
        case Chemisorption::Type::Avrami:
          requireExactKeys(params, {"rateCoefficient", "order", "maximumLoading", "heatOfChemisorption", "Isotherm"},
                           parametersContext);
          break;
        case Chemisorption::Type::General:
          requireExactKeys(params,
                           {"maximumLoading", "heatOfChemisorption", "adsorptionRateCoefficient",
                            "adsorptionActivationEnergy", "desorptionRateCoefficient", "desorptionActivationEnergy",
                            "poreConcentrationOrder", "capacityOrder", "desorptionOrder", "filmMassTransferCoefficient",
                            "poreDiffusivity", "usePoreSurfaceTransport", "Isotherm"},
                           parametersContext);
          break;
        case Chemisorption::Type::Elovich:
          requireExactKeys(params,
                           {"maximumLoading", "heatOfChemisorption", "alpha", "beta", "filmMassTransferCoefficient",
                            "poreDiffusivity", "usePoreSurfaceTransport", "Isotherm"},
                           parametersContext);
          break;
      }

      readOptionalNumber<double>(params, "rateCoefficient", siteChemisorption.rateCoefficient);
      readOptionalNonNegativeInteger(params, "order", siteChemisorption.order);
      readOptionalNumber<double>(params, "maximumLoading", siteChemisorption.maximumLoading);
      readOptionalNumber<double>(params, "heatOfChemisorption", siteChemisorption.heatOfChemisorption);
      readOptionalNumber<double>(params, "adsorptionRateCoefficient", siteChemisorption.adsorptionRateCoefficient);
      readOptionalNumber<double>(params, "adsorptionActivationEnergy", siteChemisorption.adsorptionActivationEnergy);
      readOptionalNumber<double>(params, "desorptionRateCoefficient", siteChemisorption.desorptionRateCoefficient);
      readOptionalNumber<double>(params, "desorptionActivationEnergy", siteChemisorption.desorptionActivationEnergy);
      readOptionalNonNegativeInteger(params, "poreConcentrationOrder", siteChemisorption.poreConcentrationOrder);
      readOptionalNonNegativeInteger(params, "capacityOrder", siteChemisorption.capacityOrder);
      readOptionalNonNegativeInteger(params, "desorptionOrder", siteChemisorption.desorptionOrder);
      readOptionalNumber<double>(params, "alpha", siteChemisorption.elovichAlpha);
      readOptionalNumber<double>(params, "beta", siteChemisorption.elovichBeta);
      readOptionalNumber<double>(params, "filmMassTransferCoefficient", siteChemisorption.filmMassTransferCoefficient);
      readOptionalNumber<double>(params, "poreDiffusivity", siteChemisorption.poreDiffusivity);
      readOptionalBool(params, "usePoreSurfaceTransport", siteChemisorption.usePoreSurfaceTransport);

      const std::string isothermContext = siteContext + " Isotherm";
      const nlohmann::json& isothermValue = requireKeyCaseInsensitive(params, "Isotherm", siteContext);
      if (!isothermValue.is_object())
      {
        throw std::runtime_error("Error: Isotherm must be an object (" + isothermContext + ")");
      }
      requireOnlyExactKeys(isothermValue, {"Type", "Parameters"}, isothermContext);

      const std::string isothermTypeString =
          getStringOrThrow(requireKeyCaseInsensitive(isothermValue, "Type", isothermContext), "Type", isothermContext);
      const IsothermSpec* spec = findIsothermSpec(isothermTypeString);
      if (spec == nullptr)
      {
        throw std::runtime_error("Error: unknown isotherm type '" + isothermTypeString + "' (" + isothermContext + ")");
      }
      if (nonIsothermal && !spec->nonIsothermalImplemented)
      {
        throw std::logic_error("Error: nonIsothermal not implemented for " + isothermTypeString);
      }

      const nlohmann::json& parameters = requireKeyCaseInsensitive(isothermValue, "Parameters", isothermContext);
      std::vector<double> values =
          requireDoubleParameterCount(getNumberListOrThrow<double>(parameters, isothermTypeString, isothermContext),
                                      spec->parameterCount, isothermTypeString, isothermContext);
      siteChemisorption.isotherm = Isotherm(spec->type, values, nonIsothermal);

      if (siteChemisorption.type == Chemisorption::Type::General)
      {
        if (siteChemisorption.adsorptionRateCoefficient < 0.0 || siteChemisorption.desorptionRateCoefficient < 0.0 ||
            siteChemisorption.filmMassTransferCoefficient < 0.0 || siteChemisorption.poreDiffusivity < 0.0)
        {
          throw std::runtime_error(
              "Error: General Chemisorption kinetic and transport coefficients must be non-negative (" + siteContext +
              ")");
        }
      }
      else if (siteChemisorption.type == Chemisorption::Type::Elovich)
      {
        if (siteChemisorption.elovichAlpha < 0.0 || siteChemisorption.elovichBeta < 0.0 ||
            siteChemisorption.filmMassTransferCoefficient < 0.0 || siteChemisorption.poreDiffusivity < 0.0)
        {
          throw std::runtime_error(
              "Error: Elovich Chemisorption kinetic and transport coefficients must be non-negative (" + siteContext +
              ")");
        }
      }
      else if (siteChemisorption.rateCoefficient < 0.0)
      {
        throw std::runtime_error("Error: Chemisorption rateCoefficient must be non-negative (" + siteContext + ")");
      }
      if (siteChemisorption.maximumLoading <= 0.0)
      {
        throw std::runtime_error("Error: Chemisorption maximumLoading must be positive (" + siteContext + ")");
      }
      chemisorption.add(siteChemisorption);
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

static const IsothermSpec* findIsothermSpec(const std::string& typeString)
{
  static const std::vector<std::pair<std::string_view, IsothermSpec>> specs{
      {"Langmuir", {Isotherm::Type::Langmuir, 2, true}},
      {"pH-Langmuir", {Isotherm::Type::Langmuir_pH, 3, true}},
      {"Anti-Langmuir", {Isotherm::Type::Anti_Langmuir, 2, false}},
      {"BET", {Isotherm::Type::BET, 3, false}},
      {"Henry", {Isotherm::Type::Henry, 1, false}},
      {"Freundlich", {Isotherm::Type::Freundlich, 2, false}},
      {"Sips", {Isotherm::Type::Sips, 3, true}},
      {"Langmuir-Freundlich", {Isotherm::Type::Langmuir_Freundlich, 3, true}},
      {"Redlich-Peterson", {Isotherm::Type::Redlich_Peterson, 3, false}},
      {"Toth", {Isotherm::Type::Toth, 3, true}},
      {"Unilan", {Isotherm::Type::Unilan, 3, false}},
      {"OBrien&Myers", {Isotherm::Type::OBrien_Myers, 3, false}},
      {"Quadratic", {Isotherm::Type::Quadratic, 3, false}},
      {"Temkin", {Isotherm::Type::Temkin, 3, false}},
      {"Bingel&Walton", {Isotherm::Type::BingelWalton, 3, false}},
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

static bool supportsSegregatedCompetitiveIsotherm(Isotherm::Type type)
{
  switch (type)
  {
    case Isotherm::Type::Langmuir:
    case Isotherm::Type::Anti_Langmuir:
    case Isotherm::Type::Sips:
    case Isotherm::Type::Langmuir_Freundlich:
    case Isotherm::Type::Redlich_Peterson:
    case Isotherm::Type::Toth:
      return true;
    default:
      return false;
  }
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
    const bool known = std::any_of(allowedKeys.begin(), allowedKeys.end(), [&](std::string_view allowedKey)
                                   { return caseInSensStringCompare(it.key(), std::string{allowedKey}); });
    if (known)
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

static void readPhysisorptionSites(Component& comp, const nlohmann::json& item, const std::string& context)
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

    if (site.is_object() && containsKeyCaseInsensitive(site, "Type") && containsKeyCaseInsensitive(site, "Parameters"))
    {
      requireOnlyKnownKeys(site, {"Type", "Parameters"}, siteContext);

      std::string typeString =
          getStringOrThrow(requireKeyCaseInsensitive(site, "Type", siteContext), "Type", siteContext);
      const nlohmann::json& params = requireKeyCaseInsensitive(site, "Parameters", siteContext);
      const IsothermSpec* spec = findIsothermSpec(typeString);
      if (spec == nullptr)
      {
        throw std::runtime_error("Error: unknown isotherm type '" + typeString + "'" +
                                 (siteContext.empty() ? "" : (" (" + siteContext + ")")));
      }

      std::vector<double> values = requireDoubleParameterCount(
          getNumberListOrThrow<double>(params, typeString, siteContext), spec->parameterCount, typeString, siteContext);

      if (comp.nonIsothermal)
      {
        if (!spec->nonIsothermalImplemented)
        {
          throw std::logic_error("Error: nonIsothermal not implemented for " + typeString);
        }
      }

      comp.isotherm.add(Isotherm(spec->type, values, comp.nonIsothermal));
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
  requireOnlyKnownKeys(item,
                       {"Name", "FileName", "CarrierGas", "GasPhaseMolFraction", "LiquidPhaseConcentration",
                        "InitialLiquidPhaseConcentration", "MassTransferCoefficient",
                        "AxialDispersionCoefficient", "MolecularWeight", "HeatOfAdsorption", "referenceTemperature",
                        "nonIsothermal", "ChemisorptionSites", "PhysisorptionSites"},
                       context);

  std::string componentName = getStringOrThrow(requireKeyCaseInsensitive(item, "Name", context), "Name", context);

  Component comp(componentId, componentName);

  readOptionalString(item, "FileName", comp.filename);
  readOptionalBool(item, "CarrierGas", comp.isCarrierGas);
  readOptionalNumber<double>(item, "GasPhaseMolFraction", comp.initialGasMoleFraction);
  readOptionalNumber<double>(item, "LiquidPhaseConcentration", comp.inletLiquidConcentration);
  readOptionalNumber<double>(item, "InitialLiquidPhaseConcentration", comp.initialLiquidConcentration);
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
  return comp;
}

static void requireConfigured(bool condition, const std::string& message)
{
  if (!condition)
  {
    throw std::runtime_error(message);
  }
}

static const nlohmann::json* findFirstKeyCaseInsensitive(const nlohmann::json& object,
                                                         std::initializer_list<std::string_view> keys)
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

static std::vector<double> readSizedDoubleList(const nlohmann::json& value, const std::string& key,
                                               const std::string& context, size_t expected)
{
  std::vector<double> values = getNumberListOrThrow<double>(value, key, context);
  if (values.size() != expected)
  {
    throw std::runtime_error("Error: " + key + " must contain " + std::to_string(expected) + " values (" + context +
                             ")");
  }
  for (double valueItem : values)
  {
    if (!std::isfinite(valueItem) || valueItem <= 0.0)
    {
      throw std::runtime_error("Error: " + key + " values must be positive finite numbers (" + context + ")");
    }
  }
  return values;
}

static std::vector<size_t> readReactionParticipants(const nlohmann::json& value,
                                                    const std::vector<Component>& components, const std::string& key,
                                                    const std::string& context)
{
  if (!value.is_array())
  {
    throw std::runtime_error("Error: Reaction " + key +
                             " must be an array of zero-based component indices or component names (" + context + ")");
  }

  std::vector<size_t> participants;
  participants.reserve(value.size());
  for (const nlohmann::json& item : value)
  {
    if (item.is_string())
    {
      const std::string name = item.get<std::string>();
      bool foundExact = false;
      for (size_t comp = 0; comp < components.size(); ++comp)
      {
        if (components[comp].name == name)
        {
          participants.push_back(comp);
          foundExact = true;
          break;
        }
      }
      if (foundExact)
      {
        continue;
      }

      size_t match = components.size();
      for (size_t comp = 0; comp < components.size(); ++comp)
      {
        if (!caseInSensStringCompare(components[comp].name, name)) continue;
        if (match != components.size())
        {
          throw std::runtime_error("Error: ambiguous component name '" + name + "' in Reaction " + key + " (" +
                                   context + ")");
        }
        match = comp;
      }

      if (match != components.size())
      {
        participants.push_back(match);
        continue;
      }

      throw std::runtime_error("Error: component name '" + name + "' in Reaction " + key + " was not found (" +
                               context + ")");
    }

    if (item.is_number_unsigned())
    {
      const size_t index = item.get<size_t>();
      if (index < components.size())
      {
        participants.push_back(index);
        continue;
      }
    }
    else if (item.is_number_integer())
    {
      const std::int64_t index = item.get<std::int64_t>();
      if (index >= 0 && static_cast<size_t>(index) < components.size())
      {
        participants.push_back(static_cast<size_t>(index));
        continue;
      }
    }

    throw std::runtime_error("Error: Reaction " + key +
                             " must contain zero-based component indices or component names (" + context + ")");
  }
  return participants;
}

static double geometryDoubleOrDefault(const nlohmann::json& object, std::initializer_list<std::string_view> keys,
                                      double defaultValue)
{
  const nlohmann::json* value = findFirstKeyCaseInsensitive(object, keys);
  return value == nullptr ? defaultValue : getNumberOrThrow<double>(*value, std::string{*keys.begin()}, "Geometry");
}

static double geometryRequiredDouble(const nlohmann::json& object, std::initializer_list<std::string_view> keys)
{
  const nlohmann::json* value = findFirstKeyCaseInsensitive(object, keys);
  if (value == nullptr)
  {
    throw std::runtime_error("Error: required Geometry key '" + std::string{*keys.begin()} + "' missing");
  }
  return getNumberOrThrow<double>(*value, std::string{*keys.begin()}, "Geometry");
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
  requireOnlyKnownKeys(parsed_data,
                       {"SimulationType",
                        "MixturePredictionMethod",
                        "MPDSettings",
                        "IASTMethod",
                        "BreakthroughIntegrator",
                        "BoundaryCondition",
                        "FluidPhase",
                        "LiquidDensity",
                        "pHMode",
                        "pHValue",
                        "pKw",
                        "pHComponent",
                        "PressureScale",
                        "ReadColumnFile",
                        "DisplayName",
                        "DebugForceMultibed",
                        "Temperature",
                        "ColumnVoidFraction",
                        "DynamicViscosity",
                        "ParticleDiameter",
                        "ParticleDensity",
                        "InletPressure",
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
                        "wallDensity",
                        "gasThermalConductivity",
                        "wallThermalConductivity",
                        "heatTransferGasWall",
                        "heatTransferGasSolid",
                        "heatTransferWallExternal",
                        "heatCapacityGas",
                        "liquidThermalConductivity",
                        "heatTransferLiquidSolid",
                        "heatTransferLiquidWall",
                        "heatCapacityLiquid",
                        "heatCapacitySolid",
                        "heatCapacityWall",
                        "energyBalance",
                        "SwingAdsorptionPhases",
                        "ColumnSections",
                        "Adsorbents",
                        "Components",
                        "Reactions"},
                       "general settings");

  // Track presence separately from value, because some numeric defaults are valid runtime values.
  const bool hasInletPressure = containsKeyCaseInsensitive(parsed_data, "InletPressure");
  const bool hasOutletPressure = containsKeyCaseInsensitive(parsed_data, "OutletPressure");
  const bool hasColumnEntranceVelocity = containsKeyCaseInsensitive(parsed_data, "ColumnEntranceVelocity");
  const bool hasPressureGradient = containsKeyCaseInsensitive(parsed_data, "PressureGradient");
  const bool hasTemperature = containsKeyCaseInsensitive(parsed_data, "Temperature");
  const bool hasInfluxTemperature = containsKeyCaseInsensitive(parsed_data, "InfluxTemperature");

  readOptionalMappedString(parsed_data, "SimulationType", simulationType,
                           {{"Breakthrough", SimulationType::Breakthrough},
                            {"MixturePrediction", SimulationType::MixturePrediction},
                            {"Fitting", SimulationType::Fitting},
                            {"SwingAdsorption", SimulationType::SwingAdsorption},
                            {"Test", SimulationType::Test}});

  readOptionalMappedString(parsed_data, "MixturePredictionMethod", mixturePredictionMethod,
                           {{"IAST", 0}, {"SIAST", 1}, {"EI", 2}, {"SEI", 3}, {"SCI", 4}, {"SPI", 5}, {"MPD", 6}});

  const nlohmann::json* mpdJson = findKeyCaseInsensitive(parsed_data, "MPDSettings");
  if (mixturePredictionMethod == 6 && mpdJson == nullptr)
  {
    throw std::runtime_error("Error: MixturePredictionMethod MPD requires MPDSettings");
  }
  if (mpdJson != nullptr)
  {
    if (!mpdJson->is_object()) throw std::runtime_error("Error: MPDSettings must be an object");
    requireOnlyKnownKeys(*mpdJson,
                         {"FileName", "ReferenceTemperature", "ReferenceFugacity", "ReferenceFrameworkMass",
                          "FrameworkMass", "ComponentBounds"},
                         "MPDSettings");

    MPDSettings settings;
    settings.fileName = getStringOrThrow(requireKeyCaseInsensitive(*mpdJson, "FileName", "MPDSettings"),
                                         "FileName", "MPDSettings");
    settings.referenceTemperature = getNumberOrThrow<double>(
        requireKeyCaseInsensitive(*mpdJson, "ReferenceTemperature", "MPDSettings"), "ReferenceTemperature",
        "MPDSettings");
    settings.referenceFugacity =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(*mpdJson, "ReferenceFugacity", "MPDSettings"),
                                 "ReferenceFugacity", "MPDSettings");

    const nlohmann::json* referenceMass = findKeyCaseInsensitive(*mpdJson, "ReferenceFrameworkMass");
    const nlohmann::json* frameworkMass = findKeyCaseInsensitive(*mpdJson, "FrameworkMass");
    if (referenceMass != nullptr && frameworkMass != nullptr)
    {
      throw std::runtime_error(
          "Error: MPDSettings must use only one of ReferenceFrameworkMass or FrameworkMass");
    }
    if (referenceMass == nullptr && frameworkMass == nullptr)
    {
      throw std::runtime_error("Error: required key 'ReferenceFrameworkMass' missing (MPDSettings)");
    }
    settings.referenceFrameworkMass = getNumberOrThrow<double>(
        referenceMass != nullptr ? *referenceMass : *frameworkMass, "ReferenceFrameworkMass", "MPDSettings");

    const nlohmann::json& boundsJson =
        requireKeyCaseInsensitive(*mpdJson, "ComponentBounds", "MPDSettings");
    if (!boundsJson.is_array() || boundsJson.empty())
    {
      throw std::runtime_error("Error: MPD ComponentBounds must be a non-empty array");
    }

    const auto readParticleNumber = [](const nlohmann::json& value, const std::string& key,
                                       const std::string& context) -> size_t
    {
      try
      {
        if (value.is_number_unsigned()) return value.get<size_t>();
        if (value.is_number_integer())
        {
          const std::int64_t integer = value.get<std::int64_t>();
          if (integer >= 0) return static_cast<size_t>(integer);
        }
      }
      catch (const nlohmann::json::exception&)
      {
      }
      throw std::runtime_error("Error: key '" + key + "' must be a non-negative integer (" + context + ")");
    };

    settings.componentBounds.reserve(boundsJson.size());
    for (size_t component = 0; component < boundsJson.size(); ++component)
    {
      const nlohmann::json& item = boundsJson[component];
      const std::string context = "MPDSettings ComponentBounds " + std::to_string(component);
      if (!item.is_object()) throw std::runtime_error("Error: MPD ComponentBounds entries must be objects");
      requireOnlyKnownKeys(item, {"Component", "NMin", "NMax", "DeltaN"}, context);

      MPDComponentBounds bounds;
      readOptionalString(item, "Component", bounds.component);
      bounds.nMin = readParticleNumber(requireKeyCaseInsensitive(item, "NMin", context), "NMin", context);
      bounds.nMax = readParticleNumber(requireKeyCaseInsensitive(item, "NMax", context), "NMax", context);
      bounds.deltaN = readParticleNumber(requireKeyCaseInsensitive(item, "DeltaN", context), "DeltaN", context);
      settings.componentBounds.push_back(std::move(bounds));
    }

    if (!(settings.referenceTemperature > 0.0))
    {
      throw std::runtime_error("Error: MPD ReferenceTemperature must be positive");
    }
    if (!(settings.referenceFugacity > 0.0))
    {
      throw std::runtime_error("Error: MPD ReferenceFugacity must be positive");
    }
    if (!(settings.referenceFrameworkMass > 0.0))
    {
      throw std::runtime_error("Error: MPD ReferenceFrameworkMass must be positive");
    }
    for (const MPDComponentBounds& bounds : settings.componentBounds)
    {
      if (bounds.deltaN == 0) throw std::runtime_error("Error: MPD DeltaN must be positive");
      if (bounds.nMax < bounds.nMin) throw std::runtime_error("Error: MPD NMax must be at least NMin");
      if ((bounds.nMax - bounds.nMin) % bounds.deltaN != 0)
      {
        throw std::runtime_error("Error: MPD [NMin, NMax] must be exactly divisible by DeltaN");
      }
    }

    std::filesystem::path distributionPath{settings.fileName};
    if (distributionPath.is_relative())
    {
      distributionPath = std::filesystem::path{fileName}.parent_path() / distributionPath;
    }
    settings.fileName = distributionPath.lexically_normal().string();
    mpdSettings = std::move(settings);

    if (!hasTemperature)
    {
      std::print(stderr,
                 "Warning: MPD ReferenceTemperature is set, but target Temperature is not; using the default target "
                 "temperature of {} K.\n",
                 temperature);
    }
  }

  readOptionalMappedString(parsed_data, "IASTMethod", IASTMethod, {{"FastIAST", 0}, {"NestedLoopBisection", 1}});

  readOptionalMappedString(parsed_data, "BreakthroughIntegrator", breakthroughIntegrator,
                           {{"RungeKutta3", 0}, {"CVODE", 1}, {"SIRK3", 2}});

  readOptionalMappedString(parsed_data, "BoundaryCondition", boundaryCondition,
                           {{"InletPressureInletVelocity", 0},
                            {"InletPressureOutletPressure", 1},
                            {"InletVelocityOutletPressure", 2},
                            {"FixedVelocity", 3},
                            {"FixedPressureInletVelocity", 4}});

  readOptionalMappedString(parsed_data, "FluidPhase", fluidPhase, {{"Gas", 0}, {"Liquid", 1}});
  readOptionalMappedString(parsed_data, "pHMode", pHMode, {{"Fixed", 0}, {"HPlus", 1}, {"OHMinus", 2}});

  readOptionalMappedString(parsed_data, "PressureScale", pressureScale, {{"Log", 0}, {"Linear", 1}});

  readOptionalString(parsed_data, "ReadColumnFile", readColumnFile);
  readOptionalString(parsed_data, "DisplayName", displayName);
  readOptionalBool(parsed_data, "DebugForceMultibed", debugForceMultibed);

  readOptionalNumber<double>(parsed_data, "Temperature", temperature);
  readOptionalNumber<double>(parsed_data, "ColumnVoidFraction", columnVoidFraction);
  readOptionalNumber<double>(parsed_data, "DynamicViscosity", dynamicViscosity);
  readOptionalNumber<double>(parsed_data, "LiquidDensity", liquidDensity);
  readOptionalNumber<double>(parsed_data, "pHValue", pHValue);
  readOptionalNumber<double>(parsed_data, "pKw", pKw);
  readOptionalNumber<double>(parsed_data, "ParticleDiameter", particleDiameter);
  readOptionalNumber<double>(parsed_data, "ParticleDensity", particleDensity);

  readOptionalNumber<double>(parsed_data, "InletPressure", inletPressure);

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
    columnDistances = getNumberListOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ColumnDistances", ""),
                                                   "ColumnDistances", "");
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

  readOptionalNumber<double>(parsed_data, "wallDensity", wallDensity);
  readOptionalNumber<double>(parsed_data, "gasThermalConductivity", gasThermalConductivity);
  readOptionalNumber<double>(parsed_data, "wallThermalConductivity", wallThermalConductivity);
  readOptionalNumber<double>(parsed_data, "heatTransferGasWall", heatTransferGasWall);
  readOptionalNumber<double>(parsed_data, "heatTransferGasSolid", heatTransferGasSolid);
  readOptionalNumber<double>(parsed_data, "heatTransferWallExternal", heatTransferWallExternal);
  readOptionalNumber<double>(parsed_data, "heatCapacityGas", heatCapacityGas);
  readOptionalNumber<double>(parsed_data, "liquidThermalConductivity", liquidThermalConductivity);
  readOptionalNumber<double>(parsed_data, "heatTransferLiquidSolid", heatTransferLiquidSolid);
  readOptionalNumber<double>(parsed_data, "heatTransferLiquidWall", heatTransferLiquidWall);
  readOptionalNumber<double>(parsed_data, "heatCapacityLiquid", heatCapacityLiquid);
  readOptionalNumber<double>(parsed_data, "heatCapacitySolid", heatCapacitySolid);
  readOptionalNumber<double>(parsed_data, "heatCapacityWall", heatCapacityWall);
  readOptionalBool(parsed_data, "energyBalance", energyBalance);

  const bool requireGeometry =
      simulationType == SimulationType::Breakthrough || simulationType == SimulationType::SwingAdsorption;
  const nlohmann::json* geometryJson = findKeyCaseInsensitive(parsed_data, "Geometry");
  if (geometryJson == nullptr)
  {
    if (requireGeometry)
    {
      throw std::runtime_error("Error: Geometry is required for Breakthrough and SwingAdsorption simulations");
    }
    geometry = makeGeometry(PackedBedTubeSpec{.voidFraction = columnVoidFraction,
                                              .particleDiameter = particleDiameter,
                                              .internalDiameter = internalDiameter,
                                              .outerDiameter = outerDiameter});
  }
  else
  {
    if (!geometryJson->is_object())
    {
      throw std::runtime_error("Error: Geometry must be a JSON object");
    }

    requireOnlyKnownKeys(
        *geometryJson,
        {"Type", "ColumnVoidFraction", "ParticleDiameter", "InternalDiameter", "OuterDiameter", "ChannelShape",
         "InternalChannelDimension", "NumberOfChannels", "WashcoatThickness", "WashcoatVolumePerChannelVolume",
         "ForchheimerCoefficient"},
        "Geometry");

    const std::string geometryType =
        getStringOrThrow(requireKeyCaseInsensitive(*geometryJson, "Type", "Geometry"), "Type", "Geometry");

    if (caseInSensStringCompare(geometryType, "PackedBed"))
    {
      const double eps = geometryDoubleOrDefault(*geometryJson, {"ColumnVoidFraction"}, columnVoidFraction);
      const double dp = geometryDoubleOrDefault(*geometryJson, {"ParticleDiameter"}, particleDiameter);
      const double di = geometryDoubleOrDefault(*geometryJson, {"InternalDiameter"}, internalDiameter);
      const double douter = geometryDoubleOrDefault(*geometryJson, {"OuterDiameter"}, outerDiameter);
      geometry = makeGeometry(PackedBedTubeSpec{
          .voidFraction = eps, .particleDiameter = dp, .internalDiameter = di, .outerDiameter = douter});
    }
    else if (caseInSensStringCompare(geometryType, "Monolith"))
    {
      const nlohmann::json* channelShapeValue = findFirstKeyCaseInsensitive(*geometryJson, {"ChannelShape"});
      const std::string channelShape =
          channelShapeValue == nullptr ? "circular" : getStringOrThrow(*channelShapeValue, "ChannelShape", "Geometry");
      const double dInt = geometryRequiredDouble(*geometryJson, {"InternalChannelDimension"});
      const double dOut = geometryRequiredDouble(*geometryJson, {"OuterDiameter"});
      const nlohmann::json* numberOfChannelsValue = findFirstKeyCaseInsensitive(*geometryJson, {"NumberOfChannels"});
      if (numberOfChannelsValue == nullptr)
      {
        throw std::runtime_error("Error: required Geometry key 'NumberOfChannels' missing");
      }

      std::size_t numberOfChannels = 0;
      try
      {
        if (numberOfChannelsValue->is_number_unsigned())
        {
          numberOfChannels = numberOfChannelsValue->get<std::size_t>();
        }
        else if (numberOfChannelsValue->is_number_integer())
        {
          const std::int64_t integer = numberOfChannelsValue->get<std::int64_t>();
          if (integer > 0)
          {
            numberOfChannels = static_cast<std::size_t>(integer);
          }
        }
      }
      catch (const nlohmann::json::exception&)
      {
      }
      if (numberOfChannels == 0)
      {
        throw std::runtime_error("Error: Geometry key 'NumberOfChannels' must be a positive integer");
      }

      const double washcoatThickness = geometryDoubleOrDefault(*geometryJson, {"WashcoatThickness"}, 0.0);
      const double forchheimerCoefficient =
          geometryDoubleOrDefault(*geometryJson, {"ForchheimerCoefficient"}, 0.0);
      const nlohmann::json* washcoatVolumeValue =
          findFirstKeyCaseInsensitive(*geometryJson, {"WashcoatVolumePerChannelVolume"});
      std::optional<double> washcoatVolume;
      if (washcoatVolumeValue != nullptr)
      {
        const double value =
            getNumberOrThrow<double>(*washcoatVolumeValue, "WashcoatVolumePerChannelVolume", "Geometry");
        if (value != -1.0) washcoatVolume = value;
      }
      geometry = makeGeometry(MonolithSpec{.channelShape = parseChannelShape(channelShape),
                                           .internalChannelDimension = dInt,
                                           .outerDiameter = dOut,
                                           .numberOfChannels = numberOfChannels,
                                           .washcoatThickness = washcoatThickness,
                                           .washcoatVolumePerChannelVolume = washcoatVolume,
                                           .forchheimerCoefficient = forchheimerCoefficient});
    }
    else
    {
      throw std::runtime_error("Error: Geometry Type must be PackedBed or Monolith");
    }
  }

  columnVoidFraction = geometry.voidFraction;
  if (geometry.dimensions.solidCharacteristicLength > 0.0)
  {
    particleDiameter = geometry.dimensions.solidCharacteristicLength;
  }
  internalDiameter = geometry.dimensions.internalDiameter;
  outerDiameter = geometry.dimensions.outerDiameter;

  swingAdsorptionPhases.clear();
  const nlohmann::json* phases = findKeyCaseInsensitive(parsed_data, "SwingAdsorptionPhases");
  if (phases != nullptr)
  {
    if (!phases->is_array())
    {
      throw std::runtime_error("Error: SwingAdsorptionPhases must be an array");
    }

    swingAdsorptionPhases.reserve(phases->size());
    for (size_t i = 0; i < phases->size(); ++i)
    {
      const nlohmann::json& value = (*phases)[i];
      const std::string context = "SwingAdsorptionPhases entry " + std::to_string(i);
      if (!value.is_object())
      {
        throw std::runtime_error("Error: each SwingAdsorptionPhases entry must be an object (" + context + ")");
      }

      requireOnlyKnownKeys(value, {"Name", "Temperature", "InletPressure", "NumberOfTimeSteps"}, context);

      SwingAdsorptionPhase phase;
      phase.name = "Phase " + std::to_string(i + 1);
      readOptionalString(value, "Name", phase.name);
      readOptionalNumber<double>(value, "Temperature", phase.temperature);
      readOptionalNumber<double>(value, "InletPressure", phase.inletPressure);

      const nlohmann::json* numberOfSteps = findKeyCaseInsensitive(value, "NumberOfTimeSteps");
      if (numberOfSteps == nullptr)
      {
        throw std::runtime_error("Error: required key 'NumberOfTimeSteps' missing (" + context + ")");
      }

      bool parsedNumberOfSteps = false;
      try
      {
        if (numberOfSteps->is_number_unsigned())
        {
          phase.numberOfSteps = numberOfSteps->get<size_t>();
          parsedNumberOfSteps = true;
        }
        else if (numberOfSteps->is_number_integer())
        {
          const std::int64_t integer = numberOfSteps->get<std::int64_t>();
          if (integer >= 0)
          {
            phase.numberOfSteps = static_cast<size_t>(integer);
            parsedNumberOfSteps = true;
          }
        }
      }
      catch (const nlohmann::json::exception&)
      {
      }
      if (!parsedNumberOfSteps)
      {
        throw std::runtime_error("Error: key 'NumberOfTimeSteps' must be a non-negative integer (" + context + ")");
      }
      if (phase.numberOfSteps == 0)
      {
        throw std::runtime_error("Error: SwingAdsorption phase NumberOfTimeSteps must be positive (" + context + ")");
      }
      swingAdsorptionPhases.push_back(phase);
    }
  }

  if (simulationType == SimulationType::SwingAdsorption && !swingAdsorptionPhases.empty())
  {
    const SwingAdsorptionPhase& firstPhase = swingAdsorptionPhases.front();
    if (!hasTemperature && firstPhase.temperature.has_value())
    {
      temperature = *firstPhase.temperature;
    }
    if (!hasInfluxTemperature && firstPhase.temperature.has_value())
    {
      influxTemperature = *firstPhase.temperature;
    }
    if (!hasInletPressure && firstPhase.inletPressure.has_value())
    {
      inletPressure = *firstPhase.inletPressure;
    }

    numberOfTimeSteps =
        std::accumulate(swingAdsorptionPhases.begin(), swingAdsorptionPhases.end(), size_t{0},
                        [](size_t total, const SwingAdsorptionPhase& phase) { return total + phase.numberOfSteps; });
    autoNumberOfTimeSteps = false;
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

  if (mpdSettings.has_value())
  {
    std::vector<const Component*> mpdComponents;
    for (const Component& component : components)
    {
      if (!component.isCarrierGas) mpdComponents.push_back(&component);
    }
    if (mpdSettings->componentBounds.size() != mpdComponents.size())
    {
      throw std::runtime_error("Error: MPD ComponentBounds must contain one entry for each non-carrier component");
    }
    for (size_t component = 0; component < mpdComponents.size(); ++component)
    {
      const MPDComponentBounds& bounds = mpdSettings->componentBounds[component];
      if (!bounds.component.empty() && bounds.component != mpdComponents[component]->name)
      {
        throw std::runtime_error("Error: MPD ComponentBounds entry " + std::to_string(component) + " names '" +
                                 bounds.component + "', but C-order dimension " + std::to_string(component) +
                                 " corresponds to component '" + mpdComponents[component]->name + "'");
      }
    }
  }

  // Gas equilibrium receives mole fractions. Liquid equilibrium receives concentrations directly.
  if (fluidPhase == 0 && simulationType != SimulationType::Fitting)
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
    if (fluidPhase == 0 && components[j].isCarrierGas)
    {
      carrierGasComponent = j;
      std::vector<double> values{1.0, 0.0};
      Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values, false);
      components[carrierGasComponent].isotherm.add(isotherm);

      ++numberOfCarrierGases;
    }
  }

  if (fluidPhase == 1)
  {
    for (const Component& component : components)
    {
      if (component.isCarrierGas)
      {
        throw std::runtime_error("Error: liquid simulations do not use CarrierGas components");
      }
      if (!std::isfinite(component.inletLiquidConcentration) || component.inletLiquidConcentration < 0.0 ||
          !std::isfinite(component.initialLiquidConcentration) || component.initialLiquidConcentration < 0.0)
      {
        throw std::runtime_error("Error: liquid concentrations must be finite and non-negative");
      }
    }

    if (!std::isfinite(liquidDensity) || liquidDensity <= 0.0)
    {
      throw std::runtime_error("Error: LiquidDensity must be a positive finite number");
    }

    if (pHMode != 0)
    {
      const nlohmann::json* pHComponentValue = findKeyCaseInsensitive(parsed_data, "pHComponent");
      if (pHComponentValue == nullptr)
      {
        throw std::runtime_error("Error: transported pH mode requires pHComponent");
      }
      const std::string pHComponentName = getStringOrThrow(*pHComponentValue, "pHComponent", "");
      const auto component = std::find_if(components.begin(), components.end(), [&](const Component& candidate)
                                          { return candidate.name == pHComponentName; });
      if (component == components.end())
      {
        throw std::runtime_error("Error: pHComponent '" + pHComponentName + "' was not found");
      }
      pHComponent = static_cast<size_t>(std::distance(components.begin(), component));
    }
  }

  reactions.clear();
  const nlohmann::json* reactionsJson = findKeyCaseInsensitive(parsed_data, "Reactions");
  if (reactionsJson != nullptr)
  {
    if (!reactionsJson->is_array())
    {
      throw std::runtime_error("Error: Reactions must be an array");
    }

    reactions.reserve(reactionsJson->size());
    for (size_t reactionId = 0; reactionId < reactionsJson->size(); ++reactionId)
    {
      const nlohmann::json& value = (*reactionsJson)[reactionId];
      const std::string context = "Reaction " + std::to_string(reactionId);
      if (!value.is_object())
      {
        throw std::runtime_error("Error: each Reactions entry must be an object (" + context + ")");
      }
      requireOnlyKnownKeys(value,
                           {"Phase", "Style", "Reactants", "Products", "Stoichiometry", "Kinetics", "Site",
                            "ForwardOrders", "BackwardOrders", "RateLimitTime"},
                           context);

      Reaction reaction;
      const std::string phase = getStringOrThrow(requireKeyCaseInsensitive(value, "Phase", context), "Phase", context);
      if (caseInSensStringCompare(phase, "Physisorbed"))
      {
        reaction.phase = Reaction::Phase::Physisorbed;
      }
      else if (caseInSensStringCompare(phase, "Chemisorbed"))
      {
        reaction.phase = Reaction::Phase::Chemisorbed;
      }
      else if (caseInSensStringCompare(phase, "PoreConcentration"))
      {
        reaction.phase = Reaction::Phase::PoreConcentration;
      }
      else
      {
        throw std::runtime_error("Error: Reaction Phase must be Physisorbed, Chemisorbed, or PoreConcentration (" +
                                 context + ")");
      }

      const std::string style = getStringOrThrow(requireKeyCaseInsensitive(value, "Style", context), "Style", context);
      if (caseInSensStringCompare(style, "GeneralPowerLaw"))
      {
        reaction.style = Reaction::Style::GeneralPowerLaw;
      }
      else if (caseInSensStringCompare(style, "LangmuirHinshelwood"))
      {
        reaction.style = Reaction::Style::LangmuirHinshelwood;
      }
      else if (caseInSensStringCompare(style, "LangmuirHinshelwoodHougenWatson"))
      {
        reaction.style = Reaction::Style::LangmuirHinshelwoodHougenWatson;
      }
      else
      {
        throw std::runtime_error(
            "Error: Reaction Style must be GeneralPowerLaw, LangmuirHinshelwood, or "
            "LangmuirHinshelwoodHougenWatson (" +
            context + ")");
      }

      reaction.reactants = readReactionParticipants(requireKeyCaseInsensitive(value, "Reactants", context), components,
                                                    "Reactants", context);
      reaction.products = readReactionParticipants(requireKeyCaseInsensitive(value, "Products", context), components,
                                                   "Products", context);

      if (reaction.reactants.empty() || reaction.products.empty())
      {
        throw std::runtime_error("Error: Reaction requires at least one reactant and one product (" + context + ")");
      }

      std::set<size_t> participants;
      for (size_t comp : reaction.reactants)
      {
        if (comp >= components.size())
        {
          throw std::runtime_error("Error: Reaction reactant component index is out of range (" + context + ")");
        }
        if (!participants.insert(comp).second)
        {
          throw std::runtime_error("Error: duplicate component index in Reaction (" + context + ")");
        }
      }
      for (size_t comp : reaction.products)
      {
        if (comp >= components.size())
        {
          throw std::runtime_error("Error: Reaction product component index is out of range (" + context + ")");
        }
        if (!participants.insert(comp).second)
        {
          throw std::runtime_error("Error: duplicate component index in Reaction (" + context + ")");
        }
      }

      const size_t participantCount = reaction.reactants.size() + reaction.products.size();
      const nlohmann::json* stoichiometryValue = findKeyCaseInsensitive(value, "Stoichiometry");
      std::vector<double> stoichiometry;
      if (stoichiometryValue == nullptr)
      {
        stoichiometry.assign(participantCount, 1.0);
      }
      else if (stoichiometryValue->is_array())
      {
        stoichiometry = readSizedDoubleList(*stoichiometryValue, "Stoichiometry", context, participantCount);
      }
      else if (stoichiometryValue->is_object())
      {
        requireOnlyKnownKeys(*stoichiometryValue, {"Reactants", "Products"}, context + " Stoichiometry");
        const std::vector<double> reactantValues =
            containsKeyCaseInsensitive(*stoichiometryValue, "Reactants")
                ? readSizedDoubleList(requireKeyCaseInsensitive(*stoichiometryValue, "Reactants", context), "Reactants",
                                      context + " Stoichiometry", reaction.reactants.size())
                : std::vector<double>(reaction.reactants.size(), 1.0);
        const std::vector<double> productValues =
            containsKeyCaseInsensitive(*stoichiometryValue, "Products")
                ? readSizedDoubleList(requireKeyCaseInsensitive(*stoichiometryValue, "Products", context), "Products",
                                      context + " Stoichiometry", reaction.products.size())
                : std::vector<double>(reaction.products.size(), 1.0);

        stoichiometry = reactantValues;
        stoichiometry.insert(stoichiometry.end(), productValues.begin(), productValues.end());
      }
      else
      {
        throw std::runtime_error("Error: Stoichiometry must be an array or object (" + context + ")");
      }
      reaction.reactantStoichiometry.assign(
          stoichiometry.begin(), stoichiometry.begin() + static_cast<std::ptrdiff_t>(reaction.reactants.size()));
      reaction.productStoichiometry.assign(
          stoichiometry.begin() + static_cast<std::ptrdiff_t>(reaction.reactants.size()), stoichiometry.end());

      reaction.forwardOrders = containsKeyCaseInsensitive(value, "ForwardOrders")
                                   ? readSizedDoubleList(requireKeyCaseInsensitive(value, "ForwardOrders", context),
                                                         "ForwardOrders", context, reaction.reactants.size())
                                   : reaction.reactantStoichiometry;
      reaction.backwardOrders = containsKeyCaseInsensitive(value, "BackwardOrders")
                                    ? readSizedDoubleList(requireKeyCaseInsensitive(value, "BackwardOrders", context),
                                                          "BackwardOrders", context, reaction.products.size())
                                    : reaction.productStoichiometry;
      readOptionalNonNegativeInteger(value, "Site", reaction.site);
      readOptionalNumber<double>(value, "RateLimitTime", reaction.rateLimitTime);

      if (!std::isfinite(reaction.rateLimitTime) || reaction.rateLimitTime < 0.0)
      {
        throw std::runtime_error("Error: RateLimitTime must be non-negative (" + context + ")");
      }
      if ((reaction.style == Reaction::Style::LangmuirHinshelwood ||
           reaction.style == Reaction::Style::LangmuirHinshelwoodHougenWatson) &&
          reaction.phase != Reaction::Phase::PoreConcentration)
      {
        throw std::runtime_error("Error: LH/LHHW reactions require Phase=PoreConcentration (" + context + ")");
      }

      if (reaction.phase == Reaction::Phase::Chemisorbed && !containsKeyCaseInsensitive(parsed_data, "Adsorbents"))
      {
        for (size_t comp : participants)
        {
          if (reaction.site >= components[comp].chemisorption.numberOfSites)
          {
            throw std::runtime_error("Error: Chemisorbed reaction component has no requested chemisorption site (" +
                                     context + ")");
          }
        }
      }

      const nlohmann::json& kinetics = requireKeyCaseInsensitive(value, "Kinetics", context);
      if (!kinetics.is_object())
      {
        throw std::runtime_error("Error: Reaction Kinetics must be an object (" + context + ")");
      }
      requireOnlyKnownKeys(kinetics,
                           {
                               "forwardRateCoefficient",
                               "forwardActivationEnergy",
                               "equilibriumConstant",
                               "gibbsFreeEnergy",
                           },
                           context + " Kinetics");

      reaction.forwardRateCoefficient =
          getNumberOrThrow<double>(requireKeyCaseInsensitive(kinetics, "forwardRateCoefficient", context + " Kinetics"),
                                   "forwardRateCoefficient", context + " Kinetics");
      reaction.forwardActivationEnergy = getNumberOrThrow<double>(
          requireKeyCaseInsensitive(kinetics, "forwardActivationEnergy", context + " Kinetics"),
          "forwardActivationEnergy", context + " Kinetics");
      reaction.equilibriumConstant =
          getNumberOrThrow<double>(requireKeyCaseInsensitive(kinetics, "equilibriumConstant", context + " Kinetics"),
                                   "equilibriumConstant", context + " Kinetics");
      reaction.gibbsFreeEnergy =
          getNumberOrThrow<double>(requireKeyCaseInsensitive(kinetics, "gibbsFreeEnergy", context + " Kinetics"),
                                   "gibbsFreeEnergy", context + " Kinetics");

      if (!std::isfinite(reaction.forwardRateCoefficient) || reaction.forwardRateCoefficient < 0.0)
      {
        throw std::runtime_error("Error: forwardRateCoefficient must be non-negative (" + context + ")");
      }
      if (!std::isfinite(reaction.forwardActivationEnergy))
      {
        throw std::runtime_error("Error: forwardActivationEnergy must be finite (" + context + ")");
      }
      if (!std::isfinite(reaction.equilibriumConstant) || reaction.equilibriumConstant <= 0.0)
      {
        throw std::runtime_error("Error: equilibriumConstant must be positive (" + context + ")");
      }
      if (!std::isfinite(reaction.gibbsFreeEnergy))
      {
        throw std::runtime_error("Error: gibbsFreeEnergy must be finite (" + context + ")");
      }
      reactions.push_back(std::move(reaction));
    }
  }

  adsorbentComponents.clear();
  adsorbentLengths.clear();
  adsorbentInterfaceLengths.clear();
  adsorbentMixFractions.clear();
  adsorbentGridPoints.clear();
  adsorbentVoidFractions.clear();
  adsorbentParticleDensities.clear();
  adsorbentParticleDiameters.clear();
  adsorbentGeometries.clear();

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
    std::vector<bool> hasMixFraction;
    hasMixFraction.reserve(adsorbentsJson.size());

    for (std::size_t ads = 0; ads < adsorbentsJson.size(); ++ads)
    {
      const nlohmann::json& adsorbent = adsorbentsJson[ads];
      const std::string context = "Adsorbent " + std::to_string(ads);
      if (!adsorbent.is_object())
      {
        throw std::runtime_error("Error: each adsorbent entry must be an object (" + context + ")");
      }
      requireOnlyKnownKeys(adsorbent,
                           {"Name", "ParticleDiameter", "ColumnVoidFraction", "ParticleDensity", "AdsorbentLength",
                            "MixFraction", "ComponentParameters", "Components"},
                           context);

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
      double adsorbentMixFraction = 0.0;
      readOptionalNumber<double>(adsorbent, "ParticleDiameter", adsorbentParticleDiameter);
      readOptionalNumber<double>(adsorbent, "ColumnVoidFraction", adsorbentVoidFraction);
      readOptionalNumber<double>(adsorbent, "ParticleDensity", adsorbentParticleDensity);
      readOptionalNumber<double>(adsorbent, "AdsorbentLength", adsorbentLength);
      const bool hasAdsorbentMixFraction = containsKeyCaseInsensitive(adsorbent, "MixFraction");
      readOptionalNumber<double>(adsorbent, "MixFraction", adsorbentMixFraction);

      adsorbentParticleDiameters.push_back(adsorbentParticleDiameter);
      adsorbentVoidFractions.push_back(adsorbentVoidFraction);
      adsorbentParticleDensities.push_back(adsorbentParticleDensity);
      adsorbentLengths.push_back(adsorbentLength);
      adsorbentMixFractions.push_back(adsorbentMixFraction);
      hasMixFraction.push_back(hasAdsorbentMixFraction);
      adsorbentGeometries.push_back(
          geometry.kind == GeometryKind::PackedBed
              ? makeGeometry(PackedBedTubeSpec{.voidFraction = adsorbentVoidFraction,
                                               .particleDiameter = adsorbentParticleDiameter,
                                               .internalDiameter = geometry.dimensions.internalDiameter,
                                               .outerDiameter = geometry.dimensions.outerDiameter})
              : geometry);

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
            throw std::runtime_error("Error: unknown component '" + it.key() + "' in ComponentParameters (" + context +
                                     ")");
          }

          Component& comp = componentsForAdsorbent[componentIndex->second];
          const nlohmann::json& params = it.value();
          const std::string componentContext = context + ", ComponentParameters " + it.key();
          if (!params.is_object())
          {
            throw std::runtime_error("Error: component parameters must be an object (" + componentContext + ")");
          }

          requireOnlyKnownKeys(params,
                               {"Name", "MassTransferCoefficient", "AxialDispersionCoefficient", "HeatOfAdsorption",
                                "referenceTemperature", "nonIsothermal", "ChemisorptionSites", "PhysisorptionSites"},
                               componentContext);

          readOptionalNumber<double>(params, "MassTransferCoefficient", comp.massTransferCoefficient);
          readOptionalNumber<double>(params, "AxialDispersionCoefficient", comp.axialDispersionCoefficient);
          readOptionalNumber<double>(params, "HeatOfAdsorption", comp.heatOfAdsorption);
          readOptionalNumber<double>(params, "referenceTemperature", comp.referenceTemperature);
          readOptionalBool(params, "nonIsothermal", comp.nonIsothermal);
          readChemisorption(comp.chemisorption, params, comp.nonIsothermal, componentContext);

          const bool hasIsothermOverride = containsKeyCaseInsensitive(params, "PhysisorptionSites");

          if (hasIsothermOverride)
          {
            comp.isotherm = MultiSiteIsotherm{};
          }
          readPhysisorptionSites(comp, params, componentContext);
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

    const bool usesUniformAdsorbentMix =
        std::any_of(hasMixFraction.begin(), hasMixFraction.end(), [](bool present) { return present; });
    if (!usesUniformAdsorbentMix)
    {
      adsorbentMixFractions.clear();
    }
    if (usesUniformAdsorbentMix)
    {
      if (!std::all_of(hasMixFraction.begin(), hasMixFraction.end(), [](bool present) { return present; }))
      {
        throw std::runtime_error("Error: MixFraction must be set for every adsorbent in a uniform adsorbent mix");
      }
      if (containsKeyCaseInsensitive(parsed_data, "ColumnSections"))
      {
        throw std::runtime_error("Error: MixFraction cannot be combined with ColumnSections");
      }

      double fractionSum = 0.0;
      for (std::size_t ads = 0; ads < adsorbentMixFractions.size(); ++ads)
      {
        const double fraction = adsorbentMixFractions[ads];
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
        {
          throw std::runtime_error("Error: MixFraction must be between 0 and 1 for adsorbent '" +
                                   adsorbentNames[ads] + "'");
        }
        if (containsKeyCaseInsensitive(adsorbentsJson[ads], "AdsorbentLength"))
        {
          throw std::runtime_error("Error: MixFraction cannot be combined with AdsorbentLength for adsorbent '" +
                                   adsorbentNames[ads] + "'");
        }
        fractionSum += fraction;
      }
      if (std::abs(fractionSum - 1.0) > 1.0e-12)
      {
        throw std::runtime_error("Error: adsorbent MixFraction values must sum to 1");
      }

      adsorbentLengths.clear();
      adsorbentInterfaceLengths.clear();
    }
    else if (containsKeyCaseInsensitive(parsed_data, "ColumnSections"))
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
        requireOnlyKnownKeys(section, {"Adsorbent", "Length", "InterfaceLength", "NumberOfGridPoints"}, context);

        const bool hasAdsorbent = containsKeyCaseInsensitive(section, "Adsorbent");
        const bool hasInterface = containsKeyCaseInsensitive(section, "InterfaceLength");
        if (hasAdsorbent == hasInterface)
        {
          throw std::runtime_error(
              "Error: ColumnSections entries must define either Adsorbent/Length or "
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
        throw std::runtime_error(
            "Error: ColumnSections must define one InterfaceLength between neighboring adsorbents");
      }

      adsorbentLengths = std::move(sectionLengths);
      adsorbentGridPoints = std::move(sectionGridPoints);
      adsorbentInterfaceLengths = std::move(sectionInterfaceLengths);
    }
    else
    {
      adsorbentInterfaceLengths.assign(adsorbentsJson.size() - 1, 0.0);
    }

    for (std::size_t ads = 0; ads < adsorbentsJson.size(); ++ads)
    {
      if (!usesUniformAdsorbentMix && adsorbentLengths[ads] <= 0.0)
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

    if (!usesUniformAdsorbentMix && adsorbentInterfaceLengths.size() != adsorbentsJson.size() - 1)
    {
      throw std::runtime_error("Error: ColumnSections must define Adsorbents size minus one interfaces");
    }

    if (!usesUniformAdsorbentMix)
    {
      columnLength = std::reduce(adsorbentLengths.begin(), adsorbentLengths.end(), 0.0);
    }

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
    adsorbentMixFractions.clear();
    adsorbentGridPoints.push_back(numberOfGridPoints);
    adsorbentVoidFractions.push_back(columnVoidFraction);
    adsorbentParticleDensities.push_back(particleDensity);
    adsorbentParticleDiameters.push_back(particleDiameter);
    adsorbentGeometries.push_back(geometry);
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
      if (mixturePredictionMethod != 6 && !comp.isCarrierGas && comp.isotherm.sites.empty())
      {
        throw std::runtime_error("Error: non-carrier component '" + comp.name +
                                 "' has no isotherm for at least one adsorbent");
      }
      for (const Isotherm& isotherm : comp.isotherm.sites)
      {
        if (isotherm.type != Isotherm::Type::Langmuir_pH) continue;
        if (fluidPhase != 1)
        {
          throw std::runtime_error("Error: pH-Langmuir is available only for liquid simulations");
        }
        if (mixturePredictionMethod != 5)
        {
          throw std::runtime_error("Error: pH-Langmuir requires MixturePredictionMethod SPI");
        }
      }
    }
  }

  if (fluidPhase == 1 && mixturePredictionMethod == 6)
  {
    throw std::runtime_error("Error: MPD mixture prediction is available only for gas simulations");
  }

  for (size_t reactionId = 0; reactionId < reactions.size(); ++reactionId)
  {
    const Reaction& reaction = reactions[reactionId];
    if (reaction.phase != Reaction::Phase::Chemisorbed) continue;

    for (size_t ads = 0; ads < adsorbentComponents.size(); ++ads)
    {
      for (size_t comp : reaction.reactants)
      {
        if (reaction.site >= adsorbentComponents[ads][comp].chemisorption.numberOfSites)
        {
          throw std::runtime_error("Error: Chemisorbed reaction " + std::to_string(reactionId) +
                                   " component has no requested chemisorption site in adsorbent " +
                                   std::to_string(ads));
        }
      }
      for (size_t comp : reaction.products)
      {
        if (reaction.site >= adsorbentComponents[ads][comp].chemisorption.numberOfSites)
        {
          throw std::runtime_error("Error: Chemisorbed reaction " + std::to_string(reactionId) +
                                   " component has no requested chemisorption site in adsorbent " +
                                   std::to_string(ads));
        }
      }
    }
  }

  if ((mixturePredictionMethod == 2) || (mixturePredictionMethod == 3))
  {
    for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
    {
      for (size_t i = 0; i < componentsForAdsorbent.size(); ++i)
      {
        for (size_t j = 0; j < componentsForAdsorbent[i].isotherm.sites.size(); ++j)
        {
          if (componentsForAdsorbent[i].isotherm.sites[j].type != Isotherm::Type::Langmuir)
          {
            throw std::runtime_error("Error: Explicit mixture prediction must use single Langmuir isotherms");
          }
        }
      }
    }
  }

  if (mixturePredictionMethod == 4)
  {
    for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
    {
      size_t maximumSites = 0;
      for (const Component& component : componentsForAdsorbent)
      {
        maximumSites = std::max(maximumSites, component.isotherm.sites.size());
      }

      for (size_t site = 0; site < maximumSites; ++site)
      {
        std::optional<Isotherm::Type> model;
        for (const Component& component : componentsForAdsorbent)
        {
          if (component.isCarrierGas || site >= component.isotherm.sites.size() ||
              !component.isotherm.sites[site].enabled())
            continue;

          const Isotherm::Type type = component.isotherm.sites[site].type;
          if (!supportsSegregatedCompetitiveIsotherm(type))
          {
            throw std::runtime_error("Error: unsupported physisorption isotherm model for SCI mixture prediction");
          }
          if (model.has_value() && *model != type)
          {
            throw std::runtime_error("Error: SCI requires one physisorption isotherm model per segregated site");
          }
          model = type;
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
            throw std::runtime_error("Error: SEI chemisorption mixture prediction requires Langmuir isotherms");
          }
        }
      }
    }
  }

  if (mixturePredictionMethod == 4)
  {
    for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
    {
      size_t maximumSites = 0;
      for (const Component& component : componentsForAdsorbent)
      {
        maximumSites = std::max(maximumSites, component.chemisorption.sites.size());
      }

      for (size_t site = 0; site < maximumSites; ++site)
      {
        std::optional<Isotherm::Type> model;
        for (const Component& component : componentsForAdsorbent)
        {
          if (site >= component.chemisorption.sites.size() ||
              !component.chemisorption.sites[site].isotherm.has_value() ||
              !component.chemisorption.sites[site].isotherm->enabled())
            continue;

          const Isotherm::Type type = component.chemisorption.sites[site].isotherm->type;
          if (!supportsSegregatedCompetitiveIsotherm(type))
          {
            throw std::runtime_error("Error: unsupported chemisorption isotherm model for SCI mixture prediction");
          }
          if (model.has_value() && *model != type)
          {
            throw std::runtime_error("Error: SCI requires one chemisorption isotherm model per segregated site");
          }
          model = type;
        }
      }
    }
  }
  if (hasChemisorption && mixturePredictionMethod == 2)
  {
    throw std::runtime_error(
        "Error: chemisorption mixture prediction supports IAST, SIAST, SEI, SCI, or SPI; EI is "
        "not supported");
  }

  maxIsothermTerms = 0;
  for (const std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
  {
    if (!componentsForAdsorbent.empty())
    {
      std::vector<Component>::const_iterator maxIsothermTermsIterator = std::max_element(
          componentsForAdsorbent.begin(), componentsForAdsorbent.end(), [](const Component& lhs, const Component& rhs)
          { return lhs.isotherm.sites.size() < rhs.isotherm.sites.size(); });
      maxIsothermTerms = std::max(maxIsothermTerms, maxIsothermTermsIterator->isotherm.sites.size());
    }
  }
  if (mixturePredictionMethod == 6) maxIsothermTerms = std::max(maxIsothermTerms, size_t{1});

  if (simulationType == SimulationType::Breakthrough || simulationType == SimulationType::SwingAdsorption)
  {
    static constexpr size_t inletPressureInletVelocity = 0;
    static constexpr size_t inletPressureOutletPressure = 1;
    static constexpr size_t inletVelocityOutletPressure = 2;
    static constexpr size_t fixedVelocity = 3;
    static constexpr size_t fixedPressureInletVelocity = 4;

    // Boundary conditions define which externally supplied quantities are mandatory.
    const bool boundaryNeedsInletPressure = boundaryCondition == inletPressureInletVelocity ||
                                            boundaryCondition == inletPressureOutletPressure ||
                                            boundaryCondition == fixedPressureInletVelocity;
    const bool boundaryNeedsOutletPressure =
        boundaryCondition == inletPressureOutletPressure || boundaryCondition == inletVelocityOutletPressure;
    const bool boundaryNeedsInletVelocity =
        boundaryCondition == inletPressureInletVelocity || boundaryCondition == inletVelocityOutletPressure ||
        boundaryCondition == fixedVelocity;
    const bool effectiveHasInletPressure =
        hasInletPressure || (simulationType == SimulationType::SwingAdsorption && !swingAdsorptionPhases.empty() &&
                             swingAdsorptionPhases.front().inletPressure.has_value());

    if (fluidPhase == 0)
    {
      requireConfigured(numberOfCarrierGases != 0, "Error: no carrier gas component present");
      requireConfigured(numberOfCarrierGases == 1,
                        "Error: multiple carrier gas component present (there can be only one)");
    }
    else
    {
      requireConfigured(numberOfCarrierGases == 0, "Error: liquid simulations do not use CarrierGas components");
      requireConfigured(std::isfinite(pHValue), "Error: pHValue must be finite");
      requireConfigured(std::isfinite(pKw), "Error: pKw must be finite");
    }

    requireConfigured(temperature >= 0.0, "Error: temperature not set (Use e.g.: 'Temperature 300')");
    requireConfigured(columnVoidFraction > 0.0,
                      "Error: void-fraction of the column not set or invalid (Use e.g.: 'ColumnVoidFraction 0.4')");
    requireConfigured(particleDensity >= 0.0, "Error: particle density not set (Use e.g.: 'ParticleDensity 1408.2')");
    if (simulationType == SimulationType::Breakthrough)
    {
      requireConfigured(numberOfTimeSteps != 0 || autoNumberOfTimeSteps,
                        "Error: number of time steps not set (Use e.g.: 'NumberOfTimeSteps 5000000')");
    }
    else
    {
      requireConfigured(!swingAdsorptionPhases.empty(), "Error: SwingAdsorption requires SwingAdsorptionPhases");
      for (size_t phaseIndex = 0; phaseIndex < swingAdsorptionPhases.size(); ++phaseIndex)
      {
        const SwingAdsorptionPhase& phase = swingAdsorptionPhases[phaseIndex];
        requireConfigured(phase.numberOfSteps != 0, "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) +
                                                        " must have a positive number of time steps");
        if (phase.temperature.has_value())
        {
          requireConfigured(
              std::isfinite(*phase.temperature) && *phase.temperature >= 0.0,
              "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) + " has invalid temperature");
        }
        if (phase.inletPressure.has_value())
        {
          requireConfigured(
              std::isfinite(*phase.inletPressure) && *phase.inletPressure > 0.0,
              "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) + " has invalid inlet pressure");
          if (boundaryCondition == fixedPressureInletVelocity)
          {
            requireConfigured(*phase.inletPressure + pressureGradient * columnLength > 0.0,
                              "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) +
                                  " fixed pressure profile becomes non-positive");
          }
        }
      }
    }
    requireConfigured(numberOfGridPoints != 0,
                      "Error: number of grid points not set (Use e.g.: 'NumberOfGridPoints 50')");
    requireConfigured(columnLength > 0.0, "Error: column length not set or invalid (Use e.g.: 'ColumnLength 0.3')");

    if (boundaryNeedsInletPressure)
    {
      requireConfigured(effectiveHasInletPressure && inletPressure >= 0.0,
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
      requireConfigured(
          (effectiveHasInletPressure && inletPressure >= 0.0) || (hasOutletPressure && outletPressure >= 0.0),
          "Error: FixedVelocity requires InletPressure or OutletPressure");
    }

    if (boundaryCondition == fixedPressureInletVelocity)
    {
      requireConfigured(hasPressureGradient, "Error: FixedPressureInletVelocity requires PressureGradient");
      requireConfigured(inletPressure + pressureGradient * columnLength > 0.0,
                        "Error: fixed pressure profile becomes non-positive");
    }
  }
}

bool InputReader::isMultibed() const
{
  return debugForceMultibed || adsorbentComponents.size() > 1;
}

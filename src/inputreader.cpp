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
#include <iterator>
#include <numeric>
#include <optional>
#include <print>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <tuple>
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
                                 const std::string& context, bool allowIsothermTypeKeys = false);

/**
 * \brief Requires that object holds exactly the keys requiredKeys, compared without regard to case as every key.
 */
static void requireExactKeys(const nlohmann::json& object, std::initializer_list<std::string_view> requiredKeys,
                             const std::string& context)
{
  requireOnlyKnownKeys(object, requiredKeys, context);
  for (const std::string_view key : requiredKeys)
  {
    if (!containsKeyCaseInsensitive(object, std::string{key}))
    {
      throw std::runtime_error("Error: required key '" + std::string{key} + "' missing (" + context + ")");
    }
  }
}

struct IsothermSpec
{
  Isotherm::Type type;
  std::size_t parameterCount;
};

static const IsothermSpec* findIsothermSpec(const std::string& typeString);

/**
 * \brief Reads a temperature-scalable number: a plain number, or a temperature law object
 *        {"Value", "Law", <coefficient key>, "ReferenceTemperature"} (TemperatureLaw).
 *
 * \param value The JSON value at the parameter position.
 * \param spec The laws accepted at this position; an empty mask accepts plain numbers only. With spec.heatLinked, a
 *        VantHoff law may omit "Energy": it is then marked energyInherited and takes the heat in linkHeat.
 * \param referenceTemperature T_ref of the enclosing object, used when the law gives none.
 * \param context Location of the value, for error messages.
 */
static TemperatureLaw readTemperatureLaw(const nlohmann::json& value, const TemperatureParameterSpec& spec,
                                         std::optional<double> referenceTemperature, const std::string& context)
{
  const std::string name{spec.name};
  if (value.is_number())
  {
    const double number = getNumberOrThrow<double>(value, name, context);
    if (!std::isfinite(number)) throw std::runtime_error("Error: the " + name + " must be finite (" + context + ")");
    return TemperatureLaw(number);
  }
  if (!value.is_object())
  {
    throw std::runtime_error("Error: the " + name + " must be a number or a temperature law object (" + context + ")");
  }
  if (spec.allowedLaws == 0)
  {
    throw std::runtime_error("Error: the " + name + " cannot depend on temperature; give a number (" + context + ")");
  }

  const std::string lawName = getStringOrThrow(requireKeyCaseInsensitive(value, "Law", context), "Law", context);
  const std::optional<TemperatureLaw::Type> type = TemperatureLaw::fromName(lawName);
  if (!type.has_value() ||
      (*type != TemperatureLaw::Type::Constant && (spec.allowedLaws & TemperatureLaw::mask(*type)) == 0))
  {
    throw std::runtime_error(std::format("Error: the {} cannot use the temperature law '{}'; use {} ({})", name,
                                         lawName, TemperatureLaw::describeMask(spec.allowedLaws), context));
  }

  const std::string coefficientKey{TemperatureLaw::coefficientKey(*type)};
  for (auto it = value.begin(); it != value.end(); ++it)
  {
    const bool known = caseInSensStringCompare(it.key(), "Value") || caseInSensStringCompare(it.key(), "Law") ||
                       caseInSensStringCompare(it.key(), "ReferenceTemperature") ||
                       (!coefficientKey.empty() && caseInSensStringCompare(it.key(), coefficientKey));
    if (!known)
    {
      throw std::runtime_error(std::format("Error: unknown key '{}' in the {} law of the {}; it takes Value{}{} and "
                                           "ReferenceTemperature ({})",
                                           it.key(), lawName, name, coefficientKey.empty() ? "" : ", ",
                                           coefficientKey, context));
    }
  }

  TemperatureLaw law;
  law.type = *type;
  law.value = getNumberOrThrow<double>(requireKeyCaseInsensitive(value, "Value", context), "Value", context);
  if (!std::isfinite(law.value))
  {
    throw std::runtime_error("Error: the Value of the " + name + " must be finite (" + context + ")");
  }
  readOptionalNumber<double>(value, "ReferenceTemperature", law.referenceTemperature);
  if (!law.referenceTemperature.has_value()) law.referenceTemperature = referenceTemperature;
  if (law.referenceTemperature.has_value() &&
      !(std::isfinite(*law.referenceTemperature) && *law.referenceTemperature > 0.0))
  {
    throw std::runtime_error("Error: the ReferenceTemperature of the " + name + " must be positive (" + context +
                             ")");
  }
  if (law.type == TemperatureLaw::Type::Constant) return law;

  const nlohmann::json* coefficient = findKeyCaseInsensitive(value, coefficientKey);
  if (coefficient != nullptr)
  {
    law.coefficient = getNumberOrThrow<double>(*coefficient, coefficientKey, context);
    if (!std::isfinite(law.coefficient))
    {
      throw std::runtime_error("Error: " + coefficientKey + " of the " + name + " must be finite (" + context + ")");
    }
  }
  else if (law.type == TemperatureLaw::Type::VantHoff && spec.heatLinked)
  {
    law.energyInherited = true;
  }
  else
  {
    throw std::runtime_error(
        std::format("Error: the {} law of the {} needs \"{}\" ({})", lawName, name, coefficientKey, context));
  }

  if (TemperatureLaw::requiresReferenceTemperature(law.type) && !law.referenceTemperature.has_value())
  {
    throw std::runtime_error(std::format("Error: the {} law of the {} needs a ReferenceTemperature, in the law or in "
                                         "the enclosing component or Kinetics ({})",
                                         lawName, name, context));
  }
  return law;
}

/**
 * \brief Parses the shared "Type"/"Parameters" isotherm pair used by both PhysisorptionSites and
 *        ChemisorptionSites entries.
 *
 * Each entry of Parameters is a number or, at the positions of Isotherm::temperatureParameters, a temperature law.
 * \param referenceTemperature The component's ReferenceTemperature, used by laws that give none.
 */
static Isotherm parseSiteIsotherm(const nlohmann::json& site, std::optional<double> referenceTemperature,
                                  const std::string& siteContext)
{
  const std::string typeString =
      getStringOrThrow(requireKeyCaseInsensitive(site, "Type", siteContext), "Type", siteContext);
  const nlohmann::json& params = requireKeyCaseInsensitive(site, "Parameters", siteContext);
  const IsothermSpec* spec = findIsothermSpec(typeString);
  if (spec == nullptr)
  {
    throw std::runtime_error("Error: unknown isotherm type '" + typeString + "' (" + siteContext + ")");
  }
  if (!params.is_array())
  {
    throw std::runtime_error("Error: key '" + typeString + "' must be an array of numbers (" + siteContext + ")");
  }

  const std::span<const TemperatureParameterSpec> temperatureParameters = Isotherm::temperatureParameters(spec->type);
  std::vector<double> values;
  std::vector<std::pair<size_t, TemperatureLaw>> temperatureLaws;
  for (size_t index = 0; index < params.size(); ++index)
  {
    const TemperatureParameterSpec parameter = index < temperatureParameters.size()
                                                   ? temperatureParameters[index]
                                                   : TemperatureParameterSpec{"parameter", 0, false};
    const std::string parameterContext = std::format("{}, {} parameter {}", siteContext, typeString, index);
    TemperatureLaw law = readTemperatureLaw(params[index], parameter, referenceTemperature, parameterContext);
    values.push_back(law.value);
    if (!law.isConstant()) temperatureLaws.emplace_back(index, law);
  }
  values = requireDoubleParameterCount(std::move(values), spec->parameterCount, typeString, siteContext);

  return Isotherm(spec->type, values, std::move(temperatureLaws));
}

/**
 * \brief Links the heat-linked VantHoff laws (the affinities) of isotherm sites to a heat of adsorption.
 *
 * A law without "Energy" (energyInherited) takes inheritedHeat; it is an error if there is none. Returns the
 * distinct energies of the other heat-linked VantHoff laws, from which the heat follows when the input gives none.
 */
static std::vector<double> linkHeat(std::span<Isotherm> sites, std::optional<double> inheritedHeat,
                                    const std::string& context)
{
  std::vector<double> energies;
  for (Isotherm& site : sites)
  {
    const std::span<const TemperatureParameterSpec> parameters = Isotherm::temperatureParameters(site.type);
    for (auto& [index, law] : site.temperatureLaws)
    {
      if (index >= parameters.size() || !parameters[index].heatLinked ||
          law.type != TemperatureLaw::Type::VantHoff)
      {
        continue;
      }
      if (law.energyInherited)
      {
        if (!inheritedHeat.has_value())
        {
          throw std::runtime_error(std::format(
              "Error: the VantHoff law of the {} has no \"Energy\", and there is no HeatOfAdsorption to take it "
              "from ({})",
              parameters[index].name, context));
        }
        law.coefficient = *inheritedHeat;
      }
      else if (std::find(energies.begin(), energies.end(), law.coefficient) == energies.end())
      {
        energies.push_back(law.coefficient);
      }
    }
  }
  return energies;
}

/**
 * \brief Reads the ChemisorptionSites of a component.
 *
 * \param defaultHeatOfAdsorption The component's HeatOfAdsorption, if given: the default heat of each site.
 * \param referenceTemperature The component's ReferenceTemperature, used by temperature laws that give none.
 */
static void readChemisorption(MultiSiteChemisorption& chemisorption, const nlohmann::json& object,
                              double defaultMassTransferCoefficient, std::optional<double> defaultHeatOfAdsorption,
                              std::optional<double> referenceTemperature, const std::string& context)
{
  constexpr unsigned arrhenius = TemperatureLaw::mask(TemperatureLaw::Type::Arrhenius);
  const auto readRateCoefficient = [&](const nlohmann::json& parameters, std::string_view key, const std::string& where)
  {
    return readTemperatureLaw(requireKeyCaseInsensitive(parameters, std::string{key}, where),
                              TemperatureParameterSpec{key, arrhenius, false}, referenceTemperature,
                              where + ", " + std::string{key});
  };

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

      requireOnlyKnownKeys(
          value,
          {"Type", "Parameters", "RateEquation", "RateEquationParameters", "HeatOfAdsorption", "MaximumLoading"},
          siteContext, false);

      Chemisorption siteChemisorption;
      siteChemisorption.isotherm = parseSiteIsotherm(value, referenceTemperature, siteContext);

      std::string rateEquationString = "FirstOrder";
      readOptionalString(value, "RateEquation", rateEquationString);

      if (caseInSensStringCompare(rateEquationString, "None"))
        siteChemisorption.type = Chemisorption::Type::None;
      else if (caseInSensStringCompare(rateEquationString, "FirstOrder"))
        siteChemisorption.type = Chemisorption::Type::FirstOrder;
      else if (caseInSensStringCompare(rateEquationString, "PseudoNth"))
        siteChemisorption.type = Chemisorption::Type::PseudoNth;
      else if (caseInSensStringCompare(rateEquationString, "Avrami"))
        siteChemisorption.type = Chemisorption::Type::Avrami;
      else if (caseInSensStringCompare(rateEquationString, "General"))
        siteChemisorption.type = Chemisorption::Type::General;
      else if (caseInSensStringCompare(rateEquationString, "Elovich"))
        siteChemisorption.type = Chemisorption::Type::Elovich;
      else
      {
        throw std::runtime_error("Error: invalid RateEquation '" + rateEquationString + "' (" + siteContext + ")");
      }

      if (siteChemisorption.type == Chemisorption::Type::None)
      {
        throw std::runtime_error("Error: RateEquation 'None' is not valid for a ChemisorptionSites entry (" +
                                 siteContext + ")");
      }

      const nlohmann::json* rateParams = findKeyCaseInsensitive(value, "RateEquationParameters");
      const std::string parametersContext = siteContext + " RateEquationParameters";
      if (rateParams != nullptr && !rateParams->is_object())
      {
        throw std::runtime_error("Error: RateEquationParameters must be an object (" + siteContext + ")");
      }

      switch (siteChemisorption.type)
      {
        case Chemisorption::Type::None:
          break;
        case Chemisorption::Type::FirstOrder:
          if (rateParams != nullptr)
          {
            requireOnlyKnownKeys(*rateParams,
                                 {"MassTransferCoefficient", "DesorptionMassTransferCoefficient",
                                  "DesorptionCoverageEnergy", "CoverageReferenceLoading"},
                                 parametersContext);
            siteChemisorption.rateCoefficient =
                readRateCoefficient(*rateParams, "MassTransferCoefficient", parametersContext);

            // Optional coverage-dependent desorption branch (used when q > q*).
            readOptionalNumber<double>(*rateParams, "DesorptionMassTransferCoefficient",
                                       siteChemisorption.desorptionMassTransferCoefficient);
            readOptionalNumber<double>(*rateParams, "DesorptionCoverageEnergy",
                                       siteChemisorption.desorptionCoverageEnergy);
            readOptionalNumber<double>(*rateParams, "CoverageReferenceLoading",
                                       siteChemisorption.coverageReferenceLoading);
            const bool hasDesorptionBranch =
                containsKeyCaseInsensitive(*rateParams, "DesorptionMassTransferCoefficient") ||
                containsKeyCaseInsensitive(*rateParams, "DesorptionCoverageEnergy") ||
                containsKeyCaseInsensitive(*rateParams, "CoverageReferenceLoading");
            if (hasDesorptionBranch)
            {
              if (!(siteChemisorption.desorptionMassTransferCoefficient > 0.0) ||
                  !std::isfinite(siteChemisorption.desorptionMassTransferCoefficient))
              {
                throw std::runtime_error(
                    "Error: DesorptionMassTransferCoefficient must be positive when a desorption branch is given (" +
                    parametersContext + ")");
              }
              if (!(siteChemisorption.coverageReferenceLoading > 0.0) ||
                  !std::isfinite(siteChemisorption.coverageReferenceLoading))
              {
                throw std::runtime_error(
                    "Error: CoverageReferenceLoading must be positive when a desorption branch is given (" +
                    parametersContext + ")");
              }
              if (!std::isfinite(siteChemisorption.desorptionCoverageEnergy))
              {
                throw std::runtime_error("Error: DesorptionCoverageEnergy must be finite (" + parametersContext + ")");
              }
            }
          }
          else
          {
            siteChemisorption.rateCoefficient = defaultMassTransferCoefficient;
          }
          break;
        case Chemisorption::Type::PseudoNth:
        case Chemisorption::Type::Avrami:
          if (rateParams == nullptr)
          {
            throw std::runtime_error("Error: RateEquationParameters is required for RateEquation '" +
                                     rateEquationString + "' (" + siteContext + ")");
          }
          requireExactKeys(*rateParams, {"RateCoefficient", "Order"}, parametersContext);
          siteChemisorption.rateCoefficient = readRateCoefficient(*rateParams, "RateCoefficient", parametersContext);
          siteChemisorption.order = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(*rateParams, "Order", parametersContext), "Order", parametersContext);
          break;
        case Chemisorption::Type::General:
          if (rateParams == nullptr)
          {
            throw std::runtime_error("Error: RateEquationParameters is required for RateEquation '" +
                                     rateEquationString + "' (" + siteContext + ")");
          }
          requireExactKeys(*rateParams,
                           {"AdsorptionRateCoefficient", "DesorptionRateCoefficient", "PoreConcentrationOrder",
                            "CapacityOrder", "DesorptionOrder", "FilmMassTransferCoefficient", "PoreDiffusivity",
                            "UsePoreSurfaceTransport"},
                           parametersContext);
          siteChemisorption.adsorptionRateCoefficient =
              readRateCoefficient(*rateParams, "AdsorptionRateCoefficient", parametersContext);
          siteChemisorption.desorptionRateCoefficient =
              readRateCoefficient(*rateParams, "DesorptionRateCoefficient", parametersContext);
          readOptionalNonNegativeInteger(*rateParams, "PoreConcentrationOrder",
                                         siteChemisorption.poreConcentrationOrder);
          readOptionalNonNegativeInteger(*rateParams, "CapacityOrder", siteChemisorption.capacityOrder);
          readOptionalNonNegativeInteger(*rateParams, "DesorptionOrder", siteChemisorption.desorptionOrder);
          siteChemisorption.filmMassTransferCoefficient = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(*rateParams, "FilmMassTransferCoefficient", parametersContext),
              "FilmMassTransferCoefficient", parametersContext);
          siteChemisorption.poreDiffusivity = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(*rateParams, "PoreDiffusivity", parametersContext), "PoreDiffusivity",
              parametersContext);
          readOptionalBool(*rateParams, "UsePoreSurfaceTransport", siteChemisorption.usePoreSurfaceTransport);
          break;
        case Chemisorption::Type::Elovich:
          if (rateParams == nullptr)
          {
            throw std::runtime_error("Error: RateEquationParameters is required for RateEquation '" +
                                     rateEquationString + "' (" + siteContext + ")");
          }
          requireExactKeys(*rateParams,
                           {"Alpha", "Beta", "FilmMassTransferCoefficient", "PoreDiffusivity",
                            "UsePoreSurfaceTransport"},
                           parametersContext);
          siteChemisorption.elovichAlpha = readRateCoefficient(*rateParams, "Alpha", parametersContext);
          siteChemisorption.elovichBeta = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(*rateParams, "Beta", parametersContext), "Beta", parametersContext);
          siteChemisorption.filmMassTransferCoefficient = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(*rateParams, "FilmMassTransferCoefficient", parametersContext),
              "FilmMassTransferCoefficient", parametersContext);
          siteChemisorption.poreDiffusivity = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(*rateParams, "PoreDiffusivity", parametersContext), "PoreDiffusivity",
              parametersContext);
          readOptionalBool(*rateParams, "UsePoreSurfaceTransport", siteChemisorption.usePoreSurfaceTransport);
          break;
      }

      if (siteChemisorption.order < 0.0)
      {
        throw std::runtime_error("Error: Order must be non-negative (" + parametersContext + ")");
      }

      // A VantHoff affinity law without Energy takes the site heat, else the component heat. The heat of the site is
      // the given site heat, else the energy of its affinity law, else the component heat.
      std::optional<double> siteHeat;
      readOptionalNumber<double>(value, "HeatOfAdsorption", siteHeat);
      const std::vector<double> energies =
          linkHeat(std::span<Isotherm>(&*siteChemisorption.isotherm, 1),
                   siteHeat.has_value() ? siteHeat : defaultHeatOfAdsorption, siteContext);
      siteChemisorption.heatOfChemisorption =
          siteHeat.value_or(energies.empty() ? defaultHeatOfAdsorption.value_or(0.0) : energies.front());

      siteChemisorption.maximumLoading = getNumberOrThrow<double>(
          requireKeyCaseInsensitive(value, "MaximumLoading", siteContext), "MaximumLoading", siteContext);

      if (siteChemisorption.type == Chemisorption::Type::General)
      {
        if (siteChemisorption.adsorptionRateCoefficient.value < 0.0 ||
            siteChemisorption.desorptionRateCoefficient.value < 0.0 ||
            siteChemisorption.filmMassTransferCoefficient < 0.0 || siteChemisorption.poreDiffusivity < 0.0)
        {
          throw std::runtime_error(
              "Error: General Chemisorption kinetic and transport coefficients must be non-negative (" + siteContext +
              ")");
        }
      }
      else if (siteChemisorption.type == Chemisorption::Type::Elovich)
      {
        if (siteChemisorption.elovichAlpha.value < 0.0 || siteChemisorption.elovichBeta < 0.0 ||
            siteChemisorption.filmMassTransferCoefficient < 0.0 || siteChemisorption.poreDiffusivity < 0.0)
        {
          throw std::runtime_error(
              "Error: Elovich Chemisorption kinetic and transport coefficients must be non-negative (" + siteContext +
              ")");
        }
      }
      else if (siteChemisorption.rateCoefficient.value < 0.0)
      {
        throw std::runtime_error("Error: Chemisorption MassTransferCoefficient/RateCoefficient must be non-negative (" +
                                 siteContext + ")");
      }
      if (siteChemisorption.maximumLoading <= 0.0)
      {
        throw std::runtime_error("Error: MaximumLoading must be positive (" + siteContext + ")");
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
      {"Langmuir", {Isotherm::Type::Langmuir, 2}},
      {"pH-Langmuir", {Isotherm::Type::Langmuir_pH, 3}},
      {"Anti-Langmuir", {Isotherm::Type::Anti_Langmuir, 2}},
      {"BET", {Isotherm::Type::BET, 3}},
      {"Henry", {Isotherm::Type::Henry, 1}},
      {"Freundlich", {Isotherm::Type::Freundlich, 2}},
      {"Sips", {Isotherm::Type::Sips, 3}},
      {"Langmuir-Freundlich", {Isotherm::Type::Langmuir_Freundlich, 3}},
      {"Redlich-Peterson", {Isotherm::Type::Redlich_Peterson, 3}},
      {"Toth", {Isotherm::Type::Toth, 3}},
      {"Unilan", {Isotherm::Type::Unilan, 3}},
      {"OBrien&Myers", {Isotherm::Type::OBrien_Myers, 3}},
      {"Quadratic", {Isotherm::Type::Quadratic, 3}},
      {"Temkin", {Isotherm::Type::Temkin, 3}},
      {"Bingel&Walton", {Isotherm::Type::BingelWalton, 3}},
      {"GAB", {Isotherm::Type::GAB, 3}},
      {"Humidity-enhanced-Langmuir", {Isotherm::Type::Langmuir_humidity, 4}},
      {"Moisture-swing-Langmuir", {Isotherm::Type::Langmuir_moisture_swing, 3}},
      {"GAB-coadsorption-Langmuir", {Isotherm::Type::GAB_Langmuir, 4}},
      {"GAB-coadsorption-Sips", {Isotherm::Type::GAB_Sips, 5}},
      {"GAB-coadsorption-Toth", {Isotherm::Type::GAB_Toth, 5}},
      {"GAB-modified-Toth", {Isotherm::Type::GAB_modified_Toth, 5}},
      {"Exponential-humidity-Toth", {Isotherm::Type::Toth_exponential_humidity, 4}},
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
                                 const std::string& context, bool allowIsothermTypeKeys)
{
  if (!object.is_object())
  {
    throw std::runtime_error("Error: input settings must be a JSON object" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  for (auto it = object.begin(); it != object.end(); ++it)
  {
    // Keys are case-insensitive, so two keys that differ only in case would give one setting twice.
    for (auto other = std::next(it); other != object.end(); ++other)
    {
      if (caseInSensStringCompare(it.key(), other.key()))
      {
        throw std::runtime_error("Error: keys '" + it.key() + "' and '" + other.key() +
                                 "' differ only in case; keys are case-insensitive" +
                                 (context.empty() ? "" : (" (" + context + ")")));
      }
    }

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

/**
 * \brief Reads the PhysisorptionSites of a component.
 *
 * \param referenceTemperature The component's ReferenceTemperature, used by temperature laws that give none.
 * \param heatOfAdsorption Set to the sites' HeatOfAdsorption if they give one (it must be equal for all sites).
 */
static void readPhysisorptionSites(Component& comp, const nlohmann::json& item,
                                   std::optional<double> referenceTemperature,
                                   std::optional<double>& heatOfAdsorption, const std::string& context)
{
  if (!containsKeyCaseInsensitive(item, "PhysisorptionSites")) return;

  const nlohmann::json& sites = requireKeyCaseInsensitive(item, "PhysisorptionSites", context);
  if (!sites.is_array())
  {
    throw std::runtime_error("Error: PhysisorptionSites must be an array (" + context + ")");
  }

  std::optional<double> massTransferCoefficientOverride;
  std::optional<double> heatOfAdsorptionOverride;

  for (std::size_t siteId = 0; siteId < sites.size(); ++siteId)
  {
    const nlohmann::json& site = sites[siteId];
    std::string siteContext = context + ", PhysisorptionSite " + std::to_string(siteId);

    if (!site.is_object())
    {
      throw std::runtime_error("Error: invalid PhysisorptionSites entry (" + siteContext + ")");
    }

    requireOnlyKnownKeys(site, {"Type", "Parameters", "RateEquation", "RateEquationParameters", "HeatOfAdsorption"},
                         siteContext);

    comp.isotherm.add(parseSiteIsotherm(site, referenceTemperature, siteContext));

    if (containsKeyCaseInsensitive(site, "RateEquation"))
    {
      const std::string rateEquation = getStringOrThrow(
          requireKeyCaseInsensitive(site, "RateEquation", siteContext), "RateEquation", siteContext);
      if (!caseInSensStringCompare(rateEquation, "FirstOrder"))
      {
        throw std::runtime_error(
            "Error: RateEquation '" + rateEquation +
            "' is not supported for PhysisorptionSites; only 'FirstOrder' is currently implemented (" + siteContext +
            "). Use ChemisorptionSites for a site with independent kinetics.");
      }
    }

    if (containsKeyCaseInsensitive(site, "RateEquationParameters"))
    {
      const nlohmann::json& params = requireKeyCaseInsensitive(site, "RateEquationParameters", siteContext);
      if (!params.is_object())
      {
        throw std::runtime_error("Error: RateEquationParameters must be an object (" + siteContext + ")");
      }
      requireOnlyKnownKeys(params, {"MassTransferCoefficient"}, siteContext + " RateEquationParameters");

      const double value = getNumberOrThrow<double>(
          requireKeyCaseInsensitive(params, "MassTransferCoefficient", siteContext), "MassTransferCoefficient",
          siteContext);
      if (massTransferCoefficientOverride.has_value() && *massTransferCoefficientOverride != value)
      {
        throw std::runtime_error(
            "Error: conflicting MassTransferCoefficient values across PhysisorptionSites of the same component (" +
            context + ")");
      }
      massTransferCoefficientOverride = value;
    }

    if (containsKeyCaseInsensitive(site, "HeatOfAdsorption"))
    {
      const double value = getNumberOrThrow<double>(
          requireKeyCaseInsensitive(site, "HeatOfAdsorption", siteContext), "HeatOfAdsorption", siteContext);
      if (heatOfAdsorptionOverride.has_value() && *heatOfAdsorptionOverride != value)
      {
        throw std::runtime_error(
            "Error: conflicting HeatOfAdsorption values across PhysisorptionSites of the same component (" + context +
            ")");
      }
      heatOfAdsorptionOverride = value;
    }
  }

  if (massTransferCoefficientOverride.has_value())
  {
    comp.massTransferCoefficient = *massTransferCoefficientOverride;
  }
  if (heatOfAdsorptionOverride.has_value())
  {
    heatOfAdsorption = *heatOfAdsorptionOverride;
  }
}

/**
 * \brief The temperature input of a component as given, from which the per-adsorbent overrides start.
 */
struct ComponentTemperatureInput
{
  std::optional<double> heatOfAdsorption;      ///< HeatOfAdsorption of the component or its physisorption sites.
  std::optional<double> referenceTemperature;  ///< ReferenceTemperature of the component.
  /// The differing energies of the physisorption affinity laws when no heat is given, as text; the energy balance
  /// then has no single heat for the component.
  std::string ambiguousEnergies;
};

/**
 * \brief Parses a component object.
 *
 * \param temperatureInput Set to the heat of adsorption and reference temperature as given.
 */
static Component parseComponentObject(std::size_t componentId, const nlohmann::json& item,
                                      ComponentTemperatureInput& temperatureInput)
{
  if (!item.is_object())
  {
    throw std::runtime_error("Error: each component entry must be an object");
  }

  std::string context = "Component " + std::to_string(componentId);
  requireOnlyKnownKeys(item,
                       {"Name", "FileName", "CarrierGas", "Water", "GasPhaseMolFraction", "LiquidPhaseConcentration",
                        "InitialLiquidPhaseConcentration", "MassTransferCoefficient",
                        "AxialDispersionCoefficient", "MolecularWeight", "HeatOfAdsorption", "ReferenceTemperature",
                        "ChemisorptionSites", "PhysisorptionSites"},
                       context);

  std::string componentName = getStringOrThrow(requireKeyCaseInsensitive(item, "Name", context), "Name", context);

  Component comp(componentId, componentName);

  readOptionalString(item, "FileName", comp.filename);
  readOptionalBool(item, "CarrierGas", comp.isCarrierGas);
  readOptionalBool(item, "Water", comp.isWater);
  readOptionalNumber<double>(item, "GasPhaseMolFraction", comp.initialGasMoleFraction);
  readOptionalNumber<double>(item, "LiquidPhaseConcentration", comp.inletLiquidConcentration);
  readOptionalNumber<double>(item, "InitialLiquidPhaseConcentration", comp.initialLiquidConcentration);
  readOptionalNumber<double>(item, "MassTransferCoefficient", comp.massTransferCoefficient);
  readOptionalNumber<double>(item, "AxialDispersionCoefficient", comp.axialDispersionCoefficient);
  readOptionalNumber<double>(item, "MolecularWeight", comp.molecularWeight);

  temperatureInput = ComponentTemperatureInput{};
  readOptionalNumber<double>(item, "HeatOfAdsorption", temperatureInput.heatOfAdsorption);
  readOptionalNumber<double>(item, "ReferenceTemperature", temperatureInput.referenceTemperature);
  if (temperatureInput.referenceTemperature.has_value() &&
      !(std::isfinite(*temperatureInput.referenceTemperature) && *temperatureInput.referenceTemperature > 0.0))
  {
    throw std::runtime_error("Error: ReferenceTemperature must be positive (" + context + ")");
  }
  readChemisorption(comp.chemisorption, item, comp.massTransferCoefficient, temperatureInput.heatOfAdsorption,
                    temperatureInput.referenceTemperature, context);
  readPhysisorptionSites(comp, item, temperatureInput.referenceTemperature, temperatureInput.heatOfAdsorption,
                         context);

  // The physisorption heat is the given one, else the common energy of the affinity laws (explicit heat wins).
  const std::vector<double> energies = linkHeat(comp.isotherm.sites, temperatureInput.heatOfAdsorption, context);
  comp.heatOfAdsorption = temperatureInput.heatOfAdsorption.value_or(energies.size() == 1 ? energies.front() : 0.0);
  if (!temperatureInput.heatOfAdsorption.has_value() && energies.size() > 1)
  {
    for (double energy : energies)
    {
      temperatureInput.ambiguousEnergies += std::format("{}{} J/mol", temperatureInput.ambiguousEnergies.empty() ? "" : ", ", energy);
    }
  }
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

/**
 * \brief Rejects the step-based time keys of RUPTURA 2 and gives their replacement in s.
 *
 * The replacement follows the same rules as convert.py, which rewrites a whole input file:
 * SimulationTime = NumberOfTimeSteps * TimeStep, WriteInterval = WriteEvery * TimeStep (WriteEvery defaults to
 * 10000), phase Duration = NumberOfTimeSteps * TimeStep (TimeStep defaults to 5e-4 s), PrintEvery is dropped, and
 * TimeStep and NumberOfInitTimeSteps are dropped for the adaptive integrators CVODE and ROS34PW2.
 */
static void rejectStepBasedTimeKeys(const nlohmann::json& input, const std::string& fileName)
{
  if (!input.is_object()) return;

  // CVODE and ROS34PW2 choose their own steps.
  const nlohmann::json* integrator = findKeyCaseInsensitive(input, "BreakthroughIntegrator");
  const std::string integratorName = integrator != nullptr && integrator->is_string() ? integrator->get<std::string>() : "";
  const std::string adaptiveName = caseInSensStringCompare(integratorName, "CVODE")      ? "CVODE"
                                   : caseInSensStringCompare(integratorName, "ROS34PW2") ? "ROS34PW2"
                                                                                         : "";
  const bool adaptive = !adaptiveName.empty();
  const nlohmann::json* timeStepValue = findKeyCaseInsensitive(input, "TimeStep");
  const nlohmann::json* initSteps = findKeyCaseInsensitive(input, "NumberOfInitTimeSteps");
  const nlohmann::json* steps = findKeyCaseInsensitive(input, "NumberOfTimeSteps");
  const nlohmann::json* writeEvery = findKeyCaseInsensitive(input, "WriteEvery");
  const nlohmann::json* printEvery = findKeyCaseInsensitive(input, "PrintEvery");
  const nlohmann::json* phases = findKeyCaseInsensitive(input, "SwingAdsorptionPhases");
  const bool hasPhases = phases != nullptr && phases->is_array();
  const bool hasPhaseSteps =
      hasPhases && std::any_of(phases->begin(), phases->end(), [](const nlohmann::json& phase)
                               { return phase.is_object() && findKeyCaseInsensitive(phase, "NumberOfTimeSteps"); });

  const std::string legacyName = "_legacy_" + std::filesystem::path{fileName}.filename().string();
  const std::string convertHint = std::format(
      "Run 'python convert.py {}' to convert the file; the original is kept as {}.", fileName, legacyName);

  if (steps == nullptr && writeEvery == nullptr && printEvery == nullptr && !hasPhaseSteps)
  {
    if (adaptive && timeStepValue != nullptr)
    {
      throw std::runtime_error(std::format(
          "Error: TimeStep is not used by {0}, which chooses its own internal steps (capped by {0}MaximumTimeStep); "
          "remove it and set the output interval with WriteInterval.\n{1}",
          adaptiveName, convertHint));
    }
    if (adaptive && initSteps != nullptr)
    {
      throw std::runtime_error(
          std::format("Error: NumberOfInitTimeSteps is not used by {}; remove it.\n{}", adaptiveName, convertHint));
    }
    return;
  }

  const double timeStep =
      timeStepValue != nullptr && timeStepValue->is_number() ? timeStepValue->get<double>() : 0.0005;
  const auto seconds = [&](const nlohmann::json& count)
  { return count.is_number() ? std::format("{:.12g}", count.get<double>() * timeStep) : count.dump(); };
  const auto append = [](std::string& list, const std::string& item)
  { list += (list.empty() ? "" : ", ") + item; };

  std::string removed;
  if (timeStepValue != nullptr) append(removed, "\"TimeStep\": " + timeStepValue->dump());
  if (adaptive && initSteps != nullptr) append(removed, "\"NumberOfInitTimeSteps\": " + initSteps->dump());
  if (steps != nullptr) append(removed, "\"NumberOfTimeSteps\": " + steps->dump());
  if (writeEvery != nullptr) append(removed, "\"WriteEvery\": " + writeEvery->dump());
  if (printEvery != nullptr) append(removed, "\"PrintEvery\": " + printEvery->dump());

  std::string replacement;
  if (timeStepValue != nullptr && !adaptive) append(replacement, "\"TimeStep\": " + timeStepValue->dump());
  if (steps != nullptr && !hasPhases)
  {
    append(replacement, steps->is_string() ? std::string{"\"SimulationTime\": \"auto\""}
                                           : "\"SimulationTime\": " + seconds(*steps));
  }
  append(replacement, "\"WriteInterval\": " + (writeEvery != nullptr ? seconds(*writeEvery) : seconds(10000)));

  std::string message = "Error: step-based time keys were replaced by times in s.\n";
  message += removed.empty() ? "  Add:     " + replacement + "\n"
                             : "  Replace: " + removed + "\n  with:    " + replacement + "\n";
  if (hasPhases)
  {
    for (size_t i = 0; i < phases->size(); ++i)
    {
      const nlohmann::json& phase = (*phases)[i];
      const nlohmann::json* phaseSteps = phase.is_object() ? findKeyCaseInsensitive(phase, "NumberOfTimeSteps") : nullptr;
      if (phaseSteps == nullptr) continue;
      message += std::format("  SwingAdsorptionPhases[{}]: replace \"NumberOfTimeSteps\": {} with \"Duration\": {}\n",
                             i, phaseSteps->dump(), seconds(*phaseSteps));
    }
  }
  if (adaptive && timeStepValue != nullptr)
  {
    message += std::format("{0} chooses its own internal steps, so TimeStep is not used with {0}.\n", adaptiveName);
  }
  if (timeStepValue == nullptr)
  {
    message += "The step counts were converted with the former default TimeStep of 5e-4 s.\n";
  }
  if (steps != nullptr && hasPhases)
  {
    message += "The phase durations set the simulated time of a SwingAdsorption run, so the root NumberOfTimeSteps "
               "has no replacement.\n";
  }
  if (printEvery != nullptr)
  {
    message += "Progress is printed at every output, so PrintEvery has no replacement.\n";
  }
  throw std::runtime_error(message + convertHint);
}

/**
 * \brief Rejects the temperature keys that temperature laws replaced, and gives their replacement.
 *
 * The replacements follow the rules of convert.py (convert_temperature_keys), which rewrites a whole input file:
 * "nonIsothermal": true becomes a VantHoff law on the affinity of each site, whose Energy is the HeatOfAdsorption;
 * the five-parameter GAB [q_m, C_0, K_0, H_C, H_K] becomes [q_m, {VantHoff C_0, Energy H_C}, {VantHoff K_0,
 * Energy H_K}]; the activation energies of General chemisorption and of reactions become Arrhenius laws of their
 * rate coefficients; and the gibbsFreeEnergy dG of a reaction becomes HeatOfReaction -dG with a VantHoff
 * EquilibriumConstant.
 */
static void rejectRemovedTemperatureKeys(const nlohmann::json& input, const std::string& fileName)
{
  if (!input.is_object()) return;

  std::string found;
  const auto checkComponent = [&](const nlohmann::json& component, const std::string& where)
  {
    if (!component.is_object()) return;
    if (const nlohmann::json* nonIsothermal = findKeyCaseInsensitive(component, "nonIsothermal"))
    {
      found += nonIsothermal->is_boolean() && nonIsothermal->get<bool>()
                   ? "  " + where +
                         ": replace \"nonIsothermal\": true by writing the affinity of each site as {\"Value\": b, "
                         "\"Law\": \"VantHoff\"} (its Energy is the HeatOfAdsorption)\n"
                   : "  " + where + ": remove \"nonIsothermal\"\n";
    }
    for (const std::string siteKey : {"PhysisorptionSites", "ChemisorptionSites"})
    {
      const nlohmann::json* sites = findKeyCaseInsensitive(component, siteKey);
      if (sites == nullptr || !sites->is_array()) continue;
      for (size_t site = 0; site < sites->size(); ++site)
      {
        const nlohmann::json& entry = (*sites)[site];
        if (!entry.is_object()) continue;
        const std::string siteWhere = std::format("{}.{}[{}]", where, siteKey, site);
        const nlohmann::json* type = findKeyCaseInsensitive(entry, "Type");
        const nlohmann::json* parameters = findKeyCaseInsensitive(entry, "Parameters");
        if (type != nullptr && type->is_string() && caseInSensStringCompare(type->get<std::string>(), "GAB") &&
            parameters != nullptr && parameters->is_array() && parameters->size() == 5)
        {
          const nlohmann::json& values = *parameters;
          found += std::format(
              "  {}: replace the GAB Parameters {} with [{}, {{\"Value\": {}, \"Law\": \"VantHoff\", \"Energy\": {}}}, "
              "{{\"Value\": {}, \"Law\": \"VantHoff\", \"Energy\": {}}}]\n",
              siteWhere, values.dump(), values[0].dump(), values[1].dump(), values[3].dump(), values[2].dump(),
              values[4].dump());
        }
        const nlohmann::json* rateParameters = findKeyCaseInsensitive(entry, "RateEquationParameters");
        if (rateParameters == nullptr || !rateParameters->is_object()) continue;
        for (const std::string direction : {"Adsorption", "Desorption"})
        {
          if (const nlohmann::json* energy = findKeyCaseInsensitive(*rateParameters, direction + "ActivationEnergy"))
          {
            found += energy->is_number() && energy->get<double>() == 0.0
                         ? std::format("  {}: remove \"{}ActivationEnergy\": {} (the rate coefficient is constant)\n",
                                       siteWhere, direction, energy->dump())
                         : std::format("  {}: replace \"{}ActivationEnergy\": {} by \"{}RateCoefficient\": {{\"Value\": k, "
                                       "\"Law\": \"Arrhenius\", \"ActivationEnergy\": {}}}\n",
                                       siteWhere, direction, energy->dump(), direction, energy->dump());
          }
        }
      }
    }
  };

  if (const nlohmann::json* components = findKeyCaseInsensitive(input, "Components"))
  {
    size_t index = 0;
    for (auto it = components->begin(); it != components->end(); ++it, ++index)
    {
      checkComponent(it.value(), std::format("Components[{}]", index));
    }
  }
  if (const nlohmann::json* adsorbents = findKeyCaseInsensitive(input, "Adsorbents"); adsorbents && adsorbents->is_array())
  {
    for (size_t adsorbent = 0; adsorbent < adsorbents->size(); ++adsorbent)
    {
      const nlohmann::json& entry = (*adsorbents)[adsorbent];
      if (!entry.is_object()) continue;
      if (const nlohmann::json* components = findKeyCaseInsensitive(entry, "Components"); components && components->is_array())
      {
        for (size_t index = 0; index < components->size(); ++index)
        {
          checkComponent((*components)[index], std::format("Adsorbents[{}].Components[{}]", adsorbent, index));
        }
      }
      if (const nlohmann::json* overrides = findKeyCaseInsensitive(entry, "ComponentParameters"); overrides && overrides->is_object())
      {
        for (auto it = overrides->begin(); it != overrides->end(); ++it)
        {
          checkComponent(it.value(), std::format("Adsorbents[{}].ComponentParameters.{}", adsorbent, it.key()));
        }
      }
    }
  }
  if (const nlohmann::json* reactions = findKeyCaseInsensitive(input, "Reactions"); reactions && reactions->is_array())
  {
    for (size_t reaction = 0; reaction < reactions->size(); ++reaction)
    {
      const nlohmann::json& entry = (*reactions)[reaction];
      const nlohmann::json* kinetics = entry.is_object() ? findKeyCaseInsensitive(entry, "Kinetics") : nullptr;
      if (kinetics == nullptr || !kinetics->is_object()) continue;
      if (const nlohmann::json* energy = findKeyCaseInsensitive(*kinetics, "forwardActivationEnergy"))
      {
        found += energy->is_number() && energy->get<double>() == 0.0
                     ? std::format("  Reactions[{}].Kinetics: remove \"forwardActivationEnergy\": {} (the rate "
                                   "coefficient is constant)\n",
                                   reaction, energy->dump())
                     : std::format("  Reactions[{}].Kinetics: replace \"forwardActivationEnergy\": {} by "
                                   "\"ForwardRateCoefficient\": {{\"Value\": k, \"Law\": \"Arrhenius\", "
                                   "\"ActivationEnergy\": {}}}\n",
                                   reaction, energy->dump(), energy->dump());
      }
      if (const nlohmann::json* gibbs = findKeyCaseInsensitive(*kinetics, "gibbsFreeEnergy"))
      {
        const std::string heat = gibbs->is_number() ? std::format("{:.12g}", -gibbs->get<double>()) : "-dG";
        found += gibbs->is_number() && gibbs->get<double>() == 0.0
                     ? std::format("  Reactions[{}].Kinetics: remove \"gibbsFreeEnergy\": {} (K is constant and the "
                                   "reaction releases no heat)\n",
                                   reaction, gibbs->dump())
                     : std::format("  Reactions[{}].Kinetics: replace \"gibbsFreeEnergy\": {} by \"HeatOfReaction\": "
                                   "{} and write \"EquilibriumConstant\": {{\"Value\": K, \"Law\": \"VantHoff\"}}\n",
                                   reaction, gibbs->dump(), heat);
      }
    }
  }
  if (found.empty()) return;

  const std::string legacyName = "_legacy_" + std::filesystem::path{fileName}.filename().string();
  throw std::runtime_error(
      "Error: temperature dependence is now written as temperature laws on the parameters.\n" + found +
      std::format("Run 'python convert.py {}' to convert the file; the original is kept as {} (or "
                  "_legacy_2_..., if that exists).",
                  fileName, legacyName));
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
  rejectStepBasedTimeKeys(parsed_data, fileName);
  rejectRemovedTemperatureKeys(parsed_data, fileName);
  requireOnlyKnownKeys(parsed_data,
                       {"SimulationType",
                        "MixturePredictionMethod",
                        "MPDSettings",
                        "IASTMethod",
                        "BreakthroughIntegrator",
                        "ConvectionScheme",
                        "InletCondition",
                        "FlowDirection",
                        "AxialDispersionReferenceTemperature",
                        "AxialDispersionReferencePressure",
                        "BoundaryCondition",
                        "FluidPhase",
                        "LiquidDensity",
                        "pHMode",
                        "pHValue",
                        "pKw",
                        "pHComponent",
                        "RelativeHumidity",
                        "WaterSaturationPressureModel",
                        "WaterSaturationPressure",
                        "PressureScale",
                        "ReadColumnFile",
                        "DisplayName",
                        "DebugForceMultibed",  // ignored: every column runs on the same solver
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
                        "CVODERelativeTolerance",
                        "CVODEAbsoluteToleranceConcentration",
                        "CVODEAbsoluteToleranceLoading",
                        "CVODEAbsoluteToleranceTemperature",
                        "CVODEMaximumTimeStep",
                        "CVODEMaximumOrder",
                        "CVODELinearSolver",
                        "CVODEKrylovDimension",
                        "CVODEPreconditioner",
                        "ROS34PW2RelativeTolerance",
                        "ROS34PW2AbsoluteToleranceConcentration",
                        "ROS34PW2AbsoluteToleranceLoading",
                        "ROS34PW2AbsoluteToleranceTemperature",
                        "ROS34PW2MaximumTimeStep",
                        "WriteInterval",
                        "ColumnLength",
                        "ColumnDistances",
                        "NumberOfGridPoints",
                        "ColumnPressure",
                        "ColumnLoading",
                        "ColumnError",
                        "FittingWeighting",
                        "SimulationTime",
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
                        "NumberOfCycles",
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
                           {{"RungeKutta3", 0}, {"CVODE", 1}, {"ROS34PW2", 2}});

  readOptionalMappedString(parsed_data, "ConvectionScheme", convectionScheme, {{"Upwind", 0}, {"WENO3", 1}});
  readOptionalMappedString(parsed_data, "InletCondition", inletCondition, {{"Dirichlet", 0}, {"Danckwerts", 1}});
  readOptionalMappedString(parsed_data, "FlowDirection", flowDirection, {{"Forward", 0}, {"Reverse", 1}});

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

  readOptionalNumber<double>(parsed_data, "Temperature", temperature);
  readOptionalNumber<double>(parsed_data, "ColumnVoidFraction", columnVoidFraction);
  readOptionalNumber<double>(parsed_data, "DynamicViscosity", dynamicViscosity);
  readOptionalNumber<double>(parsed_data, "AxialDispersionReferenceTemperature", axialDispersionReferenceTemperature);
  readOptionalNumber<double>(parsed_data, "AxialDispersionReferencePressure", axialDispersionReferencePressure);
  readOptionalNumber<double>(parsed_data, "LiquidDensity", liquidDensity);
  readOptionalNumber<double>(parsed_data, "pHValue", pHValue);
  readOptionalNumber<double>(parsed_data, "pKw", pKw);
  readOptionalNumber<double>(parsed_data, "RelativeHumidity", relativeHumidity);
  // A WaterSaturationPressure without a model is the Constant model.
  const bool hasConstantSaturationPressure = containsKeyCaseInsensitive(parsed_data, "WaterSaturationPressure");
  size_t saturationPressureModel = hasConstantSaturationPressure ? 2 : 0;
  readOptionalMappedString(parsed_data, "WaterSaturationPressureModel", saturationPressureModel,
                           {{"Antoine", 0}, {"Wagner", 1}, {"Constant", 2}});
  waterSaturationPressure.model = WaterSaturationPressure::Model(saturationPressureModel);
  readOptionalNumber<double>(parsed_data, "WaterSaturationPressure", waterSaturationPressure.constantPressure);
  if (waterSaturationPressure.model == WaterSaturationPressure::Model::Constant &&
      !(std::isfinite(waterSaturationPressure.constantPressure) && waterSaturationPressure.constantPressure > 0.0))
  {
    throw std::runtime_error(
        "Error: the Constant WaterSaturationPressureModel needs a positive WaterSaturationPressure in Pa");
  }
  if (waterSaturationPressure.model != WaterSaturationPressure::Model::Constant && hasConstantSaturationPressure)
  {
    throw std::runtime_error("Error: WaterSaturationPressure is the pressure of the Constant "
                             "WaterSaturationPressureModel");
  }
  if (relativeHumidity.has_value() && !(std::isfinite(*relativeHumidity) && *relativeHumidity >= 0.0))
  {
    throw std::runtime_error("Error: RelativeHumidity must be a finite non-negative number");
  }
  readOptionalNumber<double>(parsed_data, "ParticleDiameter", particleDiameter);
  readOptionalNumber<double>(parsed_data, "ParticleDensity", particleDensity);

  readOptionalNumber<double>(parsed_data, "InletPressure", inletPressure);

  readOptionalNumber<double>(parsed_data, "OutletPressure", outletPressure);
  readOptionalNumber<double>(parsed_data, "PressureStart", pressureStart);
  readOptionalNumber<double>(parsed_data, "PressureEnd", pressureEnd);
  readOptionalNonNegativeInteger(parsed_data, "NumberOfPressurePoints", numberOfPressurePoints);
  readOptionalNumber<double>(parsed_data, "PressureGradient", pressureGradient);
  readOptionalNumber<double>(parsed_data, "ColumnEntranceVelocity", columnEntranceVelocity);
  readOptionalNonNegativeInteger(parsed_data, "NumberOfInitTimeSteps", numberOfInitTimeSteps);
  readOptionalNumber<double>(parsed_data, "TimeStep", timeStep);
  readOptionalNumber<double>(parsed_data, "CVODERelativeTolerance", cvodeRelativeTolerance);
  readOptionalNumber<double>(parsed_data, "CVODEAbsoluteToleranceConcentration", cvodeAbsoluteToleranceConcentration);
  readOptionalNumber<double>(parsed_data, "CVODEAbsoluteToleranceLoading", cvodeAbsoluteToleranceLoading);
  readOptionalNumber<double>(parsed_data, "CVODEAbsoluteToleranceTemperature", cvodeAbsoluteToleranceTemperature);
  readOptionalNumber<double>(parsed_data, "CVODEMaximumTimeStep", cvodeMaximumTimeStep);
  readOptionalNonNegativeInteger(parsed_data, "CVODEMaximumOrder", cvodeMaximumOrder);
  if (cvodeMaximumOrder < 1 || cvodeMaximumOrder > 5)
  {
    throw std::runtime_error("Error: CVODEMaximumOrder must be between 1 and 5 (BDF orders)");
  }
  readOptionalMappedString(parsed_data, "CVODELinearSolver", cvodeLinearSolver,
                           {{"Dense", 0}, {"SPGMR", 1}, {"KLU", 2}});
  readOptionalNonNegativeInteger(parsed_data, "CVODEKrylovDimension", cvodeKrylovDimension);
  readOptionalMappedString(parsed_data, "CVODEPreconditioner", cvodePreconditioner,
                           {{"None", 0}, {"SorptionJacobian", 1}});
  readOptionalNumber<double>(parsed_data, "ROS34PW2RelativeTolerance", ros34pw2RelativeTolerance);
  readOptionalNumber<double>(parsed_data, "ROS34PW2AbsoluteToleranceConcentration",
                             ros34pw2AbsoluteToleranceConcentration);
  readOptionalNumber<double>(parsed_data, "ROS34PW2AbsoluteToleranceLoading", ros34pw2AbsoluteToleranceLoading);
  readOptionalNumber<double>(parsed_data, "ROS34PW2AbsoluteToleranceTemperature",
                             ros34pw2AbsoluteToleranceTemperature);
  readOptionalNumber<double>(parsed_data, "ROS34PW2MaximumTimeStep", ros34pw2MaximumTimeStep);
  if (!(ros34pw2RelativeTolerance > 0.0) || !(ros34pw2AbsoluteToleranceConcentration >= 0.0) ||
      !(ros34pw2AbsoluteToleranceLoading >= 0.0) || !(ros34pw2AbsoluteToleranceTemperature >= 0.0) ||
      !(ros34pw2MaximumTimeStep >= 0.0))
  {
    throw std::runtime_error(
        "Error: ROS34PW2RelativeTolerance must be positive, and the ROS34PW2 absolute tolerances and "
        "ROS34PW2MaximumTimeStep (0 = unlimited) non-negative");
  }
  if (cvodePreconditioner != 0 && cvodeLinearSolver != 1)
  {
    throw std::runtime_error("Error: CVODEPreconditioner needs CVODELinearSolver SPGMR");
  }
  readOptionalNumber<double>(parsed_data, "WriteInterval", writeInterval);
  readOptionalNumber<double>(parsed_data, "ColumnLength", columnLength);
  readOptionalNonNegativeInteger(parsed_data, "NumberOfGridPoints", numberOfGridPoints);
  if (containsKeyCaseInsensitive(parsed_data, "ColumnDistances"))
  {
    columnDistances = getNumberListOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ColumnDistances", ""),
                                                   "ColumnDistances", "");
    if (columnDistances.size() < 2)
    {
      throw std::runtime_error("Error: ColumnDistances must contain at least two node positions");
    }
  }
  readOptionalNonNegativeInteger(parsed_data, "ColumnPressure", columnPressure);
  readOptionalNonNegativeInteger(parsed_data, "ColumnLoading", columnLoading);
  readOptionalNonNegativeInteger(parsed_data, "ColumnError", columnError);
  if (columnPressure == 0 || columnLoading == 0 || columnError == 0)
  {
    throw std::runtime_error("Error: ColumnPressure, ColumnLoading, and ColumnError are one-based column numbers");
  }
  readOptionalMappedString(parsed_data, "FittingWeighting", fittingWeighting, {{"Absolute", 0}, {"Relative", 1}});

  const bool hasSimulationTime = containsKeyCaseInsensitive(parsed_data, "SimulationTime");
  if (hasSimulationTime)
  {
    const nlohmann::json& v = requireKeyCaseInsensitive(parsed_data, "SimulationTime", "");
    if (v.is_string())
    {
      if (!caseInSensStringCompare(v.get<std::string>(), "auto"))
      {
        throw std::runtime_error("Error: SimulationTime must be a time in s or 'auto'");
      }
      autoSimulationTime = true;
    }
    else
    {
      simulationTime = getNumberOrThrow<double>(v, "SimulationTime", "");
      autoSimulationTime = false;
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
         "ForchheimerCoefficient", "AxialDispersivity", "TaylorArisDispersion"},
        "Geometry");

    const std::string geometryType =
        getStringOrThrow(requireKeyCaseInsensitive(*geometryJson, "Type", "Geometry"), "Type", "Geometry");
    const double axialDispersivity = geometryDoubleOrDefault(*geometryJson, {"AxialDispersivity"}, 0.0);
    bool taylorArisDispersion = false;
    readOptionalBool(*geometryJson, "TaylorArisDispersion", taylorArisDispersion);

    if (caseInSensStringCompare(geometryType, "PackedBed"))
    {
      if (taylorArisDispersion)
      {
        throw std::runtime_error("Error: Geometry key 'TaylorArisDispersion' applies only to a Monolith");
      }
      const double eps = geometryDoubleOrDefault(*geometryJson, {"ColumnVoidFraction"}, columnVoidFraction);
      const double dp = geometryDoubleOrDefault(*geometryJson, {"ParticleDiameter"}, particleDiameter);
      const double di = geometryDoubleOrDefault(*geometryJson, {"InternalDiameter"}, internalDiameter);
      const double douter = geometryDoubleOrDefault(*geometryJson, {"OuterDiameter"}, outerDiameter);
      geometry = makeGeometry(PackedBedTubeSpec{.voidFraction = eps,
                                                .particleDiameter = dp,
                                                .internalDiameter = di,
                                                .outerDiameter = douter,
                                                .axialDispersivity = axialDispersivity});
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
                                           .forchheimerCoefficient = forchheimerCoefficient,
                                           .axialDispersivity = axialDispersivity,
                                           .taylorArisDispersion = taylorArisDispersion});
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
  // Per-phase feed and product-end compositions by component name; resolved to component indices once Components are
  // parsed.
  std::vector<std::vector<std::pair<std::string, double>>> phaseFeedByName;
  std::vector<std::vector<std::pair<std::string, double>>> phaseProductEndByName;
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

      requireOnlyKnownKeys(value,
                           {"Name", "Temperature", "InletPressure", "Duration", "GasPhaseMolFractions",
                            "ColumnEntranceVelocity", "FlowDirection", "ProductEndGasPhaseMolFractions", "ClosedEnd",
                            "PressureRamp"},
                           context);

      SwingAdsorptionPhase phase;
      phase.name = "Phase " + std::to_string(i + 1);
      readOptionalString(value, "Name", phase.name);
      readOptionalNumber<double>(value, "Temperature", phase.temperature);
      readOptionalNumber<double>(value, "InletPressure", phase.inletPressure);
      readOptionalNumber<double>(value, "ColumnEntranceVelocity", phase.columnEntranceVelocity);
      if (findKeyCaseInsensitive(value, "FlowDirection") != nullptr)
      {
        size_t direction = 0;
        readOptionalMappedString(value, "FlowDirection", direction, {{"Forward", 0}, {"Reverse", 1}});
        phase.flowDirection = direction;
      }
      if (findKeyCaseInsensitive(value, "ClosedEnd") != nullptr)
      {
        // Gas enters at the open end, so that ClosedEnd also sets the direction.
        if (phase.flowDirection.has_value())
        {
          throw std::runtime_error("Error: ClosedEnd sets the open end; remove FlowDirection (" + context + ")");
        }
        size_t direction = 0;
        readOptionalMappedString(value, "ClosedEnd", direction, {{"ProductEnd", 0}, {"FeedEnd", 1}});
        phase.flowDirection = direction;
        phase.closedEnd = true;
      }
      if (const nlohmann::json* ramp = findKeyCaseInsensitive(value, "PressureRamp"))
      {
        if (!ramp->is_object())
        {
          throw std::runtime_error("Error: PressureRamp must be an object with Shape and TimeConstant (" + context + ")");
        }
        requireOnlyKnownKeys(*ramp, {"Shape", "TimeConstant"}, context);
        requireKeyCaseInsensitive(*ramp, "Shape", context);
        readOptionalMappedString(*ramp, "Shape", phase.pressureRampShape,
                                 {{"Exponential", size_t{1}}, {"Linear", size_t{2}}});
        phase.pressureRampTimeConstant =
            getNumberOrThrow<double>(requireKeyCaseInsensitive(*ramp, "TimeConstant", context), "TimeConstant", context);
        if (!std::isfinite(phase.pressureRampTimeConstant) || phase.pressureRampTimeConstant <= 0.0)
        {
          throw std::runtime_error("Error: PressureRamp TimeConstant must be a positive time in s (" + context + ")");
        }
      }

      for (const auto& [key, byName] : {std::pair{"GasPhaseMolFractions", &phaseFeedByName},
                                        std::pair{"ProductEndGasPhaseMolFractions", &phaseProductEndByName}})
      {
        std::vector<std::pair<std::string, double>> fractionsByName;
        if (const nlohmann::json* fractions = findKeyCaseInsensitive(value, key))
        {
          if (!fractions->is_object() || fractions->empty())
          {
            throw std::runtime_error(std::string("Error: ") + key +
                                     " must be a non-empty object of component name to mole fraction (" + context +
                                     ")");
          }
          for (auto it = fractions->begin(); it != fractions->end(); ++it)
          {
            fractionsByName.emplace_back(it.key(), getNumberOrThrow<double>(it.value(), it.key(), context));
          }
        }
        byName->push_back(std::move(fractionsByName));
      }

      const nlohmann::json* duration = findKeyCaseInsensitive(value, "Duration");
      if (duration == nullptr)
      {
        throw std::runtime_error("Error: required key 'Duration' missing (" + context + ")");
      }
      phase.duration = getNumberOrThrow<double>(*duration, "Duration", context);
      if (!std::isfinite(phase.duration) || phase.duration <= 0.0)
      {
        throw std::runtime_error("Error: SwingAdsorption phase Duration must be a positive time in s (" + context +
                                 ")");
      }
      swingAdsorptionPhases.push_back(phase);
    }
  }

  // NumberOfCycles repeats the whole SwingAdsorptionPhases sequence; each repetition is one cycle.
  if (const nlohmann::json* cycles = findKeyCaseInsensitive(parsed_data, "NumberOfCycles"))
  {
    if (!cycles->is_number_integer() || cycles->get<std::int64_t>() < 1)
    {
      throw std::runtime_error("Error: NumberOfCycles must be a positive integer");
    }
    if (swingAdsorptionPhases.empty())
    {
      throw std::runtime_error("Error: NumberOfCycles requires SwingAdsorptionPhases");
    }
    const size_t numberOfCycles = static_cast<size_t>(cycles->get<std::int64_t>());
    if (numberOfCycles > 1)
    {
      const std::vector<SwingAdsorptionPhase> cyclePhases = swingAdsorptionPhases;
      const std::vector<std::vector<std::pair<std::string, double>>> cycleFeeds = phaseFeedByName;
      const std::vector<std::vector<std::pair<std::string, double>>> cycleProductEnds = phaseProductEndByName;
      swingAdsorptionPhases.clear();
      phaseFeedByName.clear();
      phaseProductEndByName.clear();
      for (size_t cycle = 1; cycle <= numberOfCycles; ++cycle)
      {
        for (size_t i = 0; i < cyclePhases.size(); ++i)
        {
          SwingAdsorptionPhase phase = cyclePhases[i];
          phase.name += " (cycle " + std::to_string(cycle) + ")";
          swingAdsorptionPhases.push_back(std::move(phase));
          phaseFeedByName.push_back(cycleFeeds[i]);
          phaseProductEndByName.push_back(cycleProductEnds[i]);
        }
      }
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

    if (hasSimulationTime)
    {
      throw std::runtime_error(
          "Error: the SwingAdsorptionPhases durations set the simulated time; remove SimulationTime");
    }
    simulationTime =
        std::accumulate(swingAdsorptionPhases.begin(), swingAdsorptionPhases.end(), 0.0,
                        [](double total, const SwingAdsorptionPhase& phase) { return total + phase.duration; });
    autoSimulationTime = false;
  }

  // Components
  // The heat and reference temperature of each component as given, from which the adsorbent overrides start, and the
  // components whose physisorption heat is ambiguous (checked once the energy balance is known).
  std::vector<ComponentTemperatureInput> componentTemperatureInputs;
  std::vector<std::string> ambiguousHeats;
  if (containsKeyCaseInsensitive(parsed_data, "Components"))
  {
    const nlohmann::json& comps = requireKeyCaseInsensitive(parsed_data, "Components", "");

    components.clear();

    // The order of the list sets the component indices (and the component_<i>_<name>.data files). A JSON object has
    // no order once parsed (nlohmann::json sorts its keys), so only the list form is accepted (B-IO-10).
    if (!comps.is_array())
    {
      throw std::runtime_error("Error: 'Components' must be a list of component objects, [{\"Name\": ...}, ...]");
    }
    components.reserve(comps.size());
    for (std::size_t componentId = 0; componentId < comps.size(); ++componentId)
    {
      components.push_back(
          parseComponentObject(componentId, comps[componentId], componentTemperatureInputs.emplace_back()));
    }
    for (size_t comp = 0; comp < components.size(); ++comp)
    {
      if (!componentTemperatureInputs[comp].ambiguousEnergies.empty())
      {
        ambiguousHeats.push_back(components[comp].name + ": " + componentTemperatureInputs[comp].ambiguousEnergies);
      }
    }
  }

  if (simulationType == SimulationType::Fitting)
  {
    for (const Component& component : components)
    {
      if (std::any_of(component.isotherm.sites.begin(), component.isotherm.sites.end(),
                      [](const Isotherm& site) { return !site.temperatureLaws.empty(); }))
      {
        throw std::runtime_error("Error: Fitting fits isotherm parameters at one temperature; give plain numbers "
                                 "instead of temperature laws (Component " +
                                 component.name + ")");
      }
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
    if (!components.empty() && !(sum > 0.0))
    {
      throw std::runtime_error("Error: GasPhaseMolFraction must not be zero for all components");
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

  // Resolve per-phase feed and product-end compositions; components a phase does not list are absent from them.
  for (size_t phaseIndex = 0; phaseIndex < phaseFeedByName.size(); ++phaseIndex)
  {
    SwingAdsorptionPhase& phase = swingAdsorptionPhases[phaseIndex];
    const std::string context = "SwingAdsorptionPhases '" + phase.name + "'";
    if ((phase.closedEnd || phase.pressureRampShape != 0) && (fluidPhase != 0 || boundaryCondition != 0))
    {
      throw std::runtime_error(
          "Error: ClosedEnd and PressureRamp need a gas and BoundaryCondition InletPressureInletVelocity (" + context +
          ")");
    }
    if (phase.closedEnd && inletCondition != 1)
    {
      throw std::runtime_error("Error: ClosedEnd needs InletCondition Danckwerts (" + context + ")");
    }
    if (phase.closedEnd && (phase.columnEntranceVelocity.has_value() || !phaseProductEndByName[phaseIndex].empty()))
    {
      throw std::runtime_error(
          "Error: a closed-end phase takes no ColumnEntranceVelocity or ProductEndGasPhaseMolFractions (" + context +
          ")");
    }
    if (!phaseProductEndByName[phaseIndex].empty() && phase.flowDirection.value_or(flowDirection) == 1)
    {
      throw std::runtime_error(
          "Error: ProductEndGasPhaseMolFractions sets the gas that flows back in at z = L in a forward phase; the "
          "feed of a reverse phase is GasPhaseMolFractions (" + context + ")");
    }
    for (const auto& [key, byName, target] :
         {std::tuple{"GasPhaseMolFractions", &phaseFeedByName[phaseIndex], &phase.gasPhaseMolFractions},
          std::tuple{"ProductEndGasPhaseMolFractions", &phaseProductEndByName[phaseIndex],
                     &phase.productEndGasPhaseMolFractions}})
    {
      if (byName->empty()) continue;
      if (fluidPhase != 0)
      {
        throw std::runtime_error(std::string("Error: ") + key + " is only supported for gas-phase columns (" +
                                 context + ")");
      }

      std::vector<double> fractions(components.size(), 0.0);
      for (const auto& [name, fraction] : *byName)
      {
        const auto component = std::find_if(components.begin(), components.end(),
                                            [&](const Component& candidate) { return candidate.name == name; });
        if (component == components.end())
        {
          throw std::runtime_error(std::string("Error: ") + key + " names unknown component '" + name + "' (" +
                                   context + ")");
        }
        if (!std::isfinite(fraction) || fraction < 0.0)
        {
          throw std::runtime_error(std::string("Error: ") + key + " must be finite and non-negative (" + context +
                                   ")");
        }
        fractions[static_cast<size_t>(component - components.begin())] = fraction;
      }
      const double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
      if (!(sum > 0.0))
      {
        throw std::runtime_error(std::string("Error: ") + key + " must not all be zero (" + context + ")");
      }
      if (std::abs(sum - 1.0) > 1e-15)
      {
        std::print("Normalizing: {} {} did not sum exactly to unity!\n\n", context, key);
        for (double& fraction : fractions) fraction /= sum;
      }
      *target = std::move(fractions);
    }
  }

  numberOfCarrierGases = 0;
  carrierGasComponent = 0;
  for (size_t j = 0; j < components.size(); ++j)
  {
    if (fluidPhase == 0 && components[j].isCarrierGas)
    {
      carrierGasComponent = j;
      if (components[j].isotherm.sites.empty())
      {
        std::vector<double> values{1.0, 0.0};
        Isotherm isotherm = Isotherm(Isotherm::Type::Langmuir, values);
        components[j].isotherm.add(isotherm);
      }

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
      readOptionalNonNegativeInteger(parsed_data, "pHComponent", pHComponent);
      if (pHComponent >= components.size())
      {
        throw std::runtime_error("Error: pHComponent must be a zero-based component identifier less than the number "
                                 "of Components");
      }
    }
  }

  reactions.clear();
  // Adsorbent names a reaction is restricted to; resolved to indices once the Adsorbents are parsed.
  std::vector<std::vector<std::string>> reactionAdsorbentNames;
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
                            "ForwardOrders", "BackwardOrders", "RateLimitTime", "Adsorbents"},
                           context);

      std::vector<std::string> adsorbentNamesOfReaction;
      if (const nlohmann::json* adsorbentsOfReaction = findKeyCaseInsensitive(value, "Adsorbents"))
      {
        if (!adsorbentsOfReaction->is_array() || adsorbentsOfReaction->empty())
        {
          throw std::runtime_error("Error: Reaction Adsorbents must be a non-empty array of adsorbent names (" +
                                   context + ")");
        }
        for (const nlohmann::json& name : *adsorbentsOfReaction)
        {
          adsorbentNamesOfReaction.push_back(getStringOrThrow(name, "Adsorbents", context));
        }
      }
      reactionAdsorbentNames.push_back(std::move(adsorbentNamesOfReaction));

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
      else if (caseInSensStringCompare(phase, "BulkConcentration"))
      {
        reaction.phase = Reaction::Phase::BulkConcentration;
      }
      else
      {
        throw std::runtime_error(
            "Error: Reaction Phase must be Physisorbed, Chemisorbed, PoreConcentration, or BulkConcentration (" +
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
          reaction.phase != Reaction::Phase::PoreConcentration && reaction.phase != Reaction::Phase::BulkConcentration)
      {
        throw std::runtime_error("Error: LH/LHHW reactions require Phase=PoreConcentration or BulkConcentration (" +
                                 context + ")");
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
      const std::string kineticsContext = context + " Kinetics";
      requireOnlyKnownKeys(kinetics,
                           {"ForwardRateCoefficient", "EquilibriumConstant", "HeatOfReaction", "ReferenceTemperature",
                            "AdsorptionConstants"},
                           kineticsContext);

      std::optional<double> referenceTemperature;
      readOptionalNumber<double>(kinetics, "ReferenceTemperature", referenceTemperature);
      if (referenceTemperature.has_value() && !(std::isfinite(*referenceTemperature) && *referenceTemperature > 0.0))
      {
        throw std::runtime_error("Error: ReferenceTemperature must be positive (" + kineticsContext + ")");
      }
      reaction.forwardRateCoefficient =
          readTemperatureLaw(requireKeyCaseInsensitive(kinetics, "ForwardRateCoefficient", kineticsContext),
                             {"ForwardRateCoefficient", TemperatureLaw::mask(TemperatureLaw::Type::Arrhenius), false},
                             referenceTemperature, kineticsContext);
      reaction.equilibriumConstant =
          readTemperatureLaw(requireKeyCaseInsensitive(kinetics, "EquilibriumConstant", kineticsContext),
                             {"EquilibriumConstant", TemperatureLaw::mask(TemperatureLaw::Type::VantHoff), true},
                             referenceTemperature, kineticsContext);

      // K(T) and the reaction heat are linked like an affinity and the heat of adsorption: a VantHoff K without
      // Energy takes HeatOfReaction, and without HeatOfReaction the heat is the energy of K (explicit heat wins).
      std::optional<double> heatOfReaction;
      readOptionalNumber<double>(kinetics, "HeatOfReaction", heatOfReaction);
      TemperatureLaw& equilibriumConstant = reaction.equilibriumConstant;
      if (equilibriumConstant.energyInherited)
      {
        if (!heatOfReaction.has_value())
        {
          throw std::runtime_error("Error: the VantHoff law of the EquilibriumConstant has no \"Energy\", and there "
                                   "is no HeatOfReaction to take it from (" +
                                   kineticsContext + ")");
        }
        equilibriumConstant.coefficient = *heatOfReaction;
      }
      reaction.heatOfReaction = heatOfReaction.value_or(
          equilibriumConstant.type == TemperatureLaw::Type::VantHoff ? equilibriumConstant.coefficient : 0.0);

      if (reaction.forwardRateCoefficient.value < 0.0)
      {
        throw std::runtime_error("Error: ForwardRateCoefficient must be non-negative (" + context + ")");
      }
      if (reaction.equilibriumConstant.value <= 0.0)
      {
        throw std::runtime_error("Error: EquilibriumConstant must be positive (" + context + ")");
      }
      if (!std::isfinite(reaction.heatOfReaction))
      {
        throw std::runtime_error("Error: HeatOfReaction must be finite (" + context + ")");
      }

      // LH/LHHW coverages take one adsorption constant per component, given by name; unnamed components do not adsorb.
      const nlohmann::json* adsorptionConstants = findKeyCaseInsensitive(kinetics, "AdsorptionConstants");
      if (reaction.style == Reaction::Style::GeneralPowerLaw)
      {
        if (adsorptionConstants != nullptr)
        {
          throw std::runtime_error("Error: AdsorptionConstants apply only to LH/LHHW reactions (" + kineticsContext +
                                   ")");
        }
      }
      else
      {
        if (adsorptionConstants == nullptr || !adsorptionConstants->is_object() || adsorptionConstants->empty())
        {
          throw std::runtime_error("Error: LH/LHHW reactions need AdsorptionConstants, an object of component names "
                                   "and adsorption constants K_j (" +
                                   kineticsContext + ")");
        }
        reaction.adsorptionConstants.assign(components.size(), TemperatureLaw(0.0));
        std::set<size_t> named;
        for (auto it = adsorptionConstants->begin(); it != adsorptionConstants->end(); ++it)
        {
          const std::string constantContext = kineticsContext + " AdsorptionConstants." + it.key();
          const size_t comp =
              readReactionParticipants(nlohmann::json::array({it.key()}), components, "AdsorptionConstants", context)
                  .front();
          if (!named.insert(comp).second)
          {
            throw std::runtime_error("Error: component '" + components[comp].name +
                                     "' appears twice in AdsorptionConstants (" + kineticsContext + ")");
          }
          TemperatureLaw& constant = reaction.adsorptionConstants[comp];
          constant = readTemperatureLaw(
              it.value(), {"AdsorptionConstant", TemperatureLaw::mask(TemperatureLaw::Type::VantHoff), false},
              referenceTemperature, constantContext);
          if (constant.value < 0.0)
          {
            throw std::runtime_error("Error: AdsorptionConstants must be non-negative (" + constantContext + ")");
          }
        }
        // A reactant without a constant has zero coverage, which would make the forward rate zero everywhere.
        for (size_t comp : reaction.reactants)
        {
          if (reaction.adsorptionConstants[comp].value <= 0.0)
          {
            throw std::runtime_error("Error: reactant '" + components[comp].name +
                                     "' needs a positive constant in AdsorptionConstants, or its coverage and the "
                                     "forward rate are zero (" +
                                     kineticsContext + ")");
          }
        }
      }
      reactions.push_back(std::move(reaction));
    }
  }

  adsorbentComponents.clear();
  adsorbentLengths.clear();
  adsorbentInterfaceLengths.clear();
  adsorbentMixFractions.clear();
  adsorbentGridPoints.clear();
  adsorbentInterfaceGridPoints.clear();
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
                            "MixFraction", "ComponentParameters", "Components", "AxialDispersivity"},
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
      double adsorbentDispersivity = geometry.dispersion.dispersivity;
      readOptionalNumber<double>(adsorbent, "ParticleDiameter", adsorbentParticleDiameter);
      readOptionalNumber<double>(adsorbent, "AxialDispersivity", adsorbentDispersivity);
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
                                               .outerDiameter = geometry.dimensions.outerDiameter,
                                               .axialDispersivity = adsorbentDispersivity})
              : geometry);
      if (adsorbentDispersivity < 0.0)
      {
        throw std::runtime_error("Error: AxialDispersivity must be non-negative (" + context + ")");
      }
      adsorbentGeometries.back().dispersion.dispersivity = adsorbentDispersivity;

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
                                "ReferenceTemperature", "ChemisorptionSites", "PhysisorptionSites"},
                               componentContext);

          readOptionalNumber<double>(params, "MassTransferCoefficient", comp.massTransferCoefficient);
          readOptionalNumber<double>(params, "AxialDispersionCoefficient", comp.axialDispersionCoefficient);
          // The override starts from the heat and reference temperature given for the component.
          ComponentTemperatureInput temperatureInput = componentTemperatureInputs[componentIndex->second];
          readOptionalNumber<double>(params, "HeatOfAdsorption", temperatureInput.heatOfAdsorption);
          readOptionalNumber<double>(params, "ReferenceTemperature", temperatureInput.referenceTemperature);
          readChemisorption(comp.chemisorption, params, comp.massTransferCoefficient,
                            temperatureInput.heatOfAdsorption, temperatureInput.referenceTemperature,
                            componentContext);

          const bool hasIsothermOverride = containsKeyCaseInsensitive(params, "PhysisorptionSites");

          if (hasIsothermOverride)
          {
            comp.isotherm = MultiSiteIsotherm{};
          }
          readPhysisorptionSites(comp, params, temperatureInput.referenceTemperature, temperatureInput.heatOfAdsorption,
                                 componentContext);

          // Affinity laws that inherited the component heat follow an overridden heat.
          const std::vector<double> energies =
              linkHeat(comp.isotherm.sites, temperatureInput.heatOfAdsorption, componentContext);
          comp.heatOfAdsorption =
              temperatureInput.heatOfAdsorption.value_or(energies.size() == 1 ? energies.front() : 0.0);
          if (!temperatureInput.heatOfAdsorption.has_value() && energies.size() > 1)
          {
            ambiguousHeats.push_back(componentContext);
          }
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
          ComponentTemperatureInput temperatureInput;
          Component adsorbentComponent =
              parseComponentObject(componentId, adsorbentComponentsJson[componentId], temperatureInput);
          if (!temperatureInput.ambiguousEnergies.empty())
          {
            ambiguousHeats.push_back(context + ", " + adsorbentComponent.name + ": " +
                                     temperatureInput.ambiguousEnergies);
          }
          auto componentIndex = componentIndexByName.find(adsorbentComponent.name);
          if (componentIndex == componentIndexByName.end())
          {
            throw std::runtime_error("Error: unknown component '" + adsorbentComponent.name +
                                     "' in Adsorbent Components (" + context + ")");
          }
          adsorbentComponent.id = componentsForAdsorbent[componentIndex->second].id;
          // Water is a property of the component, set in the top-level Components.
          if (containsKeyCaseInsensitive(adsorbentComponentsJson[componentId], "Water") &&
              adsorbentComponent.isWater != componentsForAdsorbent[componentIndex->second].isWater)
          {
            throw std::runtime_error("Error: Water of component '" + adsorbentComponent.name +
                                     "' differs from the top-level Components (" + context + ")");
          }
          adsorbentComponent.isWater = componentsForAdsorbent[componentIndex->second].isWater;
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
      std::vector<size_t> sectionInterfaceGridPoints;
      bool hasSectionGridPoints = false;
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

        size_t sectionIntervals = 0;
        if (containsKeyCaseInsensitive(section, "NumberOfGridPoints"))
        {
          hasSectionGridPoints = true;
          const nlohmann::json& count = requireKeyCaseInsensitive(section, "NumberOfGridPoints", context);
          if (!count.is_number_integer() || count < 0)
          {
            throw std::runtime_error("Error: ColumnSections NumberOfGridPoints must be a non-negative integer (" +
                                     context + ")");
          }
          sectionIntervals = getNumberOrThrow<size_t>(count, "NumberOfGridPoints", context);
        }

        if (hasInterface)
        {
          if (nextAdsorbentInOrder == 0 || nextAdsorbentInOrder == adsorbentsJson.size() ||
              sectionInterfaceLengths.size() + 1 != nextAdsorbentInOrder)
          {
            throw std::runtime_error("Error: each interface must occur between neighboring adsorbent sections (" +
                                     context + ")");
          }
          const double interfaceLength = getNumberOrThrow<double>(
              requireKeyCaseInsensitive(section, "InterfaceLength", context), "InterfaceLength", context);
          if (interfaceLength < 0.0)
          {
            throw std::runtime_error("Error: InterfaceLength must be non-negative (" + context + ")");
          }
          if (interfaceLength == 0.0 && sectionIntervals != 0)
          {
            throw std::runtime_error("Error: a zero-length interface must have zero NumberOfGridPoints (" + context +
                                     ")");
          }
          sectionInterfaceLengths.push_back(interfaceLength);
          sectionInterfaceGridPoints.push_back(sectionIntervals);
          continue;
        }

        if (sectionInterfaceLengths.size() != nextAdsorbentInOrder)
        {
          throw std::runtime_error("Error: ColumnSections must define an interface between neighboring adsorbents (" +
                                   context + ")");
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
        sectionGridPoints[adsorbentIndex->second] = sectionIntervals;
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

      if (hasSectionGridPoints)
      {
        if (std::ranges::any_of(sectionGridPoints, [](size_t n) { return n == 0; }))
        {
          throw std::runtime_error("Error: ColumnSections NumberOfGridPoints must be set for every adsorbent section");
        }
        for (size_t interface = 0; interface < sectionInterfaceLengths.size(); ++interface)
        {
          if (sectionInterfaceLengths[interface] > 0.0 && sectionInterfaceGridPoints[interface] == 0)
          {
            throw std::runtime_error(
                "Error: ColumnSections NumberOfGridPoints must be set for every positive-length interface");
          }
        }
      }

      adsorbentLengths = std::move(sectionLengths);
      adsorbentGridPoints = std::move(sectionGridPoints);
      adsorbentInterfaceLengths = std::move(sectionInterfaceLengths);
      adsorbentInterfaceGridPoints = std::move(sectionInterfaceGridPoints);
    }
    else
    {
      adsorbentInterfaceLengths.assign(adsorbentsJson.size() - 1, 0.0);
      adsorbentInterfaceGridPoints.assign(adsorbentsJson.size() - 1, 0);
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

    for (size_t reactionId = 0; reactionId < reactionAdsorbentNames.size(); ++reactionId)
    {
      for (const std::string& name : reactionAdsorbentNames[reactionId])
      {
        const auto adsorbentIndex = adsorbentIndexByName.find(name);
        if (adsorbentIndex == adsorbentIndexByName.end())
        {
          throw std::runtime_error("Error: Reaction " + std::to_string(reactionId) + " names unknown adsorbent '" +
                                   name + "'");
        }
        reactions[reactionId].adsorbents.push_back(adsorbentIndex->second);
      }
    }

    if (!usesUniformAdsorbentMix && adsorbentInterfaceLengths.size() != adsorbentsJson.size() - 1)
    {
      throw std::runtime_error("Error: ColumnSections must define Adsorbents size minus one interfaces");
    }

    if (!usesUniformAdsorbentMix)
    {
      columnLength = std::reduce(adsorbentLengths.begin(), adsorbentLengths.end(), 0.0) +
                     std::reduce(adsorbentInterfaceLengths.begin(), adsorbentInterfaceLengths.end(), 0.0);
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
          z += adsorbentLengths[ads];
          columnDistances.back() = z;
          if (ads < adsorbentInterfaceLengths.size() && adsorbentInterfaceLengths[ads] > 0.0)
          {
            const double interfaceDz =
                adsorbentInterfaceLengths[ads] / static_cast<double>(adsorbentInterfaceGridPoints[ads]);
            for (size_t grid = 1; grid <= adsorbentInterfaceGridPoints[ads]; ++grid)
            {
              columnDistances.push_back(z + static_cast<double>(grid) * interfaceDz);
            }
            z += adsorbentInterfaceLengths[ads];
            columnDistances.back() = z;
          }
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
    if (std::any_of(reactionAdsorbentNames.begin(), reactionAdsorbentNames.end(),
                    [](const std::vector<std::string>& names) { return !names.empty(); }))
    {
      throw std::runtime_error("Error: Reaction Adsorbents requires a multibed column (an Adsorbents list)");
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

  // The component marked Water sets RH = p_H2O / p_sat(T) and supplies the water loading of the GAB-coadsorption
  // sites; a RelativeHumidity fixes RH instead.
  const size_t numberOfWaterComponents = static_cast<size_t>(std::count_if(
      components.begin(), components.end(), [](const Component& component) { return component.isWater; }));
  if (numberOfWaterComponents > 1)
  {
    throw std::runtime_error("Error: at most one component can be marked Water");
  }
  for (const Component& component : components)
  {
    if (component.isWater && component.isCarrierGas)
    {
      throw std::runtime_error("Error: the Water component '" + component.name + "' cannot be the CarrierGas");
    }
  }
  if (fluidPhase != 0 && (numberOfWaterComponents > 0 || relativeHumidity.has_value()))
  {
    throw std::runtime_error("Error: Water and RelativeHumidity are available only for gas simulations");
  }

  for (std::vector<Component>& componentsForAdsorbent : adsorbentComponents)
  {
    for (Component& comp : componentsForAdsorbent)
    {
      for (Isotherm& isotherm : comp.isotherm.sites) isotherm.waterSaturationPressure = waterSaturationPressure;
      for (Chemisorption& site : comp.chemisorption.sites)
      {
        if (site.isotherm.has_value()) site.isotherm->waterSaturationPressure = waterSaturationPressure;
      }
    }
  }
  for (Component& comp : components)
  {
    for (Isotherm& isotherm : comp.isotherm.sites) isotherm.waterSaturationPressure = waterSaturationPressure;
    for (Chemisorption& site : comp.chemisorption.sites)
    {
      if (site.isotherm.has_value()) site.isotherm->waterSaturationPressure = waterSaturationPressure;
    }
  }

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
        if (isotherm.type == Isotherm::Type::GAB && mixturePredictionMethod != 5)
        {
          throw std::runtime_error("Error: GAB requires MixturePredictionMethod SPI");
        }
        if (isotherm.type == Isotherm::Type::GAB && fluidPhase != 0)
        {
          throw std::runtime_error("Error: GAB is available only for gas simulations");
        }
        if (isotherm.type == Isotherm::Type::Langmuir_pH && fluidPhase != 1)
        {
          throw std::runtime_error("Error: pH-Langmuir is available only for liquid simulations");
        }
        if (isotherm.type == Isotherm::Type::Langmuir_pH && mixturePredictionMethod != 5)
        {
          throw std::runtime_error("Error: pH-Langmuir requires MixturePredictionMethod SPI");
        }
        if (Isotherm::isHumidityDependent(isotherm.type))
        {
          if (mixturePredictionMethod != 5)
          {
            throw std::runtime_error("Error: humidity-dependent isotherms require MixturePredictionMethod SPI");
          }
          if (fluidPhase != 0)
          {
            throw std::runtime_error("Error: humidity-dependent isotherms are available only for gas simulations");
          }
          if (comp.isWater)
          {
            throw std::runtime_error("Error: the Water component '" + comp.name +
                                     "' cannot have a humidity-dependent isotherm");
          }
          if (numberOfWaterComponents == 0 && !relativeHumidity.has_value())
          {
            throw std::runtime_error("Error: the humidity-dependent isotherm of '" + comp.name +
                                     "' needs a component with \"Water\": true or a RelativeHumidity");
          }
          if (Isotherm::usesWaterLoading(isotherm.type) && numberOfWaterComponents == 0)
          {
            throw std::runtime_error("Error: the GAB-coadsorption and GAB-modified isotherms of '" + comp.name +
                                     "' take the water loading from a component with \"Water\": true");
          }
        }
      }
      for (const Chemisorption& site : comp.chemisorption.sites)
      {
        if (site.isotherm.has_value() && Isotherm::isHumidityDependent(site.isotherm->type))
        {
          throw std::runtime_error("Error: humidity-dependent isotherms are available only as PhysisorptionSites (" +
                                   comp.name + ")");
        }
      }
    }
  }

  if (fluidPhase == 1 && mixturePredictionMethod == 6)
  {
    throw std::runtime_error("Error: MPD mixture prediction is available only for gas simulations");
  }

  // The pore reservoir has site 0 and one site per chemisorption site of any component and adsorbent
  // (ColumnData::maximumChemisorptionSites).
  size_t poreSites = 1;
  for (const std::vector<Component>& adsorbent : adsorbentComponents)
  {
    for (const Component& component : adsorbent) poreSites = std::max(poreSites, component.chemisorption.numberOfSites);
  }
  for (size_t reactionId = 0; reactionId < reactions.size(); ++reactionId)
  {
    const Reaction& reaction = reactions[reactionId];
    if (reaction.phase == Reaction::Phase::PoreConcentration && reaction.site >= poreSites)
    {
      throw std::runtime_error(std::format("Error: PoreConcentration reaction {} uses Site {}, but the pore reservoir "
                                           "has sites 0 to {} (one per chemisorption site)",
                                           reactionId, reaction.site, poreSites - 1));
    }
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
          if (site >= component.isotherm.sites.size() ||
              !component.isotherm.sites[site].enabled())
            continue;

          if (component.isCarrierGas && component.isotherm.sites[site].type == Isotherm::Type::Langmuir &&
              component.isotherm.sites[site].parameters[1] == 0.0)
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

    // Physical properties of the column, the bed, and the energy balance that the input leaves out take a default
    // (B-IO-11); a forgotten key must not pass unnoticed. A property counts as set when the top level, the Geometry,
    // or every Adsorbents entry gives it, or when it follows from other input (the column length of a stacked
    // multibed, the void fraction of a monolith, the temperature of the first swing phase).
    {
      const nlohmann::json* adsorbentList = findKeyCaseInsensitive(parsed_data, "Adsorbents");
      const bool hasAdsorbentList = adsorbentList != nullptr && adsorbentList->is_array() && !adsorbentList->empty();
      const auto everyAdsorbentSets = [&](const std::string& key)
      {
        return hasAdsorbentList &&
               std::all_of(adsorbentList->begin(), adsorbentList->end(), [&](const nlohmann::json& adsorbent)
                           { return containsKeyCaseInsensitive(adsorbent, key); });
      };
      const auto warnDefault = [](bool set, const std::string& key, double value, const std::string& unit)
      {
        if (!set) std::print(stderr, "Warning: {} is not set; using the default {:g}{}.\n", key, value, unit);
      };
      const nlohmann::json* geometryEntry = findKeyCaseInsensitive(parsed_data, "Geometry");
      const bool stackedMultibed = hasAdsorbentList && !everyAdsorbentSets("MixFraction");
      warnDefault(hasTemperature || (simulationType == SimulationType::SwingAdsorption &&
                                     !swingAdsorptionPhases.empty() &&
                                     swingAdsorptionPhases.front().temperature.has_value()),
                  "Temperature", temperature, " K");
      warnDefault(geometry.kind != GeometryKind::PackedBed ||
                      containsKeyCaseInsensitive(parsed_data, "ColumnVoidFraction") ||
                      (geometryEntry != nullptr && containsKeyCaseInsensitive(*geometryEntry, "ColumnVoidFraction")) ||
                      everyAdsorbentSets("ColumnVoidFraction"),
                  "ColumnVoidFraction", columnVoidFraction, "");
      warnDefault(containsKeyCaseInsensitive(parsed_data, "ParticleDensity") || everyAdsorbentSets("ParticleDensity"),
                  "ParticleDensity", particleDensity, " kg/m^3");
      warnDefault(containsKeyCaseInsensitive(parsed_data, "ColumnLength") || stackedMultibed, "ColumnLength",
                  columnLength, " m");
      if (energyBalance)
      {
        if (fluidPhase == 0)
        {
          warnDefault(containsKeyCaseInsensitive(parsed_data, "heatCapacityGas"), "heatCapacityGas", heatCapacityGas,
                      " J/(kg K)");
        }
        else
        {
          warnDefault(containsKeyCaseInsensitive(parsed_data, "heatCapacityLiquid"), "heatCapacityLiquid",
                      heatCapacityLiquid, " J/(kg K)");
        }
        warnDefault(containsKeyCaseInsensitive(parsed_data, "heatCapacitySolid"), "heatCapacitySolid",
                    heatCapacitySolid, " J/(kg K)");
        warnDefault(containsKeyCaseInsensitive(parsed_data, "heatCapacityWall"), "heatCapacityWall", heatCapacityWall,
                    " J/(kg K)");
      }
    }

    requireConfigured(temperature >= 0.0, "Error: temperature not set (Use e.g.: 'Temperature 300')");
    requireConfigured(columnVoidFraction > 0.0,
                      "Error: void-fraction of the column not set or invalid (Use e.g.: 'ColumnVoidFraction 0.4')");
    requireConfigured(particleDensity >= 0.0, "Error: particle density not set (Use e.g.: 'ParticleDensity 1408.2')");
    requireConfigured(std::isfinite(writeInterval) && writeInterval > 0.0,
                      "Error: WriteInterval must be a positive time in s (Use e.g.: 'WriteInterval 5')");
    if (breakthroughIntegrator == 0)
    {
      requireConfigured(std::isfinite(timeStep) && timeStep > 0.0,
                        "Error: TimeStep must be a positive time in s (Use e.g.: 'TimeStep 0.0005')");
      // The integrators land on every output time and phase end, so an interval that is not a whole number of
      // steps ends with a shortened step.
      const auto warnIfNotMultipleOfTimeStep = [&](double interval, const std::string& name)
      {
        const double steps = interval / timeStep;
        if (std::abs(steps - std::round(steps)) > 1.0e-6)
        {
          std::print(stderr,
                     "Warning: {} ({:g} s) is not a multiple of TimeStep ({:g} s); the last step of each such interval is "
                     "shortened.\n",
                     name, interval, timeStep);
        }
      };
      warnIfNotMultipleOfTimeStep(writeInterval, "WriteInterval");
      if (simulationType == SimulationType::Breakthrough && !autoSimulationTime)
      {
        warnIfNotMultipleOfTimeStep(simulationTime, "SimulationTime");
      }
      for (size_t phaseIndex = 0; phaseIndex < swingAdsorptionPhases.size(); ++phaseIndex)
      {
        warnIfNotMultipleOfTimeStep(swingAdsorptionPhases[phaseIndex].duration,
                                    "The Duration of SwingAdsorption phase " + std::to_string(phaseIndex + 1));
      }
    }
    if (simulationType == SimulationType::Breakthrough)
    {
      requireConfigured(autoSimulationTime || (std::isfinite(simulationTime) && simulationTime > 0.0),
                        "Error: SimulationTime must be a positive time in s or 'auto' (Use e.g.: 'SimulationTime 3600')");
    }
    else
    {
      requireConfigured(!swingAdsorptionPhases.empty(), "Error: SwingAdsorption requires SwingAdsorptionPhases");
      for (size_t phaseIndex = 0; phaseIndex < swingAdsorptionPhases.size(); ++phaseIndex)
      {
        const SwingAdsorptionPhase& phase = swingAdsorptionPhases[phaseIndex];
        requireConfigured(std::isfinite(phase.duration) && phase.duration > 0.0,
                          "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) +
                              " must have a positive Duration in s");
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
        if (phase.columnEntranceVelocity.has_value())
        {
          requireConfigured(
              std::isfinite(*phase.columnEntranceVelocity) && *phase.columnEntranceVelocity > 0.0,
              "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) + " has invalid ColumnEntranceVelocity");
          // These boundary conditions compute the velocity from the pressures, which would overwrite a phase value.
          requireConfigured(boundaryCondition != inletPressureOutletPressure &&
                                boundaryCondition != fixedPressureInletVelocity,
                            "Error: SwingAdsorption phase " + std::to_string(phaseIndex + 1) +
                                " sets ColumnEntranceVelocity, but the boundary condition derives the velocity from "
                                "the pressures");
        }
      }
    }
    requireConfigured(numberOfGridPoints != 0,
                      "Error: number of grid points not set (Use e.g.: 'NumberOfGridPoints 50')");
    requireConfigured(columnLength > 0.0, "Error: column length not set or invalid (Use e.g.: 'ColumnLength 0.3')");
    // The wall balance divides by the wall heat capacity per volume.
    requireConfigured(!energyBalance || (std::isfinite(wallDensity) && wallDensity > 0.0),
                      "Error: the energy balance requires wallDensity > 0 (Use e.g.: 'wallDensity 7800')");
    requireConfigured(!energyBalance || (std::isfinite(heatCapacityWall) && heatCapacityWall > 0.0),
                      "Error: the energy balance requires heatCapacityWall > 0 (Use e.g.: 'heatCapacityWall 500')");
    // The energy balance releases one physisorption heat per component. Without a HeatOfAdsorption it follows
    // from the energy of the affinity laws, which must then agree.
    if (energyBalance && !ambiguousHeats.empty())
    {
      std::string message =
          "Error: the energy balance needs one heat of adsorption per component, but the VantHoff affinity laws of "
          "these components have different energies and no HeatOfAdsorption is given:";
      for (const std::string& item : ambiguousHeats) message += "\n  " + item;
      throw std::runtime_error(message + "\nSet HeatOfAdsorption on the component; the laws keep their own energies.");
    }

    const bool hasBulkReaction =
        std::any_of(reactions.begin(), reactions.end(), [](const Reaction& reaction)
                    { return reaction.phase == Reaction::Phase::BulkConcentration; });
    // A bulk reaction releases its heat into the solid, per kg of it (computeReactionDerivatives).
    requireConfigured(!hasBulkReaction || !energyBalance ||
                          std::all_of(adsorbentParticleDensities.begin(), adsorbentParticleDensities.end(),
                                      [](double density) { return density > 0.0; }),
                      "Error: BulkConcentration reactions release their heat into the solid, so the energy balance "
                      "needs ParticleDensity > 0");

    // The axial-dispersion reference state enables the gas-phase (Fuller) scaling of the molecular part.
    const bool hasDispersionReference =
        axialDispersionReferenceTemperature != 0.0 || axialDispersionReferencePressure != 0.0;
    requireConfigured(!hasDispersionReference ||
                          (std::isfinite(axialDispersionReferenceTemperature) &&
                           axialDispersionReferenceTemperature > 0.0 &&
                           std::isfinite(axialDispersionReferencePressure) && axialDispersionReferencePressure > 0.0),
                      "Error: AxialDispersionReferenceTemperature and AxialDispersionReferencePressure must be given "
                      "together and be positive");
    requireConfigured(!hasDispersionReference || fluidPhase == 0,
                      "Error: the axial-dispersion reference state scales gas diffusivities; it does not apply to "
                      "FluidPhase Liquid");
    // The Taylor-Aris term divides by the molecular part, which is then the molecular diffusivity of each component.
    for (size_t ads = 0; ads < adsorbentGeometries.size(); ++ads)
    {
      if (adsorbentGeometries[ads].dispersion.taylorArisCoefficient <= 0.0) continue;
      for (const Component& component : adsorbentComponents[ads])
      {
        requireConfigured(component.axialDispersionCoefficient > 0.0,
                          "Error: TaylorArisDispersion needs a positive AxialDispersionCoefficient (the molecular "
                          "diffusivity) for every component; '" +
                              component.name + "' has none");
      }
    }

    const bool hasChemisorption = std::any_of(components.begin(), components.end(), [](const Component& component)
                                              { return component.chemisorption.numberOfSites > 0; });
    requireConfigured(mixturePredictionMethod != 6 || cvodePreconditioner == 0,
                      "Error: CVODEPreconditioner SorptionJacobian does not support MPD mixture prediction");
    // The KLU Jacobian covers the concentrations and physisorbed loadings only.
    const bool kluJacobian = breakthroughIntegrator == 1 && cvodeLinearSolver == 2;
    requireConfigured(!kluJacobian || (reactions.empty() && !hasChemisorption && !energyBalance),
                      "Error: CVODELinearSolver KLU does not support reactions, chemisorption, or the energy balance; "
                      "use SPGMR or Dense");
    requireConfigured(!kluJacobian || mixturePredictionMethod != 6,
                      "Error: CVODELinearSolver KLU does not support MPD mixture prediction");
    requireConfigured(adsorbentComponents.size() <= 1 || mixturePredictionMethod != 6,
                      "Error: MPD mixture prediction (MPDSettings) supports a single adsorbent");

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
    else if (!hasColumnEntranceVelocity)
    {
      // The output files scale the time with the input velocity, tau = t v_in / L (B-IO-9). These boundary conditions
      // compute the inlet velocity, which changes during the run, so tau needs a reference velocity from the input.
      std::print(stderr,
                 "Warning: ColumnEntranceVelocity is not set; the dimensionless time tau = t v_in / L (column 1 of the "
                 "output files) has no reference velocity and is not meaningful. Set ColumnEntranceVelocity to the "
                 "velocity that tau should refer to.\n");
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

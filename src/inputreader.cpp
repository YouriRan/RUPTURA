#include "inputreader.h"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "json.h"

bool caseInSensStringCompare(const std::string& str1, const std::string& str2)
{
  return str1.size() == str2.size() && std::equal(str1.begin(), str1.end(), str2.begin(),
                                                  [](int a, int b) { return std::tolower(a) == std::tolower(b); });
}

bool startsWith(const std::string& str, const std::string& prefix)
{
  return str.size() >= prefix.size() && str.substr(0, prefix.size()) == prefix;
}

std::string trim(const std::string& s)
{
  auto start = s.begin();
  while (start != s.end() && std::isspace(*start))
  {
    start++;
  }

  auto end = s.end();
  do
  {
    end--;
  } while (std::distance(start, end) > 0 && std::isspace(*end));

  return std::string(start, end + 1);
}

template <class T>
T parse(const std::string& arguments, [[maybe_unused]] const std::string& keyword, [[maybe_unused]] size_t lineNumber)
{
  T value;

  std::string str;
  std::istringstream ss(arguments);

  ss >> value;

  return value;
}

template <typename T>
std::vector<T> parseListOfSystemValues(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  std::vector<T> list{};

  std::string str;
  std::istringstream ss(arguments);

  std::string errorString =
      "No values could be read for keyword '" + keyword + "' at line: " + std::to_string(lineNumber) + "\n";

  while (ss >> str)
  {
    if (trim(str).rfind("//", 0) == 0)
    {
      if (list.empty())
      {
        throw std::runtime_error(errorString);
      }
      return list;
    }
    T value;
    std::istringstream s(str);
    if (s >> value)
    {
      list.push_back(value);
    }
    else
    {
      if (list.empty())
      {
        throw std::runtime_error(errorString);
      }
      return list;
    }
  };

  if (list.empty())
  {
    throw std::runtime_error(errorString);
  }
  return list;
}

double parseDouble(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  double value{};

  std::string str;
  std::istringstream ss(arguments);

  if (ss >> value)
  {
    return value;
  };

  std::string errorString =
      "Numbers could not be read for keyword '" + keyword + "' at line: " + std::to_string(lineNumber) + "\n";
  throw std::runtime_error(errorString);
}

int parseBoolean(const std::string& arguments, const std::string& keyword, size_t lineNumber)
{
  bool value{};

  std::istringstream ss(arguments);

  if (ss >> std::boolalpha >> value)
  {
    return value;
  };

  std::string str;
  std::istringstream ss2(arguments);
  if (ss2 >> str)
  {
    if (caseInSensStringCompare(str, "yes")) return true;
    if (caseInSensStringCompare(str, "no")) return false;
  };

  std::string errorString =
      "Booleands could not be read for keyword '" + keyword + "' at line: " + std::to_string(lineNumber) + "\n";
  throw std::runtime_error(errorString);
}

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

  // Accept legacy-like encodings: "yes"/"no".
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

static void addIsothermSiteFromJson(Component& component, const std::string& typeString, const nlohmann::json& params,
                                    const std::string& context)
{
  if (!params.is_array())
  {
    throw std::runtime_error("Error: isotherm parameters must be an array" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  std::vector<double> values{};
  try
  {
    values = params.get<std::vector<double>>();
  }
  catch (const nlohmann::json::exception&)
  {
    throw std::runtime_error("Error: isotherm parameters must be numeric array" +
                             (context.empty() ? "" : (" (" + context + ")")));
  }

  auto requireAtLeast = [&](std::size_t n, const std::string& name)
  {
    if (values.size() < n)
    {
      throw std::runtime_error("Error: " + name + " requires " + std::to_string(n) + " parameters" +
                               (context.empty() ? "" : (" (" + context + ")")));
    }
    values.resize(n);
  };

  if (caseInSensStringCompare(typeString, "Langmuir"))
  {
    requireAtLeast(2, "Langmuir");
    component.isotherm.add(Isotherm(Isotherm::Type::Langmuir, values, 2));
    return;
  }
  if (caseInSensStringCompare(typeString, "Anti-Langmuir") || caseInSensStringCompare(typeString, "Anti_Langmuir"))
  {
    requireAtLeast(2, "Anti-Langmuir");
    component.isotherm.add(Isotherm(Isotherm::Type::Anti_Langmuir, values, 2));
    return;
  }
  if (caseInSensStringCompare(typeString, "BET"))
  {
    requireAtLeast(3, "BET");
    component.isotherm.add(Isotherm(Isotherm::Type::BET, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Henry"))
  {
    requireAtLeast(1, "Henry");
    component.isotherm.add(Isotherm(Isotherm::Type::Henry, values, 1));
    return;
  }
  if (caseInSensStringCompare(typeString, "Freundlich"))
  {
    requireAtLeast(2, "Freundlich");
    component.isotherm.add(Isotherm(Isotherm::Type::Freundlich, values, 2));
    return;
  }
  if (caseInSensStringCompare(typeString, "Sips"))
  {
    requireAtLeast(3, "Sips");
    component.isotherm.add(Isotherm(Isotherm::Type::Sips, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Langmuir-Freundlich") ||
      caseInSensStringCompare(typeString, "Langmuir_Freundlich"))
  {
    requireAtLeast(3, "Langmuir-Freundlich");
    component.isotherm.add(Isotherm(Isotherm::Type::Langmuir_Freundlich, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Redlich-Peterson") ||
      caseInSensStringCompare(typeString, "Redlich_Peterson"))
  {
    requireAtLeast(3, "Redlich-Peterson");
    component.isotherm.add(Isotherm(Isotherm::Type::Redlich_Peterson, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Toth"))
  {
    requireAtLeast(3, "Toth");
    component.isotherm.add(Isotherm(Isotherm::Type::Toth, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Unilan"))
  {
    requireAtLeast(3, "Unilan");
    component.isotherm.add(Isotherm(Isotherm::Type::Unilan, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "O'Brian&Myers") || caseInSensStringCompare(typeString, "OBrien_Myers") ||
      caseInSensStringCompare(typeString, "OBrien&Myers"))
  {
    requireAtLeast(3, "O'Brien&Myers");
    component.isotherm.add(Isotherm(Isotherm::Type::OBrien_Myers, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Quadratic"))
  {
    requireAtLeast(3, "Quadratic");
    component.isotherm.add(Isotherm(Isotherm::Type::Quadratic, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Temkin"))
  {
    requireAtLeast(3, "Temkin");
    component.isotherm.add(Isotherm(Isotherm::Type::Temkin, values, 3));
    return;
  }
  if (caseInSensStringCompare(typeString, "Bingel&Walton") || caseInSensStringCompare(typeString, "BingelWalton"))
  {
    requireAtLeast(3, "Bingel&Walton");
    component.isotherm.add(Isotherm(Isotherm::Type::BingelWalton, values, 3));
    return;
  }

  throw std::runtime_error("Error: unknown isotherm type '" + typeString + "'" +
                           (context.empty() ? "" : (" (" + context + ")")));
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

  nlohmann::json parsed_data{};
  try
  {
    parsed_data = nlohmann::json::parse(fileInput);
  }
  catch (const nlohmann::json::parse_error& ex)
  {
    std::cerr << "parse error at byte " << ex.byte << std::endl;
    throw;
  }

  // Top-level options
  if (containsKeyCaseInsensitive(parsed_data, "SimulationType"))
  {
    std::string str =
        getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "SimulationType", ""), "SimulationType", "");
    if (caseInSensStringCompare(str, "Breakthrough"))
      simulationType = SimulationType::Breakthrough;
    else if (caseInSensStringCompare(str, "MixturePrediction"))
      simulationType = SimulationType::MixturePrediction;
    else if (caseInSensStringCompare(str, "Fitting"))
      simulationType = SimulationType::Fitting;
    else if (caseInSensStringCompare(str, "Test"))
      simulationType = SimulationType::Test;
    else
      throw std::runtime_error("Error: invalid SimulationType '" + str + "'");
  }

  if (containsKeyCaseInsensitive(parsed_data, "MixturePredictionMethod"))
  {
    std::string str = getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "MixturePredictionMethod", ""),
                                       "MixturePredictionMethod", "");
    if (caseInSensStringCompare(str, "IAST"))
      mixturePredictionMethod = 0;
    else if (caseInSensStringCompare(str, "SIAST"))
      mixturePredictionMethod = 1;
    else if (caseInSensStringCompare(str, "EI"))
      mixturePredictionMethod = 2;
    else if (caseInSensStringCompare(str, "SEI"))
      mixturePredictionMethod = 3;
    else
      throw std::runtime_error("Error: invalid MixturePredictionMethod '" + str + "'");
  }

  if (containsKeyCaseInsensitive(parsed_data, "IASTMethod"))
  {
    std::string str = getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "IASTMethod", ""), "IASTMethod", "");
    if (caseInSensStringCompare(str, "FastIAS"))
      IASTMethod = 0;
    else if (caseInSensStringCompare(str, "Bisection"))
      IASTMethod = 1;
    else
      throw std::runtime_error("Error: invalid IASTMethod '" + str + "'");
  }

  if (containsKeyCaseInsensitive(parsed_data, "BreakthroughIntegrator"))
  {
    std::string str = getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "BreakthroughIntegrator", ""),
                                       "BreakthroughIntegrator", "");
    if (caseInSensStringCompare(str, "RungeKutta3"))
      breakthroughIntegrator = 0;
    else if (caseInSensStringCompare(str, "CVODE"))
      breakthroughIntegrator = 1;
    else if (caseInSensStringCompare(str, "SIRK3"))
      breakthroughIntegrator = 3;
    else
      throw std::runtime_error("Error: invalid BreakthroughIntegrator '" + str + "'");
  }

  if (containsKeyCaseInsensitive(parsed_data, "VelocityProfile"))
  {
    std::string str =
        getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "VelocityProfile", ""), "VelocityProfile", "");
    if (caseInSensStringCompare(str, "FixedPressureGradient"))
      velocityProfile = 0;
    else if (caseInSensStringCompare(str, "Ergun"))
      velocityProfile = 1;
    else if (caseInSensStringCompare(str, "FixedVelocity"))
      velocityProfile = 2;
    else
      throw std::runtime_error("Error: invalid VelocityProfile '" + str + "'");
  }

  if (containsKeyCaseInsensitive(parsed_data, "ReadColumnFile"))
  {
    readColumnFile =
        getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "ReadColumnFile", ""), "ReadColumnFile", "");
  }

  if (containsKeyCaseInsensitive(parsed_data, "DisplayName"))
  {
    displayName = getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "DisplayName", ""), "DisplayName", "");
  }

  if (containsKeyCaseInsensitive(parsed_data, "Temperature"))
  {
    temperature =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "Temperature", ""), "Temperature", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ColumnVoidFraction"))
  {
    columnVoidFraction = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ColumnVoidFraction", ""),
                                                  "ColumnVoidFraction", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "DynamicViscosity"))
  {
    dynamicViscosity = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "DynamicViscosity", ""),
                                                "DynamicViscosity", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ParticleDiameter"))
  {
    particleDiameter = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ParticleDiameter", ""),
                                                "ParticleDiameter", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ParticleDensity"))
  {
    particleDensity =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ParticleDensity", ""), "ParticleDensity", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "TotalPressure"))
  {
    totalPressure =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "TotalPressure", ""), "TotalPressure", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "PressureStart"))
  {
    pressureStart =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "PressureStart", ""), "PressureStart", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "PressureEnd"))
  {
    pressureEnd =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "PressureEnd", ""), "PressureEnd", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "NumberOfPressurePoints"))
  {
    numberOfPressurePoints = getNumberOrThrow<std::size_t>(
        requireKeyCaseInsensitive(parsed_data, "NumberOfPressurePoints", ""), "NumberOfPressurePoints", "");
  }

  if (containsKeyCaseInsensitive(parsed_data, "PressureScale"))
  {
    std::string str =
        getStringOrThrow(requireKeyCaseInsensitive(parsed_data, "PressureScale", ""), "PressureScale", "");
    if (caseInSensStringCompare(str, "Log"))
      pressureScale = 0;
    else if (caseInSensStringCompare(str, "Linear") || caseInSensStringCompare(str, "Normal"))
      pressureScale = 1;
    else
      throw std::runtime_error("Error: invalid PressureScale '" + str + "'");
  }

  if (containsKeyCaseInsensitive(parsed_data, "PressureGradient"))
  {
    pressureGradient = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "PressureGradient", ""),
                                                "PressureGradient", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ColumnEntranceVelocity"))
  {
    columnEntranceVelocity = getNumberOrThrow<double>(
        requireKeyCaseInsensitive(parsed_data, "ColumnEntranceVelocity", ""), "ColumnEntranceVelocity", "");
  }

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

  if (containsKeyCaseInsensitive(parsed_data, "NumberOfInitTimeSteps"))
  {
    numberOfInitTimeSteps = getNumberOrThrow<std::size_t>(
        requireKeyCaseInsensitive(parsed_data, "NumberOfInitTimeSteps", ""), "NumberOfInitTimeSteps", "");
  }

  if (containsKeyCaseInsensitive(parsed_data, "PulseBreakthrough"))
  {
    pulseBreakthrough =
        getBoolOrThrow(requireKeyCaseInsensitive(parsed_data, "PulseBreakthrough", ""), "PulseBreakthrough", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "PulseTime"))
  {
    pulseTime = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "PulseTime", ""), "PulseTime", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "TimeStep"))
  {
    timeStep = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "TimeStep", ""), "TimeStep", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "PrintEvery"))
  {
    printEvery =
        getNumberOrThrow<std::size_t>(requireKeyCaseInsensitive(parsed_data, "PrintEvery", ""), "PrintEvery", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "WriteEvery"))
  {
    writeEvery =
        getNumberOrThrow<std::size_t>(requireKeyCaseInsensitive(parsed_data, "WriteEvery", ""), "WriteEvery", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ColumnLength"))
  {
    columnLength =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "ColumnLength", ""), "ColumnLength", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "NumberOfGridPoints"))
  {
    numberOfGridPoints = getNumberOrThrow<std::size_t>(requireKeyCaseInsensitive(parsed_data, "NumberOfGridPoints", ""),
                                                       "NumberOfGridPoints", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ColumnPressure"))
  {
    columnPressure = getNumberOrThrow<std::size_t>(requireKeyCaseInsensitive(parsed_data, "ColumnPressure", ""),
                                                   "ColumnPressure", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ColumnLoading"))
  {
    columnLoading =
        getNumberOrThrow<std::size_t>(requireKeyCaseInsensitive(parsed_data, "ColumnLoading", ""), "ColumnLoading", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "ColumnError"))
  {
    columnError =
        getNumberOrThrow<std::size_t>(requireKeyCaseInsensitive(parsed_data, "ColumnError", ""), "ColumnError", "");
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

  if (containsKeyCaseInsensitive(parsed_data, "internalDiameter"))
  {
    internalDiameter = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "internalDiameter", ""),
                                                "internalDiameter", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "outerDiameter"))
  {
    outerDiameter =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "outerDiameter", ""), "outerDiameter", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "wallDensity"))
  {
    wallDensity =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "wallDensity", ""), "wallDensity", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "gasThermalConductivity"))
  {
    gasThermalConductivity = getNumberOrThrow<double>(
        requireKeyCaseInsensitive(parsed_data, "gasThermalConductivity", ""), "gasThermalConductivity", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "wallThermalConductivity"))
  {
    wallThermalConductivity = getNumberOrThrow<double>(
        requireKeyCaseInsensitive(parsed_data, "wallThermalConductivity", ""), "wallThermalConductivity", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "heatTransferGasSolid"))
  {
    heatTransferGasSolid = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "heatTransferGasSolid", ""),
                                                    "heatTransferGasSolid", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "heatTransferGasWall"))
  {
    heatTransferGasWall = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "heatTransferGasWall", ""),
                                                   "heatTransferGasWall", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "heatTransferWallExternal"))
  {
    heatTransferWallExternal = getNumberOrThrow<double>(
        requireKeyCaseInsensitive(parsed_data, "heatTransferWallExternal", ""), "heatTransferWallExternal", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "heatCapacityGas"))
  {
    heatCapacityGas =
        getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "heatCapacityGas", ""), "heatCapacityGas", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "heatCapacitySolid"))
  {
    heatCapacitySolid = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "heatCapacitySolid", ""),
                                                 "heatCapacitySolid", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "heatCapacityWall"))
  {
    heatCapacityWall = getNumberOrThrow<double>(requireKeyCaseInsensitive(parsed_data, "heatCapacityWall", ""),
                                                "heatCapacityWall", "");
  }
  if (containsKeyCaseInsensitive(parsed_data, "energyBalance"))
  {
    energyBalance = getBoolOrThrow(requireKeyCaseInsensitive(parsed_data, "energyBalance", ""), "energyBalance", "");
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

      // Use 'Name' for the component name.
      std::string componentName = getStringOrThrow(requireKeyCaseInsensitive(item, "Name", context), "Name", context);

      // Preserve legacy behavior: componentId is the position in the list.
      if (componentId == components.size())
      {
        components.push_back(Component(componentId, componentName));
      }
      else
      {
        components[componentId] = Component(componentId, componentName);
      }

      if (containsKeyCaseInsensitive(item, "FileName"))
      {
        components[componentId].filename =
            getStringOrThrow(requireKeyCaseInsensitive(item, "FileName", context), "FileName", context);
      }
      if (containsKeyCaseInsensitive(item, "CarrierGas"))
      {
        components[componentId].isCarrierGas =
            getBoolOrThrow(requireKeyCaseInsensitive(item, "CarrierGas", context), "CarrierGas", context);
      }
      if (containsKeyCaseInsensitive(item, "GasPhaseMolFraction"))
      {
        components[componentId].Yi0 = getNumberOrThrow<double>(
            requireKeyCaseInsensitive(item, "GasPhaseMolFraction", context), "GasPhaseMolFraction", context);
      }
      if (containsKeyCaseInsensitive(item, "MassTransferCoefficient"))
      {
        components[componentId].Kl = getNumberOrThrow<double>(
            requireKeyCaseInsensitive(item, "MassTransferCoefficient", context), "MassTransferCoefficient", context);
      }
      if (containsKeyCaseInsensitive(item, "AxialDispersionCoefficient"))
      {
        components[componentId].D =
            getNumberOrThrow<double>(requireKeyCaseInsensitive(item, "AxialDispersionCoefficient", context),
                                     "AxialDispersionCoefficient", context);
      }
      if (containsKeyCaseInsensitive(item, "MolecularWeight"))
      {
        components[componentId].molecularWeight = getNumberOrThrow<double>(
            requireKeyCaseInsensitive(item, "MolecularWeight", context), "MolecularWeight", context);
      }
      if (containsKeyCaseInsensitive(item, "HeatOfAdsorption"))
      {
        components[componentId].heatOfAdsorption = getNumberOrThrow<double>(
            requireKeyCaseInsensitive(item, "HeatOfAdsorption", context), "HeatOfAdsorption", context);
      }

      if (containsKeyCaseInsensitive(item, "NumberOfIsothermSites"))
      {
        components[componentId].isotherm.numberOfSites = getNumberOrThrow<std::size_t>(
            requireKeyCaseInsensitive(item, "NumberOfIsothermSites", context), "NumberOfIsothermSites", context);
      }

      // Preferred JSON format: "IsothermSites": [ {"Type":"Langmuir", "Parameters":[...]}, ... ]
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
            addIsothermSiteFromJson(components[componentId], typeString, params, siteContext);
          }
          else if (site.is_object() && site.size() == 1)
          {
            // Alternative compact encoding: {"Langmuir": [..]}
            auto it = site.begin();
            addIsothermSiteFromJson(components[componentId], it.key(), it.value(), siteContext);
          }
          else
          {
            throw std::runtime_error("Error: invalid IsothermSites entry (" + siteContext + ")");
          }
        }
      }
      else
      {
        // Backward-compatible JSON encoding: allow direct isotherm keys at the component level.
        for (auto it = item.begin(); it != item.end(); ++it)
        {
          const std::string& k = it.key();
          if (caseInSensStringCompare(k, "Langmuir") || caseInSensStringCompare(k, "Anti-Langmuir") ||
              caseInSensStringCompare(k, "Anti_Langmuir") || caseInSensStringCompare(k, "BET") ||
              caseInSensStringCompare(k, "Henry") || caseInSensStringCompare(k, "Freundlich") ||
              caseInSensStringCompare(k, "Sips") || caseInSensStringCompare(k, "Langmuir-Freundlich") ||
              caseInSensStringCompare(k, "Langmuir_Freundlich") || caseInSensStringCompare(k, "Redlich-Peterson") ||
              caseInSensStringCompare(k, "Redlich_Peterson") || caseInSensStringCompare(k, "Toth") ||
              caseInSensStringCompare(k, "Unilan") || caseInSensStringCompare(k, "O'Brian&Myers") ||
              caseInSensStringCompare(k, "OBrien_Myers") || caseInSensStringCompare(k, "OBrien&Myers") ||
              caseInSensStringCompare(k, "Quadratic") || caseInSensStringCompare(k, "Temkin") ||
              caseInSensStringCompare(k, "Bingel&Walton") || caseInSensStringCompare(k, "BingelWalton"))
          {
            addIsothermSiteFromJson(components[componentId], k, it.value(), context);
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
      // Allow object form { "0": {...}, "1": {...} } but keep IDs in insertion order.
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

  // normalize gas-phase mol-fractions to unity
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
    if (numberOfCarrierGases == 0)
    {
      throw std::runtime_error("Error: no carrier gas component present");
    }
    if (numberOfCarrierGases > 1)
    {
      throw std::runtime_error("Error: multiple carrier gas component present (there can be only one)");
    }
    if (temperature < 0.0)
    {
      throw std::runtime_error("Error: temperature not set (Use e.g.: 'Temperature 300'");
    }
    if (columnVoidFraction < 0.0)
    {
      throw std::runtime_error("Error: void-fraction of the colum not set (Use e.g.: 'ColumnVoidFraction 0.4'");
    }
    if (particleDensity < 0.0)
    {
      throw std::runtime_error("Error: particle density not set (Use e.g.: 'ParticleDensity 1408.2'");
    }
    if (totalPressure < 0.0)
    {
      throw std::runtime_error("Error: total pressure bot set (Use e.g.: 'TotalPressure 1e5'");
    }
    if (columnEntranceVelocity < 0.0)
    {
      throw std::runtime_error("Error: column entrance velocity not set (Use e.g.: 'columnEntranceVelocity 300'");
    }

    if ((numberOfTimeSteps == 0) && (!autoNumberOfTimeSteps))
    {
      throw std::runtime_error("Error: number of time steps not set (Use e.g.: 'NumberOfTimeSteps 5000000'");
    }
    if (numberOfGridPoints == 0)
    {
      throw std::runtime_error("Error: number of grid points not set (Use e.g.: 'NumberOfGridPoints 50'");
    }
    if (columnLength < 0)
    {
      throw std::runtime_error("Error: column length not set (Use e.g.: 'ColumnLength 0.3'");
    }
  }
}

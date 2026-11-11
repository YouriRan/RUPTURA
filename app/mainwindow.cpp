#include "mainwindow.h"

#include <QAbstractItemView>
#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QFontMetrics>
#include <QFormLayout>
#include <QFrame>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSizePolicy>
#include <QSplitter>
#include <QStatusBar>
#include <QStringList>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QTextCursor>
#include <QTimer>
#include <QVBoxLayout>
#include <algorithm>
#include <limits>

#include "simulationinput.h"
#include "inputreader.h"
#include "isotherm.h"
#include "temperature_law.h"

namespace
{
constexpr const char* stateFormat = "RupturaLabState";
constexpr qsizetype savedLogLimit = 12000;

QJsonObject selectSettings(const QJsonObject& object, const QStringList& keys)
{
  QJsonObject result;
  for (const auto& key : keys)
    if (object.contains(key)) result.insert(key, object.value(key));
  return result;
}

const QStringList columnSettingKeys{
    "FluidPhase", "LiquidDensity", "pHMode", "pHValue", "pKw", "RelativeHumidity", "WaterSaturationPressureModel",
    "WaterSaturationPressure", "liquidThermalConductivity",
    "heatTransferLiquidSolid", "heatTransferLiquidWall", "heatCapacityLiquid", "AxialDispersionReferenceTemperature",
    "AxialDispersionReferencePressure"};
const QStringList simulationSettingKeys{
    "IASTMethod", "ConvectionScheme", "InletCondition", "FlowDirection", "CVODERelativeTolerance",
    "CVODEAbsoluteToleranceConcentration",
    "CVODEAbsoluteToleranceLoading",
    "CVODEAbsoluteToleranceTemperature", "CVODEMaximumTimeStep", "CVODEMaximumOrder", "CVODELinearSolver",
    "CVODEKrylovDimension", "ROS34PW2RelativeTolerance", "ROS34PW2AbsoluteToleranceConcentration",
    "ROS34PW2AbsoluteToleranceLoading", "ROS34PW2AbsoluteToleranceTemperature", "ROS34PW2MaximumTimeStep",
    "ColumnDistances", "ReadColumnFile", "DebugForceMultibed", "ColumnPressure",
    "ColumnLoading", "ColumnError", "FittingWeighting"};
const QStringList componentSettingKeys{"LiquidPhaseConcentration", "InitialLiquidPhaseConcentration", "FileName"};

void mergeSettings(QJsonObject& object, const QJsonObject& settings)
{
  for (auto it = settings.begin(); it != settings.end(); ++it) object.insert(it.key(), it.value());
}

QString choice(const QString& value, const QStringList& values)
{
  for (const auto& item : values)
    if (value.compare(item, Qt::CaseInsensitive) == 0) return item;
  return value;
}

QString numberText(double value) { return QString::number(value, 'g', 12); }

QString optionalText(const std::optional<double>& value) { return value ? numberText(*value) : QString{}; }

QString isoDate(const QDateTime& value)
{
  return value.isValid() ? value.toUTC().toString(Qt::ISODateWithMs) : QString{};
}

QDateTime dateFromIso(const QString& value)
{
  if (value.isEmpty())
  {
    return {};
  }
  QDateTime date = QDateTime::fromString(value, Qt::ISODateWithMs);
  if (!date.isValid())
  {
    date = QDateTime::fromString(value, Qt::ISODate);
  }
  return date;
}

QVariantMap defaultsFor(const QList<ParameterDefinition>& definitions, const QVariantMap& values = {})
{
  QVariantMap result;
  for (const ParameterDefinition& definition : definitions)
  {
    if (values.contains(definition.key))
    {
      result.insert(definition.key, values.value(definition.key));
    }
    else if (definition.boolean)
    {
      result.insert(definition.key, definition.defaultValue != 0.0);
    }
    else
    {
      result.insert(definition.key, definition.defaultValue);
    }
  }
  return result;
}

// A temperature-scalable parameter holds a number or a law map {"Value", "Law", <coefficient>, "ReferenceTemperature"}.
double parameterNumber(const QVariantMap& parameters, const ParameterDefinition& definition)
{
  const QVariant value = parameters.value(definition.key, definition.defaultValue);
  return value.typeId() == QMetaType::QVariantMap ? value.toMap().value("Value").toDouble() : value.toDouble();
}

bool parameterBool(const QVariantMap& parameters, const ParameterDefinition& definition)
{
  return parameters.value(definition.key, definition.defaultValue != 0.0).toBool();
}

QJsonValue parameterJson(const QVariantMap& parameters, const ParameterDefinition& definition)
{
  if (definition.boolean) return parameterBool(parameters, definition);
  const QVariant value = parameters.value(definition.key, definition.defaultValue);
  if (value.typeId() == QMetaType::QVariantMap) return QJsonObject::fromVariantMap(value.toMap());
  return value.toDouble();
}

void clearLayout(QLayout* layout)
{
  while (layout->count() > 0)
  {
    QLayoutItem* item = layout->takeAt(0);
    if (item == nullptr)
    {
      break;
    }
    if (QWidget* widget = item->widget())
    {
      widget->deleteLater();
    }
    if (QLayout* childLayout = item->layout())
    {
      clearLayout(childLayout);
      childLayout->deleteLater();
    }
    delete item;
  }
}

QList<QPair<QString, QString>> valueOptions(const QStringList& values)
{
  QList<QPair<QString, QString>> options;
  for (const QString& value : values)
  {
    options.push_back({value, value});
  }
  return options;
}

void sizeComboForOptions(QComboBox* combo, const QList<QPair<QString, QString>>& options)
{
  int maxTextWidth = 0;
  const QFontMetrics metrics(combo->font());
  for (const auto& option : options)
  {
    maxTextWidth = std::max(maxTextWidth, metrics.horizontalAdvance(option.first));
    maxTextWidth = std::max(maxTextWidth, metrics.horizontalAdvance(option.second));
  }

  const int popupWidth = std::max(260, maxTextWidth + 56);
  const int contentChars = std::clamp((maxTextWidth / std::max(1, metrics.averageCharWidth())) + 4, 18, 42);
  combo->setMinimumContentsLength(contentChars);
  combo->setMinimumWidth(std::min(360, popupWidth));
  combo->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  combo->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
  combo->setMaxVisibleItems(16);
  combo->view()->setMinimumWidth(popupWidth);
}

QList<QPair<QString, QString>> physisorptionTypeOptions()
{
  return valueOptions({"Langmuir", "Anti-Langmuir", "BET", "Henry", "Freundlich", "Sips", "Langmuir-Freundlich",
                       "Redlich-Peterson", "Toth", "Unilan", "OBrien&Myers", "Quadratic", "Temkin",
                       "Bingel&Walton", "pH-Langmuir", "GAB", "Humidity-enhanced-Langmuir",
                       "Moisture-swing-Langmuir", "GAB-coadsorption-Langmuir", "GAB-coadsorption-Sips",
                       "GAB-coadsorption-Toth", "GAB-modified-Toth", "Exponential-humidity-Toth"});
}

QList<QPair<QString, QString>> chemisorptionTypeOptions()
{
  return valueOptions({"FirstOrder", "PseudoNth", "Avrami", "General", "Elovich"});
}

QList<QPair<QString, QString>> reactionPhaseOptions()
{
  return valueOptions({"Physisorbed", "Chemisorbed", "PoreConcentration", "BulkConcentration"});
}

QList<QPair<QString, QString>> reactionStyleOptions()
{
  return valueOptions({"GeneralPowerLaw", "LangmuirHinshelwood", "LangmuirHinshelwoodHougenWatson"});
}

QList<ParameterDefinition> phyDefs(const QString& type)
{
  QList<ParameterDefinition> definitions = [&type]() -> QList<ParameterDefinition>
  {
    if (type == "pH-Langmuir")
      return {{"qmax", "Maximum loading [mol/kg]", 1.0}, {"affinityConstant", "Affinity [m3/mol]", 0.5},
              {"pKa", "pKa [-]", 7.0}};
    if (type == "Bingel&Walton")
      return {{"qmax", "Saturation loading [mol/kg]", 1.0}, {"a", "a [1/Pa]", 1.0e-5}, {"b", "b [1/Pa]", 1.0e-5}};
    if (type == "GAB")
      return {{"qmax", "Monolayer loading [mol/kg]", 1.0}, {"c0", "C [-]", 1.0}, {"k0", "K [-]", 0.5}};
    // The humidity-dependent types: RH comes from the Water component or the fixed RelativeHumidity, the water
    // loading q_w from the isotherm of the Water component.
    if (type == "Humidity-enhanced-Langmuir")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Dry affinity b [1/Pa]", 1.0e-5},
              {"humidityEnhancement", "Humidity enhancement alpha [-]", 1.0},
              {"waterAffinity", "Water affinity k_w [-]", 1.0}};
    }
    if (type == "Moisture-swing-Langmuir")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity at RH = 1, b [1/Pa]", 1.0e-5},
              {"humidityExponent", "Humidity exponent n_w [-]", 1.0}};
    }
    if (type == "GAB-coadsorption-Langmuir")
    {
      return {{"saturationLoading", "Dry saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Dry affinity b [1/Pa]", 1.0e-5},
              {"waterEnhancement", "Water enhancement alpha [kg/mol]", 0.0},
              {"waterBlocking", "Water blocking omega [kg/mol]", 0.0}};
    }
    if (type == "GAB-coadsorption-Sips")
    {
      return {{"saturationLoading", "Dry saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Dry affinity b [1/Pa]", 1.0e-5},
              {"heterogeneityExponent", "Heterogeneity exponent n [-]", 1.0},
              {"waterEnhancement", "Water enhancement alpha [kg/mol]", 0.0},
              {"waterBlocking", "Water blocking omega [kg/mol]", 0.0}};
    }
    if (type == "GAB-coadsorption-Toth")
    {
      return {{"saturationLoading", "Dry saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Dry affinity b [1/Pa]", 1.0e-5},
              {"heterogeneityParameter", "Toth heterogeneity t [-]", 1.0},
              {"waterEnhancement", "Water enhancement alpha [kg/mol]", 0.0},
              {"waterBlocking", "Water blocking omega [kg/mol]", 0.0}};
    }
    if (type == "GAB-modified-Toth")
    {
      return {{"saturationLoading", "Dry saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Dry affinity b [1/Pa]", 1.0e-5},
              {"heterogeneityParameter", "Toth heterogeneity t [-]", 1.0},
              {"waterCapacityCoefficient", "Capacity coefficient gamma [kg/mol]", 0.0},
              {"waterAffinityCoefficient", "Affinity coefficient beta [kg/mol]", 0.0}};
    }
    if (type == "Exponential-humidity-Toth")
    {
      return {{"saturationLoading", "Dry saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity b [1/Pa]", 1.0e-5},
              {"heterogeneityParameter", "Toth heterogeneity t [-]", 1.0},
              {"humidityCapacityCoefficient", "Humidity coefficient a [-]", 0.0}};
    }
    if (type == "Anti-Langmuir")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5}};
    }
    if (type == "BET")
    {
      return {{"monolayerLoading", "Monolayer loading q_m [mol/kg]", 1.0},
              {"cParameter", "BET C parameter [-]", 10.0},
              {"saturationPressureFactor", "Saturation pressure factor [-]", 1.0e-5}};
    }
    if (type == "Henry")
    {
      return {{"henryConstant", "Henry constant K_H [mol/kg/Pa]", 1.0e-5}};
    }
    if (type == "Freundlich")
    {
      return {{"freundlichConstant", "Freundlich constant K_F", 1.0e-5},
              {"heterogeneityExponent", "Heterogeneity exponent n [-]", 1.0}};
    }
    if (type == "Sips")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5},
              {"heterogeneityExponent", "Heterogeneity exponent n [-]", 1.0}};
    }
    if (type == "Langmuir-Freundlich")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5},
              {"heterogeneityExponent", "Heterogeneity exponent nu [-]", 1.0}};
    }
    if (type == "Redlich-Peterson")
    {
      return {{"linearConstant", "Linear constant K_R [mol/kg/Pa]", 1.0e-5},
              {"affinityConstant", "Affinity constant a_R", 1.0e-5},
              {"exponent", "Exponent g [-]", 1.0}};
    }
    if (type == "Toth")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5},
              {"heterogeneityParameter", "Toth heterogeneity t [-]", 1.0}};
    }
    if (type == "Unilan")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5},
              {"energyDistributionWidth", "Energy distribution width s [-]", 1.0}};
    }
    if (type == "OBrien&Myers")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5},
              {"variance", "Energy variance [-]", 1.0}};
    }
    if (type == "Quadratic")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"linearAffinity", "Linear affinity b [1/Pa]", 1.0e-5},
              {"quadraticAffinity", "Quadratic affinity c [1/Pa^2]", 1.0e-10}};
    }
    if (type == "Temkin")
    {
      return {{"saturationLoading", "Saturation loading q_sat [mol/kg]", 1.0},
              {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5},
              {"interactionParameter", "Interaction parameter theta [-]", 0.1}};
    }
    return {{"qmax", "Maximum loading qmax [mol/kg]", 1.0}, {"affinityConstant", "Affinity constant b [1/Pa]", 1.0e-5}};
  }();

  // Which parameters may depend on temperature comes from the table the input reader uses.
  static const QList<QPair<QString, Isotherm::Type>> isothermTypes{
      {"Langmuir", Isotherm::Type::Langmuir},
      {"pH-Langmuir", Isotherm::Type::Langmuir_pH},
      {"Henry", Isotherm::Type::Henry},
      {"Freundlich", Isotherm::Type::Freundlich},
      {"Sips", Isotherm::Type::Sips},
      {"Langmuir-Freundlich", Isotherm::Type::Langmuir_Freundlich},
      {"Redlich-Peterson", Isotherm::Type::Redlich_Peterson},
      {"Toth", Isotherm::Type::Toth},
      {"GAB", Isotherm::Type::GAB},
      {"Humidity-enhanced-Langmuir", Isotherm::Type::Langmuir_humidity},
      {"Moisture-swing-Langmuir", Isotherm::Type::Langmuir_moisture_swing},
      {"GAB-coadsorption-Langmuir", Isotherm::Type::GAB_Langmuir},
      {"GAB-coadsorption-Sips", Isotherm::Type::GAB_Sips},
      {"GAB-coadsorption-Toth", Isotherm::Type::GAB_Toth},
      {"GAB-modified-Toth", Isotherm::Type::GAB_modified_Toth},
      {"Exponential-humidity-Toth", Isotherm::Type::Toth_exponential_humidity}};
  for (const auto& [name, isothermType] : isothermTypes)
  {
    if (name.compare(type, Qt::CaseInsensitive) != 0) continue;
    const std::span<const TemperatureParameterSpec> laws = Isotherm::temperatureParameters(isothermType);
    for (qsizetype i = 0; i < definitions.size() && static_cast<size_t>(i) < laws.size(); ++i)
    {
      definitions[i].temperatureLaws = laws[static_cast<size_t>(i)].allowedLaws;
      definitions[i].heatLinked = laws[static_cast<size_t>(i)].heatLinked;
    }
  }
  return definitions;
}

QList<ParameterDefinition> chemDefs(const QString& type)
{
  constexpr unsigned arrhenius = TemperatureLaw::mask(TemperatureLaw::Type::Arrhenius);
  if (type == "PseudoNth")
  {
    return {{"rateCoefficient", "Rate coefficient k", 0.04, false, false, arrhenius},
            {"order", "Reaction order n [-]", 1.0},
            {"maximumLoading", "Maximum chemisorption loading q_max [mol/kg]", 0.8},
            {"heatOfChemisorption", "Heat of chemisorption Delta H [J/mol]", 42000.0}};
  }
  if (type == "Avrami")
  {
    return {{"rateCoefficient", "Rate coefficient k", 0.04, false, false, arrhenius},
            {"order", "Avrami order n [-]", 1.0},
            {"maximumLoading", "Maximum chemisorption loading q_max [mol/kg]", 0.8},
            {"heatOfChemisorption", "Heat of chemisorption Delta H [J/mol]", 42000.0}};
  }
  if (type == "General")
  {
    return {{"maximumLoading", "Maximum chemisorption loading q_max [mol/kg]", 0.55},
            {"heatOfChemisorption", "Heat of chemisorption Delta H [J/mol]", 52000.0},
            {"adsorptionRateCoefficient", "Adsorption rate coefficient k_ads", 2.5e-5, false, false, arrhenius},
            {"desorptionRateCoefficient", "Desorption rate coefficient k_des", 1.0e-6, false, false, arrhenius},
            {"poreConcentrationOrder", "Pore concentration order [-]", 1.0, true},
            {"capacityOrder", "Capacity order [-]", 1.0, true},
            {"desorptionOrder", "Desorption order [-]", 1.0, true},
            {"filmMassTransferCoefficient", "Film mass transfer coefficient [m/s]", 1.5e-3},
            {"poreDiffusivity", "Pore diffusivity [m2/s]", 8.0e-11},
            {"usePoreSurfaceTransport", "Use pore surface transport", 1.0, false, true}};
  }
  if (type == "Elovich")
  {
    return {{"maximumLoading", "Maximum chemisorption loading q_max [mol/kg]", 0.25},
            {"heatOfChemisorption", "Heat of chemisorption Delta H [J/mol]", 46000.0},
            {"alpha", "Elovich alpha", 0.02, false, false, arrhenius},
            {"beta", "Elovich beta", 3.0},
            {"filmMassTransferCoefficient", "Film mass transfer coefficient [m/s]", 1.0e-3},
            {"poreDiffusivity", "Pore diffusivity [m2/s]", 5.0e-11},
            {"usePoreSurfaceTransport", "Use pore surface transport", 1.0, false, true}};
  }
  return {{"rateCoefficient", "Rate coefficient k [1/s]", 0.04, false, false, arrhenius},
          {"maximumLoading", "Maximum chemisorption loading q_max [mol/kg]", 0.8},
          {"heatOfChemisorption", "Heat of chemisorption Delta H [J/mol]", 42000.0}};
}

void addOptionalNumber(QJsonObject& object, const QString& key, const std::optional<double>& value)
{
  if (value)
  {
    object.insert(key, *value);
  }
}

bool hasNumber(const QJsonObject& object, const QString& key)
{
  const QJsonValue value = object.value(key);
  return value.isDouble();
}

double readNumber(const QJsonObject& object, const QString& key, double fallback)
{
  const QJsonValue value = object.value(key);
  return value.isDouble() ? value.toDouble() : fallback;
}

int readInt(const QJsonObject& object, const QString& key, int fallback)
{
  const QJsonValue value = object.value(key);
  return value.isDouble() ? std::max(0, value.toInt()) : fallback;
}

bool readBool(const QJsonObject& object, const QString& key, bool fallback)
{
  const QJsonValue value = object.value(key);
  return value.isBool() ? value.toBool() : fallback;
}

QString readString(const QJsonObject& object, const QString& key, const QString& fallback)
{
  const QJsonValue value = object.value(key);
  return value.isString() ? value.toString() : fallback;
}

std::optional<double> readOptionalNumber(const QJsonObject& object, const QString& key)
{
  if (!hasNumber(object, key))
  {
    return std::nullopt;
  }
  return object.value(key).toDouble();
}

QVariantMap parametersFromArray(const QJsonArray& array, const QList<ParameterDefinition>& definitions)
{
  QVariantMap result;
  for (qsizetype i = 0; i < definitions.size(); ++i)
  {
    const ParameterDefinition& definition = definitions[i];
    if (i < array.size() && array[i].isDouble())
    {
      result.insert(definition.key, definition.integer ? static_cast<int>(array[i].toInt()) : array[i].toDouble());
    }
    else if (i < array.size() && array[i].isObject() && definition.temperatureLaws != 0)
    {
      result.insert(definition.key, array[i].toObject().toVariantMap());
    }
  }
  return defaultsFor(definitions, result);
}

QVariantMap parametersFromObject(const QJsonObject& object, const QList<ParameterDefinition>& definitions)
{
  QVariantMap result;
  for (const ParameterDefinition& definition : definitions)
  {
    const QJsonValue value = object.value(definition.key);
    if (definition.boolean && value.isBool())
    {
      result.insert(definition.key, value.toBool());
    }
    else if (value.isDouble())
    {
      result.insert(definition.key, definition.integer ? static_cast<int>(value.toInt()) : value.toDouble());
    }
    else if (value.isObject() && definition.temperatureLaws != 0)
    {
      result.insert(definition.key, value.toObject().toVariantMap());
    }
  }
  return defaultsFor(definitions, result);
}

QJsonObject parametersToObject(const QVariantMap& parameters, const QList<ParameterDefinition>& definitions)
{
  QJsonObject object;
  for (const ParameterDefinition& definition : definitions)
  {
    object.insert(definition.key, parameterJson(parameters, definition));
  }
  return object;
}

PhysisorptionSiteConfig phySiteFromObject(const QJsonObject& object, const QString& fallbackId)
{
  PhysisorptionSiteConfig site;
  site.id = readString(object, "Id", fallbackId);
  site.type = readString(object, "Type", "Langmuir");
  for (const auto& option : physisorptionTypeOptions())
    if (site.type.compare(option.second, Qt::CaseInsensitive) == 0) site.type = option.second;
  site.settings = object.value("Settings").toObject();
  for (const QString& key : {QString("RateEquation"), QString("RateEquationParameters"), QString("HeatOfAdsorption")})
    if (object.contains(key)) site.settings.insert(key, object.value(key));
  const QJsonValue parameters = object.value("Parameters");
  if (parameters.isArray())
  {
    site.parameters = parametersFromArray(parameters.toArray(), phyDefs(site.type));
  }
  else if (parameters.isObject())
  {
    site.parameters = parametersFromObject(parameters.toObject(), phyDefs(site.type));
  }
  else
  {
    site.parameters = defaultsFor(phyDefs(site.type));
  }
  if (site.type == "GAB")
  {
    // A workspace saved before temperature laws keeps the GAB energies as dHC and dHK; they are now the VantHoff laws
    // of C and K. (An imported simulation.json has passed the input reader, which rejects the old form.)
    const QJsonObject named = parameters.toObject();
    for (const auto& [key, energy] :
         {qMakePair(QString("c0"), named.value("dHC")), qMakePair(QString("k0"), named.value("dHK"))})
    {
      if (energy.isDouble() && energy.toDouble() != 0.0 && site.parameters.value(key).typeId() != QMetaType::QVariantMap)
      {
        site.parameters.insert(
            key, QVariantMap{{"Value", site.parameters.value(key).toDouble()}, {"Law", "VantHoff"}, {"Energy", energy.toDouble()}});
      }
    }
  }
  return site;
}

ChemisorptionSiteConfig chemSiteFromObject(const QJsonObject& object, const QString& fallbackId,
                                           const QString& fallbackIsothermId, double massTransfer = 0.0,
                                           double heatOfAdsorption = 0.0)
{
  ChemisorptionSiteConfig site;
  site.id = readString(object, "Id", fallbackId);
  if (object.value("Parameters").isArray())
  {
    site.type = readString(object, "RateEquation", "FirstOrder");
    for (const auto& option : chemisorptionTypeOptions())
      if (site.type.compare(option.second, Qt::CaseInsensitive) == 0) site.type = option.second;
    site.parameters = defaultsFor(chemDefs(site.type));
    if (site.type == "FirstOrder") site.parameters.insert("rateCoefficient", massTransfer);
    const auto kinetics = object.value("RateEquationParameters").toObject();
    for (const auto& definition : chemDefs(site.type))
    {
      QString key = definition.key;
      if (key == "rateCoefficient" && site.type == "FirstOrder") key = "MassTransferCoefficient";
      else key[0] = key[0].toUpper();
      if (kinetics.contains(key)) site.parameters.insert(definition.key, kinetics.value(key).toVariant());
    }

    site.parameters.insert("maximumLoading", readNumber(object, "MaximumLoading", 0.8));
    // An absent site heat inherits the component heat at import.
    site.parameters.insert("heatOfChemisorption", readNumber(object, "HeatOfAdsorption", heatOfAdsorption));
    site.isotherm = phySiteFromObject(object, fallbackIsothermId);
    site.isotherm.settings = {};
    return site;
  }
  site.type = readString(object, "Type", "FirstOrder");
  const QJsonObject parameters = object.value("Parameters").toObject();
  site.parameters = parametersFromObject(parameters, chemDefs(site.type));
  // A workspace saved before temperature laws keeps the General activation energies apart.
  for (const QString& direction : {QString("adsorption"), QString("desorption")})
  {
    const double energy = readNumber(parameters, direction + "ActivationEnergy", 0.0);
    const QString key = direction + "RateCoefficient";
    if (energy != 0.0 && site.parameters.value(key).typeId() != QMetaType::QVariantMap)
      site.parameters.insert(key, QVariantMap{{"Value", site.parameters.value(key).toDouble()}, {"Law", "Arrhenius"},
                                              {"ActivationEnergy", energy}});
  }

  const QJsonValue isothermValue =
      parameters.contains("Isotherm") ? parameters.value("Isotherm") : object.value("Isotherm");
  site.isotherm = isothermValue.isObject()
                      ? phySiteFromObject(isothermValue.toObject(), fallbackIsothermId)
                      : PhysisorptionSiteConfig{
                            fallbackIsothermId, "Langmuir",
                            defaultsFor(phyDefs("Langmuir"), {{"qmax", site.parameters.value("maximumLoading", 0.8)},
                                                              {"affinityConstant", 2.5e-5}})};
  return site;
}

QJsonArray phySitesStateJson(const QList<PhysisorptionSiteConfig>& sites)
{
  QJsonArray array;
  for (const PhysisorptionSiteConfig& site : sites)
  {
    array.append(QJsonObject{
        {"Id", site.id}, {"Type", site.type}, {"Parameters", parametersToObject(site.parameters, phyDefs(site.type))},
        {"Settings", site.settings}});
  }
  return array;
}

QJsonArray chemSitesStateJson(const QList<ChemisorptionSiteConfig>& sites)
{
  QJsonArray array;
  for (const ChemisorptionSiteConfig& site : sites)
  {
    array.append(QJsonObject{
        {"Id", site.id},
        {"Type", site.type},
        {"Parameters", parametersToObject(site.parameters, chemDefs(site.type))},
        {"Isotherm",
         QJsonObject{{"Id", site.isotherm.id},
                     {"Type", site.isotherm.type},
                     {"Parameters", parametersToObject(site.isotherm.parameters, phyDefs(site.isotherm.type))}}}});
  }
  return array;
}

QJsonArray swingPhasesStateJson(const QList<SwingPhaseConfig>& phases)
{
  QJsonArray array;
  for (const SwingPhaseConfig& phase : phases)
  {
    QJsonObject object{{"Id", phase.id},
                       {"Name", phase.name},
                       {"Temperature", phase.temperature},
                       {"InletPressure", phase.inletPressure},
                       {"Duration", phase.duration}};
    if (!phase.flowDirection.isEmpty()) object.insert("FlowDirection", phase.flowDirection);
    for (auto key = phase.otherKeys.begin(); key != phase.otherKeys.end(); ++key) object.insert(key.key(), key.value());
    array.append(object);
  }
  return array;
}

SwingPhaseConfig swingPhaseFromObject(const QJsonObject& object, const QString& fallbackId, double fallbackTemperature,
                                      double fallbackInletPressure, double fallbackDuration)
{
  SwingPhaseConfig phase;
  phase.id = readString(object, "Id", fallbackId);
  phase.name = readString(object, "Name", "Phase");
  phase.temperature = readNumber(object, "Temperature", fallbackTemperature);
  phase.inletPressure = readNumber(object, "InletPressure", fallbackInletPressure);
  const double duration = readNumber(object, "Duration", fallbackDuration);
  phase.duration = duration > 0.0 ? duration : fallbackDuration;
  phase.flowDirection = readString(object, "FlowDirection", "");
  phase.otherKeys = object;
  for (const QString& key : {"Id", "Name", "Temperature", "InletPressure", "Duration", "FlowDirection"})
  {
    phase.otherKeys.remove(key);
  }
  return phase;
}

QStringList csvTokens(QString text)
{
  text = text.trimmed();
  if (text.startsWith('[') && text.endsWith(']'))
  {
    text = text.mid(1, text.size() - 2);
  }

  QStringList tokens;
  for (QString token : text.split(',', Qt::SkipEmptyParts))
  {
    token = token.trimmed();
    if ((token.startsWith('"') && token.endsWith('"')) || (token.startsWith('\'') && token.endsWith('\'')))
    {
      token = token.mid(1, token.size() - 2).trimmed();
    }
    if (!token.isEmpty())
    {
      tokens.push_back(token);
    }
  }
  return tokens;
}

QList<double> numbersFromText(const QString& text)
{
  QList<double> values;
  for (const QString& token : csvTokens(text))
  {
    bool ok = false;
    const double value = token.toDouble(&ok);
    if (ok)
    {
      values.push_back(value);
    }
  }
  return values;
}

QString componentNameForId(const QList<ComponentConfig>& components, const QString& componentId)
{
  auto it = std::find_if(components.cbegin(), components.cend(),
                         [&componentId](const ComponentConfig& component) { return component.id == componentId; });
  return it == components.cend() ? componentId : it->name;
}

QString componentIdFromToken(const QString& token, const QList<ComponentConfig>& components)
{
  const QString trimmed = token.trimmed();
  if (trimmed.isEmpty())
  {
    return {};
  }
  auto idIt = std::find_if(components.cbegin(), components.cend(),
                           [&trimmed](const ComponentConfig& component) { return component.id == trimmed; });
  if (idIt != components.cend())
  {
    return idIt->id;
  }
  auto nameIt = std::find_if(components.cbegin(), components.cend(),
                             [&trimmed](const ComponentConfig& component) { return component.name == trimmed; });
  if (nameIt != components.cend())
  {
    return nameIt->id;
  }
  bool ok = false;
  const int index = trimmed.toInt(&ok);
  if (ok && index >= 0 && index < static_cast<int>(components.size()))
  {
    return components[index].id;
  }
  return trimmed;
}

QString componentIdFromValue(const QJsonValue& value, const QList<ComponentConfig>& components)
{
  if (value.isDouble())
  {
    const int index = value.toInt();
    if (index >= 0 && index < static_cast<int>(components.size()))
    {
      return components[index].id;
    }
    return QString::number(index);
  }
  if (value.isString())
  {
    return componentIdFromToken(value.toString(), components);
  }
  return {};
}

double stoichiometryAt(const QJsonArray& values, qsizetype index, double fallback = 1.0)
{
  return index >= 0 && index < values.size() && values[index].isDouble() ? values[index].toDouble() : fallback;
}

qsizetype participantCountFromValue(const QJsonValue& value)
{
  if (value.isArray())
  {
    return value.toArray().size();
  }
  if (value.isString())
  {
    return csvTokens(value.toString()).size();
  }
  return 0;
}

ReactionParticipantConfig reactionParticipantFromObject(const QJsonObject& object,
                                                        const QList<ComponentConfig>& components,
                                                        double fallbackStoichiometry)
{
  QString token =
      readString(object, "ComponentId", readString(object, "Component", readString(object, "Name", QString{})));
  if (token.isEmpty() && object.value("Index").isDouble())
  {
    token = QString::number(object.value("Index").toInt());
  }
  return ReactionParticipantConfig{componentIdFromToken(token, components),
                                   readNumber(object, "Stoichiometry", fallbackStoichiometry),
                                   readOptionalNumber(object, "Order")};
}

QList<ReactionParticipantConfig> reactionParticipantsFromValue(const QJsonValue& value,
                                                               const QString& stoichiometryText,
                                                               const QJsonArray& stoichiometryArray,
                                                               const QList<ComponentConfig>& components,
                                                               int fallbackIndex)
{
  QList<ReactionParticipantConfig> participants;
  const QList<double> textStoichiometry = numbersFromText(stoichiometryText);
  auto fallbackStoichiometry = [&](qsizetype index)
  {
    if (index < textStoichiometry.size())
    {
      return textStoichiometry[index];
    }
    return stoichiometryAt(stoichiometryArray, index);
  };

  if (value.isArray())
  {
    const QJsonArray array = value.toArray();
    for (qsizetype i = 0; i < array.size(); ++i)
    {
      if (array[i].isObject())
      {
        participants.push_back(
            reactionParticipantFromObject(array[i].toObject(), components, fallbackStoichiometry(i)));
      }
      else
      {
        participants.push_back(
            ReactionParticipantConfig{componentIdFromValue(array[i], components), fallbackStoichiometry(i)});
      }
    }
  }
  else
  {
    qsizetype i = 0;
    for (const QString& token : csvTokens(value.toString()))
    {
      participants.push_back(
          ReactionParticipantConfig{componentIdFromToken(token, components), fallbackStoichiometry(i)});
      ++i;
    }
  }

  if (participants.isEmpty() && fallbackIndex >= 0 && fallbackIndex < static_cast<int>(components.size()))
  {
    participants.push_back(ReactionParticipantConfig{components[fallbackIndex].id, 1.0});
  }
  return participants;
}

QJsonArray reactionParticipantsStateJson(const QList<ReactionParticipantConfig>& participants,
                                         const QList<ComponentConfig>& components)
{
  QJsonArray array;
  for (const ReactionParticipantConfig& participant : participants)
  {
    array.append(QJsonObject{{"ComponentId", participant.componentId},
                             {"Component", componentNameForId(components, participant.componentId)},
                             {"Stoichiometry", participant.stoichiometry},
                             {"Order", participant.order ? QJsonValue(*participant.order) : QJsonValue()}});
  }
  return array;
}

QString participantSummary(const QList<ReactionParticipantConfig>& participants,
                           const QList<ComponentConfig>& components)
{
  QStringList names;
  for (const ReactionParticipantConfig& participant : participants)
  {
    names.push_back(componentNameForId(components, participant.componentId));
  }
  return names.isEmpty() ? QString("none") : names.join(" + ");
}

QList<QPair<QString, QString>> componentDropdownOptions(const QList<ComponentConfig>& components,
                                                        const QString& extraComponentId = {})
{
  QList<QPair<QString, QString>> options;
  bool hasExtra = extraComponentId.isEmpty();
  for (qsizetype i = 0; i < components.size(); ++i)
  {
    const ComponentConfig& component = components[i];
    options.push_back({QString("%1: %2").arg(i).arg(component.name), component.id});
    hasExtra = hasExtra || component.id == extraComponentId;
  }
  if (!hasExtra)
  {
    options.push_back({QString("Unknown: %1").arg(extraComponentId), extraComponentId});
  }
  return options;
}

ReactionConfig reactionFromObject(const QJsonObject& object, const QString& fallbackId, const QString& fallbackName,
                                  const QList<ComponentConfig>& components)
{
  ReactionConfig reaction;
  reaction.id = readString(object, "Id", fallbackId);
  reaction.name = readString(object, "Name", fallbackName);
  reaction.phase = readString(object, "Phase", reaction.phase);
  reaction.style = readString(object, "Style", reaction.style);
  if (!QStringList{"Physisorbed", "Chemisorbed", "PoreConcentration", "BulkConcentration"}.contains(reaction.phase))
  {
    reaction.phase = "Physisorbed";
  }
  if (!QStringList{"GeneralPowerLaw", "LangmuirHinshelwood", "LangmuirHinshelwoodHougenWatson"}.contains(
          reaction.style))
  {
    reaction.style = "GeneralPowerLaw";
  }
  if (reaction.style != "GeneralPowerLaw" && reaction.phase != "BulkConcentration")
  {
    reaction.phase = "PoreConcentration";
  }

  reaction.site = readInt(object, "Site", reaction.site);
  QJsonArray reactantStoichiometryArray;
  QJsonArray productStoichiometryArray;
  const QJsonValue stoichiometryValue = object.value("Stoichiometry");
  if (stoichiometryValue.isObject())
  {
    const QJsonObject stoichiometry = stoichiometryValue.toObject();
    reactantStoichiometryArray = stoichiometry.value("Reactants").toArray();
    productStoichiometryArray = stoichiometry.value("Products").toArray();
  }
  else if (stoichiometryValue.isArray())
  {
    const QJsonArray stoichiometry = stoichiometryValue.toArray();
    const qsizetype reactantCount = participantCountFromValue(object.value("Reactants"));
    for (qsizetype i = 0; i < stoichiometry.size(); ++i)
    {
      (i < reactantCount ? reactantStoichiometryArray : productStoichiometryArray).append(stoichiometry[i]);
    }
  }

  reaction.reactants =
      reactionParticipantsFromValue(object.value("Reactants"), readString(object, "ReactantStoichiometry", {}),
                                    reactantStoichiometryArray, components, 0);
  reaction.products =
      reactionParticipantsFromValue(object.value("Products"), readString(object, "ProductStoichiometry", {}),
                                    productStoichiometryArray, components, components.size() > 1 ? 1 : 0);
  for (const auto& pair : {qMakePair(QString("ForwardOrders"), &reaction.reactants),
                           qMakePair(QString("BackwardOrders"), &reaction.products)})
  {
    const auto orders = object.value(pair.first).toArray();
    for (qsizetype i = 0; i < orders.size() && i < pair.second->size(); ++i)
      (*pair.second)[i].order = orders[i].toDouble();
  }

  const QJsonObject kinetics = object.value("Kinetics").toObject();
  const QJsonObject& numeric = kinetics.isEmpty() ? object : kinetics;
  // Each coefficient is a number or a temperature law, in the Kinetics of simulation.json or at the top level of a
  // saved workspace (keys are case-insensitive in simulation.json).
  for (const auto& [key, target] : {qMakePair(QString("ForwardRateCoefficient"), &reaction.forwardRateCoefficient),
                                    qMakePair(QString("EquilibriumConstant"), &reaction.equilibriumConstant)})
  {
    for (const QString& candidate : {key, key.at(0).toLower() + key.mid(1)})
    {
      const QJsonValue value = numeric.value(candidate);
      if (value.isDouble() || value.isObject()) *target = value.toVariant();
    }
  }
  if (hasNumber(numeric, "rateCoefficient")) reaction.forwardRateCoefficient = numeric.value("rateCoefficient").toDouble();
  reaction.heatOfReaction = readOptionalNumber(numeric, "HeatOfReaction");
  // LH/LHHW adsorption constants: keyed by component id in a saved workspace and by (case-insensitive) component
  // name in simulation.json.
  const QJsonObject adsorptionConstants = numeric.value("AdsorptionConstants").toObject();
  for (auto it = adsorptionConstants.begin(); it != adsorptionConstants.end(); ++it)
  {
    const auto component =
        std::find_if(components.cbegin(), components.cend(), [&it](const ComponentConfig& item)
                     { return item.id == it.key() || item.name.compare(it.key(), Qt::CaseInsensitive) == 0; });
    reaction.adsorptionConstants.insert(component == components.cend() ? it.key() : component->id,
                                        it.value().toVariant());
  }
  // A workspace saved before temperature laws keeps the activation energy apart, and GibbsFreeEnergy dG set K(T) and
  // the heat: they are an Arrhenius law of k+, and a VantHoff K that takes HeatOfReaction = -dG.
  const double activationEnergy = readNumber(numeric, "ForwardActivationEnergy", 0.0);
  if (activationEnergy != 0.0 && reaction.forwardRateCoefficient.typeId() != QMetaType::QVariantMap)
  {
    reaction.forwardRateCoefficient = QVariantMap{
        {"Value", reaction.forwardRateCoefficient.toDouble()}, {"Law", "Arrhenius"}, {"ActivationEnergy", activationEnergy}};
  }
  const double gibbsFreeEnergy = readNumber(numeric, "GibbsFreeEnergy", 0.0);
  if (gibbsFreeEnergy != 0.0 && !reaction.heatOfReaction)
  {
    reaction.heatOfReaction = -gibbsFreeEnergy;
    if (reaction.equilibriumConstant.typeId() != QMetaType::QVariantMap)
    {
      reaction.equilibriumConstant = QVariantMap{{"Value", reaction.equilibriumConstant.toDouble()}, {"Law", "VantHoff"}};
    }
  }
  reaction.rateLimitTime = readNumber(object, "RateLimitTime", reaction.rateLimitTime);
  return reaction;
}

QJsonArray reactionsStateJson(const QList<ReactionConfig>& reactions, const QList<ComponentConfig>& components)
{
  QJsonArray array;
  for (const ReactionConfig& reaction : reactions)
  {
    QJsonObject json{{"Id", reaction.id},
                     {"Name", reaction.name},
                     {"Phase", reaction.phase},
                     {"Style", reaction.style},
                     {"Site", reaction.site},
                     {"Reactants", reactionParticipantsStateJson(reaction.reactants, components)},
                     {"Products", reactionParticipantsStateJson(reaction.products, components)},
                     {"ForwardRateCoefficient", QJsonValue::fromVariant(reaction.forwardRateCoefficient)},
                     {"EquilibriumConstant", QJsonValue::fromVariant(reaction.equilibriumConstant)},
                     {"RateLimitTime", reaction.rateLimitTime}};
    addOptionalNumber(json, "HeatOfReaction", reaction.heatOfReaction);
    if (!reaction.adsorptionConstants.isEmpty())
    {
      json.insert("AdsorptionConstants", QJsonObject::fromVariantMap(reaction.adsorptionConstants));
    }
    array.append(json);
  }
  return array;
}

QJsonObject runRecordStateJson(const RupturaRunner::RunRecord& record, const QString& logText)
{
  if (record.id.isEmpty())
  {
    return {};
  }
  QString cappedLog = logText;
  if (cappedLog.size() > savedLogLimit)
  {
    cappedLog = cappedLog.right(savedLogLimit);
  }
  return QJsonObject{{"Id", record.id},
                     {"Directory", record.directory},
                     {"SimulationPath", record.simulationPath},
                     {"NotebookPath", record.notebookPath},
                     {"LogPath", record.logPath},
                     {"Status", record.status},
                     {"ErrorMessage", record.errorMessage},
                     {"StartedAt", isoDate(record.startedAt)},
                     {"FinishedAt", isoDate(record.finishedAt)},
                     {"ExitCode", record.exitCode},
                     {"LogTail", cappedLog}};
}

RupturaRunner::RunRecord runRecordFromObject(const QJsonObject& object)
{
  RupturaRunner::RunRecord record;
  record.id = readString(object, "Id", {});
  record.directory = readString(object, "Directory", {});
  record.simulationPath = readString(object, "SimulationPath", {});
  record.notebookPath = readString(object, "NotebookPath", {});
  record.logPath = readString(object, "LogPath", {});
  record.status = readString(object, "Status", record.id.isEmpty() ? "created" : "completed");
  record.errorMessage = readString(object, "ErrorMessage", {});
  record.startedAt = dateFromIso(readString(object, "StartedAt", {}));
  record.finishedAt = dateFromIso(readString(object, "FinishedAt", {}));
  const QJsonValue exitCode = object.value("ExitCode");
  record.exitCode = exitCode.isDouble() ? exitCode.toInt() : record.exitCode;
  return record;
}

bool stringSetContains(const QList<ColumnConfig>& columns, const QString& id)
{
  return std::any_of(columns.cbegin(), columns.cend(), [&id](const ColumnConfig& column) { return column.id == id; });
}

bool stringSetContains(const QList<ComponentConfig>& components, const QString& id)
{
  return std::any_of(components.cbegin(), components.cend(),
                     [&id](const ComponentConfig& component) { return component.id == id; });
}

int numericIdSuffix(const QString& id)
{
  const qsizetype dash = id.lastIndexOf('-');
  if (dash < 0 || dash + 1 >= id.size())
  {
    return 0;
  }
  bool ok = false;
  const int value = id.mid(dash + 1).toInt(&ok);
  return ok ? value : 0;
}
}  // namespace

template <>
ColumnConfig& MainWindow::selected<ColumnConfig>()
{
  return columns_[selectedColumnIndex_];
}

template <>
SimulationConfig& MainWindow::selected<SimulationConfig>()
{
  return simulations_[selectedSimulationIndex_];
}

template <>
ComponentConfig& MainWindow::selected<ComponentConfig>()
{
  return components_[selectedComponentIndex_];
}

template <>
ReactionConfig& MainWindow::selected<ReactionConfig>()
{
  return reactions_[selectedReactionIndex_];
}

template <typename Config>
QLineEdit* MainWindow::addNumberField(QFormLayout* form, const QString& label, double Config::*member)
{
  return addNumberField(form, label, selected<Config>().*member,
                        [this, member](double value) { selected<Config>().*member = value; });
}

template <typename Config>
QLineEdit* MainWindow::addOptionalNumberField(QFormLayout* form, const QString& label,
                                              std::optional<double> Config::*member)
{
  return addOptionalNumberField(form, label, selected<Config>().*member,
                                [this, member](std::optional<double> value) { selected<Config>().*member = value; });
}

template <typename Config>
QComboBox* MainWindow::addComboField(QFormLayout* form, const QString& label, QString Config::*member,
                                     const QList<QPair<QString, QString>>& options)
{
  return addComboField(form, label, selected<Config>().*member, options,
                       [this, member](const QString& value) { selected<Config>().*member = value; });
}

MainWindow::MainWindow(QWidget* parent) : QMainWindow(parent)
{
  seedDefaults();
  setupUi();
  refreshAll();

  connect(&runner_, &RupturaRunner::runChanged, this,
          [this](const RupturaRunner::RunRecord& record)
          {
            SimulationConfig* simulation = simulationForRunId(record.id);
            if (simulation == nullptr && !pendingRunSimulationId_.isEmpty())
            {
              simulation = simulationById(pendingRunSimulationId_);
            }
            if (simulation != nullptr)
            {
              simulation->lastRun = record;
              if (simulation->lastRunLog.isEmpty())
              {
                simulation->lastRunLog = runner_.logTail(record);
              }
              refreshSimulationList();
            }
            refreshRunPanel();
          });
  connect(&runner_, &RupturaRunner::logChanged, this,
          [this](const QString& runId, const QString& logText)
          {
            if (SimulationConfig* simulation = simulationForRunId(runId))
            {
              simulation->lastRunLog = logText;
              if (activeSimulation() == simulation)
              {
                runLog_->setPlainText(logText);
                runLog_->moveCursor(QTextCursor::End);
              }
            }
          });
  connect(&runner_, &RupturaRunner::statusMessage, this,
          [this](const QString& message) { statusBar()->showMessage(message); });

  statusBar()->showMessage(QString("Simulation directory: %1").arg(QDir::current().filePath("simulations")));
}

void MainWindow::seedDefaults()
{
  ComponentConfig helium;
  helium.id = "component-1";
  helium.name = "Helium";
  helium.carrier = true;
  helium.molecularWeight = 0.004;
  components_.push_back(helium);

  ComponentConfig co2;
  co2.id = "component-2";
  co2.name = "CO2";
  co2.massTransfer = 0.06;
  co2.axialDispersion = 0.0;
  co2.molecularWeight = 0.044;
  co2.physisorptionSites.push_back(makePhySite("Langmuir", {{"qmax", 4.4}, {"affinityConstant", 0.000291}}));
  co2.physisorptionSites.push_back(makePhySite("Langmuir", {{"qmax", 10.0}, {"affinityConstant", 0.000000696}}));
  components_.push_back(co2);

  ComponentConfig propane;
  propane.id = "component-3";
  propane.name = "C3H8";
  propane.massTransfer = 0.06;
  propane.axialDispersion = 0.0;
  propane.physisorptionSites.push_back(makePhySite("Langmuir", {{"qmax", 4.97}, {"affinityConstant", 0.0000000093}}));
  propane.physisorptionSites.push_back(makePhySite("Langmuir", {{"qmax", 2.99}, {"affinityConstant", 0.000279}}));
  components_.push_back(propane);

  columns_.push_back(makeColumn("Column 1", "column-1"));

  QList<FeedConfig> feeds{{"component-1", 0.9}, {"component-2", 0.05}, {"component-3", 0.05}};
  simulations_.push_back(makeSimulation("MOR breakthrough", "simulation-1", feeds));
}

void MainWindow::setupUi()
{
  setWindowTitle("Ruptura Lab");
  resize(1500, 940);

  auto* central = new QWidget(this);
  auto* root = new QVBoxLayout(central);
  root->setContentsMargins(10, 10, 10, 10);
  root->setSpacing(10);

  auto* topBar = new QHBoxLayout;
  auto* title = new QLabel("Ruptura Lab");
  QFont titleFont = title->font();
  titleFont.setPointSize(titleFont.pointSize() + 5);
  titleFont.setBold(true);
  title->setFont(titleFont);
  runBadge_ = new QLabel("Idle");
  runBadge_->setObjectName("runBadge");
  topBar->addWidget(title);
  topBar->addWidget(runBadge_);
  loadButton_ = new QPushButton("Load JSON");
  saveButton_ = new QPushButton("Save State");
  topBar->addWidget(loadButton_);
  topBar->addWidget(saveButton_);
  topBar->addStretch(1);

  runTopButton_ = new QPushButton("Run");
  cancelTopButton_ = new QPushButton("Cancel");
  analysisTopButton_ = new QPushButton("Analysis");
  topBar->addWidget(runTopButton_);
  topBar->addWidget(cancelTopButton_);
  topBar->addWidget(analysisTopButton_);
  root->addLayout(topBar);

  auto* splitter = new QSplitter(Qt::Horizontal);
  root->addWidget(splitter, 1);

  componentList_ = new QListWidget;
  auto* componentControls = new QWidget;
  auto* componentButtons = new QHBoxLayout(componentControls);
  componentButtons->setContentsMargins(0, 0, 0, 0);
  auto* addComponentButton = new QPushButton("Add");
  auto* removeComponentButton = new QPushButton("Remove");
  componentButtons->addWidget(addComponentButton);
  componentButtons->addWidget(removeComponentButton);
  splitter->addWidget(
      makePanel("Components", new QLabel("Library"), componentControls, componentList_, &componentEditorLayout_));

  columnList_ = new QListWidget;
  activeColumnLabel_ = new QLabel("Column 1");
  auto* columnControls = new QWidget;
  auto* columnButtons = new QHBoxLayout(columnControls);
  columnButtons->setContentsMargins(0, 0, 0, 0);
  auto* addColumnButton = new QPushButton("Add");
  auto* removeColumnButton = new QPushButton("Remove");
  columnButtons->addWidget(addColumnButton);
  columnButtons->addWidget(removeColumnButton);
  splitter->addWidget(makePanel("Columns", activeColumnLabel_, columnControls, columnList_, &columnEditorLayout_));

  reactionList_ = new QListWidget;
  reactionsMetric_ = new QLabel("0");
  auto* reactionControls = new QWidget;
  auto* reactionButtons = new QHBoxLayout(reactionControls);
  reactionButtons->setContentsMargins(0, 0, 0, 0);
  auto* addReactionButton = new QPushButton("Add");
  auto* removeReactionButton = new QPushButton("Remove");
  reactionButtons->addWidget(addReactionButton);
  reactionButtons->addWidget(removeReactionButton);
  splitter->addWidget(
      makePanel("Reactions", reactionsMetric_, reactionControls, reactionList_, &reactionEditorLayout_));

  simulationList_ = new QListWidget;
  fractionTotalLabel_ = new QLabel("0.000");
  auto* simulationControls = new QWidget;
  auto* simulationButtons = new QHBoxLayout(simulationControls);
  simulationButtons->setContentsMargins(0, 0, 0, 0);
  auto* addSimulationButton = new QPushButton("Add");
  auto* copySimulationButton = new QPushButton("Copy");
  auto* removeSimulationButton = new QPushButton("Remove");
  simulationButtons->addWidget(addSimulationButton);
  simulationButtons->addWidget(copySimulationButton);
  simulationButtons->addWidget(removeSimulationButton);
  splitter->addWidget(
      makePanel("Simulations", fractionTotalLabel_, simulationControls, simulationList_, &simulationEditorLayout_));

  auto* experiment = new QFrame;
  experiment->setFrameShape(QFrame::StyledPanel);
  auto* experimentLayout = new QVBoxLayout(experiment);
  auto* experimentHead = new QHBoxLayout;
  auto* experimentTitle = new QLabel("Experiment");
  QFont headFont = experimentTitle->font();
  headFont.setBold(true);
  experimentTitle->setFont(headFont);
  runIdLabel_ = new QLabel("No run");
  experimentHead->addWidget(experimentTitle);
  experimentHead->addStretch(1);
  experimentHead->addWidget(runIdLabel_);
  experimentLayout->addLayout(experimentHead);

  auto* commandRow = new QHBoxLayout;
  runButton_ = new QPushButton("Run");
  cancelButton_ = new QPushButton("Cancel");
  analysisButton_ = new QPushButton("Analysis");
  commandRow->addWidget(runButton_);
  commandRow->addWidget(cancelButton_);
  commandRow->addWidget(analysisButton_);
  experimentLayout->addLayout(commandRow);

  experimentLayout->addWidget(new QLabel("Simulation JSON"));
  jsonPreview_ = new QPlainTextEdit;
  jsonPreview_->setReadOnly(true);
  jsonPreview_->setMinimumHeight(300);
  experimentLayout->addWidget(jsonPreview_, 2);
  experimentLayout->addWidget(new QLabel("Run log"));
  runLog_ = new QPlainTextEdit;
  runLog_->setReadOnly(true);
  runLog_->setMinimumHeight(220);
  runLog_->setPlaceholderText("No run log for the selected simulation.");
  experimentLayout->addWidget(runLog_, 1);
  splitter->addWidget(experiment);
  splitter->setSizes({330, 330, 330, 370, 420});

  connect(addComponentButton, &QPushButton::clicked, this, &MainWindow::addComponent);
  connect(removeComponentButton, &QPushButton::clicked, this, &MainWindow::removeComponent);
  connect(addColumnButton, &QPushButton::clicked, this, &MainWindow::addColumn);
  connect(removeColumnButton, &QPushButton::clicked, this, &MainWindow::removeColumn);
  connect(addReactionButton, &QPushButton::clicked, this, &MainWindow::addReaction);
  connect(removeReactionButton, &QPushButton::clicked, this, &MainWindow::removeReaction);
  connect(addSimulationButton, &QPushButton::clicked, this, &MainWindow::addSimulation);
  connect(copySimulationButton, &QPushButton::clicked, this, &MainWindow::copySimulation);
  connect(removeSimulationButton, &QPushButton::clicked, this, &MainWindow::removeSimulation);
  connect(loadButton_, &QPushButton::clicked, this, &MainWindow::loadJsonFile);
  connect(saveButton_, &QPushButton::clicked, this, &MainWindow::saveStateFile);
  connect(runButton_, &QPushButton::clicked, this, &MainWindow::runSimulation);
  connect(runTopButton_, &QPushButton::clicked, this, &MainWindow::runSimulation);
  connect(cancelButton_, &QPushButton::clicked, this, &MainWindow::cancelSimulation);
  connect(cancelTopButton_, &QPushButton::clicked, this, &MainWindow::cancelSimulation);
  connect(analysisButton_, &QPushButton::clicked, this, &MainWindow::openAnalysis);
  connect(analysisTopButton_, &QPushButton::clicked, this, &MainWindow::openAnalysis);

  connect(componentList_, &QListWidget::currentRowChanged, this,
          [this](int row)
          {
            if (rebuilding_ || row < 0) return;
            selectedComponentIndex_ = row;
            rebuildComponentEditor();
          });
  connect(columnList_, &QListWidget::currentRowChanged, this,
          [this](int row)
          {
            if (rebuilding_ || row < 0) return;
            selectedColumnIndex_ = row;
            rebuildColumnEditor();
          });
  connect(reactionList_, &QListWidget::currentRowChanged, this,
          [this](int row)
          {
            if (rebuilding_ || row < 0) return;
            selectedReactionIndex_ = row;
            rebuildReactionEditor();
          });
  connect(simulationList_, &QListWidget::currentRowChanged, this,
          [this](int row)
          {
            if (rebuilding_ || row < 0) return;
            selectedSimulationIndex_ = row;
            refreshAll();
          });

  setCentralWidget(central);
  setStyleSheet(
      "QFrame { border: 1px solid #cfd6dc; border-radius: 6px; }"
      "QFrame QLabel, QFrame QLineEdit, QFrame QComboBox, QFrame QListWidget, QFrame QTableWidget, "
      "QFrame QPlainTextEdit, QFrame QGroupBox { border: 0; }"
      "QLabel#runBadge { padding: 3px 8px; border: 1px solid #cfd6dc; border-radius: 5px; }"
      "QGroupBox { font-weight: 700; margin-top: 10px; padding-top: 12px; }"
      "QGroupBox::title { subcontrol-origin: margin; left: 8px; padding: 0 3px; }"
      "QPushButton { min-height: 28px; padding: 4px 10px; border: 1px solid #8894a3; border-radius: 5px; "
      "background: #f4f7fa; color: #1f2933; }"
      "QPushButton:hover { background: #e8eef5; border-color: #66758a; }"
      "QPushButton:pressed { background: #dce5ee; border-color: #4d5e73; padding-top: 5px; padding-bottom: 3px; }"
      "QPushButton:disabled { color: #8b949e; background: #eef1f4; border-color: #c4ccd4; }"
      "QPlainTextEdit { font-family: Menlo, Consolas, monospace; font-size: 12px; }");

  updateRunControls("idle");
}

QWidget* MainWindow::makePanel(const QString& title, QLabel* metric, QWidget* controls, QListWidget* list,
                               QVBoxLayout** editorLayout)
{
  auto* panel = new QFrame;
  panel->setFrameShape(QFrame::StyledPanel);
  auto* layout = new QVBoxLayout(panel);

  auto* head = new QHBoxLayout;
  auto* titleLabel = new QLabel(title);
  QFont font = titleLabel->font();
  font.setBold(true);
  titleLabel->setFont(font);
  head->addWidget(titleLabel);
  head->addStretch(1);
  head->addWidget(metric);
  layout->addLayout(head);
  layout->addWidget(controls);
  layout->addWidget(list, 0);

  auto* editor = new QWidget;
  *editorLayout = new QVBoxLayout(editor);
  (*editorLayout)->setContentsMargins(0, 0, 0, 0);
  (*editorLayout)->setSpacing(8);
  auto* scroll = new QScrollArea;
  scroll->setWidget(editor);
  scroll->setWidgetResizable(true);
  scroll->setFrameShape(QFrame::NoFrame);
  layout->addWidget(scroll, 1);
  return panel;
}

void MainWindow::refreshAll()
{
  rebuilding_ = true;
  refreshComponentList();
  refreshColumnList();
  refreshReactionList();
  refreshSimulationList();
  rebuilding_ = false;
  rebuildComponentEditor();
  rebuildColumnEditor();
  rebuildReactionEditor();
  rebuildSimulationEditor();
  updatePreview();
  updateFractionTotal();
  refreshRunPanel();
}

void MainWindow::refreshComponentList()
{
  QSignalBlocker blocker(componentList_);
  componentList_->clear();
  for (const ComponentConfig& component : components_)
  {
    const QString kind = component.carrier ? "carrier"
                                           : QString("%1 phys, %2 chem")
                                                 .arg(static_cast<int>(component.physisorptionSites.size()))
                                                 .arg(static_cast<int>(component.chemisorptionSites.size()));
    componentList_->addItem(QString("%1\n%2").arg(component.name, kind));
  }
  const int count = static_cast<int>(components_.size());
  selectedComponentIndex_ = std::clamp(selectedComponentIndex_, 0, std::max(0, count - 1));
  if (!components_.isEmpty())
  {
    componentList_->setCurrentRow(selectedComponentIndex_);
  }
}

void MainWindow::refreshColumnList()
{
  QSignalBlocker blocker(columnList_);
  columnList_->clear();
  const SimulationConfig* simulation = activeSimulation();
  for (const ColumnConfig& column : columns_)
  {
    const bool active = simulation != nullptr && simulation->selectedColumnId == column.id;
    columnList_->addItem(
        QString("%1%2\n%3, %4, T=%5 K")
            .arg(active ? "* " : "", column.name, column.mode, column.geometryType, numberText(column.temperature)));
  }
  const int count = static_cast<int>(columns_.size());
  selectedColumnIndex_ = std::clamp(selectedColumnIndex_, 0, std::max(0, count - 1));
  if (!columns_.isEmpty())
  {
    columnList_->setCurrentRow(selectedColumnIndex_);
  }
  const ColumnConfig* selected = simulation == nullptr ? nullptr : selectedColumnFor(*simulation);
  activeColumnLabel_->setText(selected == nullptr ? "No column" : selected->name);
}

void MainWindow::refreshReactionList()
{
  QSignalBlocker blocker(reactionList_);
  reactionList_->clear();
  for (const ReactionConfig& reaction : reactions_)
  {
    const QString detail = QString("%1 -> %2, %3, site %4")
                               .arg(participantSummary(reaction.reactants, components_),
                                    participantSummary(reaction.products, components_), reaction.phase)
                               .arg(reaction.site);
    reactionList_->addItem(QString("%1\n%2").arg(reaction.name, detail));
  }
  const int count = static_cast<int>(reactions_.size());
  selectedReactionIndex_ = std::clamp(selectedReactionIndex_, 0, std::max(0, count - 1));
  if (!reactions_.isEmpty())
  {
    reactionList_->setCurrentRow(selectedReactionIndex_);
  }
  reactionsMetric_->setText(QString::number(count));
}

void MainWindow::refreshSimulationList()
{
  QSignalBlocker blocker(simulationList_);
  simulationList_->clear();
  for (const SimulationConfig& simulation : simulations_)
  {
    QString detail = QString("%1, %2 components, sum y=%3")
                         .arg(simulation.type)
                         .arg(static_cast<int>(componentsJsonFor(simulation).size()))
                         .arg(numberText(feedTotal(simulation)));
    if (simulation.type == "SwingAdsorption")
    {
      double totalDuration = 0.0;
      for (const SwingPhaseConfig& phase : simulation.swingPhases)
      {
        totalDuration += phase.duration;
      }
      detail += QString(", %1 phases, %2 s")
                    .arg(static_cast<int>(simulation.swingPhases.size()))
                    .arg(numberText(totalDuration));
    }
    if (!simulation.lastRun.id.isEmpty())
    {
      detail += QString(", last run: %1").arg(simulation.lastRun.status);
    }
    simulationList_->addItem(QString("%1\n%2").arg(simulation.displayName, detail));
  }
  const int count = static_cast<int>(simulations_.size());
  selectedSimulationIndex_ = std::clamp(selectedSimulationIndex_, 0, std::max(0, count - 1));
  if (!simulations_.isEmpty())
  {
    simulationList_->setCurrentRow(selectedSimulationIndex_);
  }
}

void MainWindow::rebuildComponentEditor()
{
  clearLayout(componentEditorLayout_);
  if (components_.isEmpty())
  {
    componentEditorLayout_->addWidget(new QLabel("No component selected."));
    return;
  }

  ComponentConfig& component = components_[selectedComponentIndex_];
  auto* form = new QFormLayout;
  componentEditorLayout_->addLayout(form);
  addTextField(form, "Name", component.name,
               [this](const QString& value)
               {
                 components_[selectedComponentIndex_].name = value;
                 refreshComponentList();
                 refreshSimulationList();
               });
  addCheckField(form, "Carrier gas", component.carrier,
                [this](bool value)
                {
                  components_[selectedComponentIndex_].carrier = value;
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                });
  addCheckField(form, "Water (sets the relative humidity)", component.water,
                [this](bool value) { components_[selectedComponentIndex_].water = value; });
  addOptionalNumberField(form, "Mass transfer coefficient [1/s]", &ComponentConfig::massTransfer);
  addOptionalNumberField(form, "Axial dispersion coefficient [m2/s]", &ComponentConfig::axialDispersion);
  addOptionalNumberField(form, "Molecular weight [kg/mol]", &ComponentConfig::molecularWeight);
  addOptionalNumberField(form, "Heat of adsorption [J/mol]", &ComponentConfig::heatOfAdsorption);
  addOptionalNumberField(form, "Reference temperature [K]", &ComponentConfig::referenceTemperature);
  addSettingNumbers(form, component.settings,
                    {{"LiquidPhaseConcentration", "Liquid inlet concentration [mol/m3]", 0.0},
                     {"InitialLiquidPhaseConcentration", "Initial liquid concentration [mol/m3]", 0.0}});

  auto* physGroup = new QGroupBox("Physisorption");
  auto* physLayout = new QVBoxLayout(physGroup);
  for (qsizetype siteIndex = 0; siteIndex < component.physisorptionSites.size(); ++siteIndex)
  {
    PhysisorptionSiteConfig& site = component.physisorptionSites[siteIndex];
    auto* siteGroup = new QGroupBox(QString("Site %1").arg(siteIndex + 1));
    auto* siteLayout = new QVBoxLayout(siteGroup);
    auto* siteForm = new QFormLayout;
    siteLayout->addLayout(siteForm);
    addComboField(siteForm, "Type", site.type, physisorptionTypeOptions(),
                  [this, siteIndex](const QString& value)
                  {
                    PhysisorptionSiteConfig& current =
                        components_[selectedComponentIndex_].physisorptionSites[siteIndex];
                    current.type = value;
                    current.parameters = defaultsFor(phyDefs(value), current.parameters);
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  });
    addParameterFields(siteForm, site.parameters, phyDefs(site.type));
    addSettingNumbers(siteForm, site.settings,
                      {{"HeatOfAdsorption", "Site heat of adsorption [J/mol]", component.heatOfAdsorption.value_or(0.0)}});
    const auto rateParameters = site.settings.value("RateEquationParameters").toObject();
    addOptionalNumberField(siteForm, "Site mass transfer [1/s]", readOptionalNumber(rateParameters, "MassTransferCoefficient"),
                           [this, siteIndex](std::optional<double> value)
                           {
                             auto& settings = components_[selectedComponentIndex_].physisorptionSites[siteIndex].settings;
                             if (value) settings.insert("RateEquationParameters", QJsonObject{{"MassTransferCoefficient", *value}});
                             else settings.remove("RateEquationParameters");
                           });
    auto* remove = new QPushButton("Remove site");
    connect(remove, &QPushButton::clicked, this,
            [this, siteIndex]()
            {
              components_[selectedComponentIndex_].physisorptionSites.removeAt(siteIndex);
              refreshAll();
            });
    siteLayout->addWidget(remove);
    physLayout->addWidget(siteGroup);
  }
  auto* addPhys = new QPushButton("Add physisorption site");
  connect(addPhys, &QPushButton::clicked, this,
          [this]()
          {
            components_[selectedComponentIndex_].physisorptionSites.push_back(makePhySite());
            refreshAll();
          });
  physLayout->addWidget(addPhys);
  componentEditorLayout_->addWidget(physGroup);

  auto* chemGroup = new QGroupBox("Chemisorption");
  auto* chemLayout = new QVBoxLayout(chemGroup);
  for (qsizetype siteIndex = 0; siteIndex < component.chemisorptionSites.size(); ++siteIndex)
  {
    ChemisorptionSiteConfig& site = component.chemisorptionSites[siteIndex];
    auto* siteGroup = new QGroupBox(QString("Site %1").arg(siteIndex + 1));
    auto* siteLayout = new QVBoxLayout(siteGroup);
    auto* siteForm = new QFormLayout;
    siteLayout->addLayout(siteForm);
    addComboField(siteForm, "Type", site.type, chemisorptionTypeOptions(),
                  [this, siteIndex](const QString& value)
                  {
                    ChemisorptionSiteConfig& current =
                        components_[selectedComponentIndex_].chemisorptionSites[siteIndex];
                    current.type = value;
                    current.parameters = defaultsFor(chemDefs(value), current.parameters);
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  });
    addParameterFields(siteForm, site.parameters, chemDefs(site.type));

    auto* isoGroup = new QGroupBox("Equilibrium isotherm");
    auto* isoForm = new QFormLayout(isoGroup);
    addComboField(isoForm, "Type", site.isotherm.type, physisorptionTypeOptions(),
                  [this, siteIndex](const QString& value)
                  {
                    PhysisorptionSiteConfig& isotherm =
                        components_[selectedComponentIndex_].chemisorptionSites[siteIndex].isotherm;
                    isotherm.type = value;
                    isotherm.parameters = defaultsFor(phyDefs(value), isotherm.parameters);
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  });
    addParameterFields(isoForm, site.isotherm.parameters, phyDefs(site.isotherm.type));
    siteLayout->addWidget(isoGroup);

    auto* remove = new QPushButton("Remove site");
    connect(remove, &QPushButton::clicked, this,
            [this, siteIndex]()
            {
              components_[selectedComponentIndex_].chemisorptionSites.removeAt(siteIndex);
              refreshAll();
            });
    siteLayout->addWidget(remove);
    chemLayout->addWidget(siteGroup);
  }
  auto* addChem = new QPushButton("Add chemisorption site");
  connect(addChem, &QPushButton::clicked, this,
          [this]()
          {
            components_[selectedComponentIndex_].chemisorptionSites.push_back(makeChemSite());
            refreshAll();
          });
  chemLayout->addWidget(addChem);
  componentEditorLayout_->addWidget(chemGroup);
  componentEditorLayout_->addStretch(1);
}

void MainWindow::rebuildColumnEditor()
{
  clearLayout(columnEditorLayout_);
  if (columns_.isEmpty())
  {
    columnEditorLayout_->addWidget(new QLabel("No column selected."));
    return;
  }

  ColumnConfig& column = columns_[selectedColumnIndex_];
  auto* form = new QFormLayout;
  columnEditorLayout_->addLayout(form);
  addSettingCombo(form, "Fluid phase", column.settings, "FluidPhase", "Gas", {"Gas", "Liquid"});
  if (column.settings.value("FluidPhase").toString().compare("Liquid", Qt::CaseInsensitive) == 0)
  {
    addSettingNumbers(form, column.settings, {{"LiquidDensity", "Liquid density [kg/m3]", 1000.0},
                                             {"pHValue", "Fixed pH [-]", 7.0}, {"pKw", "pKw [-]", 14.0}});
    addSettingCombo(form, "pH mode", column.settings, "pHMode", "Fixed", {"Fixed", "HPlus", "OHMinus"});
  }
  else
  {
    // The axial dispersion coefficients are then given at this state and scale as (T/T_ref)^1.75 (P_ref/P).
    addSettingNumbers(form, column.settings,
                      {{"AxialDispersionReferenceTemperature", "Dispersion reference temperature [K] (0 = none)", 0.0},
                       {"AxialDispersionReferencePressure", "Dispersion reference pressure [Pa] (0 = none)", 0.0}});
    // The humidity-dependent isotherms take RH from the Water component unless RelativeHumidity fixes it.
    const QStringList humiditySources{"Water component", "Fixed"};
    addComboField(form, "Relative humidity",
                  QString(column.settings.contains("RelativeHumidity") ? "Fixed" : "Water component"),
                  valueOptions(humiditySources),
                  [this, &column](const QString& value)
                  {
                    if (value == "Fixed")
                      column.settings.insert("RelativeHumidity", column.settings.value("RelativeHumidity").toDouble(0.5));
                    else
                      column.settings.remove("RelativeHumidity");
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  });
    if (column.settings.contains("RelativeHumidity"))
    {
      addSettingNumbers(form, column.settings, {{"RelativeHumidity", "Fixed relative humidity [-]", 0.5}});
    }
    const QStringList saturationModels{"Antoine", "Wagner", "Constant"};
    addComboField(form, "Water saturation pressure",
                  choice(column.settings.value("WaterSaturationPressureModel").toString("Antoine"), saturationModels),
                  valueOptions(saturationModels),
                  [this, &column](const QString& value)
                  {
                    column.settings.insert("WaterSaturationPressureModel", value);
                    if (value != "Constant") column.settings.remove("WaterSaturationPressure");
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  });
    if (column.settings.value("WaterSaturationPressureModel").toString() == "Constant")
    {
      addSettingNumbers(form, column.settings, {{"WaterSaturationPressure", "Water saturation pressure [Pa]", 3169.0}});
    }
  }
  addTextField(form, "Name", column.name,
               [this](const QString& value)
               {
                 columns_[selectedColumnIndex_].name = value;
                 refreshColumnList();
                 refreshSimulationList();
               });
  addComboField(form, "Column mode", column.mode,
                {{"Single packed bed", "single"},
                 {"Multi-bed column", "multibed"},
                 {"Uniform adsorbent mix", "mixed"}},
                [this](const QString& value)
                {
                  columns_[selectedColumnIndex_].mode = value;
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                });
  addComboField(form, "Boundary condition", &ColumnConfig::boundaryCondition,
                {{"Inlet pressure + velocity", "InletPressureInletVelocity"},
                 {"Inlet + outlet pressure", "InletPressureOutletPressure"},
                 {"Velocity + outlet pressure", "InletVelocityOutletPressure"},
                 {"Fixed velocity", "FixedVelocity"},
                 {"Fixed pressure + velocity", "FixedPressureInletVelocity"}});
  addNumberField(form, "Temperature [K]", &ColumnConfig::temperature);
  addNumberField(form, "Dynamic viscosity [Pa s]", &ColumnConfig::dynamicViscosity);
  addNumberField(form, "Void fraction [-]", &ColumnConfig::voidFraction);
  addNumberField(form, "Particle density [kg/m3]", &ColumnConfig::particleDensity);
  addNumberField(form, "Particle diameter [m]", &ColumnConfig::particleDiameter);
  if (column.mode != "multibed")
  {
    addNumberField(form, "Column length [m]", &ColumnConfig::length);
  }
  QLabel* totalLengthLabel = nullptr;
  if (column.mode == "multibed")
  {
    double totalLength = 0.0;
    for (qsizetype i = 0; i < column.beds.size(); ++i)
    {
      totalLength += column.beds[i].length;
      if (i + 1 < column.beds.size()) totalLength += column.beds[i].interfaceLengthAfter;
    }
    totalLengthLabel = new QLabel(numberText(totalLength));
    form->addRow("Total column length [m]", totalLengthLabel);
  }
  addNumberField(form, "Inlet pressure [Pa]", &ColumnConfig::inletPressure);
  addOptionalNumberField(form, "Outlet pressure [Pa]", &ColumnConfig::outletPressure);
  addNumberField(form, "Pressure gradient [Pa/m]", &ColumnConfig::pressureGradient);
  addNumberField(form, "Column entrance velocity [m/s]", &ColumnConfig::velocity);
  addCheckField(form, "Energy balance", column.energyBalance,
                [this](bool value)
                {
                  columns_[selectedColumnIndex_].energyBalance = value;
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                });

  auto* geometry = new QGroupBox("Geometry");
  auto* geometryForm = new QFormLayout(geometry);
  addComboField(geometryForm, "Type", column.geometryType, {{"Packed bed", "PackedBed"}, {"Monolith", "Monolith"}},
                [this](const QString& value)
                {
                  columns_[selectedColumnIndex_].geometryType = value;
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                });
  if (column.geometryType == "Monolith")
  {
    addComboField(
        geometryForm, "Channel shape", &ColumnConfig::channelShape,
        {{"Triangular", "triangular"}, {"Square", "square"}, {"Hexagonal", "hexagonal"}, {"Circular", "circular"}});
    addNumberField(geometryForm, "Internal channel dimension [m]", &ColumnConfig::internalChannelDimension);
    addNumberField(geometryForm, "Outer diameter [m]", &ColumnConfig::outerDiameter);
    addNumberField(geometryForm, "Number of channels [-]", static_cast<double>(column.numberOfChannels),
                   [this](double value)
                   { columns_[selectedColumnIndex_].numberOfChannels = std::max(1, static_cast<int>(value)); });
    addNumberField(geometryForm, "Washcoat thickness [m]", &ColumnConfig::washcoatThickness);
    addOptionalNumberField(geometryForm, "Washcoat volume / channel volume [-]",
                           &ColumnConfig::washcoatVolumePerChannelVolume);
  }
  else
  {
    addNumberField(geometryForm, "Internal diameter [m]", &ColumnConfig::internalDiameter);
    addNumberField(geometryForm, "Outer diameter [m]", &ColumnConfig::outerDiameter);
  }
  addSettingNumbers(geometryForm, column.geometrySettings, {{"AxialDispersivity", "Axial dispersivity [m]", 0.0}});
  if (column.geometryType == "Monolith")
  {
    addSettingNumbers(geometryForm, column.geometrySettings,
                      {{"ForchheimerCoefficient", "Forchheimer coefficient [1/m]", 0.0}});
    addCheckField(geometryForm, "Taylor-Aris dispersion", column.geometrySettings.value("TaylorArisDispersion").toBool(),
                  [this](bool value) { columns_[selectedColumnIndex_].geometrySettings.insert("TaylorArisDispersion", value); });
  }
  columnEditorLayout_->addWidget(geometry);

  if (column.energyBalance)
  {
    auto* thermal = new QGroupBox("Thermal");
    auto* thermalForm = new QFormLayout(thermal);
    addNumberField(thermalForm, "Influx temperature [K]", &ColumnConfig::influxTemperature);
    addNumberField(thermalForm, "Wall density [kg/m3]", &ColumnConfig::wallDensity);
    addNumberField(thermalForm, "Gas thermal conductivity [W/m/K]", &ColumnConfig::gasThermalConductivity);
    addNumberField(thermalForm, "Wall thermal conductivity [W/m/K]", &ColumnConfig::wallThermalConductivity);
    addNumberField(thermalForm, "Gas-solid heat transfer [W/m2/K]", &ColumnConfig::heatTransferGasSolid);
    addNumberField(thermalForm, "Gas-wall heat transfer [W/m2/K]", &ColumnConfig::heatTransferGasWall);
    addNumberField(thermalForm, "Wall-external heat transfer [W/m2/K]", &ColumnConfig::heatTransferWallExternal);
    addNumberField(thermalForm, "Gas heat capacity [J/kg/K]", &ColumnConfig::heatCapacityGas);
    addNumberField(thermalForm, "Solid heat capacity [J/kg/K]", &ColumnConfig::heatCapacitySolid);
    addNumberField(thermalForm, "Wall heat capacity [J/kg/K]", &ColumnConfig::heatCapacityWall);
    if (column.settings.value("FluidPhase").toString().compare("Liquid", Qt::CaseInsensitive) == 0)
      addSettingNumbers(thermalForm, column.settings,
                        {{"liquidThermalConductivity", "Liquid thermal conductivity [W/m/K]", 0.6},
                         {"heatTransferLiquidSolid", "Liquid-solid heat transfer [W/m2/K]", 0.0},
                         {"heatTransferLiquidWall", "Liquid-wall heat transfer [W/m2/K]", 0.0},
                         {"heatCapacityLiquid", "Liquid heat capacity [J/kg/K]", 4180.0}});
    columnEditorLayout_->addWidget(thermal);
  }

  if (column.mode != "single")
  {
    auto* bedsGroup = new QGroupBox("Beds");
    auto* bedsLayout = new QVBoxLayout(bedsGroup);
    if (column.mode == "multibed")
    {
      auto* gridHint =
          new QLabel("Set grid intervals for every bed and finite interface, or leave all zero for a global grid.");
      gridHint->setWordWrap(true);
      bedsLayout->addWidget(gridHint);
    }
    for (qsizetype bedIndex = 0; bedIndex < column.beds.size(); ++bedIndex)
    {
      BedConfig& bed = column.beds[bedIndex];
      auto* bedGroup = new QGroupBox(bed.name.isEmpty() ? QString("Bed %1").arg(bedIndex + 1) : bed.name);
      auto* bedLayout = new QVBoxLayout(bedGroup);
      auto* bedForm = new QFormLayout;
      bedLayout->addLayout(bedForm);
      addTextField(bedForm, "Adsorbent name", bed.name,
                   [this, bedIndex](const QString& value)
                   {
                     columns_[selectedColumnIndex_].beds[bedIndex].name = value;
                     refreshColumnList();
                   });
      if (column.mode == "mixed")
      {
        addNumberField(bedForm, "Mix fraction [-]", bed.mixFraction, [this, bedIndex](double value)
                       { columns_[selectedColumnIndex_].beds[bedIndex].mixFraction = value; });
      }
      else
      {
        addNumberField(bedForm, "Length [m]", bed.length,
                       [this, bedIndex, totalLengthLabel](double value)
                       {
                         ColumnConfig& selected = columns_[selectedColumnIndex_];
                         selected.beds[bedIndex].length = value;
                         double totalLength = 0.0;
                         for (qsizetype i = 0; i < selected.beds.size(); ++i)
                         {
                           totalLength += selected.beds[i].length;
                           if (i + 1 < selected.beds.size()) totalLength += selected.beds[i].interfaceLengthAfter;
                         }
                         totalLengthLabel->setText(numberText(totalLength));
                       });
        addNumberField(
            bedForm, "Grid intervals (0 = global) [-]", static_cast<double>(bed.gridPoints),
            [this, bedIndex](double value)
            { columns_[selectedColumnIndex_].beds[bedIndex].gridPoints = std::max(0, static_cast<int>(value)); });
      }
      addNumberField(bedForm, "Void fraction [-]", bed.voidFraction, [this, bedIndex](double value)
                     { columns_[selectedColumnIndex_].beds[bedIndex].voidFraction = value; });
      addNumberField(bedForm, "Particle density [kg/m3]", bed.particleDensity, [this, bedIndex](double value)
                     { columns_[selectedColumnIndex_].beds[bedIndex].particleDensity = value; });
      addNumberField(bedForm, "Particle diameter [m]", bed.particleDiameter, [this, bedIndex](double value)
                     { columns_[selectedColumnIndex_].beds[bedIndex].particleDiameter = value; });
      addJsonField(bedForm, "Overrides (ComponentParameters / Components / AxialDispersivity)", bed.settings,
                   [this, bedIndex](const QJsonValue& value)
                   { columns_[selectedColumnIndex_].beds[bedIndex].settings = value.toObject(); });
      if (column.mode == "multibed" && bedIndex + 1 < column.beds.size())
      {
        addNumberField(bedForm, "Interface length after [m]", bed.interfaceLengthAfter,
                       [this, bedIndex, totalLengthLabel](double value)
                       {
                         ColumnConfig& selected = columns_[selectedColumnIndex_];
                         selected.beds[bedIndex].interfaceLengthAfter = std::max(0.0, value);
                         double totalLength = 0.0;
                         for (qsizetype i = 0; i < selected.beds.size(); ++i)
                         {
                           totalLength += selected.beds[i].length;
                           if (i + 1 < selected.beds.size()) totalLength += selected.beds[i].interfaceLengthAfter;
                         }
                         totalLengthLabel->setText(numberText(totalLength));
                       });
        addNumberField(bedForm, "Interface grid intervals (0 = global) [-]",
                       static_cast<double>(bed.interfaceGridPointsAfter),
                       [this, bedIndex](double value)
                       {
                         columns_[selectedColumnIndex_].beds[bedIndex].interfaceGridPointsAfter =
                             std::max(0, static_cast<int>(value));
                       });
      }
      auto* remove = new QPushButton("Remove bed");
      connect(remove, &QPushButton::clicked, this,
              [this, bedIndex]()
              {
                columns_[selectedColumnIndex_].beds.removeAt(bedIndex);
                if (columns_[selectedColumnIndex_].beds.isEmpty())
                {
                  columns_[selectedColumnIndex_].beds.push_back(makeBed("Bed 1"));
                }
                refreshAll();
              });
      bedLayout->addWidget(remove);
      bedsLayout->addWidget(bedGroup);
    }
    auto* addBed = new QPushButton("Add bed");
    connect(addBed, &QPushButton::clicked, this,
            [this]()
            {
              columns_[selectedColumnIndex_].beds.push_back(
                  makeBed(QString("Bed %1").arg(columns_[selectedColumnIndex_].beds.size() + 1)));
              refreshAll();
            });
    bedsLayout->addWidget(addBed);
    columnEditorLayout_->addWidget(bedsGroup);
  }
  columnEditorLayout_->addStretch(1);
}

void MainWindow::rebuildReactionEditor()
{
  clearLayout(reactionEditorLayout_);
  if (reactions_.isEmpty())
  {
    reactionEditorLayout_->addWidget(new QLabel("No reaction selected."));
    reactionEditorLayout_->addStretch(1);
    return;
  }

  selectedReactionIndex_ = std::clamp(selectedReactionIndex_, 0, static_cast<int>(reactions_.size()) - 1);
  ReactionConfig& reaction = reactions_[selectedReactionIndex_];

  auto* form = new QFormLayout;
  reactionEditorLayout_->addLayout(form);
  addTextField(form, "Name", reaction.name,
               [this](const QString& value)
               {
                 ReactionConfig& current = reactions_[selectedReactionIndex_];
                 current.name = value.trimmed().isEmpty() ? QString("Reaction %1").arg(selectedReactionIndex_ + 1)
                                                          : value.trimmed();
                 refreshReactionList();
               });
  addComboField(form, "Phase", reaction.phase, reactionPhaseOptions(),
                [this](const QString& value)
                {
                  ReactionConfig& current = reactions_[selectedReactionIndex_];
                  current.phase = value;
                  if (value != "PoreConcentration" && value != "BulkConcentration" &&
                      current.style != "GeneralPowerLaw")
                  {
                    current.style = "GeneralPowerLaw";
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  }
                  refreshReactionList();
                });
  addComboField(form, "Style", reaction.style, reactionStyleOptions(),
                [this](const QString& value)
                {
                  ReactionConfig& current = reactions_[selectedReactionIndex_];
                  current.style = value;
                  if (value != "GeneralPowerLaw" && current.phase != "BulkConcentration")
                  {
                    current.phase = "PoreConcentration";
                  }
                  // LH/LHHW show the adsorption constants, so the editor is rebuilt on every style change.
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  refreshReactionList();
                });
  addNumberField(form, "Site [-]", static_cast<double>(reaction.site),
                 [this](double value)
                 {
                   reactions_[selectedReactionIndex_].site = std::max(0, static_cast<int>(value));
                   refreshReactionList();
                 });

  auto defaultParticipantId = [this](const QList<ReactionParticipantConfig>& participants)
  {
    for (const ComponentConfig& component : components_)
    {
      if (!component.carrier && std::none_of(participants.cbegin(), participants.cend(),
                                             [&component](const ReactionParticipantConfig& participant)
                                             { return participant.componentId == component.id; }))
      {
        return component.id;
      }
    }
    return components_.isEmpty() ? QString{} : components_.front().id;
  };

  auto addParticipantEditor = [&](const QString& title, QList<ReactionParticipantConfig>& participants,
                                  const QString& addLabel, const QString& removeLabel)
  {
    auto* group = new QGroupBox(title);
    auto* layout = new QVBoxLayout(group);
    auto* table = new QTableWidget(static_cast<int>(participants.size()), 3);
    table->setHorizontalHeaderLabels({"Component", "Stoichiometry", "Order (blank = stoichiometry)"});
    table->setSelectionBehavior(QAbstractItemView::SelectRows);
    table->setSelectionMode(QAbstractItemView::SingleSelection);
    table->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    table->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    table->setMinimumHeight(std::clamp(76 + 30 * static_cast<int>(participants.size()), 110, 220));

    for (qsizetype row = 0; row < participants.size(); ++row)
    {
      auto* combo = new QComboBox;
      const QList<QPair<QString, QString>> options =
          componentDropdownOptions(components_, participants[row].componentId);
      sizeComboForOptions(combo, options);
      for (const auto& option : options)
      {
        combo->addItem(option.first, option.second);
      }
      const int index = combo->findData(participants[row].componentId);
      if (index >= 0)
      {
        combo->setCurrentIndex(index);
      }
      connect(combo, qOverload<int>(&QComboBox::currentIndexChanged), this,
              [this, combo, &participants, row](int)
              {
                if (rebuilding_ || row < 0 || row >= participants.size()) return;
                participants[row].componentId = combo->currentData().toString();
                refreshReactionList();
                updatePreview();
              });
      table->setCellWidget(static_cast<int>(row), 0, combo);

      auto* stoichiometry = new QTableWidgetItem(numberText(participants[row].stoichiometry));
      stoichiometry->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable);
      table->setItem(static_cast<int>(row), 1, stoichiometry);
      table->setItem(static_cast<int>(row), 2, new QTableWidgetItem(optionalText(participants[row].order)));
    }

    connect(table, &QTableWidget::itemChanged, this,
            [this, &participants](QTableWidgetItem* item)
            {
              if (rebuilding_ || item == nullptr || item->column() < 1) return;
              const int row = item->row();
              if (row < 0 || row >= static_cast<int>(participants.size())) return;
              bool ok = false;
              const double value = item->text().toDouble(&ok);
              if (item->column() == 2 && item->text().trimmed().isEmpty())
                participants[row].order.reset();
              else if (ok && item->column() == 2)
                participants[row].order = value;
              else if (ok)
              {
                participants[row].stoichiometry = value;
              }
              updatePreview();
            });

    layout->addWidget(table);
    auto* buttons = new QHBoxLayout;
    auto* addButton = new QPushButton(addLabel);
    auto* removeButton = new QPushButton(removeLabel);
    removeButton->setEnabled(!participants.isEmpty());
    buttons->addWidget(addButton);
    buttons->addWidget(removeButton);
    layout->addLayout(buttons);
    connect(addButton, &QPushButton::clicked, this,
            [this, &participants, defaultParticipantId]()
            {
              participants.push_back(ReactionParticipantConfig{defaultParticipantId(participants), 1.0});
              refreshAll();
            });
    connect(removeButton, &QPushButton::clicked, this,
            [this, table, &participants]()
            {
              if (participants.isEmpty()) return;
              const int currentRow = table->currentRow();
              const int row = currentRow >= 0 ? currentRow : static_cast<int>(participants.size()) - 1;
              if (row >= 0 && row < static_cast<int>(participants.size()))
              {
                participants.removeAt(row);
                refreshAll();
              }
            });
    reactionEditorLayout_->addWidget(group);
  };

  addParticipantEditor("Reactants", reaction.reactants, "Add reactant", "Remove reactant");
  addParticipantEditor("Products", reaction.products, "Add product", "Remove product");

  auto* kinetics = new QGroupBox("Kinetics");
  auto* kineticsForm = new QFormLayout(kinetics);
  addTemperatureLawField(kineticsForm, "k+", reaction.forwardRateCoefficient,
                         TemperatureLaw::mask(TemperatureLaw::Type::Arrhenius), false,
                         [this](const QVariant& value) { reactions_[selectedReactionIndex_].forwardRateCoefficient = value; });
  addTemperatureLawField(kineticsForm, "Keq", reaction.equilibriumConstant,
                         TemperatureLaw::mask(TemperatureLaw::Type::VantHoff), true,
                         [this](const QVariant& value) { reactions_[selectedReactionIndex_].equilibriumConstant = value; });
  addOptionalNumberField(kineticsForm, "Heat of reaction [J/mol]", &ReactionConfig::heatOfReaction);
  addNumberField(kineticsForm, "Limit [s]", &ReactionConfig::rateLimitTime);
  reactionEditorLayout_->addWidget(kinetics);
  if (reaction.style != "GeneralPowerLaw")
  {
    // theta_j = K_j c_j / (1 + sum_i K_i c_i), one K_j per component.
    auto* adsorption = new QGroupBox("Adsorption constants K_j [m^3/mol] (0: does not adsorb)");
    auto* adsorptionForm = new QFormLayout(adsorption);
    for (const ComponentConfig& component : components_)
    {
      addTemperatureLawField(adsorptionForm, component.name, reaction.adsorptionConstants.value(component.id, 0.0),
                             TemperatureLaw::mask(TemperatureLaw::Type::VantHoff), false,
                             [this, id = component.id](const QVariant& value)
                             { reactions_[selectedReactionIndex_].adsorptionConstants.insert(id, value); });
    }
    reactionEditorLayout_->addWidget(adsorption);
  }
  reactionEditorLayout_->addStretch(1);
}

void MainWindow::rebuildSimulationEditor()
{
  clearLayout(simulationEditorLayout_);
  feedTable_ = nullptr;
  if (simulations_.isEmpty())
  {
    simulationEditorLayout_->addWidget(new QLabel("No simulation selected."));
    return;
  }

  SimulationConfig& simulation = simulations_[selectedSimulationIndex_];
  auto* form = new QFormLayout;
  simulationEditorLayout_->addLayout(form);
  addTextField(form, "Display name", simulation.displayName,
               [this](const QString& value)
               {
                 simulations_[selectedSimulationIndex_].displayName = value;
                 refreshSimulationList();
               });
  addComboField(form, "Simulation type", simulation.type,
                valueOptions({"Breakthrough", "SwingAdsorption", "MixturePrediction"}),
                [this](const QString& value)
                {
                  simulations_[selectedSimulationIndex_].type = value;
                  if (value == "SwingAdsorption")
                  {
                    ensureSwingPhases(simulations_[selectedSimulationIndex_]);
                  }
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                });

  QList<QPair<QString, QString>> columnOptions;
  for (const ColumnConfig& column : columns_)
  {
    columnOptions.push_back({column.name, column.id});
  }
  addComboField(form, "Column", simulation.selectedColumnId, columnOptions,
                [this](const QString& value)
                {
                  simulations_[selectedSimulationIndex_].selectedColumnId = value;
                  refreshColumnList();
                });
  addComboField(form, "Mixture prediction method", &SimulationConfig::mixtureMethod,
                valueOptions({"SIAST", "IAST", "EI", "SEI", "SCI", "SPI"}));
  addSettingCombo(form, "IAST solver", simulation.settings, "IASTMethod", "FastIAST",
                  {"FastIAST", "NestedLoopBisection"});
  const auto* feedColumn = selectedColumnFor(simulation);
  const bool liquid = feedColumn && feedColumn->settings.value("FluidPhase").toString().compare("Liquid", Qt::CaseInsensitive) == 0;
  if (liquid)
  {
    auto pHOptions = componentDropdownOptions(components_, simulation.pHComponentId);
    pHOptions.prepend({"Select component", QString{}});
    addComboField(form, "Transported pH component", &SimulationConfig::pHComponentId, pHOptions);
  }

  auto* feedGroup = new QGroupBox("Feed components");
  auto* feedLayout = new QVBoxLayout(feedGroup);
  feedTable_ = new QTableWidget(static_cast<int>(components_.size()), liquid ? 4 : 3);
  feedTable_->setHorizontalHeaderLabels(liquid ? QStringList{"Use", "Component", "Inlet [mol/m3]", "Initial [mol/m3]"}
                                              : QStringList{"Use", "Component", "Gas mole fraction y_i"});
  feedTable_->horizontalHeader()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
  feedTable_->horizontalHeader()->setSectionResizeMode(1, QHeaderView::Stretch);
  feedTable_->horizontalHeader()->setSectionResizeMode(2, QHeaderView::ResizeToContents);
  for (qsizetype row = 0; row < components_.size(); ++row)
  {
    const int tableRow = static_cast<int>(row);
    const ComponentConfig& component = components_[row];
    const int feed = feedIndex(simulation, component.id);
    auto* useItem = new QTableWidgetItem;
    useItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
    useItem->setCheckState(feed >= 0 ? Qt::Checked : Qt::Unchecked);
    feedTable_->setItem(tableRow, 0, useItem);

    auto* nameItem = new QTableWidgetItem(component.name);
    nameItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    feedTable_->setItem(tableRow, 1, nameItem);

    auto* fractionItem = new QTableWidgetItem(feed >= 0 ? numberText(liquid ? component.settings.value("LiquidPhaseConcentration").toDouble()
                                                         : simulation.componentFeeds[feed].gasFraction) : "");
    fractionItem->setFlags(feed >= 0 ? Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable
                                     : Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    feedTable_->setItem(tableRow, 2, fractionItem);
    if (liquid)
    {
      auto* initial = new QTableWidgetItem(numberText(component.settings.value("InitialLiquidPhaseConcentration").toDouble()));
      initial->setFlags(fractionItem->flags());
      feedTable_->setItem(tableRow, 3, initial);
    }
  }
  connect(feedTable_, &QTableWidget::itemChanged, this,
          [this, liquid](QTableWidgetItem* item)
          {
            if (rebuilding_ || item == nullptr || selectedSimulationIndex_ >= static_cast<int>(simulations_.size()))
              return;
            SimulationConfig& current = simulations_[selectedSimulationIndex_];
            const int row = item->row();
            if (row < 0 || row >= static_cast<int>(components_.size())) return;
            const QString componentId = components_[row].id;
            if (item->column() == 0)
            {
              const int feed = feedIndex(current, componentId);
              if (item->checkState() == Qt::Checked && feed < 0)
              {
                current.componentFeeds.push_back({componentId, 0.0});
              }
              else if (item->checkState() != Qt::Checked && feed >= 0)
              {
                current.componentFeeds.removeAt(feed);
              }
              QTimer::singleShot(0, this, [this]() { refreshAll(); });
              return;
            }
            if (liquid && item->column() >= 2)
            {
              bool ok = false;
              const double value = item->text().toDouble(&ok);
              if (ok)
              {
                components_[row].settings.insert(item->column() == 2 ? "LiquidPhaseConcentration"
                                                                     : "InitialLiquidPhaseConcentration", value);
                updatePreview();
              }
              return;
            }
            if (item->column() == 2)
            {
              bool ok = false;
              const double value = item->text().toDouble(&ok);
              const int feed = feedIndex(current, componentId);
              if (ok && feed >= 0)
              {
                current.componentFeeds[feed].gasFraction = value;
                refreshSimulationList();
                updatePreview();
                updateFractionTotal();
              }
            }
          });
  feedLayout->addWidget(feedTable_);

  auto* feedButtons = new QHBoxLayout;
  auto* selectAll = new QPushButton("Select all");
  auto* clear = new QPushButton("Clear");
  auto* normalize = new QPushButton("Normalize");
  normalize->setVisible(!liquid);
  feedButtons->addWidget(selectAll);
  feedButtons->addWidget(clear);
  feedButtons->addWidget(normalize);
  feedLayout->addLayout(feedButtons);
  connect(selectAll, &QPushButton::clicked, this,
          [this]()
          {
            SimulationConfig& current = simulations_[selectedSimulationIndex_];
            for (const ComponentConfig& component : components_)
            {
              if (feedIndex(current, component.id) < 0)
              {
                current.componentFeeds.push_back({component.id, 0.0});
              }
            }
            refreshAll();
          });
  connect(clear, &QPushButton::clicked, this,
          [this]()
          {
            simulations_[selectedSimulationIndex_].componentFeeds.clear();
            refreshAll();
          });
  connect(normalize, &QPushButton::clicked, this,
          [this]()
          {
            SimulationConfig& current = simulations_[selectedSimulationIndex_];
            const double total = feedTotal(current);
            if (total > 0.0)
            {
              for (FeedConfig& feed : current.componentFeeds)
              {
                feed.gasFraction /= total;
              }
            }
            refreshAll();
          });
  simulationEditorLayout_->addWidget(feedGroup);

  if (simulation.type == "MixturePrediction")
  {
    auto* mixture = new QGroupBox("Pressure range");
    auto* mixtureForm = new QFormLayout(mixture);
    addNumberField(mixtureForm, "Pressure start [Pa]", &SimulationConfig::pressureStart);
    addNumberField(mixtureForm, "Pressure end [Pa]", &SimulationConfig::pressureEnd);
    addNumberField(mixtureForm, "Pressure points [-]", static_cast<double>(simulation.pressurePoints),
                   [this](double value)
                   { simulations_[selectedSimulationIndex_].pressurePoints = std::max(1, static_cast<int>(value)); });
    addComboField(mixtureForm, "Pressure scale", &SimulationConfig::pressureScale,
                  {{"Log", "log"}, {"Linear", "linear"}});
    simulationEditorLayout_->addWidget(mixture);
  }
  else
  {
    const bool isSwing = simulation.type == "SwingAdsorption";
    auto* breakthrough = new QGroupBox(isSwing ? "Shared column-step settings" : "Breakthrough settings");
    auto* breakthroughForm = new QFormLayout(breakthrough);
    addComboField(breakthroughForm, "Integrator", simulation.integrator,
                  valueOptions({
#if BUILD_SUNDIALS
                      "CVODE",
#endif
                      "RungeKutta3", "ROS34PW2"}),
                  [this](const QString& value)
                  {
                    simulations_[selectedSimulationIndex_].integrator = value;
                    QTimer::singleShot(0, this, [this]() { refreshAll(); });
                  });
    addSettingCombo(breakthroughForm, "Convection scheme", simulation.settings, "ConvectionScheme", "WENO3",
                    {"WENO3", "Upwind"});
    addSettingCombo(breakthroughForm, "Inlet condition", simulation.settings, "InletCondition", "Danckwerts",
                    {"Danckwerts", "Dirichlet"});
    addSettingCombo(breakthroughForm, "Flow direction (feed end)", simulation.settings, "FlowDirection", "Forward",
                    {"Forward", "Reverse"});
    if (simulation.integrator == "CVODE")
    {
      addSettingNumbers(breakthroughForm, simulation.settings,
                        {{"CVODERelativeTolerance", "Relative tolerance [-]", 1.0e-5},
                         {"CVODEAbsoluteToleranceConcentration", "Concentration absolute tolerance [mol/m3]", 1.0e-6},
                         {"CVODEAbsoluteToleranceLoading", "Loading absolute tolerance [mol/kg]", 1.0e-7},
                         {"CVODEAbsoluteToleranceTemperature", "Temperature absolute tolerance [K]", 1.0e-4},
                         {"CVODEMaximumTimeStep", "Maximum internal step [s] (-1 = auto, 0 = unlimited)", -1.0},
                         {"CVODEMaximumOrder", "Maximum BDF order [-] (1-5)", 5.0, true},
                         {"CVODEKrylovDimension", "Krylov dimension [-]", 30.0, true}});
      addSettingCombo(breakthroughForm, "CVODE linear solver", simulation.settings, "CVODELinearSolver", "SPGMR",
                      {"Dense", "SPGMR"});
    }
    if (simulation.integrator == "ROS34PW2")
    {
      addSettingNumbers(breakthroughForm, simulation.settings,
                        {{"ROS34PW2RelativeTolerance", "Relative tolerance [-]", 1.0e-6},
                         {"ROS34PW2AbsoluteToleranceConcentration", "Concentration absolute tolerance [mol/m3]", 1.0e-7},
                         {"ROS34PW2AbsoluteToleranceLoading", "Loading absolute tolerance [mol/kg]", 1.0e-8},
                         {"ROS34PW2AbsoluteToleranceTemperature", "Temperature absolute tolerance [K]", 1.0e-5},
                         {"ROS34PW2MaximumTimeStep", "Maximum step [s] (0 = unlimited)", 0.0}});
    }
    addJsonField(breakthroughForm, "Multibed grid positions [m] (empty = automatic)",
                 simulation.settings.value("ColumnDistances").toArray(),
                 [this](const QJsonValue& value)
                 {
                   auto& settings = simulations_[selectedSimulationIndex_].settings;
                   if (value.toArray().isEmpty()) settings.remove("ColumnDistances");
                   else settings.insert("ColumnDistances", value);
                 });
    addTextField(breakthroughForm, "Initial column state file", simulation.settings.value("ReadColumnFile").toString(),
                 [this](const QString& value)
                 {
                   auto& settings = simulations_[selectedSimulationIndex_].settings;
                   if (value.trimmed().isEmpty()) settings.remove("ReadColumnFile");
                   else settings.insert("ReadColumnFile", value.trimmed());
                 });
    addNumberField(breakthroughForm, "Grid points [-]", static_cast<double>(simulation.gridPoints), [this](double value)
                   { simulations_[selectedSimulationIndex_].gridPoints = std::max(1, static_cast<int>(value)); });
    if (simulation.integrator == "RungeKutta3")  // CVODE and ROS34PW2 choose their own internal steps
    {
      addNumberField(breakthroughForm, "Time step [s]", &SimulationConfig::timeStep);
      addNumberField(breakthroughForm, "Initialization steps [-]", static_cast<double>(simulation.initSteps),
                     [this](double value)
                     { simulations_[selectedSimulationIndex_].initSteps = std::max(0, static_cast<int>(value)); });
    }
    addNumberField(breakthroughForm, "Write interval [s]", &SimulationConfig::writeInterval);
    if (!isSwing)
    {
      addCheckField(breakthroughForm, "Automatic simulation time", simulation.autoTime,
                    [this](bool value)
                    {
                      simulations_[selectedSimulationIndex_].autoTime = value;
                      QTimer::singleShot(0, this, [this]() { refreshAll(); });
                    });
      if (!simulation.autoTime)
      {
        addNumberField(breakthroughForm, "Simulation time [s]", &SimulationConfig::simulationTime);
      }
    }
    simulationEditorLayout_->addWidget(breakthrough);

    if (isSwing)
    {
      ensureSwingPhases(simulation);
      auto* phasesGroup = new QGroupBox("Swing phases");
      auto* phasesLayout = new QVBoxLayout(phasesGroup);
      for (qsizetype phaseIndex = 0; phaseIndex < simulation.swingPhases.size(); ++phaseIndex)
      {
        SwingPhaseConfig& phase = simulation.swingPhases[phaseIndex];
        auto* phaseGroup = new QGroupBox(phase.name.isEmpty() ? QString("Phase %1").arg(phaseIndex + 1) : phase.name);
        auto* phaseLayout = new QVBoxLayout(phaseGroup);
        auto* phaseForm = new QFormLayout;
        phaseLayout->addLayout(phaseForm);
        addTextField(phaseForm, "Name", phase.name,
                     [this, phaseIndex](const QString& value)
                     {
                       simulations_[selectedSimulationIndex_].swingPhases[phaseIndex].name = value;
                       refreshSimulationList();
                     });
        addNumberField(phaseForm, "Temperature [K]", phase.temperature, [this, phaseIndex](double value)
                       { simulations_[selectedSimulationIndex_].swingPhases[phaseIndex].temperature = value; });
        addNumberField(phaseForm, "Inlet pressure [Pa]", phase.inletPressure, [this, phaseIndex](double value)
                       { simulations_[selectedSimulationIndex_].swingPhases[phaseIndex].inletPressure = value; });
        addNumberField(phaseForm, "Duration [s]", phase.duration, [this, phaseIndex](double value)
                       { simulations_[selectedSimulationIndex_].swingPhases[phaseIndex].duration = value; });
        // Reverse feeds the phase at z = L; empty keeps the FlowDirection of the simulation.
        addComboField(phaseForm, "Flow direction", phase.flowDirection,
                      {{"Simulation setting", ""}, {"Forward (feed at z = 0)", "Forward"},
                       {"Reverse (feed at z = L)", "Reverse"}},
                      [this, phaseIndex](const QString& value)
                      { simulations_[selectedSimulationIndex_].swingPhases[phaseIndex].flowDirection = value; });

        auto* remove = new QPushButton("Remove phase");
        connect(remove, &QPushButton::clicked, this,
                [this, phaseIndex]()
                {
                  SimulationConfig& current = simulations_[selectedSimulationIndex_];
                  current.swingPhases.removeAt(phaseIndex);
                  ensureSwingPhases(current);
                  refreshAll();
                });
        phaseLayout->addWidget(remove);
        phasesLayout->addWidget(phaseGroup);
      }

      auto* phaseButtons = new QHBoxLayout;
      auto* addPressurePhase = new QPushButton("Add pressure phase");
      auto* addTemperaturePhase = new QPushButton("Add temperature phase");
      phaseButtons->addWidget(addPressurePhase);
      phaseButtons->addWidget(addTemperaturePhase);
      phasesLayout->addLayout(phaseButtons);
      connect(addPressurePhase, &QPushButton::clicked, this,
              [this]()
              {
                SimulationConfig& current = simulations_[selectedSimulationIndex_];
                const SwingPhaseConfig base =
                    current.swingPhases.isEmpty() ? SwingPhaseConfig{} : current.swingPhases.last();
                current.swingPhases.push_back(
                    makeSwingPhase(QString("Pressure phase %1").arg(current.swingPhases.size() + 1), base.temperature,
                                   std::max(1000.0, 0.1 * base.inletPressure), base.duration));
                refreshAll();
              });
      connect(addTemperaturePhase, &QPushButton::clicked, this,
              [this]()
              {
                SimulationConfig& current = simulations_[selectedSimulationIndex_];
                const SwingPhaseConfig base =
                    current.swingPhases.isEmpty() ? SwingPhaseConfig{} : current.swingPhases.last();
                current.swingPhases.push_back(
                    makeSwingPhase(QString("Temperature phase %1").arg(current.swingPhases.size() + 1),
                                   base.temperature + 50.0, base.inletPressure, base.duration));
                refreshAll();
              });
      simulationEditorLayout_->addWidget(phasesGroup);
    }
  }
  simulationEditorLayout_->addStretch(1);
  updateFractionTotal();
}

void MainWindow::updatePreview()
{
  if (jsonPreview_ == nullptr)
  {
    return;
  }
  jsonPreview_->setPlainText(QString::fromUtf8(QJsonDocument(buildSimulationJson()).toJson(QJsonDocument::Indented)));
}

void MainWindow::updateFractionTotal()
{
  const SimulationConfig* simulation = activeSimulation();
  const double total = simulation == nullptr ? 0.0 : feedTotal(*simulation);
  const auto* column = simulation ? selectedColumnFor(*simulation) : nullptr;
  if (column && column->settings.value("FluidPhase").toString().compare("Liquid", Qt::CaseInsensitive) == 0)
  {
    fractionTotalLabel_->setText("Liquid concentrations");
    fractionTotalLabel_->setStyleSheet({});
    return;
  }
  fractionTotalLabel_->setText(QString::number(total, 'f', 3));
  fractionTotalLabel_->setStyleSheet(std::abs(total - 1.0) < 0.001 ? "color: #276749;" : "color: #b7791f;");
}

void MainWindow::refreshRunPanel()
{
  if (runIdLabel_ == nullptr || runLog_ == nullptr)
  {
    return;
  }

  const SimulationConfig* simulation = activeSimulation();
  if (simulation == nullptr || simulation->lastRun.id.isEmpty())
  {
    runIdLabel_->setText("No run");
    if (!runLog_->toPlainText().isEmpty())
    {
      runLog_->clear();
    }
    updateRunBadge("idle");
    updateRunControls("idle");
    return;
  }

  runIdLabel_->setText(simulation->lastRun.id);
  if (runLog_->toPlainText() != simulation->lastRunLog)
  {
    runLog_->setPlainText(simulation->lastRunLog);
    runLog_->moveCursor(QTextCursor::End);
  }
  updateRunBadge(simulation->lastRun.status);
  updateRunControls(simulation->lastRun.status);
}

void MainWindow::updateRunControls(const QString& status)
{
  if (runButton_ == nullptr)
  {
    return;
  }

  const bool selectedBusy = status == "running" || status == "canceling";
  const bool runnerBusy = runner_.isRunning();
  runButton_->setDisabled(runnerBusy);
  runTopButton_->setDisabled(runnerBusy);
  cancelButton_->setEnabled(status == "running" && runnerBusy);
  cancelTopButton_->setEnabled(status == "running" && runnerBusy);
  const SimulationConfig* simulation = activeSimulation();
  const bool hasRun = simulation != nullptr && !simulation->lastRun.id.isEmpty();
  analysisButton_->setEnabled(hasRun && !selectedBusy);
  analysisTopButton_->setEnabled(hasRun && !selectedBusy);
}

void MainWindow::updateRunBadge(const QString& status)
{
  const QString label = status.isEmpty() ? "Idle" : status.left(1).toUpper() + status.mid(1);
  runBadge_->setText(label);
  QString color = "#5e6872";
  if (status == "running")
  {
    color = "#096b72";
  }
  else if (status == "completed")
  {
    color = "#276749";
  }
  else if (status == "failed" || status == "canceled")
  {
    color = "#a43d36";
  }
  runBadge_->setStyleSheet(
      QString("padding: 3px 8px; border: 1px solid #cfd6dc; border-radius: 5px; color: %1;").arg(color));
}

QString MainWindow::newId(const QString& prefix)
{
  ++idCounter_;
  return QString("%1-%2").arg(prefix).arg(idCounter_);
}

PhysisorptionSiteConfig MainWindow::makePhySite(const QString& type, const QVariantMap& values)
{
  PhysisorptionSiteConfig site;
  site.id = newId("phys");
  site.type = type;
  site.parameters = defaultsFor(phyDefs(type), values);
  return site;
}

ChemisorptionSiteConfig MainWindow::makeChemSite(const QString& type, const QVariantMap& values)
{
  ChemisorptionSiteConfig site;
  site.id = newId("chem");
  site.type = type;
  site.parameters = defaultsFor(chemDefs(type), values);
  site.isotherm =
      makePhySite("Langmuir", {{"qmax", values.value("maximumLoading", 0.8)}, {"affinityConstant", 2.5e-5}});
  return site;
}

BedConfig MainWindow::makeBed(const QString& name)
{
  BedConfig bed;
  bed.id = newId("bed");
  bed.name = name;
  return bed;
}

ColumnConfig MainWindow::makeColumn(const QString& name, const QString& id)
{
  ColumnConfig column;
  column.id = id.isEmpty() ? newId("column") : id;
  column.name = name;
  column.beds.push_back(makeBed("Bed 1"));
  column.beds.push_back(makeBed("Bed 2"));
  return column;
}

SwingPhaseConfig MainWindow::makeSwingPhase(const QString& name, double temperature, double inletPressure,
                                            double duration)
{
  SwingPhaseConfig phase;
  phase.id = newId("phase");
  phase.name = name;
  phase.temperature = temperature;
  phase.inletPressure = inletPressure;
  phase.duration = duration;
  return phase;
}

ReactionConfig MainWindow::makeReaction(const QString& name)
{
  ReactionConfig reaction;
  reaction.id = newId("reaction");
  reaction.name = name.isEmpty() ? "Reaction" : name;

  QList<QString> participants;
  for (const ComponentConfig& component : components_)
  {
    if (!component.carrier)
    {
      participants.push_back(component.id);
    }
  }
  if (participants.isEmpty() && !components_.isEmpty())
  {
    participants.push_back(components_.front().id);
  }
  if (!participants.isEmpty())
  {
    reaction.reactants.push_back(ReactionParticipantConfig{participants.front(), 1.0});
  }
  if (participants.size() >= 2)
  {
    reaction.products.push_back(ReactionParticipantConfig{participants[1], 1.0});
  }
  else if (!participants.isEmpty())
  {
    reaction.products.push_back(ReactionParticipantConfig{participants.front(), 1.0});
  }
  return reaction;
}

SimulationConfig MainWindow::makeSimulation(const QString& name, const QString& id, const QList<FeedConfig>& feeds)
{
  SimulationConfig simulation;
  simulation.id = id.isEmpty() ? newId("simulation") : id;
  simulation.displayName = name;
  simulation.selectedColumnId = columns_.isEmpty() ? "column-1" : columns_.front().id;
  simulation.componentFeeds = feeds;
  ensureSwingPhases(simulation);
  return simulation;
}

void MainWindow::ensureSwingPhases(SimulationConfig& simulation)
{
  if (!simulation.swingPhases.isEmpty())
  {
    return;
  }

  const ColumnConfig* column = selectedColumnFor(simulation);
  const double temperature = column == nullptr ? 300.0 : column->temperature;
  const double inletPressure = column == nullptr ? 100000.0 : column->inletPressure;
  simulation.swingPhases.push_back(makeSwingPhase("Adsorption", temperature, inletPressure, 100000));
  simulation.swingPhases.push_back(
      makeSwingPhase("Regeneration", temperature, std::max(1000.0, 0.1 * inletPressure), 100000));
}

SimulationConfig* MainWindow::activeSimulation()
{
  if (simulations_.isEmpty())
  {
    return nullptr;
  }
  selectedSimulationIndex_ = std::clamp(selectedSimulationIndex_, 0, static_cast<int>(simulations_.size()) - 1);
  return &simulations_[selectedSimulationIndex_];
}

const SimulationConfig* MainWindow::activeSimulation() const
{
  if (simulations_.isEmpty())
  {
    return nullptr;
  }
  const int index = std::clamp(selectedSimulationIndex_, 0, static_cast<int>(simulations_.size()) - 1);
  return &simulations_[index];
}

ColumnConfig* MainWindow::selectedColumnFor(SimulationConfig& simulation)
{
  auto it = std::find_if(columns_.begin(), columns_.end(), [&simulation](const ColumnConfig& column)
                         { return column.id == simulation.selectedColumnId; });
  if (it == columns_.end())
  {
    return columns_.isEmpty() ? nullptr : &columns_.front();
  }
  return &(*it);
}

const ColumnConfig* MainWindow::selectedColumnFor(const SimulationConfig& simulation) const
{
  auto it = std::find_if(columns_.cbegin(), columns_.cend(), [&simulation](const ColumnConfig& column)
                         { return column.id == simulation.selectedColumnId; });
  if (it == columns_.cend())
  {
    return columns_.isEmpty() ? nullptr : &columns_.front();
  }
  return &(*it);
}

SimulationConfig* MainWindow::simulationById(const QString& simulationId)
{
  auto it = std::find_if(simulations_.begin(), simulations_.end(),
                         [&simulationId](const SimulationConfig& simulation) { return simulation.id == simulationId; });
  return it == simulations_.end() ? nullptr : &(*it);
}

SimulationConfig* MainWindow::simulationForRunId(const QString& runId)
{
  if (runId.isEmpty())
  {
    return nullptr;
  }
  auto it = std::find_if(simulations_.begin(), simulations_.end(),
                         [&runId](const SimulationConfig& simulation) { return simulation.lastRun.id == runId; });
  return it == simulations_.end() ? nullptr : &(*it);
}

int MainWindow::feedIndex(const SimulationConfig& simulation, const QString& componentId) const
{
  for (int i = 0; i < simulation.componentFeeds.size(); ++i)
  {
    if (simulation.componentFeeds[i].componentId == componentId)
    {
      return i;
    }
  }
  return -1;
}

double MainWindow::feedTotal(const SimulationConfig& simulation) const
{
  double total = 0.0;
  for (const FeedConfig& feed : simulation.componentFeeds)
  {
    total += feed.gasFraction;
  }
  return total;
}

QJsonObject MainWindow::buildSimulationJson() const
{
  const auto* simulation = activeSimulation();
  const auto modeled = buildModeledSimulationJson();
  if (simulation == nullptr || simulation->originalInput.isEmpty()) return modeled;
  auto result = LabInput::applyChanges(simulation->originalInput, simulation->originalModel, modeled);
  QHash<QString, QString> names;
  const auto originalComponents = simulation->originalInput.value("Components").toArray();
  for (qsizetype i = 0; i < originalComponents.size(); ++i)
  {
    const auto id = QString("component-%1").arg(i + 1);
    const auto component = std::find_if(components_.begin(), components_.end(),
                                      [&id](const auto& item) { return item.id == id; });
    const auto oldName = originalComponents[i].toObject().value("Name").toString();
    if (component != components_.end() && component->name != oldName) names.insert(oldName, component->name);
  }
  return names.isEmpty() ? result : LabInput::renameComponentReferences(result, names);
}

QJsonObject MainWindow::buildModeledSimulationJson() const
{
  const SimulationConfig* simulation = activeSimulation();
  if (simulation == nullptr)
  {
    return {};
  }
  const ColumnConfig* column = selectedColumnFor(*simulation);
  const double temperature = column == nullptr ? 300.0 : column->temperature;
  if (simulation->type == "MixturePrediction")
  {
    QJsonObject json{
        {"SimulationType", "MixturePrediction"},
        {"DisplayName", simulation->displayName},
        {"Temperature", temperature},
        {"PressureStart", simulation->pressureStart},
        {"PressureEnd", simulation->pressureEnd},
        {"NumberOfPressurePoints", simulation->pressurePoints},
        {"PressureScale", simulation->pressureScale},
        {"MixturePredictionMethod", simulation->mixtureMethod},
        {"Components", componentsJsonFor(*simulation)},
    };
    if (column) mergeSettings(json, column->settings);
    mergeSettings(json, simulation->settings);
    if (!simulation->pHComponentId.isEmpty())
      json.insert("pHComponent", static_cast<int>(componentExportOrderFor(*simulation).indexOf(simulation->pHComponentId)));
    return json;
  }

  const bool isSwing = simulation->type == "SwingAdsorption";

  const QString integrator =
      (simulation->integrator == "RungeKutta3" || simulation->integrator == "ROS34PW2") ? simulation->integrator : "CVODE";
  QJsonObject json{
      {"SimulationType", isSwing ? "SwingAdsorption" : "Breakthrough"},
      {"DisplayName", simulation->displayName},
      {"Temperature", temperature},
      {"BreakthroughIntegrator", integrator},
      {"WriteInterval", simulation->writeInterval},
      {"MixturePredictionMethod", simulation->mixtureMethod},
      {"Components", componentsJsonFor(*simulation)},
  };
  // The phase durations set the simulated time of a swing run.
  if (!isSwing)
  {
    json.insert("SimulationTime", simulation->autoTime ? QJsonValue("auto") : QJsonValue(simulation->simulationTime));
  }
  // CVODE and ROS34PW2 choose their own internal steps and take no TimeStep.
  if (integrator == "RungeKutta3")
  {
    json.insert("TimeStep", simulation->timeStep);
    json.insert("NumberOfInitTimeSteps", simulation->initSteps);
  }
  const QJsonArray reactions = reactionsJsonFor(*simulation);
  if (!reactions.isEmpty())
  {
    json.insert("Reactions", reactions);
  }
  if (isSwing)
  {
    json.insert("SwingAdsorptionPhases", swingPhasesToJson(*simulation));
  }
  if (column != nullptr)
  {
    applyColumnToBreakthrough(json, *column, *simulation);
  }
  mergeSettings(json, simulation->settings);
  if (!simulation->pHComponentId.isEmpty())
    json.insert("pHComponent", static_cast<int>(componentExportOrderFor(*simulation).indexOf(simulation->pHComponentId)));
  return json;
}

QList<QString> MainWindow::componentExportOrderFor(const SimulationConfig& simulation) const
{
  QList<QString> componentIds;
  auto appendComponent = [&](const QString& componentId)
  {
    if (!componentId.isEmpty() && !componentIds.contains(componentId) && stringSetContains(components_, componentId))
    {
      componentIds.push_back(componentId);
    }
  };

  for (const FeedConfig& feed : simulation.componentFeeds)
  {
    appendComponent(feed.componentId);
  }
  if (simulation.type != "MixturePrediction")
  {
    for (const ReactionConfig& reaction : reactions_)
    {
      for (const ReactionParticipantConfig& participant : reaction.reactants)
      {
        appendComponent(participant.componentId);
      }
      for (const ReactionParticipantConfig& participant : reaction.products)
      {
        appendComponent(participant.componentId);
      }
    }
  }
  return componentIds;
}

QJsonArray MainWindow::componentsJsonFor(const SimulationConfig& simulation) const
{
  QJsonArray components;
  for (const QString& componentId : componentExportOrderFor(simulation))
  {
    auto it = std::find_if(components_.cbegin(), components_.cend(),
                           [&componentId](const ComponentConfig& component) { return component.id == componentId; });
    if (it != components_.cend())
    {
      const int feed = feedIndex(simulation, componentId);
      components.append(componentToJson(*it, feed >= 0 ? simulation.componentFeeds[feed].gasFraction : 0.0));
    }
  }
  return components;
}

QJsonArray MainWindow::reactionsJsonFor(const SimulationConfig& simulation) const
{
  QJsonArray reactions;
  if (simulation.type == "MixturePrediction")
  {
    return reactions;
  }
  for (const ReactionConfig& reaction : reactions_)
  {
    reactions.append(reactionToJson(reaction, simulation));
  }
  return reactions;
}

QJsonObject MainWindow::componentToJson(const ComponentConfig& component, double gasFraction) const
{
  QJsonObject json{{"Name", component.name}, {"GasPhaseMolFraction", gasFraction}};
  mergeSettings(json, component.settings);
  if (component.carrier)
  {
    json.insert("CarrierGas", true);
  }
  if (component.water)
  {
    json.insert("Water", true);
  }
  addOptionalNumber(json, "MassTransferCoefficient", component.massTransfer);
  addOptionalNumber(json, "AxialDispersionCoefficient", component.axialDispersion);
  addOptionalNumber(json, "MolecularWeight", component.molecularWeight);
  addOptionalNumber(json, "HeatOfAdsorption", component.heatOfAdsorption);
  addOptionalNumber(json, "ReferenceTemperature", component.referenceTemperature);
  if (!component.carrier && !component.physisorptionSites.isEmpty())
  {
    QJsonArray sites;
    for (const PhysisorptionSiteConfig& site : component.physisorptionSites)
    {
      sites.append(phySiteToJson(site));
    }
    json.insert("PhysisorptionSites", sites);
  }
  if (!component.carrier && !component.chemisorptionSites.isEmpty())
  {
    QJsonArray sites;
    for (const ChemisorptionSiteConfig& site : component.chemisorptionSites)
    {
      sites.append(chemSiteToJson(site));
    }
    json.insert("ChemisorptionSites", sites);
  }
  return json;
}

QJsonObject MainWindow::reactionToJson(const ReactionConfig& reaction, const SimulationConfig& simulation) const
{
  const QList<QString> componentIds = componentExportOrderFor(simulation);
  auto participantsToJson = [&](const QList<ReactionParticipantConfig>& participants)
  {
    QJsonArray array;
    for (const ReactionParticipantConfig& participant : participants)
    {
      const qsizetype componentIndex = componentIds.indexOf(participant.componentId);
      if (componentIndex >= 0)
      {
        array.append(static_cast<int>(componentIndex));
        continue;
      }
      bool ok = false;
      const int rawIndex = participant.componentId.toInt(&ok);
      if (ok && rawIndex >= 0)
      {
        array.append(rawIndex);
      }
      else
      {
        array.append(componentNameForId(components_, participant.componentId));
      }
    }
    return array;
  };
  auto stoichiometryToJson = [](const QList<ReactionParticipantConfig>& participants)
  {
    QJsonArray array;
    for (const ReactionParticipantConfig& participant : participants)
    {
      array.append(participant.stoichiometry);
    }
    return array;
  };

  const QJsonArray reactants = participantsToJson(reaction.reactants);
  const QJsonArray products = participantsToJson(reaction.products);
  QJsonObject json{{"Phase", reaction.phase},
                     {"Style", reaction.style},
                     {"Site", std::max(0, reaction.site)},
                     {"Reactants", reactants},
                     {"Products", products},
                     {"Stoichiometry", QJsonObject{{"Reactants", stoichiometryToJson(reaction.reactants)},
                                                   {"Products", stoichiometryToJson(reaction.products)}}},
                     {"RateLimitTime", reaction.rateLimitTime}};
  QJsonObject kinetics{{"ForwardRateCoefficient", QJsonValue::fromVariant(reaction.forwardRateCoefficient)},
                       {"EquilibriumConstant", QJsonValue::fromVariant(reaction.equilibriumConstant)}};
  addOptionalNumber(kinetics, "HeatOfReaction", reaction.heatOfReaction);
  if (reaction.style != "GeneralPowerLaw")
  {
    // By name, for the components this simulation exports; a constant of 0 is the same as none.
    QJsonObject adsorptionConstants;
    for (auto it = reaction.adsorptionConstants.cbegin(); it != reaction.adsorptionConstants.cend(); ++it)
    {
      const bool zero = it.value().typeId() != QMetaType::QVariantMap && it.value().toDouble() == 0.0;
      if (zero || !componentIds.contains(it.key())) continue;
      adsorptionConstants.insert(componentNameForId(components_, it.key()), QJsonValue::fromVariant(it.value()));
    }
    kinetics.insert("AdsorptionConstants", adsorptionConstants);
  }
  json.insert("Kinetics", kinetics);
  for (const auto& pair : {qMakePair(QString("ForwardOrders"), &reaction.reactants),
                           qMakePair(QString("BackwardOrders"), &reaction.products)})
  {
    if (std::any_of(pair.second->begin(), pair.second->end(), [](const auto& item) { return item.order.has_value(); }))
    {
      QJsonArray orders;
      for (const auto& item : *pair.second) orders.append(item.order.value_or(item.stoichiometry));
      json.insert(pair.first, orders);
    }
  }
  return json;
}

QJsonObject MainWindow::phySiteToJson(const PhysisorptionSiteConfig& site) const
{
  QJsonArray parameters;
  for (const ParameterDefinition& definition : phyDefs(site.type))
  {
    parameters.append(parameterJson(site.parameters, definition));
  }
  QJsonObject json = site.settings;
  json.insert("Type", site.type);
  json.insert("Parameters", parameters);
  return json;
}

QJsonObject MainWindow::chemSiteToJson(const ChemisorptionSiteConfig& site) const
{
  QJsonObject json = phySiteToJson(site.isotherm);
  QJsonObject kinetics;
  for (const ParameterDefinition& definition : chemDefs(site.type))
  {
    if (definition.key == "maximumLoading" || definition.key == "heatOfChemisorption") continue;
    QString key = definition.key;
    if (key == "rateCoefficient" && site.type == "FirstOrder") key = "MassTransferCoefficient";
    else key[0] = key[0].toUpper();
    kinetics.insert(key, parameterJson(site.parameters, definition));
  }
  json.insert("RateEquation", site.type);
  json.insert("RateEquationParameters", kinetics);
  json.insert("MaximumLoading", site.parameters.value("maximumLoading", 0.8).toDouble());
  json.insert("HeatOfAdsorption", site.parameters.value("heatOfChemisorption", 0.0).toDouble());
  return json;
}

QJsonArray MainWindow::swingPhasesToJson(const SimulationConfig& simulation) const
{
  QJsonArray phases;
  for (const SwingPhaseConfig& phase : simulation.swingPhases)
  {
    QJsonObject object{{"Name", phase.name.isEmpty() ? "Phase" : phase.name},
                       {"Temperature", phase.temperature},
                       {"InletPressure", phase.inletPressure},
                       {"Duration", phase.duration}};
    if (!phase.flowDirection.isEmpty()) object.insert("FlowDirection", phase.flowDirection);
    for (auto key = phase.otherKeys.begin(); key != phase.otherKeys.end(); ++key) object.insert(key.key(), key.value());
    phases.append(object);
  }
  return phases;
}

void MainWindow::applyColumnToBreakthrough(QJsonObject& simulationJson, const ColumnConfig& column,
                                           const SimulationConfig& simulation) const
{
  mergeSettings(simulationJson, column.settings);
  simulationJson.insert("BoundaryCondition", column.boundaryCondition);
  simulationJson.insert("DynamicViscosity", column.dynamicViscosity);
  simulationJson.insert("ColumnVoidFraction", column.voidFraction);
  simulationJson.insert("ParticleDensity", column.particleDensity);
  simulationJson.insert("ParticleDiameter", column.particleDiameter);
  simulationJson.insert("InletPressure", column.inletPressure);
  simulationJson.insert("PressureGradient", column.pressureGradient);
  simulationJson.insert("ColumnEntranceVelocity", column.velocity);
  addOptionalNumber(simulationJson, "OutletPressure", column.outletPressure);
  QJsonObject geometry;
  if (column.geometryType == "Monolith")
  {
    geometry.insert("Type", "Monolith");
    geometry.insert("ChannelShape", column.channelShape);
    geometry.insert("InternalChannelDimension", column.internalChannelDimension);
    geometry.insert("OuterDiameter", column.outerDiameter);
    geometry.insert("NumberOfChannels", std::max(1, column.numberOfChannels));
    geometry.insert("WashcoatThickness", column.washcoatThickness);
    addOptionalNumber(geometry, "WashcoatVolumePerChannelVolume", column.washcoatVolumePerChannelVolume);
  }
  else
  {
    geometry.insert("Type", "PackedBed");
    geometry.insert("ColumnVoidFraction", column.voidFraction);
    geometry.insert("ParticleDiameter", column.particleDiameter);
    geometry.insert("InternalDiameter", column.internalDiameter);
    geometry.insert("OuterDiameter", column.outerDiameter);
  }
  mergeSettings(geometry, column.geometrySettings);
  if (column.geometryType != "Monolith") geometry.remove("TaylorArisDispersion");
  simulationJson.insert("Geometry", geometry);

  if ((column.mode == "multibed" || column.mode == "mixed") && !column.beds.isEmpty())
  {
    double length = 0.0;
    QJsonArray adsorbents;
    QJsonArray sections;
    for (qsizetype i = 0; i < column.beds.size(); ++i)
    {
      const BedConfig& bed = column.beds[i];
      const QString name = bed.name.isEmpty() ? QString("Bed %1").arg(i + 1) : bed.name;
      QJsonObject adsorbent{{"Name", name},
                            {"ColumnVoidFraction", bed.voidFraction},
                            {"ParticleDensity", bed.particleDensity},
                            {"ParticleDiameter", bed.particleDiameter}};
      mergeSettings(adsorbent, bed.settings);
      if (column.mode == "mixed")
      {
        adsorbent.insert("MixFraction", bed.mixFraction);
      }
      else
      {
        length += bed.length;
        QJsonObject section{{"Adsorbent", name}, {"Length", bed.length}};
        if (bed.gridPoints > 0) section.insert("NumberOfGridPoints", bed.gridPoints);
        sections.append(section);
        if (i + 1 < column.beds.size())
        {
          length += bed.interfaceLengthAfter;
          QJsonObject interface{{"InterfaceLength", bed.interfaceLengthAfter}};
          if (bed.interfaceLengthAfter > 0.0 && bed.interfaceGridPointsAfter > 0)
          {
            interface.insert("NumberOfGridPoints", bed.interfaceGridPointsAfter);
          }
          sections.append(interface);
        }
      }
      adsorbents.append(adsorbent);
    }
    simulationJson.insert("ColumnLength", column.mode == "mixed" ? column.length : length);
    simulationJson.insert("NumberOfGridPoints", std::max(1, simulation.gridPoints));
    simulationJson.insert("Adsorbents", adsorbents);
    if (column.mode == "multibed")
    {
      simulationJson.insert("ColumnSections", sections);
    }
  }
  else
  {
    simulationJson.insert("ColumnLength", column.length);
    simulationJson.insert("NumberOfGridPoints", std::max(1, simulation.gridPoints));
  }

  if (column.energyBalance)
  {
    simulationJson.insert("energyBalance", true);
    simulationJson.insert("InfluxTemperature", column.influxTemperature);
    simulationJson.insert("wallDensity", column.wallDensity);
    simulationJson.insert("gasThermalConductivity", column.gasThermalConductivity);
    simulationJson.insert("wallThermalConductivity", column.wallThermalConductivity);
    simulationJson.insert("heatTransferGasSolid", column.heatTransferGasSolid);
    simulationJson.insert("heatTransferGasWall", column.heatTransferGasWall);
    simulationJson.insert("heatTransferWallExternal", column.heatTransferWallExternal);
    simulationJson.insert("heatCapacityGas", column.heatCapacityGas);
    simulationJson.insert("heatCapacitySolid", column.heatCapacitySolid);
    simulationJson.insert("heatCapacityWall", column.heatCapacityWall);
  }
}

QJsonObject MainWindow::buildStateJson() const
{
  QJsonArray componentArray;
  for (const ComponentConfig& component : components_)
  {
    QJsonObject object{{"Id", component.id},
                       {"Name", component.name},
                       {"CarrierGas", component.carrier},
                       {"Water", component.water},
                       {"Settings", component.settings},
                       {"PhysisorptionSites", phySitesStateJson(component.physisorptionSites)},
                       {"ChemisorptionSites", chemSitesStateJson(component.chemisorptionSites)}};
    addOptionalNumber(object, "MassTransferCoefficient", component.massTransfer);
    addOptionalNumber(object, "AxialDispersionCoefficient", component.axialDispersion);
    addOptionalNumber(object, "MolecularWeight", component.molecularWeight);
    addOptionalNumber(object, "HeatOfAdsorption", component.heatOfAdsorption);
    addOptionalNumber(object, "ReferenceTemperature", component.referenceTemperature);
    componentArray.append(object);
  }

  QJsonArray columnArray;
  for (const ColumnConfig& column : columns_)
  {
    QJsonArray beds;
    for (const BedConfig& bed : column.beds)
    {
      beds.append(QJsonObject{{"Id", bed.id},
                              {"Name", bed.name},
                              {"Length", bed.length},
                              {"GridPoints", bed.gridPoints},
                              {"MixFraction", bed.mixFraction},
                              {"VoidFraction", bed.voidFraction},
                              {"ParticleDensity", bed.particleDensity},
                              {"ParticleDiameter", bed.particleDiameter},
                              {"InterfaceLengthAfter", bed.interfaceLengthAfter},
                              {"InterfaceGridPointsAfter", bed.interfaceGridPointsAfter},
                              {"Settings", bed.settings}});
    }
    QJsonObject object{{"Id", column.id},
                       {"Name", column.name},
                       {"Mode", column.mode},
                       {"BoundaryCondition", column.boundaryCondition},
                       {"Temperature", column.temperature},
                       {"DynamicViscosity", column.dynamicViscosity},
                       {"VoidFraction", column.voidFraction},
                       {"ParticleDensity", column.particleDensity},
                       {"ParticleDiameter", column.particleDiameter},
                       {"InletPressure", column.inletPressure},
                       {"PressureGradient", column.pressureGradient},
                       {"Velocity", column.velocity},
                       {"Length", column.length},
                       {"EnergyBalance", column.energyBalance},
                       {"InfluxTemperature", column.influxTemperature},
                       {"GeometryType", column.geometryType},
                       {"ChannelShape", column.channelShape},
                       {"InternalChannelDimension", column.internalChannelDimension},
                       {"InternalDiameter", column.internalDiameter},
                       {"OuterDiameter", column.outerDiameter},
                       {"NumberOfChannels", column.numberOfChannels},
                       {"WashcoatThickness", column.washcoatThickness},
                       {"WallDensity", column.wallDensity},
                       {"GasThermalConductivity", column.gasThermalConductivity},
                       {"WallThermalConductivity", column.wallThermalConductivity},
                       {"HeatTransferGasSolid", column.heatTransferGasSolid},
                       {"HeatTransferGasWall", column.heatTransferGasWall},
                       {"HeatTransferWallExternal", column.heatTransferWallExternal},
                       {"HeatCapacityGas", column.heatCapacityGas},
                       {"HeatCapacitySolid", column.heatCapacitySolid},
                       {"HeatCapacityWall", column.heatCapacityWall},
                       {"Beds", beds},
                       {"Settings", column.settings}, {"GeometrySettings", column.geometrySettings}};
    addOptionalNumber(object, "OutletPressure", column.outletPressure);
    addOptionalNumber(object, "WashcoatVolumePerChannelVolume", column.washcoatVolumePerChannelVolume);
    columnArray.append(object);
  }

  QJsonArray simulationArray;
  for (const SimulationConfig& simulation : simulations_)
  {
    QJsonArray feeds;
    for (const FeedConfig& feed : simulation.componentFeeds)
    {
      feeds.append(QJsonObject{{"ComponentId", feed.componentId}, {"GasPhaseMolFraction", feed.gasFraction}});
    }
    QJsonObject simulationObject{{"Id", simulation.id},
                                 {"Type", simulation.type},
                                 {"DisplayName", simulation.displayName},
                                 {"MixturePredictionMethod", simulation.mixtureMethod},
                                 {"PressureStart", simulation.pressureStart},
                                 {"PressureEnd", simulation.pressureEnd},
                                 {"PressurePoints", simulation.pressurePoints},
                                 {"PressureScale", simulation.pressureScale},
                                 {"Integrator", simulation.integrator},
                                 {"AutoTime", simulation.autoTime},
                                 {"SimulationTime", simulation.simulationTime},
                                 {"WriteInterval", simulation.writeInterval},
                                 {"InitSteps", simulation.initSteps},
                                 {"TimeStep", simulation.timeStep},
                                 {"GridPoints", simulation.gridPoints},
                                 {"SelectedColumnId", simulation.selectedColumnId},
                                 {"PHComponentId", simulation.pHComponentId},
                                 {"Settings", simulation.settings},
                                 {"OriginalInput", simulation.originalInput},
                                 {"OriginalModel", simulation.originalModel},
                                 {"SwingAdsorptionPhases", swingPhasesStateJson(simulation.swingPhases)},
                                 {"ComponentFeeds", feeds}};
    if (!simulation.lastRun.id.isEmpty())
    {
      simulationObject.insert("LastRun", runRecordStateJson(simulation.lastRun, simulation.lastRunLog));
    }
    simulationArray.append(simulationObject);
  }

  return QJsonObject{{"Format", stateFormat},
                     {"Version", 4},
                     {"SelectedComponentIndex", selectedComponentIndex_},
                     {"SelectedColumnIndex", selectedColumnIndex_},
                     {"SelectedReactionIndex", selectedReactionIndex_},
                     {"SelectedSimulationIndex", selectedSimulationIndex_},
                     {"Components", componentArray},
                     {"Columns", columnArray},
                     {"Reactions", reactionsStateJson(reactions_, components_)},
                     {"Simulations", simulationArray}};
}

bool MainWindow::loadStateJson(const QJsonObject& state, QString* error)
{
  const QJsonArray componentArray = state.value("Components").toArray();
  const QJsonArray columnArray = state.value("Columns").toArray();
  const QJsonArray reactionArray = state.value("Reactions").toArray();
  const QJsonArray simulationArray = state.value("Simulations").toArray();
  if (componentArray.isEmpty())
  {
    if (error != nullptr) *error = "The saved state does not contain any components.";
    return false;
  }
  if (columnArray.isEmpty())
  {
    if (error != nullptr) *error = "The saved state does not contain any columns.";
    return false;
  }

  QList<ComponentConfig> newComponents;
  for (qsizetype i = 0; i < componentArray.size(); ++i)
  {
    if (!componentArray[i].isObject())
    {
      if (error != nullptr) *error = QString("Component %1 is not a JSON object.").arg(i + 1);
      return false;
    }
    const QJsonObject object = componentArray[i].toObject();
    ComponentConfig component;
    component.id = readString(object, "Id", QString("component-%1").arg(i + 1));
    component.name = readString(object, "Name", QString("Component %1").arg(i + 1));
    component.carrier = readBool(object, "CarrierGas", false);
    component.water = readBool(object, "Water", false);
    component.massTransfer = readOptionalNumber(object, "MassTransferCoefficient");
    component.axialDispersion = readOptionalNumber(object, "AxialDispersionCoefficient");
    component.molecularWeight = readOptionalNumber(object, "MolecularWeight");
    component.heatOfAdsorption = readOptionalNumber(object, "HeatOfAdsorption");
    component.referenceTemperature = readOptionalNumber(object, "ReferenceTemperature");
    const bool nonIsothermal = readBool(object, "NonIsothermal", false);
    component.settings = object.value("Settings").toObject();

    const QJsonArray physisorptionSites = object.value("PhysisorptionSites").toArray();
    for (qsizetype siteIndex = 0; siteIndex < physisorptionSites.size(); ++siteIndex)
    {
      if (physisorptionSites[siteIndex].isObject())
      {
        component.physisorptionSites.push_back(phySiteFromObject(physisorptionSites[siteIndex].toObject(),
                                                                 QString("phys-%1-%2").arg(i + 1).arg(siteIndex + 1)));
      }
    }
    const QJsonArray chemisorptionSites = object.value("ChemisorptionSites").toArray();
    for (qsizetype siteIndex = 0; siteIndex < chemisorptionSites.size(); ++siteIndex)
    {
      if (chemisorptionSites[siteIndex].isObject())
      {
        component.chemisorptionSites.push_back(chemSiteFromObject(
            chemisorptionSites[siteIndex].toObject(), QString("chem-%1-%2").arg(i + 1).arg(siteIndex + 1),
            QString("chemiso-%1-%2").arg(i + 1).arg(siteIndex + 1)));
      }
    }

    // Saved before temperature laws: nonIsothermal scaled the affinity with the heat of adsorption. The affinity is now
    // a VantHoff law that inherits that heat.
    if (nonIsothermal)
    {
      QList<PhysisorptionSiteConfig*> isotherms;
      for (PhysisorptionSiteConfig& site : component.physisorptionSites) isotherms.push_back(&site);
      for (ChemisorptionSiteConfig& site : component.chemisorptionSites) isotherms.push_back(&site.isotherm);
      for (PhysisorptionSiteConfig* site : isotherms)
      {
        const QList<ParameterDefinition> definitions = phyDefs(site->type);
        if (definitions.isEmpty() || parameterNumber(site->parameters, definitions.front()) == 0.0) continue;
        for (const ParameterDefinition& definition : definitions)
        {
          const double affinity = parameterNumber(site->parameters, definition);
          if (definition.heatLinked && affinity != 0.0 &&
              site->parameters.value(definition.key).typeId() != QMetaType::QVariantMap)
          {
            site->parameters.insert(definition.key, QVariantMap{{"Value", affinity}, {"Law", "VantHoff"}});
          }
        }
      }
    }
    newComponents.push_back(component);
  }

  QList<ColumnConfig> newColumns;
  for (qsizetype i = 0; i < columnArray.size(); ++i)
  {
    if (!columnArray[i].isObject())
    {
      if (error != nullptr) *error = QString("Column %1 is not a JSON object.").arg(i + 1);
      return false;
    }
    const QJsonObject object = columnArray[i].toObject();
    ColumnConfig column;
    column.id = readString(object, "Id", QString("column-%1").arg(i + 1));
    column.name = readString(object, "Name", QString("Column %1").arg(i + 1));
    column.mode = readString(object, "Mode", "single");
    column.settings = object.value("Settings").toObject();
    column.geometrySettings = object.value("GeometrySettings").toObject();
    column.boundaryCondition = readString(object, "BoundaryCondition", "InletPressureInletVelocity");
    column.temperature = readNumber(object, "Temperature", column.temperature);
    column.dynamicViscosity = readNumber(object, "DynamicViscosity", column.dynamicViscosity);
    column.voidFraction = readNumber(object, "VoidFraction", column.voidFraction);
    column.particleDensity = readNumber(object, "ParticleDensity", column.particleDensity);
    column.particleDiameter = readNumber(object, "ParticleDiameter", column.particleDiameter);
    column.inletPressure = readNumber(object, "InletPressure", column.inletPressure);
    column.outletPressure = readOptionalNumber(object, "OutletPressure");
    column.pressureGradient = readNumber(object, "PressureGradient", column.pressureGradient);
    column.velocity = readNumber(object, "Velocity", column.velocity);
    column.length = readNumber(object, "Length", column.length);
    column.energyBalance = readBool(object, "EnergyBalance", false);
    column.influxTemperature = readNumber(object, "InfluxTemperature", column.influxTemperature);
    const QString geometryType = readString(object, "GeometryType", column.geometryType);
    if (geometryType.compare("PackedBed", Qt::CaseInsensitive) != 0 &&
        geometryType.compare("Monolith", Qt::CaseInsensitive) != 0)
    {
      if (error != nullptr) *error = QString("Column %1 GeometryType must be PackedBed or Monolith.").arg(i + 1);
      return false;
    }
    column.geometryType = geometryType.compare("Monolith", Qt::CaseInsensitive) == 0 ? "Monolith" : "PackedBed";
    column.channelShape = readString(object, "ChannelShape", column.channelShape).toLower();
    if (column.channelShape != "triangular" && column.channelShape != "square" && column.channelShape != "hexagonal" &&
        column.channelShape != "circular")
    {
      column.channelShape = "square";
    }
    column.internalChannelDimension = readNumber(object, "InternalChannelDimension", column.internalChannelDimension);
    column.internalDiameter = readNumber(object, "InternalDiameter", column.internalDiameter);
    column.outerDiameter = readNumber(object, "OuterDiameter", column.outerDiameter);
    column.numberOfChannels = std::max(1, readInt(object, "NumberOfChannels", column.numberOfChannels));
    column.washcoatThickness = readNumber(object, "WashcoatThickness", column.washcoatThickness);
    column.washcoatVolumePerChannelVolume = readOptionalNumber(object, "WashcoatVolumePerChannelVolume");
    column.wallDensity = readNumber(object, "WallDensity", column.wallDensity);
    column.gasThermalConductivity = readNumber(object, "GasThermalConductivity", column.gasThermalConductivity);
    column.wallThermalConductivity = readNumber(object, "WallThermalConductivity", column.wallThermalConductivity);
    column.heatTransferGasSolid = readNumber(object, "HeatTransferGasSolid", column.heatTransferGasSolid);
    column.heatTransferGasWall = readNumber(object, "HeatTransferGasWall", column.heatTransferGasWall);
    column.heatTransferWallExternal = readNumber(object, "HeatTransferWallExternal", column.heatTransferWallExternal);
    column.heatCapacityGas = readNumber(object, "HeatCapacityGas", column.heatCapacityGas);
    column.heatCapacitySolid = readNumber(object, "HeatCapacitySolid", column.heatCapacitySolid);
    column.heatCapacityWall = readNumber(object, "HeatCapacityWall", column.heatCapacityWall);

    const QJsonArray bedArray = object.value("Beds").toArray();
    for (qsizetype bedIndex = 0; bedIndex < bedArray.size(); ++bedIndex)
    {
      if (!bedArray[bedIndex].isObject())
      {
        continue;
      }
      const QJsonObject bedObject = bedArray[bedIndex].toObject();
      BedConfig bed;
      bed.id = readString(bedObject, "Id", QString("bed-%1-%2").arg(i + 1).arg(bedIndex + 1));
      bed.name = readString(bedObject, "Name", QString("Bed %1").arg(bedIndex + 1));
      bed.settings = bedObject.value("Settings").toObject();
      bed.length = readNumber(bedObject, "Length", bed.length);
      bed.gridPoints = readInt(bedObject, "GridPoints", bed.gridPoints);
      bed.mixFraction = readNumber(bedObject, "MixFraction", bed.mixFraction);
      bed.voidFraction = readNumber(bedObject, "VoidFraction", bed.voidFraction);
      bed.particleDensity = readNumber(bedObject, "ParticleDensity", bed.particleDensity);
      bed.particleDiameter = readNumber(bedObject, "ParticleDiameter", bed.particleDiameter);
      bed.interfaceLengthAfter = readNumber(bedObject, "InterfaceLengthAfter", bed.interfaceLengthAfter);
      bed.interfaceGridPointsAfter = readInt(bedObject, "InterfaceGridPointsAfter", bed.interfaceGridPointsAfter);
      column.beds.push_back(bed);
    }
    if (column.beds.isEmpty())
    {
      BedConfig bed;
      bed.id = QString("bed-%1-1").arg(i + 1);
      bed.name = "Bed 1";
      column.beds.push_back(bed);
    }
    if (column.mode != "single" && column.mode != "multibed" && column.mode != "mixed")
    {
      column.mode = column.beds.size() > 1 ? "multibed" : "single";
    }
    newColumns.push_back(column);
  }

  QList<ReactionConfig> newReactions;
  for (qsizetype reactionIndex = 0; reactionIndex < reactionArray.size(); ++reactionIndex)
  {
    if (reactionArray[reactionIndex].isObject())
    {
      newReactions.push_back(reactionFromObject(reactionArray[reactionIndex].toObject(),
                                                QString("reaction-%1").arg(reactionIndex + 1),
                                                QString("Reaction %1").arg(reactionIndex + 1), newComponents));
    }
  }

  QList<SimulationConfig> newSimulations;
  for (qsizetype i = 0; i < simulationArray.size(); ++i)
  {
    if (!simulationArray[i].isObject())
    {
      continue;
    }
    const QJsonObject object = simulationArray[i].toObject();
    SimulationConfig simulation;
    simulation.id = readString(object, "Id", QString("simulation-%1").arg(i + 1));
    simulation.type = readString(object, "Type", "Breakthrough");
    simulation.settings = object.value("Settings").toObject();
    simulation.pHComponentId = readString(object, "PHComponentId", {});
    simulation.originalInput = object.value("OriginalInput").toObject();
    simulation.originalModel = object.value("OriginalModel").toObject();
    const QJsonObject supportCheck{{"SimulationType", simulation.type},
                                  {"MixturePredictionMethod", object.value("MixturePredictionMethod")}};
    if (!LabInput::supported(supportCheck, error) ||
        (!simulation.originalInput.isEmpty() && !LabInput::supported(simulation.originalInput, error))) return false;
    if (simulation.type != "Breakthrough" && simulation.type != "SwingAdsorption" &&
        simulation.type != "MixturePrediction")
    {
      simulation.type = "Breakthrough";
    }
    simulation.displayName = readString(object, "DisplayName", QString("Simulation %1").arg(i + 1));
    simulation.mixtureMethod = readString(object, "MixturePredictionMethod", simulation.mixtureMethod);
    simulation.pressureStart = readNumber(object, "PressureStart", simulation.pressureStart);
    simulation.pressureEnd = readNumber(object, "PressureEnd", simulation.pressureEnd);
    simulation.pressurePoints = std::max(1, readInt(object, "PressurePoints", simulation.pressurePoints));
    simulation.pressureScale = readString(object, "PressureScale", simulation.pressureScale);
    simulation.integrator = readString(object, "Integrator", simulation.integrator);
    simulation.autoTime = readBool(object, "AutoTime", simulation.autoTime);
    simulation.simulationTime = readNumber(object, "SimulationTime", simulation.simulationTime);
    simulation.writeInterval = readNumber(object, "WriteInterval", simulation.writeInterval);
    simulation.initSteps = readInt(object, "InitSteps", simulation.initSteps);
    simulation.timeStep = readNumber(object, "TimeStep", simulation.timeStep);
    simulation.gridPoints = std::max(1, readInt(object, "GridPoints", simulation.gridPoints));
    simulation.selectedColumnId = readString(object, "SelectedColumnId", newColumns.front().id);
    if (!stringSetContains(newColumns, simulation.selectedColumnId))
    {
      simulation.selectedColumnId = newColumns.front().id;
    }
    const auto columnIt = std::find_if(newColumns.cbegin(), newColumns.cend(), [&simulation](const ColumnConfig& column)
                                       { return column.id == simulation.selectedColumnId; });
    const double fallbackTemperature = columnIt == newColumns.cend() ? 300.0 : columnIt->temperature;
    const double fallbackInletPressure = columnIt == newColumns.cend() ? 100000.0 : columnIt->inletPressure;

    const QJsonArray phaseArray = object.value("SwingAdsorptionPhases").toArray();
    for (qsizetype phaseIndex = 0; phaseIndex < phaseArray.size(); ++phaseIndex)
    {
      if (phaseArray[phaseIndex].isObject())
      {
        simulation.swingPhases.push_back(swingPhaseFromObject(
            phaseArray[phaseIndex].toObject(), QString("phase-%1-%2").arg(i + 1).arg(phaseIndex + 1),
            fallbackTemperature, fallbackInletPressure, SwingPhaseConfig{}.duration));
      }
    }
    if (simulation.type == "SwingAdsorption" && simulation.swingPhases.isEmpty())
    {
      simulation.swingPhases.push_back(SwingPhaseConfig{QString("phase-%1-1").arg(i + 1), "Adsorption",
                                                        fallbackTemperature, fallbackInletPressure});
      simulation.swingPhases.push_back(SwingPhaseConfig{QString("phase-%1-2").arg(i + 1), "Regeneration",
                                                        fallbackTemperature,
                                                        std::max(1000.0, 0.1 * fallbackInletPressure)});
    }

    const QJsonValue lastRunValue = object.value("LastRun");
    if (lastRunValue.isObject())
    {
      const QJsonObject lastRunObject = lastRunValue.toObject();
      simulation.lastRun = runRecordFromObject(lastRunObject);
      simulation.lastRunLog = readString(lastRunObject, "LogTail", {});
    }

    const QJsonArray feedArray = object.value("ComponentFeeds").toArray();
    for (const QJsonValue& feedValue : feedArray)
    {
      if (!feedValue.isObject())
      {
        continue;
      }
      const QJsonObject feedObject = feedValue.toObject();
      const QString componentId = readString(feedObject, "ComponentId", {});
      if (!componentId.isEmpty() && stringSetContains(newComponents, componentId))
      {
        simulation.componentFeeds.push_back({componentId, readNumber(feedObject, "GasPhaseMolFraction", 0.0)});
      }
    }
    newSimulations.push_back(simulation);
  }
  if (newSimulations.isEmpty())
  {
    SimulationConfig simulation;
    simulation.id = "simulation-1";
    simulation.displayName = "Simulation 1";
    simulation.selectedColumnId = newColumns.front().id;
    newSimulations.push_back(simulation);
  }

  components_ = newComponents;
  columns_ = newColumns;
  reactions_ = newReactions;
  simulations_ = newSimulations;
  for (SimulationConfig& simulation : simulations_)
  {
    if (!simulation.lastRun.id.isEmpty())
    {
      runner_.rememberRun(simulation.lastRun);
      if (simulation.lastRunLog.isEmpty())
      {
        simulation.lastRunLog = runner_.logTail(simulation.lastRun);
      }
    }
  }
  selectedComponentIndex_ =
      std::clamp(readInt(state, "SelectedComponentIndex", 0), 0, static_cast<int>(components_.size()) - 1);
  selectedColumnIndex_ = std::clamp(readInt(state, "SelectedColumnIndex", 0), 0, static_cast<int>(columns_.size()) - 1);
  selectedReactionIndex_ =
      std::clamp(readInt(state, "SelectedReactionIndex", 0), 0, std::max(0, static_cast<int>(reactions_.size()) - 1));
  selectedSimulationIndex_ =
      std::clamp(readInt(state, "SelectedSimulationIndex", 0), 0, static_cast<int>(simulations_.size()) - 1);
  resetIdCounterFromState();
  refreshAll();
  return true;
}

bool MainWindow::importSimulationJson(const QJsonObject& input, QString* error, const QString& sourceDirectory)
{
  QJsonObject simulationJson;
  try
  {
    simulationJson = LabInput::resolveFiles(LabInput::canonicalize(input), sourceDirectory);
  }
  catch (const std::exception& exception)
  {
    if (error) *error = QString::fromUtf8(exception.what());
    return false;
  }
  const auto parsed = LabInput::read(simulationJson, error);
  if (!parsed) return false;
  const QJsonArray componentArray = simulationJson.value("Components").toArray();
  if (componentArray.isEmpty())
  {
    if (error != nullptr) *error = "The simulation JSON does not contain a Components array.";
    return false;
  }

  QList<ComponentConfig> newComponents;
  QList<FeedConfig> feeds;
  for (qsizetype i = 0; i < componentArray.size(); ++i)
  {
    if (!componentArray[i].isObject())
    {
      if (error != nullptr) *error = QString("Component %1 is not a JSON object.").arg(i + 1);
      return false;
    }
    const QJsonObject object = componentArray[i].toObject();
    ComponentConfig component;
    component.id = QString("component-%1").arg(i + 1);
    component.name = readString(object, "Name", QString("Component %1").arg(i + 1));
    component.carrier = readBool(object, "CarrierGas", false);
    component.water = readBool(object, "Water", false);
    component.massTransfer = readOptionalNumber(object, "MassTransferCoefficient");
    component.axialDispersion = readOptionalNumber(object, "AxialDispersionCoefficient");
    component.molecularWeight = readOptionalNumber(object, "MolecularWeight");
    component.heatOfAdsorption = readOptionalNumber(object, "HeatOfAdsorption");
    component.referenceTemperature = readOptionalNumber(object, "referenceTemperature");
    if (!component.referenceTemperature)
    {
      component.referenceTemperature = readOptionalNumber(object, "ReferenceTemperature");
    }
    component.settings = selectSettings(object, componentSettingKeys);

    const QJsonArray physisorptionSites = object.value("PhysisorptionSites").toArray();
    for (qsizetype siteIndex = 0; siteIndex < physisorptionSites.size(); ++siteIndex)
    {
      if (physisorptionSites[siteIndex].isObject())
      {
        component.physisorptionSites.push_back(phySiteFromObject(physisorptionSites[siteIndex].toObject(),
                                                                 QString("phys-%1-%2").arg(i + 1).arg(siteIndex + 1)));
      }
    }

    const QJsonArray chemisorptionSites = object.value("ChemisorptionSites").toArray();
    for (qsizetype siteIndex = 0; siteIndex < chemisorptionSites.size(); ++siteIndex)
    {
      if (chemisorptionSites[siteIndex].isObject())
      {
        component.chemisorptionSites.push_back(chemSiteFromObject(
            chemisorptionSites[siteIndex].toObject(), QString("chem-%1-%2").arg(i + 1).arg(siteIndex + 1),
            QString("chemiso-%1-%2").arg(i + 1).arg(siteIndex + 1),
            component.massTransfer.value_or(0.0), component.heatOfAdsorption.value_or(0.0)));
      }
    }

    newComponents.push_back(component);
    feeds.push_back({component.id, readNumber(object, "GasPhaseMolFraction", 0.0)});
  }

  ColumnConfig column;
  column.id = "column-1";
  column.temperature = parsed->temperature;
  column.dynamicViscosity = parsed->dynamicViscosity;
  column.voidFraction = parsed->columnVoidFraction;
  column.particleDensity = parsed->particleDensity;
  column.particleDiameter = parsed->particleDiameter;
  column.inletPressure = parsed->inletPressure;
  column.pressureGradient = parsed->pressureGradient;
  column.velocity = parsed->columnEntranceVelocity;
  column.length = parsed->columnLength;
  column.influxTemperature = parsed->influxTemperature;
  column.internalDiameter = parsed->internalDiameter;
  column.outerDiameter = parsed->outerDiameter;
  column.wallDensity = parsed->wallDensity;
  column.gasThermalConductivity = parsed->gasThermalConductivity;
  column.wallThermalConductivity = parsed->wallThermalConductivity;
  column.heatTransferGasSolid = parsed->heatTransferGasSolid;
  column.heatTransferGasWall = parsed->heatTransferGasWall;
  column.heatTransferWallExternal = parsed->heatTransferWallExternal;
  column.heatCapacityGas = parsed->heatCapacityGas;
  column.heatCapacitySolid = parsed->heatCapacitySolid;
  column.heatCapacityWall = parsed->heatCapacityWall;
  column.name = readString(simulationJson, "ColumnName", readString(simulationJson, "DisplayName", "Imported column"));
  column.mode = "single";
  column.settings = selectSettings(simulationJson, columnSettingKeys);
  column.boundaryCondition = choice(readString(simulationJson, "BoundaryCondition", column.boundaryCondition),
      {"InletPressureInletVelocity", "InletPressureOutletPressure", "InletVelocityOutletPressure", "FixedVelocity",
       "FixedPressureInletVelocity"});
  column.temperature = readNumber(simulationJson, "Temperature", column.temperature);
  column.dynamicViscosity = readNumber(simulationJson, "DynamicViscosity", column.dynamicViscosity);
  column.voidFraction = readNumber(simulationJson, "ColumnVoidFraction", column.voidFraction);
  column.particleDensity = readNumber(simulationJson, "ParticleDensity", column.particleDensity);
  column.particleDiameter = readNumber(simulationJson, "ParticleDiameter", column.particleDiameter);
  column.inletPressure = readNumber(simulationJson, "InletPressure", column.inletPressure);
  column.outletPressure = readOptionalNumber(simulationJson, "OutletPressure");
  column.pressureGradient = readNumber(simulationJson, "PressureGradient", column.pressureGradient);
  column.velocity = readNumber(simulationJson, "ColumnEntranceVelocity", column.velocity);
  column.length = readNumber(simulationJson, "ColumnLength", column.length);
  column.energyBalance = readBool(simulationJson, "energyBalance", readBool(simulationJson, "EnergyBalance", false));
  column.influxTemperature = readNumber(simulationJson, "InfluxTemperature", column.influxTemperature);
  const QJsonObject geometry = simulationJson.value("Geometry").toObject();
  column.geometrySettings =
      selectSettings(geometry, {"ForchheimerCoefficient", "AxialDispersivity", "TaylorArisDispersion"});
  if (!geometry.isEmpty())
  {
    const QString geometryType = readString(geometry, "Type", column.geometryType);
    if (geometryType.compare("PackedBed", Qt::CaseInsensitive) != 0 &&
        geometryType.compare("Monolith", Qt::CaseInsensitive) != 0)
    {
      if (error != nullptr) *error = "Geometry Type must be PackedBed or Monolith.";
      return false;
    }
    column.geometryType = geometryType.compare("Monolith", Qt::CaseInsensitive) == 0 ? "Monolith" : "PackedBed";
    if (column.geometryType == "Monolith")
    {
      column.channelShape = readString(geometry, "ChannelShape", column.channelShape).toLower();
      column.internalChannelDimension =
          readNumber(geometry, "InternalChannelDimension", column.internalChannelDimension);
      column.outerDiameter = readNumber(geometry, "OuterDiameter", column.outerDiameter);
      column.numberOfChannels = std::max(1, readInt(geometry, "NumberOfChannels", column.numberOfChannels));
      column.washcoatThickness = readNumber(geometry, "WashcoatThickness", column.washcoatThickness);
      column.washcoatVolumePerChannelVolume = readOptionalNumber(geometry, "WashcoatVolumePerChannelVolume");
    }
    else
    {
      column.geometryType = "PackedBed";
      column.voidFraction = readNumber(geometry, "ColumnVoidFraction", column.voidFraction);
      column.particleDiameter = readNumber(geometry, "ParticleDiameter", column.particleDiameter);
      column.internalDiameter = readNumber(geometry, "InternalDiameter", column.internalDiameter);
      column.outerDiameter = readNumber(geometry, "OuterDiameter", column.outerDiameter);
    }
  }
  column.wallDensity = readNumber(simulationJson, "wallDensity", column.wallDensity);
  column.gasThermalConductivity = readNumber(simulationJson, "gasThermalConductivity", column.gasThermalConductivity);
  column.wallThermalConductivity =
      readNumber(simulationJson, "wallThermalConductivity", column.wallThermalConductivity);
  column.heatTransferGasSolid = readNumber(simulationJson, "heatTransferGasSolid", column.heatTransferGasSolid);
  column.heatTransferGasWall = readNumber(simulationJson, "heatTransferGasWall", column.heatTransferGasWall);
  column.heatTransferWallExternal =
      readNumber(simulationJson, "heatTransferWallExternal", column.heatTransferWallExternal);
  column.heatCapacityGas = readNumber(simulationJson, "heatCapacityGas", column.heatCapacityGas);
  column.heatCapacitySolid = readNumber(simulationJson, "heatCapacitySolid", column.heatCapacitySolid);
  column.heatCapacityWall = readNumber(simulationJson, "heatCapacityWall", column.heatCapacityWall);

  const QJsonArray adsorbents = simulationJson.value("Adsorbents").toArray();
  const QJsonArray sections = simulationJson.value("ColumnSections").toArray();
  const bool isUniformMix =
      !adsorbents.isEmpty() && sections.isEmpty() &&
      std::all_of(adsorbents.cbegin(), adsorbents.cend(), [](const QJsonValue& value)
                  { return value.isObject() && value.toObject().contains("MixFraction"); });
  if (!adsorbents.isEmpty() && (!sections.isEmpty() || isUniformMix))
  {
    column.beds.clear();
    column.mode = isUniformMix ? "mixed" : "multibed";
    if (isUniformMix)
    {
      for (const QJsonValue& adsorbentValue : adsorbents)
      {
        const QJsonObject adsorbent = adsorbentValue.toObject();
        BedConfig bed;
        bed.id = QString("bed-1-%1").arg(column.beds.size() + 1);
        bed.name = readString(adsorbent, "Name", QString("Bed %1").arg(column.beds.size() + 1));
        bed.mixFraction = readNumber(adsorbent, "MixFraction", bed.mixFraction);
        bed.voidFraction = readNumber(adsorbent, "ColumnVoidFraction", column.voidFraction);
        bed.particleDensity = readNumber(adsorbent, "ParticleDensity", column.particleDensity);
        bed.settings = selectSettings(adsorbent, {"ComponentParameters", "Components", "AxialDispersivity"});
        bed.particleDiameter = readNumber(adsorbent, "ParticleDiameter", column.particleDiameter);
        column.beds.push_back(bed);
      }
    }
    else
    {
      for (const QJsonValue& sectionValue : sections)
      {
        if (!sectionValue.isObject())
        {
          continue;
        }
        const QJsonObject section = sectionValue.toObject();
        if (section.contains("InterfaceLength"))
        {
          if (!column.beds.isEmpty())
          {
            column.beds.last().interfaceLengthAfter = readNumber(section, "InterfaceLength", 0.0);
            column.beds.last().interfaceGridPointsAfter = readInt(section, "NumberOfGridPoints", 0);
          }
          continue;
        }
        const QString adsorbentName =
            readString(section, "Adsorbent", QString("Bed %1").arg(column.beds.size() + 1));
        QJsonObject adsorbent;
        for (const QJsonValue& adsorbentValue : adsorbents)
        {
          if (adsorbentValue.isObject() && readString(adsorbentValue.toObject(), "Name", {}) == adsorbentName)
          {
            adsorbent = adsorbentValue.toObject();
            break;
          }
        }
        BedConfig bed;
        bed.id = QString("bed-1-%1").arg(column.beds.size() + 1);
        bed.name = adsorbentName;
        bed.length = readNumber(section, "Length", column.length);
        bed.gridPoints = readInt(section, "NumberOfGridPoints", 0);
        bed.voidFraction = readNumber(adsorbent, "ColumnVoidFraction", column.voidFraction);
        bed.particleDensity = readNumber(adsorbent, "ParticleDensity", column.particleDensity);
        bed.settings = selectSettings(adsorbent, {"ComponentParameters", "Components", "AxialDispersivity"});
        bed.particleDiameter = readNumber(adsorbent, "ParticleDiameter", column.particleDiameter);
        column.beds.push_back(bed);
      }
    }

  }

  if (!adsorbents.isEmpty() && sections.isEmpty() && !isUniformMix)
  {
    column.mode = "multibed";
    for (const auto& value : adsorbents)
    {
      const auto adsorbent = value.toObject();
      BedConfig bed;
      bed.id = QString("bed-1-%1").arg(column.beds.size() + 1);
      bed.name = readString(adsorbent, "Name", QString::number(column.beds.size()));
      bed.length = readNumber(adsorbent, "AdsorbentLength", column.length);
      bed.voidFraction = readNumber(adsorbent, "ColumnVoidFraction", column.voidFraction);
      bed.particleDensity = readNumber(adsorbent, "ParticleDensity", column.particleDensity);
      bed.particleDiameter = readNumber(adsorbent, "ParticleDiameter", column.particleDiameter);
      bed.settings = selectSettings(adsorbent, {"ComponentParameters", "Components", "AxialDispersivity"});
      column.beds.push_back(bed);
    }
  }

  SimulationConfig simulation;
  simulation.id = "simulation-1";
  const QString type = readString(simulationJson, "SimulationType", "Breakthrough");
  simulation.type = choice(type, {"MixturePrediction", "SwingAdsorption", "Breakthrough"});
  if (simulation.type.compare("Test", Qt::CaseInsensitive) == 0) simulation.type = "Breakthrough";
  simulation.settings = selectSettings(simulationJson, simulationSettingKeys);
  if (simulation.type != "MixturePrediction" && geometry.isEmpty())
  {
    if (error != nullptr) *error = "Breakthrough and SwingAdsorption simulation JSON must contain Geometry.";
    return false;
  }
  simulation.displayName = readString(simulationJson, "DisplayName", "Imported simulation");
  simulation.mixtureMethod = choice(readString(simulationJson, "MixturePredictionMethod", "IAST"),
                                    {"IAST", "SIAST", "EI", "SEI", "SCI", "SPI"});
  simulation.pressureStart = readNumber(simulationJson, "PressureStart", simulation.pressureStart);
  simulation.pressureEnd = readNumber(simulationJson, "PressureEnd", simulation.pressureEnd);
  simulation.pressurePoints = std::max(1, readInt(simulationJson, "NumberOfPressurePoints", simulation.pressurePoints));
  simulation.pressureScale = readString(simulationJson, "PressureScale", simulation.pressureScale);
  simulation.integrator = choice(readString(simulationJson, "BreakthroughIntegrator", "RungeKutta3"),
                                 {"CVODE", "RungeKutta3", "ROS34PW2"});
  if (simulation.integrator != "CVODE" && simulation.integrator != "RungeKutta3" && simulation.integrator != "ROS34PW2")
  {
    if (error != nullptr) *error = "BreakthroughIntegrator must be CVODE, RungeKutta3, or ROS34PW2.";
    return false;
  }
  const QJsonValue simulationTime = simulationJson.value("SimulationTime");
  simulation.autoTime = !simulationTime.isDouble();
  if (simulationTime.isDouble()) simulation.simulationTime = simulationTime.toDouble();
  simulation.writeInterval = readNumber(simulationJson, "WriteInterval", simulation.writeInterval);
  simulation.initSteps = readInt(simulationJson, "NumberOfInitTimeSteps", simulation.initSteps);
  simulation.timeStep = readNumber(simulationJson, "TimeStep", simulation.timeStep);
  simulation.gridPoints = std::max(1, readInt(simulationJson, "NumberOfGridPoints", simulation.gridPoints));
  simulation.selectedColumnId = column.id;
  simulation.componentFeeds = feeds;
  const int pHIndex = simulationJson.value("pHComponent").toInt(-1);
  if (pHIndex >= 0 && pHIndex < newComponents.size()) simulation.pHComponentId = newComponents[pHIndex].id;

  QList<ReactionConfig> newReactions;
  const QJsonArray reactionArray = simulationJson.value("Reactions").toArray();
  for (qsizetype reactionIndex = 0; reactionIndex < reactionArray.size(); ++reactionIndex)
  {
    if (reactionArray[reactionIndex].isObject())
    {
      auto reaction = reactionFromObject(reactionArray[reactionIndex].toObject(),
                                                QString("reaction-1-%1").arg(reactionIndex + 1),
                                                QString("Reaction %1").arg(reactionIndex + 1), newComponents);
      const auto& resolved = parsed->reactions[static_cast<size_t>(reactionIndex)];
      for (qsizetype i = 0; i < reaction.reactants.size(); ++i)
        reaction.reactants[i].componentId = newComponents[static_cast<qsizetype>(resolved.reactants[static_cast<size_t>(i)])].id;
      for (qsizetype i = 0; i < reaction.products.size(); ++i)
        reaction.products[i].componentId = newComponents[static_cast<qsizetype>(resolved.products[static_cast<size_t>(i)])].id;
      newReactions.push_back(reaction);
    }
  }
  if (simulation.type == "SwingAdsorption")
  {
    QJsonArray phaseArray = simulationJson.value("SwingAdsorptionPhases").toArray();
    for (qsizetype phaseIndex = 0; phaseIndex < phaseArray.size(); ++phaseIndex)
    {
      if (phaseArray[phaseIndex].isObject())
      {
        simulation.swingPhases.push_back(
            swingPhaseFromObject(phaseArray[phaseIndex].toObject(), QString("phase-1-%1").arg(phaseIndex + 1),
                                 column.temperature, column.inletPressure, SwingPhaseConfig{}.duration));
      }
    }

    if (simulation.swingPhases.isEmpty())
    {
      if (error != nullptr) *error = "SwingAdsorption simulation JSON must contain SwingAdsorptionPhases.";
      return false;
    }
  }

  components_ = newComponents;
  columns_ = {column};
  reactions_ = newReactions;
  simulations_ = {simulation};
  selectedComponentIndex_ = 0;
  selectedColumnIndex_ = 0;
  selectedReactionIndex_ = 0;
  selectedSimulationIndex_ = 0;
  simulations_[0].originalModel = buildModeledSimulationJson();
  simulations_[0].originalInput = simulationJson;
  resetIdCounterFromState();
  refreshAll();
  return true;
}

void MainWindow::resetIdCounterFromState()
{
  int maxId = 10;
  for (const ComponentConfig& component : components_)
  {
    maxId = std::max(maxId, numericIdSuffix(component.id));
    for (const PhysisorptionSiteConfig& site : component.physisorptionSites)
    {
      maxId = std::max(maxId, numericIdSuffix(site.id));
    }
    for (const ChemisorptionSiteConfig& site : component.chemisorptionSites)
    {
      maxId = std::max(maxId, numericIdSuffix(site.id));
      maxId = std::max(maxId, numericIdSuffix(site.isotherm.id));
    }
  }
  for (const ColumnConfig& column : columns_)
  {
    maxId = std::max(maxId, numericIdSuffix(column.id));
    for (const BedConfig& bed : column.beds)
    {
      maxId = std::max(maxId, numericIdSuffix(bed.id));
    }
  }
  for (const SimulationConfig& simulation : simulations_)
  {
    maxId = std::max(maxId, numericIdSuffix(simulation.id));
    for (const SwingPhaseConfig& phase : simulation.swingPhases)
    {
      maxId = std::max(maxId, numericIdSuffix(phase.id));
    }
  }
  for (const ReactionConfig& reaction : reactions_)
  {
    maxId = std::max(maxId, numericIdSuffix(reaction.id));
  }
  idCounter_ = maxId;
}

void MainWindow::loadJsonFile()
{
  const QString path =
      QFileDialog::getOpenFileName(this, "Load JSON", QDir::currentPath(), "JSON files (*.json);;All files (*)");
  if (path.isEmpty())
  {
    return;
  }

  QFile file(path);
  if (!file.open(QIODevice::ReadOnly))
  {
    QMessageBox::warning(this, "Ruptura Lab", QString("Could not open %1.").arg(path));
    return;
  }

  QJsonParseError parseError;
  const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &parseError);
  if (parseError.error != QJsonParseError::NoError || !document.isObject())
  {
    QMessageBox::warning(this, "Ruptura Lab",
                         QString("Could not parse %1 as a JSON object: %2").arg(path, parseError.errorString()));
    return;
  }

  QString error;
  const QJsonObject root = document.object();
  const bool isState =
      readString(root, "Format", {}) == stateFormat || (root.contains("Simulations") && root.contains("Columns"));
  const bool ok = isState ? loadStateJson(root, &error) : importSimulationJson(root, &error, QFileInfo(path).absolutePath());
  if (!ok)
  {
    QMessageBox::warning(this, "Ruptura Lab", error);
    return;
  }

  statusBar()->showMessage(isState ? QString("Loaded Ruptura Lab state from %1").arg(path)
                                   : QString("Imported simulation JSON from %1").arg(path));
}

void MainWindow::saveStateFile()
{
  QString path =
      QFileDialog::getSaveFileName(this, "Save Ruptura Lab State", QDir::current().filePath("ruptura-lab-state.json"),
                                   "JSON files (*.json);;All files (*)");
  if (path.isEmpty())
  {
    return;
  }
  if (QFileInfo(path).suffix().isEmpty())
  {
    path += ".json";
  }

  QFile file(path);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate))
  {
    QMessageBox::warning(this, "Ruptura Lab", QString("Could not write %1.").arg(path));
    return;
  }
  file.write(QJsonDocument(buildStateJson()).toJson(QJsonDocument::Indented));
  statusBar()->showMessage(QString("Saved Ruptura Lab state to %1").arg(path));
}

void MainWindow::addComponent()
{
  ComponentConfig component;
  component.id = newId("component");
  component.name = QString("Component %1").arg(static_cast<int>(components_.size()) + 1);
  component.massTransfer = 0.06;
  component.axialDispersion = 0.0;
  component.physisorptionSites.push_back(makePhySite());
  components_.push_back(component);
  selectedComponentIndex_ = static_cast<int>(components_.size()) - 1;
  refreshAll();
}

void MainWindow::removeComponent()
{
  if (components_.isEmpty())
  {
    return;
  }
  const QString componentId = components_[selectedComponentIndex_].id;
  components_.removeAt(selectedComponentIndex_);
  for (SimulationConfig& simulation : simulations_)
  {
    for (int i = static_cast<int>(simulation.componentFeeds.size()) - 1; i >= 0; --i)
    {
      if (simulation.componentFeeds[i].componentId == componentId)
      {
        simulation.componentFeeds.removeAt(i);
      }
    }
  }
  if (components_.isEmpty())
  {
    ComponentConfig component;
    component.id = newId("component");
    component.name = "Component 1";
    component.massTransfer = 0.06;
    component.axialDispersion = 0.0;
    component.physisorptionSites.push_back(makePhySite());
    components_.push_back(component);
  }
  const QString fallbackComponentId = components_.front().id;
  auto cleanParticipants = [&](QList<ReactionParticipantConfig>& participants)
  {
    for (int i = static_cast<int>(participants.size()) - 1; i >= 0; --i)
    {
      if (participants[i].componentId == componentId)
      {
        participants.removeAt(i);
      }
    }
    if (participants.isEmpty())
    {
      participants.push_back(ReactionParticipantConfig{fallbackComponentId, 1.0});
    }
  };
  for (ReactionConfig& reaction : reactions_)
  {
    cleanParticipants(reaction.reactants);
    cleanParticipants(reaction.products);
  }
  selectedComponentIndex_ = std::clamp(selectedComponentIndex_, 0, static_cast<int>(components_.size()) - 1);
  refreshAll();
}

void MainWindow::addColumn()
{
  ColumnConfig column = makeColumn(QString("Column %1").arg(static_cast<int>(columns_.size()) + 1));
  columns_.push_back(column);
  selectedColumnIndex_ = static_cast<int>(columns_.size()) - 1;
  if (SimulationConfig* simulation = activeSimulation())
  {
    simulation->selectedColumnId = column.id;
  }
  refreshAll();
}

void MainWindow::removeColumn()
{
  if (columns_.isEmpty())
  {
    return;
  }
  const QString columnId = columns_[selectedColumnIndex_].id;
  columns_.removeAt(selectedColumnIndex_);
  if (columns_.isEmpty())
  {
    columns_.push_back(makeColumn("Column 1", "column-1"));
  }
  for (SimulationConfig& simulation : simulations_)
  {
    if (simulation.selectedColumnId == columnId)
    {
      simulation.selectedColumnId = columns_.front().id;
    }
  }
  selectedColumnIndex_ = std::clamp(selectedColumnIndex_, 0, static_cast<int>(columns_.size()) - 1);
  refreshAll();
}

void MainWindow::addReaction()
{
  reactions_.push_back(makeReaction(QString("Reaction %1").arg(static_cast<int>(reactions_.size()) + 1)));
  selectedReactionIndex_ = static_cast<int>(reactions_.size()) - 1;
  refreshAll();
}

void MainWindow::removeReaction()
{
  if (reactions_.isEmpty())
  {
    return;
  }
  const int row = reactionList_ == nullptr ? selectedReactionIndex_ : reactionList_->currentRow();
  const int index = row >= 0 ? row : selectedReactionIndex_;
  if (index < 0 || index >= static_cast<int>(reactions_.size()))
  {
    return;
  }
  reactions_.removeAt(index);
  selectedReactionIndex_ = std::clamp(index, 0, std::max(0, static_cast<int>(reactions_.size()) - 1));
  refreshAll();
}

void MainWindow::addSimulation()
{
  QList<FeedConfig> feeds;
  if (const SimulationConfig* simulation = activeSimulation())
  {
    feeds = simulation->componentFeeds;
  }
  SimulationConfig simulation =
      makeSimulation(QString("Simulation %1").arg(static_cast<int>(simulations_.size()) + 1), {}, feeds);
  if (const SimulationConfig* current = activeSimulation())
  {
    simulation.selectedColumnId = current->selectedColumnId;
  }
  simulations_.push_back(simulation);
  selectedSimulationIndex_ = static_cast<int>(simulations_.size()) - 1;
  refreshAll();
}

void MainWindow::copySimulation()
{
  if (const SimulationConfig* current = activeSimulation())
  {
    SimulationConfig copy = *current;
    copy.id = newId("simulation");
    for (SwingPhaseConfig& phase : copy.swingPhases)
    {
      phase.id = newId("phase");
    }
    copy.displayName += " copy";
    copy.lastRun = {};
    copy.lastRunLog.clear();
    simulations_.insert(selectedSimulationIndex_ + 1, copy);
    selectedSimulationIndex_ += 1;
    refreshAll();
  }
}

void MainWindow::removeSimulation()
{
  if (simulations_.isEmpty())
  {
    return;
  }
  simulations_.removeAt(selectedSimulationIndex_);
  if (simulations_.isEmpty())
  {
    simulations_.push_back(makeSimulation("Simulation 1", "simulation-1"));
  }
  selectedSimulationIndex_ = std::clamp(selectedSimulationIndex_, 0, static_cast<int>(simulations_.size()) - 1);
  refreshAll();
}

void MainWindow::runSimulation()
{
  SimulationConfig* simulation = activeSimulation();
  if (simulation == nullptr)
  {
    return;
  }
  updatePreview();
  QString error;
  pendingRunSimulationId_ = simulation->id;
  if (!runner_.startRun(buildSimulationJson(), &error))
  {
    pendingRunSimulationId_.clear();
    QMessageBox::warning(this, "Ruptura Lab", error);
    refreshRunPanel();
    return;
  }
  simulation->lastRun = runner_.run(runner_.latestRunId());
  simulation->lastRunLog = runner_.logTail(simulation->lastRun);
  pendingRunSimulationId_.clear();
  refreshSimulationList();
  refreshRunPanel();
}

void MainWindow::cancelSimulation()
{
  QString error;
  if (!runner_.cancelRun(&error))
  {
    QMessageBox::information(this, "Ruptura Lab", error);
  }
}

void MainWindow::openAnalysis()
{
  const SimulationConfig* simulation = activeSimulation();
  if (simulation == nullptr || simulation->lastRun.id.isEmpty())
  {
    QMessageBox::information(this, "Ruptura Lab", "No run is available for the selected simulation.");
    return;
  }

  QString error;
  if (!runner_.openAnalysis(simulation->lastRun.id, &error))
  {
    QMessageBox::information(this, "Ruptura Lab", error);
  }
}

QLineEdit* MainWindow::addTextField(QFormLayout* form, const QString& label, const QString& value,
                                    const std::function<void(const QString&)>& setter)
{
  auto* edit = new QLineEdit(value);
  connect(edit, &QLineEdit::textEdited, this,
          [this, setter](const QString& text)
          {
            setter(text);
            updatePreview();
          });
  form->addRow(label, edit);
  return edit;
}

QLineEdit* MainWindow::addNumberField(QFormLayout* form, const QString& label, double value,
                                      const std::function<void(double)>& setter)
{
  auto* edit = new QLineEdit(numberText(value));
  connect(edit, &QLineEdit::textEdited, this,
          [this, setter](const QString& text)
          {
            bool ok = false;
            const double value = text.toDouble(&ok);
            if (ok)
            {
              setter(value);
              updatePreview();
              updateFractionTotal();
            }
          });
  form->addRow(label, edit);
  return edit;
}

QLineEdit* MainWindow::addOptionalNumberField(QFormLayout* form, const QString& label,
                                              const std::optional<double>& value,
                                              const std::function<void(std::optional<double>)>& setter)
{
  auto* edit = new QLineEdit(optionalText(value));
  connect(edit, &QLineEdit::textEdited, this,
          [this, setter](const QString& text)
          {
            const QString trimmed = text.trimmed();
            if (trimmed.isEmpty())
            {
              setter(std::nullopt);
              updatePreview();
              return;
            }
            bool ok = false;
            const double value = trimmed.toDouble(&ok);
            if (ok)
            {
              setter(value);
              updatePreview();
            }
          });
  form->addRow(label, edit);
  return edit;
}

QCheckBox* MainWindow::addCheckField(QFormLayout* form, const QString& label, bool checked,
                                     const std::function<void(bool)>& setter)
{
  auto* check = new QCheckBox;
  check->setChecked(checked);
  connect(check, &QCheckBox::toggled, this,
          [this, setter](bool value)
          {
            setter(value);
            updatePreview();
            updateFractionTotal();
          });
  form->addRow(label, check);
  return check;
}

QComboBox* MainWindow::addComboField(QFormLayout* form, const QString& label, const QString& value,
                                     const QList<QPair<QString, QString>>& options,
                                     const std::function<void(const QString&)>& setter)
{
  auto* combo = new QComboBox;
  sizeComboForOptions(combo, options);
  for (const auto& option : options)
  {
    combo->addItem(option.first, option.second);
  }
  const int index = combo->findData(value);
  if (index >= 0)
  {
    combo->setCurrentIndex(index);
  }
  connect(combo, qOverload<int>(&QComboBox::currentIndexChanged), this,
          [this, combo, setter](int)
          {
            setter(combo->currentData().toString());
            updatePreview();
            updateFractionTotal();
          });
  form->addRow(label, combo);
  return combo;
}

void MainWindow::addParameterFields(QFormLayout* form, QVariantMap& parameters,
                                    const QList<ParameterDefinition>& definitions)
{
  parameters = defaultsFor(definitions, parameters);
  for (const ParameterDefinition& definition : definitions)
  {
    if (definition.boolean)
    {
      auto* check = new QCheckBox;
      check->setChecked(parameterBool(parameters, definition));
      connect(check, &QCheckBox::toggled, this,
              [this, &parameters, definition](bool value)
              {
                parameters.insert(definition.key, value);
                updatePreview();
              });
      form->addRow(definition.label, check);
      continue;
    }
    if (definition.temperatureLaws != 0)
    {
      addTemperatureLawField(form, definition.label, parameters.value(definition.key), definition.temperatureLaws,
                             definition.heatLinked, [&parameters, definition](const QVariant& value)
                             { parameters.insert(definition.key, value); });
      continue;
    }

    auto* edit = new QLineEdit(numberText(parameterNumber(parameters, definition)));
    connect(edit, &QLineEdit::textEdited, this,
            [this, &parameters, definition](const QString& text)
            {
              bool ok = false;
              double value = text.toDouble(&ok);
              if (ok)
              {
                if (definition.integer)
                {
                  value = static_cast<double>(static_cast<int>(value));
                }
                parameters.insert(definition.key, value);
                updatePreview();
              }
            });
    form->addRow(definition.label, edit);
  }
}

void MainWindow::addTemperatureLawField(QFormLayout* form, const QString& label, const QVariant& value,
                                        unsigned allowedLaws, bool heatLinked,
                                        const std::function<void(const QVariant&)>& setter)
{
  // The value is a number (constant) or a law {"Value", "Law", <coefficient>, "ReferenceTemperature"}, written to
  // simulation.json as is. An empty coefficient of a heat-linked VantHoff law takes the heat of adsorption (or of
  // reaction), and an empty T_ref the ReferenceTemperature of the component (or Kinetics).
  const QVariantMap law = value.typeId() == QMetaType::QVariantMap ? value.toMap() : QVariantMap{};
  const double number = law.isEmpty() ? value.toDouble() : law.value("Value").toDouble();
  const QString lawName = law.value("Law", "Constant").toString();

  auto* container = new QWidget;
  auto* layout = new QHBoxLayout(container);
  layout->setContentsMargins(0, 0, 0, 0);
  auto* valueEdit = new QLineEdit(numberText(number));
  auto* lawCombo = new QComboBox;
  lawCombo->addItem("Constant", "Constant");
  for (const TemperatureLaw::Type type : {TemperatureLaw::Type::VantHoff, TemperatureLaw::Type::Arrhenius,
                                          TemperatureLaw::Type::Exponential, TemperatureLaw::Type::Linear,
                                          TemperatureLaw::Type::Inverse})
  {
    if ((allowedLaws & TemperatureLaw::mask(type)) == 0) continue;
    const QString name = QString::fromUtf8(TemperatureLaw::name(type).data(),
                                           static_cast<qsizetype>(TemperatureLaw::name(type).size()));
    lawCombo->addItem(name, name);
  }
  lawCombo->setCurrentIndex(std::max(0, lawCombo->findData(lawName, Qt::UserRole, Qt::MatchFixedString)));
  auto* coefficientEdit = new QLineEdit;
  auto* referenceEdit = new QLineEdit(law.contains("ReferenceTemperature")
                                          ? numberText(law.value("ReferenceTemperature").toDouble())
                                          : QString{});
  referenceEdit->setPlaceholderText("T_ref [K]");
  layout->addWidget(valueEdit, 2);
  layout->addWidget(lawCombo, 2);
  layout->addWidget(coefficientEdit, 2);
  layout->addWidget(referenceEdit, 1);

  const auto coefficientKey = [lawCombo]()
  {
    const std::optional<TemperatureLaw::Type> type =
        TemperatureLaw::fromName(lawCombo->currentData().toString().toStdString());
    const std::string_view key = TemperatureLaw::coefficientKey(type.value_or(TemperatureLaw::Type::Constant));
    return QString::fromUtf8(key.data(), static_cast<qsizetype>(key.size()));
  };
  const auto showLaw = [=]()
  {
    const QString key = coefficientKey();
    coefficientEdit->setVisible(!key.isEmpty());
    referenceEdit->setVisible(!key.isEmpty());
    const QString units = key == "Energy" || key == "ActivationEnergy" ? " [J/mol]" : key == "Slope" ? " [1/K]" : "";
    coefficientEdit->setPlaceholderText(key == "Energy" && heatLinked ? "Energy [J/mol] (empty: heat)" : key + units);
  };
  coefficientEdit->setText(law.contains(coefficientKey()) ? numberText(law.value(coefficientKey()).toDouble())
                                                          : QString{});
  showLaw();

  const auto update = [this, valueEdit, lawCombo, coefficientEdit, referenceEdit, coefficientKey, setter]()
  {
    bool ok = false;
    const double valueNumber = valueEdit->text().toDouble(&ok);
    if (!ok) return;
    const QString key = coefficientKey();
    if (key.isEmpty())
    {
      setter(valueNumber);
    }
    else
    {
      QVariantMap result{{"Value", valueNumber}, {"Law", lawCombo->currentData().toString()}};
      const double coefficient = coefficientEdit->text().toDouble(&ok);
      if (ok) result.insert(key, coefficient);
      const double referenceTemperature = referenceEdit->text().toDouble(&ok);
      if (ok) result.insert("ReferenceTemperature", referenceTemperature);
      setter(result);
    }
    updatePreview();
  };
  connect(valueEdit, &QLineEdit::textEdited, this, update);
  connect(coefficientEdit, &QLineEdit::textEdited, this, update);
  connect(referenceEdit, &QLineEdit::textEdited, this, update);
  connect(lawCombo, qOverload<int>(&QComboBox::currentIndexChanged), this,
          [showLaw, update](int)
          {
            showLaw();
            update();
          });
  form->addRow(label, container);
}

void MainWindow::addSettingNumbers(QFormLayout* form, QJsonObject& settings,
                                   const QList<ParameterDefinition>& definitions)
{
  for (const auto& definition : definitions)
    addNumberField(form, definition.label, readNumber(settings, definition.key, definition.defaultValue),
                   [&settings, definition](double value)
                   {
                     settings.insert(definition.key, definition.integer ? QJsonValue(static_cast<int>(value))
                                                                        : QJsonValue(value));
                   });
}

void MainWindow::addSettingCombo(QFormLayout* form, const QString& label, QJsonObject& settings, const QString& key,
                                 const QString& fallback, const QStringList& values)
{
  addComboField(form, label, choice(settings.value(key).toString(fallback), values), valueOptions(values),
                [this, &settings, key](const QString& value)
                {
                  settings.insert(key, value);
                  QTimer::singleShot(0, this, [this]() { refreshAll(); });
                });
}

void MainWindow::addJsonField(QFormLayout* form, const QString& label, const QJsonValue& value,
                              const std::function<void(const QJsonValue&)>& setter)
{
  auto* container = new QWidget;
  auto* layout = new QVBoxLayout(container);
  layout->setContentsMargins(0, 0, 0, 0);
  auto* editor = new QPlainTextEdit;
  editor->setPlainText(QString::fromUtf8(value.isArray() ? QJsonDocument(value.toArray()).toJson()
                                                        : QJsonDocument(value.toObject()).toJson()));
  editor->setMaximumHeight(140);
  auto* apply = new QPushButton("Apply");
  auto* errorLabel = new QLabel;
  errorLabel->setWordWrap(true);
  layout->addWidget(editor);
  layout->addWidget(errorLabel);
  layout->addWidget(apply);
  connect(apply, &QPushButton::clicked, this, [this, editor, errorLabel, setter, array = value.isArray()]()
          {
            QJsonParseError error;
            const auto document = QJsonDocument::fromJson(editor->toPlainText().toUtf8(), &error);
            if (error.error != QJsonParseError::NoError || (array ? !document.isArray() : !document.isObject()))
            {
              errorLabel->setText(error.error == QJsonParseError::NoError
                                      ? (array ? "Enter a JSON array." : "Enter a JSON object.") : error.errorString());
              return;
            }
            setter(array ? QJsonValue(document.array()) : QJsonValue(document.object()));
            errorLabel->setText("Applied");
            updatePreview();
          });
  form->addRow(label, container);
}

#include "simulationinput.h"

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QHash>
#include <QJsonArray>
#include <QJsonDocument>
#include <QStandardPaths>
#include <QTemporaryDir>
#include <QStringList>
#include <stdexcept>

#include "inputreader.h"

namespace
{
QString canonicalKey(const QString& key)
{
  static const QStringList keys = QString::fromLatin1(
      "SimulationType MixturePredictionMethod MPDSettings IASTMethod BreakthroughIntegrator BoundaryCondition "
      "FluidPhase LiquidDensity pHMode pHValue pKw pHComponent PressureScale ReadColumnFile DisplayName "
      "DebugForceMultibed Temperature ColumnVoidFraction DynamicViscosity ParticleDiameter ParticleDensity "
      "InletPressure OutletPressure PressureStart PressureEnd NumberOfPressurePoints PressureGradient "
      "ColumnEntranceVelocity NumberOfInitTimeSteps TimeStep CVODERelativeTolerance "
      "CVODEAbsoluteToleranceConcentration CVODEAbsoluteToleranceLoading CVODEAbsoluteToleranceTemperature "
      "CVODEMaximumTimeStep CVODEMaximumOrder CVODELinearSolver CVODEKrylovDimension PrintEvery WriteEvery ColumnLength "
      "ColumnDistances NumberOfGridPoints ColumnPressure ColumnLoading ColumnError NumberOfTimeSteps Geometry "
      "InfluxTemperature wallDensity gasThermalConductivity wallThermalConductivity heatTransferGasWall "
      "heatTransferGasSolid heatTransferWallExternal heatCapacityGas liquidThermalConductivity "
      "heatTransferLiquidSolid heatTransferLiquidWall heatCapacityLiquid heatCapacitySolid heatCapacityWall "
      "energyBalance SwingAdsorptionPhases ColumnSections Adsorbents Components Reactions Name FileName "
      "CarrierGas GasPhaseMolFraction LiquidPhaseConcentration InitialLiquidPhaseConcentration "
      "MassTransferCoefficient AxialDispersionCoefficient MolecularWeight HeatOfAdsorption referenceTemperature "
      "nonIsothermal ChemisorptionSites PhysisorptionSites Type Parameters RateEquation RateEquationParameters "
      "MaximumLoading ChannelShape InternalChannelDimension InternalDiameter OuterDiameter NumberOfChannels "
      "WashcoatThickness WashcoatVolumePerChannelVolume ForchheimerCoefficient AdsorbentLength MixFraction ComponentParameters "
      "Adsorbent Length InterfaceLength Phase Style Reactants Products Stoichiometry Kinetics Site "
      "ForwardOrders BackwardOrders RateLimitTime forwardRateCoefficient forwardActivationEnergy "
      "equilibriumConstant gibbsFreeEnergy").split(' ');
  for (const QString& candidate : keys)
  {
    if (candidate.compare(key, Qt::CaseInsensitive) == 0) return candidate;
  }
  return key;
}

QJsonValue normalizeEnum(const QString& key, const QJsonValue& value)
{
  if (!value.isString()) return value;
  static const QHash<QString, QStringList> options{
      {"SimulationType", {"Breakthrough", "SwingAdsorption", "MixturePrediction", "Fitting", "Test"}},
      {"MixturePredictionMethod", {"IAST", "SIAST", "EI", "SEI", "SCI", "SPI", "MPD"}},
      {"IASTMethod", {"FastIAST", "NestedLoopBisection"}},
      {"BreakthroughIntegrator", {"RungeKutta3", "CVODE", "SIRK3"}},
      {"CVODELinearSolver", {"Dense", "SPGMR"}},
      {"FluidPhase", {"Gas", "Liquid"}}, {"pHMode", {"Fixed", "HPlus", "OHMinus"}},
      {"PressureScale", {"Log", "Linear"}},
      {"BoundaryCondition", {"InletPressureInletVelocity", "InletPressureOutletPressure", "InletVelocityOutletPressure",
                              "FixedVelocity", "FixedPressureInletVelocity"}},
      {"Type", {"PackedBed", "Monolith", "Langmuir", "pH-Langmuir", "Anti-Langmuir", "BET", "Henry", "Freundlich",
                "Sips", "Langmuir-Freundlich", "Redlich-Peterson", "Toth", "Unilan", "OBrien&Myers", "Quadratic",
                "Temkin", "Bingel&Walton", "GAB"}},
      {"ChannelShape", {"triangular", "square", "hexagonal", "circular"}},
      {"RateEquation", {"FirstOrder", "PseudoNth", "Avrami", "General", "Elovich", "None"}},
      {"Phase", {"Physisorbed", "Chemisorbed", "PoreConcentration"}},
      {"Style", {"GeneralPowerLaw", "LangmuirHinshelwood", "LangmuirHinshelwoodHougenWatson"}}};
  for (const auto& candidate : options.value(key))
    if (candidate.compare(value.toString(), Qt::CaseInsensitive) == 0) return candidate;
  return value;
}

QJsonValue normalize(const QJsonValue& value, bool names = false)
{
  if (value.isArray())
  {
    QJsonArray result;
    for (const auto& item : value.toArray()) result.append(normalize(item));
    return result;
  }
  if (!value.isObject()) return value;
  QJsonObject result;
  const auto object = value.toObject();
  for (auto it = object.begin(); it != object.end(); ++it)
  {
    const QString key = names ? it.key() : canonicalKey(it.key());
    if (result.contains(key)) throw std::runtime_error("Duplicate input key: " + key.toStdString());
    // RateEquationParameters has an intentionally case-sensitive CLI schema.
    result.insert(key, key == "RateEquationParameters" ? it.value()
                       : normalizeEnum(key, normalize(it.value(), key == "ComponentParameters" || key == "Components")));
  }
  return result;
}

QJsonValue resolve(const QJsonValue& value, const QDir& directory)
{
  if (value.isArray())
  {
    QJsonArray result;
    for (const auto& item : value.toArray()) result.append(resolve(item, directory));
    return result;
  }
  if (!value.isObject()) return value;
  QJsonObject object = value.toObject();
  for (auto it = object.begin(); it != object.end(); ++it)
  {
    if ((it.key() == "ReadColumnFile" || it.key() == "FileName") && it.value().isString() &&
        !it.value().toString().isEmpty())
      it.value() = directory.absoluteFilePath(it.value().toString());
    else
      it.value() = resolve(it.value(), directory);
  }
  return object;
}

QJsonValue patch(const QJsonValue& original, const QJsonValue& before, const QJsonValue& after)
{
  if (before == after) return original;
  if (before.isObject() && after.isObject())
  {
    QJsonObject result = original.toObject();
    const auto oldObject = before.toObject();
    const auto newObject = after.toObject();
    for (auto it = oldObject.begin(); it != oldObject.end(); ++it)
      if (!newObject.contains(it.key())) result.remove(it.key());
    for (auto it = newObject.begin(); it != newObject.end(); ++it)
    {
      if (oldObject.value(it.key()) != it.value())
        result.insert(it.key(), patch(result.value(it.key()), oldObject.value(it.key()), it.value()));
    }
    return result;
  }
  if (before.isArray() && after.isArray() && original.isArray())
  {
    auto result = original.toArray();
    const auto oldArray = before.toArray();
    const auto newArray = after.toArray();
    if (result.size() == oldArray.size() && oldArray.size() == newArray.size())
    {
      for (qsizetype i = 0; i < result.size(); ++i) result[i] = patch(result[i], oldArray[i], newArray[i]);
      return result;
    }
  }
  return after;
}
}

QJsonObject LabInput::canonicalize(const QJsonObject& input)
{
  auto result = normalize(input).toObject();
  if (result.value("Components").isObject())
  {
    QJsonArray components;
    const auto object = result.value("Components").toObject();
    for (const auto& component : object) components.append(component);
    result.insert("Components", components);
  }
  return result;
}

QJsonObject LabInput::resolveFiles(const QJsonObject& input, const QString& directory)
{
  return resolve(input, QDir(directory.isEmpty() ? QDir::currentPath() : directory)).toObject();
}

QJsonObject LabInput::applyChanges(const QJsonObject& original, const QJsonObject& before, const QJsonObject& after)
{
  return patch(original, before, after).toObject();
}

QJsonObject LabInput::renameComponentReferences(const QJsonObject& input, const QHash<QString, QString>& names)
{
  auto result = input;
  auto adsorbents = result.value("Adsorbents").toArray();
  for (qsizetype i = 0; i < adsorbents.size(); ++i)
  {
    auto bed = adsorbents[i].toObject();
    auto parameters = bed.value("ComponentParameters").toObject();
    QJsonObject renamed;
    for (auto it = parameters.begin(); it != parameters.end(); ++it)
      renamed.insert(names.value(it.key(), it.key()), it.value());
    if (bed.contains("ComponentParameters")) bed.insert("ComponentParameters", renamed);
    auto components = bed.value("Components").toArray();
    for (qsizetype j = 0; j < components.size(); ++j)
    {
      auto component = components[j].toObject();
      const auto name = component.value("Name").toString();
      if (names.contains(name)) component.insert("Name", names.value(name));
      components[j] = component;
    }
    if (bed.contains("Components")) bed.insert("Components", components);
    adsorbents[i] = bed;
  }
  if (result.contains("Adsorbents")) result.insert("Adsorbents", adsorbents);
  auto reactions = result.value("Reactions").toArray();
  for (qsizetype i = 0; i < reactions.size(); ++i)
  {
    auto reaction = reactions[i].toObject();
    for (const auto& key : {QString("Reactants"), QString("Products")})
    {
      auto participants = reaction.value(key).toArray();
      for (qsizetype j = 0; j < participants.size(); ++j)
      {
        if (!participants[j].isString()) continue;
        const auto name = participants[j].toString();
        if (names.contains(name)) participants[j] = names.value(name);
        else
          for (auto it = names.begin(); it != names.end(); ++it)
            if (it.key().compare(name, Qt::CaseInsensitive) == 0) { participants[j] = it.value(); break; }
      }
      reaction.insert(key, participants);
    }
    reactions[i] = reaction;
  }
  if (result.contains("Reactions")) result.insert("Reactions", reactions);
  return result;
}

bool LabInput::supported(const QJsonObject& input, QString* error)
{
  const QString type = input.value("SimulationType").toString("Breakthrough");
  const bool knownType = type.compare("Breakthrough", Qt::CaseInsensitive) == 0 ||
                         type.compare("SwingAdsorption", Qt::CaseInsensitive) == 0 ||
                         type.compare("MixturePrediction", Qt::CaseInsensitive) == 0 ||
                         type.compare("Test", Qt::CaseInsensitive) == 0;
  if (!knownType || input.value("MixturePredictionMethod").toString().compare("MPD", Qt::CaseInsensitive) == 0 ||
      input.contains("MPDSettings"))
  {
    if (error) *error = "Ruptura Lab supports breakthrough, swing adsorption, and mixture prediction; MPD and fitting are excluded.";
    return false;
  }
#if !BUILD_SUNDIALS
  if (input.value("BreakthroughIntegrator").toString().compare("CVODE", Qt::CaseInsensitive) == 0)
  {
    if (error) *error = "CVODE is unavailable in this build. Select RungeKutta3 or SIRK3, or rebuild with BUILD_SUNDIALS=ON.";
    return false;
  }
#endif
  return true;
}

bool LabInput::validate(const QJsonObject& input, QString* error)
{
  return read(input, error) != nullptr;
}

std::unique_ptr<InputReader> LabInput::read(const QJsonObject& input, QString* error)
{
  if (!supported(input, error)) return {};
  try
  {
    QTemporaryDir directory;
    QFile file(directory.filePath("simulation.json"));
    if (!directory.isValid() || !file.open(QIODevice::WriteOnly))
      throw std::runtime_error("Could not create a temporary simulation input for validation.");
    const QByteArray data = QJsonDocument(input).toJson();
    if (file.write(data) != data.size()) throw std::runtime_error("Could not write the simulation input for validation.");
    file.close();
    return std::make_unique<InputReader>(file.fileName().toStdString());
  }
  catch (const std::exception& exception)
  {
    if (error) *error = QString::fromUtf8(exception.what());
    return {};
  }
}

QString LabInput::findCli(const QString& applicationDirectory, const QString& buildExecutable)
{
  const QString override = qEnvironmentVariable("RUPTURA_CLI_EXECUTABLE");
  if (!override.isEmpty()) return QStandardPaths::findExecutable(override);
  const QString sibling = QDir(applicationDirectory).filePath("ruptura");
  for (const QString& path : {sibling, buildExecutable})
    if (QFileInfo(path).isFile() && QFileInfo(path).isExecutable()) return path;
  return QStandardPaths::findExecutable("ruptura");
}

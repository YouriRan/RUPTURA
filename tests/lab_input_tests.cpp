#include <gtest/gtest.h>

#include <algorithm>
#include <QApplication>
#include <QDir>
#include <QEventLoop>
#include <QFile>
#include <QJsonDocument>
#include <QProcess>
#include <QTemporaryDir>
#include <QTimer>

#include "mainwindow.h"
#include "simulationinput.h"

struct MainWindowTestAccess
{
  static bool import(MainWindow& window, const QJsonObject& input, QString* error, const QString& directory = {})
  { return window.importSimulationJson(input, error, directory); }
  static QJsonObject output(const MainWindow& window) { return window.buildSimulationJson(); }
  static QJsonObject modeled(const MainWindow& window) { return window.buildModeledSimulationJson(); }
  static QJsonObject state(const MainWindow& window) { return window.buildStateJson(); }
  static bool restore(MainWindow& window, const QJsonObject& state, QString* error)
  { return window.loadStateJson(state, error); }
  static SimulationConfig& simulation(MainWindow& window) { return *window.activeSimulation(); }
  static QList<ComponentConfig>& components(MainWindow& window) { return window.components_; }
  static ColumnConfig& column(MainWindow& window) { return window.columns_.front(); }
};

namespace
{
QJsonObject example(const QString& path)
{
  QFile file(QDir(QString::fromUtf8(RUPTURA_SOURCE_DIR)).filePath(path));
  if (!file.open(QIODevice::ReadOnly)) throw std::runtime_error("Missing example: " + path.toStdString());
  auto input = QJsonDocument::fromJson(file.readAll()).object();
#if !BUILD_SUNDIALS
  if (input.value("BreakthroughIntegrator").toString() == "CVODE") input.insert("BreakthroughIntegrator", "RungeKutta3");
#endif
  return input;
}

QJsonObject mixture()
{
  return QJsonDocument::fromJson(R"({
    "SimulationType": "MixturePrediction", "MixturePredictionMethod": "SPI",
    "Temperature": 300, "PressureStart": 1, "PressureEnd": 10, "NumberOfPressurePoints": 2,
    "Components": [{"Name": "CO2", "GasPhaseMolFraction": 1, "MassTransferCoefficient": 0.04,
                    "PhysisorptionSites": [{"Type": "Langmuir", "Parameters": [1, 0.00001]}]}]
  })").object();
}

void writeJson(const QString& path, const QJsonObject& input)
{
  QFile file(path);
  ASSERT_TRUE(file.open(QIODevice::WriteOnly));
  const auto bytes = QJsonDocument(input).toJson();
  ASSERT_EQ(file.write(bytes), bytes.size());
}

class LabExamples : public testing::TestWithParam<const char*> {};

TEST_P(LabExamples, PreservesImportedInputAndSavedState)
{
  const auto input = example(QString::fromUtf8(GetParam()));
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  EXPECT_EQ(MainWindowTestAccess::output(window), LabInput::canonicalize(input));
  // Exercise the exporter as well as the lossless import path.
  EXPECT_TRUE(LabInput::validate(MainWindowTestAccess::modeled(window), &error)) << error.toStdString();
  MainWindow restored;
  ASSERT_TRUE(MainWindowTestAccess::restore(restored, MainWindowTestAccess::state(window), &error)) << error.toStdString();
  EXPECT_EQ(MainWindowTestAccess::output(restored), MainWindowTestAccess::output(window));
  MainWindowTestAccess::simulation(restored).displayName = "Edited in Lab";
  auto expected = LabInput::canonicalize(input);
  expected.insert("DisplayName", "Edited in Lab");
  EXPECT_EQ(MainWindowTestAccess::output(restored), expected);
}

INSTANTIATE_TEST_SUITE_P(SupportedExamples, LabExamples, testing::Values(
    "examples/BEA-alkanes-C7/psa/simulation.json",
    "examples/BEA-alkanes-C7/tsa/simulation.json",
    "examples/Chemisorption-CO2-N2/monolith/simulation.json",
    "examples/MOR-CO2-C3H8/siast/simulation.json",
    "examples/Liquid/pH-dependent/simulation.json",
    "examples/Liquid/multibed-pH-front/simulation.json",
    "examples/Chemisorption-CO2/breakthrough/simulation.json",
    "examples/Chemisorption-CO2-N2/breakthrough/simulation.json",
    "examples/Multibed/CO2-N2-polishing/simulation.json",
    "examples/Multibed/CO2-N2-AB-mixture/simulation.json",
    "examples/Multibed/CO2-C3H8-gradient/simulation.json",
    "examples/Multibed/C6-alkanes-layered/simulation.json",
    "examples/Reactions/physisorbed-power-law/simulation.json",
    "examples/Reactions/pore-langmuir-hinshelwood/simulation.json"));

TEST(LabInput, ExcludesMpdAndFittingWithoutReplacingState)
{
  MainWindow window;
  const auto before = MainWindowTestAccess::state(window);
  QString error;
  for (const QString& type : {QString("Fitting"), QString("fitting")})
  {
    auto input = mixture();
    input.insert("SimulationType", type);
    EXPECT_FALSE(MainWindowTestAccess::import(window, input, &error));
    EXPECT_TRUE(error.contains("excluded"));
    EXPECT_EQ(MainWindowTestAccess::state(window), before);
  }
  auto input = mixture();
  input.insert("MixturePredictionMethod", "mpd");
  EXPECT_FALSE(MainWindowTestAccess::import(window, input, &error));
  EXPECT_EQ(MainWindowTestAccess::state(window), before);
  auto state = before;
  auto simulations = state.value("Simulations").toArray();
  auto simulation = simulations[0].toObject();
  simulation.insert("Type", "Fitting");
  simulations[0] = simulation;
  state.insert("Simulations", simulations);
  EXPECT_FALSE(MainWindowTestAccess::restore(window, state, &error));
}

TEST(LabInput, LegacySavedChemisorptionExportsCurrentCliSchema)
{
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, mixture(), &error));
  auto state = MainWindowTestAccess::state(window);
  auto components = state.value("Components").toArray();
  auto component = components[0].toObject();
  component.insert("ChemisorptionSites", QJsonArray{QJsonObject{
      {"Type", "FirstOrder"}, {"Parameters", QJsonObject{{"rateCoefficient", 0.04},
          {"maximumLoading", 0.8}, {"heatOfChemisorption", 42000.0}}},
      {"Isotherm", QJsonObject{{"Type", "Langmuir"}, {"Parameters", QJsonArray{1.0, 1.0e-5}}}}}});
  components[0] = component;
  state.insert("Components", components);
  ASSERT_TRUE(MainWindowTestAccess::restore(window, state, &error));
  const auto input = MainWindowTestAccess::modeled(window);
  ASSERT_TRUE(LabInput::validate(input, &error)) << error.toStdString();
  const auto site = input.value("Components").toArray()[0].toObject().value("ChemisorptionSites").toArray()[0].toObject();
  EXPECT_EQ(site.value("Type").toString(), "Langmuir");
  EXPECT_EQ(site.value("RateEquation").toString(), "FirstOrder");
  EXPECT_DOUBLE_EQ(site.value("MaximumLoading").toDouble(), 0.8);
}

TEST(LabInput, SolverSettingsAndExplicitGridSurviveEdits)
{
  auto input = example("examples/Multibed/CO2-C3H8-gradient/simulation.json");
  input.insert("CVODERelativeTolerance", 2.0e-7);
  input.insert("CVODEAbsoluteToleranceConcentration", 3.0e-9);
  input.insert("CVODEMaximumTimeStep", 0.01);
  input.insert("CVODELinearSolver", "SPGMR");
  input.insert("CVODEKrylovDimension", 12);
  input.insert("IASTMethod", "NestedLoopBisection");
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  MainWindowTestAccess::simulation(window).settings.insert("CVODERelativeTolerance", 5.0e-8);
  input.insert("CVODERelativeTolerance", 5.0e-8);
  EXPECT_EQ(MainWindowTestAccess::output(window), input);
  EXPECT_TRUE(LabInput::validate(MainWindowTestAccess::output(window), &error)) << error.toStdString();
}

TEST(LabInput, TransportedPhUsesComponentIdentityWhenFeedsReorder)
{
  const auto input = example("examples/Liquid/multibed-pH-front/simulation.json");
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  auto& simulation = MainWindowTestAccess::simulation(window);
  ASSERT_FALSE(simulation.pHComponentId.isEmpty());
  std::reverse(simulation.componentFeeds.begin(), simulation.componentFeeds.end());
  const auto output = MainWindowTestAccess::output(window);
  const auto originalComponents = input.value("Components").toArray();
  const QString expectedName = originalComponents[input.value("pHComponent").toInt()].toObject().value("Name").toString();
  EXPECT_EQ(output.value("Components").toArray()[output.value("pHComponent").toInt()].toObject().value("Name").toString(), expectedName);
  EXPECT_TRUE(LabInput::validate(output, &error)) << error.toStdString();
}

TEST(LabInput, ReactionOrdersSurviveStateAndComponentEdits)
{
  auto input = example("examples/Reactions/pore-langmuir-hinshelwood/simulation.json");
  auto reactions = input.value("Reactions").toArray();
  auto reaction = reactions[0].toObject();
  QJsonArray forward, backward;
  for (qsizetype i = 0; i < reaction.value("Reactants").toArray().size(); ++i) forward.append(0.5);
  for (qsizetype i = 0; i < reaction.value("Products").toArray().size(); ++i) backward.append(1.5);
  reaction.insert("ForwardOrders", forward);
  reaction.insert("BackwardOrders", backward);
  reactions[0] = reaction;
  input.insert("Reactions", reactions);
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  MainWindow restored;
  ASSERT_TRUE(MainWindowTestAccess::restore(restored, MainWindowTestAccess::state(window), &error));
  const auto output = MainWindowTestAccess::modeled(restored);
  EXPECT_EQ(output.value("Reactions").toArray()[0].toObject().value("ForwardOrders").toArray(), forward);
  EXPECT_EQ(output.value("Reactions").toArray()[0].toObject().value("BackwardOrders").toArray(), backward);
  EXPECT_TRUE(LabInput::validate(output, &error)) << error.toStdString();
}

TEST(LabInput, CaseInsensitiveFieldsAndObjectComponents)
{
  auto input = mixture();
  auto component = input.take("Components").toArray()[0].toObject();
  component.insert("name", component.take("Name"));
  input.insert("components", QJsonObject{{"0", component}});
  input.insert("simulationtype", input.take("SimulationType").toString().toLower());
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  EXPECT_EQ(MainWindowTestAccess::output(window).value("Components").toArray().size(), 1);
  EXPECT_TRUE(LabInput::validate(MainWindowTestAccess::output(window), &error)) << error.toStdString();
}

TEST(LabInput, RelativeRestartPathsResolveAtImport)
{
  QTemporaryDir source;
  auto input = example("examples/Liquid/pH-dependent/simulation.json");
  input.insert("ReadColumnFile", "initial.json");
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error, source.path())) << error.toStdString();
  EXPECT_EQ(MainWindowTestAccess::output(window).value("ReadColumnFile").toString(), source.filePath("initial.json"));
}

TEST(LabInput, AllChemisorptionRateEquationsExportValidIndependentKinetics)
{
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, mixture(), &error));
  auto& component = MainWindowTestAccess::components(window).front();
  for (const QString& type : {QString("FirstOrder"), QString("PseudoNth"), QString("Avrami"),
                              QString("General"), QString("Elovich")})
  {
    SCOPED_TRACE(type.toStdString());
    ChemisorptionSiteConfig site;
    site.type = type;
    site.parameters.insert("maximumLoading", 0.6);
    component.chemisorptionSites = {site};
    auto output = MainWindowTestAccess::modeled(window);
    EXPECT_TRUE(LabInput::validate(output, &error)) << error.toStdString();
    MainWindow restored;
    ASSERT_TRUE(MainWindowTestAccess::import(restored, output, &error)) << error.toStdString();
    EXPECT_EQ(MainWindowTestAccess::modeled(restored).value("Components"), output.value("Components"));
  }
}

TEST(LabInput, AddedIsothermsRetainEveryParameter)
{
  for (const auto& pair : {qMakePair(QString("GAB"), QJsonArray{1.0, 1.0, 0.5, 0.0, 0.0}),
                           qMakePair(QString("Bingel&Walton"), QJsonArray{1.0, 1.0e-5, 2.0e-5})})
  {
    auto input = mixture();
    auto components = input.value("Components").toArray();
    auto component = components[0].toObject();
    component.insert("PhysisorptionSites", QJsonArray{QJsonObject{{"Type", pair.first}, {"Parameters", pair.second}}});
    components[0] = component;
    input.insert("Components", components);
    MainWindow window;
    QString error;
    ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
    EXPECT_EQ(MainWindowTestAccess::modeled(window).value("Components"), input.value("Components"));
    EXPECT_TRUE(LabInput::validate(MainWindowTestAccess::modeled(window), &error)) << error.toStdString();
  }
}

TEST(LabInput, MultibedLengthEditsKeepDistinctAdsorbentIsotherms)
{
  const auto input = example("examples/Multibed/CO2-N2-polishing/simulation.json");
  MainWindow window;
  QString error;
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  MainWindowTestAccess::column(window).beds[0].length += 0.01;
  const auto output = MainWindowTestAccess::output(window);
  const auto beforeBeds = input.value("Adsorbents").toArray();
  const auto afterBeds = output.value("Adsorbents").toArray();
  ASSERT_EQ(afterBeds.size(), beforeBeds.size());
  for (qsizetype i = 0; i < beforeBeds.size(); ++i)
    EXPECT_EQ(afterBeds[i].toObject().value("ComponentParameters"), beforeBeds[i].toObject().value("ComponentParameters"));
  EXPECT_TRUE(LabInput::validate(output, &error)) << error.toStdString();
}

TEST(LabInput, RenamingComponentsUpdatesNamedReactionAndAdsorbentReferences)
{
  MainWindow window;
  QString error;
  auto input = example("examples/Reactions/pore-langmuir-hinshelwood/simulation.json");
  auto reactions = input.value("Reactions").toArray();
  const auto components = input.value("Components").toArray();
  for (qsizetype i = 0; i < reactions.size(); ++i)
  {
    auto reaction = reactions[i].toObject();
    for (const auto& key : {QString("Reactants"), QString("Products")})
    {
      auto participants = reaction.value(key).toArray();
      for (qsizetype j = 0; j < participants.size(); ++j)
        if (participants[j].isDouble()) participants[j] = components[participants[j].toInt()].toObject().value("Name").toString().toLower();
      reaction.insert(key, participants);
    }
    reactions[i] = reaction;
  }
  input.insert("Reactions", reactions);
  ASSERT_TRUE(MainWindowTestAccess::import(window, input, &error)) << error.toStdString();
  MainWindowTestAccess::components(window)[1].name = "Renamed reactant";
  EXPECT_TRUE(LabInput::validate(MainWindowTestAccess::output(window), &error)) << error.toStdString();
  EXPECT_TRUE(LabInput::validate(MainWindowTestAccess::modeled(window), &error)) << error.toStdString();
  ASSERT_TRUE(MainWindowTestAccess::import(window, example("examples/Multibed/CO2-N2-polishing/simulation.json"), &error));
  MainWindowTestAccess::components(window)[1].name = "Renamed CO2";
  const auto output = MainWindowTestAccess::output(window);
  EXPECT_TRUE(output.value("Adsorbents").toArray()[0].toObject().value("ComponentParameters").toObject().contains("Renamed CO2"));
  EXPECT_TRUE(LabInput::validate(output, &error)) << error.toStdString();
}

TEST(LabInput, FindsInstalledSiblingBeforeStaleBuildPath)
{
  QTemporaryDir directory;
  QFile executable(directory.filePath("ruptura"));
  ASSERT_TRUE(executable.open(QIODevice::WriteOnly));
  executable.write("#!/bin/sh\nexit 0\n");
  executable.close();
  ASSERT_TRUE(executable.setPermissions(QFileDevice::ReadOwner | QFileDevice::WriteOwner | QFileDevice::ExeOwner));
  EXPECT_EQ(LabInput::findCli(directory.path(), "/missing/build/ruptura"), executable.fileName());
}

#ifdef LAB_CLI_EXECUTABLE
TEST(LabInput, ExportedGasLiquidAndMultibedInputsRunInCli)
{
  for (const QString& path : {QString("examples/Liquid/pH-dependent/simulation.json"),
       QString("examples/Chemisorption-CO2/breakthrough/simulation.json"),
       QString("examples/Multibed/CO2-N2-polishing/simulation.json")})
  {
    SCOPED_TRACE(path.toStdString());
    MainWindow window;
    QString error;
    ASSERT_TRUE(MainWindowTestAccess::import(window, example(path), &error)) << error.toStdString();
    auto input = MainWindowTestAccess::modeled(window);
    input.insert("NumberOfTimeSteps", 2);
    input.insert("TimeStep", 1.0e-6);
    input.insert("PrintEvery", 1);
    input.insert("WriteEvery", 1);
    QTemporaryDir directory;
    writeJson(directory.filePath("simulation.json"), input);
    QProcess process;
    process.setWorkingDirectory(directory.path());
    process.start(QString::fromUtf8(LAB_CLI_EXECUTABLE), QStringList{});
    ASSERT_TRUE(process.waitForFinished(20000));
    EXPECT_EQ(process.exitStatus(), QProcess::NormalExit);
    EXPECT_EQ(process.exitCode(), 0) << process.readAllStandardError().toStdString();
    EXPECT_TRUE(QFile::exists(directory.filePath("analysis.ipynb")));
  }
}

TEST(LabInput, RunnerValidatesAndRunsActualExecutable)
{
  RupturaRunner runner;
  QString error;
  auto bad = mixture();
  bad.insert("SimulationType", "Fitting");
  EXPECT_FALSE(runner.startRun(bad, &error));
  EXPECT_FALSE(runner.isRunning());
  QTemporaryDir directory;
  const auto previous = QDir::currentPath();
  ASSERT_TRUE(QDir::setCurrent(directory.path()));
  QEventLoop loop;
  QTimer timer;
  timer.setSingleShot(true);
  QObject::connect(&runner, &RupturaRunner::runChanged, &loop, [&](const auto& record)
      { if (record.status == "completed" || record.status == "failed") loop.quit(); });
  QObject::connect(&timer, &QTimer::timeout, &loop, &QEventLoop::quit);
  const bool started = runner.startRun(mixture(), &error);
  if (started) { timer.start(15000); loop.exec(); }
  EXPECT_TRUE(QDir::setCurrent(previous));
  ASSERT_TRUE(started) << error.toStdString();
  EXPECT_EQ(runner.latestRun().status, "completed") << runner.logTail(runner.latestRun()).toStdString();
  EXPECT_EQ(runner.latestRun().exitCode, 0);
}
#endif
}

int main(int argc, char** argv)
{
  qputenv("QT_QPA_PLATFORM", "offscreen");
  QApplication application(argc, argv);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

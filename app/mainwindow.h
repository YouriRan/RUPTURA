#pragma once

#include <QJsonArray>
#include <QJsonObject>
#include <QList>
#include <QMainWindow>
#include <QString>
#include <QVariant>
#include <functional>
#include <optional>

#include "runner.h"

class QCheckBox;
class QComboBox;
class QFormLayout;
class QGroupBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QPlainTextEdit;
class QPushButton;
class QTableWidget;
class QVBoxLayout;

struct ParameterDefinition
{
  QString key;
  QString label;
  double defaultValue{0.0};
  bool integer{false};
  bool boolean{false};
};

struct PhysisorptionSiteConfig
{
  QString id;
  QString type{"Langmuir"};
  QVariantMap parameters;
};

struct ChemisorptionSiteConfig
{
  QString id;
  QString type{"FirstOrder"};
  QVariantMap parameters;
  PhysisorptionSiteConfig isotherm;
};

struct ComponentConfig
{
  QString id;
  QString name;
  bool carrier{false};
  std::optional<double> massTransfer;
  std::optional<double> axialDispersion;
  std::optional<double> molecularWeight;
  std::optional<double> heatOfAdsorption;
  std::optional<double> referenceTemperature;
  bool nonIsothermal{false};
  QList<PhysisorptionSiteConfig> physisorptionSites;
  QList<ChemisorptionSiteConfig> chemisorptionSites;
};

struct BedConfig
{
  QString id;
  QString name;
  double length{0.15};
  double mixFraction{0.5};
  double voidFraction{0.4};
  double particleDensity{1000.0};
  double particleDiameter{0.001};
  double interfaceLengthAfter{0.0};
};

struct ColumnConfig
{
  QString id;
  QString name;
  QString mode{"single"};
  QString boundaryCondition{"InletPressureInletVelocity"};
  double temperature{300.0};
  double dynamicViscosity{1.0e-5};
  double voidFraction{0.4};
  double particleDensity{1693.89};
  double particleDiameter{0.001};
  double inletPressure{2500000.0};
  std::optional<double> outletPressure;
  double pressureGradient{0.0};
  double velocity{0.1};
  double length{0.3};
  bool energyBalance{false};
  double influxTemperature{300.0};
  QString geometryType{"PackedBed"};
  QString channelShape{"square"};
  double internalChannelDimension{0.001};
  double internalDiameter{0.04};
  double outerDiameter{0.045};
  int numberOfChannels{10};
  double washcoatThickness{0.0001};
  std::optional<double> washcoatVolumePerChannelVolume;
  double wallDensity{7850.0};
  double gasThermalConductivity{0.09};
  double wallThermalConductivity{16.0};
  double heatTransferGasSolid{220.0};
  double heatTransferGasWall{20.0};
  double heatTransferWallExternal{10.0};
  double heatCapacityGas{5000.0};
  double heatCapacitySolid{900.0};
  double heatCapacityWall{400.0};
  QList<BedConfig> beds;
};

struct FeedConfig
{
  QString componentId;
  double gasFraction{0.0};
};

struct SwingPhaseConfig
{
  QString id;
  QString name;
  double temperature{300.0};
  double inletPressure{100000.0};
  int numberOfSteps{100000};
};

struct ReactionParticipantConfig
{
  QString componentId;
  double stoichiometry{1.0};
};

struct ReactionConfig
{
  QString id;
  QString name;
  QString phase{"Physisorbed"};
  QString style{"GeneralPowerLaw"};
  int site{0};
  QList<ReactionParticipantConfig> reactants;
  QList<ReactionParticipantConfig> products;
  double forwardRateCoefficient{1.0};
  double forwardActivationEnergy{0.0};
  double equilibriumConstant{1.0};
  double gibbsFreeEnergy{0.0};
  double rateLimitTime{1.0e-4};
};

struct SimulationConfig
{
  QString id;
  QString type{"Breakthrough"};
  QString displayName;
  QString mixtureMethod{"SIAST"};
  double pressureStart{0.0001};
  double pressureEnd{10000000.0};
  int pressurePoints{100};
  QString pressureScale{"log"};
  QString integrator{"CVODE"};
  bool autoSteps{true};
  int timeSteps{500000};
  int initSteps{0};
  int printEvery{10000};
  int writeEvery{10000};
  double timeStep{0.0005};
  int gridPoints{100};
  QString selectedColumnId{"column-1"};
  QList<FeedConfig> componentFeeds;
  QList<SwingPhaseConfig> swingPhases;
  RupturaRunner::RunRecord lastRun;
  QString lastRunLog;
};

class MainWindow : public QMainWindow
{
  Q_OBJECT

 public:
  explicit MainWindow(QWidget* parent = nullptr);

 private:
  RupturaRunner runner_;
  int idCounter_{10};
  bool rebuilding_{false};
  int selectedComponentIndex_{0};
  int selectedColumnIndex_{0};
  int selectedReactionIndex_{0};
  int selectedSimulationIndex_{0};

  QList<ComponentConfig> components_;
  QList<ColumnConfig> columns_;
  QList<ReactionConfig> reactions_;
  QList<SimulationConfig> simulations_;

  QLabel* runBadge_{nullptr};
  QLabel* activeColumnLabel_{nullptr};
  QLabel* reactionsMetric_{nullptr};
  QLabel* fractionTotalLabel_{nullptr};
  QLabel* runIdLabel_{nullptr};
  QListWidget* componentList_{nullptr};
  QListWidget* columnList_{nullptr};
  QListWidget* reactionList_{nullptr};
  QListWidget* simulationList_{nullptr};
  QVBoxLayout* componentEditorLayout_{nullptr};
  QVBoxLayout* columnEditorLayout_{nullptr};
  QVBoxLayout* reactionEditorLayout_{nullptr};
  QVBoxLayout* simulationEditorLayout_{nullptr};
  QTableWidget* feedTable_{nullptr};
  QPlainTextEdit* jsonPreview_{nullptr};
  QPlainTextEdit* runLog_{nullptr};
  QPushButton* runButton_{nullptr};
  QPushButton* runTopButton_{nullptr};
  QPushButton* cancelButton_{nullptr};
  QPushButton* cancelTopButton_{nullptr};
  QPushButton* analysisButton_{nullptr};
  QPushButton* analysisTopButton_{nullptr};
  QPushButton* loadButton_{nullptr};
  QPushButton* saveButton_{nullptr};
  QString pendingRunSimulationId_;

  void seedDefaults();
  void setupUi();
  void refreshAll();
  void refreshComponentList();
  void refreshColumnList();
  void refreshReactionList();
  void refreshSimulationList();
  void rebuildComponentEditor();
  void rebuildColumnEditor();
  void rebuildReactionEditor();
  void rebuildSimulationEditor();
  void updatePreview();
  void updateFractionTotal();
  void refreshRunPanel();
  void updateRunControls(const QString& status);
  void updateRunBadge(const QString& status);

  QString newId(const QString& prefix);
  PhysisorptionSiteConfig makePhySite(const QString& type = "Langmuir", const QVariantMap& values = {});
  ChemisorptionSiteConfig makeChemSite(const QString& type = "FirstOrder", const QVariantMap& values = {});
  BedConfig makeBed(const QString& name);
  ColumnConfig makeColumn(const QString& name, const QString& id = {});
  SwingPhaseConfig makeSwingPhase(const QString& name, double temperature, double inletPressure, int steps = 100000);
  ReactionConfig makeReaction(const QString& name = {});
  SimulationConfig makeSimulation(const QString& name, const QString& id = {}, const QList<FeedConfig>& feeds = {});
  void ensureSwingPhases(SimulationConfig& simulation);

  SimulationConfig* activeSimulation();
  const SimulationConfig* activeSimulation() const;
  ColumnConfig* selectedColumnFor(SimulationConfig& simulation);
  const ColumnConfig* selectedColumnFor(const SimulationConfig& simulation) const;
  SimulationConfig* simulationById(const QString& simulationId);
  SimulationConfig* simulationForRunId(const QString& runId);
  int feedIndex(const SimulationConfig& simulation, const QString& componentId) const;
  double feedTotal(const SimulationConfig& simulation) const;

  QJsonObject buildSimulationJson() const;
  QJsonObject buildStateJson() const;
  QList<QString> componentExportOrderFor(const SimulationConfig& simulation) const;
  QJsonArray componentsJsonFor(const SimulationConfig& simulation) const;
  QJsonArray reactionsJsonFor(const SimulationConfig& simulation) const;
  QJsonObject componentToJson(const ComponentConfig& component, double gasFraction) const;
  QJsonObject reactionToJson(const ReactionConfig& reaction, const SimulationConfig& simulation) const;
  QJsonObject phySiteToJson(const PhysisorptionSiteConfig& site) const;
  QJsonObject chemSiteToJson(const ChemisorptionSiteConfig& site) const;
  QJsonArray swingPhasesToJson(const SimulationConfig& simulation) const;
  void applyColumnToBreakthrough(QJsonObject& simulationJson, const ColumnConfig& column,
                                 const SimulationConfig& simulation) const;
  bool loadStateJson(const QJsonObject& state, QString* error);
  bool importSimulationJson(const QJsonObject& simulationJson, QString* error);
  void resetIdCounterFromState();

  void loadJsonFile();
  void saveStateFile();
  void addComponent();
  void removeComponent();
  void addColumn();
  void removeColumn();
  void addReaction();
  void removeReaction();
  void addSimulation();
  void copySimulation();
  void removeSimulation();
  void runSimulation();
  void cancelSimulation();
  void openAnalysis();

  QWidget* makePanel(const QString& title, QLabel* metric, QWidget* controls, QListWidget* list,
                     QVBoxLayout** editorLayout);
  QLineEdit* addTextField(QFormLayout* form, const QString& label, const QString& value,
                          const std::function<void(const QString&)>& setter);
  QLineEdit* addNumberField(QFormLayout* form, const QString& label, double value,
                            const std::function<void(double)>& setter);
  QLineEdit* addOptionalNumberField(QFormLayout* form, const QString& label, const std::optional<double>& value,
                                    const std::function<void(std::optional<double>)>& setter);
  QCheckBox* addCheckField(QFormLayout* form, const QString& label, bool checked,
                           const std::function<void(bool)>& setter);
  QComboBox* addComboField(QFormLayout* form, const QString& label, const QString& value,
                           const QList<QPair<QString, QString>>& options,
                           const std::function<void(const QString&)>& setter);
  void addParameterFields(QFormLayout* form, QVariantMap& parameters, const QList<ParameterDefinition>& definitions);
};

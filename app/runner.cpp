#include "runner.h"

#include <QCoreApplication>
#include <QDesktopServices>
#include <QDir>
#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QPointer>
#include <QRegularExpression>
#include <QStandardPaths>
#include <QTcpServer>
#include <QTimer>
#include <QUrl>
#include <algorithm>

#ifndef RUPTURA_CLI_EXECUTABLE
#define RUPTURA_CLI_EXECUTABLE "ruptura"
#endif

#ifndef RUPTURA_SOURCE_DIR
#define RUPTURA_SOURCE_DIR ""
#endif

namespace
{
QString isoDate(const QDateTime& value)
{
  return value.isValid() ? value.toUTC().toString(Qt::ISODateWithMs) : QString{};
}

QJsonArray notebookSource(const QString& source)
{
  QJsonArray lines;
  const QStringList split = source.split('\n');
  for (qsizetype i = 0; i < split.size(); ++i)
  {
    QString line = split.at(i);
    if (i + 1 < split.size())
    {
      line += '\n';
    }
    lines.append(line);
  }
  return lines;
}

QJsonObject markdownCell(const QString& source)
{
  return QJsonObject{
      {"cell_type", "markdown"},
      {"metadata", QJsonObject{}},
      {"source", notebookSource(source)},
  };
}

QJsonObject codeCell(const QString& source)
{
  return QJsonObject{
      {"cell_type", "code"},
      {"execution_count", QJsonValue(QJsonValue::Null)},
      {"metadata", QJsonObject{}},
      {"outputs", QJsonArray{}},
      {"source", notebookSource(source)},
  };
}

QJsonObject runRecordJson(const RupturaRunner::RunRecord& record)
{
  return QJsonObject{
      {"id", record.id},
      {"directory", record.directory},
      {"simulationPath", record.simulationPath},
      {"notebookPath", record.notebookPath},
      {"logPath", record.logPath},
      {"status", record.status},
      {"startedAt", isoDate(record.startedAt)},
      {"finishedAt", isoDate(record.finishedAt)},
      {"exitCode", record.exitCode},
      {"errorMessage", record.errorMessage},
  };
}
}  // namespace

RupturaRunner::RupturaRunner(QObject* parent) : QObject(parent)
{
  qRegisterMetaType<RupturaRunner::RunRecord>("RupturaRunner::RunRecord");
  repoRoot_ = QString::fromUtf8(RUPTURA_SOURCE_DIR);
  cliExecutable_ = QString::fromUtf8(RUPTURA_CLI_EXECUTABLE);

  if (repoRoot_.isEmpty())
  {
    repoRoot_ = QDir::currentPath();
  }
  if (cliExecutable_.isEmpty())
  {
    cliExecutable_ = QStandardPaths::findExecutable("ruptura");
  }
}

bool RupturaRunner::isRunning() const { return runProcess_ != nullptr && runProcess_->state() != QProcess::NotRunning; }

QString RupturaRunner::latestRunId() const { return latestRunId_; }

RupturaRunner::RunRecord RupturaRunner::latestRun() const
{
  return latestRunId_.isEmpty() || !runs_.contains(latestRunId_) ? RunRecord{} : runs_.value(latestRunId_);
}

RupturaRunner::RunRecord RupturaRunner::run(const QString& runId) const
{
  return runId.isEmpty() || !runs_.contains(runId) ? RunRecord{} : runs_.value(runId);
}

void RupturaRunner::rememberRun(const RunRecord& record)
{
  if (record.id.isEmpty())
  {
    return;
  }
  runs_.insert(record.id, record);
  if (!runOrder_.contains(record.id))
  {
    runOrder_.append(record.id);
  }
  if (latestRunId_.isEmpty())
  {
    latestRunId_ = record.id;
  }
}

bool RupturaRunner::startRun(const QJsonObject& simulation, QString* errorMessage)
{
  if (isRunning())
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "A simulation is already running.";
    }
    return false;
  }

  const QString displayName = simulation.value("DisplayName").toString("experiment");
  const QString runId = nextRunId(displayName);
  const QString runRoot = QDir(QDir::currentPath()).filePath("simulations");
  const QString runDirectory = QDir(runRoot).filePath(runId);

  if (!QDir().mkpath(runDirectory))
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "Could not create the simulation directory.";
    }
    return false;
  }

  RunRecord record;
  record.id = runId;
  record.directory = runDirectory;
  record.simulationPath = QDir(runDirectory).filePath("simulation.json");
  record.notebookPath = QDir(runDirectory).filePath("analysis.ipynb");
  record.logPath = QDir(runDirectory).filePath("run.log");
  record.status = "running";
  record.startedAt = QDateTime::currentDateTimeUtc();

  QFile simulationFile(record.simulationPath);
  if (!simulationFile.open(QIODevice::WriteOnly | QIODevice::Truncate))
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "Could not write simulation.json.";
    }
    return false;
  }
  simulationFile.write(QJsonDocument(simulation).toJson(QJsonDocument::Indented));
  simulationFile.close();

  QFile logFile(record.logPath);
  if (logFile.open(QIODevice::WriteOnly | QIODevice::Truncate))
  {
    logFile.write(QString("Ruptura CLI: %1\nWorking directory: %2\n\n").arg(cliExecutable_, runDirectory).toUtf8());
  }

  runs_.insert(runId, record);
  runOrder_.append(runId);
  latestRunId_ = runId;
  activeRunId_ = runId;
  writeManifest(record);
  emitRecord(record);

  auto* process = new QProcess(this);
  runProcess_ = process;
  process->setProgram(cliExecutable_);
  process->setWorkingDirectory(runDirectory);
  process->setProcessChannelMode(QProcess::MergedChannels);

  connect(process, &QProcess::readyReadStandardOutput, this,
          [this, process, runId]() { appendProcessOutput(runId, process->readAllStandardOutput()); });

  connect(process, &QProcess::errorOccurred, this,
          [this, process, runId](QProcess::ProcessError error)
          {
            RunRecord& run = runs_[runId];
            if (run.status == "canceling")
            {
              return;
            }
            run.errorMessage = QString("Process error %1: %2").arg(static_cast<int>(error)).arg(process->errorString());
            if (error == QProcess::FailedToStart)
            {
              run.status = "failed";
              run.finishedAt = QDateTime::currentDateTimeUtc();
              activeRunId_.clear();
              writeManifest(run);
              emitRecord(run);
              emit statusMessage(run.errorMessage);
            }
          });

  connect(process, qOverload<int, QProcess::ExitStatus>(&QProcess::finished), this,
          [this, process, runId](int exitCode, QProcess::ExitStatus exitStatus)
          {
            appendProcessOutput(runId, process->readAllStandardOutput());

            RunRecord& run = runs_[runId];
            const bool wasCanceled = run.status == "canceling" || run.status == "canceled";
            run.exitCode = exitCode;
            run.finishedAt = QDateTime::currentDateTimeUtc();
            if (wasCanceled)
            {
              run.status = "canceled";
              if (run.errorMessage.isEmpty() || run.errorMessage == "Cancellation requested by user.")
              {
                run.errorMessage = "Canceled by user.";
              }
            }
            else
            {
              run.status = (exitStatus == QProcess::NormalExit && exitCode == 0) ? "completed" : "failed";
              if (run.status == "failed" && run.errorMessage.isEmpty())
              {
                run.errorMessage = QString("ruptura exited with code %1").arg(exitCode);
              }
            }
            writeManifest(run);

            if (runProcess_ == process)
            {
              runProcess_ = nullptr;
            }
            activeRunId_.clear();
            emitRecord(run);
            emit statusMessage(QString("Simulation %1 %2").arg(runId, run.status));
            process->deleteLater();
          });

  process->start();
  if (!process->waitForStarted(1000))
  {
    RunRecord& run = runs_[runId];
    run.status = "failed";
    run.finishedAt = QDateTime::currentDateTimeUtc();
    run.errorMessage = process->errorString();
    writeManifest(run);
    runProcess_ = nullptr;
    activeRunId_.clear();
    process->deleteLater();
    emitRecord(run);
    if (errorMessage != nullptr)
    {
      *errorMessage = run.errorMessage;
    }
    return false;
  }

  emit statusMessage(QString("Simulation %1 running").arg(runId));
  return true;
}

bool RupturaRunner::cancelRun(QString* errorMessage)
{
  if (!isRunning() || activeRunId_.isEmpty() || !runs_.contains(activeRunId_))
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "No simulation is running.";
    }
    return false;
  }

  RunRecord& run = runs_[activeRunId_];
  run.status = "canceling";
  run.errorMessage = "Cancellation requested by user.";
  writeManifest(run);
  appendProcessOutput(run.id, "\nCancel requested by user.\n");
  emitRecord(run);
  emit statusMessage(QString("Canceling simulation %1").arg(run.id));

  QPointer<QProcess> process(runProcess_);
  process->terminate();
  QTimer::singleShot(1500, this,
                     [process]()
                     {
                       if (process != nullptr && process->state() != QProcess::NotRunning)
                       {
                         process->kill();
                       }
                     });
  return true;
}

bool RupturaRunner::openAnalysis(QString* errorMessage) { return openAnalysis(latestRunId_, errorMessage); }

bool RupturaRunner::openAnalysis(const QString& runId, QString* errorMessage)
{
  if (runId.isEmpty() || !runs_.contains(runId))
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "No simulation is available for analysis.";
    }
    return false;
  }

  RunRecord& record = runs_[runId];
  if (!writeAnalysisNotebook(record, errorMessage))
  {
    emitRecord(record);
    return false;
  }

  const QString url = startJupyter(record, errorMessage);
  emitRecord(record);
  return !url.isEmpty();
}

QString RupturaRunner::logTail(const RunRecord& record, qint64 maxBytes) const
{
  QFile file(record.logPath);
  if (!file.open(QIODevice::ReadOnly))
  {
    return {};
  }

  if (file.size() > maxBytes)
  {
    file.seek(file.size() - maxBytes);
  }
  return QString::fromUtf8(file.readAll());
}

bool RupturaRunner::writeAnalysisNotebook(RunRecord& record, QString* errorMessage)
{
  const QString escapedRepo = QString(repoRoot_).replace('\\', "\\\\").replace('\'', "\\'");
  const QString escapedRun = QString(record.directory).replace('\\', "\\\\").replace('\'', "\\'");

  QJsonArray cells;
  cells.append(markdownCell(QString("# Ruptura analysis\n\nSimulation `%1`").arg(record.id)));
  cells.append(codeCell(QString("from pathlib import Path\n"
                                "import json\n"
                                "import sys\n\n"
                                "repo_root = Path(r'%1').resolve()\n"
                                "run_dir = Path(r'%2').resolve()\n"
                                "if str(repo_root) not in sys.path:\n"
                                "    sys.path.insert(0, str(repo_root))\n\n"
                                "simulation_json = run_dir / 'simulation.json'\n"
                                "data_dir = run_dir\n\n"
                                "import ruptura\n"
                                "from IPython.display import display\n\n"
                                "simulation_config = json.loads(simulation_json.read_text())\n"
                                "simulation_type = simulation_config.get('SimulationType', 'Breakthrough')\n"
                                "simulation_type\n")
                            .arg(escapedRepo, escapedRun)));
  cells.append(
      codeCell("if simulation_type == 'MixturePrediction':\n"
               "    from ruptura.plot_mixture import MixturePredictionPlotly\n"
               "    plotter = MixturePredictionPlotly.from_simulation_json(simulation_json, data_dir=data_dir)\n"
               "    display(plotter.pure_components(include_carrier_gas=False))\n"
               "    display(plotter.mixture_loading(include_carrier_gas=False))\n"
               "    display(plotter.mixture_adsorbed_molfractions(include_carrier_gas=False))\n"
               "else:\n"
               "    from ruptura.plot_breakthrough import BreakthroughPlotly\n"
               "    plotter = BreakthroughPlotly.from_simulation_json(simulation_json, data_dir=data_dir)\n"
               "    display(plotter.breakthrough(\n"
               "        x_units='min',\n"
               "        y_units='normalized concentration',\n"
               "        include_carrier_gas=False,\n"
               "        show_markers=True,\n"
               "    ))\n"));
  cells.append(
      codeCell("if simulation_type in ('Breakthrough', 'SwingAdsorption'):\n"
               "    try:\n"
               "        display(plotter.temperature_triplet_time(x_units='min', y_units='kelvin', show_markers=True))\n"
               "    except Exception as exc:\n"
               "        print(f'Temperature plot unavailable: {exc}')\n"
               "else:\n"
               "    print('Temperature widgets are specific to column runs.')\n"));
  cells.append(
      codeCell("if simulation_type in ('Breakthrough', 'SwingAdsorption'):\n"
               "    explorer = plotter.explorer()\n"
               "    explorer.display()\n"
               "else:\n"
               "    print('Mixture prediction figures are loaded above.')\n"));
  cells.append(codeCell(
      "if simulation_type == 'SwingAdsorption':\n"
      "    print('Swing adsorption writes file output for plotting; in-memory compute is not available yet.')\n"
      "else:\n"
      "    try:\n"
      "        result = ruptura.run(simulation_json)\n"
      "        print(result.kind, result.shape)\n"
      "    except Exception as exc:\n"
      "        print(f'In-memory ruptura compute unavailable: {exc}')\n"));

  const QJsonObject metadata{
      {"kernelspec", QJsonObject{{"display_name", "Python 3"}, {"language", "python"}, {"name", "python3"}}},
      {"language_info", QJsonObject{{"name", "python"}, {"pygments_lexer", "ipython3"}}},
  };
  const QJsonObject notebook{
      {"cells", cells},
      {"metadata", metadata},
      {"nbformat", 4},
      {"nbformat_minor", 5},
  };

  QFile file(record.notebookPath);
  if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate))
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "Could not write the analysis notebook.";
    }
    return false;
  }
  file.write(QJsonDocument(notebook).toJson(QJsonDocument::Indented));
  file.close();

  writeManifest(record);
  return true;
}

QString RupturaRunner::startJupyter(const RunRecord& record, QString* errorMessage)
{
  stopJupyter();

  QString program = QStandardPaths::findExecutable("jupyter");
  QStringList args;
  if (!program.isEmpty())
  {
    args << "lab";
  }
  else
  {
    program = QStandardPaths::findExecutable("python3");
    if (program.isEmpty())
    {
      program = QStandardPaths::findExecutable("python");
    }
    if (program.isEmpty())
    {
      if (errorMessage != nullptr)
      {
        *errorMessage = "Notebook generated, but neither jupyter nor python was found on PATH.";
      }
      return {};
    }
    args << "-m" << "jupyter" << "lab";
  }

  const quint16 port = reservePort();
  if (port == 0)
  {
    if (errorMessage != nullptr)
    {
      *errorMessage = "Notebook generated, but no local port was available for Jupyter.";
    }
    return {};
  }

  args << "--no-browser"
       << "--notebook-dir" << record.directory << "--port" << QString::number(port) << "--ServerApp.token="
       << "--NotebookApp.token="
       << "--ServerApp.password="
       << "--NotebookApp.password=";

  jupyterProcess_ = new QProcess(this);
  jupyterProcess_->setProgram(program);
  jupyterProcess_->setArguments(args);
  jupyterProcess_->setWorkingDirectory(record.directory);
  jupyterProcess_->setProcessChannelMode(QProcess::MergedChannels);

  connect(jupyterProcess_, &QProcess::readyReadStandardOutput, this,
          [this]()
          {
            const QString output = QString::fromUtf8(jupyterProcess_->readAllStandardOutput()).trimmed();
            if (!output.isEmpty())
            {
              emit statusMessage(output);
            }
          });

  connect(jupyterProcess_, qOverload<int, QProcess::ExitStatus>(&QProcess::finished), this,
          [this](int, QProcess::ExitStatus)
          {
            jupyterUrl_.clear();
            if (jupyterProcess_ != nullptr)
            {
              jupyterProcess_->deleteLater();
              jupyterProcess_ = nullptr;
            }
          });

  jupyterProcess_->start();
  if (!jupyterProcess_->waitForStarted(2000))
  {
    if (errorMessage != nullptr)
    {
      *errorMessage =
          QString("Notebook generated, but Jupyter could not start: %1").arg(jupyterProcess_->errorString());
    }
    jupyterProcess_->deleteLater();
    jupyterProcess_ = nullptr;
    return {};
  }

  jupyterUrl_ = QString("http://127.0.0.1:%1/lab/tree/analysis.ipynb").arg(port);
  QDesktopServices::openUrl(QUrl(jupyterUrl_));
  emit statusMessage(QString("Opened analysis notebook for %1").arg(record.id));
  return jupyterUrl_;
}

void RupturaRunner::stopJupyter()
{
  if (jupyterProcess_ == nullptr)
  {
    return;
  }

  if (jupyterProcess_->state() != QProcess::NotRunning)
  {
    jupyterProcess_->terminate();
    if (!jupyterProcess_->waitForFinished(1500))
    {
      jupyterProcess_->kill();
      jupyterProcess_->waitForFinished(1500);
    }
  }
  jupyterProcess_->deleteLater();
  jupyterProcess_ = nullptr;
  jupyterUrl_.clear();
}

QString RupturaRunner::nextRunId(const QString& displayName) const
{
  const QString stamp = QDateTime::currentDateTimeUtc().toString("yyyyMMdd-HHmmss-zzz");
  return stamp + "-" + safeFileName(displayName, "simulation").left(36);
}

QString RupturaRunner::safeFileName(const QString& value, const QString& fallback) const
{
  QString safe = value.trimmed();
  safe.replace(QRegularExpression("[^A-Za-z0-9._-]+"), "-");
  safe.replace(QRegularExpression("-+"), "-");
  safe = safe.mid(0, 64).trimmed();
  if (safe.isEmpty() || safe == "-" || safe == ".")
  {
    return fallback;
  }
  return safe;
}

quint16 RupturaRunner::reservePort() const
{
  QTcpServer portServer;
  if (!portServer.listen(QHostAddress::LocalHost, 0))
  {
    return 0;
  }
  const quint16 port = portServer.serverPort();
  portServer.close();
  return port;
}

void RupturaRunner::appendProcessOutput(const QString& runId, const QByteArray& output)
{
  if (output.isEmpty() || !runs_.contains(runId))
  {
    return;
  }

  QFile logFile(runs_[runId].logPath);
  if (logFile.open(QIODevice::WriteOnly | QIODevice::Append))
  {
    logFile.write(output);
  }
  emit logChanged(runId, logTail(runs_[runId]));
}

void RupturaRunner::writeManifest(const RunRecord& record) const
{
  QFile file(QDir(record.directory).filePath("manifest.json"));
  if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate))
  {
    return;
  }

  QJsonObject manifest = runRecordJson(record);
  manifest.insert("repoRoot", repoRoot_);
  manifest.insert("cliExecutable", cliExecutable_);
  file.write(QJsonDocument(manifest).toJson(QJsonDocument::Indented));
}

void RupturaRunner::emitRecord(const RunRecord& record)
{
  emit runChanged(record);
  emit logChanged(record.id, logTail(record));
}

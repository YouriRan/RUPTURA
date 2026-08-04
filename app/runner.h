#pragma once

#include <QDateTime>
#include <QHash>
#include <QJsonObject>
#include <QObject>
#include <QProcess>
#include <QString>
#include <QStringList>

class RupturaRunner : public QObject
{
  Q_OBJECT

 public:
  struct RunRecord
  {
    QString id;
    QString directory;
    QString simulationPath;
    QString notebookPath;
    QString logPath;
    QString status{"created"};
    QString errorMessage;
    QDateTime startedAt;
    QDateTime finishedAt;
    int exitCode{-1};
  };

  explicit RupturaRunner(QObject* parent = nullptr);

  bool isRunning() const;
  QString latestRunId() const;
  RunRecord latestRun() const;
  RunRecord run(const QString& runId) const;
  QString logTail(const RunRecord& record, qint64 maxBytes = 12000) const;
  void rememberRun(const RunRecord& record);

  bool startRun(const QJsonObject& simulation, QString* errorMessage = nullptr);
  bool cancelRun(QString* errorMessage = nullptr);
  bool openAnalysis(QString* errorMessage = nullptr);
  bool openAnalysis(const QString& runId, QString* errorMessage = nullptr);

 signals:
  void runChanged(const RupturaRunner::RunRecord& record);
  void logChanged(const QString& runId, const QString& logText);
  void statusMessage(const QString& message);

 private:
  QHash<QString, RunRecord> runs_;
  QStringList runOrder_;
  QProcess* runProcess_{nullptr};
  QProcess* jupyterProcess_{nullptr};
  QString activeRunId_;
  QString latestRunId_;
  QString repoRoot_;
  QString cliExecutable_;
  QString jupyterUrl_;

  bool writeAnalysisNotebook(RunRecord& record, QString* errorMessage);
  QString startJupyter(const RunRecord& record, QString* errorMessage);
  void stopJupyter();

  QString nextRunId(const QString& displayName) const;
  QString safeFileName(const QString& value, const QString& fallback) const;
  quint16 reservePort() const;
  void appendProcessOutput(const QString& runId, const QByteArray& output);
  void writeManifest(const RunRecord& record) const;
  void emitRecord(const RunRecord& record);
};

Q_DECLARE_METATYPE(RupturaRunner::RunRecord)

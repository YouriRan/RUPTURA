#pragma once

#include <QJsonObject>
#include <QHash>
#include <QString>
#include <memory>

struct InputReader;

namespace LabInput
{
// Normalize CLI field spelling, without changing component/adsorbent names.
QJsonObject canonicalize(const QJsonObject& input);
QJsonObject resolveFiles(const QJsonObject& input, const QString& directory);
QJsonObject applyChanges(const QJsonObject& original, const QJsonObject& before, const QJsonObject& after);
QJsonObject renameComponentReferences(const QJsonObject& input, const QHash<QString, QString>& names);
bool supported(const QJsonObject& input, QString* error);
bool validate(const QJsonObject& input, QString* error);
std::unique_ptr<InputReader> read(const QJsonObject& input, QString* error);
QString findCli(const QString& applicationDirectory, const QString& buildExecutable);
}

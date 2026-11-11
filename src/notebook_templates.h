#pragma once

#include <filesystem>
#include <string>
#include <string_view>

#include "json.h"

namespace RupturaNotebooks
{
using Notebook = nlohmann::ordered_json;

enum class SimulationType
{
  Breakthrough,
  MixturePrediction,
  Fitting,
};

enum class WriteResult
{
  Created,
  AlreadyExists,
};

SimulationType simulationTypeFromName(std::string_view name);
// The notebooks hold no paths: they read the results from the folder they are opened in (the run directory) and
// import the installed ruptura package.
Notebook defaultMixturePredictionNotebook(std::string_view runId);
Notebook defaultBreakthroughNotebook(std::string_view runId);
Notebook defaultFittingNotebook(std::string_view runId);
Notebook defaultNotebookForSimulation(SimulationType simulationType, std::string_view runId);
WriteResult writeDefaultNotebook(const std::filesystem::path& notebookPath, SimulationType simulationType,
                                 std::string_view runId);
}  // namespace RupturaNotebooks

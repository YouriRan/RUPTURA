#pragma once

#include <filesystem>
#include <string>
#include <string_view>

#include "inputreader.h"
#include "json.h"

namespace RupturaNotebooks
{
using Notebook = nlohmann::ordered_json;

enum class WriteResult
{
  Created,
  AlreadyExists,
};

// The notebooks hold no paths: they read the results from the folder they are opened in (the run directory) and
// import the installed ruptura package.
Notebook defaultMixturePredictionNotebook(std::string_view runId);
Notebook defaultBreakthroughNotebook(std::string_view runId);
Notebook defaultFittingNotebook(std::string_view runId);
/// The notebook of a SwingAdsorption run is that of a Breakthrough run.
Notebook defaultNotebookForSimulation(InputReader::SimulationType simulationType, std::string_view runId);
WriteResult writeDefaultNotebook(const std::filesystem::path& notebookPath, InputReader::SimulationType simulationType,
                                 std::string_view runId);
}  // namespace RupturaNotebooks

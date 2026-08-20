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
Notebook defaultMixturePredictionNotebook(std::string_view runId, const std::filesystem::path& repoRoot,
                                          const std::filesystem::path& runDirectory);
Notebook defaultBreakthroughNotebook(std::string_view runId, const std::filesystem::path& repoRoot,
                                     const std::filesystem::path& runDirectory);
Notebook defaultFittingNotebook(std::string_view runId, const std::filesystem::path& repoRoot,
                                const std::filesystem::path& runDirectory);
Notebook defaultNotebookForSimulation(SimulationType simulationType, std::string_view runId,
                                      const std::filesystem::path& repoRoot, const std::filesystem::path& runDirectory);
WriteResult writeDefaultNotebook(const std::filesystem::path& notebookPath, SimulationType simulationType,
                                 std::string_view runId, const std::filesystem::path& repoRoot,
                                 const std::filesystem::path& runDirectory);
}  // namespace RupturaNotebooks

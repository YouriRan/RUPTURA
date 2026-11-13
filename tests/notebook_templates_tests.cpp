#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "json.h"
#include "notebook_templates.h"

namespace
{
using Notebook = RupturaNotebooks::Notebook;

std::string notebookSource(const Notebook& notebook)
{
  std::string source;
  for (const Notebook& cell : notebook.at("cells"))
  {
    for (const Notebook& line : cell.at("source"))
    {
      source += line.get<std::string>();
    }
  }
  return source;
}

void expectValidNotebook(const Notebook& notebook)
{
  EXPECT_EQ(notebook.at("nbformat"), 4);
  EXPECT_EQ(notebook.at("nbformat_minor"), 5);
  EXPECT_FALSE(notebook.at("cells").empty());
  EXPECT_EQ(notebook.at("metadata").at("kernelspec").at("name"), "python3");
}

class TemporaryDirectory
{
 public:
  TemporaryDirectory()
  {
    static std::atomic<unsigned long long> sequence{0};
    const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
    path_ = std::filesystem::temp_directory_path() /
            ("ruptura-notebook-tests-" + std::to_string(timestamp) + "-" + std::to_string(sequence++));
    std::filesystem::create_directory(path_);
  }

  ~TemporaryDirectory()
  {
    std::error_code error;
    std::filesystem::remove_all(path_, error);
  }

  const std::filesystem::path& path() const { return path_; }

 private:
  std::filesystem::path path_;
};

std::string readFile(const std::filesystem::path& path)
{
  std::ifstream file(path, std::ios::binary);
  std::ostringstream contents;
  contents << file.rdbuf();
  return contents.str();
}
}  // namespace

TEST(NotebookTemplates, SelectsMixturePredictionTemplate)
{
  const Notebook notebook = RupturaNotebooks::defaultNotebookForSimulation(
      RupturaNotebooks::SimulationType::MixturePrediction, "run-1", "/repo", "/runs/run-1");
  const std::string source = notebookSource(notebook);

  expectValidNotebook(notebook);
  EXPECT_NE(source.find("MixturePredictionPlotly"), std::string::npos);
  EXPECT_NE(source.find("mixture_loading"), std::string::npos);
  EXPECT_EQ(source.find("BreakthroughPlotly"), std::string::npos);
}

TEST(NotebookTemplates, SelectsBreakthroughTemplate)
{
  const Notebook notebook = RupturaNotebooks::defaultNotebookForSimulation(
      RupturaNotebooks::SimulationType::Breakthrough, "run-2", "/repo", "/runs/run-2");
  const std::string source = notebookSource(notebook);

  expectValidNotebook(notebook);
  EXPECT_NE(source.find("BreakthroughPlotly"), std::string::npos);
  EXPECT_NE(source.find("plotter.explorer()"), std::string::npos);
  EXPECT_EQ(source.find("FittingPlotly"), std::string::npos);
}

TEST(NotebookTemplates, SelectsFittingTemplate)
{
  const Notebook notebook = RupturaNotebooks::defaultNotebookForSimulation(RupturaNotebooks::SimulationType::Fitting,
                                                                           "run-3", "/repo", "/runs/run-3");
  const std::string source = notebookSource(notebook);

  expectValidNotebook(notebook);
  EXPECT_NE(source.find("FittingPlotly"), std::string::npos);
  EXPECT_NE(source.find("_fitted.json"), std::string::npos);
  EXPECT_NE(source.find("all_fit_figure"), std::string::npos);
  EXPECT_EQ(source.find("BreakthroughPlotly"), std::string::npos);
}

TEST(NotebookTemplates, MapsSwingAdsorptionToBreakthroughTemplate)
{
  EXPECT_EQ(RupturaNotebooks::simulationTypeFromName("SwingAdsorption"),
            RupturaNotebooks::SimulationType::Breakthrough);
}

TEST(NotebookTemplates, WritesSelectedNotebookToRunDirectory)
{
  TemporaryDirectory runDirectory;
  const std::filesystem::path notebookPath = runDirectory.path() / "analysis.ipynb";

  EXPECT_EQ(RupturaNotebooks::writeDefaultNotebook(notebookPath, RupturaNotebooks::SimulationType::Fitting, "run-5",
                                                   "/repo", runDirectory.path()),
            RupturaNotebooks::WriteResult::Created);
  EXPECT_NE(notebookSource(Notebook::parse(readFile(notebookPath))).find("FittingPlotly"), std::string::npos);
}

TEST(NotebookTemplates, NeverOverwritesAnExistingNotebook)
{
  TemporaryDirectory runDirectory;
  const std::filesystem::path notebookPath = runDirectory.path() / "analysis.ipynb";
  const std::string userNotebook = "user-edited notebook\n";
  {
    std::ofstream file(notebookPath, std::ios::binary);
    file << userNotebook;
  }

  EXPECT_EQ(RupturaNotebooks::writeDefaultNotebook(notebookPath, RupturaNotebooks::SimulationType::Breakthrough,
                                                   "run-6", "/repo", runDirectory.path()),
            RupturaNotebooks::WriteResult::AlreadyExists);
  EXPECT_EQ(readFile(notebookPath), userNotebook);
}

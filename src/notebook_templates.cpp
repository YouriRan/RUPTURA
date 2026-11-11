#include "notebook_templates.h"

#include <fstream>
#include <stdexcept>
#include <vector>

namespace
{
using RupturaNotebooks::Notebook;

std::vector<std::string> notebookSource(const std::string& source)
{
  std::vector<std::string> lines;
  size_t start = 0;
  while (start < source.size())
  {
    const size_t newline = source.find('\n', start);
    if (newline == std::string::npos)
    {
      lines.push_back(source.substr(start));
      return lines;
    }
    lines.push_back(source.substr(start, newline - start + 1));
    start = newline + 1;
  }
  if (source.empty() || source.back() == '\n')
  {
    lines.emplace_back();
  }
  return lines;
}

Notebook markdownCell(const std::string& source)
{
  return Notebook{
      {"cell_type", "markdown"},
      {"metadata", Notebook::object()},
      {"source", notebookSource(source)},
  };
}

Notebook codeCell(const std::string& source)
{
  return Notebook{
      {"cell_type", "code"},          {"execution_count", nullptr},       {"metadata", Notebook::object()},
      {"outputs", Notebook::array()}, {"source", notebookSource(source)},
  };
}

Notebook notebook(const Notebook& cells)
{
  return Notebook{
      {"cells", cells},
      {"metadata",
       Notebook{{"kernelspec", Notebook{{"display_name", "Python 3"}, {"language", "python"}, {"name", "python3"}}},
                {"language_info", Notebook{{"name", "python"}, {"pygments_lexer", "ipython3"}}}}},
      {"nbformat", 4},
      {"nbformat_minor", 5},
  };
}

// The notebook holds no paths: Jupyter starts its kernel in the folder of the notebook, which is the run directory,
// and the plotting code comes from the installed ruptura package. A run directory can then be moved or shared (B-PY-4).
Notebook setupCell()
{
  return codeCell(
      "from pathlib import Path\n"
      "import json\n\n"
      "run_dir = Path.cwd()  # the folder of this notebook, which holds the simulation results\n"
      "simulation_json = run_dir / 'simulation.json'\n"
      "data_dir = run_dir\n"
      "simulation_config = json.loads(simulation_json.read_text())\n\n"
      "from IPython.display import display\n"
      "simulation_config.get('SimulationType')\n");
}

Notebook initialCells(const std::string& title, const std::string& description, std::string_view runId)
{
  Notebook cells = Notebook::array();
  cells.push_back(markdownCell("# " + title + "\n\n" + description + "\n\nSimulation `" + std::string(runId) + "`"));
  cells.push_back(setupCell());
  return cells;
}
}  // namespace

namespace RupturaNotebooks
{
SimulationType simulationTypeFromName(std::string_view name)
{
  if (name == "MixturePrediction") return SimulationType::MixturePrediction;
  if (name == "Fitting") return SimulationType::Fitting;
  return SimulationType::Breakthrough;
}

Notebook defaultMixturePredictionNotebook(std::string_view runId)
{
  Notebook cells = initialCells("Mixture prediction analysis",
                                "Run the cells after the simulation has written its component data files.", runId);
  cells.push_back(markdownCell("## Load the mixture-prediction results"));
  cells.push_back(
      codeCell("from ruptura.plot_mixture import MixturePredictionPlotly\n\n"
               "plotter = MixturePredictionPlotly.from_simulation_json(\n"
               "    simulation_json,\n"
               "    data_dir=data_dir,\n"
               ")\n"
               "plotter\n"));
  cells.push_back(markdownCell("## Pure-component isotherms"));
  cells.push_back(codeCell("display(plotter.pure_components(include_carrier_gas=False))\n"));
  cells.push_back(markdownCell("## Mixture loadings"));
  cells.push_back(codeCell("display(plotter.mixture_loading(include_carrier_gas=False))\n"));
  cells.push_back(markdownCell("## Adsorbed-phase mole fractions"));
  cells.push_back(codeCell("display(plotter.mixture_adsorbed_molfractions(include_carrier_gas=False))\n"));
  cells.push_back(markdownCell("## Chemisorption and total loadings"));
  cells.push_back(
      codeCell("if plotter.has_chemisorption():\n"
               "    display(plotter.mixture_chemisorption_loading(include_carrier_gas=False))\n"
               "    display(plotter.mixture_total_loading(include_carrier_gas=False))\n"
               "else:\n"
               "    print('This simulation has no chemisorption equilibrium isotherm.')\n"));
  return notebook(cells);
}

Notebook defaultBreakthroughNotebook(std::string_view runId)
{
  Notebook cells = initialCells("Breakthrough analysis",
                                "Run the cells after the simulation has written `column.data` and its component data "
                                "files.",
                                runId);
  cells.push_back(markdownCell("## Load the breakthrough results"));
  cells.push_back(
      codeCell("from ruptura.plot_breakthrough import BreakthroughPlotly\n\n"
               "plotter = BreakthroughPlotly.from_simulation_json(\n"
               "    simulation_json,\n"
               "    data_dir=data_dir,\n"
               ")\n"
               "plotter\n"));
  cells.push_back(markdownCell("## Breakthrough curves"));
  cells.push_back(
      codeCell("display(plotter.breakthrough(\n"
               "    x_units='min',\n"
               "    y_units='normalized_concentration',\n"
               "    include_carrier_gas=False,\n"
               "    show_markers=True,\n"
               "))\n"));
  cells.push_back(markdownCell("## Temperature history"));
  cells.push_back(
      codeCell("try:\n"
               "    display(plotter.temperature_triplet_time(\n"
               "        x_units='min',\n"
               "        y_units='kelvin',\n"
               "        show_markers=True,\n"
               "    ))\n"
               "except Exception as exc:\n"
               "    print(f'Temperature plot unavailable: {exc}')\n"));
  cells.push_back(markdownCell("## Interactive column explorer"));
  cells.push_back(
      codeCell("explorer = plotter.explorer()\n"
               "explorer.display()\n"));
  return notebook(cells);
}

Notebook defaultFittingNotebook(std::string_view runId)
{
  Notebook cells =
      initialCells("Fitting analysis", "Run the cells after fitting has written the fitted-components JSON file.",
                   runId);
  cells.push_back(markdownCell("## Load the fitting results"));
  cells.push_back(
      codeCell("from ruptura.plot_fitting import FittingPlotly\n\n"
               "display_name = simulation_config.get('DisplayName', '')\n"
               "fitted_json = data_dir / (f'{display_name}_fitted.json' if display_name else "
               "'fitted.json')\n"
               "plotter = FittingPlotly.from_simulation_json(\n"
               "    simulation_json,\n"
               "    fit_json=fitted_json,\n"
               "    data_dir=data_dir,\n"
               ")\n"
               "plotter\n"));
  cells.push_back(markdownCell("## Raw data and fitted isotherms"));
  cells.push_back(codeCell("display(plotter.all_fit_figure(show_fit=True))\n"));
  cells.push_back(markdownCell("Set `show_fit=False` to inspect the raw fitting data by itself."));
  cells.push_back(codeCell("display(plotter.all_fit_figure(show_fit=False))\n"));
  return notebook(cells);
}

Notebook defaultNotebookForSimulation(SimulationType simulationType, std::string_view runId)
{
  switch (simulationType)
  {
    case SimulationType::MixturePrediction:
      return defaultMixturePredictionNotebook(runId);
    case SimulationType::Fitting:
      return defaultFittingNotebook(runId);
    case SimulationType::Breakthrough:
    default:
      return defaultBreakthroughNotebook(runId);
  }
}

WriteResult writeDefaultNotebook(const std::filesystem::path& notebookPath, SimulationType simulationType,
                                 std::string_view runId)
{
  std::error_code existsError;
  if (std::filesystem::exists(notebookPath, existsError)) return WriteResult::AlreadyExists;
  if (existsError)
  {
    throw std::runtime_error("Could not inspect the default analysis notebook path '" + notebookPath.string() +
                             "': " + existsError.message());
  }

  std::ofstream file(notebookPath, std::ios::out | std::ios::binary | std::ios::noreplace);
  if (!file)
  {
    if (std::filesystem::exists(notebookPath)) return WriteResult::AlreadyExists;
    throw std::runtime_error("Could not create the default analysis notebook '" + notebookPath.string() + "'");
  }

  file << defaultNotebookForSimulation(simulationType, runId).dump(2) << '\n';
  file.close();
  if (!file)
  {
    throw std::runtime_error("Could not finish writing the default analysis notebook '" + notebookPath.string() + "'");
  }
  return WriteResult::Created;
}
}  // namespace RupturaNotebooks

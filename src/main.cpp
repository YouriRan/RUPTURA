#include <exception>
#include <filesystem>
#include <print>

#include "breakthrough.h"
#include "fitting.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "notebook_templates.h"
#include "special_functions.h"
#include "swing_adsorption.h"

#ifndef RUPTURA_SOURCE_DIR
#define RUPTURA_SOURCE_DIR ""
#endif

namespace
{
template <typename Simulation>
void runSimulation(const InputReader& reader)
{
  Simulation simulation(reader);
  simulation.print();
  simulation.run();
}

RupturaNotebooks::SimulationType notebookType(InputReader::SimulationType simulationType)
{
  switch (simulationType)
  {
    case InputReader::SimulationType::MixturePrediction:
      return RupturaNotebooks::SimulationType::MixturePrediction;
    case InputReader::SimulationType::Fitting:
      return RupturaNotebooks::SimulationType::Fitting;
    case InputReader::SimulationType::SwingAdsorption:
    case InputReader::SimulationType::Breakthrough:
    default:
      return RupturaNotebooks::SimulationType::Breakthrough;
  }
}

void ensureDefaultNotebook(InputReader::SimulationType simulationType)
{
  const std::filesystem::path runDirectory = std::filesystem::current_path();
  const std::filesystem::path notebookPath = runDirectory / "analysis.ipynb";
  const std::string runId = runDirectory.filename().string();
  const auto result = RupturaNotebooks::writeDefaultNotebook(notebookPath, notebookType(simulationType), runId,
                                                             RUPTURA_SOURCE_DIR, runDirectory);
  if (result == RupturaNotebooks::WriteResult::Created)
  {
    std::println("Wrote default analysis notebook: {}", notebookPath.string());
  }
}
}  // namespace

int main(void)
{
  try
  {
    InputReader reader("simulation.json");
    ensureDefaultNotebook(reader.simulationType);

    switch (reader.simulationType)
    {
      case InputReader::SimulationType::Breakthrough:
      default:
      {
        if (reader.isMultibed())
        {
          runSimulation<Breakthrough<MultibedColumn>>(reader);
        }
        else
        {
          runSimulation<Breakthrough<Column>>(reader);
        }
        break;
      }
      case InputReader::SimulationType::MixturePrediction:
      {
        MixturePrediction mixture(reader);

        mixture.print();
        mixture.run();
        mixture.print();
        break;
      }
      case InputReader::SimulationType::Fitting:
      {
        Fitting fitting(reader);

        fitting.run();
        break;
      }
      case InputReader::SimulationType::SwingAdsorption:
      {
        if (reader.isMultibed())
        {
          runSimulation<SwingAdsorption<MultibedColumn>>(reader);
        }
        else
        {
          runSimulation<SwingAdsorption<Column>>(reader);
        }
        break;
      }
    }
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what();
    exit(-1);
  }

  return 0;
}

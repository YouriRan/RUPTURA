#include <exception>

#include "breakthrough.h"
#include "fitting.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "special_functions.h"
#include "swing_adsorption.h"

namespace
{
template <typename Simulation>
void runSimulation(const InputReader& reader)
{
  Simulation simulation(reader);
  simulation.print();
  simulation.run();
}
}  // namespace

int main(void)
{
  try
  {
    InputReader reader("simulation.json");

    switch (reader.simulationType)
    {
      case InputReader::SimulationType::Breakthrough:
      default:
      {
        if (reader.adsorbentComponents.size() > 1)
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
        if (reader.adsorbentComponents.size() > 1)
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

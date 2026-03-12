#include <exception>

#include "breakthrough.h"
#include "fitting.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "special_functions.h"

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
        Breakthrough breakthrough(reader);

        breakthrough.print();
        breakthrough.run();
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
    }
  }
  catch (std::exception const& e)
  {
    std::cerr << e.what();
    exit(-1);
  }

  return 0;
}

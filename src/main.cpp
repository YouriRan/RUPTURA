#include <exception>
#include <cstdio>
#include <filesystem>
#include <print>

#include "breakthrough.h"
#include "fitting.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "notebook_templates.h"
#include "special_functions.h"
#include "swing_adsorption.h"
#if RUPTURA_BUILD_PROCESS
#include "cycle/main_input.hpp"
#include "cycle/plotting.hpp"
#include "cycle/runner.hpp"
#endif

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

int main(int argc, char* argv[])
{
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  try
  {
    if (argc == 2 && std::string_view(argv[1]) == "--capabilities")
    {
      nlohmann::json capabilities={{"sundials",static_cast<bool>(BUILD_SUNDIALS)},
        {"input_validation",true},{"shared_fvm",static_cast<bool>(RUPTURA_BUILD_PROCESS)},
        {"shared_column",static_cast<bool>(RUPTURA_BUILD_PROCESS)},{"version","2.0"}};
#if RUPTURA_BUILD_PROCESS
      capabilities["process"]=nlohmann::json::parse(cycle::process_capabilities().dump());
#endif
      std::println("{}",capabilities.dump());
      return 0;
    }
    bool canonicalize=false,validate=false,restart=false,no_plots=false,plot_only=false;int cycles=0;
    std::filesystem::path input_path="simulation.json";bool path_given=false;
    for(int i=1;i<argc;++i){std::string arg=argv[i];
      if(arg=="--help"||arg=="-h"){std::println("Usage: ruptura [simulation.json] [--validate] [--canonicalize] [--restart] [--cycles N] [--no-plots] [--plot-only]\nWithout arguments, run simulation.json. All times are seconds. Process controls use the shared column with SpatialDiscretization FDM or FVM.");return 0;}
      if(arg=="--canonicalize")canonicalize=true;
      else if(arg=="--validate")validate=true;
      else if(arg=="--restart")restart=true;
      else if(arg=="--no-plots")no_plots=true;
      else if(arg=="--plot-only")plot_only=true;
      else if(arg=="--cycles"){if(++i>=argc)throw std::runtime_error("--cycles needs a positive integer");std::size_t used=0;cycles=std::stoi(argv[i],&used);if(used!=std::string(argv[i]).size()||cycles<1)throw std::runtime_error("--cycles needs a positive integer");}
      else if(arg=="--config"){if(++i>=argc)throw std::runtime_error("--config needs an input path");input_path=argv[i];path_given=true;}
      else if(!arg.empty()&&arg.front()!='-'&&!path_given){input_path=arg;path_given=true;}
      else throw std::runtime_error("Unknown argument: "+arg+". Use ruptura --help.");
    }
    input_path=std::filesystem::absolute(input_path);
#if RUPTURA_BUILD_PROCESS
    if(canonicalize){std::cout<<cycle::canonicalize_chemisorption_input(cycle::read_input_document(input_path)).dump(2)<<'\n';return 0;}
    if(cycle::is_process_input(cycle::read_input_document(input_path))) {
      auto config=cycle::load_main_input(input_path);
      if(restart)config["RESTART"]=true;
      if(cycles){config["RUN_MODE"]="fixed_cycles";config["N_CYCLES"]=cycles;}
      if(no_plots)config["MAKE_PLOTS"]=false;
      if(plot_only){auto metadata=cycle::read_input_document(std::filesystem::path(config["OUTPUT_DIR"].get<std::string>())/"metadata.json");auto output=config["OUTPUT_DIR"];config=cycle::normalize_main_input(metadata.at("input_document"),input_path.parent_path());config["OUTPUT_DIR"]=output;config["_PLOT_ONLY"]=true;}
      cycle::CycleRunner simulation(config);
      if(validate){std::cout<<simulation.validate().dump(2)<<'\n';return 0;}
      if(plot_only){simulation.load_plot_progress();simulation.write_plots();return 0;}
      std::filesystem::current_path(input_path.parent_path());
      ensureDefaultNotebook(InputReader::SimulationType::Breakthrough);
      const auto summary=simulation.run();
      cycle::write_main_outputs(simulation.cfg,input_path.parent_path());
      if(config.value("MAKE_PLOTS",true))simulation.write_plots();
      std::cout<<summary.dump(2)<<'\n';return 0;
    }
#endif
#if RUPTURA_BUILD_PROCESS
    if(plot_only){cycle::write_legacy_breakthrough_plots(input_path);return 0;}
#endif
    if(restart||plot_only||cycles)throw std::runtime_error("These process controls require SpatialDiscretization FDM or FVM. Existing main runs retain ReadColumnFile and NumberOfCycles.");
    if(validate){InputReader reader(input_path.string());std::println("Input is valid.");return 0;}
    std::filesystem::current_path(input_path.parent_path());
    InputReader reader(input_path.string());
    ensureDefaultNotebook(reader.simulationType);

    switch (reader.simulationType)
    {
      case InputReader::SimulationType::Breakthrough:
      default:
      {
        runSimulation<Breakthrough>(reader);
#if RUPTURA_BUILD_PROCESS
        const auto document = cycle::read_input_document(input_path);
        const auto output = document.value("Output", cycle::Json::object());
        if (!no_plots && output.value("Plots", false)) cycle::write_legacy_breakthrough_plots(input_path);
#endif
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
        runSimulation<SwingAdsorption>(reader);
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

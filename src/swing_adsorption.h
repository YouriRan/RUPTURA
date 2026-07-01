#pragma once

#include <cstddef>
#include <string>
#include <tuple>
#include <vector>

#include "breakthrough.h"
#include "column.h"
#include "component.h"
#include "cvode.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "rk3.h"
#include "rk3_si.h"
#include "timing.h"

/**
 * \brief Coordinates repeated breakthrough sub-stages for swing adsorption.
 *
 * Stores the base breakthrough simulation and the temperature/pressure stages
 * used to cycle the column.
 */
struct SwingAdsorption
{
  /**
   * \brief Operating conditions for one swing adsorption sub-stage.
   */
  struct SubStage
  {
    double temperature;    ///< Sub-stage temperature in K.
    double pressure;       ///< Sub-stage pressure in Pa.
    size_t numberOfSteps;  ///< Number of time steps in this sub-stage.
  };

  /**
   * \brief Constructs a swing adsorption simulation from parsed input data.
   */
  SwingAdsorption(const InputReader& inputReader);

  /**
   * \brief Runs all configured swing adsorption sub-stages.
   */
  void run();

  /**
   * \brief Prints the simulation representation to standard output.
   */
  void print() const;

  /**
   * \brief Returns a compact text representation of the swing adsorption setup.
   */
  std::string repr() const;

  Breakthrough breakthrough;        ///< Breakthrough simulation reused across sub-stages.
  std::vector<SubStage> subStages;  ///< Ordered list of swing adsorption sub-stages.
  Timing timings;                   ///< Accumulated timing diagnostics.
};

#include <cstddef>
#include <tuple>
#include <vector>

#include "column.h"
#include "component.h"
#include "cvode.h"
#include "inputreader.h"
#include "mixture_prediction.h"
#include "rk3.h"
#include "rk3_si.h"

#ifdef PYBUILD
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
namespace py = pybind11;
#endif  // PYBUILD

/**
 * \brief Simulates a breakthrough process in an adsorption column.
 *
 * The Breakthrough struct encapsulates the parameters and methods required to simulate
 * the breakthrough of gases in an adsorption column. It handles the initialization of
 * simulation parameters, computation of time steps, and generation of output scripts
 * for plotting and visualization.
 */
struct Breakthrough
{
  enum class IntegrationScheme
  {
    SSP_RK = 0,     ///< Strong Stability Preserving Runge-Kutta method.
    CVODE = 1,      ///< CVODE integration
    Iterative = 2,  ///< Iterative integration scheme.
    SIRK3 = 3,
  };

  /**
   * \brief Constructs a Breakthrough simulation using an InputReader.
   *
   * Initializes a Breakthrough object with parameters specified in the InputReader.
   *
   * \param inputreader Reference to an InputReader containing simulation parameters.
   */
  Breakthrough(const InputReader& inputreader);

  /**
   * \brief Constructs a Breakthrough simulation with specified parameters.
   *
   * Initializes a Breakthrough object with the provided simulation parameters.
   *
   * \param _displayName Name of the simulation for display purposes.
   * \param _components Vector of components involved in the simulation.
   * \param _carrierGasComponent Index of the carrier gas component.
   * \param _numberOfGridPoints Number of grid points in the column.
   * \param _printEvery Frequency of printing time steps to the screen.
   * \param _writeEvery Frequency of writing data to files.
   * \param _temperature Simulation temperature in Kelvin.
   * \param _p_total Total pressure in the column in Pascals.
   * \param _columnVoidFraction Void fraction of the column.
   * \param _pressureGradient Pressure gradient in the column.
   * \param _particleDensity Particle density in kg/m³.
   * \param _columnEntranceVelocity Interstitial velocity at the beginning of the column in m/s.
   * \param _columnLength Length of the column in meters.
   * \param _timeStep Time step for the simulation.
   * \param _numberOfTimeSteps Total number of time steps.
   * \param _autoSteps Flag to use automatic number of steps.
   * \param _pulse Flag to indicate pulsed inlet condition.
   * \param _pulseTime Pulse time.
   * \param _mixture MixturePrediction object for mixture predictions.
   * \param _breakthroughIntegrator Type of integrator used for solving PDE.
   */
  Breakthrough(std::string _displayName, std::vector<Component> _components, size_t _carrierGasComponent,
               size_t _numberOfGridPoints, size_t _printEvery, size_t _writeEvery, double _temperature, double _p_total,
               double _columnVoidFraction, double _pressureGradient, double _particleDensity,
               double _columnEntranceVelocity, double _columnLength, double _timeStep, size_t _numberOfTimeSteps,
               bool _autoSteps, bool _pulse, double _pulseTime, const MixturePrediction _mixture,
               size_t _breakthroughIntegrator);

  /**
   * \brief Prints the representation of the Breakthrough object to the console.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the Breakthrough object.
   *
   * \return A string representing the Breakthrough object.
   */
  std::string repr() const;

  /**
   * \brief Runs the Breakthrough simulation.
   *
   * Executes the simulation over the specified number of time steps.
   */
  void run();

  /**
   * \brief Creates a Gnuplot script for plotting breakthrough curves.
   *
   * Generates a Gnuplot script to visualize the simulation results.
   */
  void createPlotScript();

  /**
   * \brief Creates scripts for generating movies of the simulation.
   *
   * Generates scripts to create movies visualizing the simulation over time.
   */
  void createMovieScripts();

#ifdef PYBUILD
  /**
   * \brief Computes the Breakthrough simulation and returns the results.
   *
   * Executes the simulation and returns a NumPy array containing the simulation data.
   *
   * \return A NumPy array of simulation results.
   */
  py::array_t<double> compute();

  /**
   * \brief Sets the component parameters for the simulation.
   *
   * Updates the mole fractions and isotherm parameters for each component.
   *
   * \param molfracs Vector of mole fractions for each component.
   * \param params Vector of isotherm parameters for the components.
   */
  void setComponentsParameters(std::vector<double> molfracs, std::vector<double> params);

  /**
   * \brief Retrieves the component parameters used in the simulation.
   *
   * Returns the isotherm parameters for each component.
   *
   * \return A vector containing the isotherm parameters for the components.
   */
  std::vector<double> getComponentsParameters();
#endif  // PYBUILD

  const std::string displayName;  ///< Name of the simulation for display purposes.
  size_t carrierGasComponent{0};  ///< Index of the carrier gas component.
  size_t Ncomp;                   ///< Number of components.
  size_t Ngrid;                   ///< Number of grid points.

  size_t printEvery;  ///< Frequency of printing time steps to the screen.
  size_t writeEvery;  ///< Frequency of writing data to files.

  double dt;                ///< Time step for integration.
  size_t Nsteps;            ///< Total number of steps.
  bool autoSteps;           ///< Flag to use automatic number of steps.
  bool pulse;               ///< Pulsed inlet condition for breakthrough.
  double tpulse;            ///< Pulse time.
  size_t maxIsothermTerms;  ///< Maximum number of isotherm terms.

  Column state;
  RungeKutta3 rk3;
  SemiImplicitRungeKutta3 sirk3;
  CVODE cvode;
  IntegrationScheme integrationScheme;

  /**
   * \brief Computes the first derivatives of concentrations and pressures.
   *
   * Calculates the derivatives Dq/dt and Dp/dt along the column.
   *
   * \param dqdt Output vector for the derivatives of Q with respect to time.
   * \param dpdt Output vector for the derivatives of P with respect to time.
   * \param q_eq Equilibrium adsorption amounts.
   * \param q Current adsorption amounts.
   * \param v Interstitial gas velocities.
   * \param p Partial pressures.
   */
  void computeFirstDerivatives(std::vector<double>& dqdt, std::vector<double>& dpdt, const std::vector<double>& q_eq,
                               const std::vector<double>& q, const std::vector<double>& v,
                               const std::vector<double>& p);

  /**
   * \brief Computes a single simulation step.
   *
   * Advances the simulation by one time step.
   *
   * \param step The current time step index.
   */
  void computeStep(size_t step);

  /**
   * \brief Computes the equilibrium loadings for the current time step.
   *
   * Calculates the equilibrium adsorption amounts based on current pressures.
   */
  void computeEquilibriumLoadings();

  /**
   * \brief Computes the interstitial gas velocities along the column.
   *
   * Updates the velocities based on current pressures and adsorption amounts.
   */
  void computeVelocity();

  /**
   * \brief Creates a script to generate a movie for the interstitial gas velocity.
   */
  void createMovieScriptColumnV();

  /**
   * \brief Creates a script to generate a movie for the total pressure along the column.
   */
  void createMovieScriptColumnPt();

  /**
   * \brief Creates a script to generate a movie for the adsorption amounts Q.
   */
  void createMovieScriptColumnQ();

  /**
   * \brief Creates a script to generate a movie for the equilibrium adsorption amounts Qeq.
   */
  void createMovieScriptColumnQeq();

  /**
   * \brief Creates a script to generate a movie for the partial pressures P.
   */
  void createMovieScriptColumnP();

  /**
   * \brief Creates a script to generate a movie for the derivatives of pressure Dp/dt.
   */
  void createMovieScriptColumnDpdt();

  /**
   * \brief Creates a script to generate a movie for the derivatives of adsorption amounts Dq/dt.
   */
  void createMovieScriptColumnDqdt();

  /**
   * \brief Creates a script to generate a movie for the normalized partial pressures.
   */
  void createMovieScriptColumnPnormalized();
};

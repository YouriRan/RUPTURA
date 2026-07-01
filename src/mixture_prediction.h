#pragma once

#include <span>
#include <tuple>
#include <vector>

#include "component.h"
#include "inputreader.h"

/**
 * \brief Class for predicting mixture adsorption isotherms.
 *
 * The MixturePrediction class provides methods to predict mixture adsorption isotherms using various methods like IAST,
 * SIAST, and explicit isotherm models. It handles the computation of adsorbed phase mole fractions and loadings based
 * on gas phase compositions and pressure.
 */
struct MixturePrediction
{
  /**
   * \brief Enum class for prediction methods.
   *
   * Specifies the method used for predicting mixture adsorption isotherms.
   */
  enum class PredictionMethod
  {
    IAST = 0,   ///< Ideal Adsorbed Solution Theory
    SIAST = 1,  ///< Segregated Ideal Adsorbed Solution Theory
    EI = 2,     ///< Explicit Isotherm
    SEI = 3     ///< Segregated Explicit Isotherm
  };

  /**
   * \brief Enum class for IAST methods.
   *
   * Specifies the method used for solving IAST equations.
   */
  enum class IASTMethod
  {
    FastIAST = 0,            ///< Fast IAST algorithm
    NestedLoopBisection = 1  ///< Nested Loop Bisection method
  };

  /**
   * \brief Constructs a MixturePrediction object from an InputReader.
   *
   * Initializes the MixturePrediction instance using the parameters provided by the InputReader.
   *
   * \param inputreader The InputReader containing the simulation parameters.
   */
  MixturePrediction(const InputReader& inputreader);

  /**
   * \brief Constructs a MixturePrediction object with specified parameters.
   *
   * Initializes the MixturePrediction instance using the provided parameters.
   *
   * \param _displayName The display name for the simulation.
   * \param _components A vector of Component objects representing the mixture components.
   * \param _numberOfCarrierGases The number of carrier gases in the mixture.
   * \param _carrierGasComponent The index of the carrier gas component.
   * \param _temperature The temperature of the system.
   * \param _pressureStart The starting pressure for the simulation.
   * \param _pressureEnd The ending pressure for the simulation.
   * \param _numberOfPressurePoints The number of pressure points in the simulation.
   * \param _pressureScale The pressure scale (0 for Log, 1 for Normal).
   * \param _predictionMethod The prediction method to use.
   * \param _iastMethod The IAST method to use.
   */
  MixturePrediction(std::string _displayName, std::vector<Component> _components, size_t _numberOfCarrierGases,
                    size_t _carrierGasComponent, double _temperature, double _pressureStart, double _pressureEnd,
                    size_t _numberOfPressurePoints, size_t _pressureScale, size_t _predictionMethod,
                    size_t _iastMethod);

  /**
   * \brief Prints the mixture prediction data.
   *
   * Outputs the mixture prediction data to the standard output.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the MixturePrediction object.
   *
   * \return A string representing the MixturePrediction object.
   */
  std::string repr() const;

  /**
   * \brief Runs the mixture prediction simulation.
   *
   * Performs the mixture prediction calculations and writes the results to output files.
   */
  void run();

  /**
   * \brief Gets the maximum number of isotherm terms.
   *
   * \return The maximum number of isotherm terms.
   */
  size_t getMaxIsothermTerms() const { return maxIsothermTerms; }

  /**
   * \brief Predicts the mixture adsorption isotherm.
   *
   * Computes the adsorbed phase mole fractions and loadings based on the gas phase compositions and pressure.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \param gasTemperature Temperature of the gas.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> predictMixture(std::span<const double> idealGasMolFractions, const double& externalPressure,
                                           std::span<double> adsorbedMolFractions, std::span<double> numberOfMolecules,
                                           double* cachedPressure, double* cachedGrandPotential,
                                           double& gasTemperature);

  std::string displayName;                  ///< The display name for the simulation.
  std::vector<Component> components;        ///< The vector of components in the mixture.
  std::vector<Component> sortedComponents;  ///< Components sorted according to specific criteria.
  size_t numberOfComponents;                ///< The total number of components.
  size_t numberOfSortedComponents;          ///< The number of sorted components.
  size_t numberOfCarrierGases;              ///< The number of carrier gases in the mixture.
  size_t carrierGasComponent;               ///< The index of the carrier gas component.
  PredictionMethod predictionMethod;        ///< The method used for predicting mixture adsorption isotherms.
  IASTMethod iastMethod;                    ///< The method used for solving IAST equations.
  size_t maxIsothermTerms;                  ///< The maximum number of isotherm terms.
  std::vector<std::vector<Component>>
      segregatedSortedComponents;  ///< Segregated and sorted components for SIAST/SEI methods.

  std::vector<double> firstExplicitIsothermAlpha;    ///< Intermediate calculation vector for explicit isotherms.
  std::vector<double> secondExplicitIsothermAlpha;   ///< Intermediate calculation vector for explicit isotherms.
  std::vector<double> explicitIsothermAlphaProduct;  ///< Product of alpha values for explicit isotherms.
  std::vector<double> adsorbedMoleFractionsScratch;  ///< Intermediate calculation vector.

  std::vector<double> hypotheticalPressure;   ///< Hypothetical pressures in IAST calculations in Pa.
  std::vector<double> reducedGrandPotential;  ///< Reduced grand potential values.
  std::vector<double> residualVector;         ///< Intermediate calculation vector in IAST.
  std::vector<double> correctionVector;       ///< Correction vector in IAST.
  std::vector<double> jacobianMatrix;         ///< Jacobian matrix in IAST calculations.

  /**
   * \brief Enum class for pressure scales.
   *
   * Specifies the scale to use for pressure in the simulation.
   */
  enum class PressureScale
  {
    Log = 0,    ///< Logarithmic pressure scale
    Normal = 1  ///< Linear pressure scale
  };
  double temperature{300.0};                        ///< The temperature of the system in K.
  double pressureStart{1e3};                        ///< The starting pressure for the simulation in Pa.
  double pressureEnd{1e8};                          ///< The ending pressure for the simulation in Pa.
  size_t numberOfPressurePoints{100};               ///< The number of pressure points in the simulation.
  PressureScale pressureScale{PressureScale::Log};  ///< The pressure scale to use.

  /**
   * \brief Initializes the pressure points for the simulation.
   *
   * Generates a vector of pressure points based on the starting and ending pressures and the pressure scale.
   *
   * \return A vector containing the pressure points for the simulation.
   */
  std::vector<double> initPressures();

  /**
   * \brief Sorts the components based on specific criteria.
   *
   * Sorts the components to optimize calculations in prediction methods.
   */
  void sortComponents();

  /**
   * \brief Computes mixture prediction using Fast IAST method.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \return A pair containing the number of IAST steps and a status code.a
   */
  std::pair<size_t, size_t> computeFastIAST(std::span<const double> idealGasMolFractions,
                                            const double& externalPressure, std::span<double> adsorbedMolFractions,
                                            std::span<double> numberOfMolecules, double* cachedPressure,
                                            double* cachedGrandPotential, double& gasTemperature);

  /**
   * \brief Computes mixture prediction using Fast SIAST method.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeFastSIAST(std::span<const double> idealGasMolFractions,
                                             const double& externalPressure, std::span<double> adsorbedMolFractions,
                                             std::span<double> numberOfMolecules, double* cachedPressure,
                                             double* cachedGrandPotential, double& gasTemperature);

  /**
   * \brief Computes mixture prediction for a specific term using Fast SIAST method.
   *
   * \param term The index of the isotherm term.
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeFastSIAST(size_t term, std::span<const double> idealGasMolFractions,
                                             const double& externalPressure, std::span<double> adsorbedMolFractions,
                                             std::span<double> numberOfMolecules, double* cachedPressure,
                                             double* cachedGrandPotential, double& gasTemperature);

  /**
   * \brief Computes mixture prediction using IAST with nested loop bisection method.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeIASTNestedLoopBisection(std::span<const double> idealGasMolFractions,
                                                           const double& externalPressure,
                                                           std::span<double> adsorbedMolFractions,
                                                           std::span<double> numberOfMolecules, double* cachedPressure,
                                                           double* cachedGrandPotential, double& gasTemperature);

  /**
   * \brief Computes mixture prediction using SIAST with nested loop bisection method.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeSIASTNestedLoopBisection(std::span<const double> idealGasMolFractions,
                                                            const double& externalPressure,
                                                            std::span<double> adsorbedMolFractions,
                                                            std::span<double> numberOfMolecules, double* cachedPressure,
                                                            double* cachedGrandPotential, double& gasTemperature);

  /**
   * \brief Computes mixture prediction for a specific term using SIAST with nested loop bisection method.
   *
   * \param term The index of the isotherm term.
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \param cachedPressure An array to cache intermediate pressure calculations.
   * \param cachedGrandPotential An array to cache intermediate reducedGrandPotential calculations.
   * \return A pair containing the number of IAST steps and a status code.
   */
  std::pair<size_t, size_t> computeSIASTNestedLoopBisection(size_t term, std::span<const double> idealGasMolFractions,
                                                            const double& externalPressure,
                                                            std::span<double> adsorbedMolFractions,
                                                            std::span<double> numberOfMolecules, double* cachedPressure,
                                                            double* cachedGrandPotential, double& gasTemperature);

  /**
   * \brief Computes mixture prediction using explicit isotherm model.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \return A pair containing the number of steps and a status code.
   */
  std::pair<size_t, size_t> computeExplicitIsotherm(std::span<const double> idealGasMolFractions,
                                                    const double& externalPressure,
                                                    std::span<double> adsorbedMolFractions,
                                                    std::span<double> numberOfMolecules, double& gasTemperature);

  /**
   * \brief Computes mixture prediction using segregated explicit isotherm model.
   *
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \return A pair containing the number of steps and a status code.
   */
  std::pair<size_t, size_t> computeSegratedExplicitIsotherm(std::span<const double> idealGasMolFractions,
                                                            const double& externalPressure,
                                                            std::span<double> adsorbedMolFractions,
                                                            std::span<double> numberOfMolecules,
                                                            double& gasTemperature);

  /**
   * \brief Computes mixture prediction for a specific term using segregated explicit isotherm model.
   *
   * \param site The index of the isotherm site.
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param P The total pressure.
   * \param adsorbedMolFractions The adsorbed phase mole fractions (output).
   * \param numberOfMolecules The number of adsorbed molecules of each component (output).
   * \return A pair containing the number of steps and a status code.
   */
  std::pair<size_t, size_t> computeSegratedExplicitIsotherm(size_t site, std::span<const double> idealGasMolFractions,
                                                            const double& externalPressure,
                                                            std::span<double> adsorbedMolFractions,
                                                            std::span<double> numberOfMolecules,
                                                            double& gasTemperature);

  /**
   * \brief Prints error status for debugging purposes.
   *
   * Outputs the current state of variables when an error occurs in IAST calculations.
   *
   * \param reducedGrandPotential The current reducedGrandPotential value.
   * \param sum The current sum of mole fractions.
   * \param P The total pressure.
   * \param idealGasMolFractions The gas phase mole fractions.
   * \param cachedPressure An array of cached pressure values.
   */
  void printErrorStatus(double reducedGrandPotential, double sum, double P,
                        std::span<const double> idealGasMolFractions, double cachedPressure[], double gasTemperature);
};

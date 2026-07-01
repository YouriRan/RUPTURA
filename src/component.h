#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "multi_site_isotherm.h"
#include "utils.h"

/**
 * \brief Represents a chemical component in the simulation.
 *
 * The Component struct encapsulates the properties and behaviors of a chemical component within the system.
 * It includes identifiers, names, isotherm data, and parameters related to mass transfer and diffusion.
 */
struct Component
{
  /**
   * \brief Constructs a Component with an id and name.
   *
   * Initializes a Component using the provided identifier and name.
   */
  Component(size_t id, std::string name);

  /**
   * \brief Constructs a Component with a vector of single-site isotherms.
   *
   * Initializes the component and constructs its MultiSiteIsotherm from the supplied isotherm sites.
   */
  Component(size_t id, std::string name, std::vector<Isotherm> isotherms, double initialGasMoleFraction,
            double massTransferCoefficient, double axialDispersionCoefficient, bool isCarrierGas = false,
            double molecularWeight = 1.0, double heatOfAdsorption = 0.0);

  /**
   * \brief Constructs a Component with an already-constructed multi-site isotherm.
   *
   * This overload is convenient from Python after constructing a MultiSiteIsotherm directly.
   */
  Component(size_t id, std::string name, MultiSiteIsotherm isotherm, double initialGasMoleFraction,
            double massTransferCoefficient, double axialDispersionCoefficient, bool isCarrierGas = false,
            double molecularWeight = 1.0, double heatOfAdsorption = 0.0);

  size_t id{0};                                              ///< Identifier of the component.
  std::string name{};                                        ///< Name of the component.
  std::string filename{};                                    ///< Filename associated with the component data.
  MultiSiteIsotherm isotherm;                                ///< Isotherm information for the component.
  double initialGasMoleFraction{0.0};                        ///< Gas-phase mole fraction [-].
  double massTransferCoefficient{0.0};                       ///< Mass transfer coefficient in 1/s.
  double axialDispersionCoefficient{0.0};                    ///< Axial dispersion coefficient in m^2/s.
  double heatOfAdsorption{0.0};                              ///< Heat of adsorption in J/mol.
  bool isCarrierGas{false};                                  ///< Flag indicating if this is the carrier gas.
  double molecularWeight{1.0};                               ///< Molecular weight in kg/mol.
  bool nonIsothermal{false};                                 ///< Enables temperature scaling when true.
  std::optional<double> referenceTemperature{std::nullopt};  ///< Reference temperature in K.

  /**
   * \brief Prints the component information to the console.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the Component.
   */
  std::string repr() const;

  /**
   * \brief Returns the isotherm scaling factor at a gas temperature.
   */
  inline double scale(double temperature) const
  {
    if (!nonIsothermal) return 1.0;

    const double beta = 1.0 / (R * temperature);
    const double betaReference = referenceTemperature.has_value() ? 1.0 / (R * referenceTemperature.value()) : 0.0;
    return std::exp((beta - betaReference) * heatOfAdsorption);
  }
};

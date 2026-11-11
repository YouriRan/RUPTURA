#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>
#define _USE_MATH_DEFINES
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <limits>
#include <print>
#include <span>
#include <string>
#include <utility>

#include "random_numbers.h"
#include "special_functions.h"
#include "temperature_law.h"
#include "utils.h"
#include "water_saturation_pressure.h"

// Langmuir:
// parameter 0: K
// parameter 1: N
//
// Langmuir-Freundlich
// parameter 0: K
// parameter 1: N
// parameter 2: power
//
// GAB
// parameter 0: q_m,mono [mol/kg]
// parameter 1: C [-] (temperature-dependent through a VantHoff law)
// parameter 2: K [-] (temperature-dependent through a VantHoff law)
//
// Humidity-dependent sites (RH is the relative humidity, q_w the equilibrium water loading in mol/kg):
// Humidity-enhanced Langmuir   [q_sat, b, alpha, k_w]   b_eff = b [1 + alpha k_w RH / (1 + k_w RH)]
// Moisture-swing Langmuir      [q_sat, b, n_w]          b_eff = b RH^(-n_w)
// GAB-coadsorption Langmuir    [q_sat, b, alpha, omega] b_eff = b (1 + alpha q_w), q_eff = q_sat max(0, 1 - omega q_w)
// GAB-coadsorption Sips, Toth  [q_sat, b, nu, alpha, omega], the same modifiers
// GAB-modified Toth            [q_sat, b, nu, gamma, beta]  q_eff = q_sat / (1 - gamma q_w), b_eff = b (1 + beta q_w)
// Exponential-humidity Toth    [q_sat, b, nu, a]        q_eff = q_sat exp(a RH)

/**
 * \brief Inverse pressure 1/p at which psiForPressure(p) equals the reduced grand potential, by bracketing and bisection.
 *
 * Used for single isotherms without an analytical inverse and for multi-site isotherms. The bracket starts at cachedP0
 * (5 when unset), widens by factors of 2, and is bisected to a relative width of 1e-15; the final pressure is stored in
 * cachedP0 for the next call.
 *
 * \param psiForPressure Reduced grand potential as a function of pressure; increasing in p.
 * \param reduced_grand_potential The reduced grand potential to invert.
 * \param cachedP0 Starting pressure on input (when positive), the solution on output.
 * \return The inverse of the pressure corresponding to the given reduced grand potential.
 */
template <typename PsiForPressure>
double inversePressureByBisection(const PsiForPressure& psiForPressure, double reduced_grand_potential,
                                  double& cachedP0)
{
  const double tiny = 1.0e-15;

  // from here on, work with pressure, and return 1.0 / pressure at the end of the routine
  double p_start;
  if (cachedP0 <= 0.0)
  {
    p_start = 5.0;
  }
  else
  {
    // use the last value of Pi0
    p_start = cachedP0;
  }

  // use bisection algorithm
  double s = psiForPressure(p_start);

  size_t nr_steps = 0;
  double left_bracket = p_start;
  double right_bracket = p_start;

  if (s < reduced_grand_potential)
  {
    // find the bracket on the right
    do
    {
      right_bracket *= 2.0;
      s = psiForPressure(right_bracket);

      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("reduced_grand_potential: {}\n", reduced_grand_potential);
        std::print("reducedGrandPotential: {}\n", s);
        std::print("p_start: {}\n", p_start);
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
      }
    } while (s < reduced_grand_potential);
  }
  else
  {
    // find the bracket on the left
    do
    {
      left_bracket *= 0.5;
      s = psiForPressure(left_bracket);

      ++nr_steps;
      if (nr_steps > 100000)
      {
        std::print("reduced_grand_potential: {}\n", reduced_grand_potential);
        std::print("reducedGrandPotential: {}\n", s);
        std::print("p_start: {}\n", p_start);
        std::print("Left bracket: {}\n", left_bracket);
        std::print("Right bracket: {}\n", right_bracket);
        throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum > 1) does NOT converge\n");
      }
    } while (s > reduced_grand_potential);
  }

  do
  {
    double middle = 0.5 * (left_bracket + right_bracket);
    s = psiForPressure(middle);

    if (s > reduced_grand_potential)
      right_bracket = middle;
    else
      left_bracket = middle;

    ++nr_steps;
    if (nr_steps > 100000)
    {
      std::print("Left bracket: {}\n", left_bracket);
      std::print("Right bracket: {}\n", right_bracket);
      throw std::runtime_error("Error (Inverse bisection): initial bracketing (for sum < 1) does NOT converge\n");
    }
  } while (std::abs(left_bracket - right_bracket) / std::abs(left_bracket + right_bracket) > tiny);

  double middle = 0.5 * (left_bracket + right_bracket);

  //  Store the last value of Pi0
  cachedP0 = middle;

  return 1.0 / middle;
}

/**
 * \brief Represents an isotherm model for adsorption processes.
 *
 * The Isotherm struct encapsulates various isotherm models used to describe the adsorption equilibrium
 * between a fluid and a solid at a constant temperature. It supports different types of isotherm models,
 * such as Langmuir, Freundlich, BET, and others.
 */
struct Isotherm
{
  /**
   * \brief Enumeration of the different types of isotherm models.
   */
  enum class Type
  {
    Langmuir = 0,             ///< Langmuir isotherm model
    Anti_Langmuir = 1,        ///< Anti-Langmuir isotherm model
    BET = 2,                  ///< Brunauer–Emmett–Teller (BET) isotherm model
    Henry = 3,                ///< Henry's law isotherm model
    Freundlich = 4,           ///< Freundlich isotherm model
    Sips = 5,                 ///< Sips isotherm model
    Langmuir_Freundlich = 6,  ///< Langmuir-Freundlich isotherm model
    Redlich_Peterson = 7,     ///< Redlich-Peterson isotherm model
    Toth = 8,                 ///< Toth isotherm model
    Unilan = 9,               ///< Unilan isotherm model
    OBrien_Myers = 10,        ///< O'Brien and Myers isotherm model
    Quadratic = 11,           ///< Quadratic isotherm model
    Temkin = 12,              ///< Temkin isotherm model
    BingelWalton = 13,        ///< Bingel and Walton isotherm model
    Langmuir_pH = 14,         ///< pH-dependent Langmuir isotherm model
    GAB = 15,                 ///< Guggenheim-Anderson-de Boer water isotherm model
    Langmuir_humidity = 16,         ///< Humidity-enhanced Langmuir isotherm model
    Langmuir_moisture_swing = 17,   ///< Moisture-swing (Wang-Lackner) Langmuir isotherm model
    GAB_Langmuir = 18,              ///< Langmuir isotherm modified by the GAB water loading
    GAB_Sips = 19,                  ///< Sips isotherm modified by the GAB water loading
    GAB_Toth = 20,                  ///< Toth isotherm modified by the GAB water loading
    GAB_modified_Toth = 21,         ///< Toth isotherm of Stampi-Bombelli et al. with water-dependent q_sat and b
    Toth_exponential_humidity = 22  ///< Toth isotherm with the capacity factor exp(a RH)
  };

  /**
   * \brief Returns whether an isotherm type depends on the relative humidity (all but GAB itself, which is the water
   *        isotherm).
   */
  static constexpr bool isHumidityDependent(Type type) noexcept
  {
    return type >= Type::Langmuir_humidity && type <= Type::Toth_exponential_humidity;
  }

  /**
   * \brief Returns whether an isotherm type depends on the water loading, which the water component supplies.
   */
  static constexpr bool usesWaterLoading(Type type) noexcept
  {
    return type >= Type::GAB_Langmuir && type <= Type::GAB_modified_Toth;
  }

  /**
   * \brief Constructs an Isotherm with specified type and parameters.
   *
   * Initializes an Isotherm object with the given isotherm type, parameter values, and the number of parameters.
   *
   * \param type The type of the isotherm model.
   * \param values A vector of parameter values for the isotherm model.
   * \param temperatureLaws Temperature-dependent parameters as (parameter index, law) pairs; the value of each law
   *        replaces the corresponding entry of values.
   */
  Isotherm(Isotherm::Type type, const std::vector<double>& values,
           std::vector<std::pair<size_t, TemperatureLaw>> temperatureLaws = {});

  /**
   * \brief Constructs an Isotherm with specified type index and parameters.
   *
   * Initializes an Isotherm object using the type index, parameter values, and the number of parameters.
   *
   * \param t The index of the isotherm type.
   * \param values A vector of parameter values for the isotherm model.
   */
  Isotherm(size_t t, const std::vector<double>& parameters);

  Isotherm::Type type;             ///< The type of the isotherm model.
  std::vector<double> parameters;  ///< Parameter values at the last setTemperature; units depend on isotherm type.
  size_t numberOfParameters;       ///< Number of active parameters in the isotherm.
  /// Temperature-dependent parameters as (index in parameters, law) pairs; the other parameters are constant.
  std::vector<std::pair<size_t, TemperatureLaw>> temperatureLaws{};
  /// Temperature in K of the last setTemperature. GAB also evaluates its saturation pressure at it.
  double temperature{std::numeric_limits<double>::quiet_NaN()};
  /// Saturation pressure of water, from which GAB forms its relative pressure p / p_sat(T).
  WaterSaturationPressure waterSaturationPressure{};

  /**
   * \brief Evaluates the temperature-dependent parameters at a temperature in K.
   *
   * Each law writes its value into parameters, so that value(), psiForPressure(), and inversePressureForPsi() describe
   * the isotherm at that temperature. Constant parameters are left unchanged.
   */
  inline void setTemperature(double newTemperature)
  {
    temperature = newTemperature;
    for (const auto& [index, law] : temperatureLaws)
    {
      parameters[index] = law(newTemperature);
    }
  }

  /**
   * \brief Returns whether the isotherm depends on temperature: through a temperature law, or GAB through p_sat.
   */
  [[nodiscard]] bool dependsOnTemperature() const noexcept { return !temperatureLaws.empty() || type == Type::GAB; }

  /**
   * \brief Returns the temperature law of a parameter, or nullptr if the parameter is constant.
   */
  [[nodiscard]] const TemperatureLaw* temperatureLaw(size_t index) const noexcept;

  /**
   * \brief Returns, for each entry of the Parameters array of an isotherm type, whether and how it may depend on
   *        temperature.
   *
   * The span is empty for types without temperature-dependent parameters. Positions with an empty law mask are
   * constant. The input reader and the graphical interfaces use this table.
   */
  static std::span<const TemperatureParameterSpec> temperatureParameters(Type type) noexcept;

  /**
   * \brief Returns whether this site can contribute adsorption loading.
   *
   * The first parameter is the loading coefficient (or saturation loading) for every supported isotherm type. A zero
   * value therefore makes the site inactive while leaving it available for multisite index alignment.
   */
  [[nodiscard]] bool enabled() const noexcept { return !parameters.empty() && parameters[0] != 0.0; }

  /**
   * \brief Prints a representation of the isotherm to standard output.
   */
  void print() const;

  /**
   * \brief Returns a string representation of the isotherm.
   *
   * \return A string representing the isotherm model and its parameters.
   */
  std::string repr() const;

  /**
   * \brief Computes the adsorption amount at a given pressure.
   *
   * Calculates the adsorption loading based on the isotherm model and parameters for the specified pressure.
   *
   * \param pressure The pressure at which to evaluate the isotherm.
   * \param pH The local pH (used by pH-Langmuir only).
   * \param relativeHumidity The local relative humidity [-] (used by the humidity-dependent types only).
   * \param waterLoading The equilibrium water loading in mol/kg (used by the GAB-coadsorption types only).
   * \return The adsorption amount at the given pressure and the last setTemperature.
   */
  inline double value(double pressure, double pH = 7.0, double relativeHumidity = 0.0, double waterLoading = 0.0) const
  {
    if (!enabled()) return 0.0;

    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double bp = parameters[1] * pressure;
        return parameters[0] * bp / (1.0 + bp);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        return parameters[0] * pressure / (1.0 - parameters[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return parameters[0] * parameters[1] * pressure /
               ((1.0 - parameters[2] * pressure) * (1.0 - parameters[2] + parameters[1] * pressure));
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return parameters[0] * std::pow(pressure, 1.0 / parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        double temp = std::pow(parameters[1] * pressure, 1.0 / parameters[2]);
        return parameters[0] * temp / (1.0 + temp);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double bp = parameters[1] * std::pow(pressure, parameters[2]);
        return parameters[0] * bp / (1.0 + bp);
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        return parameters[0] * pressure / (1.0 + parameters[1] * std::pow(pressure, parameters[2]));
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        return parameters[0] * temp / std::pow(1.0 + std::pow(temp, parameters[2]), 1.0 / parameters[2]);
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] *
               (temp1 / temp2 + 0.5 * parameters[2] * parameters[2] * temp1 * (1.0 - temp1) / (temp2 * temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = parameters[2] * pressure * pressure;
        return parameters[0] * (temp1 + 2.0 * temp2) / (1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (temp1 + parameters[2] * temp1 * temp1 * (temp1 - 1.0));
      }
      case Isotherm::Type::BingelWalton:
      {
        return parameters[0] * (1.0 - std::exp(-(parameters[1] + parameters[2]) * pressure)) /
               (1.0 + (parameters[2] / parameters[1]) * std::exp(-(parameters[1] + parameters[2]) * pressure));
      }
      default:
        return valueSlowPath(pressure, pH, relativeHumidity, waterLoading);
    }
  }

  /**
   * \brief Computes the slope dq/dp of the adsorption amount at a given pressure.
   *
   * The analytical derivative of value() with respect to pressure, at the last setTemperature. Freundlich, Sips, and
   * Langmuir-Freundlich with a power of p below one (parameters[1] > 1, parameters[2] > 1, parameters[2] < 1) have an
   * unbounded slope as p -> 0; callers pass a floored pressure.
   *
   * \param pressure The pressure at which to evaluate the slope.
   * \param pH The local pH (used by pH-Langmuir only; its pH dependence is not differentiated).
   * \param relativeHumidity The local relative humidity [-], held fixed (humidity-dependent types only).
   * \param waterLoading The equilibrium water loading in mol/kg, held fixed (GAB-coadsorption types only).
   * \return The slope of the adsorption amount at the given pressure, in loading units per pressure unit.
   */
  inline double slope(double pressure, double pH = 7.0, double relativeHumidity = 0.0, double waterLoading = 0.0) const
  {
    if (!enabled()) return 0.0;

    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double denominator = 1.0 + parameters[1] * pressure;
        return parameters[0] * parameters[1] / (denominator * denominator);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        double denominator = 1.0 - parameters[1] * pressure;
        return parameters[0] / (denominator * denominator);
      }
      case Isotherm::Type::BET:
      {
        double denominator = (1.0 - parameters[2] * pressure) * (1.0 - parameters[2] + parameters[1] * pressure);
        double denominatorSlope = -parameters[2] * (1.0 - parameters[2] + parameters[1] * pressure) +
                                  parameters[1] * (1.0 - parameters[2] * pressure);
        return parameters[0] * parameters[1] * (denominator - pressure * denominatorSlope) /
               (denominator * denominator);
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0];
      }
      case Isotherm::Type::Freundlich:
      {
        return (parameters[0] / parameters[1]) * std::pow(pressure, 1.0 / parameters[1] - 1.0);
      }
      case Isotherm::Type::Sips:
      {
        double temp = std::pow(parameters[1] * pressure, 1.0 / parameters[2]);
        double tempSlope = std::pow(parameters[1], 1.0 / parameters[2]) *
                           std::pow(pressure, 1.0 / parameters[2] - 1.0) / parameters[2];
        return parameters[0] * tempSlope / ((1.0 + temp) * (1.0 + temp));
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double bp = parameters[1] * std::pow(pressure, parameters[2]);
        double bpSlope = parameters[1] * parameters[2] * std::pow(pressure, parameters[2] - 1.0);
        return parameters[0] * bpSlope / ((1.0 + bp) * (1.0 + bp));
      }
      case Isotherm::Type::Redlich_Peterson:
      {
        double bp = parameters[1] * std::pow(pressure, parameters[2]);
        return parameters[0] * (1.0 + (1.0 - parameters[2]) * bp) / ((1.0 + bp) * (1.0 + bp));
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        return parameters[0] * parameters[1] *
               std::pow(1.0 + std::pow(temp, parameters[2]), -(1.0 + parameters[2]) / parameters[2]);
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] * parameters[1] *
               (1.0 / (temp2 * temp2) + 0.5 * parameters[2] * parameters[2] * (1.0 - 4.0 * temp1 + temp1 * temp1) /
                                            (temp2 * temp2 * temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double numerator = parameters[1] * pressure + 2.0 * parameters[2] * pressure * pressure;
        double denominator = 1.0 + parameters[1] * pressure + parameters[2] * pressure * pressure;
        double numeratorSlope = parameters[1] + 4.0 * parameters[2] * pressure;
        double denominatorSlope = parameters[1] + 2.0 * parameters[2] * pressure;
        return parameters[0] * (numeratorSlope * denominator - numerator * denominatorSlope) /
               (denominator * denominator);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (1.0 + parameters[2] * (3.0 * temp1 * temp1 - 2.0 * temp1)) * parameters[1] /
               ((1.0 + temp) * (1.0 + temp));
      }
      case Isotherm::Type::BingelWalton:
      {
        double rate = parameters[1] + parameters[2];
        double ratio = parameters[2] / parameters[1];
        double decay = std::exp(-rate * pressure);
        return parameters[0] * rate * decay * (1.0 + ratio) / ((1.0 + ratio * decay) * (1.0 + ratio * decay));
      }
      default:
        return slopeSlowPath(pressure, pH, relativeHumidity, waterLoading);
    }
  }

  /**
   * \brief Computes the reduced grand potential (spreading pressure) at a given pressure.
   *
   * Calculates the reduced grand potential reducedGrandPotential based on the isotherm model and parameters for the
   * specified pressure.
   *
   * \param pressure The pressure at which to evaluate the reduced grand potential.
   * \return The reduced grand potential reducedGrandPotential at the given pressure.
   */
  inline double psiForPressure(double pressure) const
  {
    if (!enabled()) return 0.0;

    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double bp = parameters[1] * pressure;
        return parameters[0] * std::log(1.0 + bp);
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        return -(parameters[0] / parameters[1]) * std::log(1.0 - parameters[1] * pressure);
      }
      case Isotherm::Type::BET:
      {
        return (parameters[0] * parameters[1]) *
               std::log((1.0 - parameters[2] + parameters[1] * pressure) /
                        ((1.0 - parameters[2]) * (1.0 - parameters[2] * pressure))) /
               (parameters[1] + parameters[2] - parameters[2] * parameters[2]);
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] * pressure;
      }
      case Isotherm::Type::Freundlich:
      {
        return parameters[0] * parameters[1] * std::pow(pressure, 1.0 / parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[2] * parameters[0] *
               std::log(1.0 + std::pow(parameters[1] * pressure, 1.0 / parameters[2]));
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double bp = parameters[1] * std::pow(pressure, parameters[2]);
        return (parameters[0] / parameters[2]) * std::log(1.0 + bp);
      }
      case Isotherm::Type::OBrien_Myers:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = 1.0 + temp1;
        return parameters[0] * (std::log(temp2) + 0.5 * parameters[2] * parameters[2] * temp1 / (temp2 * temp2));
      }
      case Isotherm::Type::Quadratic:
      {
        double temp1 = parameters[1] * pressure;
        double temp2 = parameters[2] * pressure * pressure;
        return parameters[0] * std::log(1.0 + temp1 + temp2);
      }
      case Isotherm::Type::Temkin:
      {
        double temp = parameters[1] * pressure;
        double temp1 = temp / (1.0 + temp);
        return parameters[0] * (std::log(1.0 + temp) - 0.5 * parameters[2] * temp1 * temp1);
      }
      default:
        return psiForPressureSlowPath(pressure);
    }
  }

  /**
   * \brief Computes the inverse pressure corresponding to a given reduced grand potential reducedGrandPotential.
   *
   * Calculates the pressure that corresponds to the specified reduced grand potential reducedGrandPotential using the
   * isotherm model. This function may cache intermediate results to improve performance in repeated calculations.
   *
   * \param reduced_grand_potential The reduced grand potential reducedGrandPotential.
   * \param cachedP0 A reference to a cached pressure value used to initialize the calculation.
   * \return The inverse of the pressure corresponding to the given reducedGrandPotential.
   */
  inline double inversePressureForPsi(double reduced_grand_potential, double& cachedP0) const
  {
    if (!enabled()) return 0.0;

    switch (type)
    {
      case Isotherm::Type::Langmuir:
      {
        double denominator = std::exp(reduced_grand_potential / parameters[0]) - 1.0;
        double b = parameters[1];
        return b / denominator;
      }
      case Isotherm::Type::Anti_Langmuir:
      {
        double denominator = 1.0 - std::exp(-parameters[1] * reduced_grand_potential / parameters[0]);
        return parameters[1] / denominator;
      }
      case Isotherm::Type::Henry:
      {
        return parameters[0] / reduced_grand_potential;
      }
      case Isotherm::Type::Freundlich:
      {
        return std::pow((parameters[0] * parameters[1]) / reduced_grand_potential, parameters[1]);
      }
      case Isotherm::Type::Sips:
      {
        return parameters[1] /
               std::pow((std::exp(reduced_grand_potential / (parameters[2] * parameters[0])) - 1.0), parameters[2]);
      }
      case Isotherm::Type::Langmuir_Freundlich:
      {
        double denominator = std::exp(reduced_grand_potential * parameters[2] / parameters[0]) - 1.0;
        double b = parameters[1];
        return std::pow(b / denominator, 1.0 / parameters[2]);
      }
      case Isotherm::Type::Langmuir_pH:
        throw std::runtime_error("Error: pH-dependent Langmuir is not compatible with spreading-pressure methods");
      default:
        return inversePressureByBisection([this](double pressure) { return psiForPressure(pressure); },
                                          reduced_grand_potential, cachedP0);
    }
  }

  /**
   * \brief Out-of-line part of value() for the isotherm types that need many registers.
   *
   * Keeping them out of the inline value() spares the common types from saving all callee-saved registers on
   * every call.
   */
  [[gnu::noinline]] double valueSlowPath(double pressure, double pH, double relativeHumidity,
                                         double waterLoading) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir_pH:
      {
        const double pHFactor = std::pow(10.0, std::clamp(parameters[2] - pH, -300.0, 300.0));
        const double bp = parameters[1] * pressure / (1.0 + pHFactor);
        return parameters[0] * bp / (1.0 + bp);
      }
      case Isotherm::Type::Unilan:
      {
        double temp1 = 1.0 + parameters[1] * std::exp(parameters[2]) * pressure;
        double temp2 = 1.0 + parameters[1] * std::exp(-parameters[2]) * pressure;
        return parameters[0] * (0.5 / parameters[2]) * std::log(temp1 / temp2);
      }
      case Isotherm::Type::Langmuir_humidity:
      case Isotherm::Type::Langmuir_moisture_swing:
      case Isotherm::Type::GAB_Langmuir:
      {
        [[maybe_unused]] const auto [saturationLoading, affinity, exponent] =
            humidityParameters(relativeHumidity, waterLoading);
        // An infinite affinity (the moisture-swing Langmuir at RH = 0) saturates the site at any positive pressure.
        if (std::isinf(affinity)) return pressure > 0.0 ? saturationLoading : 0.0;
        const double bp = affinity * pressure;
        return saturationLoading * bp / (1.0 + bp);
      }
      case Isotherm::Type::GAB_Sips:
      {
        const auto [saturationLoading, affinity, exponent] = humidityParameters(relativeHumidity, waterLoading);
        const double temp = std::pow(affinity * pressure, 1.0 / exponent);
        return saturationLoading * temp / (1.0 + temp);
      }
      case Isotherm::Type::GAB_Toth:
      case Isotherm::Type::GAB_modified_Toth:
      case Isotherm::Type::Toth_exponential_humidity:
      {
        const auto [saturationLoading, affinity, exponent] = humidityParameters(relativeHumidity, waterLoading);
        const double temp = affinity * pressure;
        return saturationLoading * temp / std::pow(1.0 + std::pow(temp, exponent), 1.0 / exponent);
      }
      case Isotherm::Type::GAB:
      {
        if (pressure <= 0.0) return 0.0;
        if (!(temperature > 0.0))
        {
          throw std::domain_error("Error: GAB needs a positive temperature (set it with setTemperature)");
        }

        const double saturationPressure = waterSaturationPressure(temperature);
        const double relativePressure = pressure / saturationPressure;
        // C and K are already evaluated at the temperature by their VantHoff laws (setTemperature).
        const double c = parameters[1];
        const double k = parameters[2];
        const double kRelativePressure = k * relativePressure;
        if (kRelativePressure >= 1.0)
        {
          throw std::domain_error("Error: GAB isotherm is undefined for K p/p_sat >= 1");
        }
        const double denominator =
            (1.0 - kRelativePressure) * (1.0 + kRelativePressure * (c - 1.0));
        return parameters[0] * c * kRelativePressure / denominator;
      }
      default:
        throw std::runtime_error("Error: unknown isotherm type");
    }
  }

  /**
   * \brief Out-of-line part of slope() for the isotherm types that need many registers.
   */
  [[gnu::noinline]] double slopeSlowPath(double pressure, double pH, double relativeHumidity,
                                         double waterLoading) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir_pH:
      {
        const double pHFactor = std::pow(10.0, std::clamp(parameters[2] - pH, -300.0, 300.0));
        const double b = parameters[1] / (1.0 + pHFactor);
        const double denominator = 1.0 + b * pressure;
        return parameters[0] * b / (denominator * denominator);
      }
      case Isotherm::Type::Unilan:
      {
        double upper = parameters[1] * std::exp(parameters[2]);
        double lower = parameters[1] * std::exp(-parameters[2]);
        return parameters[0] * (0.5 / parameters[2]) *
               (upper / (1.0 + upper * pressure) - lower / (1.0 + lower * pressure));
      }
      case Isotherm::Type::Langmuir_humidity:
      case Isotherm::Type::Langmuir_moisture_swing:
      case Isotherm::Type::GAB_Langmuir:
      {
        [[maybe_unused]] const auto [saturationLoading, affinity, exponent] =
            humidityParameters(relativeHumidity, waterLoading);
        if (std::isinf(affinity)) return 0.0;
        const double denominator = 1.0 + affinity * pressure;
        return saturationLoading * affinity / (denominator * denominator);
      }
      case Isotherm::Type::GAB_Sips:
      {
        const auto [saturationLoading, affinity, exponent] = humidityParameters(relativeHumidity, waterLoading);
        const double temp = std::pow(affinity * pressure, 1.0 / exponent);
        const double tempSlope =
            std::pow(affinity, 1.0 / exponent) * std::pow(pressure, 1.0 / exponent - 1.0) / exponent;
        return saturationLoading * tempSlope / ((1.0 + temp) * (1.0 + temp));
      }
      case Isotherm::Type::GAB_Toth:
      case Isotherm::Type::GAB_modified_Toth:
      case Isotherm::Type::Toth_exponential_humidity:
      {
        const auto [saturationLoading, affinity, exponent] = humidityParameters(relativeHumidity, waterLoading);
        const double temp = affinity * pressure;
        return saturationLoading * affinity * std::pow(1.0 + std::pow(temp, exponent), -(1.0 + exponent) / exponent);
      }
      case Isotherm::Type::GAB:
      {
        if (!(temperature > 0.0))
        {
          throw std::domain_error("Error: GAB needs a positive temperature (set it with setTemperature)");
        }

        // The same saturation pressure as valueSlowPath; below p = 0 the slope is its limit at p -> 0.
        const double saturationPressure = waterSaturationPressure(temperature);
        const double c = parameters[1];
        const double k = parameters[2];
        const double kRelativePressure = k * std::max(pressure, 0.0) / saturationPressure;
        if (kRelativePressure >= 1.0)
        {
          throw std::domain_error("Error: GAB isotherm is undefined for K p/p_sat >= 1");
        }
        const double denominator = (1.0 - kRelativePressure) * (1.0 + kRelativePressure * (c - 1.0));
        const double denominatorSlope = (c - 2.0) - 2.0 * (c - 1.0) * kRelativePressure;
        return parameters[0] * c * (denominator - kRelativePressure * denominatorSlope) / (denominator * denominator) *
               k / saturationPressure;
      }
      default:
        throw std::runtime_error("Error: unknown isotherm type");
    }
  }

  /**
   * \brief Effective saturation loading, affinity, and exponent of a humidity-dependent site at a relative humidity and
   *        water loading.
   *
   * The exponent is nu for the Sips and Toth types and 1 for the Langmuir types. The moisture-swing affinity is infinite
   * at RH = 0.
   */
  [[nodiscard]] std::array<double, 3> humidityParameters(double relativeHumidity, double waterLoading) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir_humidity:
      {
        const double waterAffinity = parameters[3] * relativeHumidity;
        return {parameters[0], parameters[1] * (1.0 + parameters[2] * waterAffinity / (1.0 + waterAffinity)), 1.0};
      }
      case Isotherm::Type::Langmuir_moisture_swing:
        return {parameters[0], parameters[1] * std::pow(relativeHumidity, -parameters[2]), 1.0};
      case Isotherm::Type::GAB_Langmuir:
        return {parameters[0] * std::max(0.0, 1.0 - parameters[3] * waterLoading),
                parameters[1] * (1.0 + parameters[2] * waterLoading), 1.0};
      case Isotherm::Type::GAB_Sips:
      case Isotherm::Type::GAB_Toth:
        return {parameters[0] * std::max(0.0, 1.0 - parameters[4] * waterLoading),
                parameters[1] * (1.0 + parameters[3] * waterLoading), parameters[2]};
      case Isotherm::Type::GAB_modified_Toth:
      {
        const double capacityReduction = parameters[3] * waterLoading;
        if (capacityReduction >= 1.0)
        {
          throw std::domain_error("Error: GAB-modified Toth is undefined for gamma q_w >= 1");
        }
        return {parameters[0] / (1.0 - capacityReduction), parameters[1] * (1.0 + parameters[4] * waterLoading),
                parameters[2]};
      }
      case Isotherm::Type::Toth_exponential_humidity:
        return {parameters[0] * std::exp(parameters[3] * relativeHumidity), parameters[1], parameters[2]};
      default:
        throw std::runtime_error("Error: the isotherm type does not depend on humidity");
    }
  }

  /**
   * \brief Out-of-line part of psiForPressure() for the isotherm types that need many registers.
   *
   * Keeping them out of the inline psiForPressure() spares the common types from saving all callee-saved registers on
   * every call.
   */
  [[gnu::noinline]] double psiForPressureSlowPath(double pressure) const
  {
    switch (type)
    {
      case Isotherm::Type::Langmuir_pH:
        throw std::runtime_error("Error: pH-dependent Langmuir is not compatible with spreading-pressure methods");
      case Isotherm::Type::Redlich_Peterson:
      {
        if (parameters[1] * std::pow(pressure, parameters[2]) < 1.0)
        {
          return parameters[0] * pressure *
                 hypergeometric2F1(1.0, 1.0 / parameters[2], 1.0 + 1.0 / parameters[2],
                                   -parameters[1] * std::pow(pressure, parameters[2]));
        }
        else
        {
          double prefactor = parameters[0] / parameters[2];
          double temp = M_PI / (std::pow(parameters[1], 1.0 / parameters[2]) * std::sin(M_PI * 1.0 / parameters[2]));

          double term1 = -1.0 / (parameters[1] * std::pow(pressure, parameters[2]));
          double numerator = 1.0;
          double sum = 0.0;
          // quickly converging sum
          for (size_t k = 1; k <= 15; k++)
          {
            numerator *= term1;
            sum += numerator / (static_cast<double>(k) * parameters[2] - 1.0);
          }
          return prefactor * (temp + pressure * parameters[2] * sum);
        }
      }
      case Isotherm::Type::Toth:
      {
        double temp = parameters[1] * pressure;
        double theta = temp / std::pow(1.0 + std::pow(temp, parameters[2]), 1.0 / parameters[2]);
        double theta_pow = std::pow(theta, parameters[2]);
        double reducedGrandPotential = parameters[0] * (theta - (theta / parameters[2]) * std::log(1.0 - theta_pow));

        // use the first 100 terms of the sum
        double temp1 = parameters[0] * theta;
        double temp2 = 0.0;
        for (size_t k = 1; k <= 100; ++k)
        {
          temp1 *= theta_pow;
          temp2 += parameters[2];
          reducedGrandPotential -= temp1 / (temp2 * (temp2 + 1.0));
        }

        return reducedGrandPotential;
      }
      case Isotherm::Type::Unilan:
      {
        return (0.5 * parameters[0] / parameters[2]) * (li2(-parameters[1] * std::exp(-parameters[2]) * pressure) -
                                                        li2(-parameters[1] * std::exp(parameters[2]) * pressure));
      }
      case Isotherm::Type::BingelWalton:
      {
        double start = 1e-14;
        size_t max_steps = 100;
        double acc = 1e-6;

        // Romberg integration: https://en.wikipedia.org/wiki/Romberg%27s_method
        std::vector<double> R1(max_steps), R2(max_steps);  // buffers
        double *Rp = &R1[0], *Rc = &R2[0];                 // Rp is previous row, Rc is current row
        double h = pressure - start;                       // step size
        Rp[0] = (value(start) / start + value(pressure) / pressure) * h * 0.5;  // first trapezoidal step

        for (size_t i = 1; i < max_steps; ++i)
        {
          h /= 2.0;
          double c = 0;
          size_t ep = size_t{1} << (i - 1);  // 2^(n-1)
          for (size_t j = 1; j <= ep; ++j)
          {
            c +=
                value(start + static_cast<double>(2 * j - 1) * h) / (start + static_cast<double>(2 * j - 1) * h);
          }
          Rc[0] = h * c + 0.5 * Rp[0];  // R(i,0)

          double n_k = 1.0;
          for (size_t j = 1; j <= i; ++j)
          {
            n_k *= 4.0;                                         // 4^j, exact in double
            Rc[j] = (n_k * Rc[j - 1] - Rp[j - 1]) / (n_k - 1);  // compute R(i,j)
          }

          if (i > 1 && std::fabs(Rp[i - 1] - Rc[i]) < acc)
          {
            return Rc[i];
          }

          double* rt = Rp;
          Rp = Rc;
          Rc = rt;
        }
        return Rp[max_steps - 1];  // return our best guess
      }
      case Isotherm::Type::GAB:
        throw std::runtime_error("Error: GAB is compatible only with pure-isotherm methods such as SPI");
      case Isotherm::Type::Langmuir_humidity:
      case Isotherm::Type::Langmuir_moisture_swing:
      case Isotherm::Type::GAB_Langmuir:
      case Isotherm::Type::GAB_Sips:
      case Isotherm::Type::GAB_Toth:
      case Isotherm::Type::GAB_modified_Toth:
      case Isotherm::Type::Toth_exponential_humidity:
        throw std::runtime_error("Error: humidity-dependent isotherms are compatible only with SPI");
      default:
        throw std::runtime_error("Error: unknown isotherm type");
    }
  }

  /**
   * \brief Returns whether psiForPressure() is defined, which it is not for pH-Langmuir, GAB, and the
   *        humidity-dependent types.
   */
  [[nodiscard]] bool definesGrandPotential() const noexcept
  {
    return type != Type::Langmuir_pH && type != Type::GAB && !isHumidityDependent(type);
  }

  /**
   * \brief Randomizes the isotherm parameters within specified bounds.
   *
   * Sets the isotherm parameters to random values within physically meaningful ranges,
   * based on the specified maximum loading.
   *
   * \param maximumLoading The maximum adsorption loading used to scale the randomized parameters.
   */
  void randomize(double maximumLoading);

  /**
   * \brief Checks if the isotherm parameters are physically meaningful.
   *
   * Determines whether the isotherm parameters are within acceptable physical ranges.
   *
   * \return True if the parameters are unphysical; false otherwise.
   */
  bool isUnphysical() const;
};

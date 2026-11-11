#pragma once

#include <cmath>
#include <stdexcept>

/**
 * \brief Saturation pressure of liquid water, which turns a water partial pressure into a relative humidity.
 *
 * The models are
 *
 * - Antoine:  10^(8.07131 - 1730.63/(T - 39.724)) mmHg, the Kelvin form of the Antoine correlation for 274-373 K
 *             (the default, and the only choice of GAB before WaterSaturationPressureModel existed);
 * - Wagner:   the IAPWS saturation-pressure equation of Wagner and Pruss (1993), from 273.16 K to the critical
 *             temperature 647.096 K;
 * - Constant: a fixed pressure in Pa.
 */
struct WaterSaturationPressure
{
  /**
   * \brief The saturation-pressure correlation.
   */
  enum class Model
  {
    Antoine = 0,  ///< Antoine correlation for 274-373 K
    Wagner = 1,   ///< IAPWS equation of Wagner and Pruss
    Constant = 2  ///< Fixed pressure
  };

  Model model{Model::Antoine};
  double constantPressure{0.0};  ///< Saturation pressure in Pa for Model::Constant.

  /**
   * \brief Returns the saturation pressure in Pa at a temperature in K.
   */
  double operator()(double temperature) const
  {
    switch (model)
    {
      case Model::Antoine:
      default:
        // The result is converted from mmHg to Pa.
        return std::pow(10.0, 8.07131 - 1730.63 / (temperature - 39.724)) * 133.32236842105263;
      case Model::Wagner:
      {
        const double criticalTemperature = 647.096;
        const double criticalPressure = 22.064e6;
        if (!(temperature > 0.0 && temperature < criticalTemperature))
        {
          throw std::domain_error("Error: the Wagner saturation pressure of water needs 0 < T < 647.096 K");
        }
        const double tau = 1.0 - temperature / criticalTemperature;
        const double sqrtTau = std::sqrt(tau);
        const double tau3 = tau * tau * tau;
        const double sum = -7.85951783 * tau + 1.84408259 * tau * sqrtTau - 11.7866497 * tau3 +
                           22.6807411 * tau3 * sqrtTau - 15.9618719 * tau3 * tau +
                           1.80122502 * tau3 * tau3 * tau * sqrtTau;
        return criticalPressure * std::exp(criticalTemperature / temperature * sum);
      }
      case Model::Constant:
        return constantPressure;
    }
  }
};

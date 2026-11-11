#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "block_banded.h"
#include "column.h"

/**
 * \brief Banded Jacobian of the fluid concentrations and physisorbed loadings ("Analytical sorption Jacobian" in the
 *        supplementary information of the paper).
 *
 * The unknowns are node-major, (c_1, ..., c_N, q_1, ..., q_N) at every node, with node 0 the prescribed inlet
 * (with a Danckwerts inlet, the rows next to the inlet couple to the inlet-face value of computeMassDerivatives).
 * The transport rows are the exact derivatives of computeMassDerivatives at the velocity of the last precompute:
 * upwind convection, or WENO3 with its nonlinear weights and, in the gas phase, the mole fractions and the upwind
 * total molar flux; the equilibrium enters through computePhysisorptionEquilibriumSensitivities at a fixed total
 * pressure (not for MPD, which has no sensitivities: there the loading rows hold the mass-transfer rate only).
 * Chemisorption, reactions, pore transport, the velocity closures, and the energy balance are not included.
 *
 * The matrix is kept in LAPACK general band storage with room for the fill-in of dgbtrf: with b sub- and
 * superdiagonals (2N with upwind convection, 5N with WENO3, whose stencil reaches two nodes upstream) and 3b + 1 rows,
 * entry (row, col) is band[2b + row - col + col (3b + 1)].
 */
struct SorptionJacobian
{
  explicit SorptionJacobian(const Column& column);

  size_t numberOfComponents;                   ///< N.
  size_t numberOfNodes;                        ///< numberOfGridPoints + 1.
  size_t size;                                 ///< Number of unknowns, 2 N numberOfNodes.
  size_t bandwidth;                            ///< Number of sub- and superdiagonals.
  size_t leadingDimension;                     ///< Rows of the band storage, 3 bandwidth + 1.
  std::vector<double> band;                    ///< Column-major band storage, see above.
  std::vector<double> equilibriumSensitivity;  ///< dq*_i/dc_l per node, column-major N x N blocks.
  std::vector<double> mixtureSensitivity;      ///< Scratch for one mixture, N x N.

  double& operator()(size_t row, size_t col) noexcept
  {
    return band[2 * bandwidth + row - col + col * leadingDimension];
  }
  double operator()(size_t row, size_t col) const noexcept
  {
    return band[2 * bandwidth + row - col + col * leadingDimension];
  }

  /**
   * \brief Entry (row, col), zero outside the band.
   */
  double entry(size_t row, size_t col) const noexcept
  {
    return (row > col ? row - col : col - row) <= bandwidth ? (*this)(row, col) : 0.0;
  }

  size_t concentrationIndex(size_t node, size_t component) const noexcept
  {
    return 2 * node * numberOfComponents + component;
  }
  size_t loadingIndex(size_t node, size_t component) const noexcept
  {
    return (2 * node + 1) * numberOfComponents + component;
  }

  /**
   * \brief Copies the field-major concentration and loading blocks of a column into a node-major vector.
   */
  void gather(std::span<const double> concentration, std::span<const double> loading,
              std::span<double> nodeMajor) const noexcept;

  /**
   * \brief Copies a node-major vector back into field-major concentration and loading blocks.
   */
  void scatter(std::span<const double> nodeMajor, std::span<double> concentration,
               std::span<double> loading) const noexcept;
};

/**
 * \brief Assembles the Jacobian at the given concentrations and gas temperatures; call after the right-hand side at
 *        that state (precompute, or the CVODE right-hand side), which leaves the velocity, total pressure, mole
 *        fractions, and the hypothetical pressures of the equilibrium sensitivities in the column.
 */
void assembleSorptionJacobian(Column& column, std::span<const double> concentration,
                              std::span<const double> gasTemperature, SorptionJacobian& jacobian);

/**
 * \brief Jacobian of the gas, solid, and wall temperatures of computeEnergyDerivatives with respect to themselves, as a
 *        block-tridiagonal matrix of 3 x 3 node blocks (T_g, T_s, T_w): the heat exchange between the three, the
 *        upwind convection and conduction of the gas temperature, and the conduction of the wall. The gas density,
 *        velocity, and conductivity are held fixed, and the coupling to the loadings (heat of adsorption, temperature
 *        dependence of the isotherms) is left out. Call after computeDerivatives, which sets the conduction
 *        coefficients of the column.
 */
void assembleThermalJacobian(const Column& column, BlockBandedMatrix& jacobian);

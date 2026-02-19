#include "rk3_si.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <vector>

#include "utils.h"

extern "C"
{
  void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
  void dgbsv_(int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
}

bool SemiImplicitRungeKutta3::propagate(Column& column, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;
  Column newcolumn(column);

  std::vector<double> solved(column.partialPressure.size());

  std::vector<double> implicitInvKLs(Ncomp);
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    implicitInvKLs[comp] = 1.0 / (1.0 + timeStep * column.components[comp].Kl);
  }

  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance = std::max(tolerance, std::abs((column.partialPressure[Ngrid * Ncomp + j] /
                                                (column.exitPressure * column.components[j].Yi0)) -
                                               1.0));
    }

    // consider 1% as being visibily indistinguishable from 'converged'
    // use a 10% longer time for display purposes
    if (tolerance < 0.01)
    {
      std::cout << "\nConvergence criteria reached, running 10% longer\n\n" << std::endl;
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoSteps = false;
    }
  }

  // SSP-RK Step 1
  // ======================================================================

  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P

  // Dqdt and Dpdt are calculated at old time step
  // make estimate for the new loadings and new gas phase partial pressures
  // first iteration is made using the Explicit Euler scheme

  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;
    newcolumn.adsorption[i] =
        (column.adsorption[i] + timeStep * kl * column.equilibriumAdsorption[i]) * implicitInvKLs[comp];
  }

  computeVelocity(newcolumn);
  computePressureUpdateMatrix(newcolumn, timeStep, solved);

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.partialPressure[i] = solved[i];
  }

  computeEquilibriumLoadings(newcolumn);

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep

  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;
    newcolumn.adsorption[i] =
        0.75 * column.adsorption[i] +
        0.25 * (newcolumn.adsorption[i] + timeStep * kl * column.equilibriumAdsorption[i]) * implicitInvKLs[comp];
  }

  computeVelocity(newcolumn);
  computePressureUpdateMatrix(newcolumn, timeStep, solved);

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.partialPressure[i] = 0.75 * column.partialPressure[i] + 0.25 * solved[i];
  }

  computeEquilibriumLoadings(newcolumn);

  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep

  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;
    newcolumn.adsorption[i] = (1.0 / 3.0) * column.adsorption[i] +
                              (2.0 / 3.0) *
                                  (newcolumn.adsorption[i] + timeStep * kl * column.equilibriumAdsorption[i]) *
                                  implicitInvKLs[comp];
  }

  computeVelocity(newcolumn);
  computePressureUpdateMatrix(newcolumn, timeStep, solved);

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.partialPressure[i] = (1.0 / 3.0) * column.partialPressure[i] + (2.0 / 3.0) * solved[i];
  }

  computeEquilibriumLoadings(newcolumn);

  // update to the new time step
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;

    newcolumn.adsorption[i] =
        (newcolumn.adsorption[i] + timeStep * timeStep * kl * kl * newcolumn.equilibriumAdsorption[i]) /
        (1.0 + timeStep * timeStep * kl * kl);
  }

  computeVelocity(newcolumn);

  computePressureUpdateMatrixFinal(newcolumn, timeStep, solved);
  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.partialPressure[i] = solved[i];
  }

  column = newcolumn;

  // pulse boundary condition
  for (size_t j = 0; j < Ncomp; ++j)
  {
    if (column.pulse && t > column.pulseTime)
    {
      if (j == column.carrierGasComponent)
      {
        column.partialPressure[0 * Ncomp + j] = column.externalPressure;
      }
      else
      {
        column.partialPressure[0 * Ncomp + j] = 0.0;
      }
    }
    else
    {
      column.partialPressure[0 * Ncomp + j] = column.externalPressure * column.components[j].Yi0;
    }
  }

  return (!autoSteps && step >= numberOfSteps - 1);
}

void computePressureFirstDerivative(size_t Ncomp, size_t Ngrid, double resolution,
                                    std::span<const double> partialPressure, std::span<double> partialPressureDot,
                                    std::span<const double> interstitialGasVelocity,
                                    const std::vector<Component>& components)
{
  double idx = 1.0 / resolution;
  double idx2 = idx * idx;

  std::mdspan<const double, std::dextents<size_t, 2>> spanPartialPressure(partialPressure.data(), Ngrid + 1, Ncomp);

  // first gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    partialPressureDot[comp] = 0.0;
  }

  // middle gridpoints
  for (size_t grid = 1; grid < Ngrid; ++grid)
  {
    for (size_t comp = 0; comp < Ncomp; ++comp)
    {
      partialPressureDot[grid * Ncomp + comp] =
          (interstitialGasVelocity[grid - 1] * spanPartialPressure[grid - 1, comp] -
           interstitialGasVelocity[grid] * spanPartialPressure[grid, comp]) *
              idx +
          components[comp].D *
              (spanPartialPressure[grid + 1, comp] - 2.0 * spanPartialPressure[grid, comp] +
               spanPartialPressure[grid - 1, comp]) *
              idx2;
    }
  }

  // last gridpoint
  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    partialPressureDot[Ngrid * Ncomp + comp] =
        (interstitialGasVelocity[Ngrid - 1] * spanPartialPressure[Ngrid - 1, comp] -
         interstitialGasVelocity[Ngrid] * spanPartialPressure[Ngrid, comp]) *
            idx +
        components[comp].D * (spanPartialPressure[Ngrid - 1, comp] - spanPartialPressure[Ngrid, comp]) * idx2;
  }
}

// u' = f(u) + g(u) u
// u^{i+1} = u + dt * f / (1 + dt * g)

// p' = f(p) + g(p) p
// f(p) = (q - qeq)
// g(p) = (D / dx) * v_i-1 * p_i-1 + v_i + v_i+1
// A = (D / dx)
// | v_i v_{i+1} 0 0 0 0 |
// | v_{i-1} v_i v_{i+1} 0 0 0 |
// | 0 v_{i-1} v_i v_{i+1} 0 0 |

// p^{i+1} = p + dt * f / (1 + dt * A)
// (1 + dt * A)^{-1} p^{i+1} = p + dt * f

void computePressureUpdateMatrix(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  int n = static_cast<int>(Ngrid + 1);
  int nrhs = 1;
  int info = 0;

  std::vector<double> upper(Ngrid);
  std::vector<double> lower(Ngrid);
  std::vector<double> diag(Ngrid + 1);
  std::vector<double> rhs(Ngrid + 1);

  // pdot = Ap - b(q_eq - q)
  // diag: Aii = -(v_i / dx + 2D / (dx)^2)
  // lower: Ai,i-1 = (v_{i-1} / dx + D / (dx)^2)
  // upper: D / (dx)^2
  // SIRK3: (1 - dt * A) p^{i+1} = p^i + dt * b * (q_eq - q)
  // Matrix: I - dt * A

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double D = column.components[comp].D;
    double b = column.prefactorMassTransfer[comp];
    for (size_t i = 0; i < Ngrid; ++i)
    {
      lower[i] = -timeStep * (column.interstitialGasVelocity[i] * idx + D * idx2);
      upper[i] = -timeStep * (D * idx2);
    }
    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      diag[i] = 1.0 + timeStep * (column.interstitialGasVelocity[i] * idx + 2.0 * D * idx2);
      rhs[i] = column.partialPressure[i * Ncomp + comp] +
               timeStep * (-b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]));
    }

    // boundary conditions
    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.partialPressure[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

void computePressureUpdateMatrixEnergyBalance(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  int n = static_cast<int>(Ngrid + 1);
  int nrhs = 1;
  int info = 0;

  std::vector<double> upper(Ngrid);
  std::vector<double> lower(Ngrid);
  std::vector<double> diag(Ngrid + 1);
  std::vector<double> rhs(Ngrid + 1);

  for (size_t grid = 0; grid < Ngrid; ++grid)
  {
    column.facePressures[grid] = 0.5 * (column.totalPressure[grid] + column.totalPressure[grid + 1]);
  }

  // pdot = Ap - b(q_eq - q)
  // diag: Aii = -(v_i / dx + 2D / (dx)^2)
  // lower: Ai,i-1 = (v_{i-1} / dx + D / (dx)^2)
  // upper: D / (dx)^2
  // SIRK3: (1 - dt * A) p^{i+1} = p^i + dt * b * (q_eq - q)
  // Matrix: I - dt * A

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double D = column.components[comp].D;
    double b = column.prefactorMassTransfer[comp];
    for (size_t i = 0; i < Ngrid; ++i)
    {
      lower[i] = -timeStep *
                 ((column.gasTemperature[i + 1] / column.gasTemperature[i]) * column.interstitialGasVelocity[i] * idx +
                  D * idx2 *
                      ((column.facePressures[i] / column.totalPressure[i]) +
                       (column.totalPressure[i + 1] / column.gasTemperature[i + 1]) *
                           (column.gasTemperature[i] - column.gasTemperature[i + 1])));

      upper[i] = -timeStep * (D * idx2 * (column.facePressures[i] / column.totalPressure[i + 1]));
    }
    for (size_t i = 1; i < Ngrid; ++i)
    {
      // FIX: gasTemperatureDot is not computed!!!
      diag[i] =
          1.0 - timeStep * (-column.interstitialGasVelocity[i] * idx +
                            (column.gasTemperatureDot[i] / column.gasTemperature[i]) +
                            D * idx2 *
                                (((column.facePressures[i - 1] + column.facePressures[i]) / column.totalPressure[i]) +
                                 (column.totalPressure[i] / column.gasTemperature[i]) *
                                     (column.gasTemperature[i - 1] - column.gasTemperature[i])));
      rhs[i] = column.partialPressure[i * Ncomp + comp] +
               timeStep * (-b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]));
    }

    diag[Ngrid] = 1.0 - timeStep * (-column.interstitialGasVelocity[Ngrid] * idx +
                                    (column.gasTemperatureDot[Ngrid] / column.gasTemperature[Ngrid]) +
                                    D * idx2 *
                                        ((column.facePressures[Ngrid - 1] / column.totalPressure[Ngrid]) +
                                         (column.totalPressure[Ngrid] / column.gasTemperature[Ngrid]) *
                                             (column.gasTemperature[Ngrid - 1] - column.gasTemperature[Ngrid])));

    // boundary conditions
    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.partialPressure[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

void computeEnergyUpdateMatrixEnergyBalance(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  double relativeVolume = ((1 - column.voidFraction) / column.voidFraction);
  double accessibleSurface = 6.0 / column.particleDiameter;  // prefactor 6.0 in python and 2.0 in eqs
  double prefactorGasSolid = relativeVolume * accessibleSurface * column.heatTransferGasSolid / column.heatCapacityGas;

  // is this extra 1/eps necessary? It's in python not in eqs
  double prefactorGasWall =
      4.0 * column.heatTransferGasWall / (column.voidFraction * column.heatCapacityGas * column.internalDiameter);
  double prefactorGasGas = -(prefactorGasSolid + prefactorGasWall);

  double coeffSolidGas =
      accessibleSurface * column.heatTransferGasSolid / (column.heatCapacitySolid * column.particleDensity);
  double coeffSolidSolid = -coeffSolidGas;

  // prefactor 4 in python, 2 in eqs
  double internalArea =
      4.0 * column.internalDiameter /
      (column.outerDiameter * column.outerDiameter - column.internalDiameter * column.internalDiameter);
  double externalArea =
      4.0 * column.outerDiameter /
      (column.outerDiameter * column.outerDiameter - column.internalDiameter * column.internalDiameter);
  double invHeatDensityWall = 1.0 / (column.heatCapacityWall * column.wallDensity);
  double coeffWallGas = column.heatTransferGasWall * internalArea * invHeatDensityWall;
  double coeffWallWall = -coeffWallGas - column.heatTransferWallExternal * externalArea * invHeatDensityWall;
  double coeffDiffusion = column.gasThermalConductivity / column.heatCapacityGas;

  int n = static_cast<int>(Ngrid + 1);
  int nrhs = 1;
  int info = 0;

  //   Eigen::SparseMatrix<double> A(9 * n * n);
  //   std::vector<Eigen::Triplet<double>> trips;

  //   for (size_t grid = 0; grid < Ngrid + 1; ++grid)
  //   {
  //     column.gasDensity[grid] = 0.0;
  //     for (size_t comp = 0; comp < Ncomp; ++comp)
  //     {
  //       column.gasDensity[grid] += column.components[comp].molecularWeight *
  //                                  std::max(0.0, column.partialPressure[grid * Ncomp + comp]) /
  //                                  (R * column.gasTemperature[grid]);
  //     }
  //   }

  //   for (size_t grid = 1; grid < Ngrid; ++grid)
  //   {
  //     // Ggg
  //     double invGasDensity = 1.0 / column.gasDensity[grid];
  //     trips.emplace_back(
  //         0 * n + grid, 0 * n + grid - 1,
  //         -timeStep * (idx2 * coeffDiffusion * invGasDensity + idx * column.interstitialGasVelocity[grid]));
  //     trips.emplace_back(0 * n + grid, 0 * n + grid,
  //                        1.0 - timeStep * (-(idx2 * coeffDiffusion - prefactorGasGas) * invGasDensity -
  //                                          idx * column.interstitialGasVelocity[grid]));
  //     trips.emplace_back(0 * n + grid, 0 * n + grid + 1, -timeStep * idx2 * coeffDiffusion * invGasDensity);

  //     // Ggs
  //     trips.emplace_back(0 * n + grid, 1 * n + grid, -timeStep * prefactorGasSolid / column.gasDensity[grid]);

  //     // Ggw
  //     trips.emplace_back(0 * n + grid, 2 * n + grid, -timeStep * prefactorGasWall / column.gasDensity[grid]);

  //     // Gsg
  //     trips.emplace_back(1 * n + grid, 0 * n + grid, -timeStep * coeffSolidGas);

  //     // Gss
  //     trips.emplace_back(1 * n + grid, 1 * n + grid, 1.0 - timeStep * coeffSolidSolid);

  //     // Gsw: 0

  //     // Gwg
  //     trips.emplace_back(2 * n + grid, 0 * n + grid, -timeStep * coeffWallGas);

  //     // Gws: 0

  //     // Gww
  //     trips.emplace_back(2 * n + grid, 2 * n + grid - 1, -timeStep * ());
  //   }
}

void computePressureUpdateMatrixFinal(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  double dt2 = timeStep * timeStep;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  std::vector<double> upper(Ngrid);
  std::vector<double> lower(Ngrid);
  std::vector<double> diag(Ngrid + 1);
  std::vector<double> preRHS(Ngrid + 1);
  std::vector<double> rhs(Ngrid + 1);

  // rows:
  // 0 empty
  // 1 empty
  // 2 second superdiagonal
  // 3 first superdiagonal
  // 4 main diagonal
  // 5 first subdiagonal
  // 6 second subdiagonal
  std::vector<double> ab(7 * (Ngrid + 1), 0.0);
  std::mdspan abS(ab.data(), 7, Ngrid + 1);
  std::vector<int> ipiv(Ngrid + 1);

  int n = static_cast<int>(Ngrid + 1);
  int info = 0;
  int ku = 2;
  int kl = 2;
  int ldab = 7;
  int nrhs = 1;

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double D = column.components[comp].D;
    double b = column.prefactorMassTransfer[comp];

    // First create tridiagonal matrix A and f
    for (size_t i = 0; i < Ngrid; ++i)
    {
      lower[i] = column.interstitialGasVelocity[i] * idx + D * idx2;
      upper[i] = D * idx2;
    }
    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      diag[i] = -(column.interstitialGasVelocity[i] * idx + 2.0 * D * idx2);
      preRHS[i] = -b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]);
    }

    // Now calculate (1 + (dt * A)^2) and (p - dt * A * f)
    for (size_t i = 0; i < Ngrid - 1; ++i)
    {
      abS[6, i] = dt2 * lower[i + 1] * lower[i];
      abS[2, i + 2] = dt2 * upper[i] * upper[i + 1];
    }
    for (size_t i = 0; i < Ngrid; ++i)
    {
      abS[5, i] = dt2 * lower[i] * (diag[i] + diag[i + 1]);
      abS[3, i + 1] = dt2 * upper[i] * (diag[i] + diag[i + 1]);
    }

    // skip the first row as these are set in boundary conditions anyway

    for (size_t i = 1; i < Ngrid; ++i)
    {
      abS[4, i] = 1.0 + dt2 * (lower[i - 1] * upper[i - 1] + diag[i] * diag[i] + lower[i] * upper[i]);
      rhs[i] = column.partialPressure[i * Ncomp + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }
    abS[4, Ngrid] = 1.0 + dt2 * (lower[Ngrid - 1] * upper[Ngrid - 1] + diag[Ngrid] * diag[Ngrid]);
    rhs[Ngrid] = column.partialPressure[Ngrid * Ncomp + comp] -
                 dt2 * (lower[Ngrid - 1] * preRHS[Ngrid - 1] + diag[Ngrid] * preRHS[Ngrid]);

    // boundary conditions
    abS[4, 0] = 1.0;
    abS[3, 0] = 0.0;
    abS[2, 0] = 0.0;

    rhs[0] = column.partialPressure[comp];

    dgbsv_(&n, &kl, &ku, &nrhs, ab.data(), &ldab, ipiv.data(), rhs.data(), &n, &info);

    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

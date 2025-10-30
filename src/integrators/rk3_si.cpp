
#include "rk3_si.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <vector>

extern "C"
{
  void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
  void dgbsv_(int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
}

bool SemiImplicitRungeKutta3::propagate(Column& state, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  size_t Ngrid = state.Ngrid;
  size_t Ncomp = state.Ncomp;
  Column newState(state);

  std::vector<double> solved(state.partialPressure.size());

  if (autoSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < Ncomp; ++j)
    {
      tolerance = std::max(
          tolerance,
          std::abs((state.partialPressure[Ngrid * Ncomp + j] / (state.exitPressure * state.components[j].Yi0)) - 1.0));
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

  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = state.components[comp].Kl;
    newState.adsorption[i] =
        (state.adsorption[i] + timeStep * kl * state.equilibriumAdsorption[i]) / (1.0 + timeStep * kl);
  }

  computeVelocity(newState);
  computePressureUpdateMatrix(newState, timeStep, solved);

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newState.partialPressure[i] = solved[i];
  }


  computeEquilibriumLoadings(newState);

  // SSP-RK Step 2
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep

  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = state.components[comp].Kl;
    newState.adsorption[i] =
        0.75 * state.adsorption[i] +
        0.25 * (newState.adsorption[i] + timeStep * kl * state.equilibriumAdsorption[i]) / (1.0 + timeStep * kl);
  }

  computeVelocity(newState);
  computePressureUpdateMatrix(newState, timeStep, solved);

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newState.partialPressure[i] = 0.75 * state.partialPressure[i] + 0.25 * solved[i];
  }

  computeEquilibriumLoadings(newState);

  // SSP-RK Step 3
  // ======================================================================

  // calculate new derivatives at new (current) timestep
  // calculate the derivatives Dq/dt and Dp/dt based on Qeq, Q, V, and P at new (current) timestep

  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = state.components[comp].Kl;
    newState.adsorption[i] =
        (1.0 / 3.0) * state.adsorption[i] +
        (2.0 / 3.0) * (newState.adsorption[i] + timeStep * kl * state.equilibriumAdsorption[i]) / (1.0 + timeStep * kl);
  }

  computeVelocity(newState);
  computePressureUpdateMatrix(newState, timeStep, solved);

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newState.partialPressure[i] = (1.0 / 3.0) * state.partialPressure[i] + (2.0 / 3.0) * solved[i];
  }


  computeEquilibriumLoadings(newState);

  // update to the new time step
  for (size_t i = 0; i < state.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = state.components[comp].Kl;

    newState.adsorption[i] =
        (newState.adsorption[i] + timeStep * timeStep * kl * kl * newState.equilibriumAdsorption[i]) /
        (1.0 + timeStep * timeStep * kl * kl);
  }

  computeVelocity(newState);

  computePressureUpdateMatrixFinal(newState, timeStep, solved);
  for (size_t i = 0; i < solved.size(); ++i)
  {
    newState.partialPressure[i] = solved[i];
  }


  state = newState;

  // pulse boundary condition
  if (state.pulse)
  {
    if (t > state.pulseTime)
    {
      for (size_t j = 0; j < Ncomp; ++j)
      {
        if (j == state.carrierGasComponent)
        {
          state.partialPressure[0 * Ncomp + j] = state.externalPressure;
        }
        else
        {
          state.partialPressure[0 * Ncomp + j] = 0.0;
        }
      }
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

void computePressureUpdateMatrix(Column& state, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / state.resolution;
  double idx2 = idx * idx;
  size_t Ngrid = state.Ngrid;
  size_t Ncomp = state.Ncomp;

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
  //

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    double D = state.components[comp].D;
    double b = state.prefactorMassTransfer[comp];
    for (size_t i = 0; i < Ngrid; ++i)
    {
      lower[i] = -timeStep * (state.interstitialGasVelocity[i] * idx + D * idx2);
      upper[i] = -timeStep * (D * idx2);
    }
    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      diag[i] = 1.0 + timeStep * (state.interstitialGasVelocity[i] * idx + 2.0 * D * idx2);
      rhs[i] = state.partialPressure[i * Ncomp + comp] +
               timeStep * (-b * (state.equilibriumAdsorption[i * Ncomp + comp] - state.adsorption[i * Ncomp + comp]));
    }

    // boundary conditions
    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = state.partialPressure[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

void computePressureUpdateMatrixFinal(Column& state, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / state.resolution;
  double idx2 = idx * idx;
  double dt2 = timeStep * timeStep;
  size_t Ngrid = state.Ngrid;
  size_t Ncomp = state.Ncomp;

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
    double D = state.components[comp].D;
    double b = state.prefactorMassTransfer[comp];

    // First create tridiagonal matrix A and f 
    for (size_t i = 0; i < Ngrid; ++i)
    {
      lower[i] = state.interstitialGasVelocity[i] * idx + D * idx2;
      upper[i] = D * idx2;
    }
    for (size_t i = 0; i < Ngrid + 1; ++i)
    {
      diag[i] = -(state.interstitialGasVelocity[i] * idx + 2.0 * D * idx2);
      preRHS[i] = -b * (state.equilibriumAdsorption[i * Ncomp + comp] - state.adsorption[i * Ncomp + comp]);
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
      rhs[i] = state.partialPressure[i * Ncomp + comp]  -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }
    abS[4, Ngrid] = 1.0 + dt2 * (lower[Ngrid - 1] * upper[Ngrid - 1] + diag[Ngrid] * diag[Ngrid]);
    rhs[Ngrid] = state.partialPressure[Ngrid * Ncomp + comp]  -
                 dt2 * (lower[Ngrid - 1] * preRHS[Ngrid - 1] + diag[Ngrid] * preRHS[Ngrid]);

    // boundary conditions
    abS[4, 0] = 1.0;
    abS[3, 0] = 0.0;
    abS[2, 0] = 0.0;

    rhs[0] = state.partialPressure[comp];

    dgbsv_(&n, &kl, &ku, &nrhs, ab.data(), &ldab, ipiv.data(), rhs.data(), &n, &info);

    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

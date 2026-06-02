#include "rk3_si.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mdspan>
#include <stdexcept>
#include <vector>

#include "compute.h"
#include "utils.h"

extern "C"
{
  void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
  void dgbsv_(int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
}

void computeConcentrationUpdateMatrix(Column& column, double timeStep, std::vector<double>& solved);
void computeConcentrationUpdateMatrixEnergyBalance(Column& column, double timeStep, std::vector<double>& solved);
void computeConcentrationUpdateMatrixFinal(Column& column, double timeStep, std::vector<double>& solved);
void computeConcentrationUpdateMatrixFinalEnergyBalance(Column& column, double timeStep, std::vector<double>& solved);

bool SemiImplicitRungeKutta3::propagate(Column& column, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;
  Column newcolumn(column);

  std::vector<double> solved(column.concentration.size());

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
      tolerance =
          std::max(tolerance, std::abs((column.moleFraction[Ngrid * Ncomp + j] / column.components[j].Yi0) - 1.0));
    }

    if (tolerance < 0.01)
    {
      std::cout << "\nConvergence criteria reached, running 10% longer\n\n" << std::endl;
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoSteps = false;
    }
  }

  // SSP-RK Step 1
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;
    newcolumn.adsorption[i] =
        (column.adsorption[i] + timeStep * kl * column.equilibriumAdsorption[i]) * implicitInvKLs[comp];
  }

  computeVelocity(newcolumn);
  if (newcolumn.energyBalance)
  {
    computeConcentrationUpdateMatrixEnergyBalance(newcolumn, timeStep, solved);
  }
  else
  {
    computeConcentrationUpdateMatrix(newcolumn, timeStep, solved);
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.concentration[i] = solved[i];
  }

  computePressure(newcolumn);
  computeEquilibriumLoadings(newcolumn);

  // SSP-RK Step 2
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;
    newcolumn.adsorption[i] =
        0.75 * column.adsorption[i] +
        0.25 * (newcolumn.adsorption[i] + timeStep * kl * column.equilibriumAdsorption[i]) * implicitInvKLs[comp];
  }

  computeVelocity(newcolumn);
  if (newcolumn.energyBalance)
  {
    computeConcentrationUpdateMatrixEnergyBalance(newcolumn, timeStep, solved);
  }
  else
  {
    computeConcentrationUpdateMatrix(newcolumn, timeStep, solved);
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.concentration[i] = 0.75 * column.concentration[i] + 0.25 * solved[i];
  }

  computePressure(newcolumn);
  computeEquilibriumLoadings(newcolumn);

  // SSP-RK Step 3
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
  if (newcolumn.energyBalance)
  {
    computeConcentrationUpdateMatrixEnergyBalance(newcolumn, timeStep, solved);
  }
  else
  {
    computeConcentrationUpdateMatrix(newcolumn, timeStep, solved);
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.concentration[i] = (1.0 / 3.0) * column.concentration[i] + (2.0 / 3.0) * solved[i];
  }

  computePressure(newcolumn);
  computeEquilibriumLoadings(newcolumn);

  // final implicit adsorption update
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % Ncomp;
    double kl = column.components[comp].Kl;

    newcolumn.adsorption[i] =
        (newcolumn.adsorption[i] + timeStep * timeStep * kl * kl * newcolumn.equilibriumAdsorption[i]) /
        (1.0 + timeStep * timeStep * kl * kl);
  }

  computeVelocity(newcolumn);

  if (newcolumn.energyBalance)
  {
    computeConcentrationUpdateMatrixFinalEnergyBalance(newcolumn, timeStep, solved);
  }
  else
  {
    computeConcentrationUpdateMatrixFinal(newcolumn, timeStep, solved);
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newcolumn.concentration[i] = solved[i];
  }

  computePressure(newcolumn);
  computeEquilibriumLoadings(newcolumn);

  column = newcolumn;
  enforceBoundaryCondition(column);
  return (!autoSteps && step >= numberOfSteps - 1);
}

void computeConcentrationUpdateMatrix(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  int n = static_cast<int>(Ngrid + 1);
  int nrhs = 1;

  std::vector<double> upper(Ngrid);
  std::vector<double> lower(Ngrid);
  std::vector<double> diag(Ngrid + 1);
  std::vector<double> rhs(Ngrid + 1);

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    int info = 0;

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
      rhs[i] = column.concentration[i * Ncomp + comp] +
               timeStep * (-b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]));
    }

    diag[Ngrid] = 1.0 + timeStep * (column.interstitialGasVelocity[Ngrid] * idx + D * idx2);

    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.concentration[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgtsv failed in computeConcentrationUpdateMatrix");
    }

    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

void computeConcentrationUpdateMatrixEnergyBalance(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

  int n = static_cast<int>(Ngrid + 1);
  int nrhs = 1;

  std::vector<double> upper(Ngrid);
  std::vector<double> lower(Ngrid);
  std::vector<double> diag(Ngrid + 1);
  std::vector<double> rhs(Ngrid + 1);

  for (size_t grid = 0; grid < Ngrid; ++grid)
  {
    column.facePressures[grid] = 0.5 * (column.totalPressure[grid] + column.totalPressure[grid + 1]) / R;
  }

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    int info = 0;

    double D = column.components[comp].D;
    double b = column.prefactorMassTransfer[comp];

    for (size_t i = 0; i < Ngrid; ++i)
    {
      double invGasTemperatureUpper = 1.0 / std::max(1e-10, column.gasTemperature[i]);
      double invTotalConcentrationUpper = 1.0 / std::max(1e-10, column.totalConcentration[i + 1]);
      double invGasTemperatureLower = 1.0 / std::max(1e-10, column.gasTemperature[i + 1]);
      double invTotalConcentrationLower = 1.0 / std::max(1e-10, column.totalConcentration[i]);

      upper[i] = -timeStep * (D * idx2 * invGasTemperatureUpper * column.facePressures[i] * invTotalConcentrationUpper);

      lower[i] = -timeStep * (column.interstitialGasVelocity[i] * idx +
                              D * idx2 * invGasTemperatureLower * column.facePressures[i] * invTotalConcentrationLower);
    }

    for (size_t i = 1; i < Ngrid; ++i)
    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[i]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[i]);

      diag[i] = 1.0 + timeStep * (column.interstitialGasVelocity[i] * idx +
                                  D * idx2 * invGasTemperature *
                                      (column.facePressures[i - 1] + column.facePressures[i]) * invTotalConcentration);

      rhs[i] = column.concentration[i * Ncomp + comp] +
               timeStep * (-b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]));
    }

    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[Ngrid]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[Ngrid]);

      diag[Ngrid] =
          1.0 + timeStep * (column.interstitialGasVelocity[Ngrid] * idx +
                            D * idx2 * invGasTemperature * column.facePressures[Ngrid - 1] * invTotalConcentration);

      rhs[Ngrid] =
          column.concentration[Ngrid * Ncomp + comp] +
          timeStep *
              (-b * (column.equilibriumAdsorption[Ngrid * Ncomp + comp] - column.adsorption[Ngrid * Ncomp + comp]));
    }

    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.concentration[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgtsv failed in computeConcentrationUpdateMatrixEnergyBalance");
    }

    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

void computeConcentrationUpdateMatrixFinal(Column& column, double timeStep, std::vector<double>& solved)
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
    info = 0;
    std::fill(ab.begin(), ab.end(), 0.0);

    double D = column.components[comp].D;
    double b = column.prefactorMassTransfer[comp];

    for (size_t i = 0; i < Ngrid; ++i)
    {
      lower[i] = column.interstitialGasVelocity[i] * idx + D * idx2;
      upper[i] = D * idx2;
    }

    upper[0] = 0.0;
    diag[0] = 0.0;
    preRHS[0] = 0.0;

    for (size_t i = 1; i < Ngrid; ++i)
    {
      diag[i] = -(column.interstitialGasVelocity[i] * idx + 2.0 * D * idx2);
      preRHS[i] = -b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]);
    }

    diag[Ngrid] = -(column.interstitialGasVelocity[Ngrid] * idx + D * idx2);
    preRHS[Ngrid] = -b * (column.equilibriumAdsorption[Ngrid * Ncomp + comp] - column.adsorption[Ngrid * Ncomp + comp]);

    for (size_t i = 0; i + 1 < Ngrid; ++i)
    {
      abS[6, i] = dt2 * lower[i + 1] * lower[i];
      abS[2, i + 2] = dt2 * upper[i] * upper[i + 1];
    }

    for (size_t i = 0; i < Ngrid; ++i)
    {
      abS[5, i] = dt2 * lower[i] * (diag[i] + diag[i + 1]);
      abS[3, i + 1] = dt2 * upper[i] * (diag[i] + diag[i + 1]);
    }

    for (size_t i = 1; i < Ngrid; ++i)
    {
      abS[4, i] = 1.0 + dt2 * (lower[i - 1] * upper[i - 1] + diag[i] * diag[i] + lower[i] * upper[i]);
      rhs[i] = column.concentration[i * Ncomp + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }

    abS[4, Ngrid] = 1.0 + dt2 * (lower[Ngrid - 1] * upper[Ngrid - 1] + diag[Ngrid] * diag[Ngrid]);
    rhs[Ngrid] = column.concentration[Ngrid * Ncomp + comp] -
                 dt2 * (lower[Ngrid - 1] * preRHS[Ngrid - 1] + diag[Ngrid] * preRHS[Ngrid]);

    abS[4, 0] = 1.0;
    rhs[0] = column.concentration[comp];

    dgbsv_(&n, &kl, &ku, &nrhs, ab.data(), &ldab, ipiv.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgbsv failed in computeConcentrationUpdateMatrixFinal");
    }

    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}

void computeConcentrationUpdateMatrixFinalEnergyBalance(Column& column, double timeStep, std::vector<double>& solved)
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

  std::vector<double> ab(7 * (Ngrid + 1), 0.0);
  std::mdspan abS(ab.data(), 7, Ngrid + 1);
  std::vector<int> ipiv(Ngrid + 1);

  int n = static_cast<int>(Ngrid + 1);
  int info = 0;
  int ku = 2;
  int kl = 2;
  int ldab = 7;
  int nrhs = 1;

  for (size_t grid = 0; grid < Ngrid; ++grid)
  {
    column.facePressures[grid] = 0.5 * (column.totalPressure[grid] + column.totalPressure[grid + 1]) / R;
  }

  for (size_t comp = 0; comp < Ncomp; ++comp)
  {
    info = 0;
    std::fill(ab.begin(), ab.end(), 0.0);

    double D = column.components[comp].D;
    double b = column.prefactorMassTransfer[comp];

    for (size_t i = 0; i < Ngrid; ++i)
    {
      double invGasTemperatureUpper = 1.0 / std::max(1e-10, column.gasTemperature[i]);
      double invTotalConcentrationUpper = 1.0 / std::max(1e-10, column.totalConcentration[i + 1]);
      double invGasTemperatureLower = 1.0 / std::max(1e-10, column.gasTemperature[i + 1]);
      double invTotalConcentrationLower = 1.0 / std::max(1e-10, column.totalConcentration[i]);

      upper[i] = D * idx2 * invGasTemperatureUpper * column.facePressures[i] * invTotalConcentrationUpper;

      lower[i] = column.interstitialGasVelocity[i] * idx +
                 D * idx2 * invGasTemperatureLower * column.facePressures[i] * invTotalConcentrationLower;
    }

    upper[0] = 0.0;
    diag[0] = 0.0;
    preRHS[0] = 0.0;

    for (size_t i = 1; i < Ngrid; ++i)
    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[i]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[i]);

      diag[i] = -(column.interstitialGasVelocity[i] * idx +
                  D * idx2 * invGasTemperature * (column.facePressures[i - 1] + column.facePressures[i]) *
                      invTotalConcentration);

      preRHS[i] = -b * (column.equilibriumAdsorption[i * Ncomp + comp] - column.adsorption[i * Ncomp + comp]);
    }

    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[Ngrid]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[Ngrid]);

      diag[Ngrid] = -(column.interstitialGasVelocity[Ngrid] * idx +
                      D * idx2 * invGasTemperature * column.facePressures[Ngrid - 1] * invTotalConcentration);

      preRHS[Ngrid] =
          -b * (column.equilibriumAdsorption[Ngrid * Ncomp + comp] - column.adsorption[Ngrid * Ncomp + comp]);
    }

    for (size_t i = 0; i + 1 < Ngrid; ++i)
    {
      abS[6, i] = dt2 * lower[i + 1] * lower[i];
      abS[2, i + 2] = dt2 * upper[i] * upper[i + 1];
    }

    for (size_t i = 0; i < Ngrid; ++i)
    {
      abS[5, i] = dt2 * lower[i] * (diag[i] + diag[i + 1]);
      abS[3, i + 1] = dt2 * upper[i] * (diag[i] + diag[i + 1]);
    }

    for (size_t i = 1; i < Ngrid; ++i)
    {
      abS[4, i] = 1.0 + dt2 * (lower[i - 1] * upper[i - 1] + diag[i] * diag[i] + lower[i] * upper[i]);
      rhs[i] = column.concentration[i * Ncomp + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }

    abS[4, Ngrid] = 1.0 + dt2 * (lower[Ngrid - 1] * upper[Ngrid - 1] + diag[Ngrid] * diag[Ngrid]);
    rhs[Ngrid] = column.concentration[Ngrid * Ncomp + comp] -
                 dt2 * (lower[Ngrid - 1] * preRHS[Ngrid - 1] + diag[Ngrid] * preRHS[Ngrid]);

    abS[4, 0] = 1.0;
    rhs[0] = column.concentration[comp];

    dgbsv_(&n, &kl, &ku, &nrhs, ab.data(), &ldab, ipiv.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgbsv failed in computeConcentrationUpdateMatrixFinalEnergyBalance");
    }

    for (size_t grid = 0; grid < Ngrid + 1; ++grid)
    {
      solved[grid * Ncomp + comp] = rhs[grid];
    }
  }
}
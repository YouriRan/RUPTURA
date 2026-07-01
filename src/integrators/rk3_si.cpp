#include "rk3_si.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mdspan>
#include <print>
#include <stdexcept>
#include <vector>

#include "compute.h"
#include "utils.h"

extern "C"
{
  void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
  void dgbsv_(int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
}

bool SemiImplicitRungeKutta3::propagate(Column& column, size_t step, Timing& timings)
{
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;
  Column newColumn(column);

  std::vector<double> solved(column.moleFraction.size());

  std::vector<double> implicitInvKLs(numberOfComponents);
  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    implicitInvKLs[comp] = 1.0 / (1.0 + timeStep * column.components[comp].massTransferCoefficient);
  }

  if (autoNumberOfSteps)
  {
    double tolerance = 0.0;
    for (size_t j = 0; j < numberOfComponents; ++j)
    {
      tolerance = std::max(tolerance, std::abs((column.moleFraction[numberOfGridPoints * numberOfComponents + j] /
                                                column.components[j].initialGasMoleFraction) -
                                               1.0));
    }

    if (tolerance < 0.01)
    {
      std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoNumberOfSteps = false;
    }
  }

  // SSP-RK Step 1
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;
    newColumn.adsorption[i] =
        (column.adsorption[i] + timeStep * kl * column.equilibriumAdsorption[i]) * implicitInvKLs[comp];
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&] { computeConcentrationUpdateMatrixEnergyBalance(newColumn, timeStep, solved); });
  }
  else
  {
    timings.measure(timings.computeDerivatives, [&] { computeConcentrationUpdateMatrix(newColumn, timeStep, solved); });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.moleFraction[i] = solved[i];
  }

  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });
  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });
  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  // SSP-RK Step 2
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;
    newColumn.adsorption[i] =
        0.75 * column.adsorption[i] +
        0.25 * (newColumn.adsorption[i] + timeStep * kl * newColumn.equilibriumAdsorption[i]) * implicitInvKLs[comp];
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&] { computeConcentrationUpdateMatrixEnergyBalance(newColumn, timeStep, solved); });
  }
  else
  {
    timings.measure(timings.computeDerivatives, [&] { computeConcentrationUpdateMatrix(newColumn, timeStep, solved); });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.moleFraction[i] = 0.75 * column.moleFraction[i] + 0.25 * solved[i];
  }

  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });
  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });
  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  // SSP-RK Step 3
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;
    newColumn.adsorption[i] = (1.0 / 3.0) * column.adsorption[i] +
                              (2.0 / 3.0) *
                                  (newColumn.adsorption[i] + timeStep * kl * newColumn.equilibriumAdsorption[i]) *
                                  implicitInvKLs[comp];
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&] { computeConcentrationUpdateMatrixEnergyBalance(newColumn, timeStep, solved); });
  }
  else
  {
    timings.measure(timings.computeDerivatives, [&] { computeConcentrationUpdateMatrix(newColumn, timeStep, solved); });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.moleFraction[i] = (1.0 / 3.0) * column.moleFraction[i] + (2.0 / 3.0) * solved[i];
  }

  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });
  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });
  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  // final implicit adsorption update
  for (size_t i = 0; i < column.adsorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;

    newColumn.adsorption[i] =
        (newColumn.adsorption[i] + timeStep * timeStep * kl * kl * newColumn.equilibriumAdsorption[i]) /
        (1.0 + timeStep * timeStep * kl * kl);
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&] { computeConcentrationUpdateMatrixEnergyBalanceFinal(newColumn, timeStep, solved); });
  }
  else
  {
    timings.measure(timings.computeDerivatives,
                    [&] { computeConcentrationUpdateMatrixFinal(newColumn, timeStep, solved); });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.moleFraction[i] = solved[i];
  }

  timings.measure(timings.computePressure, [&] { computePressure(newColumn); });
  timings.measure(timings.computeEquilibriumLoadings, [&] { computeEquilibriumLoadings(newColumn); });
  timings.measure(timings.computeVelocity, [&] { computeVelocity(newColumn); });

  column = newColumn;
  enforceBoundaryCondition(column);
  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}

void computeConcentrationUpdateMatrix(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  int n = static_cast<int>(numberOfGridPoints + 1);
  int nrhs = 1;

  std::vector<double> upper(numberOfGridPoints);
  std::vector<double> lower(numberOfGridPoints);
  std::vector<double> diag(numberOfGridPoints + 1);
  std::vector<double> rhs(numberOfGridPoints + 1);

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    int info = 0;

    double axialDispersionCoefficient = column.components[comp].axialDispersionCoefficient;
    double b = column.prefactorMassTransfer[comp];

    std::fill(lower.begin(), lower.end(), 0.0);
    std::fill(upper.begin(), upper.end(), 0.0);

    for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
    {
      const double invCT = 1.0 / std::max(1e-10, column.totalConcentration[i]);
      const double dctDz = i == 0 ? 0.0 : (column.totalConcentration[i] - column.totalConcentration[i - 1]) * idx;
      const double d2ctDz2 =
          i == 0 ? 0.0
          : i < numberOfGridPoints
              ? (column.totalConcentration[i + 1] - 2.0 * column.totalConcentration[i] +
                 column.totalConcentration[i - 1]) *
                    idx2
              : (column.totalConcentration[numberOfGridPoints - 1] - column.totalConcentration[numberOfGridPoints]) *
                    idx2;
      const double lowerOp = i == 0 ? 0.0
                                    : column.interstitialGasVelocity[i] * idx +
                                          axialDispersionCoefficient * (idx2 - 2.0 * dctDz * invCT * idx);
      const double upperOp = (i > 0 && i < numberOfGridPoints) ? axialDispersionCoefficient * idx2 : 0.0;
      const double diagOp = -column.interstitialGasVelocity[i] * idx +
                            axialDispersionCoefficient * (2.0 * dctDz * invCT * idx - 2.0 * idx2 + d2ctDz2 * invCT);
      const double adsorptionSource = -b * (column.equilibriumAdsorption[i * numberOfComponents + comp] -
                                            column.adsorption[i * numberOfComponents + comp]);
      diag[i] = 1.0 - timeStep * diagOp;
      rhs[i] = column.moleFraction[i * numberOfComponents + comp] + timeStep * adsorptionSource * invCT;
      if (i > 0) lower[i - 1] = -timeStep * lowerOp;
      if (i < numberOfGridPoints) upper[i] = -timeStep * upperOp;
    }

    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.moleFraction[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgtsv failed in computeConcentrationUpdateMatrix");
    }

    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      solved[grid * numberOfComponents + comp] = rhs[grid];
    }
  }
}

void computeConcentrationUpdateMatrixEnergyBalance(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  int n = static_cast<int>(numberOfGridPoints + 1);
  int nrhs = 1;

  std::vector<double> upper(numberOfGridPoints);
  std::vector<double> lower(numberOfGridPoints);
  std::vector<double> diag(numberOfGridPoints + 1);
  std::vector<double> rhs(numberOfGridPoints + 1);

  for (size_t grid = 0; grid < numberOfGridPoints; ++grid)
  {
    column.facePressures[grid] = 0.5 * (column.totalPressure[grid] + column.totalPressure[grid + 1]) / R;
  }

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    int info = 0;

    double axialDispersionCoefficient = column.components[comp].axialDispersionCoefficient;
    double b = column.prefactorMassTransfer[comp];

    std::fill(lower.begin(), lower.end(), 0.0);
    std::fill(upper.begin(), upper.end(), 0.0);

    for (size_t i = 1; i < numberOfGridPoints; ++i)
    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[i]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[i]);
      double dctDz = (column.totalConcentration[i] - column.totalConcentration[i - 1]) * idx;
      double dTgDz = (column.gasTemperature[i] - column.gasTemperature[i - 1]) * idx;
      double gamma = invTotalConcentration * dctDz + invGasTemperature * dTgDz;

      double lowerOp = column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * (idx2 - gamma * idx);
      double diagOp =
          -column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * (gamma * idx - 2.0 * idx2);
      double upperOp = axialDispersionCoefficient * idx2;
      double adsorptionSource = -b * (column.equilibriumAdsorption[i * numberOfComponents + comp] -
                                      column.adsorption[i * numberOfComponents + comp]);

      lower[i - 1] = -timeStep * lowerOp;
      diag[i] = 1.0 - timeStep * diagOp;
      upper[i] = -timeStep * upperOp;
      rhs[i] = column.moleFraction[i * numberOfComponents + comp] + timeStep * adsorptionSource * invTotalConcentration;
    }

    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[numberOfGridPoints]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[numberOfGridPoints]);
      double dctDz =
          (column.totalConcentration[numberOfGridPoints] - column.totalConcentration[numberOfGridPoints - 1]) * idx;
      double dTgDz = (column.gasTemperature[numberOfGridPoints] - column.gasTemperature[numberOfGridPoints - 1]) * idx;
      double gamma = invTotalConcentration * dctDz + invGasTemperature * dTgDz;
      double lowerOp =
          column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * (idx2 - gamma * idx);
      double diagOp =
          -column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * (gamma * idx - idx2);
      double adsorptionSource = -b * (column.equilibriumAdsorption[numberOfGridPoints * numberOfComponents + comp] -
                                      column.adsorption[numberOfGridPoints * numberOfComponents + comp]);

      lower[numberOfGridPoints - 1] = -timeStep * lowerOp;
      diag[numberOfGridPoints] = 1.0 - timeStep * diagOp;
      rhs[numberOfGridPoints] = column.moleFraction[numberOfGridPoints * numberOfComponents + comp] +
                                timeStep * adsorptionSource * invTotalConcentration;
    }

    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.moleFraction[comp];

    dgtsv_(&n, &nrhs, lower.data(), diag.data(), upper.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgtsv failed in computeConcentrationUpdateMatrixEnergyBalance");
    }

    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      solved[grid * numberOfComponents + comp] = rhs[grid];
    }
  }
}

void computeConcentrationUpdateMatrixFinal(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  double dt2 = timeStep * timeStep;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  std::vector<double> upper(numberOfGridPoints);
  std::vector<double> lower(numberOfGridPoints);
  std::vector<double> diag(numberOfGridPoints + 1);
  std::vector<double> preRHS(numberOfGridPoints + 1);
  std::vector<double> rhs(numberOfGridPoints + 1);

  std::vector<double> ab(7 * (numberOfGridPoints + 1), 0.0);
  std::mdspan<double, std::dextents<size_t, 2>, std::layout_left> abS(ab.data(), 7, numberOfGridPoints + 1);
  std::vector<int> ipiv(numberOfGridPoints + 1);

  int n = static_cast<int>(numberOfGridPoints + 1);
  int info = 0;
  int ku = 2;
  int kl = 2;
  int ldab = 7;
  int nrhs = 1;

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    info = 0;
    std::fill(ab.begin(), ab.end(), 0.0);

    double axialDispersionCoefficient = column.components[comp].axialDispersionCoefficient;
    double b = column.prefactorMassTransfer[comp];

    for (size_t i = 0; i < numberOfGridPoints; ++i)
    {
      lower[i] = 0.0;
      upper[i] = 0.0;
    }

    upper[0] = 0.0;
    diag[0] = 0.0;
    preRHS[0] = 0.0;

    for (size_t i = 1; i < numberOfGridPoints; ++i)
    {
      const double invCT = 1.0 / std::max(1e-10, column.totalConcentration[i]);
      const double dctDz = (column.totalConcentration[i] - column.totalConcentration[i - 1]) * idx;
      const double d2ctDz2 =
          (column.totalConcentration[i + 1] - 2.0 * column.totalConcentration[i] + column.totalConcentration[i - 1]) *
          idx2;
      lower[i - 1] =
          column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * (idx2 - 2.0 * dctDz * invCT * idx);
      upper[i] = axialDispersionCoefficient * idx2;
      diag[i] = -column.interstitialGasVelocity[i] * idx +
                axialDispersionCoefficient * (2.0 * dctDz * invCT * idx - 2.0 * idx2 + d2ctDz2 * invCT);
      preRHS[i] = -b *
                  (column.equilibriumAdsorption[i * numberOfComponents + comp] -
                   column.adsorption[i * numberOfComponents + comp]) *
                  invCT;
    }

    {
      const double invCT = 1.0 / std::max(1e-10, column.totalConcentration[numberOfGridPoints]);
      const double dctDz =
          (column.totalConcentration[numberOfGridPoints] - column.totalConcentration[numberOfGridPoints - 1]) * idx;
      const double d2ctDz2 =
          (column.totalConcentration[numberOfGridPoints - 1] - column.totalConcentration[numberOfGridPoints]) * idx2;
      lower[numberOfGridPoints - 1] = column.interstitialGasVelocity[numberOfGridPoints] * idx +
                                      axialDispersionCoefficient * (idx2 - 2.0 * dctDz * invCT * idx);
      diag[numberOfGridPoints] = -column.interstitialGasVelocity[numberOfGridPoints] * idx +
                                 axialDispersionCoefficient * (2.0 * dctDz * invCT * idx - idx2 + d2ctDz2 * invCT);
      preRHS[numberOfGridPoints] = -b *
                                   (column.equilibriumAdsorption[numberOfGridPoints * numberOfComponents + comp] -
                                    column.adsorption[numberOfGridPoints * numberOfComponents + comp]) *
                                   invCT;
    }

    for (size_t i = 0; i + 1 < numberOfGridPoints; ++i)
    {
      abS[6, i] = dt2 * lower[i + 1] * lower[i];
      abS[2, i + 2] = dt2 * upper[i] * upper[i + 1];
    }

    for (size_t i = 0; i < numberOfGridPoints; ++i)
    {
      abS[5, i] = dt2 * lower[i] * (diag[i] + diag[i + 1]);
      abS[3, i + 1] = dt2 * upper[i] * (diag[i] + diag[i + 1]);
    }

    for (size_t i = 1; i < numberOfGridPoints; ++i)
    {
      abS[4, i] = 1.0 + dt2 * (lower[i - 1] * upper[i - 1] + diag[i] * diag[i] + lower[i] * upper[i]);
      rhs[i] = column.moleFraction[i * numberOfComponents + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }

    abS[4, numberOfGridPoints] = 1.0 + dt2 * (lower[numberOfGridPoints - 1] * upper[numberOfGridPoints - 1] +
                                              diag[numberOfGridPoints] * diag[numberOfGridPoints]);
    rhs[numberOfGridPoints] = column.moleFraction[numberOfGridPoints * numberOfComponents + comp] -
                              dt2 * (lower[numberOfGridPoints - 1] * preRHS[numberOfGridPoints - 1] +
                                     diag[numberOfGridPoints] * preRHS[numberOfGridPoints]);

    abS[4, 0] = 1.0;
    rhs[0] = column.moleFraction[comp];

    dgbsv_(&n, &kl, &ku, &nrhs, ab.data(), &ldab, ipiv.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgbsv failed in computeConcentrationUpdateMatrixFinal");
    }

    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      solved[grid * numberOfComponents + comp] = rhs[grid];
    }
  }
}

void computeConcentrationUpdateMatrixEnergyBalanceFinal(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  double dt2 = timeStep * timeStep;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;

  std::vector<double> upper(numberOfGridPoints);
  std::vector<double> lower(numberOfGridPoints);
  std::vector<double> diag(numberOfGridPoints + 1);
  std::vector<double> preRHS(numberOfGridPoints + 1);
  std::vector<double> rhs(numberOfGridPoints + 1);

  std::vector<double> ab(7 * (numberOfGridPoints + 1), 0.0);
  std::mdspan abS(ab.data(), 7, numberOfGridPoints + 1);
  std::vector<int> ipiv(numberOfGridPoints + 1);

  int n = static_cast<int>(numberOfGridPoints + 1);
  int info = 0;
  int ku = 2;
  int kl = 2;
  int ldab = 7;
  int nrhs = 1;

  for (size_t grid = 0; grid < numberOfGridPoints; ++grid)
  {
    column.facePressures[grid] = 0.5 * (column.totalPressure[grid] + column.totalPressure[grid + 1]) / R;
  }

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    info = 0;
    std::fill(ab.begin(), ab.end(), 0.0);

    double axialDispersionCoefficient = column.components[comp].axialDispersionCoefficient;
    double b = column.prefactorMassTransfer[comp];

    for (size_t i = 0; i < numberOfGridPoints; ++i)
    {
      lower[i] = 0.0;
      upper[i] = 0.0;
    }

    upper[0] = 0.0;
    diag[0] = 0.0;
    preRHS[0] = 0.0;

    for (size_t i = 1; i < numberOfGridPoints; ++i)
    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[i]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[i]);
      double dctDz = (column.totalConcentration[i] - column.totalConcentration[i - 1]) * idx;
      double dTgDz = (column.gasTemperature[i] - column.gasTemperature[i - 1]) * idx;
      double gamma = invTotalConcentration * dctDz + invGasTemperature * dTgDz;

      lower[i - 1] = column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * (idx2 - gamma * idx);
      diag[i] = -column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * (gamma * idx - 2.0 * idx2);
      upper[i] = axialDispersionCoefficient * idx2;

      preRHS[i] = -b *
                  (column.equilibriumAdsorption[i * numberOfComponents + comp] -
                   column.adsorption[i * numberOfComponents + comp]) *
                  invTotalConcentration;
    }

    {
      double invGasTemperature = 1.0 / std::max(1e-10, column.gasTemperature[numberOfGridPoints]);
      double invTotalConcentration = 1.0 / std::max(1e-10, column.totalConcentration[numberOfGridPoints]);
      double dctDz =
          (column.totalConcentration[numberOfGridPoints] - column.totalConcentration[numberOfGridPoints - 1]) * idx;
      double dTgDz = (column.gasTemperature[numberOfGridPoints] - column.gasTemperature[numberOfGridPoints - 1]) * idx;
      double gamma = invTotalConcentration * dctDz + invGasTemperature * dTgDz;

      lower[numberOfGridPoints - 1] =
          column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * (idx2 - gamma * idx);
      diag[numberOfGridPoints] =
          -column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * (gamma * idx - idx2);

      preRHS[numberOfGridPoints] = -b *
                                   (column.equilibriumAdsorption[numberOfGridPoints * numberOfComponents + comp] -
                                    column.adsorption[numberOfGridPoints * numberOfComponents + comp]) *
                                   invTotalConcentration;
    }

    for (size_t i = 0; i + 1 < numberOfGridPoints; ++i)
    {
      abS[6, i] = dt2 * lower[i + 1] * lower[i];
      abS[2, i + 2] = dt2 * upper[i] * upper[i + 1];
    }

    for (size_t i = 0; i < numberOfGridPoints; ++i)
    {
      abS[5, i] = dt2 * lower[i] * (diag[i] + diag[i + 1]);
      abS[3, i + 1] = dt2 * upper[i] * (diag[i] + diag[i + 1]);
    }

    for (size_t i = 1; i < numberOfGridPoints; ++i)
    {
      abS[4, i] = 1.0 + dt2 * (lower[i - 1] * upper[i - 1] + diag[i] * diag[i] + lower[i] * upper[i]);
      rhs[i] = column.moleFraction[i * numberOfComponents + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }

    abS[4, numberOfGridPoints] = 1.0 + dt2 * (lower[numberOfGridPoints - 1] * upper[numberOfGridPoints - 1] +
                                              diag[numberOfGridPoints] * diag[numberOfGridPoints]);
    rhs[numberOfGridPoints] = column.moleFraction[numberOfGridPoints * numberOfComponents + comp] -
                              dt2 * (lower[numberOfGridPoints - 1] * preRHS[numberOfGridPoints - 1] +
                                     diag[numberOfGridPoints] * preRHS[numberOfGridPoints]);

    abS[4, 0] = 1.0;
    rhs[0] = column.moleFraction[comp];

    dgbsv_(&n, &kl, &ku, &nrhs, ab.data(), &ldab, ipiv.data(), rhs.data(), &n, &info);
    if (info != 0)
    {
      throw std::runtime_error("dgbsv failed in computeConcentrationUpdateMatrixFinalEnergyBalance");
    }

    for (size_t grid = 0; grid < numberOfGridPoints + 1; ++grid)
    {
      solved[grid * numberOfComponents + comp] = rhs[grid];
    }
  }
}

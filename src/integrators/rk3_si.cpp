#include "rk3_si.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <mdspan>
#include <print>
#include <stdexcept>
#include <vector>

#include "compute.h"
#include "rk3.h"
#include "sorption.h"
#include "transport.h"
#include "utils.h"

extern "C"
{
  void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
  void dgbsv_(int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
}

namespace
{
double sorptionDotAt(const Column& column, size_t index)
{
  double sorptionDot = column.physisorptionDot[index];
  const size_t componentBlockSize = column.physisorptionDot.size();
  for (size_t site = 0; site < column.maxChemisorptionSites; ++site)
  {
    sorptionDot += column.chemisorptionDot[site * componentBlockSize + index];
  }
  return sorptionDot;
}
}  // namespace

bool SemiImplicitRungeKutta3::propagate(Column& column, size_t step, Timing& timings)
{
  if (column.surfacePoreTransportEnabled)
  {
    throw std::runtime_error(
        "SemiImplicitRungeKutta3 does not support General chemisorption with surface-pore transport; use RungeKutta3 "
        "or CVODE");
  }

  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;
  const double startTime = static_cast<double>(step) * timeStep;
  Column newColumn(column);

  std::vector<double> solved(column.concentration.size());

  auto finalizeStage = [&](Column& stage)
  {
    timings.measure(timings.updateVelocityAndPressure,
                    [&]
                    {
                      computeBulkSpeciesSink(stage.components, stage.numberOfGridPoints, stage.numberOfComponents,
                                             stage.maxChemisorptionSites, stage.geometry, stage.particleDensity,
                                             stage.concentration, stage.physisorptionDot, stage.chemisorptionDot,
                                             stage.surfaceConcentration, stage.bulkSpeciesSink,
                                             stage.reactionPhysisorptionSource, stage.reactionChemisorptionSource);
                      updateVelocityAndPressure(
                          stage.components, stage.boundaryCondition, stage.numberOfGridPoints, stage.numberOfComponents,
                          stage.inletPressure, stage.outletPressure, stage.pressureGradient, stage.columnLength,
                          stage.geometry, stage.columnEntranceVelocity, stage.dynamicViscosity, stage.resolution,
                          stage.interstitialGasVelocity, stage.gasDensity, stage.totalConcentration,
                          stage.totalPressure, stage.concentration, stage.partialPressure, stage.moleFraction,
                          stage.bulkSpeciesSink, stage.gasTemperature, stage.fluidPhase, stage.liquidDensity,
                          stage.pHMode, stage.pHValue, stage.pKw, stage.pHComponent, stage.pH);
                    });
    timings.measure(timings.computeEquilibriumLoadings,
                    [&]
                    {
                      computePhysisorptionEquilibriumLoadings(
                          stage.physisorptionMixture, stage.numberOfGridPoints, stage.numberOfComponents,
                          stage.maxIsothermTerms, stage.iastPerformance, stage.idealGasMolFractions,
                          stage.adsorbedMolFractions, stage.numberOfMolecules, stage.totalPressure,
                          stage.equilibriumPhysisorption, stage.cachedPressure, stage.cachedGrandPotential,
                          stage.moleFraction, stage.gasTemperature,
                          stage.fluidPhase == Column::FluidPhase::Gas
                              ? MixturePrediction::DrivingForceInput::MoleFraction
                              : MixturePrediction::DrivingForceInput::Concentration,
                          stage.concentration, stage.pH);
                      computeChemisorptionEquilibriumLoadings(
                          stage.chemisorptionMixture, stage.numberOfGridPoints, stage.numberOfComponents,
                          stage.maxChemisorptionSites, stage.iastPerformance, stage.idealGasMolFractions,
                          stage.adsorbedMolFractions, stage.numberOfMolecules, stage.totalPressure,
                          stage.equilibriumChemisorption, stage.cachedChemisorptionPressure,
                          stage.cachedChemisorptionGrandPotential, stage.moleFraction, stage.gasTemperature,
                          stage.fluidPhase == Column::FluidPhase::Gas
                              ? MixturePrediction::DrivingForceInput::MoleFraction
                              : MixturePrediction::DrivingForceInput::Concentration,
                          stage.concentration, stage.pH);
                    });
  };

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
      const bool gasPhase = column.fluidPhase == Column::FluidPhase::Gas;
      const double feed = gasPhase ? column.components[j].initialGasMoleFraction
                                   : column.components[j].inletLiquidConcentration;
      if (feed <= 0.0) continue;

      const size_t outlet = numberOfGridPoints * numberOfComponents + j;
      const double outletValue = gasPhase ? column.moleFraction[outlet] : column.concentration[outlet];
      tolerance = std::max(tolerance, std::abs((outletValue / feed) - 1.0));
    }

    if (tolerance < 0.01)
    {
      std::print("\nConvergence criteria reached, running 10% longer\n\n\n");
      numberOfSteps = static_cast<size_t>(1.1 * static_cast<double>(step));
      autoNumberOfSteps = false;
    }
  }

  // SSP-RK Step 1
  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;
    newColumn.physisorption[i] =
        (column.physisorption[i] + timeStep * kl * column.equilibriumPhysisorption[i]) * implicitInvKLs[comp];
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + timeStep);
                      computeConcentrationUpdateMatrixEnergyBalance(newColumn, timeStep, solved);
                    });
  }
  else
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + timeStep);
                      computeConcentrationUpdateMatrix(newColumn, timeStep, solved);
                    });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.concentration[i] = solved[i];
  }

  finalizeStage(newColumn);

  // SSP-RK Step 2
  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;
    newColumn.physisorption[i] =
        0.75 * column.physisorption[i] +
        0.25 * (newColumn.physisorption[i] + timeStep * kl * newColumn.equilibriumPhysisorption[i]) *
            implicitInvKLs[comp];
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + 0.5 * timeStep);
                      computeConcentrationUpdateMatrixEnergyBalance(newColumn, timeStep, solved);
                    });
  }
  else
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + 0.5 * timeStep);
                      computeConcentrationUpdateMatrix(newColumn, timeStep, solved);
                    });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.concentration[i] = 0.75 * column.concentration[i] + 0.25 * solved[i];
  }

  finalizeStage(newColumn);

  // SSP-RK Step 3
  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;
    newColumn.physisorption[i] =
        (1.0 / 3.0) * column.physisorption[i] +
        (2.0 / 3.0) * (newColumn.physisorption[i] + timeStep * kl * newColumn.equilibriumPhysisorption[i]) *
            implicitInvKLs[comp];
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + timeStep);
                      computeConcentrationUpdateMatrixEnergyBalance(newColumn, timeStep, solved);
                    });
  }
  else
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + timeStep);
                      computeConcentrationUpdateMatrix(newColumn, timeStep, solved);
                    });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.concentration[i] = (1.0 / 3.0) * column.concentration[i] + (2.0 / 3.0) * solved[i];
  }

  finalizeStage(newColumn);

  // final implicit physisorption update
  for (size_t i = 0; i < column.physisorption.size(); ++i)
  {
    size_t comp = i % numberOfComponents;
    double kl = column.components[comp].massTransferCoefficient;

    newColumn.physisorption[i] =
        (newColumn.physisorption[i] + timeStep * timeStep * kl * kl * newColumn.equilibriumPhysisorption[i]) /
        (1.0 + timeStep * timeStep * kl * kl);
  }

  if (newColumn.energyBalance)
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + timeStep);
                      computeConcentrationUpdateMatrixEnergyBalanceFinal(newColumn, timeStep, solved);
                    });
  }
  else
  {
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(newColumn, startTime + timeStep);
                      computeConcentrationUpdateMatrixFinal(newColumn, timeStep, solved);
                    });
  }

  for (size_t i = 0; i < solved.size(); ++i)
  {
    newColumn.concentration[i] = solved[i];
  }
  clampNonnegative(newColumn.state);

  finalizeStage(newColumn);

  column = newColumn;
  return (!autoNumberOfSteps && step >= numberOfSteps - 1);
}

void computeConcentrationUpdateMatrix(Column& column, double timeStep, std::vector<double>& solved)
{
  double idx = 1.0 / column.resolution;
  double idx2 = idx * idx;
  size_t numberOfGridPoints = column.numberOfGridPoints;
  size_t numberOfComponents = column.numberOfComponents;
  const double adsorptionPrefactor = column.geometry.loadingPrefactor(column.particleDensity);

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

    std::fill(lower.begin(), lower.end(), 0.0);
    std::fill(upper.begin(), upper.end(), 0.0);

    for (size_t i = 0; i < numberOfGridPoints + 1; ++i)
    {
      const double lowerOp = i == 0 ? 0.0 : column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * idx2;
      const double upperOp = (i > 0 && i < numberOfGridPoints) ? axialDispersionCoefficient * idx2 : 0.0;
      const double diagOp = -column.interstitialGasVelocity[i] * idx -
                            axialDispersionCoefficient * (i < numberOfGridPoints ? 2.0 : 1.0) * idx2;
      const double adsorptionSource = -adsorptionPrefactor * sorptionDotAt(column, i * numberOfComponents + comp);
      diag[i] = 1.0 - timeStep * diagOp;
      rhs[i] = column.concentration[i * numberOfComponents + comp] + timeStep * adsorptionSource;
      if (i > 0) lower[i - 1] = -timeStep * lowerOp;
      if (i < numberOfGridPoints) upper[i] = -timeStep * upperOp;
    }

    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.concentration[comp];

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
  const double adsorptionPrefactor = column.geometry.loadingPrefactor(column.particleDensity);

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

    std::fill(lower.begin(), lower.end(), 0.0);
    std::fill(upper.begin(), upper.end(), 0.0);

    for (size_t i = 1; i < numberOfGridPoints; ++i)
    {
      double lowerOp = column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * idx2;
      double diagOp = -column.interstitialGasVelocity[i] * idx - 2.0 * axialDispersionCoefficient * idx2;
      double upperOp = axialDispersionCoefficient * idx2;
      double adsorptionSource = -adsorptionPrefactor * sorptionDotAt(column, i * numberOfComponents + comp);

      lower[i - 1] = -timeStep * lowerOp;
      diag[i] = 1.0 - timeStep * diagOp;
      upper[i] = -timeStep * upperOp;
      rhs[i] = column.concentration[i * numberOfComponents + comp] + timeStep * adsorptionSource;
    }

    {
      double lowerOp = column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * idx2;
      double diagOp = -column.interstitialGasVelocity[numberOfGridPoints] * idx - axialDispersionCoefficient * idx2;
      double adsorptionSource =
          -adsorptionPrefactor * sorptionDotAt(column, numberOfGridPoints * numberOfComponents + comp);

      lower[numberOfGridPoints - 1] = -timeStep * lowerOp;
      diag[numberOfGridPoints] = 1.0 - timeStep * diagOp;
      rhs[numberOfGridPoints] =
          column.concentration[numberOfGridPoints * numberOfComponents + comp] + timeStep * adsorptionSource;
    }

    diag[0] = 1.0;
    upper[0] = 0.0;
    rhs[0] = column.concentration[comp];

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
  const double adsorptionPrefactor = column.geometry.loadingPrefactor(column.particleDensity);

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
      lower[i - 1] = column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * idx2;
      upper[i] = axialDispersionCoefficient * idx2;
      diag[i] = -column.interstitialGasVelocity[i] * idx - 2.0 * axialDispersionCoefficient * idx2;
      preRHS[i] = -adsorptionPrefactor * sorptionDotAt(column, i * numberOfComponents + comp);
    }

    {
      lower[numberOfGridPoints - 1] =
          column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * idx2;
      diag[numberOfGridPoints] =
          -column.interstitialGasVelocity[numberOfGridPoints] * idx + -axialDispersionCoefficient * idx2;
      preRHS[numberOfGridPoints] =
          -adsorptionPrefactor * sorptionDotAt(column, numberOfGridPoints * numberOfComponents + comp);
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
      rhs[i] = column.concentration[i * numberOfComponents + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }

    abS[4, numberOfGridPoints] = 1.0 + dt2 * (lower[numberOfGridPoints - 1] * upper[numberOfGridPoints - 1] +
                                              diag[numberOfGridPoints] * diag[numberOfGridPoints]);
    rhs[numberOfGridPoints] = column.concentration[numberOfGridPoints * numberOfComponents + comp] -
                              dt2 * (lower[numberOfGridPoints - 1] * preRHS[numberOfGridPoints - 1] +
                                     diag[numberOfGridPoints] * preRHS[numberOfGridPoints]);

    abS[4, 0] = 1.0;
    rhs[0] = column.concentration[comp];

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
  const double adsorptionPrefactor = column.geometry.loadingPrefactor(column.particleDensity);

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
      lower[i - 1] = column.interstitialGasVelocity[i] * idx + axialDispersionCoefficient * idx2;
      diag[i] = -column.interstitialGasVelocity[i] * idx - 2.0 * axialDispersionCoefficient * idx2;
      upper[i] = axialDispersionCoefficient * idx2;

      preRHS[i] = -adsorptionPrefactor * sorptionDotAt(column, i * numberOfComponents + comp);
    }

    {
      lower[numberOfGridPoints - 1] =
          column.interstitialGasVelocity[numberOfGridPoints] * idx + axialDispersionCoefficient * idx2;
      diag[numberOfGridPoints] =
          -column.interstitialGasVelocity[numberOfGridPoints] * idx - axialDispersionCoefficient * idx2;

      preRHS[numberOfGridPoints] =
          -adsorptionPrefactor * sorptionDotAt(column, numberOfGridPoints * numberOfComponents + comp);
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
      rhs[i] = column.concentration[i * numberOfComponents + comp] -
               dt2 * (lower[i - 1] * preRHS[i - 1] + diag[i] * preRHS[i] + upper[i] * preRHS[i + 1]);
    }

    abS[4, numberOfGridPoints] = 1.0 + dt2 * (lower[numberOfGridPoints - 1] * upper[numberOfGridPoints - 1] +
                                              diag[numberOfGridPoints] * diag[numberOfGridPoints]);
    rhs[numberOfGridPoints] = column.concentration[numberOfGridPoints * numberOfComponents + comp] -
                              dt2 * (lower[numberOfGridPoints - 1] * preRHS[numberOfGridPoints - 1] +
                                     diag[numberOfGridPoints] * preRHS[numberOfGridPoints]);

    abS[4, 0] = 1.0;
    rhs[0] = column.concentration[comp];

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

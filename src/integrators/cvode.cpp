#include "cvode.h"

#include <cmath>
#include <iostream>
#include <mdspan>
#include <span>
#include <vector>

// static inline sunindextype idxP(size_t grid, size_t comp, size_t Ncomp)
// {
//   return static_cast<sunindextype>(grid * Ncomp + comp);
// }

// static inline sunindextype idxQ(size_t grid, size_t comp, size_t Ncomp, size_t offset)
// {
//   return static_cast<sunindextype>(offset + grid * Ncomp + comp);
// }

// struct CscPattern
// {
//   std::vector<sunindextype> colptr;  // size N+1
//   std::vector<sunindextype> rowind;  // size nnz
// };

// static void buildPatternCSC(SUNMatrix A, size_t Ngrid, size_t Ncomp)
// {
//   const size_t nNodes = Ngrid + 1, offset = nNodes * Ncomp, N = 2 * offset;

//   auto* colptr = SM_INDEXPTRS_S(A);
//   auto* rowind = SM_INDEXVALS_S(A);

//   sunindextype nnz = 0;
//   for (size_t col = 0; col < N; ++col)
//   {
//     colptr[col] = nnz;

//     const bool isPcol = (col < offset);
//     const size_t local = isPcol ? col : (col - offset);
//     const size_t grid = local / Ncomp;
//     const size_t comp = local % Ncomp;

//     if (isPcol)
//     {
//       // rows: pDot(grid), pDot(grid+1), pDot(grid-1), with exclusions (pDot row 0 is zero)
//       if (grid >= 1) rowind[nnz++] = P(grid, comp, Ncomp);
//       if (grid + 1 <= Ngrid) rowind[nnz++] = P(grid + 1, comp, Ncomp);
//       if (grid >= 2) rowind[nnz++] = P(grid - 1, comp, Ncomp);

//       // sort+unique within this column (KLU likes sorted row indices)
//       auto* beg = rowind + colptr[col];
//       auto* end = rowind + nnz;
//       std::sort(beg, end);
//       end = std::unique(beg, end);
//       nnz = static_cast<sunindextype>(end - rowind);
//     }
//     else
//     {
//       // rows: qDot(grid) always; pDot(grid) if grid>=1
//       rowind[nnz++] = Q(grid, comp, Ncomp, offset);
//       if (grid >= 1) rowind[nnz++] = P(grid, comp, Ncomp);
//       if (grid >= 1 && rowind[nnz - 2] > rowind[nnz - 1]) std::swap(rowind[nnz - 2], rowind[nnz - 1]);
//     }
//   }
//   colptr[N] = nnz;
// }

bool CVODE::propagate(Column& column, size_t step)
{
  double t = static_cast<double>(step) * timeStep;
  double tNext = static_cast<double>(step + 1) * timeStep;
  size_t Ngrid = column.Ngrid;
  size_t Ncomp = column.Ncomp;

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

  sunrealtype tReturn = t;

  // maybe not necessary? Are there external things influencing these values?
  copyFromcolumn(column, u);
  copyFromcolumnDot(column, uDot);

  CVode(cvodeMem, tNext, u, &tReturn, CV_NORMAL);

  copyIntocolumn(column, u);
  copyIntocolumnDot(column, uDot);

  return (!autoSteps && step >= numberOfSteps - 1);
}

void CVODE::initialize(Column& column)
{
  // initialize logger and context
  SUNContext_Create(SUN_COMM_NULL, &sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &sunLogger);
  SUNContext_SetLogger(sunContext, sunLogger);

  // create vector that cvode will operate on
  // adsorption, partialPressure, gasTemperature, solidTemperature, wallTemperature
  const sunindextype totalSize = (column.energyBalance)
                                     ? static_cast<sunindextype>(2 * (column.Ncomp + 3) * (column.Ngrid + 1))
                                     : static_cast<sunindextype>(2 * column.Ncomp * (column.Ngrid + 1));
  u = N_VNew_Serial(totalSize, sunContext);
  copyFromcolumn(column, u);
  uDot = N_VNew_Serial(totalSize, sunContext);
  copyFromcolumnDot(column, uDot);

  cvodeMem = CVodeCreate(CV_BDF, sunContext);
  CVodeSetMaxNumSteps(cvodeMem, 1e6);

  CVodeSetUserData(cvodeMem, &column);

  const sunrealtype t0 = 0.0;
  CVodeInit(cvodeMem, f, t0, u);

  CVodeSStolerances(cvodeMem, relativeTolerance, absoluteTolerance);
  solver = SUNNonlinSol_Newton(u, sunContext);
  CVodeSetNonlinearSolver(cvodeMem, solver);

  A = SUNDenseMatrix(totalSize, totalSize, sunContext);
  // CscPattern pat = buildJacobianPatternCSC(column.Ngrid, column.Ncomp);
  // A = SUNSparseMatrix(totalSize, totalSize, static_cast<sunindextype>(pat.rowind.size()), CSC_MAT, sunContext);
  // auto* colptr = SM_INDEXPTRS_S(A);
  // auto* rowind = SM_INDEXVALS_S(A);
  // std::copy(pat.colptr.begin(), pat.colptr.end(), colptr);
  // std::copy(pat.rowind.begin(), pat.rowind.end(), rowind);

  // linSolver = SUNLinSol_Dense(u, A, sunContext);
  linSolver = SUNLinSol_KLU(u, A, sunContext);
  CVodeSetLinearSolver(cvodeMem, linSolver, A);

  CVodeSetJacFn(cvodeMem, nullptr);
  // CVodeSetJacFn(cvodeMem, JacSparse);
}

inline std::span<double> getPartialPressureSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 0 * small + 0 * big, big};
}

inline std::span<double> getAdsorptionSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 0 * small + 1 * big, big};
}

inline std::span<double> getGasTemperatureSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 0 * small + 2 * big, small};
}

inline std::span<double> getSolidTemperatureSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 1 * small + 1 * big, small};
}

inline std::span<double> getWallTemperatureSpan(N_Vector u, size_t Ngrid, size_t Ncomp)
{
  double* base = static_cast<double*>(N_VGetArrayPointer(u));
  size_t small = Ngrid + 1;
  size_t big = small * Ncomp;
  return {base + 2 * small + 1 * big, small};
}

inline void copyFromcolumn(const Column& column, N_Vector u)
{
  std::copy(column.partialPressure.begin(), column.partialPressure.end(),
            getPartialPressureSpan(u, column.Ngrid, column.Ncomp).begin());
  std::copy(column.adsorption.begin(), column.adsorption.end(),
            getAdsorptionSpan(u, column.Ngrid, column.Ncomp).begin());
  if (column.energyBalance)
  {
    std::copy(column.gasTemperature.begin(), column.gasTemperature.end(),
              getGasTemperatureSpan(u, column.Ngrid, column.Ncomp).begin());
    std::copy(column.solidTemperature.begin(), column.solidTemperature.end(),
              getSolidTemperatureSpan(u, column.Ngrid, column.Ncomp).begin());
    std::copy(column.wallTemperature.begin(), column.wallTemperature.end(),
              getWallTemperatureSpan(u, column.Ngrid, column.Ncomp).begin());
  }
}

inline void copyFromcolumnDot(const Column& column, N_Vector uDot)
{
  std::copy(column.partialPressureDot.begin(), column.partialPressureDot.end(),
            getPartialPressureSpan(uDot, column.Ngrid, column.Ncomp).begin());
  std::copy(column.adsorptionDot.begin(), column.adsorptionDot.end(),
            getAdsorptionSpan(uDot, column.Ngrid, column.Ncomp).begin());
  if (column.energyBalance)
  {
    std::copy(column.gasTemperatureDot.begin(), column.gasTemperatureDot.end(),
              getGasTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin());
    std::copy(column.solidTemperatureDot.begin(), column.solidTemperatureDot.end(),
              getSolidTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin());
    std::copy(column.wallTemperatureDot.begin(), column.wallTemperatureDot.end(),
              getWallTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin());
  }
}

inline void copyIntocolumn(Column& column, N_Vector u)
{
  std::copy(getPartialPressureSpan(u, column.Ngrid, column.Ncomp).begin(),
            getPartialPressureSpan(u, column.Ngrid, column.Ncomp).end(), column.partialPressure.begin());
  std::copy(getAdsorptionSpan(u, column.Ngrid, column.Ncomp).begin(),
            getAdsorptionSpan(u, column.Ngrid, column.Ncomp).end(), column.adsorption.begin());
  if (column.energyBalance)
  {
    std::copy(getGasTemperatureSpan(u, column.Ngrid, column.Ncomp).begin(),
              getGasTemperatureSpan(u, column.Ngrid, column.Ncomp).end(), column.gasTemperature.begin());
    std::copy(getSolidTemperatureSpan(u, column.Ngrid, column.Ncomp).begin(),
              getSolidTemperatureSpan(u, column.Ngrid, column.Ncomp).end(), column.solidTemperature.begin());
    std::copy(getWallTemperatureSpan(u, column.Ngrid, column.Ncomp).begin(),
              getWallTemperatureSpan(u, column.Ngrid, column.Ncomp).end(), column.wallTemperature.begin());
  }
}

inline void copyIntocolumnDot(Column& column, N_Vector uDot)
{
  std::copy(getPartialPressureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
            getPartialPressureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.partialPressureDot.begin());
  std::copy(getAdsorptionSpan(uDot, column.Ngrid, column.Ncomp).begin(),
            getAdsorptionSpan(uDot, column.Ngrid, column.Ncomp).end(), column.adsorptionDot.begin());
  if (column.energyBalance)
  {
    std::copy(getGasTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
              getGasTemperatureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.gasTemperatureDot.begin());
    std::copy(getSolidTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
              getSolidTemperatureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.solidTemperatureDot.begin());
    std::copy(getWallTemperatureSpan(uDot, column.Ngrid, column.Ncomp).begin(),
              getWallTemperatureSpan(uDot, column.Ngrid, column.Ncomp).end(), column.wallTemperatureDot.begin());
  }
}
static int f(sunrealtype t, N_Vector u, N_Vector uDot, void* user_data)
{
  auto* column = reinterpret_cast<Column*>(user_data);

  auto spanPartialPressure = getPartialPressureSpan(u, column->Ngrid, column->Ncomp);
  auto spanAdsorption = getAdsorptionSpan(u, column->Ngrid, column->Ncomp);

  auto spanAdsorptionDot = getAdsorptionSpan(uDot, column->Ngrid, column->Ncomp);
  auto spanPartialPressureDot = getPartialPressureSpan(uDot, column->Ngrid, column->Ncomp);

  computeEquilibriumLoadings(column->mixture, column->Ngrid, column->Ncomp, column->maxIsothermTerms,
                             column->iastPerformance, column->pressureGradient, column->columnLength,
                             column->totalPressure, spanPartialPressure, column->idealGasMolFractions,
                             column->adsorbedMolFractions, column->numberOfMolecules, column->equilibriumAdsorption,
                             column->cachedPressure, column->cachedGrandPotential);

  switch (column->velocityProfile)
  {
    case Column::VelocityProfile::FixedPressureGradient:
    {
      computeVelocityFixedGradient(column->components, column->Ngrid, column->Ncomp, column->pressureGradient,
                                   column->columnEntranceVelocity, column->resolution, column->prefactorMassTransfer,
                                   column->interstitialGasVelocity, column->totalPressure, spanAdsorption,
                                   column->equilibriumAdsorption, spanPartialPressure);
      break;
    }
    case Column::VelocityProfile::Ergun:
    {
      computeVelocityErgun(column->components, column->Ngrid, column->Ncomp, column->externalTemperature,
                           column->voidFraction, column->columnEntranceVelocity, column->dynamicViscosity,
                           column->particleDiameter, column->resolution, column->interstitialGasVelocity,
                           column->totalPressure, spanPartialPressure);
      break;
    }
    case Column::VelocityProfile::FixedVelocity:
      break;
    default:
      break;
  }

  if (column->energyBalance)
  {
    auto spanGasTemperature = getGasTemperatureSpan(u, column->Ngrid, column->Ncomp);
    auto spanSolidTemperature = getSolidTemperatureSpan(u, column->Ngrid, column->Ncomp);
    auto spanWallTemperature = getWallTemperatureSpan(u, column->Ngrid, column->Ncomp);

    auto spanGasTemperatureDot = getGasTemperatureSpan(uDot, column->Ngrid, column->Ncomp);
    auto spanSolidTemperatureDot = getSolidTemperatureSpan(uDot, column->Ngrid, column->Ncomp);
    auto spanWallTemperatureDot = getWallTemperatureSpan(uDot, column->Ngrid, column->Ncomp);

    computeFirstDerivativesEnergyBalance(
        column->components, column->Ngrid, column->Ncomp, column->externalTemperature, column->voidFraction,
        column->particleDensity, column->particleDiameter, column->internalDiameter, column->outerDiameter,
        column->wallDensity, column->gasThermalConductivity, column->wallThermalConductivity,
        column->heatTransferGasSolid, column->heatTransferGasWall, column->heatTransferWallExternal,
        column->heatCapacityGas, column->heatCapacitySolid, column->heatCapacityWall, column->resolution,
        column->prefactorMassTransfer, column->interstitialGasVelocity, column->totalPressure, spanGasTemperature,
        spanGasTemperatureDot, spanSolidTemperature, spanSolidTemperatureDot, spanWallTemperature,
        spanWallTemperatureDot, spanPartialPressure, spanPartialPressureDot, spanAdsorption, spanAdsorptionDot,
        column->equilibriumAdsorption, column->moleFraction, column->gasDensity, column->coeffGasGas,
        column->coeffGasSolid, column->coeffGasWall, column->coeffDiffusion, column->facePressures, column->energyFlow);
  }
  else
  {
    computeFirstDerivatives(column->components, column->Ngrid, column->Ncomp, column->resolution,
                            column->prefactorMassTransfer, column->interstitialGasVelocity, spanPartialPressure,
                            spanPartialPressureDot, spanAdsorption, spanAdsorptionDot, column->equilibriumAdsorption,
                            column->totalPressureDot);
  }

  return 0;
}

// static int JacSparseCSC(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix Jac, void* user_data, N_Vector tmp1,
//                         N_Vector tmp2, N_Vector tmp3)
// {
//   (void)t;
//   (void)y;
//   (void)fy;
//   (void)tmp1;
//   (void)tmp2;
//   (void)tmp3;

//   Column& column = *static_cast<Column*>(user_data);

//   const std::vector<Component>& components = column.components;
//   const size_t Ngrid = column.Ngrid, Ncomp = column.Ncomp, nNodes = Ngrid + 1, offset = nNodes * Ncomp, N = 2 *
//   offset; const double idx = 1.0 / column.resolution, idx2 = idx * idx;

//   auto* colptr = SM_INDEXPTRS_S(Jac);
//   auto* rowind = SM_INDEXVALS_S(Jac);
//   auto* data = SM_DATA_S(Jac);

//   for (sunindextype k = 0, nnz = SM_NNZ_S(Jac); k < nnz; ++k) data[k] = 0.0;

//   for (size_t col = 0; col < N; ++col)
//   {
//     const size_t local = (col < offset) ? col : (col - offset);
//     const size_t grid = local / Ncomp;
//     const size_t comp = local % Ncomp;

//     for (sunindextype k = colptr[col]; k < colptr[col + 1]; ++k)
//     {
//       const sunindextype row = rowind[k];
//       const bool isProw = (row < static_cast<sunindextype>(offset));
//       const size_t rlocal = isProw ? static_cast<size_t>(row) : static_cast<size_t>(row) - offset;
//       const size_t rgrid = rlocal / Ncomp;
//       const size_t rcomp = rlocal % Ncomp;

//       if (rcomp != comp) continue;

//       if (col < offset)
//       {
//         if (!isProw || rgrid == 0) continue;

//         if (rgrid < Ngrid)
//         {
//           if (grid + 1 == rgrid)
//             data[k] = interstitialGasVelocity[rgrid - 1] * idx + components[comp].D * idx2;
//           else if (grid == rgrid)
//             data[k] = -interstitialGasVelocity[rgrid] * idx - 2.0 * components[comp].D * idx2;
//           else if (grid == rgrid + 1)
//             data[k] = components[comp].D * idx2;
//         }
//         else
//         {
//           if (grid + 1 == Ngrid)
//             data[k] = interstitialGasVelocity[Ngrid - 1] * idx + components[comp].D * idx2;
//           else if (grid == Ngrid)
//             data[k] = -interstitialGasVelocity[Ngrid] * idx - components[comp].D * idx2;
//         }
//       }
//       else
//       {
//         if (isProw)
//           data[k] = (rgrid >= 1 && rgrid == grid) ? prefactorMassTransfer[comp] : 0.0;
//         else
//           data[k] = (rgrid == grid) ? -components[comp].Kl : 0.0;
//       }
//     }
//   }

//   return 0;
// }
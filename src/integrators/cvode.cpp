// cvode.cpp
#include "cvode.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <limits>
#include <print>
#include <span>
#include <string>
#include <type_traits>

#include "column.h"
#include "compute.h"
#include "sorption.h"
#include "transport.h"

#if BUILD_SUNDIALS

CVODE::~CVODE()
{
  if (linSolver) SUNLinSolFree(linSolver);
  if (linearMatrix) SUNMatDestroy(linearMatrix);
  if (solver) SUNNonlinSolFree(solver);
  if (absoluteToleranceVector) N_VDestroy(absoluteToleranceVector);
  if (stateDerivativeVector) N_VDestroy(stateDerivativeVector);
  if (stateVector) N_VDestroy(stateVector);
  if (cvodeMem) CVodeFree(&cvodeMem);
  if (sunLogger) SUNLogger_Destroy(&sunLogger);
  if (sunContext) SUNContext_Free(&sunContext);
}

namespace
{
void advanceCVODE(CVODE& integrator, Column& column, double endTime, Timing& timings)
{
  sunrealtype tReturn = integrator.currentTime;
  auto timer = timings.scoped(timings.total);

  int flag = CVode(integrator.cvodeMem, endTime, integrator.stateVector, &tReturn, CV_NORMAL);
  if (flag < 0)
  {
    throw std::runtime_error("CVode failed during propagation");
  }
  integrator.currentTime = tReturn;
  clampNonnegative(column.state);

  // Refreshes the derived column fields (pressure, mole fractions, rates) for the interpolated state.
  flag = CVODE::evaluateDerivatives(tReturn, integrator.stateVector, integrator.stateDerivativeVector, &column);
  if (flag != 0)
  {
    throw std::runtime_error("CVODE derivative evaluation failed after propagation");
  }
}

/**
 * \brief Smallest spacing between neighbouring grid nodes in m.
 */
double minimumGridSpacing(const Column& column)
{
  double spacing = std::numeric_limits<double>::max();
  for (size_t grid = 1; grid < column.columnDistances.size(); ++grid)
  {
    const double delta = column.columnDistances[grid] - column.columnDistances[grid - 1];
    if (delta > 0.0) spacing = std::min(spacing, delta);
  }
  return spacing == std::numeric_limits<double>::max() ? column.resolution : spacing;
}

/**
 * \brief Largest advective speed available at initialization in m/s.
 */
double characteristicVelocity(const Column& column)
{
  double velocity = std::abs(column.columnEntranceVelocity);
  for (const double v : column.interstitialGasVelocity) velocity = std::max(velocity, std::abs(v));
  return velocity;
}

/**
 * \brief Resolves the step cap handed to CVODE, in s. Returns 0.0 when the integrator should stay uncapped.
 *
 * A negative CVODE::maximumTimeStep selects the automatic value: the time the front needs to traverse a
 * single grid cell. Without a cap CVODE grows its step freely while the outlet is quiescent, and the first
 * step that meets the breakthrough front is then taken at high order from a history built entirely before
 * the front existed, which shows up as a non-monotonic hook in the breakthrough curve.
 */
double resolveMaximumTimeStep(const CVODE& integrator, const Column& column)
{
  if (integrator.maximumTimeStep > 0.0) return integrator.maximumTimeStep;
  if (integrator.maximumTimeStep == 0.0) return 0.0;

  const double spacing = minimumGridSpacing(column);
  const double velocity = characteristicVelocity(column);
  if (!(spacing > 0.0) || !(velocity > 0.0)) return 0.0;

  return spacing / velocity;
}

/**
 * \brief Fills the per-component absolute tolerance vector from the column's state layout.
 *
 * The state vector mixes mol/m^3, mol/kg and K, so a single scalar absolute tolerance is meaningless for at
 * least one of the blocks. Each block gets the tolerance appropriate to its own units instead.
 */
void fillAbsoluteTolerances(const CVODE& integrator, const Column& column, N_Vector tolerances)
{
  // Concentration, surface-concentration and pore-concentration blocks all share mol/m^3, so seed the whole
  // vector with that value and overwrite only the blocks carrying different units.
  N_VConst(integrator.absoluteToleranceConcentration, tolerances);

  double* base = static_cast<double*>(N_VGetArrayPointer(tolerances));
  const auto layout = column.stateLayout();

  const auto fill = [](std::span<double> block, double value)
  { std::fill(block.begin(), block.end(), value); };

  fill(layout.physisorption(base), integrator.absoluteToleranceLoading);
  fill(layout.chemisorption(base), integrator.absoluteToleranceLoading);
  fill(layout.gasTemperature(base), integrator.absoluteToleranceTemperature);
  fill(layout.solidTemperature(base), integrator.absoluteToleranceTemperature);
  fill(layout.wallTemperature(base), integrator.absoluteToleranceTemperature);
}

void initializeCVODE(CVODE& integrator, Column& column)
{
  static_assert(std::is_same_v<sunrealtype, double>, "CVODE currently requires SUNDIALS sunrealtype to be double");

  SUNContext_Create(SUN_COMM_NULL, &integrator.sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &integrator.sunLogger);
  SUNContext_SetLogger(integrator.sunContext, integrator.sunLogger);

  const sunindextype totalSize = static_cast<sunindextype>(column.state.size());
  integrator.stateVector = N_VMake_Serial(totalSize, column.state.data(), integrator.sunContext);
  integrator.stateDerivativeVector = N_VMake_Serial(totalSize, column.stateDot.data(), integrator.sunContext);
  if (!integrator.stateVector || !integrator.stateDerivativeVector)
  {
    throw std::runtime_error("Failed to allocate CVODE state vectors");
  }

  integrator.cvodeMem = CVodeCreate(CV_BDF, integrator.sunContext);
  CVodeSetMaxNumSteps(integrator.cvodeMem, 100000);
  CVodeSetUserData(integrator.cvodeMem, &column);
  integrator.currentTime = 0.0;

  int flag = CVodeInit(integrator.cvodeMem, CVODE::evaluateDerivatives, integrator.currentTime, integrator.stateVector);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeInit failed");

  integrator.absoluteToleranceVector = N_VNew_Serial(totalSize, integrator.sunContext);
  if (!integrator.absoluteToleranceVector)
  {
    throw std::runtime_error("Failed to allocate CVODE absolute-tolerance vector");
  }
  fillAbsoluteTolerances(integrator, column, integrator.absoluteToleranceVector);

  flag = CVodeSVtolerances(integrator.cvodeMem, integrator.relativeTolerance, integrator.absoluteToleranceVector);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSVtolerances failed");

  flag = CVodeSetMaxOrd(integrator.cvodeMem, integrator.maximumOrder);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetMaxOrd failed");

  integrator.appliedMaximumTimeStep = resolveMaximumTimeStep(integrator, column);
  if (integrator.appliedMaximumTimeStep > 0.0)
  {
    flag = CVodeSetMaxStep(integrator.cvodeMem, integrator.appliedMaximumTimeStep);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetMaxStep failed");
  }

  std::print("CVODE: {} solver, rtol {:.3e}, atol {:.3e} [mol/m^3] / {:.3e} [mol/kg] / {:.3e} [K], max step {}\n",
             integrator.linearSolverType == CVODE::LinearSolverType::SPGMR
                 ? std::format("SPGMR (matrix-free, Krylov dim {})", integrator.krylovDimension)
                 : std::string("dense"),
             integrator.relativeTolerance, integrator.absoluteToleranceConcentration,
             integrator.absoluteToleranceLoading, integrator.absoluteToleranceTemperature,
             integrator.appliedMaximumTimeStep > 0.0 ? std::format("{:.3e} s", integrator.appliedMaximumTimeStep)
                                                     : std::string("unlimited"));

  integrator.solver = SUNNonlinSol_Newton(integrator.stateVector, integrator.sunContext);
  flag = CVodeSetNonlinearSolver(integrator.cvodeMem, integrator.solver);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetNonlinearSolver failed");

  if (integrator.linearSolverType == CVODE::LinearSolverType::SPGMR)
  {
    // Matrix-free: no SUNMatrix is created and none is handed to CVODE, which then forms Jacobian-vector
    // products by difference quotients (one right-hand-side evaluation each) instead of assembling and
    // factoring J. Left unpreconditioned for now; a node-local block preconditioner goes here next.
    integrator.linSolver =
        SUNLinSol_SPGMR(integrator.stateVector, SUN_PREC_NONE, integrator.krylovDimension, integrator.sunContext);
    if (!integrator.linSolver) throw std::runtime_error("SUNLinSol_SPGMR failed");

    if (SUNLinSol_SPGMRSetGSType(integrator.linSolver, SUN_MODIFIED_GS) != SUN_SUCCESS)
    {
      throw std::runtime_error("SUNLinSol_SPGMRSetGSType failed");
    }

    flag = CVodeSetLinearSolver(integrator.cvodeMem, integrator.linSolver, nullptr);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");
  }
  else
  {
    integrator.linearMatrix = SUNDenseMatrix(totalSize, totalSize, integrator.sunContext);
    integrator.linSolver = SUNLinSol_Dense(integrator.stateVector, integrator.linearMatrix, integrator.sunContext);
    flag = CVodeSetLinearSolver(integrator.cvodeMem, integrator.linSolver, integrator.linearMatrix);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");

    CVodeSetJacFn(integrator.cvodeMem, nullptr);
  }
}
}  // namespace

void CVODE::advance(Column& column, double endTime, Timing& timings)
{
  advanceCVODE(*this, column, endTime, timings);
}

void CVODE::initialize(Column& column) { initializeCVODE(*this, column); }

void CVODE::reinitialize()
{
  const int flag = CVodeReInit(cvodeMem, currentTime, stateVector);
  if (flag != CV_SUCCESS) throw std::runtime_error("CVodeReInit failed");
}

namespace
{
void reinitializeCVODE(CVODE& integrator, const Column& column)
{
  // Only the automatic cap (negative maximumTimeStep) follows the velocity; a user-given cap is kept.
  if (integrator.maximumTimeStep < 0.0)
  {
    integrator.appliedMaximumTimeStep = resolveMaximumTimeStep(integrator, column);
    // A zero cap restores CVODE's default, an unlimited step.
    const int flag = CVodeSetMaxStep(integrator.cvodeMem, integrator.appliedMaximumTimeStep);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetMaxStep failed");
  }
  integrator.reinitialize();
}
}  // namespace

void CVODE::reinitialize(const Column& column) { reinitializeCVODE(*this, column); }

void CVODE::printStatistics() const
{
  if (!cvodeMem) return;

  // Every getter returns a CVODE flag; a counter that is unavailable for the active linear solver reports
  // as -1 rather than aborting the report.
  // Deduced rather than typed as int(*)(void*, long int*): the SUNDIALS getters have C language linkage,
  // and a pointer to them is formally a distinct type from the C++-linkage equivalent.
  const auto counter = [this](auto getter)
  {
    long int value = 0;
    return getter(cvodeMem, &value) == CV_SUCCESS ? value : -1L;
  };

  int lastOrder = 0;
  CVodeGetLastOrder(cvodeMem, &lastOrder);
  sunrealtype lastStep = 0.0;
  CVodeGetLastStep(cvodeMem, &lastStep);

  std::print("CVODE statistics:\n");
  std::print("  internal steps:            {}\n", counter(CVodeGetNumSteps));
  std::print("  rhs evaluations:           {} (+{} inside the linear solver)\n", counter(CVodeGetNumRhsEvals),
             counter(CVodeGetNumLinRhsEvals));
  std::print("  linear solver setups:      {}\n", counter(CVodeGetNumLinSolvSetups));
  std::print("  newton iterations:         {} (convergence failures: {})\n", counter(CVodeGetNumNonlinSolvIters),
             counter(CVodeGetNumNonlinSolvConvFails));
  std::print("  error test failures:       {}\n", counter(CVodeGetNumErrTestFails));
  std::print("  last order / step:         {} / {:.6e} s\n", lastOrder, lastStep);

  if (linearSolverType == LinearSolverType::SPGMR)
  {
    const long int iterations = counter(CVodeGetNumLinIters);
    const long int setups = counter(CVodeGetNumLinSolvSetups);
    std::print("  krylov iterations:         {} (convergence failures: {})\n", iterations,
               counter(CVodeGetNumLinConvFails));
    if (setups > 0 && iterations >= 0)
    {
      std::print("  krylov iterations / setup: {:.2f}\n",
                 static_cast<double>(iterations) / static_cast<double>(setups));
    }
    std::print("  preconditioner evals:      {}\n", counter(CVodeGetNumPrecEvals));
  }
  else
  {
    std::print("  jacobian evaluations:      {}\n", counter(CVodeGetNumJacEvals));
  }
}

/*
 * CVODE evaluates intermediate trial states in storage owned by SUNDIALS. The
 * callbacks below therefore bind layout views to the supplied N_Vector rather
 * than using the column's accepted-state spans directly.
 */
int CVODE::evaluateDerivatives(sunrealtype t, N_Vector stateVector, N_Vector stateDerivativeVector, void* user_data)
{
  auto* column = static_cast<Column*>(user_data);
  N_VConst(0.0, stateDerivativeVector);

  const ColumnStateLayout layout = column->stateLayout();
  double* stateBase = static_cast<double*>(N_VGetArrayPointer(stateVector));
  double* derivativeBase = static_cast<double*>(N_VGetArrayPointer(stateDerivativeVector));

  auto spanConcentration = layout.concentration(stateBase);
  auto spanPhysisorption = layout.physisorption(stateBase);
  auto spanChemisorption = layout.chemisorption(stateBase);
  const bool usePoreSurfaceTransport = column->surfacePoreTransportEnabled;
  auto spanSurfaceConcentration =
      usePoreSurfaceTransport ? layout.surfaceConcentration(stateBase) : column->surfaceConcentration;
  auto spanPoreConcentration =
      usePoreSurfaceTransport ? layout.poreConcentration(stateBase) : column->poreConcentration;

  auto spanConcentrationDot = layout.concentration(derivativeBase);
  auto spanPhysisorptionDot = layout.physisorption(derivativeBase);
  auto spanChemisorptionDot = layout.chemisorption(derivativeBase);
  auto spanSurfaceConcentrationDot =
      usePoreSurfaceTransport ? layout.surfaceConcentration(derivativeBase) : column->surfaceConcentrationDot;
  auto spanPoreConcentrationDot =
      usePoreSurfaceTransport ? layout.poreConcentration(derivativeBase) : column->poreConcentrationDot;

  auto spanGasTemperature = column->energyBalance ? layout.gasTemperature(stateBase) : column->gasTemperature;
  auto spanSolidTemperature = column->energyBalance ? layout.solidTemperature(stateBase) : column->solidTemperature;
  auto spanWallTemperature = column->energyBalance ? layout.wallTemperature(stateBase) : column->wallTemperature;
  auto spanGasTemperatureDot =
      column->energyBalance ? layout.gasTemperature(derivativeBase) : column->gasTemperatureDot;
  auto spanSolidTemperatureDot =
      column->energyBalance ? layout.solidTemperature(derivativeBase) : column->solidTemperatureDot;
  auto spanWallTemperatureDot =
      column->energyBalance ? layout.wallTemperature(derivativeBase) : column->wallTemperatureDot;

  std::fill(column->reactionPhysisorptionSource.begin(), column->reactionPhysisorptionSource.end(), 0.0);
  std::fill(column->reactionChemisorptionSource.begin(), column->reactionChemisorptionSource.end(), 0.0);
  std::fill(column->reactionPoreConcentrationSource.begin(), column->reactionPoreConcentrationSource.end(), 0.0);
  std::fill(column->reactionBulkSource.begin(), column->reactionBulkSource.end(), 0.0);
  std::fill(column->reactionHeat.begin(), column->reactionHeat.end(), 0.0);

  // The second pass feeds the uptake sink of the first back into the velocity. With a fixed velocity, a fixed
  // pressure profile, or a liquid, the velocity and pressure do not depend on the sink, and one pass suffices.
  const bool sinkSetsVelocity = column->fluidPhase == Column::FluidPhase::Gas &&
                                column->boundaryCondition != Column::BoundaryCondition::FixedVelocity &&
                                column->boundaryCondition != Column::BoundaryCondition::FixedPressureInletVelocity;
  for (size_t iteration = 0; iteration < (sinkSetsVelocity ? 2 : 1); ++iteration)
  {
    computeBulkSpeciesSink(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                           column->numberOfAdsorbents, column->maxChemisorptionSites, column->geometries,
                           column->particleDensities, column->fractionOfAdsorbent, column->totalVoidFraction,
                           spanConcentration, spanPhysisorptionDot, spanChemisorptionDot, spanSurfaceConcentration,
                           column->bulkSpeciesSink, column->reactionPhysisorptionSource,
                           column->reactionChemisorptionSource, column->reactionBulkSource);

    updateVelocityAndPressure(
        column->components, column->boundaryCondition, column->numberOfGridPoints, column->numberOfComponents,
        column->inletPressure, column->outletPressure, column->pressureGradient, column->columnLength,
        column->numberOfAdsorbents, column->columnEntranceVelocity, column->dynamicViscosity, column->columnDistances,
        column->fractionOfAdsorbent, column->geometries, column->totalVoidFraction, column->interstitialGasVelocity,
        column->gasDensity, column->totalConcentration, column->totalPressure, spanConcentration,
        column->partialPressure, column->moleFraction, column->bulkSpeciesSink, spanGasTemperature, column->fluidPhase,
        column->liquidDensity, column->pHMode, column->pHValue, column->pKw, column->pHComponent, column->pH);

    computePhysisorptionEquilibriumLoadings(
        column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->fractionOfAdsorbent, column->hasAdsorbentOfType, column->maxIsothermTerms,
        column->iastPerformance, column->idealGasMolFractions, column->adsorbedMolFractions, column->numberOfMolecules,
        column->totalPressure, column->equilibriumPhysisorption, column->cachedPressure, column->cachedGrandPotential,
        column->moleFraction, spanGasTemperature,
        column->fluidPhase == Column::FluidPhase::Gas
            ? MixturePrediction::DrivingForceInput::MoleFraction
            : MixturePrediction::DrivingForceInput::Concentration,
        spanConcentration, column->pH);

    computeChemisorptionEquilibriumLoadings(
        column->chemisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->fractionOfAdsorbent, column->hasAdsorbentOfType,
        column->maxChemisorptionSites, column->iastPerformance, column->idealGasMolFractions,
        column->adsorbedMolFractions, column->numberOfMolecules, column->totalPressure,
        column->equilibriumChemisorption, column->cachedChemisorptionPressure,
        column->cachedChemisorptionGrandPotential, column->moleFraction, spanGasTemperature,
        column->fluidPhase == Column::FluidPhase::Gas
            ? MixturePrediction::DrivingForceInput::MoleFraction
            : MixturePrediction::DrivingForceInput::Concentration,
        spanConcentration, column->pH);

    computePhysisorption(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                         column->numberOfAdsorbents, column->fractionOfAdsorbent, column->equilibriumPhysisorption,
                         spanPhysisorption, spanPhysisorptionDot);
    computeChemisorption(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                         column->numberOfAdsorbents, column->maxChemisorptionSites, column->externalTemperature, t,
                         column->fractionOfAdsorbent, column->particleDensities, column->equilibriumChemisorption,
                         spanConcentration, spanChemisorption, spanChemisorptionDot, spanPoreConcentration,
                         spanSolidTemperature);
    computeChemisorptionTransportDerivatives(
        column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->maxChemisorptionSites, column->fractionOfAdsorbent, column->geometries,
        column->particleDensities, spanConcentration, spanChemisorptionDot, spanSurfaceConcentration,
        spanSurfaceConcentrationDot, spanPoreConcentration, spanPoreConcentrationDot);

    computeReactionDerivatives(
        column->physisorptionMixtures.front().components, column->reactions, column->numberOfGridPoints,
        column->numberOfComponents, column->maxChemisorptionSites, column->externalTemperature, spanPhysisorption,
        spanPhysisorptionDot, spanChemisorption, spanChemisorptionDot, spanPoreConcentration, spanPoreConcentrationDot,
        spanSolidTemperature, column->reactionPhysisorptionSource, column->reactionChemisorptionSource,
        column->reactionPoreConcentrationSource, column->reactionHeat, spanConcentration, column->reactionBulkSource,
        column->fractionOfAdsorbent, column->numberOfAdsorbents);

    computeBulkSpeciesSink(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                           column->numberOfAdsorbents, column->maxChemisorptionSites, column->geometries,
                           column->particleDensities, column->fractionOfAdsorbent, column->totalVoidFraction,
                           spanConcentration, spanPhysisorptionDot, spanChemisorptionDot, spanSurfaceConcentration,
                           column->bulkSpeciesSink, column->reactionPhysisorptionSource,
                           column->reactionChemisorptionSource, column->reactionBulkSource);
  }

  computeMassDerivatives(column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
                         column->numberOfAdsorbents, column->columnDistances, column->fractionOfAdsorbent,
                         column->interstitialGasVelocity, spanConcentration, spanConcentrationDot,
                         column->bulkSpeciesSink, column->fluidPhase, column->totalVoidFraction,
                         column->convectionScheme);

  if (column->energyBalance)
  {
    computeEnergyDerivatives(
        column->physisorptionMixtures, column->numberOfGridPoints, column->numberOfComponents,
        column->numberOfAdsorbents, column->externalTemperature, column->geometries, column->particleDensities,
        column->fractionOfAdsorbent, column->wallDensity, column->gasThermalConductivity,
        column->wallThermalConductivity, column->heatTransferGasSolid, column->heatTransferGasWall,
        column->heatTransferWallExternal, column->heatCapacityGas, column->heatCapacitySolid, column->heatCapacityWall,
        column->columnDistances, column->interstitialGasVelocity, column->gasDensity, column->coeffDiffusion,
        column->maxChemisorptionSites, spanPhysisorptionDot, spanChemisorptionDot, spanGasTemperature,
        spanGasTemperatureDot, spanSolidTemperature, spanSolidTemperatureDot, spanWallTemperature,
        spanWallTemperatureDot, column->reactionPhysisorptionSource, column->reactionChemisorptionSource,
        column->reactionHeat);
  }

  return 0;
}

#else

CVODE::~CVODE() = default;

void CVODE::advance(Column&, double, Timing&)
{
  throw std::runtime_error("CVODE::advance() called, but this build was compiled without SUNDIALS support");
}

void CVODE::initialize(Column&) {}

void CVODE::reinitialize() {}

void CVODE::reinitialize(const Column&) {}

void CVODE::printStatistics() const {}

#endif

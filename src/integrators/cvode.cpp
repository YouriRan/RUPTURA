// cvode.cpp
#include "cvode.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <print>
#include <span>
#include <string>
#include <type_traits>

#include "column.h"
#include "rk3.h"

#if BUILD_SUNDIALS

extern "C"
{
  void dgbtrf_(int* m, int* n, int* kl, int* ku, double* ab, int* ldab, int* ipiv, int* info);
  void dgbtrs_(char* trans, int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb,
               int* info);
}

CVODE::~CVODE()
{
  if (preconditionerScratch) N_VDestroy(preconditionerScratch);
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
// Internal steps of one CVode call (CVodeSetMaxNumSteps); advanceCVODE continues a call that reaches the limit.
constexpr long maximumStepsPerCall = 100000;

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
  return spacing;
}

void advanceCVODE(CVODE& integrator, Column& column, double endTime, Timing& timings)
{
  sunrealtype tReturn = integrator.currentTime;
  auto timer = timings.scoped(timings.total);
  integrator.userData.timings = &timings;

  integrator.userData.lastError.clear();
  // A call that reaches maximumStepsPerCall is continued: an output interval that is long against the step can need
  // more steps (B-ITG-19). It stops when those steps advanced the time by less than 1e-3 of the time in which the
  // fluid crosses the smallest cell, or of the step cap if that is smaller: CVODE then crawls, as at steep fronts on
  // an empty bed at BDF orders > 1 (B-ITG-11). Without flow there is no crossing time, and only CVODE's own failures
  // stop the run.
  int flag = CV_TOO_MUCH_WORK;
  while (flag == CV_TOO_MUCH_WORK)
  {
    const double callStartTime = tReturn;
    flag = CVode(integrator.cvodeMem, endTime, integrator.stateVector, &tReturn, CV_NORMAL);
    if (flag != CV_TOO_MUCH_WORK) break;

    double fastest = 0.0;
    for (const double velocity : column.interstitialGasVelocity) fastest = std::max(fastest, std::abs(velocity));
    if (fastest <= 0.0) continue;
    double referenceTime = minimumGridSpacing(column) / fastest;
    if (integrator.appliedMaximumTimeStep > 0.0)
    {
      referenceTime = std::min(referenceTime, integrator.appliedMaximumTimeStep);
    }
    const double meanStep = (tReturn - callStartTime) / static_cast<double>(maximumStepsPerCall);
    if (meanStep < 1.0e-3 * referenceTime)
    {
      throw std::runtime_error(std::format(
          "CVODE crawls at t = {} s: {} steps advanced the time by {:.3e} s (mean step {:.3e} s, {:.1e} of the "
          "{:.3e} s in which the fluid crosses the smallest cell or of the step cap). At steep fronts on an empty "
          "bed, CVODEMaximumOrder 1 or the ROS34PW2 integrator can help.",
          tReturn, maximumStepsPerCall, tReturn - callStartTime, meanStep, meanStep / referenceTime, referenceTime));
    }
  }
  if (flag < 0)
  {
    char* name = CVodeGetReturnFlagName(flag);
    std::string message = std::format("CVode failed during propagation at t = {} s ({})", tReturn, name);
    std::free(name);
    if (!integrator.userData.lastError.empty())
    {
      message += "; last right-hand-side error: " + integrator.userData.lastError;
    }
    throw std::runtime_error(message);
  }
  integrator.currentTime = tReturn;
  clampNonnegative(column.state);

  // Refreshes the derived column fields (pressure, mole fractions, rates) for the interpolated state.
  flag = CVODE::evaluateDerivatives(tReturn, integrator.stateVector, integrator.stateDerivativeVector,
                                    &integrator.userData);
  if (flag != 0)
  {
    throw std::runtime_error("CVODE derivative evaluation failed after propagation: " + integrator.userData.lastError);
  }
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
 * A negative CVODE::maximumTimeStep selects the automatic value: the time the front needs to traverse four
 * grid cells. Without a cap CVODE grows its step freely while the outlet is quiescent, and the first
 * step that meets the breakthrough front is then taken at high order from a history built entirely before
 * the front existed, which shows up as a non-monotonic hook in the breakthrough curve. Up to four cells the
 * error stays at that of a one-cell cap (benchmark/cvode-max-step); from eight cells on it grows in the monolith.
 */
double resolveMaximumTimeStep(const CVODE& integrator, const Column& column)
{
  if (integrator.maximumTimeStep > 0.0) return integrator.maximumTimeStep;
  if (integrator.maximumTimeStep == 0.0) return 0.0;

  const double spacing = minimumGridSpacing(column);
  const double velocity = characteristicVelocity(column);
  if (!(spacing > 0.0) || !(velocity > 0.0)) return 0.0;

  return 4.0 * spacing / velocity;
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

/**
 * \brief KLU: the compressed-sparse-column pattern of the Jacobian in the field-major CVODE state.
 *
 * A concentration column has the concentrations and loadings of its node (sorption) and the concentrations that its
 * transport reaches: the same component at the neighbouring nodes with upwind convection; with WENO3 the nodes from
 * one upstream to two downstream, in the gas phase every component; a loading column has the concentration and
 * loading of its component and node; every other state has its diagonal only, so that CVODE can form I - gamma J in
 * place. Each nonzero records its index in the band storage of SorptionJacobian, or SIZE_MAX where the Jacobian has
 * no entry.
 */
void buildSparsePattern(CVODE& integrator, const Column& column)
{
  const SorptionJacobian& jacobian = *integrator.sorptionJacobian;
  const size_t n = column.numberOfComponents;
  const size_t nodes = column.numberOfGridPoints + 1;
  const size_t blockSize = nodes * n;
  const size_t concentrationOffset = static_cast<size_t>(column.concentration.data() - column.state.data());
  const size_t loadingOffset = static_cast<size_t>(column.physisorption.data() - column.state.data());
  // Node-major index of a field-major state index, or SIZE_MAX outside the concentration and loading blocks.
  const auto nodeMajor = [&](size_t index)
  {
    if (index >= concentrationOffset && index < concentrationOffset + blockSize)
    {
      return jacobian.concentrationIndex((index - concentrationOffset) / n, (index - concentrationOffset) % n);
    }
    if (index >= loadingOffset && index < loadingOffset + blockSize)
    {
      return jacobian.loadingIndex((index - loadingOffset) / n, (index - loadingOffset) % n);
    }
    return std::numeric_limits<size_t>::max();
  };

  integrator.sparseColumnStarts.assign(1, 0);
  integrator.sparseRowIndices.clear();
  integrator.sparseBandIndex.clear();
  std::vector<size_t> rows;
  for (size_t col = 0; col < column.state.size(); ++col)
  {
    rows.assign(1, col);
    if (col >= concentrationOffset && col < concentrationOffset + blockSize)
    {
      const size_t node = (col - concentrationOffset) / n;
      const size_t comp = (col - concentrationOffset) % n;
      for (size_t other = 0; other < n; ++other)
      {
        rows.push_back(concentrationOffset + node * n + other);
        rows.push_back(loadingOffset + node * n + other);
      }
      // A WENO3 face flux reaches two nodes upstream, on either side with the sign of the velocity.
      const bool weno = column.convectionScheme == Column::ConvectionScheme::WENO3;
      const bool allComponents = weno && column.fluidPhase == Column::FluidPhase::Gas;
      const size_t reach = weno ? 2 : 1;
      for (size_t row = node > reach ? node - reach : 0; row < std::min(nodes, node + reach + 1); ++row)
      {
        for (size_t other = 0; other < n; ++other)
        {
          if (allComponents || other == comp) rows.push_back(concentrationOffset + row * n + other);
        }
      }
    }
    else if (col >= loadingOffset && col < loadingOffset + blockSize)
    {
      rows.push_back(concentrationOffset + (col - loadingOffset));
    }
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

    const size_t columnNodeMajor = nodeMajor(col);
    for (size_t row : rows)
    {
      const size_t rowNodeMajor = nodeMajor(row);
      integrator.sparseRowIndices.push_back(static_cast<sunindextype>(row));
      integrator.sparseBandIndex.push_back(
          rowNodeMajor == std::numeric_limits<size_t>::max() || columnNodeMajor == std::numeric_limits<size_t>::max()
              ? std::numeric_limits<size_t>::max()
              : 2 * jacobian.bandwidth + rowNodeMajor - columnNodeMajor + columnNodeMajor * jacobian.leadingDimension);
    }
    integrator.sparseColumnStarts.push_back(static_cast<sunindextype>(integrator.sparseRowIndices.size()));
  }
}

void initializeCVODE(CVODE& integrator, Column& column)
{
  static_assert(std::is_same_v<sunrealtype, double>, "CVODE currently requires SUNDIALS sunrealtype to be double");

  SUNContext_Create(SUN_COMM_NULL, &integrator.sunContext);
  SUNLogger_Create(SUN_COMM_NULL, 0, &integrator.sunLogger);
  SUNContext_SetLogger(integrator.sunContext, integrator.sunLogger);
  // The default handler logs every error; reaching the step limit of one call is not one (advanceCVODE continues).
  SUNContext_ClearErrHandlers(integrator.sunContext);
  SUNContext_PushErrHandler(
      integrator.sunContext,
      [](int line, const char* function, const char* file, const char* message, SUNErrCode code, void* data,
         SUNContext context)
      {
        if (code != CV_TOO_MUCH_WORK) SUNLogErrHandlerFn(line, function, file, message, code, data, context);
      },
      nullptr);

  const sunindextype totalSize = static_cast<sunindextype>(column.state.size());
  integrator.stateVector = N_VMake_Serial(totalSize, column.state.data(), integrator.sunContext);
  integrator.stateDerivativeVector = N_VMake_Serial(totalSize, column.stateDot.data(), integrator.sunContext);
  if (!integrator.stateVector || !integrator.stateDerivativeVector)
  {
    throw std::runtime_error("Failed to allocate CVODE state vectors");
  }

  integrator.cvodeMem = CVodeCreate(CV_BDF, integrator.sunContext);
  CVodeSetMaxNumSteps(integrator.cvodeMem, maximumStepsPerCall);
  integrator.userData = {&column, &integrator};
  CVodeSetUserData(integrator.cvodeMem, &integrator.userData);
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
                 ? std::format("SPGMR (matrix-free, Krylov dim {}{})", integrator.krylovDimension,
                               integrator.preconditionerType == CVODE::PreconditionerType::SorptionJacobian
                                   ? ", sorption-Jacobian preconditioner"
                                   : "")
                 : std::string(integrator.linearSolverType == CVODE::LinearSolverType::KLU
                                   ? "KLU (sparse direct, sorption Jacobian)"
                                   : "dense"),
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
    // factoring J. The optional left preconditioner factorizes the banded sorption Jacobian.
    const bool preconditioned = integrator.preconditionerType == CVODE::PreconditionerType::SorptionJacobian;
    integrator.linSolver = SUNLinSol_SPGMR(integrator.stateVector, preconditioned ? SUN_PREC_LEFT : SUN_PREC_NONE,
                                           integrator.krylovDimension, integrator.sunContext);
    if (!integrator.linSolver) throw std::runtime_error("SUNLinSol_SPGMR failed");

    if (SUNLinSol_SPGMRSetGSType(integrator.linSolver, SUN_MODIFIED_GS) != SUN_SUCCESS)
    {
      throw std::runtime_error("SUNLinSol_SPGMRSetGSType failed");
    }

    flag = CVodeSetLinearSolver(integrator.cvodeMem, integrator.linSolver, nullptr);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");

    if (preconditioned)
    {
      integrator.sorptionJacobian.emplace(column);
      integrator.sorptionJacobianAssembled = false;
      integrator.preconditionerBand.resize(integrator.sorptionJacobian->band.size());
      integrator.preconditionerPivots.resize(integrator.sorptionJacobian->size);
      integrator.preconditionerVector.resize(integrator.sorptionJacobian->size);
      integrator.preconditionerScratch = N_VNew_Serial(totalSize, integrator.sunContext);
      if (!integrator.preconditionerScratch) throw std::runtime_error("Failed to allocate the preconditioner vector");
      flag = CVodeSetPreconditioner(integrator.cvodeMem, CVODE::setupPreconditioner, CVODE::solvePreconditioner);
      if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetPreconditioner failed");
    }
  }
  else if (integrator.linearSolverType == CVODE::LinearSolverType::KLU)
  {
    integrator.sorptionJacobian.emplace(column);
    buildSparsePattern(integrator, column);
    integrator.linearMatrix =
        SUNSparseMatrix(totalSize, totalSize, static_cast<sunindextype>(integrator.sparseRowIndices.size()), CSC_MAT,
                        integrator.sunContext);
    integrator.linSolver = SUNLinSol_KLU(integrator.stateVector, integrator.linearMatrix, integrator.sunContext);
    if (!integrator.linearMatrix || !integrator.linSolver) throw std::runtime_error("SUNLinSol_KLU failed");
    flag = CVodeSetLinearSolver(integrator.cvodeMem, integrator.linSolver, integrator.linearMatrix);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetLinearSolver failed");
    flag = CVodeSetJacFn(integrator.cvodeMem, CVODE::evaluateJacobian);
    if (flag != CV_SUCCESS) throw std::runtime_error("CVodeSetJacFn failed");
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
    std::print("  preconditioner evals:      {} (solves: {})\n", counter(CVodeGetNumPrecEvals),
               counter(CVodeGetNumPrecSolves));
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
  auto* data = static_cast<CVODE::UserData*>(user_data);
  Column& column = *data->column;
  Timing& timings = *data->timings;
  N_VConst(0.0, stateDerivativeVector);
  column.advancePressureRamp(t);

  // No exception may unwind through the SUNDIALS frames. A trial state for which the closures or the mixture
  // prediction fail (negative outlet pressure, no velocity bracket, failed IAST inversion) is a recoverable error:
  // CVODE retries with a smaller step (B-ITG-18). Other exceptions stop the integration; advanceCVODE reports the
  // message.
  try
  {
    ColumnStateViews views =
        column.stateLayout().views(static_cast<double*>(N_VGetArrayPointer(stateVector)),
                                   static_cast<double*>(N_VGetArrayPointer(stateDerivativeVector)));
    // Without the energy balance the temperatures are the fixed fields of the column.
    if (!column.energyBalance)
    {
      views.gasTemperature = column.gasTemperature;
      views.solidTemperature = column.solidTemperature;
      views.wallTemperature = column.wallTemperature;
    }

    // The first sink starts from zero rates, the reaction sources included.
    std::fill(column.reactionPhysisorptionSource.begin(), column.reactionPhysisorptionSource.end(), 0.0);
    std::fill(column.reactionChemisorptionSource.begin(), column.reactionChemisorptionSource.end(), 0.0);
    std::fill(column.reactionPoreConcentrationSource.begin(), column.reactionPoreConcentrationSource.end(), 0.0);
    std::fill(column.reactionBulkSource.begin(), column.reactionBulkSource.end(), 0.0);
    std::fill(column.reactionHeat.begin(), column.reactionHeat.end(), 0.0);

    // A first pass of the sorption rates feeds its uptake sink into the velocity and the equilibrium loadings of the
    // second. With a fixed velocity, a fixed pressure profile, or a liquid, the velocity and pressure do not depend on
    // the sink, and one pass suffices.
    const bool sinkSetsVelocity = column.fluidPhase == Column::FluidPhase::Gas &&
                                  column.boundaryCondition != Column::BoundaryCondition::FixedVelocity &&
                                  column.boundaryCondition != Column::BoundaryCondition::FixedPressureInletVelocity;
    precompute(column, views, timings);
    if (sinkSetsVelocity)
    {
      timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column, views, true); });
      precompute(column, views, timings);
    }
    timings.measure(timings.computeDerivatives, [&] { computeDerivatives(column, views); });
  }
  catch (const std::runtime_error& error)
  {
    data->lastError = error.what();
    return 1;
  }
  catch (const std::exception& error)
  {
    data->lastError = error.what();
    return -1;
  }

  return 0;
}

int CVODE::setupPreconditioner(sunrealtype t, N_Vector stateVector, N_Vector, sunbooleantype jacobianOk,
                               sunbooleantype* jacobianCurrent, sunrealtype gamma, void* user_data)
{
  auto* data = static_cast<CVODE::UserData*>(user_data);
  CVODE& integrator = *data->integrator;
  Column& column = *data->column;
  SorptionJacobian& jacobian = *integrator.sorptionJacobian;
  try
  {
    if (!jacobianOk || !integrator.sorptionJacobianAssembled)
    {
      // The right-hand side at the trial state leaves its velocity, pressure, mole fractions, and hypothetical
      // pressures in the column, from which the Jacobian is assembled.
      if (const int flag = evaluateDerivatives(t, stateVector, integrator.preconditionerScratch, user_data); flag != 0)
      {
        return flag;
      }
      const ColumnStateLayout layout = column.stateLayout();
      double* base = N_VGetArrayPointer(stateVector);
      assembleSorptionJacobian(column, layout.concentration(base),
                               column.energyBalance ? layout.gasTemperature(base) : column.gasTemperature, jacobian);
      integrator.sorptionJacobianAssembled = true;
      *jacobianCurrent = SUNTRUE;
    }
    else
    {
      *jacobianCurrent = SUNFALSE;
    }

    std::vector<double>& band = integrator.preconditionerBand;
    for (size_t index = 0; index < band.size(); ++index) band[index] = -gamma * jacobian.band[index];
    for (size_t row = 0; row < jacobian.size; ++row)
      band[2 * jacobian.bandwidth + row * jacobian.leadingDimension] += 1.0;
    int n = static_cast<int>(jacobian.size);
    int bands = static_cast<int>(jacobian.bandwidth);
    int leadingDimension = static_cast<int>(jacobian.leadingDimension);
    int info = 0;
    dgbtrf_(&n, &n, &bands, &bands, band.data(), &leadingDimension, integrator.preconditionerPivots.data(), &info);
    // A singular matrix is recoverable: CVODE retries with a fresh Jacobian or a smaller step.
    return info == 0 ? 0 : 1;
  }
  catch (const std::runtime_error& error)
  {
    data->lastError = error.what();
    return 1;
  }
  catch (const std::exception& error)
  {
    data->lastError = error.what();
    return -1;
  }
}

int CVODE::evaluateJacobian(sunrealtype t, N_Vector stateVector, N_Vector, SUNMatrix matrix, void* user_data,
                             N_Vector scratch, N_Vector, N_Vector)
{
  auto* data = static_cast<CVODE::UserData*>(user_data);
  CVODE& integrator = *data->integrator;
  Column& column = *data->column;
  SorptionJacobian& jacobian = *integrator.sorptionJacobian;
  try
  {
    // The right-hand side at the trial state leaves its velocity, pressure, mole fractions, and hypothetical
    // pressures in the column, from which the Jacobian is assembled.
    if (const int flag = evaluateDerivatives(t, stateVector, scratch, user_data); flag != 0) return flag;
    const ColumnStateLayout layout = column.stateLayout();
    double* base = N_VGetArrayPointer(stateVector);
    assembleSorptionJacobian(column, layout.concentration(base),
                             column.energyBalance ? layout.gasTemperature(base) : column.gasTemperature, jacobian);

    // SUNMatZero, which CVODE calls before this function, also clears the pattern.
    std::copy(integrator.sparseColumnStarts.begin(), integrator.sparseColumnStarts.end(),
              SUNSparseMatrix_IndexPointers(matrix));
    std::copy(integrator.sparseRowIndices.begin(), integrator.sparseRowIndices.end(),
              SUNSparseMatrix_IndexValues(matrix));
    sunrealtype* values = SUNSparseMatrix_Data(matrix);
    for (size_t nonzero = 0; nonzero < integrator.sparseBandIndex.size(); ++nonzero)
    {
      const size_t index = integrator.sparseBandIndex[nonzero];
      values[nonzero] = index == std::numeric_limits<size_t>::max() ? 0.0 : jacobian.band[index];
    }
    return 0;
  }
  catch (const std::runtime_error& error)
  {
    data->lastError = error.what();
    return 1;
  }
  catch (const std::exception& error)
  {
    data->lastError = error.what();
    return -1;
  }
}

int CVODE::solvePreconditioner(sunrealtype, N_Vector, N_Vector, N_Vector r, N_Vector z, sunrealtype, sunrealtype, int,
                               void* user_data)
{
  auto* data = static_cast<CVODE::UserData*>(user_data);
  CVODE& integrator = *data->integrator;
  const SorptionJacobian& jacobian = *integrator.sorptionJacobian;

  N_VScale(1.0, r, z);
  const ColumnStateLayout layout = data->column->stateLayout();
  double* base = N_VGetArrayPointer(z);
  jacobian.gather(layout.concentration(base), layout.physisorption(base), integrator.preconditionerVector);
  char transpose = 'N';
  int n = static_cast<int>(jacobian.size);
  int bands = static_cast<int>(jacobian.bandwidth);
  int leadingDimension = static_cast<int>(jacobian.leadingDimension);
  int rightHandSides = 1;
  int info = 0;
  dgbtrs_(&transpose, &n, &bands, &bands, &rightHandSides, integrator.preconditionerBand.data(), &leadingDimension,
          integrator.preconditionerPivots.data(), integrator.preconditionerVector.data(), &n, &info);
  if (info != 0) return -1;
  jacobian.scatter(integrator.preconditionerVector, layout.concentration(base), layout.physisorption(base));
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

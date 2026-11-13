// Optional SUNDIALS CVODE backend. CVODE performs the stiff integration;
// RUPTURA retains its existing accepted-step, event and output interfaces.
// API reference:
// https://sundials.readthedocs.io/en/v7.4.0/cvode/Usage/index.html
#include <algorithm>
#include <climits>
#include <cmath>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "cycle/integrator.hpp"

#ifndef RUPTURA_HAS_SUNDIALS
#define RUPTURA_HAS_SUNDIALS 0
#endif

#if RUPTURA_HAS_SUNDIALS
#include <cvode/cvode.h>
#include <cvode/cvode_bandpre.h>
#include <cvode/cvode_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#endif

namespace cycle
{

bool cvode_available() noexcept { return RUPTURA_HAS_SUNDIALS != 0; }

std::string cvode_version()
{
#if RUPTURA_HAS_SUNDIALS
  return SUNDIALS_VERSION;
#else
  return "unavailable (built without SUNDIALS)";
#endif
}

#if RUPTURA_HAS_SUNDIALS
namespace
{
constexpr double eps = std::numeric_limits<double>::epsilon();

bool finite(const OdeState& state)
{
  return std::all_of(state.begin(), state.end(), [](double x) { return std::isfinite(x); });
}

template <class Flag>
void checked(Flag code, const char* operation)
{
  if (static_cast<long>(code) < 0)
    throw std::runtime_error(std::string(operation) + " failed with SUNDIALS flag " +
                             std::to_string(static_cast<long>(code)));
}

struct Handles
{
  SUNContext context = nullptr;
  N_Vector state = nullptr, tolerance = nullptr, interpolated = nullptr, constraints = nullptr;
  SUNMatrix matrix = nullptr;
  SUNLinearSolver linear = nullptr;
  void* memory = nullptr;
  ~Handles()
  {
    if (memory) CVodeFree(&memory);
    if (linear) SUNLinSolFree(linear);
    if (matrix) SUNMatDestroy(matrix);
    if (constraints) N_VDestroy(constraints);
    if (interpolated) N_VDestroy(interpolated);
    if (tolerance) N_VDestroy(tolerance);
    if (state) N_VDestroy(state);
    if (context) SUNContext_Free(&context);
  }
  Handles() = default;
  Handles(const Handles&) = delete;
  Handles& operator=(const Handles&) = delete;
};

// All user operations called from a SUNDIALS C callback are guarded here.
// A thrown exception is saved and reported after CVode returns to C++.
struct CallbackData
{
  const OdeRhs& rhs;
  const BdfOptions& options;
  const OdeState& absolute;
  const std::size_t n, physical;
  std::size_t calls = 0;
  std::exception_ptr pending;
  std::string recoverable_error;
  OdeState state, derivative, factors;
  std::vector<std::vector<std::size_t>> rows, groups;

  CallbackData(const OdeRhs& fun, const BdfOptions& opts, const OdeState& atol, std::size_t dimension)
      : rhs(fun),
        options(opts),
        absolute(atol),
        n(dimension),
        physical(opts.physical_dimension ? opts.physical_dimension : dimension),
        state(n),
        derivative(n),
        factors(physical, std::sqrt(eps))
  {
    if (opts.jacobian_rows.empty())
    {
      rows.resize(physical);
      groups.resize(physical);
      for (std::size_t col = 0; col < physical; ++col)
      {
        rows[col].resize(n);
        std::iota(rows[col].begin(), rows[col].end(), 0);
        groups[col].push_back(col);
      }
      return;
    }
    rows = opts.jacobian_rows;
    if (rows.size() != physical) throw std::invalid_argument("jacobian_rows must have one entry per physical state");
    std::vector<std::vector<bool>> occupied;
    for (std::size_t col = 0; col < physical; ++col)
    {
      auto& rs = rows[col];
      std::sort(rs.begin(), rs.end());
      rs.erase(std::unique(rs.begin(), rs.end()), rs.end());
      if (!rs.empty() && rs.back() >= n) throw std::invalid_argument("Jacobian sparsity row out of range");
      std::size_t group = 0;
      for (; group < groups.size(); ++group)
      {
        bool conflict = false;
        for (auto row : rs)
          if (occupied[group][row])
          {
            conflict = true;
            break;
          }
        if (!conflict) break;
      }
      if (group == groups.size())
      {
        groups.push_back({});
        occupied.emplace_back(n, false);
      }
      groups[group].push_back(col);
      for (auto row : rs) occupied[group][row] = true;
    }
  }

  // Non-finite trial derivatives are recoverable; state/dimension errors and
  // user exceptions are unrecoverable. Neither case crosses the C boundary.
  bool evaluate(double time, const OdeState& y, OdeState& out)
  {
    if (options.maximum_rhs_evaluations && calls >= options.maximum_rhs_evaluations)
      throw std::runtime_error("CVODE maximum RHS evaluation count exceeded");
    ++calls;
    if (!finite(y))
    {
      recoverable_error = "ODE trial state is non-finite";
      return false;
    }
    out.assign(n, 0.0);
    try
    {
      rhs(time, y, out);
    }
    catch (const std::overflow_error& e)
    {
      recoverable_error = e.what();
      return false;
    }
    catch (const std::invalid_argument& e)
    {
      const std::string message = e.what();
      if (message.find("finite") != std::string::npos)
      {
        recoverable_error = message;
        return false;
      }
      throw;
    }
    if (out.size() != n) throw std::runtime_error("ODE RHS returned incorrect dimension");
    if (!finite(out))
    {
      recoverable_error = "ODE RHS produced non-finite trial derivatives";
      return false;
    }
    return true;
  }

  static int rhs_callback(sunrealtype time, N_Vector y, N_Vector f, void* user) noexcept
  {
    auto& self = *static_cast<CallbackData*>(user);
    try
    {
      const double* yp = N_VGetArrayPointer(y);
      double* fp = N_VGetArrayPointer(f);
      if (!yp || !fp) throw std::runtime_error("CVODE serial vector has no array storage");
      std::copy(yp, yp + self.n, self.state.begin());
      if (!self.evaluate(time, self.state, self.derivative)) return 1;
      std::copy(self.derivative.begin(), self.derivative.end(), fp);
      return 0;
    }
    catch (...)
    {
      self.pending = std::current_exception();
      return -1;
    }
  }

  static int jacobian_callback(sunrealtype time, N_Vector y, N_Vector fy, SUNMatrix matrix, void* user, N_Vector,
                               N_Vector, N_Vector) noexcept
  {
    auto& self = *static_cast<CallbackData*>(user);
    try
    {
      const double* yp = N_VGetArrayPointer(y);
      const double* fp = N_VGetArrayPointer(fy);
      if (!yp || !fp) throw std::runtime_error("CVODE Jacobian vector lacks array storage");
      const OdeState base(yp, yp + self.n), f(fp, fp + self.n);
      OdeState increment(self.physical), yscale(self.physical), perturbed, output;
      checked(SUNMatZero(matrix), "SUNMatZero");
      for (std::size_t col = 0; col < self.physical; ++col)
      {
        yscale[col] = (f[col] >= 0 ? 1.0 : -1.0) * std::max(self.absolute[col], std::abs(base[col]));
        increment[col] = (base[col] + self.factors[col] * yscale[col]) - base[col];
        while (increment[col] == 0.0)
        {
          self.factors[col] *= 10.0;
          if (!std::isfinite(self.factors[col]))
            throw std::runtime_error("Cannot select CVODE finite-difference increment");
          increment[col] = (base[col] + self.factors[col] * yscale[col]) - base[col];
        }
      }
      const double reject = std::pow(eps, .875), small = std::pow(eps, .75), big = std::pow(eps, .25);
      for (const auto& group : self.groups)
      {
        perturbed = base;
        for (auto col : group) perturbed[col] += increment[col];
        if (!self.evaluate(time, perturbed, output)) return 1;
        for (auto col : group)
        {
          if (self.rows[col].empty()) continue;
          double maximum = -1, scale = 0;
          for (auto row : self.rows[col])
          {
            double diff = output[row] - f[row];
            SM_ELEMENT_D(matrix, row, col) = diff;
            if (std::abs(diff) > maximum)
            {
              maximum = std::abs(diff);
              scale = std::max(std::abs(output[row]), std::abs(f[row]));
            }
          }
          if (maximum < reject * scale)
          {
            double factor = std::min(self.factors[col] * 10,1e-2);
            double step = (base[col] + factor * yscale[col]) - base[col];
            OdeState y2(base), f2;
            y2[col] += step;
            if (!self.evaluate(time, y2, f2)) return 1;
            double max2 = -1, scale2 = 0;
            for (auto row : self.rows[col])
              if (std::abs(f2[row] - f[row]) > max2)
              {
                max2 = std::abs(f2[row] - f[row]);
                scale2 = std::max(std::abs(f2[row]), std::abs(f[row]));
              }
            if (maximum * scale2 < max2 * scale)
            {
              self.factors[col] = factor;
              increment[col] = step;
              maximum = max2;
              scale = scale2;
              for (auto row : self.rows[col]) SM_ELEMENT_D(matrix, row, col) = f2[row] - f[row];
            }
          }
          for (auto row : self.rows[col])
          {
            double& value = SM_ELEMENT_D(matrix, row, col);
            value /= increment[col];
            if (!std::isfinite(value))
            {
              self.recoverable_error = "Non-finite CVODE numerical Jacobian";
              return 1;
            }
          }
          if (maximum < small * scale) self.factors[col] *= 10;
          if (maximum > big * scale) self.factors[col] *= .1;
          self.factors[col] = std::clamp(self.factors[col],1000*eps,1e-2);
        }
      }
      // Trailing ledger columns remain zero. CVODE still solves the complete
      // augmented Newton system, retaining physical-to-ledger coupling.
      return 0;
    }
    catch (...)
    {
      self.pending = std::current_exception();
      return -1;
    }
  }
};

void statistics(Handles& handles, const CallbackData& data, BdfStatistics& stats)
{
  long steps = 0, errors = 0, nonlinear_failures = 0, jacobians = 0, setups = 0;
  int order = 1;
  (void)CVodeGetNumSteps(handles.memory, &steps);
  (void)CVodeGetNumErrTestFails(handles.memory, &errors);
  (void)CVodeGetNumNonlinSolvConvFails(handles.memory, &nonlinear_failures);
  (void)CVodeGetNumJacEvals(handles.memory, &jacobians);
  (void)CVodeGetNumLinSolvSetups(handles.memory, &setups);
  (void)CVodeGetLastOrder(handles.memory, &order);
  stats.accepted_steps = static_cast<std::size_t>(std::max(steps, 0L));
  stats.rejected_steps = static_cast<std::size_t>(std::max(errors + nonlinear_failures, 0L));
  stats.rhs_evaluations = data.calls;  // includes custom finite-difference evaluations
  stats.jacobian_evaluations = static_cast<std::size_t>(std::max(jacobians, 0L));
  stats.lu_factorizations = static_cast<std::size_t>(std::max(setups, 0L));
  stats.largest_order = std::max(stats.largest_order, order);
}

void validate(const OdeRhs& rhs, double start, const OdeState& state, double end, const BdfOptions& options,
              const std::vector<OdeEvent>& events)
{
  if (!rhs) throw std::invalid_argument("ODE RHS is empty");
  if (!std::isfinite(start) || !std::isfinite(end) || end < start)
    throw std::invalid_argument("CVODE requires a finite forward time interval");
  if (state.empty() || !finite(state) || state.size() > static_cast<std::size_t>(LONG_MAX))
    throw std::invalid_argument("Invalid CVODE initial state");
  if (!(options.relative_tolerance > 0) || !std::isfinite(options.relative_tolerance))
    throw std::invalid_argument("Relative tolerance must be finite and positive");
  if (options.absolute_tolerance.size() != 1 && options.absolute_tolerance.size() != state.size())
    throw std::invalid_argument("Absolute tolerance must be scalar or state-sized");
  for (double value : options.absolute_tolerance)
    if (!(value > 0) || !std::isfinite(value))
      throw std::invalid_argument("Absolute tolerances must be finite and positive");
  if (!(options.maximum_step > 0) || std::isnan(options.maximum_step) || options.first_step < 0 ||
      !std::isfinite(options.first_step) || !options.maximum_steps || options.physical_dimension > state.size())
    throw std::invalid_argument("Invalid CVODE step/dimension limit");
  double previous = start;
  for (double time : options.sample_times)
  {
    if (!std::isfinite(time) || time < previous || time > end)
      throw std::invalid_argument("Sample times must be ordered within integration interval");
    previous = time;
  }
  for (const auto& event : events)
    if (!event.value) throw std::invalid_argument("Empty event function");
}
}  // namespace

BdfResult integrate_cvode(const OdeRhs& rhs, double start, const OdeState& initial, double end,
                          const BdfOptions& options, const OdeAcceptedCallback& accepted,
                          const std::vector<OdeEvent>& events, const OdeSampleCallback& sample)
{
  static_assert(std::is_same_v<sunrealtype, double>, "RUPTURA-Cycles requires double-precision SUNDIALS");
  validate(rhs, start, initial, end, options, events);
  if (!options.state_permutation.empty() || !options.constraint_offsets.empty())
  {
    const auto n = initial.size();
    auto permutation = options.state_permutation;
    if (permutation.empty())
    {
      permutation.resize(n);
      std::iota(permutation.begin(), permutation.end(), 0);
    }
    auto offsets = options.constraint_offsets;
    if (offsets.empty()) offsets.assign(n, 0.);
    if (offsets.size() != n || !finite(offsets))
      throw std::invalid_argument("CVODE constraint offsets have incompatible dimensions");
    if (permutation.size() != n) throw std::invalid_argument("CVODE state permutation has incorrect size");
    std::vector<std::size_t> inverse(n, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      auto j = permutation[i];
      if (j >= n || inverse[j] != n) throw std::invalid_argument("CVODE state permutation is not bijective");
      inverse[j] = i;
    }
    const auto physical = options.physical_dimension ? options.physical_dimension : n;
    for (std::size_t i = 0; i < n; ++i)
      if ((i < physical) != (permutation[i] < physical))
        throw std::invalid_argument("CVODE state permutation must retain the physical/ledger partition");
    auto pack = [&](const OdeState& source)
    {
      if (source.size() != n) throw std::runtime_error("ODE RHS returned incorrect dimension during permutation");
      OdeState out(n);
      for (std::size_t i = 0; i < n; ++i) out[i] = source[permutation[i]];
      return out;
    };
    auto unpack = [&](const OdeState& source)
    {
      OdeState out(n);
      for (std::size_t i = 0; i < n; ++i) out[permutation[i]] = source[i] - offsets[permutation[i]];
      return out;
    };
    auto packed_initial = pack(initial);
    for (std::size_t i = 0; i < n; ++i) packed_initial[i] += offsets[permutation[i]];
    BdfOptions local = options;
    local.constraint_offsets.clear();
    local.state_permutation.clear();
    if (local.absolute_tolerance.size() != 1) local.absolute_tolerance = pack(local.absolute_tolerance);
    for (auto& i : local.nonnegative_indices) i = inverse.at(i);
    if (!local.jacobian_rows.empty())
    {
      auto original = local.jacobian_rows;
      for (std::size_t i = 0; i < physical; ++i)
      {
        auto& rows = local.jacobian_rows[i];
        rows = original[permutation[i]];
        for (auto& row : rows) row = inverse.at(row);
      }
    }
    OdeRhs ordered = [&](double t, const OdeState& y, OdeState& dy)
    {
      auto public_state = unpack(y);
      OdeState derivative(n);
      rhs(t, public_state, derivative);
      dy = pack(derivative);
    };
    std::vector<OdeEvent> ordered_events;
    for (const auto& event : events)
      ordered_events.push_back({event.name, [&, event](double t, const OdeState& y)
                                { return event.value(t, unpack(y)); }, event.direction, event.terminal});
    OdeAcceptedCallback ordered_accepted;
    if (accepted)
      ordered_accepted = [&](double t, const OdeState& y, const BdfStatistics& stats)
      { return accepted(t, unpack(y), stats); };
    OdeSampleCallback ordered_sample;
    if (sample) ordered_sample = [&](double t, const OdeState& y) { sample(t, unpack(y)); };
    auto result =
        integrate_cvode(ordered, start, packed_initial, end, local, ordered_accepted, ordered_events, ordered_sample);
    result.state = unpack(result.state);
    for (auto& state : result.states) state = unpack(state);
    for (auto& event : result.events) event.state = unpack(event.state);
    result.statistics.linear_solver +=
        options.state_permutation.empty() ? "; precision-bounded state" : "; node-major precision-bounded state";
    return result;
  }
  const auto n = initial.size();
  OdeState state(initial), absolute(n), interpolated(n);
  for (std::size_t i = 0; i < n; ++i)
    absolute[i] = options.absolute_tolerance[options.absolute_tolerance.size() == 1 ? 0 : i];
  CallbackData data(rhs, options, absolute, n);
  Handles handles;
  BdfResult result;
  result.time = start;
  result.state = initial;
  result.statistics.linear_solver = options.cvode_spgmr
                                        ? "SUNDIALS CVODE BDF/Newton; matrix-free SPGMR with band preconditioner"
                                        : "SUNDIALS CVODE BDF/Newton; dense direct LU (full augmented state)";
  auto store = [&](double time, const OdeState& y)
  {
    if (!result.times.empty() && result.times.back() == time)
      result.states.back() = y;
    else
    {
      result.times.push_back(time);
      result.states.push_back(y);
    }
  };
  auto emit = [&](double time, const OdeState& y)
  {
    store(time, y);
    if (sample) sample(time, y);
  };
  store(start, initial);
  std::size_t sample_index = 0;
  try
  {
    while (sample_index < options.sample_times.size() && options.sample_times[sample_index] == start)
    {
      emit(start, initial);
      ++sample_index;
    }
  }
  catch (const std::exception& error)
  {
    result.message = error.what();
    return result;
  }
  catch (...)
  {
    result.message = "Unknown exception from initial CVODE output handler";
    return result;
  }
  if (start == end)
  {
    result.success = true;
    result.message = "Reached integration endpoint";
    return result;
  }
  try
  {
    checked(SUNContext_Create(SUN_COMM_NULL, &handles.context), "SUNContext_Create");
    const auto dimension = static_cast<sunindextype>(n);
    if (dimension <= 0 || static_cast<std::size_t>(dimension) != n)
      throw std::invalid_argument("CVODE state exceeds SUNDIALS index range");
    handles.state = N_VMake_Serial(dimension, state.data(), handles.context);
    handles.tolerance = N_VMake_Serial(dimension, absolute.data(), handles.context);
    handles.interpolated = N_VMake_Serial(dimension, interpolated.data(), handles.context);
    if (!handles.state || !handles.tolerance || !handles.interpolated)
      throw std::runtime_error("Failed to allocate CVODE serial vectors");
    handles.memory = CVodeCreate(CV_BDF, handles.context);
    if (!handles.memory) throw std::runtime_error("CVodeCreate failed");
    checked(CVodeInit(handles.memory, CallbackData::rhs_callback, start, handles.state), "CVodeInit");
    checked(CVodeSVtolerances(handles.memory, std::max(options.relative_tolerance, 100 * eps), handles.tolerance),
            "CVodeSVtolerances");
    checked(CVodeSetUserData(handles.memory, &data), "CVodeSetUserData");
    if (options.apply_nonnegative_constraints && !options.nonnegative_indices.empty())
    {
      handles.constraints = N_VClone(handles.state);
      if (!handles.constraints) throw std::runtime_error("Cannot allocate CVODE constraints");
      N_VConst(0., handles.constraints);
      auto* constraints = N_VGetArrayPointer(handles.constraints);
      for (auto i : options.nonnegative_indices)
      {
        if (i >= n) throw std::invalid_argument("Nonnegative state index out of range");
        constraints[i] = 1.;
      }
      checked(CVodeSetConstraints(handles.memory, handles.constraints), "CVodeSetConstraints");
    }
    if (options.cvode_spgmr)
    {
      if (options.cvode_krylov_dimension < 1) throw std::invalid_argument("CVODEKrylovDimension must be positive");
      handles.linear = SUNLinSol_SPGMR(handles.state, SUN_PREC_RIGHT, options.cvode_krylov_dimension, handles.context);
      if (!handles.linear) throw std::runtime_error("Cannot create CVODE SPGMR solver");
    }
    else
    {
      handles.matrix = SUNDenseMatrix(dimension, dimension, handles.context);
      if (!handles.matrix) throw std::runtime_error("Failed to allocate CVODE dense matrix");
      handles.linear = SUNLinSol_Dense(handles.state, handles.matrix, handles.context);
      if (!handles.linear) throw std::runtime_error("Failed to allocate CVODE dense linear solver");
    }
    if (options.cvode_maximum_order < 1 || options.cvode_maximum_order > 5)
      throw std::invalid_argument("CVODEMaximumOrder must lie between 1 and 5");
    checked(CVodeSetMaxOrd(handles.memory, options.cvode_maximum_order), "CVodeSetMaxOrd");
    checked(CVodeSetLinearSolver(handles.memory, handles.linear, handles.matrix), "CVodeSetLinearSolver");
    if (options.cvode_spgmr)
    {
      auto bandwidth = static_cast<sunindextype>(std::min(n - 1, options.cvode_preconditioner_bandwidth));
      checked(CVBandPrecInit(handles.memory, dimension, bandwidth, bandwidth), "CVBandPrecInit");
    }
    if (!options.cvode_spgmr) checked(CVodeSetJacFn(handles.memory, CallbackData::jacobian_callback), "CVodeSetJacFn");
    checked(CVodeSetMaxNumSteps(handles.memory,
                                static_cast<long>(std::min(options.maximum_steps, static_cast<std::size_t>(LONG_MAX)))),
            "CVodeSetMaxNumSteps");
    checked(CVodeSetMaxStep(handles.memory, std::isfinite(options.maximum_step) ? options.maximum_step : 0.0),
            "CVodeSetMaxStep");
    if (options.first_step > 0)
      checked(CVodeSetInitStep(handles.memory, std::min({options.first_step, options.maximum_step, end - start})),
              "CVodeSetInitStep");
    // Every phase/checkpoint segment is a hard boundary: no RHS evaluation
    // beyond it is permitted, even when requesting one internal step at a time.
    checked(CVodeSetStopTime(handles.memory, end), "CVodeSetStopTime");
    checked(CVodeSetInterpolateStopTime(handles.memory, SUNFALSE), "CVodeSetInterpolateStopTime");
    std::vector<double> previous(events.size());
    for (std::size_t i = 0; i < events.size(); ++i)
    {
      previous[i] = events[i].value(start, initial);
      if (!std::isfinite(previous[i])) throw std::runtime_error("Non-finite initial event value");
    }
    double interpolation_left_time = start, interpolation_right_time = start;
    OdeState interpolation_left = initial;
    auto interpolate = [&](double time)
    {
      checked(CVodeGetDky(handles.memory, time, 0, handles.interpolated), "CVodeGetDky");
      if (!finite(interpolated)) throw std::runtime_error("CVODE dense output is non-finite");
      bool outside = false;
      for (auto i : options.nonnegative_indices) outside |= interpolated[i] < 0.;
      if (outside && interpolation_right_time > interpolation_left_time)
      {
        // CVODE constrains accepted states, but not its dense polynomial.
        // Use a convex interpolation of the two accepted states rather than
        // clipping sampled states or returning an invalid event restart.
        const double w =
            std::clamp((time - interpolation_left_time) / (interpolation_right_time - interpolation_left_time), 0., 1.);
        for (std::size_t i = 0; i < n; ++i) interpolated[i] = (1. - w) * interpolation_left[i] + w * state[i];
      }
      return interpolated;
    };
    double time = start;
    while (time < end)
    {
      if (result.statistics.accepted_steps >= options.maximum_steps)
        throw std::runtime_error("CVODE maximum accepted-step count exceeded");
      const double old_time = time;
      interpolation_left_time = time;
      interpolation_left = state;
      sunrealtype reached = time;
      const int flag = CVode(handles.memory, end, handles.state, &reached, CV_ONE_STEP);
      statistics(handles, data, result.statistics);
      if (data.pending) std::rethrow_exception(std::exchange(data.pending, {}));
      if (flag < 0)
      {
        std::string message = "CVode failed with SUNDIALS flag " + std::to_string(flag);
        if (!data.recoverable_error.empty()) message += ": " + data.recoverable_error;
        throw std::runtime_error(message);
      }
      if (flag != CV_SUCCESS && flag != CV_TSTOP_RETURN)
        throw std::runtime_error("Unexpected CVODE return flag " + std::to_string(flag));
      time = reached;
      interpolation_right_time = time;
      if (!(time > old_time) || time > end || !finite(state))
        throw std::runtime_error("CVODE returned a non-finite state or invalid accepted time");
      result.time = time;
      result.state = state;

      // Use the same initial-zero, direction and nonterminal-event semantics
      // as native BDF. CVODE supplies the accepted polynomial; its built-in
      // rootfinder deliberately ignores initial roots and is not used here.
      double stop = time;
      std::vector<OdeEventHit> hits;
      for (std::size_t i = 0; i < events.size(); ++i)
      {
        double g0 = previous[i], g1 = events[i].value(time, state);
        if (!std::isfinite(g1)) throw std::runtime_error("Non-finite event value");
        previous[i] = g1;
        bool up = g0 <= 0 && g1 >= 0, down = g0 >= 0 && g1 <= 0;
        bool crossing = (events[i].direction == 0 && (up || down)) || (events[i].direction > 0 && up) ||
                        (events[i].direction < 0 && down);
        if (!crossing) continue;
        double left = old_time, right = time, gl = g0;
        if (g0 == 0)
          right = left;
        else if (g1 == 0)
          left = right;
        else
          for (int k = 0; k < 100; ++k)
          {
            if (right - left <= 4 * eps * std::max({1.0, std::abs(left), std::abs(right)})) break;
            double middle = left + .5 * (right - left);
            double gm = events[i].value(middle, interpolate(middle));
            if (!std::isfinite(gm)) throw std::runtime_error("Non-finite event root value");
            if (gm == 0)
            {
              left = right = middle;
              break;
            }
            if ((gl <= 0 && gm >= 0) || (gl >= 0 && gm <= 0))
              right = middle;
            else
            {
              left = middle;
              gl = gm;
            }
          }
        double root = .5 * (left + right);
        hits.push_back({i, events[i].name, root, interpolate(root)});
        if (events[i].terminal) stop = std::min(stop, root);
      }
      std::sort(hits.begin(), hits.end(), [](const auto& a, const auto& b) { return a.time < b.time; });
      for (const auto& hit : hits)
        if (hit.time <= stop)
        {
          result.events.push_back(hit);
          if (events[hit.index].terminal) result.stopped_by_event = true;
        }
      while (sample_index < options.sample_times.size() && options.sample_times[sample_index] <= stop)
      {
        double requested = options.sample_times[sample_index++];
        emit(requested, interpolate(requested));
      }
      if (result.stopped_by_event)
      {
        result.time = stop;
        result.state = interpolate(stop);
      }
      if (options.store_accepted_steps) store(result.time, result.state);
      if (accepted && !accepted(result.time, result.state, result.statistics)) result.stopped_by_callback = true;
      if (result.stopped_by_event || result.stopped_by_callback) break;
    }
    result.success = true;
    result.message = result.stopped_by_event      ? "Stopped at terminal event"
                     : result.stopped_by_callback ? "Stopped by accepted-step callback"
                                                  : "Reached integration endpoint";
  }
  catch (const std::exception& error)
  {
    result.success = false;
    result.message = error.what();
    if (handles.memory) statistics(handles, data, result.statistics);
  }
  catch (...)
  {
    result.success = false;
    result.message = "Unknown exception from CVODE callback or output handler";
    if (handles.memory) statistics(handles, data, result.statistics);
  }
  store(result.time, result.state);
  return result;
}
#else
BdfResult integrate_cvode(const OdeRhs&, double, const OdeState&, double, const BdfOptions&, const OdeAcceptedCallback&,
                          const std::vector<OdeEvent>&, const OdeSampleCallback&)
{
  throw std::runtime_error(
      "CVODE was requested, but this build has no SUNDIALS support. "
      "Rebuild with RUPTURA_CVODE=ON and SUNDIALS 7, or select BDF "
      "explicitly.");
}
#endif
}  // namespace cycle

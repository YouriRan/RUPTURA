#pragma once
#include <cstddef>
#include <functional>
#include <limits>
#include <string>
#include <vector>
namespace cycle
{
using OdeState = std::vector<double>;
using OdeRhs = std::function<void(double, const OdeState&, OdeState&)>;
struct BdfOptions
{
  double relative_tolerance = 1e-6;
  OdeState absolute_tolerance{1e-9};
  double maximum_step = std::numeric_limits<double>::infinity(), first_step = 0;
  std::size_t maximum_steps = 1000000, maximum_rhs_evaluations = 0, physical_dimension = 0;
  // Rows per physical-state column; include every possible nonzero, including
  // passive ledger rows. Trailing states must not affect RHS when eliminated.
  std::vector<std::vector<std::size_t>> jacobian_rows;
  bool use_lapack_if_available = true;
  bool cvode_spgmr = false;
  bool apply_nonnegative_constraints = true;
  std::size_t cvode_preconditioner_bandwidth = 0;
  // Solver index -> public state index. Physical states remain before ledgers.
  std::vector<std::size_t> state_permutation;
  // Internal z=y+offset permits a numerical lower bound without state clipping.
  OdeState constraint_offsets;
  std::vector<std::size_t> nonnegative_indices;
  int cvode_krylov_dimension = 30, cvode_maximum_order = 5;
  std::vector<double> sample_times;
  bool store_accepted_steps = false;
  int maximum_order = 5;
  bool ndf = true, adaptive_rk = true;
  std::size_t initialization_steps = 0;
  double event_time_tolerance = 1e-6;
};
struct OdeEvent
{
  std::string name;
  std::function<double(double, const OdeState&)> value;
  int direction = 0;
  bool terminal = true;
};
struct OdeEventHit
{
  std::size_t index = 0;
  std::string name;
  double time = 0;
  OdeState state;
};
struct BdfStatistics
{
  std::size_t accepted_steps = 0, rejected_steps = 0, rhs_evaluations = 0, jacobian_evaluations = 0,
              lu_factorizations = 0;
  int largest_order = 1;
  std::string linear_solver;
};
struct BdfResult
{
  bool success = false, stopped_by_event = false, stopped_by_callback = false;
  std::string message;
  double time = 0;
  OdeState state;
  std::vector<double> times;
  std::vector<OdeState> states;
  std::vector<OdeEventHit> events;
  BdfStatistics statistics;
};
// Return true to continue, false to stop cleanly at the accepted state.
using OdeAcceptedCallback = std::function<bool(double, const OdeState&, const BdfStatistics&)>;
using OdeSampleCallback = std::function<void(double, const OdeState&)>;
// Native forward-time variable-order BDF/NDF (1--5). No Python/CVODE runtime.
// Initial/final states are retained. Requested samples use dense interpolation.
// Invalid input throws; numerical failure returns the last accepted state.
BdfResult integrate_bdf(const OdeRhs&, double, const OdeState&, double, const BdfOptions& = {},
                        const OdeAcceptedCallback& = {}, const std::vector<OdeEvent>& = {},
                        const OdeSampleCallback& = {});
// Optional SUNDIALS CVODE BDF/Newton backend. Uses the same state/tolerance,
// accepted-step, dense-sample and terminal-event interface. There is no silent
// fallback: requesting CVODE in a build without it throws a clear exception.
bool cvode_available() noexcept;
std::string cvode_version();
BdfResult integrate_cvode(const OdeRhs&, double, const OdeState&, double, const BdfOptions& = {},
                          const OdeAcceptedCallback& = {}, const std::vector<OdeEvent>& = {},
                          const OdeSampleCallback& = {});
// Backward Euler uses the native Newton/error-control engine at order one,
// with the NDF modification disabled.
BdfResult integrate_bdf1(const OdeRhs&, double, const OdeState&, double, const BdfOptions& = {},
                         const OdeAcceptedCallback& = {}, const std::vector<OdeEvent>& = {},
                         const OdeSampleCallback& = {});
// Third-order SSP explicit RK and Pareschi--Russo SSP3(4,3,3) additive RK.
// The latter accepts the explicit part separately; implicit = full - explicit.
BdfResult integrate_ssprk3(const OdeRhs&, double, const OdeState&, double, const BdfOptions& = {},
                           const OdeAcceptedCallback& = {}, const std::vector<OdeEvent>& = {},
                           const OdeSampleCallback& = {});
BdfResult integrate_semi_implicit_ssprk3(const OdeRhs&, const OdeRhs& explicit_rhs,
    double, const OdeState&, double, const BdfOptions& = {}, const OdeAcceptedCallback& = {},
    const std::vector<OdeEvent>& = {}, const OdeSampleCallback& = {});
}  // namespace cycle

#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <type_traits>
#include <utility>

/**
 * \brief Accumulates wall-clock timings for major simulation phases.
 */
struct Timing
{
  /**
   * \brief Monotonic clock used for timing scopes.
   */
  using Clock = std::chrono::steady_clock;

  double total = 0.0;                       ///< Total elapsed time in s.
  double computePressure = 0.0;             ///< Time spent computing pressure fields in s.
  double computeEquilibriumLoadings = 0.0;  ///< Time spent computing equilibrium loadings in s.
  double computeVelocity = 0.0;             ///< Time spent computing velocity fields in s.
  double computeDerivatives = 0.0;          ///< Time spent computing derivatives in s.

  /**
   * \brief RAII timer that adds elapsed time to a referenced accumulator.
   */
  struct Scope
  {
    double& value;            ///< Accumulated elapsed time in s.
    Clock::time_point start;  ///< Start time for this scope.
    bool active = true;       ///< Indicates whether this scope still owns an active timing interval.

    /**
     * \brief Starts timing into the supplied accumulator.
     */
    explicit Scope(double& value) : value(value), start(Clock::now()) {}

    /**
     * \brief Copy construction is disabled to avoid double-counting.
     */
    Scope(const Scope&) = delete;

    /**
     * \brief Copy assignment is disabled to avoid double-counting.
     */
    Scope& operator=(const Scope&) = delete;

    /**
     * \brief Moves an active timing scope without stopping it.
     */
    Scope(Scope&& other) noexcept : value(other.value), start(other.start), active(other.active)
    {
      other.active = false;
    }

    /**
     * \brief Move assignment is disabled because the accumulator reference cannot be rebound.
     */
    Scope& operator=(Scope&&) = delete;

    /**
     * \brief Stops the timer and adds elapsed time to the accumulator.
     */
    ~Scope() { stop(); }

    /**
     * \brief Stops timing if this scope is still active.
     */
    void stop()
    {
      if (!active) return;
      value += std::chrono::duration<double>(Clock::now() - start).count();
      active = false;
    }
  };

  /**
   * \brief Resets all accumulated timing values to zero.
   */
  void clear()
  {
    total = 0.0;
    computePressure = 0.0;
    computeEquilibriumLoadings = 0.0;
    computeVelocity = 0.0;
    computeDerivatives = 0.0;
  }

  /**
   * \brief Creates an RAII scope that accumulates elapsed time into the supplied value.
   */
  Scope scoped(double& value) { return Scope(value); }

  /**
   * \brief Runs a callable while accumulating its elapsed time.
   */
  template <class Callable>
  decltype(auto) measure(double& value, Callable&& callable)
  {
    auto timer = scoped(value);
    using Result = std::invoke_result_t<Callable&>;
    if constexpr (std::is_void_v<Result>)
    {
      std::forward<Callable>(callable)();
      return;
    }
    else
    {
      return std::forward<Callable>(callable)();
    }
  }

  /**
   * \brief Prints accumulated timings to the provided C file handle.
   */
  void print(std::FILE* output = stdout) const
  {
    std::print(output, "Total time:                    {:.6f} s\n", total);
    std::print(output, "  computePressure:             {:.6f} s\n", computePressure);
    std::print(output, "  computeEquilibriumLoadings:  {:.6f} s\n", computeEquilibriumLoadings);
    std::print(output, "  computeVelocity:             {:.6f} s\n", computeVelocity);
    std::print(output, "  computeDerivatives:          {:.6f} s\n", computeDerivatives);
  }
};

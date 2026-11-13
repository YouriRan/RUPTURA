// Native C++ adaptation of scipy.integrate._ivp.bdf/common, SciPy 1.17.0.
// Difference-history transformations, order/step control, Newton stopping and
// adaptive finite differences are derived from SciPy (BSD-3-Clause).
// Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
// All rights reserved. Complete redistribution terms and disclaimer are in
// third_party/SCIPY_LICENSE.txt, which must accompany binary distributions.
#include "cycle/integrator.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#if defined(__unix__) || defined(__APPLE__)
#include <dlfcn.h>
#endif
namespace cycle {
namespace {
constexpr double eps = std::numeric_limits<double>::epsilon();
using V = OdeState;
using History = std::array<V, 8>;
bool finite(const V& v) {
  return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
}
double norm(const V& v, const V& s) {
  double max = 0, sum = 0;
  for (std::size_t i = 0; i < v.size(); ++i) {
    double a = std::abs(v[i] / s[i]);
    if (!std::isfinite(a)) return INFINITY;
    if (!a) continue;
    if (a > max) {
      sum = 1 + sum * (max / a) * (max / a);
      max = a;
    } else
      sum += (a / max) * (a / max);
  }
  return max ? max * std::sqrt(sum / v.size()) : 0;
}
void change_D(History& d, int order, double factor) {
  double r[6][6]{}, u[6][6]{}, p[6][6]{};
  for (int j = 0; j <= order; ++j) {
    r[0][j] = u[0][j] = 1;
    for (int i = 1; i <= order; ++i) {
      r[i][j] = r[i - 1][j] * (i - 1.0 - factor * j) / i;
      u[i][j] = u[i - 1][j] * (i - 1.0 - j) / i;
    }
  }
  for (int i = 0; i <= order; ++i)
    for (int j = 0; j <= order; ++j)
      for (int k = 0; k <= order; ++k) p[i][j] += r[i][k] * u[k][j];
  for (std::size_t k = 0; k < d[0].size(); ++k) {
    double next[6]{};
    for (int i = 0; i <= order; ++i)
      for (int j = 0; j <= order; ++j) next[i] += p[j][i] * d[j][k];
    for (int i = 0; i <= order; ++i) d[i][k] = next[i];
  }
}
struct Lapack {
  using Getrf = void (*)(const int*, const int*, double*, const int*, int*, int*);
  using Getrs = void (*)(const char*, const int*, const int*, const double*, const int*, const int*,
                         double*, const int*, int*);
  Getrf getrf = nullptr;
  Getrs getrs = nullptr;
  void* handle = nullptr;
  Lapack() {
#if defined(__unix__) || defined(__APPLE__)
    // Conventional LP64 symbols only; ILP64 LAPACK has a different integer ABI.
    for (const char* name : {"liblapack.so.3", "liblapack.so", "liblapack.dylib"}) {
      handle = dlopen(name, RTLD_LAZY | RTLD_LOCAL);
      if (!handle) continue;
      void* a = dlsym(handle, "dgetrf_");
      void* b = dlsym(handle, "dgetrs_");
      if (a && b) {
        static_assert(sizeof(a) == sizeof(getrf));
        std::memcpy(&getrf, &a, sizeof(a));
        std::memcpy(&getrs, &b, sizeof(b));
        break;
      }
      dlclose(handle);
      handle = nullptr;
    }
#endif
  }
  ~Lapack() {
#if defined(__unix__) || defined(__APPLE__)
    if (handle) dlclose(handle);
#endif
  }
};
struct Linear {
  std::size_t n, m;
  const Lapack* lapack;
  V lu, ledger;
  std::vector<int> piv;
  double c = 0;
  bool valid = false;
  Linear(std::size_t n_, std::size_t m_, const Lapack* l)
      : n(n_), m(m_), lapack(l), lu(m * m), ledger((n - m) * m), piv(m) {}
  bool factor(const V& jac, double c_) {
    c = c_;
    for (std::size_t j = 0; j < m; ++j) {
      for (std::size_t i = 0; i < m; ++i) lu[j * m + i] = (i == j ? 1 : 0) - c * jac[j * n + i];
      for (std::size_t i = m; i < n; ++i) ledger[j * (n - m) + i - m] = jac[j * n + i];
    }
    if (lapack) {
      int dim = static_cast<int>(m), info = 0;
      lapack->getrf(&dim, &dim, lu.data(), &dim, piv.data(), &info);
      return valid = info == 0 && finite(lu);
    }
    for (std::size_t k = 0; k < m; ++k) {
      std::size_t p = k;
      for (std::size_t i = k + 1; i < m; ++i)
        if (std::abs(lu[k * m + i]) > std::abs(lu[k * m + p])) p = i;
      piv[k] = static_cast<int>(p);
      if (!std::isfinite(lu[k * m + p]) || lu[k * m + p] == 0) return valid = false;
      if (p != k)
        for (std::size_t j = 0; j < m; ++j) std::swap(lu[j * m + k], lu[j * m + p]);
      for (std::size_t i = k + 1; i < m; ++i) lu[k * m + i] /= lu[k * m + k];
      for (std::size_t j = k + 1; j < m; ++j) {
        double a = lu[j * m + k];
        if (!a) continue;
        for (std::size_t i = k + 1; i < m; ++i) lu[j * m + i] -= lu[k * m + i] * a;
      }
    }
    return valid = finite(lu);
  }
  bool solve(V& b) const {
    if (!valid) return false;
    if (lapack) {
      int dim = static_cast<int>(m), one = 1, info = 0;
      const char trans = 'N';
      lapack->getrs(&trans, &dim, &one, lu.data(), &dim, piv.data(), b.data(), &dim, &info);
      if (info) return false;
    } else {
      for (std::size_t k = 0; k < m; ++k)
        if (static_cast<std::size_t>(piv[k]) != k) std::swap(b[k], b[piv[k]]);
      for (std::size_t k = 0; k < m; ++k)
        for (std::size_t i = k + 1; i < m; ++i) b[i] -= lu[k * m + i] * b[k];
      for (std::size_t k = m; k-- > 0;) {
        b[k] /= lu[k * m + k];
        for (std::size_t i = 0; i < k; ++i) b[i] -= lu[k * m + i] * b[k];
      }
    }
    for (std::size_t j = 0; j < m; ++j)
      for (std::size_t i = m; i < n; ++i) b[i] += c * ledger[j * (n - m) + i - m] * b[j];
    return finite(b);
  }
};
class Bdf {
  const OdeRhs& rhs;
  const BdfOptions& o;
  BdfResult result;
  std::size_t n, m;
  V atol, factors, jac;
  std::vector<std::vector<std::size_t>> rows, groups;
  History d;
  std::array<double, 6> gamma{}, alpha{}, ec{};
  double rtol, newton_tol, t, bound, h = 0;
  int order = 1, equal_steps = 0;
  std::string trial_error;
  Linear lu;
  void evaluate(double time, const V& y, V& f) {
    if (o.maximum_rhs_evaluations && result.statistics.rhs_evaluations >= o.maximum_rhs_evaluations)
      throw std::runtime_error("BDF maximum RHS evaluation count exceeded");
    ++result.statistics.rhs_evaluations;
    f.assign(n, 0);
    rhs(time, y, f);
    if (f.size() != n) throw std::runtime_error("ODE RHS returned incorrect dimension");
    if (!finite(f)) throw std::runtime_error("ODE RHS produced non-finite values");
  }
  bool trial(double time, const V& y, V& f) {
    try {
      evaluate(time, y, f);
      return true;
    } catch (const std::exception& e) {
      trial_error = e.what();
      if (o.maximum_rhs_evaluations &&
          result.statistics.rhs_evaluations >= o.maximum_rhs_evaluations)
        throw;
      return false;
    }
  }
  V scale(const V& y) const {
    V s(n);
    for (std::size_t i = 0; i < n; ++i) s[i] = atol[i] + rtol * std::abs(y[i]);
    return s;
  }
  void coloring() {
    if (o.jacobian_rows.empty()) {
      rows.resize(m);
      groups.resize(m);
      for (std::size_t j = 0; j < m; ++j) {
        rows[j].resize(n);
        std::iota(rows[j].begin(), rows[j].end(), 0);
        groups[j].push_back(j);
      }
      return;
    }
    rows = o.jacobian_rows;
    if (rows.size() != m)
      throw std::invalid_argument("jacobian_rows must have one entry per physical state");
    std::vector<std::vector<bool>> used;
    for (std::size_t j = 0; j < m; ++j) {
      auto& rs = rows[j];
      std::sort(rs.begin(), rs.end());
      rs.erase(std::unique(rs.begin(), rs.end()), rs.end());
      if (!rs.empty() && rs.back() >= n)
        throw std::invalid_argument("Jacobian sparsity row out of range");
      std::size_t g = 0;
      for (; g < groups.size(); ++g) {
        bool conflict = false;
        for (auto i : rs)
          if (used[g][i]) {
            conflict = true;
            break;
          }
        if (!conflict) break;
      }
      if (g == groups.size()) {
        groups.push_back({});
        used.emplace_back(n, false);
      }
      groups[g].push_back(j);
      for (auto i : rs) used[g][i] = true;
    }
  }
  void jacobian(double time, const V& y, const V* known = nullptr) {
    ++result.statistics.jacobian_evaluations;
    V f, yp, fp, inc(m), ys(m);
    if (known)
      f = *known;
    else
      evaluate(time, y, f);
    std::fill(jac.begin(), jac.end(), 0);
    for (std::size_t j = 0; j < m; ++j) {
      ys[j] = (f[j] >= 0 ? 1 : -1) * std::max(atol[j], std::abs(y[j]));
      inc[j] = (y[j] + factors[j] * ys[j]) - y[j];
      while (inc[j] == 0) {
        factors[j] *= 10;
        inc[j] = (y[j] + factors[j] * ys[j]) - y[j];
        if (!std::isfinite(factors[j]))
          throw std::runtime_error("Cannot select finite-difference increment");
      }
    }
    const double reject = std::pow(eps, .875), small = std::pow(eps, .75), big = std::pow(eps, .25);
    for (const auto& group : groups) {
      yp = y;
      for (auto j : group) yp[j] += inc[j];
      evaluate(time, yp, fp);
      for (auto j : group) {
        if (rows[j].empty()) continue;
        double md = -1, fs = 0;
        for (auto i : rows[j]) {
          double df = fp[i] - f[i];
          jac[j * n + i] = df;
          if (std::abs(df) > md) {
            md = std::abs(df);
            fs = std::max(std::abs(fp[i]), std::abs(f[i]));
          }
        }
        if (md < reject * fs) {
          double factor = std::min(factors[j] * 10,1e-2), hi = (y[j] + factor * ys[j]) - y[j];
          V y2(y), f2;
          y2[j] += hi;
          evaluate(time, y2, f2);
          double md2 = -1, fs2 = 0;
          for (auto i : rows[j])
            if (std::abs(f2[i] - f[i]) > md2) {
              md2 = std::abs(f2[i] - f[i]);
              fs2 = std::max(std::abs(f2[i]), std::abs(f[i]));
            }
          if (md * fs2 < md2 * fs) {
            factors[j] = factor;
            inc[j] = hi;
            md = md2;
            fs = fs2;
            for (auto i : rows[j]) jac[j * n + i] = f2[i] - f[i];
          }
        }
        for (auto i : rows[j]) jac[j * n + i] /= inc[j];
        if (md < small * fs) factors[j] *= 10;
        if (md > big * fs) factors[j] *= .1;
        factors[j] = std::clamp(factors[j],1000*eps,1e-2);
      }
    }
    if (!finite(jac)) throw std::runtime_error("Non-finite BDF numerical Jacobian");
  }
  double initial_step(const V& y, const V& f) {
    auto s = scale(y);
    double d0 = norm(y, s), d1 = norm(f, s), h0 = d0 < 1e-5 || d1 < 1e-5 ? 1e-6 : .01 * d0 / d1;
    h0 = std::min(h0, bound - t);
    V y1(n), f1, df(n);
    for (int k = 0; k < 30; ++k) {
      for (std::size_t i = 0; i < n; ++i) y1[i] = y[i] + h0 * f[i];
      if (trial(t + h0, y1, f1)) break;
      h0 *= .1;
      if (k == 29) throw std::runtime_error("Unable to evaluate initial BDF probe: " + trial_error);
    }
    for (std::size_t i = 0; i < n; ++i) df[i] = f1[i] - f[i];
    double d2 = norm(df, s) / h0, h1 = d1 <= 1e-15 && d2 <= 1e-15
                                           ? std::max(1e-6, h0 * .001)
                                           : std::sqrt(.01 / std::max(d1, d2));
    return std::min({100 * h0, h1, bound - t, o.maximum_step});
  }
  bool newton(double time, const V& prediction, double c, const V& psi, const V& s, V& next,
              V& correction, int& iterations) {
    next = prediction;
    correction.assign(n, 0);
    V f, delta(n);
    double old = -1;
    for (int k = 0; k < 4; ++k) {
      iterations = k + 1;
      if (!trial(time, next, f)) return false;
      for (std::size_t i = 0; i < n; ++i) delta[i] = c * f[i] - psi[i] - correction[i];
      if (!lu.solve(delta)) return false;
      double dn = norm(delta, s), rate = old < 0 ? -1 : dn / old;
      if (rate >= 0 && (rate >= 1 || std::pow(rate, 4 - k) / (1 - rate) * dn > newton_tol))
        return false;
      for (std::size_t i = 0; i < n; ++i) {
        next[i] += delta[i];
        correction[i] += delta[i];
      }
      if (dn == 0 || (rate >= 0 && rate / (1 - rate) * dn < newton_tol)) return finite(next);
      old = dn;
    }
    return false;
  }
  void step() {
    double min_h = 10 * std::abs(std::nextafter(t, INFINITY) - t);
    if (h > o.maximum_step) {
      change_D(d, order, o.maximum_step / h);
      h = o.maximum_step;
      equal_steps = 0;
    } else if (h < min_h) {
      change_D(d, order, min_h / h);
      h = min_h;
      equal_steps = 0;
    }
    bool current_jac = false;
    V predict(n), psi(n), next, correction, s;
    double time = t, error = 0, safety = .9;
    int iterations = 0;
    for (std::size_t attempt = 0;; ++attempt) {
      if (h < min_h || h <= 0 || !std::isfinite(h) || attempt > 1000)
        throw std::runtime_error("BDF step too small or repeatedly rejected" +
                                 (trial_error.empty() ? std::string() : ": " + trial_error));
      time = t + h;
      if (time > bound) {
        time = bound;
        change_D(d, order, (time - t) / h);
        equal_steps = 0;
        lu.valid = false;
      }
      h = time - t;
      for (std::size_t i = 0; i < n; ++i) {
        predict[i] = psi[i] = 0;
        for (int j = 0; j <= order; ++j) predict[i] += d[j][i];
        for (int j = 1; j <= order; ++j) psi[i] += d[j][i] * gamma[j] / alpha[order];
      }
      s = scale(predict);
      double c = h / alpha[order];
      bool converged = false;
      for (;;) {
        if (!lu.valid) {
          ++result.statistics.lu_factorizations;
          lu.factor(jac, c);
        }
        converged = lu.valid && newton(time, predict, c, psi, s, next, correction, iterations);
        if (converged || current_jac) break;
        jacobian(time, predict);
        lu.valid = false;
        current_jac = true;
      }
      if (!converged) {
        h *= .5;
        change_D(d, order, .5);
        equal_steps = 0;
        lu.valid = false;
        ++result.statistics.rejected_steps;
        continue;
      }
      safety = .9 * 9 / (8 + iterations);
      s = scale(next);
      error = ec[order] * norm(correction, s);
      if (!std::isfinite(error) || error > 1) {
        double factor = std::max(.2, safety * std::pow(error, -1.0 / (order + 1)));
        h *= factor;
        change_D(d, order, factor);
        equal_steps = 0;
        ++result.statistics.rejected_steps;
        continue;
      }
      break;
    }
    ++equal_steps;
    ++result.statistics.accepted_steps;
    t = time;
    result.time = t;
    result.state = next;
    for (std::size_t i = 0; i < n; ++i) {
      d[order + 2][i] = correction[i] - d[order + 1][i];
      d[order + 1][i] = correction[i];
      for (int j = order; j >= 0; --j) d[j][i] += d[j + 1][i];
    }
    if (equal_steps < order + 1) return;
    double lower = order > 1 ? ec[order - 1] * norm(d[order], s) : INFINITY,
           upper = order < o.maximum_order ? ec[order + 1] * norm(d[order + 2], s) : INFINITY;
    double f[] = {std::pow(lower, -1.0 / order), std::pow(error, -1.0 / (order + 1)),
                  std::pow(upper, -1.0 / (order + 2))};
    int best = 0;
    for (int j = 1; j < 3; ++j)
      if (f[j] > f[best]) best = j;
    order += best - 1;
    result.statistics.largest_order = std::max(result.statistics.largest_order, order);
    double factor = std::min(10.0, safety * f[best]);
    h *= factor;
    change_D(d, order, factor);
    equal_steps = 0;
    lu.valid = false;
  }
  V interpolate(double time) const {
    V y(d[0]);
    double p = 1;
    for (int j = 0; j < order; ++j) {
      p *= (time - (t - h * j)) / (h * (j + 1));
      for (std::size_t i = 0; i < n; ++i) y[i] += p * d[j + 1][i];
    }
    return y;
  }
  void store(double time, const V& y) {
    if (!result.times.empty() && result.times.back() == time) {
      result.states.back() = y;
      return;
    }
    result.times.push_back(time);
    result.states.push_back(y);
  }

 public:
  Bdf(const OdeRhs& rhs_, double t0, const V& y0, double end, const BdfOptions& options,
      const Lapack* lapack)
      : rhs(rhs_),
        o(options),
        n(y0.size()),
        m(o.physical_dimension ? o.physical_dimension : n),
        atol(n),
        factors(m, std::sqrt(eps)),
        jac(n * m),
        rtol(std::max(o.relative_tolerance, 100 * eps)),
        newton_tol(std::max(10 * eps / rtol, std::min(.03, std::sqrt(rtol)))),
        t(t0),
        bound(end),
        lu(n, m, lapack) {
    result.time = t;
    result.state = y0;
    result.statistics.linear_solver = lapack
                                          ? "LAPACK dense LU (passive ledgers eliminated)"
                                          : "native pivoted dense LU (passive ledgers eliminated)";
    for (std::size_t i = 0; i < n; ++i)
      atol[i] = o.absolute_tolerance[o.absolute_tolerance.size() == 1 ? 0 : i];
    for (auto& row : d) row.assign(n, 0);
    double kappa[] = {0, -.185, -1.0 / 9, -.0823, -.0415, 0};
    if (!o.ndf) for (double& value : kappa) value = 0.;
    for (int j = 1; j <= 5; ++j) gamma[j] = gamma[j - 1] + 1.0 / j;
    for (int j = 0; j <= 5; ++j) {
      alpha[j] = (1 - kappa[j]) * gamma[j];
      ec[j] = kappa[j] * gamma[j] + 1.0 / (j + 1);
    }
    coloring();
  }
  BdfResult run(const OdeAcceptedCallback& accepted, const std::vector<OdeEvent>& events,
                const OdeSampleCallback& sample) {
    store(t, result.state);
    std::size_t si = 0;
    auto emit = [&](double time, const V& y) {
      store(time, y);
      if (sample) sample(time, y);
    };
    while (si < o.sample_times.size() && o.sample_times[si] == t) {
      emit(t, result.state);
      ++si;
    }
    if (t == bound) {
      result.success = true;
      result.message = "Reached integration endpoint";
      return result;
    }
    try {
      V f;
      evaluate(t, result.state, f);
      h = o.first_step > 0 ? std::min(o.first_step, bound - t) : initial_step(result.state, f);
      d[0] = result.state;
      for (std::size_t i = 0; i < n; ++i) d[1][i] = h * f[i];
      jacobian(t, result.state, &f);
      std::vector<double> previous(events.size());
      for (std::size_t i = 0; i < events.size(); ++i) {
        previous[i] = events[i].value(t, result.state);
        if (!std::isfinite(previous[i])) throw std::runtime_error("Non-finite initial event value");
      }
      while (t < bound) {
        if (result.statistics.accepted_steps >= o.maximum_steps)
          throw std::runtime_error("BDF maximum accepted-step count exceeded");
        double old = t;
        step();
        double stop = t;
        std::vector<OdeEventHit> hits;
        for (std::size_t i = 0; i < events.size(); ++i) {
          double g0 = previous[i], g1 = events[i].value(t, result.state);
          if (!std::isfinite(g1)) throw std::runtime_error("Non-finite event value");
          bool up = g0 <= 0 && g1 >= 0, down = g0 >= 0 && g1 <= 0;
          bool crossing = (events[i].direction == 0 && (up || down)) ||
                          (events[i].direction > 0 && up) || (events[i].direction < 0 && down);
          previous[i] = g1;
          if (!crossing) continue;
          double left = old, right = t, gl = g0;
          if (g0 == 0)
            right = left;
          else if (g1 == 0)
            left = right;
          else
            for (int k = 0; k < 100; ++k) {
              if (right - left <= 4 * eps * std::max({1.0, std::abs(left), std::abs(right)})) break;
              double mid = left + .5 * (right - left), gm = events[i].value(mid, interpolate(mid));
              if (!std::isfinite(gm)) throw std::runtime_error("Non-finite event root value");
              if (gm == 0) {
                left = right = mid;
                break;
              }
              if ((gl <= 0 && gm >= 0) || (gl >= 0 && gm <= 0))
                right = mid;
              else {
                left = mid;
                gl = gm;
              }
            }
          double root = .5 * (left + right);
          hits.push_back({i, events[i].name, root, interpolate(root)});
          if (events[i].terminal) stop = std::min(stop, root);
        }
        std::sort(hits.begin(), hits.end(),
                  [](const auto& a, const auto& b) { return a.time < b.time; });
        for (const auto& hit : hits)
          if (hit.time <= stop) {
            result.events.push_back(hit);
            if (events[hit.index].terminal) result.stopped_by_event = true;
          }
        while (si < o.sample_times.size() && o.sample_times[si] <= stop) {
          double time = o.sample_times[si++];
          emit(time, interpolate(time));
        }
        if (result.stopped_by_event) {
          result.time = stop;
          result.state = interpolate(stop);
        }
        if (o.store_accepted_steps) store(result.time, result.state);
        if (accepted && !accepted(result.time, result.state, result.statistics))
          result.stopped_by_callback = true;
        if (result.stopped_by_event || result.stopped_by_callback) break;
      }
      result.success = true;
      result.message = result.stopped_by_event      ? "Stopped at terminal event"
                       : result.stopped_by_callback ? "Stopped by accepted-step callback"
                                                    : "Reached integration endpoint";
    } catch (const std::exception& e) {
      result.success = false;
      result.message = e.what();
    }
    store(result.time, result.state);
    return result;
  }
};
}  // namespace
BdfResult integrate_bdf(const OdeRhs& rhs, double t0, const OdeState& y0, double bound,
                        const BdfOptions& o, const OdeAcceptedCallback& accepted,
                        const std::vector<OdeEvent>& events, const OdeSampleCallback& sample) {
  if (!rhs) throw std::invalid_argument("ODE RHS is empty");
  if(o.maximum_order<1||o.maximum_order>5)throw std::invalid_argument("Native BDF maximum order must be 1 to 5");
  if (!std::isfinite(t0) || !std::isfinite(bound) || bound < t0)
    throw std::invalid_argument("BDF requires finite forward time interval");
  if (y0.empty() || !finite(y0) || y0.size() > static_cast<std::size_t>(INT_MAX))
    throw std::invalid_argument("Invalid BDF initial state");
  if (!(o.relative_tolerance > 0) || !std::isfinite(o.relative_tolerance))
    throw std::invalid_argument("Relative tolerance must be finite and positive");
  if (o.absolute_tolerance.size() != 1 && o.absolute_tolerance.size() != y0.size())
    throw std::invalid_argument("Absolute tolerance must be scalar or state-sized");
  for (double a : o.absolute_tolerance)
    if (!(a > 0) || !std::isfinite(a))
      throw std::invalid_argument("Absolute tolerances must be finite and positive");
  if (!(o.maximum_step > 0) || std::isnan(o.maximum_step) || o.first_step < 0 ||
      !std::isfinite(o.first_step) || !o.maximum_steps || o.physical_dimension > y0.size())
    throw std::invalid_argument("Invalid BDF step/dimension limit");
  double last = t0;
  for (double time : o.sample_times) {
    if (!std::isfinite(time) || time < last || time > bound)
      throw std::invalid_argument("Sample times must be ordered within interval");
    last = time;
  }
  for (const auto& event : events)
    if (!event.value) throw std::invalid_argument("Empty event function");
  static const Lapack lapack;
  Bdf b(rhs, t0, y0, bound, o,
        o.use_lapack_if_available && lapack.getrf && lapack.getrs ? &lapack : nullptr);
  return b.run(accepted, events, sample);
}
}  // namespace cycle

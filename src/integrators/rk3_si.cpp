#include "rk3_si.h"

#include <algorithm>
#include <mdspan>
#include <stdexcept>
#include <vector>

#include "utils.h"

extern "C"
{
  void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
  void dgbsv_(int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info);
}

namespace
{
/**
 * \brief Solves the implicit concentration updates of every component; node 0 is the fixed inlet.
 *
 * L is the upwind convection and axial dispersion operator with a zero-gradient outlet; s = -sink is the sorption
 * source. A first-order stage solves (I - dt L) c = c_old + dt s (tridiagonal), the final stage
 * (I + dt^2 L^2) c = c_old - dt^2 L s (pentadiagonal, LAPACK band storage, column-major).
 */
void solveConcentrations(const Column& column, double dt, bool finalStage, std::vector<double>& solved)
{
  const size_t nodes = column.numberOfGridPoints + 1;
  const size_t last = column.numberOfGridPoints;
  const size_t numberOfComponents = column.numberOfComponents;
  const double idx = 1.0 / column.resolution;
  const double idx2 = idx * idx;
  const std::vector<Component>& components = column.physisorptionMixtures.front().components;

  // Row i of L: lower[i - 1] c_{i-1} + diag[i] c_i + upper[i] c_{i+1}; row 0 is zero (fixed inlet).
  std::vector<double> lower(last), diag(nodes), upper(last), source(nodes), rhs(nodes);
  std::vector<double> subdiagonal(last), diagonal(nodes), superdiagonal(last), band(7 * nodes);
  std::vector<int> pivots(nodes);
  std::mdspan<double, std::dextents<size_t, 2>, std::layout_left> ab(band.data(), 7, nodes);
  int n = static_cast<int>(nodes);
  int nrhs = 1;
  int info = 0;

  for (size_t comp = 0; comp < numberOfComponents; ++comp)
  {
    const double dispersion = components[comp].axialDispersionCoefficient;
    diag[0] = 0.0;
    upper[0] = 0.0;
    source[0] = 0.0;
    for (size_t i = 1; i < nodes; ++i)
    {
      const double velocity = column.interstitialGasVelocity[i];
      lower[i - 1] = velocity * idx + dispersion * idx2;
      diag[i] = -velocity * idx - (i < last ? 2.0 : 1.0) * dispersion * idx2;
      if (i < last) upper[i] = dispersion * idx2;
      source[i] = -column.bulkSpeciesSink[i * numberOfComponents + comp];
    }

    rhs[0] = column.concentration[comp];
    if (!finalStage)
    {
      for (size_t i = 1; i < nodes; ++i) rhs[i] = column.concentration[i * numberOfComponents + comp] + dt * source[i];
      diagonal[0] = 1.0;
      superdiagonal[0] = 0.0;
      for (size_t i = 1; i < nodes; ++i)
      {
        subdiagonal[i - 1] = -dt * lower[i - 1];
        diagonal[i] = 1.0 - dt * diag[i];
        if (i < last) superdiagonal[i] = -dt * upper[i];
      }
      dgtsv_(&n, &nrhs, subdiagonal.data(), diagonal.data(), superdiagonal.data(), rhs.data(), &n, &info);
      if (info != 0) throw std::runtime_error("dgtsv failed in the SIRK3 concentration update");
    }
    else
    {
      // A = I + dt^2 L^2 in band storage: ab[kl + ku + i - j, j] = A[i, j] with kl = ku = 2.
      const double dt2 = dt * dt;
      std::fill(band.begin(), band.end(), 0.0);
      for (size_t i = 0; i + 2 < nodes; ++i)
      {
        ab[6, i] = dt2 * lower[i + 1] * lower[i];
        ab[2, i + 2] = dt2 * upper[i] * upper[i + 1];
      }
      for (size_t i = 0; i < last; ++i)
      {
        ab[5, i] = dt2 * lower[i] * (diag[i] + diag[i + 1]);
        ab[3, i + 1] = dt2 * upper[i] * (diag[i] + diag[i + 1]);
      }
      ab[4, 0] = 1.0;
      for (size_t i = 1; i < nodes; ++i)
      {
        const double couplingAbove = i < last ? lower[i] * upper[i] : 0.0;
        const double sourceAbove = i < last ? upper[i] * source[i + 1] : 0.0;
        ab[4, i] = 1.0 + dt2 * (lower[i - 1] * upper[i - 1] + diag[i] * diag[i] + couplingAbove);
        rhs[i] = column.concentration[i * numberOfComponents + comp] -
                 dt2 * (lower[i - 1] * source[i - 1] + diag[i] * source[i] + sourceAbove);
      }
      int kl = 2;
      int ku = 2;
      int ldab = 7;
      dgbsv_(&n, &kl, &ku, &nrhs, band.data(), &ldab, pivots.data(), rhs.data(), &n, &info);
      if (info != 0) throw std::runtime_error("dgbsv failed in the SIRK3 concentration update");
    }

    for (size_t grid = 0; grid < nodes; ++grid) solved[grid * numberOfComponents + comp] = rhs[grid];
  }
}
}  // namespace

void SemiImplicitRungeKutta3::step(Column& column, double startTime, double dt, Timing& timings)
{
  if (column.convectionScheme != Column::ConvectionScheme::Upwind || column.numberOfAdsorbents != 1)
  {
    throw std::runtime_error("SemiImplicitRungeKutta3 supports a single adsorbent with ConvectionScheme Upwind only");
  }
  auto totalTimer = timings.scoped(timings.total);

  const size_t numberOfComponents = column.numberOfComponents;
  const std::vector<Component>& components = column.physisorptionMixtures.front().components;
  const std::vector<double> initialPhysisorption(column.physisorption.begin(), column.physisorption.end());
  const std::vector<double> initialConcentration(column.concentration.begin(), column.concentration.end());
  std::vector<double> solved(column.concentration.size());

  // One SSP-RK3 stage: the LDF loading update (implicit in its own rate), then the concentration solve at the
  // updated loadings, each blended with the initial state as alpha * initial + beta * stage.
  auto stage = [&](double alpha, double beta, double elapsedTime, bool finalStage)
  {
    for (size_t i = 0; i < column.physisorption.size(); ++i)
    {
      const double rate = components[i % numberOfComponents].massTransferCoefficient;
      if (finalStage)
      {
        const double rate2 = dt * dt * rate * rate;
        column.physisorption[i] = (column.physisorption[i] + rate2 * column.equilibriumPhysisorption[i]) / (1.0 + rate2);
      }
      else
      {
        column.physisorption[i] =
            alpha * initialPhysisorption[i] +
            beta * (column.physisorption[i] + dt * rate * column.equilibriumPhysisorption[i]) / (1.0 + dt * rate);
      }
    }
    timings.measure(timings.computeDerivatives,
                    [&]
                    {
                      computeDerivatives(column, elapsedTime);
                      solveConcentrations(column, dt, finalStage, solved);
                    });
    for (size_t i = 0; i < solved.size(); ++i)
    {
      column.concentration[i] = finalStage ? solved[i] : alpha * initialConcentration[i] + beta * solved[i];
    }
    if (finalStage) clampNonnegative(column.state);
    precompute(column, timings);
  };

  stage(0.0, 1.0, startTime + dt, false);
  stage(0.75, 0.25, startTime + 0.5 * dt, false);
  stage(1.0 / 3.0, 2.0 / 3.0, startTime + dt, false);
  stage(0.0, 1.0, startTime + dt, true);
}

#include "block_banded.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace
{
/**
 * \brief Inverts the column-major n x n matrix in place by Gauss-Jordan elimination with partial pivoting; false if it
 *        is singular. work holds n^2 values.
 */
bool invert(std::span<double> matrix, size_t n, std::span<double> work)
{
  // work starts as the identity and becomes the inverse while matrix is reduced to the identity.
  std::fill(work.begin(), work.end(), 0.0);
  for (size_t i = 0; i < n; ++i) work[i * n + i] = 1.0;
  for (size_t col = 0; col < n; ++col)
  {
    size_t pivot = col;
    for (size_t row = col + 1; row < n; ++row)
    {
      if (std::abs(matrix[col * n + row]) > std::abs(matrix[col * n + pivot])) pivot = row;
    }
    const double value = matrix[col * n + pivot];
    if (value == 0.0 || !std::isfinite(value)) return false;
    if (pivot != col)
    {
      for (size_t k = 0; k < n; ++k)
      {
        std::swap(matrix[k * n + pivot], matrix[k * n + col]);
        std::swap(work[k * n + pivot], work[k * n + col]);
      }
    }
    const double inverse = 1.0 / value;
    for (size_t k = 0; k < n; ++k)
    {
      matrix[k * n + col] *= inverse;
      work[k * n + col] *= inverse;
    }
    for (size_t row = 0; row < n; ++row)
    {
      const double factor = matrix[col * n + row];
      if (row == col || factor == 0.0) continue;
      for (size_t k = 0; k < n; ++k)
      {
        matrix[k * n + row] -= factor * matrix[k * n + col];
        work[k * n + row] -= factor * work[k * n + col];
      }
    }
  }
  std::copy(work.begin(), work.end(), matrix.begin());
  return true;
}

/**
 * \brief result -= left right, for column-major n x n blocks.
 */
void subtractProduct(std::span<double> result, std::span<const double> left, std::span<const double> right, size_t n)
{
  for (size_t col = 0; col < n; ++col)
  {
    for (size_t k = 0; k < n; ++k)
    {
      const double factor = right[col * n + k];
      if (factor == 0.0) continue;
      for (size_t row = 0; row < n; ++row) result[col * n + row] -= left[k * n + row] * factor;
    }
  }
}
}  // namespace

BlockBandedMatrix::BlockBandedMatrix(size_t numberOfBlocks, size_t blockSize, size_t lowerBlocks, size_t upperBlocks)
    : numberOfBlocks(numberOfBlocks),
      blockSize(blockSize),
      lowerBlocks(lowerBlocks),
      upperBlocks(upperBlocks),
      blocks(numberOfBlocks * (lowerBlocks + upperBlocks + 1) * blockSize * blockSize, 0.0),
      work(2 * blockSize * blockSize, 0.0)
{
}

void BlockBandedMatrix::setZero() noexcept { std::fill(blocks.begin(), blocks.end(), 0.0); }

void BlockBandedMatrix::scaleAndShift(double factor, double shift) noexcept
{
  for (double& value : blocks) value *= factor;
  for (size_t block = 0; block < numberOfBlocks; ++block)
  {
    for (size_t i = 0; i < blockSize; ++i) (*this)(block, block, i, i) += shift;
  }
}

bool BlockBandedMatrix::factorize(size_t nonzeroUpperBlocks)
{
  const size_t upper = std::min(upperBlocks, nonzeroUpperBlocks);
  factorizedUpperBlocks = upper;
  const size_t area = blockSize * blockSize;
  const auto block = [&](size_t blockRow, size_t blockColumn)
  { return std::span<double>(blocks.data() + blockOffset(blockRow, blockColumn), area); };
  const std::span<double> product(work.data(), area);
  const std::span<double> scratch(work.data() + area, area);
  for (size_t k = 0; k < numberOfBlocks; ++k)
  {
    // U_kj = A_kj - sum_{m < k} L_km U_mj for the blocks of row k on and above the diagonal, over the m for which both
    // factors lie in their bands; with a single block above the diagonal, U_{k,k+1} is that of A.
    for (size_t j = k; j < std::min(numberOfBlocks, k + upper + 1); ++j)
    {
      for (size_t m = std::max(k > lowerBlocks ? k - lowerBlocks : 0, j > upper ? j - upper : 0); m < k; ++m)
      {
        subtractProduct(block(k, j), block(k, m), block(m, j), blockSize);
      }
    }
    if (!invert(block(k, k), blockSize, scratch)) return false;
    for (size_t i = k + 1; i < std::min(numberOfBlocks, k + lowerBlocks + 1); ++i)
    {
      // L_ik = (A_ik - sum_{m < k} L_im U_mk) U_kk^-1
      for (size_t m = std::max(i > lowerBlocks ? i - lowerBlocks : 0, k > upper ? k - upper : 0); m < k; ++m)
      {
        subtractProduct(block(i, k), block(i, m), block(m, k), blockSize);
      }
      const double* lower = blocks.data() + blockOffset(i, k);
      const double* inverse = blocks.data() + blockOffset(k, k);
      std::fill(product.begin(), product.end(), 0.0);
      for (size_t col = 0; col < blockSize; ++col)
      {
        for (size_t j = 0; j < blockSize; ++j)
        {
          const double factor = inverse[col * blockSize + j];
          for (size_t row = 0; row < blockSize; ++row) product[col * blockSize + row] += lower[j * blockSize + row] * factor;
        }
      }
      std::copy(product.begin(), product.end(), block(i, k).begin());
    }
  }
  return true;
}

void BlockBandedMatrix::solve(std::span<double> x) const
{
  // L y = x: y_k = x_k - sum_j L_kj y_j.
  for (size_t k = 1; k < numberOfBlocks; ++k)
  {
    double* target = x.data() + k * blockSize;
    for (size_t j = k > lowerBlocks ? k - lowerBlocks : 0; j < k; ++j)
    {
      const double* lower = blocks.data() + blockOffset(k, j);
      const double* source = x.data() + j * blockSize;
      for (size_t col = 0; col < blockSize; ++col)
      {
        for (size_t row = 0; row < blockSize; ++row) target[row] -= lower[col * blockSize + row] * source[col];
      }
    }
  }
  // U x = y: x_k = U_kk^-1 (y_k - sum_j U_kj x_j).
  for (size_t k = numberOfBlocks; k-- > 0;)
  {
    double* target = x.data() + k * blockSize;
    for (size_t j = k + 1; j < std::min(numberOfBlocks, k + factorizedUpperBlocks + 1); ++j)
    {
      const double* upper = blocks.data() + blockOffset(k, j);
      const double* next = x.data() + j * blockSize;
      for (size_t col = 0; col < blockSize; ++col)
      {
        for (size_t row = 0; row < blockSize; ++row) target[row] -= upper[col * blockSize + row] * next[col];
      }
    }
    const double* inverse = blocks.data() + blockOffset(k, k);
    for (size_t row = 0; row < blockSize; ++row)
    {
      double sum = 0.0;
      for (size_t col = 0; col < blockSize; ++col) sum += inverse[col * blockSize + row] * target[col];
      work[row] = sum;
    }
    std::copy(work.begin(), work.begin() + static_cast<std::ptrdiff_t>(blockSize), target);
  }
}

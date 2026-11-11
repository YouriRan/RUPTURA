#pragma once

#include <cstddef>
#include <limits>
#include <span>
#include <vector>

/**
 * \brief Square matrix of dense blocks that is nonzero from lowerBlocks blocks below the diagonal to upperBlocks blocks
 *        above it, with an in-place block LU factorization.
 *
 * factorize() replaces the matrix by A = L U, with L block-lower with identity diagonal blocks and U block-upper. It
 * keeps the inverses of the diagonal blocks of U and does not exchange rows between blocks, which suits the
 * diagonally dominant stage matrices I - gamma h J of ROS34PW2; L and U then keep the bands of A. Vectors are
 * block-major: entry i of block k is x[k blockSize + i]. With blocks of size b, factorize() costs about
 * (2 lowerBlocks upperBlocks + 2) b^3 operations per block row and solve() about (lowerBlocks + upperBlocks + 1) 2 b^2.
 */
struct BlockBandedMatrix
{
  BlockBandedMatrix() = default;
  BlockBandedMatrix(size_t numberOfBlocks, size_t blockSize, size_t lowerBlocks, size_t upperBlocks = 1);

  /**
   * \brief Entry (row, col) of block (blockRow, blockColumn); blockColumn runs from blockRow - lowerBlocks to
   *        blockRow + upperBlocks.
   */
  double& operator()(size_t blockRow, size_t blockColumn, size_t row, size_t col) noexcept
  {
    return blocks[blockOffset(blockRow, blockColumn) + col * blockSize + row];
  }
  double operator()(size_t blockRow, size_t blockColumn, size_t row, size_t col) const noexcept
  {
    return blocks[blockOffset(blockRow, blockColumn) + col * blockSize + row];
  }

  void setZero() noexcept;

  /**
   * \brief Multiplies every entry by factor and adds shift to the diagonal: forms I - gamma h J from J.
   */
  void scaleAndShift(double factor, double shift) noexcept;

  /**
   * \brief Factorizes the matrix in place; false if a diagonal block of U is singular. The blocks more than
   *        nonzeroUpperBlocks above the diagonal must be zero; factorize() and solve() skip them.
   */
  bool factorize(size_t nonzeroUpperBlocks = std::numeric_limits<size_t>::max());

  /**
   * \brief Overwrites x with A^-1 x, after factorize().
   */
  void solve(std::span<double> x) const;

  size_t numberOfBlocks{0};  ///< Block rows and columns.
  size_t blockSize{0};       ///< Rows and columns of a block.
  size_t lowerBlocks{0};     ///< Blocks below the diagonal.
  size_t upperBlocks{1};     ///< Blocks above the diagonal.

 private:
  /// Start of block (blockRow, blockColumn): each block row stores lowerBlocks + upperBlocks + 1 column-major blocks,
  /// from block column blockRow - lowerBlocks to blockRow + upperBlocks.
  size_t blockOffset(size_t blockRow, size_t blockColumn) const noexcept
  {
    return (blockRow * (lowerBlocks + upperBlocks + 1) + blockColumn + lowerBlocks - blockRow) * blockSize * blockSize;
  }

  std::vector<double> blocks;        ///< The blocks; after factorize() L below, inverted U on, and U above the diagonal.
  size_t factorizedUpperBlocks{1};   ///< Blocks above the diagonal that the last factorize() processed.
  mutable std::vector<double> work;  ///< Two blocks of scratch.
};

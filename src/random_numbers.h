#pragma once

#include <random>

/**
 * \brief Singleton random-number generator used by optimization routines.
 */
class RandomNumber
{
 public:
  /**
   * \brief Returns a uniform random number in the range [0, 1).
   */
  static double Uniform() { return getInstance().uniformDistribution(getInstance().randomEngine); }

  /**
   * \brief Returns a standard-normal random number.
   */
  static double Gaussian() { return getInstance().normalDistribution(getInstance().randomEngine); }

  /**
   * \brief Returns a uniform integer in the inclusive range [i, j].
   */
  static size_t Integer(size_t i, size_t j)
  {
    return i + static_cast<size_t>(static_cast<double>(j + 1 - i) * Uniform());
  }

  /**
   * \brief Returns a uniform random 64-bit unsigned integer.
   */
  static uint64_t UInt64() { return getInstance().uniformUInt64Distribution(getInstance().randomEngine); }

 private:
  /**
   * \brief Seeds the random-number engine from std::random_device.
   */
  RandomNumber()
  {
    std::random_device rd;
    randomEngine = std::mt19937_64(rd());
    // randomEngine = std::mt19937_64(10);
    uniformDistribution = std::uniform_real_distribution<double>(0.0, 1.0);
    uniformUInt64Distribution = std::uniform_int_distribution<uint64_t>();
    normalDistribution = std::normal_distribution<double>();
  }

  /**
   * \brief Destroys the singleton random-number generator.
   */
  ~RandomNumber() {}

  /**
   * \brief Returns the process-wide random-number generator instance.
   */
  static RandomNumber& getInstance()
  {
    static RandomNumber s;
    return s;
  }

  /**
   * \brief Copy construction is disabled for the singleton.
   */
  RandomNumber(RandomNumber const&) = delete;

  /**
   * \brief Copy assignment is disabled for the singleton.
   */
  RandomNumber& operator=(RandomNumber const&) = delete;

  std::mt19937_64 randomEngine;                                       ///< Pseudo-random number engine.
  std::uniform_real_distribution<double> uniformDistribution;         ///< Uniform real distribution on [0, 1).
  std::normal_distribution<double> normalDistribution;                ///< Standard-normal distribution.
  std::uniform_int_distribution<uint64_t> uniformUInt64Distribution;  ///< Uniform 64-bit integer distribution.
};

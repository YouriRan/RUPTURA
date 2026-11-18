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
   * \brief Returns a uniform integer in the inclusive range [i, j].
   */
  static size_t Integer(size_t i, size_t j)
  {
    return i + static_cast<size_t>(static_cast<double>(j + 1 - i) * Uniform());
  }

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
};

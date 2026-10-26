#pragma once

#include <chrono>
#include <filesystem>
#include <random>
#include <string>
#include <system_error>

/**
 * \brief Directory for the temporary input files of this test process.
 *
 * The tests write their inputs under fixed file names. One directory per process keeps two ruptura_tests runs
 * (two build trees or two sessions) from overwriting each other's files. It is removed when the process exits.
 */
inline const std::filesystem::path& testTemporaryDirectory()
{
  struct Directory
  {
    std::filesystem::path path;

    Directory()
    {
      const auto timestamp = std::chrono::steady_clock::now().time_since_epoch().count();
      path = std::filesystem::temp_directory_path() /
             ("ruptura-tests-" + std::to_string(timestamp) + "-" + std::to_string(std::random_device{}()));
      std::filesystem::create_directories(path);
    }

    ~Directory()
    {
      std::error_code error;
      std::filesystem::remove_all(path, error);
    }
  };

  static const Directory directory;
  return directory.path;
}

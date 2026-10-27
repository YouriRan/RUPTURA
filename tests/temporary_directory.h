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

/**
 * \brief Changes the working directory for the lifetime of the guard, e.g. for a simulation that writes its output
 * files to the current directory.
 */
class WorkingDirectoryGuard
{
 public:
  explicit WorkingDirectoryGuard(const std::filesystem::path& path) : originalPath(std::filesystem::current_path())
  {
    std::filesystem::create_directories(path);
    std::filesystem::current_path(path);
  }

  ~WorkingDirectoryGuard()
  {
    std::error_code error;
    std::filesystem::current_path(originalPath, error);
  }

 private:
  std::filesystem::path originalPath;
};

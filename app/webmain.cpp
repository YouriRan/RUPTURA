#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

#ifndef RUPTURA_SOURCE_DIR
#define RUPTURA_SOURCE_DIR ""
#endif
#ifndef RUPTURA_LAB_PYTHON_EXECUTABLE
#define RUPTURA_LAB_PYTHON_EXECUTABLE "python3"
#endif

int main(int argc, char* argv[])
{
  namespace fs = std::filesystem;
  fs::path executable = argv[0];
  if (!executable.has_parent_path())
  {
    const char* paths = std::getenv("PATH");
    std::string path = paths ? paths : "";
    size_t start = 0;
    while (start <= path.size())
    {
      const auto end = path.find(':', start);
      const auto candidate = fs::path(path.substr(start, end - start)) / executable;
      if (fs::is_regular_file(candidate)) { executable = candidate; break; }
      if (end == std::string::npos) break;
      start = end + 1;
    }
  }
  executable = fs::weakly_canonical(fs::absolute(executable));
  const auto installed = executable.parent_path().parent_path() / "share/ruptura/python";
  const fs::path packageRoot = fs::is_directory(installed / "ruptura") ? installed : fs::path(RUPTURA_SOURCE_DIR);
  if (fs::is_directory(packageRoot / "ruptura"))
  {
    const char* existing = std::getenv("PYTHONPATH");
    const std::string paths = packageRoot.string() + (existing ? ":" + std::string(existing) : "");
    setenv("PYTHONPATH", paths.c_str(), 1);
  }
  if (!std::getenv("RUPTURA_CLI_EXECUTABLE"))
  {
    const auto sibling = executable.parent_path() / "ruptura";
    if (fs::is_regular_file(sibling)) setenv("RUPTURA_CLI_EXECUTABLE", sibling.c_str(), 1);
#ifdef RUPTURA_CLI_EXECUTABLE
    else if (fs::is_regular_file(RUPTURA_CLI_EXECUTABLE))
      setenv("RUPTURA_CLI_EXECUTABLE", RUPTURA_CLI_EXECUTABLE, 1);
#endif
  }
  const char* configured = std::getenv("RUPTURA_PYTHON_EXECUTABLE");
  std::string python = configured && *configured ? configured : RUPTURA_LAB_PYTHON_EXECUTABLE;
  std::string launch = "import sys\nsys.argv[0] = 'ruptura-lab'\ntry:\n from ruptura.webapp.server import main\n"
                       "except ImportError as error:\n print('Could not load Ruptura Lab: ' + str(error) + "
                       "'\\nSet RUPTURA_PYTHON_EXECUTABLE to an environment containing ruptura[web].', file=sys.stderr)\n"
                       " raise SystemExit(1)\nmain()\n";
  std::vector<char*> arguments{python.data(), const_cast<char*>("-c"), launch.data()};
  for (int i = 1; i < argc; ++i) arguments.push_back(argv[i]);
  arguments.push_back(nullptr);
  execvp(python.c_str(), arguments.data());
  std::cerr << "Could not start Python. Set RUPTURA_PYTHON_EXECUTABLE to an interpreter with ruptura[web] installed.\n";
  return 1;
}

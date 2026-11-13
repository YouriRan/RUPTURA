#!/usr/bin/env bash
# Portable CPU build of the CLI, shared FDM/FVM column, cycle GUI and bindings.
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"
python_build=OFF
gui_build=ON
tests_build=OFF
cvode_build=ON
portable_toolchain=auto
sundials_prefix="${SUNDIALS_ROOT:-}"
build_dir="${RUPTURA_BUILD_DIR:-$project_dir/build}"
jobs="${RUPTURA_BUILD_JOBS:-4}"
while (($#)); do
  case "$1" in
    --python) python_build=ON ;;
    --gui) gui_build=ON ;;
    --no-gui) gui_build=OFF ;;
    --tests) tests_build=ON ;;
    --toolchain) portable_toolchain=yes ;;
    --no-cvode) cvode_build=OFF ;;
    --sundials-prefix) shift; sundials_prefix="${1:?Supply the SUNDIALS install prefix}" ;;
    --help|-h)
      cat <<'HELP'
./install.sh [--python] [--gui | --no-gui] [--tests] [--toolchain] [--no-cvode]
             [--sundials-prefix /path/to/sundials]

Builds build/ruptura and build/app/ruptura_lab with CPU-specific flags disabled.
The browser GUI is included by default; --no-gui builds the CLI without it.
Python packages are installed in a project virtual environment. A compiler with
C++23 <print> and <mdspan>, LAPACK/BLAS, zlib headers and a C compiler are needed.
If the system C++ library is too old, a pinned Zig/Clang toolchain is installed
in that environment. --toolchain selects this tested fallback explicitly.
SUNDIALS 7.5 is downloaded by CMake unless an installed prefix is supplied.
--no-cvode explicitly selects the built-in BDF solver; CVODE inputs then fail
with an explanation. --tests also builds the inherited main and process tests.
Set CXX, CC, RUPTURA_BUILD_JOBS or RUPTURA_BUILD_DIR to override build choices.
No root privileges or system package changes are used.
HELP
      exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
  shift
done
command -v python3 >/dev/null || { echo 'Python 3.11 or newer is required.' >&2; exit 1; }
python3 -m venv "$project_dir/.venv-ruptura"
py="$project_dir/.venv-ruptura/bin/python"
"$py" -m pip install --upgrade 'cmake>=3.28' 'nanobind>=2.2' numpy matplotlib
if [[ "$python_build" == ON || "$gui_build" == ON ]]; then
  "$py" -m pip install pandas 'plotly>=5.24'
fi
if [[ "$tests_build" == ON ]]; then
  "$py" -m pip install pytest
fi
cmake_bin="$project_dir/.venv-ruptura/bin/cmake"
probe_dir="$(mktemp -d)"
trap 'rm -rf -- "$probe_dir"' EXIT
cat > "$probe_dir/probe.cpp" <<'CPP'
#include <print>
#include <mdspan>
int main(){std::println("C++23 standard library available");}
CPP
compiler="${CXX:-c++}"
if [[ "$portable_toolchain" == yes ]] || ! "$compiler" -std=c++23 "$probe_dir/probe.cpp" -o "$probe_dir/probe" > "$probe_dir/compiler.log" 2>&1; then
  echo 'Using the pinned portable Clang/libc++ toolchain.'
  "$py" -m pip install 'ziglang==0.17.1'
  zig_bin="$("$py" -c 'import pathlib,ziglang;print(pathlib.Path(ziglang.__file__).parent/"zig")')"
  mkdir -p "$project_dir/.build-tools"
  printf '#!/usr/bin/env bash\nexec %q c++ -mcpu=baseline "$@"\n' "$zig_bin" > "$project_dir/.build-tools/cxx"
  chmod +x "$project_dir/.build-tools/cxx"
  compiler="$project_dir/.build-tools/cxx"
fi
cmake_args=(-S "$project_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release
  "-DCMAKE_CXX_COMPILER=$compiler" -DRUPTURA_NATIVE_ARCH=OFF -DBUILD_DOCS=OFF
  -DCMAKE_CXX_SCAN_FOR_MODULES=OFF "-DBUILD_APP=$gui_build"
  -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON
  "-DRUPTURA_LAB_PYTHON_EXECUTABLE=$py"
  "-DBUILD_SUNDIALS=$cvode_build" "-DBUILD_PYTHON=$python_build"
  "-DRUPTURA_PYTHON_EXECUTABLE=$py" "-DBUILD_TESTING=$tests_build")
if [[ -n "$sundials_prefix" ]]; then
  cmake_args+=("-DSUNDIALS_ROOT=$sundials_prefix" -DRUPTURA_USE_SYSTEM_SUNDIALS=ON)
fi
"$cmake_bin" "${cmake_args[@]}"
"$cmake_bin" --build "$build_dir" --parallel "$jobs"
if [[ "$tests_build" == ON ]]; then
  "$project_dir/.venv-ruptura/bin/ctest" --test-dir "$build_dir" --output-on-failure --parallel "$jobs"
fi
printf '\nBuilt: %s/ruptura\n' "$build_dir"
printf 'Check an input: %q/ruptura /path/to/simulation.json --validate\n' "$build_dir"
if [[ "$gui_build" == ON ]]; then
  printf 'Start the cycle GUI: %q/app/ruptura_lab\n' "$build_dir"
fi
if [[ "$python_build" == ON ]]; then
  printf 'Python/Lab: source %q/.venv-ruptura/bin/activate\n' "$project_dir"
  printf 'Then: PYTHONPATH=%q python -m ruptura.webapp --help\n' "$project_dir"
fi

cmake_minimum_required(VERSION 4.3)
project(ruptura LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(LAPACK)

option(BUILD_SUNDIALS "Build with SUNDIALS-based CVODE integrator support" ON)
option(BUILD_PYTHON "Build nanobind-based Python module" OFF)

set(SUNDIALS_ROOT "" CACHE PATH
  "Prefix of a preinstalled SUNDIALS (expects include/ and lib/ under this directory)"
)

set(NANOBIND_SOURCE_DIR "" CACHE PATH
  "Path to a nanobind source checkout. Leave empty to use find_package()/pip-installed nanobind."
)

set(RUPTURA_PYTHON_MODULE_NAME "_ruptura" CACHE STRING
  "Python import name for the nanobind module"
)

if(SUNDIALS_ROOT)
  set(_ruptura_default_use_system_sundials ON)
else()
  set(_ruptura_default_use_system_sundials OFF)
endif()

option(RUPTURA_USE_SYSTEM_SUNDIALS
  "Use a preinstalled SUNDIALS instead of FetchContent"
  ${_ruptura_default_use_system_sundials}
)

if(BUILD_SUNDIALS)

  # ---- SuiteSparse / KLU (system dependency) ----
  find_path(KLU_INCLUDE_DIR
    NAMES klu.h
    HINTS
      ${KLU_ROOT}
      /opt/homebrew
      /usr/local
      /usr
    PATH_SUFFIXES include/suitesparse suitesparse
  )

  find_library(KLU_LIBRARY
    NAMES klu
    HINTS
      ${KLU_ROOT}
      /opt/homebrew
      /usr/local
      /usr
    PATH_SUFFIXES lib lib64
  )

  if(NOT KLU_INCLUDE_DIR OR NOT KLU_LIBRARY)
    message(FATAL_ERROR
      "Could not find SuiteSparse KLU (klu.h / libklu). "
      "Install SuiteSparse and/or set KLU_ROOT/KLU_INCLUDE_DIR/KLU_LIBRARY."
    )
  endif()

  get_filename_component(KLU_LIBRARY_DIR "${KLU_LIBRARY}" DIRECTORY)

  if(RUPTURA_USE_SYSTEM_SUNDIALS)
    if(NOT SUNDIALS_ROOT)
      message(FATAL_ERROR
        "RUPTURA_USE_SYSTEM_SUNDIALS=ON requires -DSUNDIALS_ROOT=/path/to/sundials/prefix"
      )
    endif()

    find_package(SUNDIALS CONFIG REQUIRED
      HINTS
        "${SUNDIALS_ROOT}"
      PATH_SUFFIXES
        lib/cmake/sundials
        lib64/cmake/sundials
        lib/cmake/SUNDIALS
        lib64/cmake/SUNDIALS
    )

  else()
    include(FetchContent)
    FetchContent_Declare(sundials
      GIT_REPOSITORY https://github.com/LLNL/sundials.git
      GIT_TAG v7.5.0
    )

    set(SUNDIALS_BUILD_CVODE ON CACHE BOOL "" FORCE)
    set(SUNDIALS_BUILD_NVEC_SERIAL ON CACHE BOOL "" FORCE)
    set(SUNDIALS_BUILD_NVEC_MANYVECTOR ON CACHE BOOL "" FORCE)
    set(EXAMPLES_ENABLE OFF CACHE BOOL "" FORCE)
    set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)

    # Enable KLU inside the vendored SUNDIALS build.
    set(ENABLE_KLU ON CACHE BOOL "" FORCE)
    set(KLU_INCLUDE_DIR "${KLU_INCLUDE_DIR}" CACHE PATH "" FORCE)
    set(KLU_LIBRARY_DIR "${KLU_LIBRARY_DIR}" CACHE PATH "" FORCE)

    FetchContent_MakeAvailable(sundials)
  endif()

endif()

if(BUILD_PYTHON)
  set(DEV_MODULE Development.Module)

  set(RUPTURA_PYTHON_EXECUTABLE "" CACHE FILEPATH
    "Exact Python interpreter to use for the Python module"
  )

  # Prefer the active virtualenv/conda environment when one is activated.
  # Use ONLY instead of FIRST if you want configuration to fail rather than
  # silently fall back to a system Python.
  set(Python_FIND_VIRTUALENV FIRST)
  set(Python_FIND_STRATEGY LOCATION)

  if(RUPTURA_PYTHON_EXECUTABLE)
    set(Python_EXECUTABLE "${RUPTURA_PYTHON_EXECUTABLE}" CACHE FILEPATH "" FORCE)
  endif()

  find_package(Python 3.12 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

  message(STATUS "Using Python executable: ${Python_EXECUTABLE}")
  message(STATUS "Using Python version: ${Python_VERSION}")

  if(NANOBIND_SOURCE_DIR)
    add_subdirectory("${NANOBIND_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/nanobind")
  endif()

  if(NOT COMMAND nanobind_add_module)
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
      OUTPUT_VARIABLE nanobind_ROOT
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )

    if(nanobind_ROOT)
      list(APPEND CMAKE_PREFIX_PATH "${nanobind_ROOT}")
      find_package(nanobind CONFIG QUIET)
    endif()
  endif()

  if(NOT COMMAND nanobind_add_module)
    find_package(nanobind CONFIG QUIET)
  endif()

  if(NOT COMMAND nanobind_add_module)
    message(FATAL_ERROR
      "Could not find nanobind. Install it with Python/pip, or configure with "
      "-DNANOBIND_SOURCE_DIR=/path/to/nanobind."
    )
  endif()
endif()

set(RUPTURA_SOURCES
  src/breakthrough.cpp
  src/column.cpp
  src/component.cpp
  src/fitting.cpp
  src/inputreader.cpp
  src/isotherm.cpp
  src/mixture_prediction.cpp
  src/multi_site_isotherm.cpp
  src/random_numbers.cpp
  src/special_functions.cpp
  src/swing_adsorption.cpp
  src/integrators/rk3.cpp
  src/integrators/rk3_si.cpp
  src/integrators/cvode.cpp
  src/integrators/compute.cpp
)

set(RUPTURA_COMPILE_OPTIONS
  -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align
  -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wno-reorder-ctor
  -Wnull-dereference -Wdouble-promotion -Wformat=2
  -O3 -march=native -ffast-math -fomit-frame-pointer -ftree-vectorize -funroll-loops
)

set(RUPTURA_INCLUDE_DIRS
  src
  src/integrators
  $<$<BOOL:${BUILD_SUNDIALS}>:${KLU_INCLUDE_DIR}>
)

set(RUPTURA_LINK_LIBRARIES
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::cvode>
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::nvecserial>
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::sunlinsolklu>
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::sunmatrixsparse>
)

# ---------------------------------------------------------------------------
# Stage 1: build the reusable core library.
#
# Keep this target free of BUILD_PYTHON so it remains a normal C++ library that can
# be linked into both the executable and the Python module. Python-specific
# wrappers should live in ruptura/pybindings.cpp.
# ---------------------------------------------------------------------------
add_library(rupturalib STATIC
  ${RUPTURA_SOURCES}
)

set_target_properties(rupturalib PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)

target_compile_options(rupturalib PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(rupturalib PUBLIC
  BUILD_SUNDIALS=$<BOOL:${BUILD_SUNDIALS}>
)

target_include_directories(rupturalib PUBLIC
  ${RUPTURA_INCLUDE_DIRS}
)

target_link_libraries(rupturalib PUBLIC
  ${RUPTURA_LINK_LIBRARIES}
)


# ---------------------------------------------------------------------------
# Stage 2: app and pybuild
# ---------------------------------------------------------------------------
option(BUILD_APP "Build Ruptura command-line executable" ON)

if(BUILD_APP)
  add_subdirectory(app)
endif()

if(BUILD_PYTHON)
  add_subdirectory(ruptura)
endif()

add_subdirectory(docs)

cmake_minimum_required(VERSION 3.14)
project(ruptura LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(LAPACK)

option(BUILD_SUNDIALS "Build with SUNDIALS-based CVODE integrator support" ON)

set(SUNDIALS_ROOT "" CACHE PATH
  "Prefix of a preinstalled SUNDIALS (expects include/ and lib/ under this directory)"
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

include_directories(src)

set(SOURCES
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
  src/integrators/rk3.cpp
  src/integrators/rk3_si.cpp
  src/integrators/cvode.cpp
  src/integrators/compute.cpp
  src/utils.h
)

add_executable(ruptura src/main.cpp ${SOURCES})

target_compile_options(ruptura PRIVATE
  -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align
  -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wno-reorder-ctor
  -Wnull-dereference -Wdouble-promotion -Wformat=2
  -O3 -march=native -ffast-math -fomit-frame-pointer -ftree-vectorize -funroll-loops
)

target_compile_definitions(ruptura PRIVATE
  BUILD_SUNDIALS=$<BOOL:${BUILD_SUNDIALS}>
)

target_include_directories(ruptura PRIVATE
  src
  src/integrators
  $<$<BOOL:${BUILD_SUNDIALS}>:${KLU_INCLUDE_DIR}>
)

target_link_libraries(ruptura PRIVATE
  ${LAPACK_LIBRARIES}
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::cvode>
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::nvecserial>
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::sunlinsolklu>
  $<$<BOOL:${BUILD_SUNDIALS}>:SUNDIALS::sunmatrixsparse>
)

add_subdirectory(docs)
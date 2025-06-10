cmake_minimum_required(VERSION 3.14)

project(ruptura LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SUNDIALS COMPONENTS cvode nvecserial nvecmanyvector
        HINTS "${SUNDIALS_ROOT}" REQUIRED)

# Include directories
include_directories(src)

# Compiler flags
set(CXX_COMPILE_FLAGS
        -g
        -O3
        -march=native
        -ffast-math
        -Wall
        -Wextra
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Woverloaded-virtual
        -Wpedantic
        -Wconversion
        -Wsign-conversion
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
        -fomit-frame-pointer
        -ftree-vectorize
        -fno-stack-check
        -funroll-loops
)

# Collect source files
set(SOURCES
        src/breakthrough.cpp
        src/breakthrough_state.cpp
        src/component.cpp
        src/fitting.cpp
        src/inputreader.cpp
        src/isotherm.cpp
        src/mixture_prediction.cpp
        src/multi_site_isotherm.cpp
        src/random_numbers.cpp
        src/special_functions.cpp
        src/integrators/rk3.cpp
        src/integrators/cvode.cpp
        src/integrators/compute.cpp
        src/utils.h
)

# -------------------------------
# Build the ruptura executable
# -------------------------------
add_executable(ruptura src/main.cpp ${SOURCES})

target_compile_options(ruptura PRIVATE ${CXX_COMPILE_FLAGS})

target_include_directories(ruptura PRIVATE src src/integrators)

target_link_libraries(ruptura PRIVATE
        SUNDIALS::cvode
        SUNDIALS::nvecserial
)

# -------------------------------
# Doxygen Documentation
# -------------------------------
add_subdirectory(docs)

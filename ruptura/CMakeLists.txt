project(ruptura_python)

set(NANOBIND_SOURCE_DIR "" CACHE PATH
  "Path to a nanobind source checkout. Leave empty to use find_package()/pip-installed nanobind."
)

set(RUPTURA_PYTHON_MODULE_NAME "_ruptura" CACHE STRING
  "Python import name for the nanobind module"
)

set(RUPTURA_PYTHON_EXECUTABLE "" CACHE FILEPATH
  "Exact Python interpreter to use for the Python module"
)

set(DEV_MODULE Development.Module)

# Prefer the active virtualenv/conda environment when one is activated.
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

nanobind_add_module(${PROJECT_NAME}
  NOMINSIZE
  NB_STATIC
  pybindings.cpp
)

set_target_properties(${PROJECT_NAME} PROPERTIES
  OUTPUT_NAME "${RUPTURA_PYTHON_MODULE_NAME}"
)

target_compile_options(${PROJECT_NAME} PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(${PROJECT_NAME} PRIVATE
  MODULE_NAME=${RUPTURA_PYTHON_MODULE_NAME}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  rupturalib
)

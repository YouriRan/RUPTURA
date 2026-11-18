find_package(GTest CONFIG QUIET)

if(NOT GTest_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
    SYSTEM  # its headers are compiled with our warning flags otherwise
  )

  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

add_executable(ruptura_tests
  convection_tests.cpp
  example_inputreader_tests.cpp
  fluid_phase_tests.cpp
  humidity_isotherm_tests.cpp
  mixture_prediction_tests.cpp
  multibed_tests.cpp
  notebook_templates_tests.cpp
  pressure_swing_tests.cpp
  ros34pw2_tests.cpp
  single_adsorbent_equivalence_tests.cpp
  sorption_jacobian_tests.cpp
  temperature_law_tests.cpp
)

# A GoogleTest in /usr/local/include is not a system header for Apple Clang: its headers gave ~200 warnings.
target_compile_options(ruptura_tests PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
  $<$<CXX_COMPILER_ID:AppleClang,Clang>:--system-header-prefix=gtest/>
)

target_compile_definitions(ruptura_tests PRIVATE
  RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_link_libraries(ruptura_tests PRIVATE
  rupturalib
  GTest::gtest_main
)

include(GoogleTest)
# The tests write output files into their working directory, which must not be the source tree.
gtest_discover_tests(ruptura_tests
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

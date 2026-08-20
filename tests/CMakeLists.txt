find_package(GTest CONFIG QUIET)

if(NOT GTest_FOUND)
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.zip
  )

  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

add_executable(ruptura_tests
  example_inputreader_tests.cpp
  fluid_phase_tests.cpp
  mixture_prediction_tests.cpp
  multibed_tests.cpp
  notebook_templates_tests.cpp
  single_adsorbent_equivalence_tests.cpp
)

target_compile_options(ruptura_tests PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(ruptura_tests PRIVATE
  RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_link_libraries(ruptura_tests PRIVATE
  rupturalib
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(ruptura_tests
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

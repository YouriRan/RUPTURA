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
  convection_tests.cpp
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

# A small method harness independent of column physics and SUNDIALS.
add_executable(ruptura_ros3_tests
  ros3_tests.cpp
  ../src/integrators/ros3.cpp
)
target_include_directories(ruptura_ros3_tests PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_compile_options(ruptura_ros3_tests PRIVATE ${RUPTURA_COMPILE_OPTIONS} -fno-fast-math)
target_link_libraries(ruptura_ros3_tests PRIVATE ${LAPACK_LIBRARIES} GTest::gtest_main)
gtest_discover_tests(ruptura_ros3_tests)

if(TARGET ruptura_lab_support)
  add_executable(ruptura_lab_tests lab_input_tests.cpp)
  target_link_libraries(ruptura_lab_tests PRIVATE ruptura_lab_support GTest::gtest)
  target_compile_definitions(ruptura_lab_tests PRIVATE RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  if(TARGET ruptura)
    target_compile_definitions(ruptura_lab_tests PRIVATE LAB_CLI_EXECUTABLE="$<TARGET_FILE:ruptura>")
  endif()
  add_test(NAME RupturaLab.InputCompatibility COMMAND ruptura_lab_tests)
  set_tests_properties(RupturaLab.InputCompatibility PROPERTIES ENVIRONMENT "QT_QPA_PLATFORM=offscreen")
endif()

find_package(benchmark CONFIG REQUIRED)

add_executable(multibed_single_adsorbent_benchmark
  benchmark_cases.h
  single_adsorbent_benchmark.cpp
)

add_executable(multibed_single_adsorbent_results
  benchmark_cases.h
  compare_results.cpp
)

target_compile_options(multibed_single_adsorbent_benchmark PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(multibed_single_adsorbent_benchmark PRIVATE
  RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_compile_options(multibed_single_adsorbent_results PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(multibed_single_adsorbent_results PRIVATE
  RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_link_libraries(multibed_single_adsorbent_benchmark PRIVATE
  rupturalib
  benchmark::benchmark
  benchmark::benchmark_main
)

target_link_libraries(multibed_single_adsorbent_results PRIVATE
  rupturalib
)

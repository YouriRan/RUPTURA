find_package(Qt6 6.5 REQUIRED COMPONENTS
  Network
  Widgets
)

set(CMAKE_AUTOMOC ON)

add_executable(ruptura_lab
  main.cpp
  mainwindow.cpp
  mainwindow.h
  runner.cpp
  runner.h
)

target_compile_options(ruptura_lab PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(ruptura_lab PRIVATE
  RUPTURA_CLI_EXECUTABLE="$<TARGET_FILE:${PROJECT_NAME}>"
  RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_link_libraries(ruptura_lab PRIVATE
  Qt6::Network
  Qt6::Widgets
)

add_dependencies(ruptura_lab ${PROJECT_NAME})

install(TARGETS ruptura_lab DESTINATION bin COMPONENT app)

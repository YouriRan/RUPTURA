if(BUILD_APP)
  set(RUPTURA_LAB_PYTHON_EXECUTABLE "" CACHE FILEPATH "Default Python interpreter for the Ruptura Lab web launcher")
  add_executable(ruptura_lab webmain.cpp)
  target_compile_definitions(ruptura_lab PRIVATE RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  if(RUPTURA_LAB_PYTHON_EXECUTABLE)
    target_compile_definitions(ruptura_lab PRIVATE RUPTURA_LAB_PYTHON_EXECUTABLE="${RUPTURA_LAB_PYTHON_EXECUTABLE}")
  elseif(BUILD_PYTHON AND Python_EXECUTABLE)
    target_compile_definitions(ruptura_lab PRIVATE RUPTURA_LAB_PYTHON_EXECUTABLE="${Python_EXECUTABLE}")
  endif()
  if(TARGET ruptura)
    target_compile_definitions(ruptura_lab PRIVATE RUPTURA_CLI_EXECUTABLE="$<TARGET_FILE:ruptura>")
    add_dependencies(ruptura_lab ruptura)
  endif()
  install(TARGETS ruptura_lab DESTINATION bin COMPONENT app)
  install(DIRECTORY "${PROJECT_SOURCE_DIR}/ruptura/" DESTINATION share/ruptura/python/ruptura
    COMPONENT app FILES_MATCHING PATTERN "*.py" PATTERN "*.html" PATTERN "*.css" PATTERN "*.js" PATTERN "LICENSE.txt"
    PATTERN "*.png" PATTERN "*.json"
    PATTERN "__pycache__" EXCLUDE)
endif()

if(BUILD_QT_APP)
find_package(Qt6 6.5 REQUIRED COMPONENTS
  Network
  Widgets
)

set(CMAKE_AUTOMOC ON)

add_library(ruptura_lab_support STATIC
  mainwindow.cpp
  mainwindow.h
  runner.cpp
  runner.h
  simulationinput.cpp
  simulationinput.h
)

target_compile_options(ruptura_lab_support PRIVATE
  ${RUPTURA_COMPILE_OPTIONS}
)

target_compile_definitions(ruptura_lab_support PRIVATE
  RUPTURA_SOURCE_DIR="${PROJECT_SOURCE_DIR}"
)

target_include_directories(ruptura_lab_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ruptura_lab_support PUBLIC
  rupturalib
  Qt6::Network
  Qt6::Widgets
)

if(TARGET ruptura)
  target_compile_definitions(ruptura_lab_support PRIVATE RUPTURA_CLI_EXECUTABLE="$<TARGET_FILE:ruptura>")
  add_dependencies(ruptura_lab_support ruptura)
endif()

add_executable(ruptura_lab_qt main.cpp)
target_link_libraries(ruptura_lab_qt PRIVATE ruptura_lab_support)

install(TARGETS ruptura_lab_qt DESTINATION bin COMPONENT app)
endif()

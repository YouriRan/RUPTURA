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

find_program(DOXYGEN_PATH doxygen)

if(NOT DOXYGEN_PATH)
    message(STATUS
        "Doxygen not found; the 'documentation' target will not be available. "
        "Install Doxygen and re-run CMake to enable it."
    )
    return()
endif()

message(STATUS "Doxygen found: ${DOXYGEN_PATH}")

foreach(_css doxygen-awesome.css doxygen-awesome-sidebar-only.css)
    if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/${_css}")
        file(DOWNLOAD
            "https://raw.githubusercontent.com/jothepro/doxygen-awesome-css/v2.3.1/${_css}"
            "${CMAKE_CURRENT_LIST_DIR}/${_css}"
            STATUS _css_status
        )
        list(GET _css_status 0 _css_code)
        if(NOT _css_code EQUAL 0)
            list(GET _css_status 1 _css_message)
            message(WARNING "Could not download ${_css}: ${_css_message}")
            file(REMOVE "${CMAKE_CURRENT_LIST_DIR}/${_css}")
        endif()
    endif()
endforeach()

add_custom_target(documentation
    COMMAND ${DOXYGEN_PATH} ${CMAKE_CURRENT_LIST_DIR}/Doxyfile
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    BYPRODUCTS ${CMAKE_BINARY_DIR}/html/index.html
)

add_custom_command(TARGET documentation POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo
            "Documentation successfully generated. You can preview at: ${CMAKE_BINARY_DIR}/html/index.html"
)

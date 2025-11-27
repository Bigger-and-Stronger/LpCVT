########################################################################################################################
# main
########################################################################################################################
add_executable(${PROJECT_NAME}_MAIN
        main.cpp
)

target_include_directories(${PROJECT_NAME}_MAIN
        PRIVATE
        "${CMAKE_SOURCE_DIR}/deps"
        "${CMAKE_SOURCE_DIR}/src"
)

target_link_libraries(${PROJECT_NAME}_MAIN
        PRIVATE
        ${PROJECT_NAME}_COMMON
        Geogram::geogram
)
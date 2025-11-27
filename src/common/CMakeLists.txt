add_library(${PROJECT_NAME}_COMMON)

target_include_directories(${PROJECT_NAME}_COMMON
        PRIVATE
        "${CMAKE_SOURCE_DIR}/src"
        "${CMAKE_SOURCE_DIR}/deps"
)

# configure sources
file(GLOB COMMON_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cxx)

# add sources to the library
target_sources(${PROJECT_NAME}_COMMON
        PRIVATE
        ${COMMON_SOURCE}
)
target_link_libraries(${PROJECT_NAME}_COMMON
        PRIVATE
)
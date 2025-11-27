# this CMakeLists.txt is copy from https://github.com/secantalpha/LpCVT/blob/main/sources/LpCVT/CMakeLists.txt, and make some modifies

########################################################################################################################
# CMake options
########################################################################################################################
cmake_minimum_required(VERSION 3.30)

project(LpCVT)

if(WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /utf-8")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

########################################################################################################################
# CMake module path
########################################################################################################################
list(PREPEND CMAKE_MODULE_PATH
        ${CMAKE_CURRENT_LIST_DIR}/cmake)

########################################################################################################################
# Set variables
########################################################################################################################
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin) # Set output path for the compiled files
set(GEOGRAM_DIR         "$ENV{GEOGRAM_DIR}" CACHE PATH "Directory where Geogram is installed")
set(OPTION_BUILD_TESTS  ON CACHE BOOL "Build tests (based on GoogleTest).")

########################################################################################################################
# Include 3rdParty
########################################################################################################################
include(Dependencies)

########################################################################################################################
# Add subdirectories
########################################################################################################################
add_subdirectory(src)
add_subdirectory(app)

message(STATUS "ENABLE_BUILD_TESTS = ${OPTION_BUILD_TESTS}")
if(OPTION_BUILD_TESTS)
    add_subdirectory(tests)
endif()
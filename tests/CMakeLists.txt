############################################################################################
# Download and install googletest
############################################################################################
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

############################################################################################
# Build tests
############################################################################################
enable_testing()

file(GLOB TESTS_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)

add_executable(${PROJECT_NAME}_TESTS
        ${TESTS_SOURCE})

target_include_directories(${PROJECT_NAME}_TESTS
        PRIVATE
        "${CMAKE_SOURCE_DIR}/src"
        "${CMAKE_SOURCE_DIR}/deps"
)

target_link_libraries(${PROJECT_NAME}_TESTS
        ${PROJECT_NAME}_BASE
        ${PROJECT_NAME}_UTILS
        ${PROJECT_NAME}_IO
        ${PROJECT_NAME}_EVAL
        ${PROJECT_NAME}_COMMON
        ${PROJECT_NAME}_CORE
        Geogram::geogram
        GTest::gtest_main
)

add_compile_definitions(TEST_DATA_PATH="${CMAKE_CURRENT_SOURCE_DIR}/testdata/")

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}_TESTS)
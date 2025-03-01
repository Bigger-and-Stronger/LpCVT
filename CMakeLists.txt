# this CMakeLists.txt is copy from https://github.com/secantalpha/LpCVT/blob/main/sources/LpCVT/CMakeLists.txt, and make some modifies

cmake_minimum_required(VERSION 3.21)
project(LpCVT)

find_package(boost REQUIRED)
find_package(CGAL 5.0 QUIET REQUIRED COMPONENTS Core)

file(GLOB_RECURSE SOURCES 
    ${PROJECT_SOURCE_DIR}/src/LpCVT/*.h
    ${PROJECT_SOURCE_DIR}/src/LpCVT/*.cpp)

add_executable(LpCVT ${SOURCES})

target_include_directories(LpCVT PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_include_directories(LpCVT PRIVATE ${GMP_INCLUDE_DIR})

target_link_libraries(LpCVT PRIVATE CGAL::CGAL)

set(GMP_DLL "${GMP_INCLUDE_DIR}/../lib/libgmp-10.dll")
add_custom_command(
  TARGET LpCVT
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GMP_DLL}
          $<TARGET_FILE_DIR:LpCVT>)
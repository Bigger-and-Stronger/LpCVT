# this CMakeLists.txt is copy from https://github.com/secantalpha/LpCVT/blob/main/sources/LpCVT/CMakeLists.txt, and make some modifies

cmake_minimum_required(VERSION 3.21)
project(LpCVT)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      "${CMAKE_CURRENT_LIST_DIR}/Release"
      CACHE PATH "default install path" FORCE)
endif()

find_package(Boost REQUIRED)
find_package(CGAL QUIET REQUIRED COMPONENTS Core)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_LIST_DIR}/src/*.h
     ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)

include_directories(${CMAKE_CURRENT_LIST_DIR}/src)

add_executable(LpCVT ${SOURCES})

cmake_path(GET CMAKE_CURRENT_LIST_DIR PARENT_PATH LP_CVT_INCLUDE_DIR)

target_include_directories(LpCVT PRIVATE ${LP_CVT_INCLUDE_DIR})

target_link_libraries(LpCVT PRIVATE CGAL::CGAL)

set(GMP_DLL "$ENV{GMP_DIR}/lib/libgmp-10.dll")

add_custom_command(
  TARGET LpCVT
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${GMP_DLL}
          $<TARGET_FILE_DIR:LpCVT>)

install(TARGETS LpCVT DESTINATION .)
install(FILES ${GMP_DLL} DESTINATION .)
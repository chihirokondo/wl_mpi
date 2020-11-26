cmake_minimum_required(VERSION 3.1)
project(EXACT CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

add_executable(exact_ferro_ising.app exact_ferro_ising.cpp)

cmake_minimum_required(VERSION 3.1)
project(REWL_MPI CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

find_package(MPI)
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
endif(MPI_FOUND)
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})

add_subdirectory(src)

add_executable(main_mpi.app main_mpi.cpp)
target_link_libraries(main_mpi.app ${MPI_CXX_LIBRARIES} rewl_lib)

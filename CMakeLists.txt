cmake_minimum_required(VERSION 3.1)
project(REWL_MPI CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

## OpenMPI
find_package(MPI)
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
endif(MPI_FOUND)
## eigen3
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})
## JSON for Modern C++
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(include/nlohmann_json_cmake_fetchcontent)

add_executable(sample.app sample.cpp)
target_link_libraries(sample.app
  ${MPI_CXX_LIBRARIES}
  nlohmann_json::nlohmann_json
)

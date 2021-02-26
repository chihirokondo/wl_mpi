cmake_minimum_required(VERSION 3.1)
project(REWL_MPI_BENCHMARKS CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

add_executable(benchmark.app benchmark.cpp)
add_executable(ferro_ising_square.app ferro_ising_square.cpp)

## OpenMPI
find_package(MPI)
if(MPI_FOUND)
  target_include_directories(benchmark.app PUBLIC ${MPI_INCLUDE_PATH})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
endif(MPI_FOUND)
## eigen3
find_package(Eigen3 REQUIRED)
include(${EIGEN3_USE_FILE})
## boost
find_package(Boost REQUIRED)
target_include_directories(benchmark.app PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(ferro_ising_square.app PRIVATE ${Boost_INCLUDE_DIRS})
## JSON for Modern C++
set(JSON_BuildTests OFF CACHE INTERNAL "")
add_subdirectory(../include/nlohmann_json_cmake_fetchcontent nlohmann_json)
## WL_MPI and model_sample
set(WL_MPI_INCLUDE_PATH "../")
target_include_directories(benchmark.app PUBLIC ${WL_MPI_INCLUDE_PATH})

target_link_libraries(benchmark.app
  ${MPI_CXX_LIBRARIES}
  nlohmann_json::nlohmann_json
  )
target_link_libraries(ferro_ising_square.app
  nlohmann_json::nlohmann_json
  )

cmake_minimum_required(VERSION 3.1)

add_library(rewl_lib STATIC
  rewl.cpp
  histo_env_manager.cpp
  mpi_setting.cpp
  stop_callback.cpp
  window.cpp
)
target_include_directories(rewl_lib
  PUBLIC ../include
)

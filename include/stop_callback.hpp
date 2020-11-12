#ifndef WANGLANDAU_STOP_CALLBACK_H_
#define WANGLANDAU_STOP_CALLBACK_H_


#include <iostream>
#include <chrono>
#include "mpi_setting.hpp"


class StopCallback {
 public:
  StopCallback(MPIV mpiv, double timelimit_secs)
      : mpiv_(mpiv),
        start_(std::chrono::system_clock::now()),
        timelimit_secs_(timelimit_secs) {}
  bool operator()() const;
 private:
  const MPIV mpiv_;
  const std::chrono::system_clock::time_point start_;
  const double timelimit_secs_;
};


#endif // WANGLANDAU_STOP_CALLBACK_H_

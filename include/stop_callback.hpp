#ifndef WANGLANDAU_STOP_CALLBACK_H_
#define WANGLANDAU_STOP_CALLBACK_H_


#include <iostream>
#include <chrono>
#include "mpi.hpp"


class StopCallback {
 public:
  StopCallback(MPIV mpiv, double timelimit)
      : mpiv_(mpiv),
        start_(std::chrono::system_clock::now()),
        timelimit_(timelimit) {}
  bool operator()() const;
 private:
  const MPIV mpiv_;
  const std::chrono::system_clock::time_point start_;
  const double timelimit_;
};


bool StopCallback::operator()() const {
  bool to_stop;
  if (mpiv_.myid() == 0) {
    std::chrono::system_clock::time_point end;
    end = std::chrono::system_clock::now();
    to_stop =
        std::chrono::duration_cast<std::chrono::seconds>(end-start_).count() >=
        timelimit_;
  }
  MPI_Bcast(&to_stop, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  return to_stop;
}


#endif // WANGLANDAU_STOP_CALLBACK_H_
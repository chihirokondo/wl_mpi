#include "stop_callback.hpp"


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
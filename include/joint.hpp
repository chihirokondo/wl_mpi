#ifndef WANGLANDAU_JOINT_H_
#define WANGLANDAU_JOINT_H_


#include <mpi.h>
#include <cmath>
#include <vector>
#include "mpi_setting.hpp"
#include "window.hpp"


inline void joint_ln_dos(std::vector<double> *ln_dos_ptr,
    const WindowManager &window, const MPIV &mpiv);
inline int search_joint_point(int imin, int imax,
    const std::vector<double> &ln_dos,
    const std::vector<double> &ln_dos_next_window);


void joint_ln_dos(std::vector<double> *ln_dos_ptr, const WindowManager &window,
    const MPIV &mpiv) {
  MPI_Status status;
  std::vector<double> &ln_dos(*ln_dos_ptr);
  if (mpiv.myid() == 0) {
    // Boss process joints individual "ln_dos"s in worker process.
    int imin_next_window;
    int imax_this_window = window.imax();
    std::vector<double> ln_dos_next_window(ln_dos.size(), 0.0);
    for (int i=1; i<mpiv.num_windows(); ++i) {
      // Recieve next window's information and "ln_dos".
      MPI_Recv(&imin_next_window, 1, MPI_INT, i*mpiv.num_walkers_window(), 101,
          MPI_COMM_WORLD, &status);
      MPI_Recv(&ln_dos_next_window[0], ln_dos.size(), MPI_DOUBLE,
          i*mpiv.num_walkers_window(), 202, MPI_COMM_WORLD, &status);
      // Search joint point.
      int ijoint = search_joint_point(imin_next_window, imax_this_window,
          ln_dos, ln_dos_next_window);
      // Joint "ln_dos"s at "ijoint".
      double delta = ln_dos[ijoint] - ln_dos_next_window[ijoint];
      MPI_Recv(&imax_this_window, 1, MPI_INT, i*mpiv.num_walkers_window(), 303,
          MPI_COMM_WORLD, &status);
      for (int i=ijoint; i<=imax_this_window; ++i) {
        if (ln_dos_next_window[i] == 0.0) {
          ln_dos[i] = 0.0;
        } else {
          ln_dos[i] = delta + ln_dos_next_window[i];
        }
      }
    }
  } else if (mpiv.myid()%mpiv.num_walkers_window() == 0) {
    // Worker processes send their information and "ln_dos" to boss process.
    int imin = window.imin();
    int imax = window.imax();
    MPI_Send(&imin, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
    MPI_Send(&ln_dos[0], ln_dos.size(), MPI_DOUBLE, 0, 202, MPI_COMM_WORLD);
    MPI_Send(&imax, 1, MPI_INT, 0, 303, MPI_COMM_WORLD);
  }
}


int search_joint_point(int imin, int imax, const std::vector<double> &ln_dos,
    const std::vector<double> &ln_dos_next_window) {
  // Search index of joint point by comparing gradient.
  int iprevious = imin;
  double min_grad_diff;
  bool first_min_candidate = true;
  int ijoint;
  for (int i=iprevious+1; i<=imax; ++i) {
    if (ln_dos[iprevious] == 0.0) {
      ++iprevious;
    } else if (ln_dos[i] != 0.0) {
      // Calculate gradient.
      double grad = (ln_dos[i]-ln_dos[iprevious])/(i-iprevious);
      double grad_next_window =
          (ln_dos_next_window[i]-ln_dos_next_window[iprevious])/(i-iprevious);
      if ((min_grad_diff>std::fabs(grad-grad_next_window)) ||
          first_min_candidate) {
        min_grad_diff = std::fabs(grad-grad_next_window);
        ijoint = iprevious;
        first_min_candidate = false;
      }
      iprevious = i;
    }
  }
  return ijoint;
}


#endif // WANGLANDAU_JOINT_H_
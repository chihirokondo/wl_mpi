#include <mpi.h>
#include <cmath>
#include <random>
#include <iostream>
#include <vector>
#include "mpi_setting.hpp"
#include "histo_env_manager.hpp"
#include "window.hpp"
#include "wl_params.hpp"
#include "rewl.hpp"


int generate_partner(std::mt19937 &engine, int exchange_pattern,
    const MPIV &mpiv) {
  std::vector<size_t> partner_list(2*mpiv.multiple());
  int partner;
  // 'head-node' in the window determines pairs of flippartners.
  if (mpiv.local_id(exchange_pattern) == 0) {
    int choose_from = mpiv.multiple();
    int select;
    std::vector<size_t> lib_re(mpiv.multiple());
    for (size_t i=0; i<mpiv.multiple(); ++i) lib_re[i] = mpiv.multiple()+i;
    for (size_t i=0; i<mpiv.multiple(); ++i) {
      std::uniform_int_distribution<> dist(0, choose_from-1);
      select = dist(engine);
      partner_list[i] = lib_re[select];
      partner_list[lib_re[select]] = i;
      --choose_from;
      for (size_t j=select; j<choose_from; ++j) lib_re[j] = lib_re[j+1];
    }
  }
  // At this point, every walker has a swap partner assigned,
  // now they must be communicated.
  if (mpiv.comm_id() != -1) {
    MPI_Scatter(&partner_list[0], 1, MPI_INT, &partner, 1, MPI_INT, 0,
        mpiv.local_comm(mpiv.comm_id()));
  } else {
    partner = -1;
  }
  return partner;
}


bool check_histoflat(const WindowManager &window,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv) {
  // Check flatness of the histogram for all walkers in the window.
  MPI_Status status;
  bool my_flat = true;
  bool other_flat;
  int num_bins = 0;
  int min_histo = histogram[window.imin()];
  double average = 0.0;
  for (size_t i=window.imin(); i<=window.imax(); ++i) {
    if (histogram[i] != 0) {
      ++num_bins;
      average += (double)histogram[i];
      if ((min_histo == 0.0) | (histogram[i]<min_histo)) {
        min_histo = histogram[i];
      }
    }
  }
  average /= num_bins;
  if ((double)min_histo < flatness*average) my_flat = false;
  // Now talk to all the other walkers in the window.
  if (mpiv.myid()%mpiv.multiple() == 0) {
    // 'root' in window, receive individual flatnesses.
    for (int i=1; i<mpiv.multiple(); ++i) {
      MPI_Recv(&other_flat, 1, MPI_CXX_BOOL, mpiv.myid()+i, 66, MPI_COMM_WORLD,
          &status);
      my_flat *= other_flat;
    }
    for (int i=1; i<mpiv.multiple(); ++i) {
      // Let everybody know.
      MPI_Send(&my_flat, 1, MPI_CXX_BOOL, mpiv.myid()+i, 88, MPI_COMM_WORLD);
    }
  } else {
    // Send individual flatness and receive 'merged' flatness.
    MPI_Send(&my_flat, 1, MPI_CXX_BOOL,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 66, MPI_COMM_WORLD);
    MPI_Recv(&other_flat, 1, MPI_CXX_BOOL,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 88, MPI_COMM_WORLD, &status);
    my_flat = other_flat;  // Replace individual flatness by merged.
  }
  return my_flat;
  // Note: By now, myflat refers to the 'collective' flatness in the window,
  //       not the flatness of an individual walker.
}


void merge_ln_dos(std::vector<double> *ln_dos_ptr, const MPIV &mpiv) {
  MPI_Status status;
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<double> ln_dos_buf = ln_dos;
  if (mpiv.myid()%mpiv.multiple() == 0) {
    // 'root' in window, receive individual g(E) and send merged g(E).
    for (int i=1; i<mpiv.multiple(); ++i) {
      MPI_Recv(&ln_dos_buf[0], ln_dos.size(), MPI_DOUBLE, mpiv.myid()+i, 77,
          MPI_COMM_WORLD, &status);
      for (size_t j=0; j<ln_dos.size(); ++j) ln_dos[j] += ln_dos_buf[j];
    }
    int num_bins = 0;
    double mean = 0.0;
    for (size_t i=0; i<ln_dos.size(); ++i) {
      ln_dos[i] /= mpiv.multiple();
      mean += ln_dos[i];
      if (ln_dos[i] != 0.0) {
        ++num_bins;
      }
    }
    mean /= num_bins;
    for (double &ln_dos_i : ln_dos) {
      if (ln_dos_i != 0.0) {
        ln_dos_i -= mean;
      }
    }
    for (int i=1; i<mpiv.multiple(); ++i) {
      MPI_Send(&ln_dos[0], ln_dos.size(), MPI_DOUBLE, mpiv.myid()+i, 99,
          MPI_COMM_WORLD);
    }
  } else {
    // Send individual g(E) and receive merged g(E).
    MPI_Send(&ln_dos[0], ln_dos.size(), MPI_DOUBLE,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 77, MPI_COMM_WORLD);
    MPI_Recv(&ln_dos[0], ln_dos.size(), MPI_DOUBLE,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 99, MPI_COMM_WORLD,
        &status);
  }
}

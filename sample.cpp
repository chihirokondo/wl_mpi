#include <mpi.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "include/mpi_setting.hpp"
#include "include/histo_env_manager.hpp"
#include "include/window.hpp"
#include "include/wl_params.hpp"
#include "include/rewl.hpp"
// Model
#include "model_sample/lattice/graph.hpp"
#include "model_sample/ferro_ising.hpp"


int main(int argc, char *argv[]) {
  int numprocs, myid, num_walkers_window;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // Check command line arguments.
  try {
    if (argc != 7) throw 0;
  }
  catch (int err_status) {
    if (myid == 0) {
      std::cerr
          << "ERROR: Unexpected number of command line arguments!\n"
          << "       Expect 6 arguments, " << argc - 1 << " were provided.\n"
          << "Syntax: ./a.out [arg1] [arg2] [arg3] [arg4] [arg5] [arg6] \n\n"
          << "Please provide the following command line arguments:\n"
          << "1. Overlap between consecutive windows."
          << " [double, 0 <= overlap <= 1]\n"
          << "2. Number of walkers per window. [integer]\n"
          << "3. Number of Monte Carlo steps between replica exchange."
          << " [integer]\n"
          << "4. Random number seed. [integer]\n"
          << "5. Time limit (secs). [double]\n"
          << "6. Should execute from the top. [integer (bool)]\n"
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  num_walkers_window = atoi(argv[2]);
  try {
    // "numprocs" must be a multiple of "num_walkers_window".
    if (numprocs%num_walkers_window != 0) throw 0;
  }
  catch (int err_status) {
    if (myid == 0) {
      std::cerr
          << "ERROR: Number of processes must be a multiple of"
          << " the second command line argument.\n"
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPIV mpiv(numprocs, myid, num_walkers_window);
  // Model dependent variables.
  int dim = 2;
  int length = 4;
  lattice::graph lat = lattice::graph::simple(dim, length);
  double condition_value = std::pow(2.0, (double)lat.num_sites());
  int sweeps = lat.num_sites();
  FerroIsing model(lat);
  HistoEnvManager histo_env(model.ene_min(), model.ene_max(), model.num_bins(),
      true);
  // Original Wang-Landau parameters.
  int check_flatness_every = 500;
  double lnf = 1.0;
  double lnfmin = 1e-8;
  double flatness = 0.95;
  std::vector<double> ln_dos(histo_env.num_bins(), 0.0);
  std::mt19937 engine(atoi(argv[4])+mpiv.myid());
  // Replica exchange Wang-Landau (REWL) parameters.
  double overlap = atof(argv[1]);
  int swap_every = atoi(argv[3]);
  WLParams wl_params(sweeps, check_flatness_every, lnf, lnfmin, flatness,
      swap_every);
  // Settings for the windows.
  WindowManager window(histo_env, mpiv, overlap);
  int width_index = window.iwidth();
  int width_index_min;
  MPI_Allreduce(&width_index, &width_index_min, 1, MPI_INT, MPI_MIN,
      MPI_COMM_WORLD);
  try {
    if (width_index_min < 1) throw 0;
  }
  catch (int err_status) {
    if (mpiv.myid() == 0) {
      std::cerr
          << "ERROR: Too many number of the windows or"
          << " too few number of the bins were provided.\n"
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // Program control.
  double timelimit_secs = atof(argv[5]);
  bool from_the_top = atoi(argv[6]);
  int running_state;
  // REWL routine.
  running_state = rewl<FerroIsing>(&ln_dos, &model, histo_env, &wl_params,
      window, &mpiv, engine, timelimit_secs, from_the_top);
  if (running_state == 1) {
    // Output.
    merge_ln_dos(&ln_dos, mpiv);
    if (mpiv.myid()%mpiv.num_walkers_window() == 0) {
      std::string filename = "./rawdata/ln_val_window" +
          std::to_string(mpiv.myid()/mpiv.num_walkers_window()) + ".dat";
      std::ofstream ofs(filename, std::ios::out);
      ofs << "# dim: " << dim << ", length: " << length << "\n";
      ofs << "# condition_type: " << "sum" << "\n";
      ofs << "# condition_value: " << condition_value << "\n";
      ofs << "# sewing_point: " << window.isew() << "\n";
      ofs << "# number_of_windows: " << mpiv.num_windows() << "\n";
      ofs << "# value \t # lngV\n";
      for (size_t i=window.imin(); i<=window.imax(); ++i) {
        ofs << histo_env.GetVal(i, "mid") << "\t"
            << std::scientific << std::setprecision(15) << ln_dos[i] << "\n";
      }
      ofs << std::endl;
    }
  } else if (running_state == -1) {
    if (mpiv.myid() == 0) {
      std::cerr
          << "ERROR: Cannot restart the experiment.\n"
          << "       Last-time job was completely finished or "
          << "some conditions have been changed."
          << std::endl;
    }
  }
  MPI_Finalize();
  return 0;
  //// result state.
}

#include <cmath>
#include <mpi.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "include/wl_mpi.hpp"
// Sample model (classical ferro magnetic ising model)
#include "model_sample/lattice/graph.hpp"
#include "model_sample/ferro_ising.hpp"


int main(int argc, char *argv[]) {
  int numprocs, myid, num_walkers_window;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // Check command line arguments.
  if (argc != 5)  {
    if (myid == 0) {
      std::cerr
          << "ERROR: Unexpected number of command line arguments!\n"
          << "       Expect 6 arguments, " << argc - 1 << " were provided.\n"
          << "Syntax: " << argv[0]
          << " [arg1] [arg2] [arg3] [arg4]  \n\n"
          << "Please provide the following command line arguments:\n"
          << "1. Number of walkers per window. [integer]\n"
          << "2. Random number seed. [integer]\n"
          << "3. Time limit (secs). [double]\n"
          << "4. Should execute from the top. [integer (bool)]\n"
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  num_walkers_window = atoi(argv[1]);
  MPIV mpiv(numprocs, myid, num_walkers_window);
  // Model dependent variables.
  int dim = 2;
  int length = 4;
  lattice::graph lat = lattice::graph::simple(dim, length);
  FerroIsing model(lat);
  HistoEnvManager histo_env(model.ene_min(), model.ene_max(), model.num_bins(),
      true);
  // Replica exchange Wang-Landau (REWL) parameters.
  int check_flatness_every = 500;
  double lnf = 1.0;
  double lnfmin = 1e-8;
  double flatness = 0.95;
  double overlap = 0.75; // 0<= overlap <= 1.
  int exch_every = 100;
  WLParams wl_params(check_flatness_every, lnf, lnfmin, flatness, overlap,
      exch_every);
  std::mt19937 engine(atoi(argv[2])+mpiv.myid());
  // Program control variables.
  double timelimit_secs = atof(argv[3]);
  bool from_the_top = atoi(argv[4]);
  // REWL routine.
  std::vector<double> ln_dos;
  RunningState running_state = rewl<FerroIsing>(&ln_dos, &model, histo_env,
      &wl_params, &mpiv, engine, timelimit_secs, from_the_top);
  int result = 0;
  if (running_state == RunningState::ERROR) {
    result = static_cast<int>(RunningState::ERROR);
  }
  if ((running_state==RunningState::ALL_FINISHED) && (mpiv.myid()==0)) {
    double ln_const = ln_dos[0] - std::log(2.0);
    for (auto &ln_dos_i : ln_dos) {
      if (ln_dos_i != 0.0) ln_dos_i -= ln_const;
    }
    std::ofstream ofs("ln_dos_jointed.dat", std::ios::out);
    ofs << "# ferro ising model\n";
    ofs << "# dim = " << dim << ", length = " << length << "\n";
    ofs << "# energy\t # log_e (density of states)\n";
    for (size_t i=0; i<ln_dos.size(); ++i) {
      ofs << std::scientific << std::setprecision(15)
          << histo_env.GetVal(i, "mid") << "\t" << ln_dos[i] << "\n";
    }
    ofs << std::endl;
  }
  MPI_Finalize();
  return result;
}

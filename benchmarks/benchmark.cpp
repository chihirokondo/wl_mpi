#include <cmath>
#include <mpi.h>
#include <random>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>
#include "boost/date_time/posix_time/posix_time.hpp"
#include <nlohmann/json.hpp>
#include "include/wl_mpi.hpp"
// Sample model (classical ferro magnetic ising model)
#include "model_sample/lattice/graph.hpp"
#include "model_sample/ferro_ising.hpp"
using json = nlohmann::json;
using namespace boost::posix_time;


int main(int argc, char *argv[]) {
  ptime now = second_clock::local_time();
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
  int length = 8;
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
  std::chrono::system_clock::time_point start =
      std::chrono::system_clock::now();
  RunningState running_state = rewl<FerroIsing>(&ln_dos, &model, histo_env,
      &wl_params, &mpiv, engine, timelimit_secs, from_the_top);
  std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
  double elapsed_time =
      std::chrono::duration_cast<std::chrono::seconds>(end-start).count();
  MPI_Allreduce(MPI_IN_PLACE, &elapsed_time, 1, MPI_DOUBLE, MPI_MAX,
      MPI_COMM_WORLD);
  int result = 0;
  if (running_state == RunningState::ERROR) {
    result = static_cast<int>(RunningState::ERROR);
  }
  if ((running_state==RunningState::ALL_FINISHED) && (mpiv.myid()==0)) {
    // Normalization.
    double ln_const = ln_dos[0] - std::log(2.0);
    for (auto &ln_dos_i : ln_dos) {
      if (ln_dos_i != 0.0) ln_dos_i -= ln_const;
    }
    json data_wl_json;
    // Mapping.
    data_wl_json["mode"] = "test";
    data_wl_json["model"] = "classical ferro ising";
    data_wl_json["lattice_type"] = "square";
    data_wl_json["lattice_dim"] = dim;
    data_wl_json["lattice_length_x"] = length;
    data_wl_json["lattice_length_y"] = length;
    data_wl_json["normalization"] = "first_bin";
    data_wl_json["elapsed_time"] = elapsed_time;
    data_wl_json["num_procs"] = mpiv.numprocs();
    data_wl_json["num_procs_per_window"] = mpiv.num_walkers_window();
    data_wl_json["num_windows"] = mpiv.num_windows();
    data_wl_json["check_flatness_every"] = check_flatness_every;
    data_wl_json["initial_lnf"] = lnf;
    data_wl_json["lnf_min"] = lnfmin;
    data_wl_json["flatness"] = flatness;
    data_wl_json["overlap"] = overlap;
    data_wl_json["exch_every"] = exch_every;
    data_wl_json["random_seed_rank0"] = atoi(argv[2]);
    json json_ln_dos(ln_dos);
    data_wl_json["ln_dos"] = json_ln_dos; // TODO: check precision.
    std::vector<double> energies(ln_dos.size(), 0.0);
    for (size_t i=0; i<energies.size(); ++i) {
      energies[i] = histo_env.GetVal(i, "mid");
    }
    json json_energies(energies);
    data_wl_json["energies"] = json_energies;
    std::string datetime = to_iso_extended_string(now);
    std::string filename = "../analysis/data_wl/" + datetime + ".json";
    std::ofstream ofs(filename, std::ios::out);
    ofs << std::setw(4) << data_wl_json << std::endl;
  }
  MPI_Finalize();
  return result;
}

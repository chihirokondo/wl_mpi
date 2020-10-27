// Ferromagnetic Ising model.
#include <mpi.h>
#include <cmath>
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "include/mpi.hpp"
#include "include/lattice/graph.hpp"
#include "include/ferro_ising.hpp"
#include "include/histo_env_manager.hpp"
#include "include/window.hpp"
#include "include/wl_params.hpp"


template <typename Model>
void rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 &engine);
int generate_partner(std::mt19937 &engine, int exchange_pattern,
    const MPIV &mpiv);
bool check_histoflat(const WindowManager &window,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv);
template <typename Model>
bool replica_exchange(double *val_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const Model &model,
    const HistoEnvManager &histo_env, const WindowManager &window,
    const MPIV &mpiv, std::mt19937 &engine);
void merge_ln_dos(std::vector<double> *ln_dos_ptr, const MPIV &mpiv);


int main(int argc, char *argv[]) {
  int numprocs, myid, multiple;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  // Check command line arguments.
  try {
    if (argc != 5) throw 0;
  }
  catch (int err_status) {
    if (myid == 0) {
      std::cerr
          << "ERROR: Unexpected number of command line arguments!\n"
          << "       Expect 4 arguments, " << argc - 1 << " were provided.\n"
          << "Syntax: ./a.out [arg1] [arg2] [arg3] [arg4] \n\n"
          << "Please provide the following command line arguments:\n"
          << "1. Overlap between consecutive windows."
          << " [double, 0 <= overlap <= 1]\n"
          << "2. Number of walkers per window. [integer]\n"
          << "3. Number of Monte Carlo steps between replica exchange."
          << " [integer]\n"
          << "4. Random number seed. [integer]\n" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // mpiv.multiple() is # of walkers per window.
  multiple = atoi(argv[2]);
  try {
    // "numprocs" must be a multiple of "multiple".
    if (numprocs%multiple != 0) throw 0;
  }
  catch (int err_status) {
    if (myid == 0) {
      std::cerr
          << "ERROR: Number of processes must be a multiple of"
          << " the second command line argument.\n" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPIV mpiv(numprocs, myid, multiple);
  if (mpiv.num_windows()>1) {
    // Create new groups and communicators for each window.
    mpiv.create_local_communicator();
    // Get the local id (in the local communicators).
    mpiv.set_local_id();
  }
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
          << " too few number of the bins were provided.\n" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // REWL.
  rewl<FerroIsing>(&ln_dos, &model, histo_env, &wl_params, window, &mpiv,
      engine);
  // Output.
  merge_ln_dos(&ln_dos, mpiv);
  if (mpiv.myid()%mpiv.multiple() == 0) {
    std::string filename = "./rawdata/lngE_proc" +
        std::to_string(mpiv.myid()/mpiv.multiple()) + ".dat";
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
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  MPI_Finalize();
  return 0;
}


template <typename Model>
void rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 &engine) {
  Model &model(*model_ptr);
  WLParams &wl_params(*wl_params_ptr);
  MPIV &mpiv(*mpiv_ptr);
  MPI_Status status;
  std::uniform_real_distribution<> uniform01_double{0., 1.};
  bool is_flat, is_exchange_accepted;
  double lnf_slowest = wl_params.lnf();
  int exchange_pattern = 0; // 0 or 1.
  int partner;
  double lnf_tmp, val_tmp;
  int swap_count_down = wl_params.swap_every();
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<int> histogram(histo_env.num_bins(), 0);
  // Initialize the configuration.
  while ((val_tmp < window.valmin()) || (window.valmax() < val_tmp)) {
    val_tmp = model.Propose(engine);
    if (std::log(uniform01_double(engine)) <
        ln_dos[histo_env.GetIndex(model.GetVal())] -
        ln_dos[histo_env.GetIndex(val_tmp)]) {
      model.Update();
    }
    ln_dos[histo_env.GetIndex(model.GetVal())] += wl_params.lnf();
  }
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  // Main Wang-Landau routine.
  while (lnf_slowest > wl_params.lnfmin()) {
    for (int i=0; i<wl_params.check_flatness_every(); ++i) {
      for (int j=0; j<wl_params.sweeps(); ++j) {
        val_tmp = model.Propose(engine);
        if ((window.valmin() <= val_tmp) && (val_tmp <= window.valmax()) &&
            (std::log(uniform01_double(engine)) <
            ln_dos[histo_env.GetIndex(model.GetVal())] -
            ln_dos[histo_env.GetIndex(val_tmp)])) {
          // Accept.
          model.Update();
        }
        ln_dos[histo_env.GetIndex(model.GetVal())] += wl_params.lnf();
        histogram[histo_env.GetIndex(model.GetVal())] += 1;
      } // End 1 sweep.
      --swap_count_down;
      // Start RE.
      if ((mpiv.num_windows()>1) && (swap_count_down==0)) {
        swap_count_down = wl_params.swap_every();
        if (mpiv.num_windows()>2) {
          exchange_pattern ^= 1;
        }
        mpiv.set_comm_id(exchange_pattern);
        // Get exchange partner.
        partner = generate_partner(engine, exchange_pattern, mpiv);
        MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
        if (partner != -1) {
          // Replica exchange.
          is_exchange_accepted = replica_exchange<Model>(&val_tmp, partner,
              exchange_pattern, ln_dos, model, histo_env, window, mpiv,
              engine);
          if (is_exchange_accepted) {
            // Exchange configuration.
            model.ExchangeConfig(partner, mpiv.local_comm(mpiv.comm_id()),
                val_tmp);
          }
        }
        // Update histograms (independently of whether RE happened or not).
        ln_dos[histo_env.GetIndex(model.GetVal())] += wl_params.lnf();
        histogram[histo_env.GetIndex(model.GetVal())] += 1;
      } // End RE.
    }
    // Check flatness.
    is_flat = check_histoflat(window, histogram, wl_params.flatness(), mpiv);
    if (is_flat) {
      wl_params.update_lnf();
      for (int &i : histogram) i = 0;
      // Merge g(E) estimators from multiple walkers in the same window.
      merge_ln_dos(&ln_dos, mpiv);
    }
    // Check progress from all other windows.
    lnf_tmp = wl_params.lnf();
    MPI_Allreduce(&lnf_tmp, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX,
        MPI_COMM_WORLD);
    ////
    if (mpiv.myid() == 1) {
      std::cout << "lnf_slowest : " << lnf_slowest << std::endl;
    }
    ////
  } // End while(lnf_slowest>lnfmin) -> this terminates the simulation.
}


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


template <typename Model>
bool replica_exchange(double *val_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const Model &model,
    const HistoEnvManager &histo_env, const WindowManager &window,
    const MPIV &mpiv, std::mt19937 &engine) {
  MPI_Status status;
  std::uniform_real_distribution<> uniform01_double{0., 1.};
  double my_frac, other_frac;
  bool is_exchange_accepted;
  // Get the "value" from my exchange partner.
  *val_partner = model.GetVal();
  MPI_Sendrecv_replace(val_partner, 1, MPI_DOUBLE, partner, 1, partner, 1,
      mpiv.local_comm(mpiv.comm_id()), &status);
  if ((*val_partner < window.valmin()) || (window.valmax() < *val_partner)) {
    my_frac = -1.0;
  } else {
    my_frac = std::exp(ln_dos[histo_env.GetIndex(*val_partner)] -
        ln_dos[histo_env.GetIndex(model.GetVal())]);
  }
  if (mpiv.local_id(exchange_pattern)<mpiv.multiple()) {
    // Receiver calculate combined exchange probability.
    // Get my partner's part of the exchange probability.
    MPI_Recv(&other_frac, 1, MPI_DOUBLE, partner, 2,
        mpiv.local_comm(mpiv.comm_id()), &status);
    // Calculate combined exchange probability and do exchange trial.
    if ((my_frac>0.0)&&(other_frac>0.0)&&
        (uniform01_double(engine)<my_frac*other_frac)) {
      // Exchange accepted.
      is_exchange_accepted = true;
    } else {
      is_exchange_accepted = false;
    }
    MPI_Send(&is_exchange_accepted, 1, MPI_CXX_BOOL, partner, 3,
        mpiv.local_comm(mpiv.comm_id()));
  } else {
    // Send my part of exchange probability and await decision.
    MPI_Send(&my_frac, 1, MPI_DOUBLE, partner, 2,
        mpiv.local_comm(mpiv.comm_id()));
    MPI_Recv(&is_exchange_accepted, 1, MPI_CXX_BOOL, partner, 3,
        mpiv.local_comm(mpiv.comm_id()), &status);
  } // Now all process know whether the replica exchange will be executed.
  return is_exchange_accepted;
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

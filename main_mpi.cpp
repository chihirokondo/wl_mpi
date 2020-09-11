// Ferromagnetic Ising model.
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "include/random.hpp"
#include "include/mpi.hpp"
#include "include/lattice/graph.hpp"
#include "include/ferro_ising.hpp"
#include "include/histo_env_manager.hpp"


int generate_partner(irandom::MTRandom &random, int exchange_pattern,
    const MPIV &mpiv);
bool check_histoflat(size_t imin, size_t imax,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv);
bool replica_exchange(double *val_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const FerroIsing &model,
    const HistoEnvManager &histo_env_manager, const MPIV &mpiv,
    irandom::MTRandom &random, double minval_window, double maxval_window);
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
  HistoEnvManager histo_env_manager(model.ene_min(), model.ene_max(),
    model.num_bins(), true);
  // Original Wang-Landau parameters.
  int check_flatness_every = 500;
  double lnf = 1.0;
  double lnfmin = 1e-8;
  double flatness = 0.95;
  std::vector<double> ln_dos(histo_env_manager.num_bins(), 0.0);
  std::vector<int> histogram(histo_env_manager.num_bins(), 0);
  // Replica exchange Wang-Landau (REWL) parameters.
  bool exchange_accepted;
  double lnf_slowest = lnf;
  double overlap = atof(argv[1]);
  int swap_every = atoi(argv[3]);
  int swap_count_down = swap_every;
  int exchange_pattern = 0; // 0 or 1.
  int partner;
  // Settings for the windows.
  double width = (histo_env_manager.maxval()-histo_env_manager.minval()) /
      (1 + (mpiv.num_windows()-1)*(1-overlap));
  double minval_window = histo_env_manager.minval() +
      (mpiv.myid()/mpiv.multiple())*(1-overlap)*width;
  double maxval_window = minval_window+width;
  size_t imin = histo_env_manager.GetIndex(minval_window);
  size_t imax = histo_env_manager.GetIndex(maxval_window);
  minval_window = histo_env_manager.GetVal(imin, "min");
  maxval_window = histo_env_manager.GetVal(imax, "max");
  if (mpiv.myid() < mpiv.multiple()) {
    imin = 0;
    minval_window = histo_env_manager.minval();
  }
  if (mpiv.myid() >= mpiv.numprocs()-mpiv.multiple()) {
    imax = histo_env_manager.num_bins()-1;
    maxval_window = histo_env_manager.maxval();
  }
  int width_index = imax-imin;
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
  size_t isew = histo_env_manager.GetIndex(minval_window+width*(1-overlap)) -
      imin;
  // For statistics.
  int try_right = 0;
  int try_left = 0;
  int exchange_right = 0;
  int exchange_left = 0;
  // Other variables.
  double val_tmp;
  std::vector<int> config_buf;
  bool is_flat;
  irandom::MTRandom random(atoi(argv[4])+mpiv.myid());
  // Initialize the configuration.
  while ((model.GetVal() < (minval_window + (maxval_window-minval_window)/3)) ||
      (model.GetVal() > (minval_window + 2*(maxval_window-minval_window)/3))) {
    ////
    val_tmp = model.Propose(random);
    model.Update();
    ////
  }
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  // Main Wang-Landau routine.
  while (lnf_slowest > lnfmin) {
    for (int i=0; i<check_flatness_every; ++i) {
      for (int j=0; j<sweeps; ++j) {
        val_tmp = model.Propose(random);
        if ((val_tmp >= minval_window) && (val_tmp <= maxval_window) &&
            (std::log(random.Random()) <
            ln_dos[histo_env_manager.GetIndex(model.GetVal())] -
            ln_dos[histo_env_manager.GetIndex(val_tmp)])) {
          // Accept.
          model.Update();
        }
        ln_dos[histo_env_manager.GetIndex(model.GetVal())] += lnf;
        histogram[histo_env_manager.GetIndex(model.GetVal())] += 1;
      } // End 1 sweep.
      --swap_count_down;
      // Start RE.
      if ((mpiv.num_windows()>1) && (swap_count_down==0)) {
        swap_count_down = swap_every;
        if (mpiv.num_windows()>2) {
          exchange_pattern ^= 1;
        }
        mpiv.set_comm_id(exchange_pattern);
        // Get exchange partner.
        partner = generate_partner(random, exchange_pattern, mpiv);
        MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
        if (partner != -1) {
          // Statistics.
          if (partner > mpiv.local_id(exchange_pattern)) ++try_right;
          else ++try_left;
          // Replica exchange.
          exchange_accepted = replica_exchange(&val_tmp, partner,
              exchange_pattern, ln_dos, model, histo_env_manager, mpiv, random,
              minval_window, maxval_window);
          if (exchange_accepted) {
            // Exchange configuration.
            model.ExchangeConfig(partner, mpiv.local_comm(mpiv.comm_id()),
                val_tmp);
            // Statistics.
            if (partner>mpiv.local_id(exchange_pattern)) ++exchange_right;
            else ++exchange_left;
          }
        }
        // Update histograms (independently of whether RE happened or not).
        ln_dos[histo_env_manager.GetIndex(model.GetVal())] += lnf;
        histogram[histo_env_manager.GetIndex(model.GetVal())] += 1;
      } // End RE.
    }
    // Check flatness.
    is_flat = check_histoflat(imin, imax, histogram, flatness, mpiv);
    if (is_flat) {
      lnf /= 2.0;
      for (int &i : histogram) i = 0;
      // Merge g(E) estimators from multiple walkers in the same window.
      merge_ln_dos(&ln_dos, mpiv);
    }
    // Check progress from all other windows.
    MPI_Allreduce(&lnf, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    ////
    if (mpiv.myid() == 1) {
      std::cout << "lnf_slowest : " << lnf_slowest << std::endl;
    }
    ////
  } // End while(lnf_slowest>lnfmin) -> this terminates the simulation.
  // Output.
  merge_ln_dos(&ln_dos, mpiv);
  if (mpiv.myid()%mpiv.multiple() == 0) {
    std::string filename = "./rawdata/lngE_proc" +
        std::to_string(mpiv.myid()/mpiv.multiple()) + ".dat";
    std::ofstream ofs(filename, std::ios::out);
    ofs << "# dim: " << dim << ", length: " << length << "\n";
    ofs << "# condition_type: " << "sum" << "\n";
    ofs << "# condition_value: " << condition_value << "\n";
    ofs << "# sewing_point: " << isew << "\n";
    ofs << "# number_of_windows: " << mpiv.num_windows() << "\n";
    ofs << "# value \t # lngV\n";
    for (size_t i=imin; i<=imax; ++i) {
      ofs << histo_env_manager.GetVal(i, "mid") << "\t"
          << std::scientific << std::setprecision(15) << ln_dos[i] << "\n";
    }
    ofs << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  MPI_Finalize();
  return 0;
}


int generate_partner(irandom::MTRandom &random, int exchange_pattern,
    const MPIV &mpiv) {
  std::vector<int> partner_list(2*mpiv.multiple());
  int partner;
  // 'head-node' in the window determines pairs of flippariners.
  if (mpiv.local_id(exchange_pattern) == 0) {
    int choose_from = mpiv.multiple();
    int select;
    std::vector<int> lib_re(mpiv.multiple());
    for (int i=0; i<mpiv.multiple(); ++i) lib_re[i] = mpiv.multiple()+i;
    for (int i=0; i<mpiv.multiple(); ++i) {
      select = random.Randrange(choose_from);
      partner_list[i] = lib_re[select];
      partner_list[lib_re[select]] = i;
      --choose_from;
      for (int j=select; j<choose_from; ++j) lib_re[j] = lib_re[j+1];
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


bool check_histoflat(size_t imin, size_t imax,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv) {
  // Check flatness of the histogram for all walkers in the window.
  MPI_Status status;
  bool my_flat = true;
  bool other_flat;
  int num_bins = 0;
  int min_histo = histogram[imin];
  double average = 0.0;
  for (size_t i=imin; i<=imax; ++i) {
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


bool replica_exchange(double *val_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const FerroIsing &model,
    const HistoEnvManager &histo_env_manager, const MPIV &mpiv,
    irandom::MTRandom &random, double minval_window, double maxval_window) {
  MPI_Status status;
  double my_frac, other_frac;
  bool exchange_accepted;
  // Get the "value" from my exchange partner.
  *val_partner = model.GetVal();
  MPI_Sendrecv_replace(val_partner, 1, MPI_DOUBLE, partner, 1, partner, 1,
      mpiv.local_comm(mpiv.comm_id()), &status);
  if ((*val_partner<minval_window) || (maxval_window<*val_partner)) {
    my_frac = -1.0;
  } else {
    my_frac = std::exp(ln_dos[histo_env_manager.GetIndex(*val_partner)] -
        ln_dos[histo_env_manager.GetIndex(model.GetVal())]);
  }
  if (mpiv.local_id(exchange_pattern)<mpiv.multiple()) {
    // Receiver calculate combined exchange probability.
    // Get my partner's part of the exchange probability.
    MPI_Recv(&other_frac, 1, MPI_DOUBLE, partner, 2,
        mpiv.local_comm(mpiv.comm_id()), &status);
    // Calculate combined exchange probability and do exchange trial.
    if ((my_frac>0.0)&&(other_frac>0.0)&&
        (random.Random()<my_frac*other_frac)) {
      // Exchange accepted.
      exchange_accepted = true;
    } else {
      exchange_accepted = false;
    }
    MPI_Send(&exchange_accepted, 1, MPI_CXX_BOOL, partner, 3,
        mpiv.local_comm(mpiv.comm_id()));
  } else {
    // Send my part of exchange probability and await decision.
    MPI_Send(&my_frac, 1, MPI_DOUBLE, partner, 2,
        mpiv.local_comm(mpiv.comm_id()));
    MPI_Recv(&exchange_accepted, 1, MPI_CXX_BOOL, partner, 3,
        mpiv.local_comm(mpiv.comm_id()), &status);
  } // Now all process know whether the replica exchange will be executed.
  return exchange_accepted;
}


void merge_ln_dos(std::vector<double> *ln_dos_ptr, const MPIV &mpiv) {
  MPI_Status status;
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<double> ln_dos_buf = ln_dos;
  if (mpiv.myid()%mpiv.multiple() == 0) {
    // 'root' in window, receive individual g(E) and send merged g(E).
    for (int i=1; i<mpiv.multiple(); ++i) {
      MPI_Recv(&ln_dos_buf[0], (int)ln_dos.size(), MPI_DOUBLE, mpiv.myid()+i,
          77, MPI_COMM_WORLD, &status);
      for (int j=0; j<ln_dos.size(); ++j) ln_dos[j] += ln_dos_buf[j];
    }
    int num_bins = 0;
    double mean = 0.0;
    for (int i=0; i<ln_dos.size(); ++i) {
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
      MPI_Send(&ln_dos[0], (int)ln_dos.size(), MPI_DOUBLE, mpiv.myid()+i, 99,
          MPI_COMM_WORLD);
    }
  } else {
    // Send individual g(E) and receive merged g(E).
    MPI_Send(&ln_dos[0], (int)ln_dos.size(), MPI_DOUBLE,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 77, MPI_COMM_WORLD);
    MPI_Recv(&ln_dos[0], (int)ln_dos.size(), MPI_DOUBLE,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 99, MPI_COMM_WORLD,
        &status);
  }
}

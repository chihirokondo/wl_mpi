// Ferromagnetic Ising model.
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "include/random.hpp"
#include "include/mpi.hpp"
#include "include/lattice/graph.hpp"
#include "include/ferro_ising.hpp"


int generate_partner(irandom::MTRandom &random, int exchange_pattern,
    const MPIV &mpiv);
bool check_histoflat(int imin, int imax, const std::vector<int> &histogram,
    double flatness, const MPIV &mpiv);
bool replica_exchange(int *energy_partner, int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const FerroIsing &model,
    const MPIV &mpiv, irandom::MTRandom &random, double energy_min_window,
    double energy_max_window);
void merge_ln_dos(std::vector<double> *ln_dos_ptr, const MPIV &mpiv);


int main(int argc, char *argv[]) {
  int numprocs, myid, multiple;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  try {
    if (numprocs <= 2) throw 0;
  }
  catch (int err_status) {
    if (myid == 0) {
      std::cerr
          << "ERROR: Few number of processes!\n"
          << "# of processes must be greater than 2."
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
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
          << "1. Overlap between consecutive windows. [double, 0 <= overlap <= 1]\n"
          << "2. Number of walkers per energy subwindow. [integer]\n"
          << "3. Number of Monte Carlo steps between replica exchange. [integer]\n"
          << "4. Random number seed. [integer]\n"
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  // mpiv.multiple() is # of walkers per energy window.
  multiple = atoi(argv[2]);
  try {
    // "numprocs" must be a multiple of "multiple".
    if (numprocs%multiple != 0) throw 0;
  }
  catch (int err_status) {
    if (myid == 0) {
      std::cerr
          << "ERROR: # of processes must be a multiple of the second command line argument.\n"
          << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPIV mpiv(numprocs, myid, multiple);
  // Create new groups and communicators for each energy window.
  mpiv.create_local_communicator();
  // Get the local id (in the local communicators).
  mpiv.set_local_id();
  // Model dependent variables.
  int dim = 2;
  int length = 4;
  lattice::graph lat = lattice::graph::simple(dim, length);
  double condition_value = std::pow(2.0, (double)lat.num_sites());
  FerroIsing model(lat);
  // Original Wang-Landau parameters.
  int check_flatness_every = 500;
  double lnf = 1.0;
  double lnfmin = 1e-8;
  double flatness = 0.95;
  std::vector<double> ln_dos(model.num_values(), 0.0);
  std::vector<int> histogram(model.num_values(), 0);
  // Replica exchange Wang-Landau (REWL) parameters.
  bool is_exchange_accepted;
  double lnf_slowest = lnf;
  double overlap = atof(argv[1]);
  int swap_every = atoi(argv[3]);
  int swap_count_down = swap_every;
  int exchange_pattern = 1; // 0 or 1.
  int partner;
  double energy_width = (model.max_value()-model.min_value())/
      (1.0 + ((double)(mpiv.numprocs()/mpiv.multiple()) - 1.0)*(1.0-overlap));
  double energy_min_window = model.min_value()+
      (mpiv.myid()/mpiv.multiple())*(1.0-overlap)*energy_width;
  double energy_max_window = energy_min_window+energy_width;
  int imin = model.get_index(energy_min_window, "ceil");
  int imax = model.get_index(energy_max_window);
  if (mpiv.myid() >= mpiv.numprocs()-mpiv.multiple()) {
    imax = model.num_values()-1;
  }
  int isew =
      model.get_index(energy_min_window+energy_width*(1-overlap), "ceil");
  isew -= imin;
  // For statistics.
  int try_right = 0;
  int try_left = 0;
  int exchange_right = 0;
  int exchange_left = 0;
  // Other variables.
  int energy_tmp;
  std::vector<int> config_buf;
  bool is_flat;
  irandom::MTRandom random(atoi(argv[4])+mpiv.myid());
  // Initiate configuration.
  while ((model.value() <
      (energy_min_window + (energy_max_window-energy_min_window)/3)) ||
      (model.value() >
      (energy_min_window + 2*(energy_max_window-energy_min_window)/3))) {
    energy_tmp = model.Propose(random);
    model.Update();
  }
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  // Main Wang-Landau routine.
  while (lnf_slowest > lnfmin) {
    for (int i=0; i<check_flatness_every; ++i) {
      for (int j=0; j<lat.num_sites(); ++j) {
        energy_tmp = model.Propose(random);
        if ((energy_tmp >= energy_min_window) &&
            (energy_tmp <= ceil(energy_max_window)) &&
            (std::log(random.Random()) <
            ln_dos[model.get_index(model.value())] -
            ln_dos[model.get_index(energy_tmp)])) {
          // Accept.
          model.Update();
        }
        ln_dos[model.get_index(model.value())] += lnf;
        histogram[model.get_index(model.value())] += 1;
      } // End 1 sweep.
      --swap_count_down;
      // Start RE.
      if (swap_count_down==0) {
        swap_count_down = swap_every;
        exchange_pattern ^= 1;
        mpiv.set_comm_id(exchange_pattern);
        // Get exchange partner.
        partner = generate_partner(random, exchange_pattern, mpiv);
        MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
        if (partner != -1) {
          // Statistics.
          if (partner > mpiv.local_id(exchange_pattern)) ++try_right;
          else ++try_left;
          // Replica exchange.
          is_exchange_accepted = replica_exchange(&energy_tmp, partner,
              exchange_pattern, ln_dos, model, mpiv, random, energy_min_window,
              energy_max_window);
          if (is_exchange_accepted) {
            // Exchange configuration.
            config_buf = model.config();
            MPI_Sendrecv_replace(&config_buf[0], config_buf.size(),
                MPI_INT,partner, 1, partner, 1, mpiv.local_comm(mpiv.comm_id()),
                &status);
            // Update.
            model.Update(energy_tmp, config_buf);
            // Statistics.
            if (partner>mpiv.local_id(exchange_pattern)) ++exchange_right;
            else ++exchange_left;
          }
        }
        ln_dos[model.get_index(model.value())] += lnf;
        histogram[model.get_index(model.value())] += 1;
      } // End RE.
    }
    // Check flatness.
    is_flat = check_histoflat(imin, imax, histogram, flatness, mpiv);
    if (is_flat) {
      lnf /= 2.0;
      for (int &i : histogram) i = 0;
      // Merge g(E) estimators from multiple walkers in the same energy window.
      merge_ln_dos(&ln_dos, mpiv);
    }
    // Check progress from all other windows.
    MPI_Allreduce(&lnf, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    ////
    {
      if (mpiv.myid() == 1) {
        std::cout << "lnf_slowest : " << lnf_slowest << std::endl;
      }
    }
    ////
  } // End while(lnf_slowest>lnfmin) -> this terminates the simulation.
  // Output.
  std::string filename =
      "./rawdata/lngE_proc" + std::to_string(mpiv.myid()) + ".dat";
  std::ofstream ofs(filename, std::ios::out);
  ofs << "# dim: " << dim << ", length: " << length << "\n";
  ofs << "# energy \t # lngE\n";
  ofs << "# condition_type: " << "sum" << "\n";
  ofs << "# condition_value: " << condition_value << "\n";
  ofs << "# sewing_point: " << isew << "\n";
  for (int i=imin; i<=imax; ++i) {
    ofs << model.values(i) << "\t" << ln_dos[i] << "\n";
  }
  ofs << std::endl;
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  MPI_Finalize();
  return 0;
}


int generate_partner(irandom::MTRandom &random, int exchange_pattern,
    const MPIV &mpiv) {
  std::vector<int> partner_list(2*mpiv.multiple());
  int partner;
  // 'head-node' in the energy window determines pairs of flippariners.
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


bool check_histoflat(int imin, int imax, const std::vector<int> &histogram,
    double flatness, const MPIV &mpiv) {
  // Check flatness of the histogram for all walkers in the energy window.
  MPI_Status status;
  bool my_flat = true;
  bool other_flat;
  int num_bins = 0;
  int min_histo = histogram[imin];
  double average = 0.0;
  for (int i=imin; i<=imax; ++i) {
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
  // Now talk to all the other walkers in the energy window.
  // (! this whole thing can be reduced to an MPI_Allreduce once there are separate communicators for energy windows !)
  if (mpiv.myid()%mpiv.multiple() == 0) {
    // 'root' in energy window, receive individual flatnesses.
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
    MPI_Send(&my_flat, 1, MPI_CXX_BOOL, mpiv.myid()-(mpiv.myid()%mpiv.multiple()),
        66, MPI_COMM_WORLD);
    MPI_Recv(&other_flat, 1, MPI_CXX_BOOL,
        mpiv.myid()-(mpiv.myid()%mpiv.multiple()), 88, MPI_COMM_WORLD, &status);
    my_flat = other_flat;  // Replace individual flatness by merged.
  }
  return my_flat;
  // Note: By now, myflat refers to the 'collective' flatness in the energy window, not the flatness of an individual walker.
}


bool replica_exchange(int *energy_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const FerroIsing &model,
    const MPIV &mpiv, irandom::MTRandom &random, double energy_min_window,
    double energy_max_window) {
  MPI_Status status;
  double my_frac, other_frac;
  bool is_exchange_accepted;
  // Get the energy from my exchange partner.
  *energy_partner = model.value();
  MPI_Sendrecv_replace(energy_partner, 1, MPI_INT, partner, 1, partner, 1,
      mpiv.local_comm(mpiv.comm_id()), &status);
  if ((*energy_partner>energy_max_window) ||
      (*energy_partner<energy_min_window)) {
    my_frac = -1.0;
  } else {
    my_frac = std::exp(ln_dos[model.get_index(*energy_partner)] -
        ln_dos[model.get_index(model.value())]);
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
    // 'root' in energy window, receive individual g(E) and send merged g(E).
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
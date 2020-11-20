#ifndef WANGLANDAU_REWL_H_
#define WANGLANDAU_REWL_H_


#include <mpi.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include "mpi_setting.hpp"
#include "stop_callback.hpp"
#include "histo_env_manager.hpp"
#include "window.hpp"
#include "wl_params.hpp"
#include "log_for_json.hpp"


template <typename Model>
int rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 *engine,
    StopCallback stop_callback, bool from_the_top);
int generate_partner(std::mt19937 *engine, int exchange_pattern,
    const MPIV &mpiv);
template <typename Model>
bool replica_exchange(double *val_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const Model &model,
    const HistoEnvManager &histo_env, const WindowManager &window,
    const MPIV &mpiv, std::mt19937 *engine);
bool check_histoflat(const WindowManager &window,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv);
void merge_ln_dos(std::vector<double> *ln_dos_ptr, const MPIV &mpiv);


template <typename Model>
int rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 *engine,
    StopCallback stop_callback, bool from_the_top) {
  Model &model(*model_ptr);
  WLParams &wl_params(*wl_params_ptr);
  MPIV &mpiv(*mpiv_ptr);
  std::mt19937 &engine(*engine);
  MPI_Status status;
  std::string log_file_name = "./log/proc" + std::to_string(mpiv.myid()) +
      ".json";
  std::ofstream ofs_log(log_file_name, std::ios::out);
  std::string model_file_name = "./log/proc" + std::to_string(mpiv.myid()) +
      "_model_state";
  std::ofstream ofs_model_log(model_file_name, std::ios::out);
  int running_state = 0;
  std::uniform_real_distribution<> uniform01_double{0., 1.};
  bool is_flat, is_exchange_accepted;
  double lnf_slowest = wl_params.lnf();
  int exchange_pattern = 0; // 0 or 1.
  int partner;
  double lnf_tmp, val_proposed;
  int swap_count_down = wl_params.swap_every();
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<int> histogram(histo_env.num_bins(), 0);
  // Initialize the configuration.
  val_proposed = model.GetVal();
  while ((val_proposed < window.valmin()) || (window.valmax() < val_proposed)) {
    val_proposed = model.Propose(&engine);
    if (std::log(uniform01_double(engine)) <
        ln_dos[histo_env.GetIndex(model.GetVal())] -
        ln_dos[histo_env.GetIndex(val_proposed)]) {
      model.Update();
    }
    ln_dos[histo_env.GetIndex(model.GetVal())] += wl_params.lnf();
  }
  ++running_state;
  // Main Wang-Landau routine.
  while (lnf_slowest > wl_params.lnfmin()) {
    // Check elapsed time.
    if (stop_callback()) {
      // Leave log files and stop.
      write_log_json(&ofs_log, running_state, mpiv, wl_params, ln_dos, engine,
          histogram, swap_count_down, exchange_pattern, lnf_slowest);
      model.WriteState(&ofs_model_log);
      return running_state;
    }
    for (int i=0; i<wl_params.check_flatness_every(); ++i) {
      for (int j=0; j<wl_params.sweeps(); ++j) {
        val_proposed = model.Propose(&engine);
        if ((window.valmin() <= val_proposed) &&
            (val_proposed <= window.valmax()) &&
            (std::log(uniform01_double(engine)) <
            ln_dos[histo_env.GetIndex(model.GetVal())] -
            ln_dos[histo_env.GetIndex(val_proposed)])) {
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
        partner = generate_partner(&engine, exchange_pattern, mpiv);
        MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
        if (partner != -1) {
          // Replica exchange.
          is_exchange_accepted = replica_exchange<Model>(&val_proposed, partner,
              exchange_pattern, ln_dos, model, histo_env, window, mpiv,
              &engine);
          if (is_exchange_accepted) {
            // Exchange configuration.
            model.ExchangeConfig(partner, mpiv.local_comm(mpiv.comm_id()),
                val_proposed);
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
  ++running_state;
  // Leave log files.
  write_log_json(&ofs_log, running_state, mpiv, wl_params, ln_dos, engine,
      histogram, swap_count_down, exchange_pattern, lnf_slowest);
  model.WriteState(&ofs_model_log);
  return running_state;
}


int generate_partner(std::mt19937 *engine, int exchange_pattern,
    const MPIV &mpiv) {
  std::mt19937 &engine(*engine);
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


template <typename Model>
bool replica_exchange(double *val_partner ,int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const Model &model,
    const HistoEnvManager &histo_env, const WindowManager &window,
    const MPIV &mpiv, std::mt19937 *engine) {
  std::mt19937 &engine(*engine);
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


#endif // WANGLANDAU_REWL_H_
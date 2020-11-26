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
#include "is_time_out.hpp"
#include "histo_env_manager.hpp"
#include "window.hpp"
#include "wl_params.hpp"
#include "log_for_json.hpp"


template <typename Model>
int rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 &engine,
    double timelimit_secs, bool from_the_top);
int generate_partner(std::mt19937 &engine, int exchange_pattern,
    const MPIV &mpiv);
template <typename Model>
void exchange_config(Model *model_ptr, int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const HistoEnvManager &histo_env,
    const WindowManager &window, const MPIV &mpiv, std::mt19937 &engine);
bool are_all_hists_flat(const WindowManager &window,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv);
void merge_ln_dos(std::vector<double> *ln_dos_ptr, const MPIV &mpiv);


template <typename Model>
int rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 &engine,
    double timelimit_secs, bool from_the_top) {
  std::vector<double> &ln_dos(*ln_dos_ptr);
  Model &model(*model_ptr);
  WLParams &wl_params(*wl_params_ptr);
  MPIV &mpiv(*mpiv_ptr);
  MPI_Status status;
  IsTimeOut is_time_out(timelimit_secs);
  // For log files.
  std::string log_file_name = "./log/proc" + std::to_string(mpiv.myid()) +
      ".json";
  std::string model_file_name = "./log/proc" + std::to_string(mpiv.myid()) +
      "_model_state";
  std::ofstream ofs_log(log_file_name, std::ios::out);
  std::ofstream ofs_model_log(model_file_name, std::ios::out);
  int running_state = 0;
  std::vector<int> histogram(histo_env.num_bins(), 0);
  int swap_count_down = wl_params.swap_every();
  int exchange_pattern = 0; // 0 or 1.
  double lnf_slowest = wl_params.lnf();
  if (from_the_top) {
    // Initialize the configuration.
    std::uniform_real_distribution<> uniform01{0., 1.};
    double val_new = model.val();
    while ((val_new < window.valmin()) || (window.valmax() < val_new)) {
      val_new = model.Propose(engine);
      size_t old_i = histo_env.GetIndex(model.val());
      size_t new_i = histo_env.GetIndex(val_new);
      if (std::log(uniform01(engine)) < ln_dos[old_i] - ln_dos[new_i]) {
        model.Update();
      }
      ln_dos[histo_env.GetIndex(model.val())] += wl_params.lnf();
    }
  } else {
    // Read log files.
    std::ifstream ifs_log(log_file_name, std::ios::in);
    bool is_consistent;
    if (mpiv.myid() == 0) {
      is_consistent = check_log_json(&ifs_log, mpiv, wl_params, ln_dos);
    }
    MPI_Bcast(&is_consistent, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    if (!is_consistent) return -1; // Error occured.
    set_from_log_json(&ifs_log, &wl_params, &ln_dos, engine, &histogram,
        &swap_count_down, &exchange_pattern, &lnf_slowest);
    // Read model log file.
    std::ifstream ifs_model_log(model_file_name, std::ios::in);
    model.SetFromLog(&ifs_model_log);
  }
  // Main Wang-Landau routine.
  while (lnf_slowest > wl_params.lnfmin()) {
    // Check elapsed time.
    bool should_stop = is_time_out();
    MPI_Bcast(&should_stop, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    if (should_stop) {
      // Leave log files.
      write_log_json(&ofs_log, running_state, mpiv, wl_params, ln_dos, engine,
          histogram, swap_count_down, exchange_pattern, lnf_slowest);
      model.WriteState(&ofs_model_log);
      return running_state;
    }
    for (int i=0; i<wl_params.check_flatness_every(); ++i) {
      for (int j=0; j<wl_params.sweeps(); ++j) {
        std::uniform_real_distribution<> uniform01{0., 1.};
        double val_new = model.Propose(engine);
        size_t old_i = histo_env.GetIndex(model.val());
        size_t new_i = histo_env.GetIndex(val_new);
        if ((window.valmin() <= val_new) && (val_new <= window.valmax()) &&
            (std::log(uniform01(engine)) < ln_dos[old_i] - ln_dos[new_i])) {
          // Accept.
          model.Update();
        }
        ln_dos[histo_env.GetIndex(model.val())] += wl_params.lnf();
        histogram[histo_env.GetIndex(model.val())] += 1;
      } // End 1 sweep.
      --swap_count_down;
      // Start RE.
      if ((mpiv.num_windows()>1) && (swap_count_down==0)) {
        swap_count_down = wl_params.swap_every();
        if (mpiv.num_windows()>2) exchange_pattern ^= 1;
        mpiv.set_comm_id(exchange_pattern);
        // Get exchange partner.
        int partner = generate_partner(engine, exchange_pattern, mpiv);
        MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
        if (partner != -1) {
          // Replica exchange.
          exchange_config<Model>(&model, partner, exchange_pattern, ln_dos,
              histo_env, window, mpiv, engine);
        }
        // Update histograms (independently of whether RE happened or not).
        ln_dos[histo_env.GetIndex(model.val())] += wl_params.lnf();
        histogram[histo_env.GetIndex(model.val())] += 1;
      } // End RE.
    }
    // Check flatness of the histograms in the window.
    if (are_all_hists_flat(window, histogram, wl_params.flatness(), mpiv)) {
      wl_params.update_lnf();
      for (int &i : histogram) i = 0;
      // Merge g(E) estimators from multiple walkers in the same window.
      merge_ln_dos(&ln_dos, mpiv);
    }
    // Check progress from all other windows.
    double lnf_tmp = wl_params.lnf();
    MPI_Allreduce(&lnf_tmp, &lnf_slowest, 1, MPI_DOUBLE, MPI_MAX,
        MPI_COMM_WORLD);
    ////
    if (mpiv.myid() == 0) {
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


int generate_partner(std::mt19937 &engine, int exchange_pattern,
    const MPIV &mpiv) {
  std::vector<size_t> partner_list(2*mpiv.multiple());
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
  if (mpiv.comm_id() == -1) return -1;
  int partner;
  MPI_Scatter(&partner_list[0], 1, MPI_INT, &partner, 1, MPI_INT, 0,
      mpiv.local_comm(mpiv.comm_id()));
  return partner;
}


template <typename Model>
void exchange_config(Model *model_ptr, int partner, int exchange_pattern,
    const std::vector<double> &ln_dos, const HistoEnvManager &histo_env,
    const WindowManager &window, const MPIV &mpiv, std::mt19937 &engine) {
  MPI_Status status;
  Model &model(*model_ptr);
  // Get the "value" from my exchange partner.
  double val_partner = model.val();
  MPI_Sendrecv_replace(&val_partner, 1, MPI_DOUBLE, partner, 1, partner, 1,
      mpiv.local_comm(mpiv.comm_id()), &status);
  double my_frac;
  if ((val_partner < window.valmin()) || (window.valmax() < val_partner)) {
    my_frac = -1.0;
  } else {
    size_t partner_i = histo_env.GetIndex(val_partner);
    size_t my_i = histo_env.GetIndex(model.val());
    my_frac = std::exp(ln_dos[partner_i] - ln_dos[my_i]);
  }
  bool is_exchange_accepted;
  if (mpiv.local_id(exchange_pattern) < mpiv.multiple()) {
    // Receiver calculate combined exchange probability.
    // Get my partner's part of the exchange probability.
    double other_frac;
    MPI_Recv(&other_frac, 1, MPI_DOUBLE, partner, 2,
        mpiv.local_comm(mpiv.comm_id()), &status);
    // Calculate combined exchange probability and do exchange trial.
    std::uniform_real_distribution<> uniform01{0., 1.};
    if ((my_frac>0.0) && (other_frac>0.0) &&
        (uniform01(engine)<my_frac*other_frac)) {
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
  if (is_exchange_accepted) {
    // Exchange configuration.
    model.set_val(val_partner);
    model.ExchangeConfig(partner, mpiv.local_comm(mpiv.comm_id()));
  }
}


// Check flatness of the histogram for all walkers in the window.
bool are_all_hists_flat(const WindowManager &window,
    const std::vector<int> &histogram, double flatness, const MPIV &mpiv) {
  MPI_Status status;
  int num_bins = 0;
  int min_histo = 0;
  double average = 0.0;
  for (size_t i=window.imin(); i<=window.imax(); ++i) {
    if (histogram[i]!=0) {
      ++num_bins;
      average += (double)histogram[i];
      if ((min_histo==0) || (histogram[i]<min_histo)) min_histo = histogram[i];
    }
  }
  average /= num_bins;
  bool my_flat = true;
  if ((double)min_histo < flatness*average) my_flat = false;
  // Now talk to all the other walkers in the window.
  bool other_flat;
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
        (mpiv.myid()/mpiv.multiple())*mpiv.multiple(), 66, MPI_COMM_WORLD);
    MPI_Recv(&other_flat, 1, MPI_CXX_BOOL,
        (mpiv.myid()/mpiv.multiple())*mpiv.multiple(), 88, MPI_COMM_WORLD,
        &status);
    my_flat = other_flat;  // Replace individual flatness by merged one.
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
      if (ln_dos[i] != 0.0) ++num_bins;
    }
    mean /= num_bins;
    for (double &ln_dos_i : ln_dos) {
      if (ln_dos_i != 0.0) ln_dos_i -= mean;
    }
    for (int i=1; i<mpiv.multiple(); ++i) {
      MPI_Send(&ln_dos[0], ln_dos.size(), MPI_DOUBLE, mpiv.myid()+i, 99,
          MPI_COMM_WORLD);
    }
  } else {
    // Send individual g(E) and receive merged g(E).
    MPI_Send(&ln_dos[0], ln_dos.size(), MPI_DOUBLE,
        (mpiv.myid()/mpiv.multiple())*mpiv.multiple(), 77, MPI_COMM_WORLD);
    MPI_Recv(&ln_dos[0], ln_dos.size(), MPI_DOUBLE,
        (mpiv.myid()/mpiv.multiple())*mpiv.multiple(), 99, MPI_COMM_WORLD,
        &status);
  }
}


#endif // WANGLANDAU_REWL_H_
#ifndef WANGLANDAU_REWL_H_
#define WANGLANDAU_REWL_H_


#include <fstream>
#include <random>
#include <string>
#include "mpi_setting.hpp"
#include "stop_callback.hpp"


template <typename Model>
void rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 &engine,
    StopCallback stop_callback);
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


template <typename Model>
void rewl(std::vector<double> *ln_dos_ptr, Model *model_ptr,
    const HistoEnvManager &histo_env, WLParams *wl_params_ptr,
    const WindowManager &window, MPIV *mpiv_ptr, std::mt19937 &engine,
    StopCallback stop_callback) {
  Model &model(*model_ptr);
  WLParams &wl_params(*wl_params_ptr);
  MPIV &mpiv(*mpiv_ptr);
  MPI_Status status;
  std::uniform_real_distribution<> uniform01_double{0., 1.};
  bool is_flat, is_exchange_accepted;
  const double LNF_SLOWEST = wl_params.lnf();
  int exchange_pattern = 0; // 0 or 1.
  int partner;
  double lnf_tmp, val_proposed;
  int swap_count_down = wl_params.swap_every();
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<int> histogram(histo_env.num_bins(), 0);
  // Initialize the configuration.
  val_proposed = model.GetVal();
  while ((val_proposed < window.valmin()) || (window.valmax() < val_proposed)) {
    val_proposed = model.Propose(engine);
    if (std::log(uniform01_double(engine)) <
        ln_dos[histo_env.GetIndex(model.GetVal())] -
        ln_dos[histo_env.GetIndex(val_proposed)]) {
      model.Update();
    }
    ln_dos[histo_env.GetIndex(model.GetVal())] += wl_params.lnf();
  }
  // Main Wang-Landau routine.
  while (LNF_SLOWEST > wl_params.lnfmin()) {
    for (int i=0; i<wl_params.check_flatness_every(); ++i) {
      for (int j=0; j<wl_params.sweeps(); ++j) {
        val_proposed = model.Propose(engine);
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
        partner = generate_partner(engine, exchange_pattern, mpiv);
        MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
        if (partner != -1) {
          // Replica exchange.
          is_exchange_accepted = replica_exchange<Model>(&val_proposed, partner,
              exchange_pattern, ln_dos, model, histo_env, window, mpiv, engine);
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
    MPI_Allreduce(&lnf_tmp, &LNF_SLOWEST, 1, MPI_DOUBLE, MPI_MAX,
        MPI_COMM_WORLD);
    ////
    if (mpiv.myid() == 1) {
      std::cout << "LNF_SLOWEST : " << LNF_SLOWEST << std::endl;
    }
    ////
    // Time check.
    if (stop_callback()) {
      //// Leave log files and stop.
      std::string filename = "./log/proc" + mpiv.myid() + ".json";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } // End while(LNF_SLOWEST>lnfmin) -> this terminates the simulation.
  //// leave log files.
}


#endif // WANGLANDAU_REWL_H_
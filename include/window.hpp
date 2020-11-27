#ifndef WANGLANDAU_WINDOW_H_
#define WANGLANDAU_WINDOW_H_


#include <iostream>
#include "mpi_setting.hpp"
#include "histo_env_manager.hpp"


class WindowManager {
 public:
  WindowManager(HistoEnvManager histo_env, MPIV mpiv, double overlap);
  // Gettor.
  double valmin() const {return valmin_;}
  double valmax() const {return valmax_;}
  size_t imin() const {return imin_;}
  size_t imax() const {return imax_;}
  size_t isew() const {return isew_;}
  size_t iwidth() const {return iwidth_;}
 private:
  double valmin_, valmax_;
  size_t imin_, imax_, isew_, iwidth_;
};


WindowManager::WindowManager(HistoEnvManager histo_env, MPIV mpiv,
    double overlap) {
  double width = (histo_env.maxval()-histo_env.minval()) /
      (1 + (mpiv.num_windows()-1)*(1-overlap));
  valmin_ = histo_env.minval() +
      (mpiv.myid()/mpiv.num_walkers_window())*(1-overlap)*width;
  valmax_ = valmin_+width;
  imin_ = histo_env.GetIndex(valmin_);
  imax_ = histo_env.GetIndex(valmax_);
  valmin_ = histo_env.GetVal(imin_, "min");
  valmax_ = histo_env.GetVal(imax_, "max");
  if (mpiv.myid() < mpiv.num_walkers_window()) {
    imin_ = 0;
    valmin_ = histo_env.minval();
  }
  if (mpiv.myid() >= mpiv.numprocs()-mpiv.num_walkers_window()) {
    imax_ = histo_env.num_bins()-1;
    valmax_ = histo_env.maxval();
  }
  iwidth_ = imax_ - imin_;
  isew_ = histo_env.GetIndex(valmin_+width*(1-overlap)) - imin_;
}


#endif // WANGLANDAU_WINDOW_H_

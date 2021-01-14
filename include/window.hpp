#ifndef WANGLANDAU_WINDOW_H_
#define WANGLANDAU_WINDOW_H_


#include <stdexcept>
#include <iostream>
#include <string>
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
 private:
  double valmin_, valmax_;
  size_t imin_, imax_;
};


inline WindowManager::WindowManager(HistoEnvManager histo_env, MPIV mpiv,
    double overlap) {
  if (overlap < 0 || overlap > 1) {
    std::string message =
        "ERROR: Provided argument \"overlap\" was out of range.";
    message += "       Must be 0 <= \"overlap\" <=1.";
    throw std::out_of_range(message);
  }
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
  if (imax_-imin_ < 1) {
    std::string message = "ERROR: Width of the window is narrow ";
    message += "compared to bin width\n";
    throw std::out_of_range(message);
  }
}


#endif // WANGLANDAU_WINDOW_H_

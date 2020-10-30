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


#endif // WANGLANDAU_WINDOW_H_

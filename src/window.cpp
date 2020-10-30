#include "window.hpp"


WindowManager::WindowManager(HistoEnvManager histo_env, MPIV mpiv,
    double overlap) {
  double width = (histo_env.maxval()-histo_env.minval()) /
    (1 + (mpiv.num_windows()-1)*(1-overlap));
  valmin_ = histo_env.minval() + (mpiv.myid()/mpiv.multiple())*(1-overlap)*width;
  valmax_ = valmin_+width;
  imin_ = histo_env.GetIndex(valmin_);
  imax_ = histo_env.GetIndex(valmax_);
  valmin_ = histo_env.GetVal(imin_, "min");
  valmax_ = histo_env.GetVal(imax_, "max");
  if (mpiv.myid() < mpiv.multiple()) {
    imin_ = 0;
    valmin_ = histo_env.minval();
  }
  if (mpiv.myid() >= mpiv.numprocs()-mpiv.multiple()) {
    imax_ = histo_env.num_bins()-1;
    valmax_ = histo_env.maxval();
  }
  iwidth_ = imax_ - imin_;
  isew_ = histo_env.GetIndex(valmin_+width*(1-overlap)) - imin_;
}

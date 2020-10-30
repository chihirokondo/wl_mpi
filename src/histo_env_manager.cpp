#include "histo_env_manager.hpp"


HistoEnvManager::HistoEnvManager(double min, double max, size_t num_bins,
    bool centering) : min_(min), max_(max), num_bins_(num_bins),
    centering_(centering) {
  if (centering_) {
    width_ = (max_-min_)/(num_bins_-1);
  } else {
    width_ = (max_-min_)/num_bins_;
  }
}


size_t HistoEnvManager::GetIndex(double val) const {
  if (centering_) {
    return (size_t)((val-min_+0.5*width_)/width_);
  } else {
    return (size_t)((val-min_)/width_);
  }
}


double HistoEnvManager::GetVal(size_t index, std::string loc) const {
  if (loc == "min") {
    if (centering_) {
      return min_ + (index-0.5)*width_;
    } else {
      return min_ + index*width_;
    }
  } else if (loc == "mid") {
    if (centering_) {
      return min_ + index*width_;
    } else {
      return min_ + (index+0.5)*width_;
    }
  } else if (loc == "max") {
    if (centering_) {
      return min_ + (index+0.5)*width_;
    } else {
      return min_ + (index+1)*width_;
    }
  } else {
    // Same as the case of "mid".
    if (centering_) {
      return min_ + index*width_;
    } else {
      return min_ + (index+0.5)*width_;
    }
  }
}
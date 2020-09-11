#ifndef WANGLANDAU_HISTO_ENV_MANAGER_H_
#define WANGLANDAU_HISTO_ENV_MANAGER_H_


#include <iostream>
#include <string>


class HistoEnvManager {
 public:
  HistoEnvManager(double min, double max, size_t num_bins,
      bool centering);
  static HistoEnvManager Construct(double min, double max, size_t num_bins,
      bool centering) {
    return HistoEnvManager(min, max, num_bins, centering);
  }
  size_t GetIndex(double val) const;
  double GetVal(size_t index, std::string loc) const;
  // Gettor.
  double minval() const {return min_;}
  double maxval() const {return max_;}
  size_t num_bins() const {return num_bins_;}
 private:
  double min_, max_, width_;
  size_t num_bins_;
  bool centering_;
};


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
  }
}


#endif // WANGLANDAU_HISTO_ENV_MANAGER_H_
#ifndef WANGLANDAU_HISTO_ENV_MANAGER_H_
#define WANGLANDAU_HISTO_ENV_MANAGER_H_


#include <iostream>
#include <string>


class HistoEnvManager {
 public:
  HistoEnvManager(double min, double max, size_t num_bins,
      bool centering) :
    min_(min),
    max_(max),
    num_bins_(num_bins),
    centering_(centering) {
    if (centering_) {
      width_ = (max_-min_)/(num_bins_-1);
    } else {
      width_ = (max_-min_)/num_bins_;
    }
  }
  static HistoEnvManager Construct(double min, double max, size_t num_bins,
      bool centering) {
    return HistoEnvManager(min, max, num_bins, centering);
  }
  size_t GetIndex(double val) const {
    if (centering_) return (size_t)((val-min_+0.5*width_)/width_);
    return (size_t)((val-min_)/width_);
  }
  double GetVal(size_t index, std::string loc) const {
    if (loc == "min") {
      if (centering_) return min_ + (index-0.5)*width_;
      return min_ + index*width_;
    }
    if (loc == "max") {
      if (centering_) return min_ + (index+0.5)*width_;
      return min_ + (index+1)*width_;
    } 
    // Return middle value.
    if (centering_) return min_ + index*width_;
    return min_ + (index+0.5)*width_;
  }
  // Gettor.
  double minval() const {return min_;}
  double maxval() const {return max_;}
  size_t num_bins() const {return num_bins_;}
 private:
  double min_, max_, width_;
  size_t num_bins_;
  bool centering_;
};


#endif // WANGLANDAU_HISTO_ENV_MANAGER_H_
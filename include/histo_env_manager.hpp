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


#endif // WANGLANDAU_HISTO_ENV_MANAGER_H_
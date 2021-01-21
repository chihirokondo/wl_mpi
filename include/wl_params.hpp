#ifndef WANGLANDAU_WL_PARAMS_H_
#define WANGLANDAU_WL_PARAMS_H_


#include <stdexcept>
#include <iostream>
#include <string>


class WLParams {
 public:
  WLParams(int check_flatness_every, double lnf, double lnfmin, double flatness,
      double overlap, int exch_every) :
    check_flatness_every_(check_flatness_every),
    lnf_(lnf),
    lnfmin_(lnfmin),
    flatness_(flatness),
    overlap_(overlap),
    exch_every_(exch_every) {
    if ((flatness_<0) || (flatness_>1)) {
      std::string message = "ERROR: Out of range argument is provided.\n";
      message += "       The 3rd argument is defined on [0, 1]";
      throw std::out_of_range(message);
    }
    if ((overlap_<0) || (overlap_>1)) {
      std::string message = "ERROR: Out of range argument is provided.\n";
      message += "       The 4th argument is defined on [0, 1]";
      throw std::out_of_range(message);
    }
  }
  void update_lnf() {lnf_ /= 2.0;}
  void set_lnf(double lnf) {lnf_ = lnf;}
  // Gettor.
  int check_flatness_every() const {return check_flatness_every_;}
  double lnf() const {return lnf_;}
  double lnfmin() const {return lnfmin_;}
  double flatness() const {return flatness_;}
  double overlap() const {return overlap_;}
  int exch_every() const {return exch_every_;}
 private:
  int check_flatness_every_, exch_every_;
  double lnf_, lnfmin_, flatness_, overlap_;
};


#endif // WANGLANDAU_WL_PARAMS_H_
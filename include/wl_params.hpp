#ifndef WANGLANDAU_WL_PARAMS_H_
#define WANGLANDAU_WL_PARAMS_H_


#include <iostream>


class WLParams {
 public:
  WLParams(int sweeps, int check_flatness_every, double lnf, double lnfmin,
      double flatness, int swap_every) :
    sweeps_(sweeps),
    check_flatness_every_(check_flatness_every),
    lnf_(lnf),
    lnfmin_(lnfmin),
    flatness_(flatness),
    swap_every_(swap_every) {}
  void update_lnf() {lnf_ /= 2.0;}
  void set_lnf(double lnf) {lnf_ = lnf;}
  // Gettor.
  int sweeps() const {return sweeps_;}
  int check_flatness_every() const {return check_flatness_every_;}
  double lnf() const {return lnf_;}
  double lnfmin() const {return lnfmin_;}
  double flatness() const {return flatness_;}
  int swap_every() const {return swap_every_;}
 private:
  int check_flatness_every_, swap_every_, sweeps_;
  double lnf_, lnfmin_, flatness_;
};


#endif // WANGLANDAU_WL_PARAMS_H_
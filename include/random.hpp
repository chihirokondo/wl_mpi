#ifndef WANGLANDAU_RANDOM_H_
#define WANGLANDAU_RANDOM_H_

#include <cmath>
#include <random>

namespace irandom {
class MTRandom {
 public:
  MTRandom(unsigned int seed=1) : seed_(seed) {
    std::mt19937 mt_gen_tmp{seed};
    mt_gen_ = mt_gen_tmp;
  }
  double Random() {
    return Uniform01_double_(mt_gen_);
  }
  int Randrange(int arg_max) {
    return (int)(Uniform01_double_(mt_gen_)*arg_max);
  }
  int Randrange(int arg_min, int arg_max) {
    return (int)(Uniform01_double_(mt_gen_)*(arg_max-arg_min)) + arg_min;
  }
  double ExponentialDistribution(double lambda) {
    return -std::log(1.0-Uniform01_double_(mt_gen_))/(lambda);
  }

 private:
  unsigned int seed_;
  std::mt19937 mt_gen_{seed_}; // Check!!!
  std::uniform_real_distribution<> Uniform01_double_{0., 1.};
};

} // end namespace Random

#endif // WANGLANDAU_RANDOM_H_

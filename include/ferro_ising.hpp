#ifndef WANGLANDAU_FERRO_ISING_H_
#define WANGLANDAU_FERRO_ISING_H_


#include <iostream>
#include <string>
#include "lattice/graph.hpp"
#include "random"


class FerroIsing {
 public:
  FerroIsing(lattice::graph lat);
  int get_index(double ene, std::string mode="floor") const {
    if (mode == "floor") {
      return (int)(ene+lat_.num_bonds())/4;
    } else if (mode=="ceil") {
      return (int)ceil((ene+lat_.num_bonds())/4);
    } else {
      return -1;
    }
  }
  int Propose(irandom::MTRandom &random) {
    energy_proposed_ = energy_;
    site_ = random.Randrange(lat_.num_sites());
    for (int i=0; i<lat_.num_neighbors(site_); ++i) {
      energy_proposed_ +=
          2*spin_config_[site_]*spin_config_[lat_.neighbor(site_, i)];
    }
    return energy_proposed_;
  }
  void Update() {
    energy_ = energy_proposed_;
    spin_config_[site_] -= 2*spin_config_[site_];
  }
  // Gettor.
  int values(int index) {return energies_[index];}
  size_t num_values() {return energies_.size();}
  int min_value() {return energies_.front();}
  int max_value() {return energies_.back();}
  std::vector<int> config() {return spin_config_;}
  int config(int index) {return spin_config_[index];}
  int value() {return energy_;}
 private:
  const lattice::graph lat_;
  int energy_proposed_, site_, energy_;
  std::vector<int> energies_;
  std::vector<int> spin_config_;
};


FerroIsing::FerroIsing(lattice::graph lat) : lat_(lat) {
  spin_config_.resize(lat_.num_sites(), 1);
  spin_config_.shrink_to_fit();
  energy_ = -2*lat_.num_sites();
  for (int i=0; i<=lat_.num_sites(); ++i) {
    energies_.push_back(energy_+4*i);
  }
}


#endif // WANGLANDAU_FERRO_ISING_H_
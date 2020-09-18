#ifndef WANGLANDAU_FERRO_ISING_H_
#define WANGLANDAU_FERRO_ISING_H_


#include <mpi.h>
#include <iostream>
#include <string>
#include "lattice/graph.hpp"
#include "random.hpp"


class FerroIsing {
 public:
  FerroIsing(lattice::graph lat);
  double Propose(irandom::MTRandom &random) {
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
  void ExchangeConfig(int partner, MPI_Comm local_comm, double energy_new);
  // Gettor.
  double GetVal() const {return energy_;}
  size_t num_bins() const {return num_bins_;}
  double ene_min() const {return ene_min_;}
  double ene_max() const {return ene_max_;}
 private:
  const lattice::graph lat_;
  int site_;
  size_t num_bins_;
  double ene_min_, ene_max_, energy_proposed_, energy_;
  std::vector<int> spin_config_;
};


FerroIsing::FerroIsing(lattice::graph lat) : lat_(lat) {
  ene_min_ = -2*(double)lat_.num_sites();
  ene_max_ = -ene_min_;
  num_bins_ = lat_.num_sites()+1;
  spin_config_.resize(lat_.num_sites(), 1);
  spin_config_.shrink_to_fit();
  energy_ = ene_min_;
}


void FerroIsing::ExchangeConfig(int partner, MPI_Comm local_comm,
    double energy_new) {
  MPI_Status status;
  MPI_Sendrecv_replace(&spin_config_[0], spin_config_.size(), MPI_INT,
      partner, 1, partner, 1, local_comm, &status);
  energy_ = energy_new;
}


#endif // WANGLANDAU_FERRO_ISING_H_
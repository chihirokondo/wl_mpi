#include "ferro_ising.hpp"


FerroIsing::FerroIsing(lattice::graph lat) : lat_(lat) {
  ene_min_ = -2*(double)lat_.num_sites();
  ene_max_ = -ene_min_;
  num_bins_ = lat_.num_sites()+1;
  spin_config_.resize(lat_.num_sites(), 1);
  spin_config_.shrink_to_fit();
  energy_ = ene_min_;
  std::uniform_int_distribution<size_t> dist(0, lat_.num_sites()-1);
  dist_ = dist;
}


void FerroIsing::ExchangeConfig(int partner, MPI_Comm local_comm,
    double energy_new) {
  MPI_Status status;
  MPI_Sendrecv_replace(&spin_config_[0], spin_config_.size(), MPI_INT,
      partner, 1, partner, 1, local_comm, &status);
  energy_ = energy_new;
}
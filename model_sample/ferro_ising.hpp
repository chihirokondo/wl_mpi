#ifndef WANGLANDAU_FERRO_ISING_H_
#define WANGLANDAU_FERRO_ISING_H_


#include <random>
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>
#include "lattice/graph.hpp"
using json = nlohmann::json;


class FerroIsing {
 public:
  FerroIsing(lattice::graph lat);
  double Propose(std::mt19937 &engine) {
    energy_proposed_ = energy_;
    site_ = dist_(engine);
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
  void Exchange(int partner, MPI_Comm local_comm);
  void StoreLog(std::ofstream *ofs_ptr);
  void SetFromLog(std::ifstream &ifs);
  // Gettor and settor.
  double val() const {return energy_;}
  size_t num_bins() const {return num_bins_;}
  size_t sweeps() const {return sweeps_;}
  double ene_min() const {return ene_min_;}
  double ene_max() const {return ene_max_;}
 private:
  const lattice::graph lat_;
  int site_;
  size_t num_bins_, sweeps_;
  double ene_min_, ene_max_, energy_proposed_, energy_;
  std::vector<int> spin_config_;
  std::uniform_int_distribution<size_t> dist_;
};


FerroIsing::FerroIsing(lattice::graph lat)
    : lat_(lat),
      sweeps_(lat.num_sites()),
      dist_(std::uniform_int_distribution<size_t>(0, lat_.num_sites()-1)) {
  ene_min_ = -2*(double)lat_.num_sites();
  ene_max_ = -ene_min_;
  num_bins_ = lat_.num_sites()+1;
  spin_config_.resize(lat_.num_sites(), 1);
  spin_config_.shrink_to_fit();
  energy_ = ene_min_;
}


void FerroIsing::Exchange(int partner, MPI_Comm local_comm) {
  MPI_Status status;
  MPI_Sendrecv_replace(&energy_, 1, MPI_DOUBLE, partner, 1, partner, 1,
      local_comm, &status);
  MPI_Sendrecv_replace(&spin_config_[0], spin_config_.size(), MPI_INT,
      partner, 2, partner, 2, local_comm, &status);
}


void FerroIsing::StoreLog(std::ofstream *ofs_ptr) {
  std::ofstream &ofs(*ofs_ptr);
  json log_json;
  // Main information.
  log_json["energy"] = energy_;
  json json_spin_config(spin_config_);
  log_json["spin_configuration"] = json_spin_config;
  // For checking.
  log_json["energy_min"] = ene_min_;
  log_json["energy_max"] = ene_max_;
  // Output with pretty printing.
  ofs << std::setw(4) << log_json << std::endl;
}


void FerroIsing::SetFromLog(std::ifstream &ifs) {
  json log_json;
  ifs >> log_json;
  energy_ = log_json["energy"].get<double>();
  spin_config_ = log_json["spin_configuration"].get<std::vector<int>>();
}


#endif // WANGLANDAU_FERRO_ISING_H_
#ifndef WANGLANDAU_LOG_FOR_JSON_H_
#define WANGLANDAU_LOG_FOR_JSON_H_


#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "mpi_setting.hpp"
#include "wl_params.hpp"
#include "flags.hpp"
using json = nlohmann::json;


inline void write_log_json(std::ofstream *ofs_ptr, RunningState running_state,
    const MPIV &mpiv, const WLParams &wl_params,
    const std::vector<double> &ln_dos, const std::mt19937 &engine,
    const std::vector<int> &histogram, int exch_count_down, double lnf_slowest);
inline bool check_log_json(const json &log_json, const MPIV &mpiv,
    const WLParams &wl_params, const std::vector<double> &ln_dos);
inline bool set_from_log_json(std::ifstream &ifs, MPIV *mpiv_ptr,
    WLParams *wl_params_ptr, std::vector<double> *ln_dos_ptr,
    std::mt19937 *engine_ptr, std::vector<int> *histogram_ptr,
    int *exch_count_down_ptr, double *lnf_slowest_ptr);


void write_log_json(std::ofstream *ofs_ptr, RunningState running_state,
    const MPIV &mpiv, const WLParams &wl_params,
    const std::vector<double> &ln_dos, const std::mt19937 &engine,
    const std::vector<int> &histogram, int exch_count_down,
    double lnf_slowest) {
  std::ofstream &ofs(*ofs_ptr);
  json log_json;
  // Mapping.
  log_json["last_time_state"] = static_cast<int>(running_state);
  log_json["mpi_setting"]["num_procs_tot"] = mpiv.numprocs();
  log_json["mpi_setting"]["num_procs_per_window"] = mpiv.num_walkers_window();
  log_json["mpi_setting"]["exch_pattern_id"] = mpiv.exch_pattern_id();
  log_json["wang_landau_params"]["check_flatness_every"] =
      wl_params.check_flatness_every();
  log_json["wang_landau_params"]["lnf"] = wl_params.lnf();
  log_json["wang_landau_params"]["lnf_min"] = wl_params.lnfmin();
  log_json["wang_landau_params"]["flatness"] = wl_params.flatness();
  log_json["wang_landau_params"]["overlap"] = wl_params.overlap();
  log_json["wang_landau_params"]["exch_every"] = wl_params.exch_every();
  json json_ln_dos(ln_dos);
  log_json["ln_dos"] = json_ln_dos;
  std::ostringstream oss;
  std::string engine_state;
  oss << engine;
  engine_state = oss.str();
  log_json["engine_state"] = engine_state;
  json json_histo(histogram);
  log_json["histogram"] = json_histo;
  log_json["exch_count_down"] = exch_count_down;
  log_json["lnf_slowest"] = lnf_slowest;
  // Output with pretty printing.
  ofs << std::setw(4) << log_json << std::endl;
}


bool check_log_json(const json &log_json, const MPIV &mpiv,
    const WLParams &wl_params, const std::vector<double> &ln_dos) {
  // Last state.
  if (log_json["last_time_state"].get<int>() ==
      static_cast<int>(RunningState::ALL_FINISHED)) return false;
  // MPI information.
  if (log_json["mpi_setting"]["num_procs_tot"].get<int>() !=
      mpiv.numprocs()) return false;
  if (log_json["mpi_setting"]["num_procs_per_window"].get<int>() !=
      mpiv.num_walkers_window()) return false;
  // WL parameters.
  if (log_json["wang_landau_params"]["check_flatness_every"].get<int>() !=
      wl_params.check_flatness_every()) return false;
  if (log_json["wang_landau_params"]["lnf_min"].get<double>() !=
      wl_params.lnfmin()) return false;
  if (log_json["wang_landau_params"]["flatness"].get<double>() !=
      wl_params.flatness()) return false;
  if (log_json["wang_landau_params"]["overlap"].get<double>() !=
      wl_params.overlap()) return false;
  if (log_json["wang_landau_params"]["exch_every"].get<int>() !=
      wl_params.exch_every()) return false;
  // Size of ln_dos.
  if (log_json["ln_dos"].get<std::vector<double>>().size() !=
      ln_dos.size()) return false;
  return true;
}


bool set_from_log_json(std::ifstream &ifs, MPIV *mpiv_ptr,
    WLParams *wl_params_ptr, std::vector<double> *ln_dos_ptr,
    std::mt19937 *engine_ptr, std::vector<int> *histogram_ptr,
    int *exch_count_down_ptr, double *lnf_slowest_ptr) {
  if (!ifs) return false;
  MPIV &mpiv(*mpiv_ptr);
  WLParams &wl_params(*wl_params_ptr);
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<int> &histogram(*histogram_ptr);
  std::mt19937 &engine(*engine_ptr);
  json log_json;
  ifs >> log_json;
  // Check consistency before setting.
  bool is_consistent;
  if (mpiv.myid() == 0) {
    is_consistent = check_log_json(log_json, mpiv, wl_params, ln_dos);
  }
  MPI_Barrier(MPI_COMM_WORLD); // Is this necessary?
  MPI_Bcast(&is_consistent, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
  if (!is_consistent) return false;
  // Set parameters from log file.
  mpiv.set_exch_pattern_id(
      log_json["mpi_setting"]["exch_pattern_id"].get<int>());
  wl_params.set_lnf(log_json["wang_landau_params"]["lnf"].get<double>());
  ln_dos = log_json["ln_dos"].get<std::vector<double>>();
  std::istringstream iss(log_json["engine_state"].get<std::string>());
  iss >> engine;
  histogram = log_json["histogram"].get<std::vector<int>>();
  *exch_count_down_ptr = log_json["exch_count_down"].get<int>();
  *lnf_slowest_ptr = log_json["lnf_slowest"].get<double>();
  return true;
}


#endif // WANGLANDAU_LOG_FOR_JSON_H_

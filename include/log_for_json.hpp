#ifndef WANGLANDAU_LOG_FOR_JSON_H_
#define WANGLANDAU_LOG_FOR_JSON_H_


#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include "mpi_setting.hpp"
#include "wl_params.hpp"
using json = nlohmann::json;


void write_log_json(std::ofstream *ofs_ptr, int running_state, const MPIV &mpiv,
    const WLParams &wl_params, const std::vector<double> &ln_dos,
    const std::mt19937 &engine, const std::vector<int> &histogram,
    int swap_count_down, int exchange_pattern, double lnf_slowest);
bool check_log_json(std::ifstream *ifs_ptr, const MPIV &mpiv,
    const WLParams &wl_params, const std::vector<double> &ln_dos);
void set_from_log_json(std::ifstream *ifs_ptr, WLParams *wl_params_ptr,
    std::vector<double> *ln_dos_ptr, std::mt19937 &engine,
    std::vector<int> *histogram_ptr, int *swap_count_down_ptr,
    int *exchange_pattern_ptr, double *lnf_slowest_ptr);


void write_log_json(std::ofstream *ofs_ptr, int running_state, const MPIV &mpiv,
    const WLParams &wl_params, const std::vector<double> &ln_dos,
    const std::mt19937 &engine, const std::vector<int> &histogram,
    int swap_count_down, int exchange_pattern, double lnf_slowest) {
  std::ofstream &ofs(*ofs_ptr);
  json log_json;
  // Mapping.
  log_json["last_time_state"] = running_state;
  log_json["mpi_setting"]["num_procs_tot"] = mpiv.numprocs();
  log_json["mpi_setting"]["num_procs_per_window"] = mpiv.multiple();
  log_json["wang_landau_params"]["sweeps"] = wl_params.sweeps();
  log_json["wang_landau_params"]["check_flatness_every"] =
      wl_params.check_flatness_every();
  log_json["wang_landau_params"]["lnf"] = wl_params.lnf();
  log_json["wang_landau_params"]["lnf_min"] = wl_params.lnfmin();
  log_json["wang_landau_params"]["flatness"] = wl_params.flatness();
  log_json["wang_landau_params"]["swap_every"] = wl_params.swap_every();
  json json_ln_dos(ln_dos);
  log_json["ln_dos"] = json_ln_dos;
  std::ostringstream oss;
  std::string engine_state;
  oss << engine;
  engine_state = oss.str();
  log_json["engine_state"] = engine_state;
  json json_histo(histogram);
  log_json["histogram"] = json_histo;
  log_json["swap_count_down"] = swap_count_down;
  log_json["exchange_pattern"] = exchange_pattern;
  log_json["lnf_slowest"] = lnf_slowest;
  // Output with pretty printing.
  ofs << std::setw(4) << log_json << std::endl;
}


bool check_log_json(std::ifstream *ifs_ptr, const MPIV &mpiv,
    const WLParams &wl_params, const std::vector<double> &ln_dos) {
  std::ifstream &ifs(*ifs_ptr);
  if (!ifs) return false;
  json log_json;
  ifs >> log_json;
  // Last state.
  if (log_json["last_time_state"].get<int>() == 1) return false;
  // MPI information.
  if (log_json["mpi_setting"]["num_procs_tot"].get<int>() !=
      mpiv.numprocs()) return false;
  if (log_json["mpi_setting"]["num_procs_per_window"].get<int>() !=
      mpiv.multiple()) return false;
  // WL parameters.
  if (log_json["wang_landau_params"]["sweeps"].get<int>() !=
      wl_params.sweeps()) return false;
  if (log_json["wang_landau_params"]["check_flatness_every"].get<int>() !=
      wl_params.check_flatness_every()) return false;
  if (log_json["wang_landau_params"]["lnf_min"].get<double>() !=
      wl_params.lnfmin()) return false;
  if (log_json["wang_landau_params"]["flatness"].get<double>() !=
      wl_params.flatness()) return false;
  if (log_json["wang_landau_params"]["swap_every"].get<int>() !=
      wl_params.swap_every()) return false;
  // Size of ln_dos.
  if (log_json["ln_dos"].get<std::vector<double>>().size() !=
      ln_dos.size()) return false;
  return true;
}


void set_from_log_json(std::ifstream *ifs_ptr, WLParams *wl_params_ptr,
    std::vector<double> *ln_dos_ptr, std::mt19937 &engine,
    std::vector<int> *histogram_ptr, int *swap_count_down_ptr,
    int *exchange_pattern_ptr, double *lnf_slowest_ptr) {
  std::ifstream &ifs(*ifs_ptr);
  WLParams &wl_params(*wl_params_ptr);
  std::vector<double> &ln_dos(*ln_dos_ptr);
  std::vector<int> &histogram(*histogram_ptr);
  json log_json;
  ifs >> log_json;
  wl_params.set_lnf(log_json["wang_landau_params"]["lnf"].get<double>());
  ln_dos = log_json["ln_dos"].get<std::vector<double>>();
  std::istringstream iss(log_json["engine_state"].get<std::string>());
  iss >> engine;
  histogram = log_json["histogram"].get<std::vector<int>>();
  *swap_count_down_ptr = log_json["swap_count_down"].get<int>();
  *exchange_pattern_ptr = log_json["exchange_pattern"].get<int>();
  *lnf_slowest_ptr = log_json["lnf_slowest"].get<double>();
}


#endif // WANGLANDAU_LOG_FOR_JSON_H_

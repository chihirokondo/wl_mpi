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


#endif // WANGLANDAU_LOG_FOR_JSON_H_

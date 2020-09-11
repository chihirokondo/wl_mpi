#include <iostream>
#include "../include/lattice/graph.hpp"
#include "../include/ferro_ising.hpp"
#include "../include/histo_env_manager.hpp"


int main() {
  double min = 0;
  double max = 4;
  size_t num_bins = 5;
  bool centering = true;
  double val = 4.3;
  HistoEnvManager test(min, max, num_bins, centering);
  std::cout << "min=" << min << ", max=" << max << ", num_bins=" << num_bins
      << ", centering=" << centering << std::endl;
  std::cout << "central values: {";
  for (int i=0; i<test.num_bins(); ++i) {
     std::cout << test.GetCentralVal(i) << ", ";
  }
  std::cout << "}" << std::endl;
  std::cout << "get index(val=" << val << "): ";
  std::cout << test.GetIndex(val) << std::endl;
  std::cout << std::endl;

  num_bins = 4;
  centering = false;
  test = HistoEnvManager::Construct(min, max, num_bins, centering);
  std::cout << "min=" << min << ", max=" << max << ", num_bins=" << num_bins
      << ", centering=" << centering << std::endl;
  std::cout << "central values : {";
  for (int i=0; i<test.num_bins(); ++i) {
      std::cout << test.GetCentralVal(i) << ", ";
  }
  std::cout << "}" << std::endl;
  val = 0.01;
  std::cout << "get index(val=" << val << "): ";
  std::cout << test.GetIndex(val) << std::endl;
  return 0;
}
#include <iostream>
#include "../include/lattice/graph.hpp"

int main() {
  lattice::graph lat = lattice::graph::simple(2, 2);
  for (int i=0; i<100; ++i) {
    std::cout << "site : " << i%4 << std::endl;
    for (int j=0; j<lat.num_neighbors(i%4); ++j) {
      std::cout << lat.neighbor(i%4, j) << std::endl;
    }
  }
  return 0;
}

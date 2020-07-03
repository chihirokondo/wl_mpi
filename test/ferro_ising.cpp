#include <iostream>
#include "../include/lattice/graph.hpp"
#include "../include/ferro_ising.hpp"


int main() {
  int dim = 2;
  int length = 4;
  lattice::graph lat = lattice::graph::simple(dim, length);
  FerroIsing model(lat);
}

// Ferro ising.
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "../lattice/graph.hpp"


int main() {
  int dim = 2;
  int length = 4;
  lattice::graph lat = lattice::graph::simple(dim, length);

  std::vector<int> dos((int)lat.num_bonds()/2+1);
  int energy, source_state, target_state;

  for (int state=0; state<(int)std::pow(2.0, (double)lat.num_sites()); ++state) {
    energy = 0;
    for (int site=0; site<lat.num_sites(); ++site) {
      source_state = (state & (int)std::pow(2.0, (double)site)) >> site;
      for(int n=0; n<lat.num_neighbors(site); ++n) {
        target_state = (state & (int)std::pow(2.0, (double)lat.neighbor(site, n))) >>
            (int)lat.neighbor(site, n);
        energy += (int)(-4.0*((double)source_state-0.5)*((double)target_state-0.5));
      }
    }
    energy /= 2;
    dos[(energy+lat.num_bonds())/4] += 1;
  }

  // output.
  std::ofstream outputfile;
  outputfile.open("exact.dat", std::ios::out);
  outputfile << "# dim = " << dim << ", length = " << length << "\n";
  outputfile << "# energy \t # dos\n";
  for (int i=0; i<dos.size(); ++i) {
    outputfile << 4*i-(int)lat.num_bonds() << "\t" << dos[i] << "\n";
  }
  outputfile << std::endl;
  outputfile.close();

  return 0;
}

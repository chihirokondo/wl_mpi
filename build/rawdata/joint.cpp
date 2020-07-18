#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include "include/rewl_data.hpp"
using namespace std;


int main() {
  // Read file.
  ifstream ifs("../lngE_proc0.dat", ios::in);
  try {
    if (!ifs) throw 0;
  }
  catch (int err_status) {
    cerr << "ERROR: Failed to open file." << endl;
    return 1;
  }
  REWLData data(ifs);
  // Routine for joint.
  string filename;
  REWLData data_append;
  for (int i=1; i<data.num_windows(); ++i) {
    filename = "../lngE_proc" + to_string(i) + ".dat";
    ifstream ifs(filename, ios::in);
    data_append = REWLData::construct_data(ifs);
    // Determine joint point.
    double min_diff =
        fabs(data.ln_dos_diff(data.diff_sew())-data_append.ln_dos_diff(0));
    int joint_diff = data.diff_sew();
    for (int j=data.diff_sew()+1; j<data.ln_dos_diff_size(); ++j) {
      if (min_diff > fabs(data.ln_dos_diff(j) -
          data_append.ln_dos_diff(j-data.diff_sew()))) {
        joint_diff = j;
      }
    }
    // Joint.
    data.Joint(joint_diff, data_append);
  }
  // ln_dos to bare dos.
  double min_ln_dos = 0.0;
  for (int i=0; i<data.ln_dos_size(); ++i) {
    if (data.ln_dos(i) != 0.0) {
      if ((min_ln_dos == 0.0) | (min_ln_dos > data.ln_dos(i))) {
        min_ln_dos = data.ln_dos(i);
      }
    }
  }
  double sum_dos_before = 0.0;
  vector<double> dos;
  for (int i=0; i<data.ln_dos_size(); ++i) {
    if (data.ln_dos(i) == 0.0) {
      dos.push_back(0.0);
    } else {
      dos.push_back(exp(data.ln_dos(i)-min_ln_dos));
      sum_dos_before += exp(data.ln_dos(i)-min_ln_dos);
    }
  }
  // Normalization.
  for (double &dos_i : dos) {
    dos_i *= data.constraint_condition_value() / sum_dos_before;
  }
  // Output.
  ofstream ofs("../final_dos.dat", ios::out);
  ofs << data.info() << "\n";
  ofs << "# energy \t # density of states\n";
  ofs << scientific << setprecision(15);
  for (int i=0; i<dos.size(); ++i) {
    ofs << data.energy(i) << "\t";
    ofs << dos[i] << "\n";
  }
  ofs << endl;
  cout << "Results:\n" << "{";
  for (const double &x : dos) cout << x << ", ";
  cout << "}" << endl;
  return 0;
}
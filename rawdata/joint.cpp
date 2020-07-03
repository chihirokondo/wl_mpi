#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "include/rewl_data.hpp"
using namespace std;


int main(int argc, char *argv[]) {
  try {
    if (argc != 2) throw 0;
  }
  catch (int err_status) {
    cerr
        << "ERROR: Unexpected number of command line arguments!\n"
        << "       Expect 1 argument, " << argc-1 << " were provided.\n"
        << "Syntax: ./a.out [arg1] \n\n"
        << "Please provide the following command line arguments:\n"
        << "1. Number of files. [integer]\n"
        << endl;
    return 0;
  }
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

  string filename;
  REWLData data_append;
  for (int i=1; i<atoi(argv[1]); ++i) {
    filename = "../lngE_proc" + to_string(i) + ".dat";
    ifstream ifs(filename, ios::in);
    data_append = REWLData::construct_data(ifs);
    // Determine joint point.
    double min_diff =
        fabs(data.ln_dos_diff_[data.diff_sew_]-data_append.ln_dos_diff_[0]);
    int joint_diff = data.diff_sew_;
    for (int j=data.diff_sew_+1; j<data.ln_dos_diff_.size(); ++j) {
      if (min_diff > fabs(data.ln_dos_diff_[j] -
          data_append.ln_dos_diff_[j-data.diff_sew_])) {
        joint_diff = j;
      }
    }
    // Joint.
    data.Joint(joint_diff, data_append);
  }

  // ln_dos to bare dos.
  double min_ln_dos = 0.0;
  for (const double &ln_dos_i : data.ln_dos_) {
    if (ln_dos_i != 0.0) {
      if ((min_ln_dos == 0.0) | (min_ln_dos > ln_dos_i)) min_ln_dos = ln_dos_i;
    }
  }
  double sum_dos_before = 0.0;
  vector<double> dos;
  for (const double &ln_dos_i : data.ln_dos_) {
    if (ln_dos_i == 0.0) {
      dos.push_back(0.0);
    } else {
      dos.push_back(exp(ln_dos_i-min_ln_dos));
      sum_dos_before += exp(ln_dos_i-min_ln_dos);
    }
  }
  // Normalization.
  for (double &dos_i : dos) {
    dos_i *= data.constraint_condition_value / sum_dos_before;
  }

  // Output.
  ofstream ofs("final_dos.dat", ios::out);
  ofs << data.info << "\n";
  ofs << "# energy \t # density of states\n";
  for (int i=0; i<dos.size(); ++i) {
    ofs << data.energy_[i] << "\t" << dos[i] << "\n";
  }
  ofs << endl;

  ////
  cout << "{";
  for (const double &x : dos) cout << x << ", ";
  cout << "}" << endl;
  ////
  return 0;
}
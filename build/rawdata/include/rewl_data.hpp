#ifndef WANGLANDAU_REWL_DATA_H_
#define WANGLANDAU_REWL_DATA_H_


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;


class REWLData {
 public:
  REWLData() {}
  REWLData(ifstream &ifs);
  static REWLData construct_data(ifstream &ifs) {return REWLData(ifs);}
  void Joint(int joint_diff, REWLData data_append);
  int diff_sew() {return diff_sew_;}
  int num_windows() {return num_windows_;}
  string info() {return info_;}
  double constraint_condition_GetVal() {return constraint_condition_value_;}
  double ln_dos(int index) {return ln_dos_[index];}
  size_t ln_dos_size() {return ln_dos_.size();}
  double ln_dos_diff(int index) {return ln_dos_diff_[index];}
  size_t ln_dos_diff_size() {return ln_dos_diff_.size();}
  double energy(int index) {return energy_[index];}
 private:
  int sew_, diff_sew_, num_windows_;
  string info_, constraint_condition_type_;
  double constraint_condition_value_;
  vector<double> energy_, ln_dos_, ln_dos_diff_;
  vector<int> to_original_i_;
};


REWLData::REWLData (ifstream &ifs) {
  energy_ = {};
  ln_dos_ = {};
  double energy_i, ln_dos_i;
  string str;
  getline(ifs, info_);
  ifs >> str >> str >> constraint_condition_type_;
  ifs >> str >> str >> constraint_condition_value_;
  ifs >> str >> str >> sew_;
  ifs >> str >> str >> num_windows_;
  ifs >> str >> str >> str >> str;
  while (ifs >> energy_i >> ln_dos_i) {
    energy_.push_back(energy_i);
    ln_dos_.push_back(ln_dos_i);
  }
  diff_sew_ = 0;
  ln_dos_diff_ = {};
  to_original_i_ = {};
  double ln_dos_old = 0.0;
  double energy_old;
  int i_old;
  for (int i=0; i<ln_dos_.size(); ++i) {
    if (ln_dos_[i] != 0.0) {
      if (ln_dos_old != 0.0) {
        if ((i>sew_) && (diff_sew_==0)) diff_sew_ = ln_dos_diff_.size();
        to_original_i_.push_back(i_old);
        ln_dos_diff_.push_back((ln_dos_[i]-ln_dos_old)/(energy_[i]-energy_old));
      }
      i_old = i;
      energy_old = energy_[i];
      ln_dos_old = ln_dos_[i];
    }
  }
}


void REWLData::Joint (int joint_diff, REWLData data_append) {
  int joint_point = to_original_i_[joint_diff];
  double modi_height =
      ln_dos_[joint_point] - data_append.ln_dos_[joint_point-sew_];
  for (int i=ln_dos_.size()-1; i>=joint_point; --i) {
    energy_.pop_back();
    ln_dos_.pop_back();
  }
  for (int i=joint_point-sew_; i<data_append.ln_dos_.size(); ++i) {
    energy_.push_back(data_append.energy_[i]);
    if (data_append.ln_dos_[i] == 0.0) {
      ln_dos_.push_back(0.0);
    } else {
      ln_dos_.push_back(data_append.ln_dos_[i] + modi_height);
    }
  }
  sew_ += data_append.sew_;

  for (int i=ln_dos_diff_.size()-1; i>=joint_diff; --i) {
    to_original_i_.pop_back();
    ln_dos_diff_.pop_back();
  }
  diff_sew_ = 0;
  double ln_dos_old = ln_dos_[joint_point];
  double energy_old = energy_[joint_point];
  int i_old = joint_point;
  for (int i=joint_point+1; i<ln_dos_.size(); ++i) {
    if (ln_dos_[i] != 0.0) {
      if ((i>sew_) && (diff_sew_==0)) diff_sew_ = ln_dos_diff_.size();
      to_original_i_.push_back(i_old);
      ln_dos_diff_.push_back((ln_dos_[i]-ln_dos_old)/(energy_[i]-energy_old));
      i_old = i;
      energy_old = energy_[i];
      ln_dos_old = ln_dos_[i];
    }
  }
}


#endif // WANGLANDAU_REWL_DATA_H_
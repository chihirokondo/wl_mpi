#ifndef WANGLANDAU_MPI_H_
#define WANGLANDAU_MPI_H_


#include <mpi.h>
#include <iostream>
#include <vector>


class MPIV {
 public:
  MPIV(int numprocs, int myid, int num_walkers_window);
  void switch_exch_pattern() {
    if (num_windows_>2) exch_pattern_id_ ^= 1;
    if ((exch_pattern_id_==0) && belong_exchblock_pattern_0_) {
      comm_id_ = 2*(myid_/(2*num_walkers_window_));
    } else if ((exch_pattern_id_==1) && belong_exchblock_pattern_1_) {
      comm_id_ = ((myid_-num_walkers_window_)/(2*num_walkers_window_))*2 + 1;
    } else {
      comm_id_ = not_belong_any_exchblock_;
    }
  }
  bool have_exch_partner() {
    if (comm_id_ == not_belong_any_exchblock_) return false;
    return true;
  }
  void set_exch_pattern_id(int exch_pattern_id) {
    exch_pattern_id_ = exch_pattern_id;
  }
  // Gettor.
  int numprocs() const {return numprocs_;}
  int myid() const {return myid_;}
  int num_walkers_window() const {return num_walkers_window_;}
  int num_windows() const {return num_windows_;}
  int exch_pattern_id() const {return exch_pattern_id_;}
  int id_in_exchblock() const {return id_in_exchblock_[exch_pattern_id_];}
  MPI_Comm mpi_comm_exchblock() const {return mpi_comm_exchblock_[comm_id_];}
 private:
  const int numprocs_, myid_, num_walkers_window_, num_windows_;
  const int not_belong_any_exchblock_ = -1;
  bool belong_exchblock_pattern_0_, belong_exchblock_pattern_1_;
  int comm_id_, exch_pattern_id_;
  std::vector<int> id_in_exchblock_;
  std::vector<MPI_Comm> mpi_comm_exchblock_;
};


inline MPIV::MPIV(int numprocs, int myid, int num_walkers_window)
    : numprocs_(numprocs),
      myid_(myid),
      num_walkers_window_(num_walkers_window),
      num_windows_(numprocs/num_walkers_window) {
  exch_pattern_id_ = 0; // 0 or 1.
  bool num_windows_is_odd = num_windows_%2 == 1;
  bool num_windows_is_even = num_windows_%2 == 0;
  bool belong_first_window = myid_ < num_walkers_window_;
  bool belong_last_window = myid_ >= numprocs_-num_walkers_window_;
  belong_exchblock_pattern_0_ = !(belong_last_window && num_windows_is_odd);
  belong_exchblock_pattern_1_ = !(belong_first_window ||
      (belong_last_window && num_windows_is_even));
  // Prepare exchange envionment.
  if (num_windows_ > 1) {
    // Create new communicators.
    MPI_Group world;
    MPI_Comm_group(MPI_COMM_WORLD, &world);
    std::vector<int> ranks(2*num_walkers_window_);
    std::vector<MPI_Group> local_group(num_windows_-1);
    mpi_comm_exchblock_.resize(num_windows_-1);
    for (int i=0; i<local_group.size(); ++i) {
      for (int j=0; j<2*num_walkers_window_;++j) {
        // Get the process rank in MPI_COMM_WORLD.
        ranks[j] = i*num_walkers_window_+j;
      }
      MPI_Group_incl(world, 2*num_walkers_window_, &ranks[0], &local_group[i]);
      MPI_Comm_create(MPI_COMM_WORLD, local_group[i], &mpi_comm_exchblock_[i]);
    }
    // Get the local id (in the local communicators).
    id_in_exchblock_ = {not_belong_any_exchblock_, not_belong_any_exchblock_};
    int comm_id;
    // CASE: "exch_pattern_id_" == 0.
    if (belong_exchblock_pattern_0_) {
      comm_id = 2*(myid_/(2*num_walkers_window_));
      MPI_Comm_rank(mpi_comm_exchblock_[comm_id], &id_in_exchblock_[0]);
    }
    // CASE: "exch_pattern_id" == 1.
    if (belong_exchblock_pattern_1_) {
      comm_id = 2*((myid_-num_walkers_window_)/(2*num_walkers_window_)) + 1;
      MPI_Comm_rank(mpi_comm_exchblock_[comm_id], &id_in_exchblock_[1]);
    }
  }
}


#endif // WANGLANDAU_MPI_H_

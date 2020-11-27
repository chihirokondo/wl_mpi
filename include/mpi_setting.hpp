#ifndef WANGLANDAU_MPI_H_
#define WANGLANDAU_MPI_H_


#include <mpi.h>
#include <iostream>
#include <vector>


class MPIV {
 public:
  MPIV(int numprocs, int myid, int num_walkers_window);
  void create_local_communicator();
  void set_id_in_exchblock();
  void set_comm_id(int exch_pattern_id);
  // Gettor.
  int numprocs() const {return numprocs_;}
  int myid() const {return myid_;}
  int num_walkers_window() const {return num_walkers_window_;}
  int num_windows() const {return num_windows_;}
  int comm_id() const {return comm_id_;}
  int id_in_exchblock(int exch_pattern_id) const {
    return id_in_exchblock_[exch_pattern_id];
  }
  MPI_Comm mpi_comm_exchblock() const {return mpi_comm_exchblock_[comm_id_];}
 private:
  const int numprocs_, myid_, num_walkers_window_, num_windows_,
      num_windows_mod2_;
  int comm_id_;
  std::vector<int> id_in_exchblock_;
  std::vector<MPI_Comm> mpi_comm_exchblock_;
};


MPIV::MPIV(int numprocs, int myid, int num_walkers_window)
    : numprocs_(numprocs),
      myid_(myid),
      num_walkers_window_(num_walkers_window),
      num_windows_(numprocs/num_walkers_window),
      num_windows_mod2_((numprocs/num_walkers_window)%2) {}


void MPIV::create_local_communicator() {
  MPI_Group world;
  MPI_Comm_group(MPI_COMM_WORLD, &world);
  std::vector<int> ranks(2*num_walkers_window_);
  std::vector<MPI_Group> local_group(num_windows_-1);
  std::vector<MPI_Comm> mpi_comm_exchblock_buf(num_windows_-1);
  for (int i=0; i<local_group.size(); ++i) {
    for (int j=0; j<2*num_walkers_window_;++j) {
      // Get the process rank in MPI_COMM_WORLD.
      ranks[j] = i*num_walkers_window_+j;
    }
    MPI_Group_incl(world, 2*num_walkers_window_, &ranks[0], &local_group[i]);
    MPI_Comm_create(MPI_COMM_WORLD, local_group[i], &mpi_comm_exchblock_buf[i]);
  }
  mpi_comm_exchblock_ = mpi_comm_exchblock_buf;
}


void MPIV::set_id_in_exchblock() {
  id_in_exchblock_ = {-1, -1};
  int comm_id_tmp;
  if (!((num_windows_mod2_==1)&&(myid_>=numprocs_-num_walkers_window_))) {
    comm_id_tmp = 2*(myid_/(2*num_walkers_window_));
    MPI_Comm_rank(mpi_comm_exchblock_[comm_id_tmp], &id_in_exchblock_[0]);
  }
  if (myid_>=num_walkers_window_) {
    if (!((num_windows_mod2_==0)&&(myid_>=numprocs_-num_walkers_window_))) {
      comm_id_tmp = 2*((myid_-num_walkers_window_)/(2*num_walkers_window_)) + 1;
      MPI_Comm_rank(mpi_comm_exchblock_[comm_id_tmp], &id_in_exchblock_[1]);
    }
  }
}


void MPIV::set_comm_id(int exch_pattern_id) {
  if ((exch_pattern_id==0) &&
      !((num_windows_mod2_==1)&&(myid_>=numprocs_-num_walkers_window_))) {
    comm_id_ = 2*(myid_/(2*num_walkers_window_));
  } else if ((exch_pattern_id==1)&&(myid_>=num_walkers_window_)&&
      !((num_windows_mod2_==0)&&(myid_>=numprocs_-num_walkers_window_))) {
    comm_id_ = ((myid_-num_walkers_window_)/(2*num_walkers_window_))*2 + 1;
  } else {
    comm_id_ = -1;
  }
}


#endif // WANGLANDAU_MPI_H_

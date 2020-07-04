#ifndef WANGLANDAU_MPI_H_
#define WANGLANDAU_MPI_H_


#include <vector>


class MPIV {
 public:
  MPIV(int numprocs, int myid, int multiple);
  void create_local_communicator();
  void set_local_id();
  void set_comm_id(int exchange_pattern);
  // Gettor.
  int numprocs() const {return numprocs_;}
  int myid() const {return myid_;}
  int multiple() const {return multiple_;}
  int comm_id() const {return comm_id_;}
  int local_id(int index) const {return local_id_[index];}
  MPI_Comm local_comm(int index) const {return local_comm_[index];}
 private:
  const int numprocs_, myid_, multiple_, num_window_mod2_;
  int comm_id_;
  std::vector<int> local_id_;
  std::vector<MPI_Comm> local_comm_;
};


MPIV::MPIV(int numprocs, int myid, int multiple)
    : numprocs_(numprocs),
      myid_(myid),
      multiple_(multiple),
      num_window_mod2_((numprocs/multiple)%2) {}


void MPIV::create_local_communicator() {
  MPI_Group world;
  MPI_Comm_group(MPI_COMM_WORLD, &world);
  std::vector<int> ranks(2*multiple_);
  std::vector<MPI_Group> local_group((numprocs_/multiple_)-1);
  std::vector<MPI_Comm> local_comm_buf((numprocs_/multiple_)-1);
  for (int i=0; i<local_group.size(); ++i) {
    for (int j=0; j<2*multiple_;++j) {
      // Get the process rank in MPI_COMM_WORLD.
      ranks[j] = i*multiple_+j;
    }
    MPI_Group_incl(world, 2*multiple_, &ranks[0], &local_group[i]);
    MPI_Comm_create(MPI_COMM_WORLD, local_group[i], &local_comm_buf[i]);
  }
  local_comm_ = local_comm_buf;
}


void MPIV::set_local_id() {
  local_id_ = {-1, -1};
  int comm_id_tmp;
  if (!((num_window_mod2_==1)&&(myid_>=numprocs_-multiple_))) {
    comm_id_tmp = 2*(myid_/(2*multiple_));
    MPI_Comm_rank(local_comm_[comm_id_tmp], &local_id_[0]);
  }
  if (myid_>=multiple_) {
    if (!((num_window_mod2_==0)&&(myid_>numprocs_-multiple_))) {
      comm_id_tmp = 2*((myid_-multiple_)/(2*multiple_)) + 1;
      MPI_Comm_rank(local_comm_[comm_id_tmp], &local_id_[1]);
    }
  }
}


void MPIV::set_comm_id(int exchange_pattern) {
  if ((exchange_pattern==0) &&
      !((num_window_mod2_==1)&&(myid_>=numprocs_-multiple_))) {
    comm_id_ = 2*(myid_/(2*multiple_));
  } else if ((exchange_pattern==1)&&(myid_>=multiple_)&&
      !((num_window_mod2_==0)&&(myid_>numprocs_-multiple_))) {
    comm_id_ = ((myid_-multiple_)/(2*multiple_))*2 + 1;
  } else {
    comm_id_ = -1;
  }
}


#endif // WANGLANDAU_MPI_H_
#ifndef WANGLANDAU_MPI_H_
#define WANGLANDAU_MPI_H_


#include <mpi.h>
#include <iostream>
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
  int num_windows() const {return num_windows_;}
  int comm_id() const {return comm_id_;}
  int local_id(int index) const {return local_id_[index];}
  MPI_Comm local_comm(int index) const {return local_comm_[index];}
 private:
  const int numprocs_, myid_, multiple_, num_windows_, num_windows_mod2_;
  int comm_id_;
  std::vector<int> local_id_;
  std::vector<MPI_Comm> local_comm_;
};


#endif // WANGLANDAU_MPI_H_

#include "domain_2.h"

SquarePartition::SquarePartition(int proc_nx, int proc_ny, MPI_Comm group_comm)
  : comm_(group_comm){
  MPI_Comm_rank(comm_, &my_rank_);
  int my_col = my_rank_ % proc_nx;
  int my_row = my_rank_ / proc_nx;

  if (proc_nx > 1) {
    neighbor_[0][0] = (my_col - 1 + proc_nx) % proc_nx + my_row * proc_nx;
    neighbor_[0][1] = (my_col + 1) % proc_nx + my_row * proc_nx;
  } else {
    neighbor_[0][0] = neighbor_[0][1] = MPI_PROC_NULL;
  }

  if (proc_ny > 1) {
    neighbor_[1][0] = my_col + (my_row - 1 + proc_ny) % proc_ny * proc_nx;
    neighbor_[1][1] = my_col + (my_row + 1) % proc_ny * proc_nx;
  } else {
    neighbor_[1][0] = neighbor_[1][1] = MPI_PROC_NULL;
  }
}

#pragma once
#include <mpi.h>

class SquarePartition {
public:
  SquarePartition(int nx, int ny, MPI_Comm group_comm);

  MPI_Comm get_communicator() const { return comm_; }

  bool is_root() const { return my_rank_ == 0; }

  int get_my_rank() const { return my_rank_; }


  template <typename T>
  void copy_to_padded_region(T* field, int Lx, int Ly,
                             bool x_padded, bool y_padded,
                             MPI_Datatype row_type, MPI_Datatype col_type) const;

private:
  MPI_Comm comm_;
  int my_rank_;
  int neighbor_[2][2]{};
};

template <typename T>
void SquarePartition::copy_to_padded_region(T* field, int Lx, int Ly,
                                            bool x_padded, bool y_padded,
                                            MPI_Datatype row_type, MPI_Datatype col_type) const {
  MPI_Request req[4];
  MPI_Status stat[4];

  if (y_padded) {
    T* row_recv0 = field + Lx * (Ly - 1);
    T* row_send0 = field + Lx * 1;
    T* row_recv1 = field + Lx * 0;
    T* row_send1 = field + Lx * (Ly - 2);

    //! Transfer data downward
    MPI_Irecv(row_recv0, 1, row_type, neighbor_[1][1], 21, comm_, &req[0]);
    MPI_Isend(row_send0, 1, row_type, neighbor_[1][0], 21, comm_, &req[1]);

    //! Transfer data upward
    MPI_Irecv(row_recv1, 1, row_type, neighbor_[1][0], 12, comm_, &req[2]);
    MPI_Isend(row_send1, 1, row_type, neighbor_[1][1], 12, comm_, &req[3]);

    MPI_Waitall(4, req, stat);

  }

  if (x_padded) {
    T* col_recv0 = field + Lx - 1;
    T* col_send0 = field + 1;
    T* col_recv1 = field + 0;
    T* col_send1 = field + Lx - 2;

    //! Transfer data leftward
    MPI_Irecv(col_recv0, 1, col_type, neighbor_[0][1], 21, comm_, &req[0]);
    MPI_Isend(col_send0, 1, col_type, neighbor_[0][0], 21, comm_, &req[1]);

    //! Transfer data rightward
    MPI_Irecv(col_recv1, 1, col_type, neighbor_[0][0], 12, comm_, &req[2]);
    MPI_Isend(col_send1, 1, col_type, neighbor_[0][1], 12, comm_, &req[3]);

    MPI_Waitall(4, req, stat);
  }
}

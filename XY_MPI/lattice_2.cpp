#include "lattice_2.h"

SquareLatticePadded_2::SquareLatticePadded_2(int gl_Lx, int gl_Ly, 
                                             int proc_nx, int proc_ny,
                                             MPI_Comm group_comm)
  : gl_Lx_(gl_Lx), gl_Ly_(gl_Ly),
    SquarePartition(proc_nx, proc_ny, group_comm),
    x_padded_(proc_nx > 1), y_padded_(proc_ny > 1) {
  int real_Lx = gl_Lx / proc_nx;
  int real_Ly = gl_Ly / proc_ny;
  Lx_ = real_Lx;
  Ly_ = real_Ly;
  real_beg_x_ = 0;
  real_beg_y_ = 0;

  if (x_padded_) {
    Lx_ += 2;
    real_beg_x_ += 1;
  }
  if (y_padded_) {
    Ly_ += 2;
    real_beg_y_ += 1;
  }

  real_end_x_ = real_beg_x_ + real_Lx;
  real_end_y_ = real_beg_y_ + real_Ly;
  n_ = Lx_ * Ly_;

  MPI_Type_vector(Ly_, 1, Lx_, MPI_DOUBLE, &col_double_type_);
  MPI_Type_vector(Lx_, 1, 1, MPI_DOUBLE, &row_double_type_);

  MPI_Type_commit(&col_double_type_);
  MPI_Type_commit(&row_double_type_);

  MPI_Datatype gl_block_type2;
  MPI_Type_vector(real_Ly, real_Lx, gl_Lx, MPI_DOUBLE, &gl_block_type2);

  //std::cout << "size of double: " << sizeof(double) << std::endl;
  MPI_Type_create_resized(gl_block_type2, 0, sizeof(double), &gl_block_type_);
  MPI_Type_commit(&gl_block_type_);

  MPI_Type_vector(real_Ly, real_Lx, Lx_, MPI_DOUBLE, &my_block_type_);
  MPI_Type_commit(&my_block_type_);

  //std::cout << real_Ly << "\t" << real_Lx << "\t" << gl_Lx << std::endl;


  int np = proc_nx * proc_ny;
  counts_ = new int[np];
  displs_ = new int[np];

  for (int row = 0; row < proc_ny; row++) {
    for (int col = 0; col < proc_nx; col++) {
      int idx = col + row * proc_nx;
      counts_[idx] = 1;
      displs_[idx] = real_Lx * col + gl_Lx * row * real_Ly;
      //std::cout << "dis: " << displs_[idx] << std::endl;
    }
  }

}

SquareLatticePadded_2::~SquareLatticePadded_2(){
  MPI_Type_free(&col_double_type_);
  MPI_Type_free(&row_double_type_);
  MPI_Type_free(&gl_block_type_);
  MPI_Type_free(&my_block_type_);
  delete[] counts_;
  delete[] displs_;
}


void SquareLatticePadded_2::scatter_fields(double* gl_f, double* my_f) const {
  int root = 0;
  int real_beg = get_idx(real_beg_x_, real_beg_y_);
  
  //if (is_root()) {
  //  for (int row = 0; row < gl_Ly_; row++) {
  //    for (int col = 0; col < gl_Lx_; col++) {
  //      int idx = col + row * gl_Lx_;
  //      std::cout << gl_f[idx] << " ";
  //    }
  //    std::cout << std::endl;
  //  }

  //  std::cout << "dis: ";
  //  for (int i = 0; i < 4; i++) {
  //    std::cout << displs_[i] << "\t";
  //  }
  //  std::cout << std::endl;
  //}

  MPI_Scatterv(gl_f, counts_, displs_, gl_block_type_,
    &my_f[real_beg], 1, my_block_type_, root, get_communicator());


  //if (get_my_rank() == 2) {
  //  for (int row = 0; row < Ly_; row++) {
  //    for (int col = 0; col < Lx_; col++) {
  //      int idx = col + row * Lx_;
  //      std::cout << my_f[idx] << " ";
  //    }
  //    std::cout << std::endl;
  //  }
  //}

}

void SquareLatticePadded_2::gather_fields(double* gl_f, double* my_f) const {
  int root = 0;
  int real_beg = get_idx(real_beg_x_, real_beg_y_);
  MPI_Gatherv(&my_f[real_beg], 1, my_block_type_,
    gl_f, counts_, displs_, gl_block_type_, root, get_communicator());
}
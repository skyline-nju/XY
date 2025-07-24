#pragma once
#include <vector>
#include <mpi.h>
#include "comn.h"
#include "domain_2.h"

class SquareLatticePadded_2: public SquarePartition {
public:
  SquareLatticePadded_2(int gl_Lx, int gl_Ly, int proc_nx, int proc_ny, MPI_Comm group_comm);

  ~SquareLatticePadded_2();

  void comm_before_cal_force(double* theta) {
    copy_to_padded_region(theta, Lx_, Ly_, x_padded_, y_padded_, row_double_type_, col_double_type_);
  }

  template <typename BiFunc>
  void for_each_pair(BiFunc pair_force) const;

  template <typename UniFunc>
  void for_each_real_site(UniFunc func) const;

  void scatter_fields(const double* gl_f, double* my_f) const;

  void gather_fields(double* gl_f, const double* my_f) const;

  int get_gl_Lx() const { return gl_Lx_; }
  int get_gl_Ly() const { return gl_Ly_; }
  int get_Lx() const { return Lx_; }
  int get_Ly() const { return Ly_; }
  int get_n() const { return n_; }
  int get_real_n() const { return (real_end_x_ - real_beg_x_) * (real_end_y_ - real_beg_y_); }
  int get_idx(int x, int y) const { return x + Lx_ * y; }

private:
  int gl_Lx_;
  int gl_Ly_;
  int Lx_;
  int Ly_;
  int n_;
  int real_beg_x_;
  int real_beg_y_;
  int real_end_x_;
  int real_end_y_;

  bool x_padded_;
  bool y_padded_;

  MPI_Datatype col_double_type_;
  MPI_Datatype row_double_type_;

  MPI_Datatype gl_block_type_;
  MPI_Datatype my_block_type_;

  int* counts_;
  int* displs_;
};


template <typename BiFunc>
void SquareLatticePadded_2::for_each_pair(BiFunc pair_force) const {
  for (int row = 0; row < real_end_y_; row++) {
    int row_upper = (row + 1) % Ly_;
    size_t row_upper_by_Lx = size_t(row_upper) * Lx_;
    size_t row_by_Lx = size_t(row) * Lx_;
    for (int col = 0; col < real_end_x_; col++) {
      int col_right = (col + 1) % Lx_;
      size_t my_i = col + row_by_Lx;
      size_t i_right = col_right + row_by_Lx;
      size_t i_upper = col + row_upper_by_Lx;
      pair_force(my_i, i_right);
      pair_force(my_i, i_upper);
    }
  }
}


template<typename UniFunc>
void SquareLatticePadded_2::for_each_real_site(UniFunc func) const {
  for (int row = real_beg_y_; row < real_end_y_; row++) {
    size_t row_Lx = row * Lx_;
    for (int col = real_beg_x_; col < real_end_x_; col++) {
      size_t idx = col + row * Lx_;
      func(idx);
    }
  }
}

/*****************************************************************************
*         Functions to initialize fields
*****************************************************************************/

template <typename T, typename TRan>
void get_rand_fields(T* f, size_t n, TRan& myran,
                     double f_a, double f_b) {
  double df = f_b - f_a;
  for (size_t i = 0; i < n; i++) {
    f[i] = f_a + df * myran.doub();
  }
}

template <typename T>
void get_uniform_fields(T* f, size_t n, T f0) {
  for (size_t i = 0; i < n; i++) {
    f[i] = f0;
  }
}

template <typename T, typename TRan>
void get_bimodal_fields(T* f, size_t n, TRan& myran, T f1, T f2) {
  size_t half_n = n / 2;
  for (size_t i = 0; i < half_n; i++) {
    f[i] = f1;
  }
  for (size_t i = half_n; i < n; i++) {
    f[i] = f2;
  }
  shuffle(f, n, myran);
}
#pragma once
#include <vector>
#include "comn.h"

class SquareLattice_2 {
public:
  SquareLattice_2(int Lx, int Ly) : Lx_(Lx), Ly_(Ly){}

  template <typename BiFunc>
  void for_each_pair(BiFunc pair_force) const;

  template <typename TriFunc>
  void for_each_pair_NR(TriFunc NR_force) const;

  int get_Lx() const { return Lx_; }
  int get_Ly() const { return Ly_; }

private:
  int Lx_;
  int Ly_;
};


template<typename BiFunc>
void SquareLattice_2::for_each_pair(BiFunc pair_force) const {
  for (int row = 0; row < Ly_; row++) {
    int row_upper = (row + 1) % Ly_;
    size_t row_upper_by_Lx = size_t(row_upper) * Lx_;
    size_t row_by_Lx = size_t(row) * Lx_;
    for (int col = 0; col < Lx_; col++) {
      int col_right = (col + 1) % Lx_;
      size_t my_i = col + row_by_Lx;
      size_t i_right = col_right + row_by_Lx;
      size_t i_upper = col + row_upper_by_Lx;
      pair_force(my_i, i_right);
      pair_force(my_i, i_upper);
    }
  }
}


template<typename TriFunc>
void SquareLattice_2::for_each_pair_NR(TriFunc NR_force) const {
  for (int row = 0; row < Ly_; row++) {
    int row_upper = (row + 1) % Ly_;
    size_t row_upper_by_Lx = size_t(row_upper) * Lx_;
    //int row_lower = (row - 1 + Ly_) % Ly_;
    //size_t row_lower_by_Lx = size_t(row_lower) * Lx_;
    size_t row_by_Lx = size_t(row) * Lx_;
    for (int col = 0; col < Lx_; col++) {
      int col_right = (col + 1) % Lx_;
      //int col_left = (col - 1 + Lx_) % Lx_;
      size_t my_i = col + row_by_Lx;
      size_t i_right = col_right + row_by_Lx;
      size_t i_upper = col + row_upper_by_Lx;
      //size_t i_left = col_left + row_by_Lx;
      //size_t i_lower = col + row_lower_by_Lx;
      NR_force(my_i, i_right, 0.);
      NR_force(my_i, i_upper, half_PI);
      //NR_force(i_right, my_i, PI);
      //NR_force(i_upper, my_i, -half_PI);
      //NR_force(my_i, i_left,  PI);
      //NR_force(my_i, i_lower, -half_PI);
    }
  }
}
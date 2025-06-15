#pragma once
#include <vector>

class SquareLattice_2 {
public:
  SquareLattice_2(int Lx, int Ly) : Lx_(Lx), Ly_(Ly){}

  template <typename BiFunc>
  void for_each_pair(BiFunc pair_force) const;

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

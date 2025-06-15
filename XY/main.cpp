#include <cstdio>
#include <iostream>
#include "lattice_2.h"
#include "xy_2.h"
#include "rand.h"


int main(int argc, char* argv[]) {
  int Lx = 128;
  int Ly = 128;

  double T = 0.9;
  double h = 0.05;

  unsigned long long seed = 1234;
  Ran myran(seed);

  int snap_dt = 1000;
  int n_step = 100000;

  XY_2<SquareLattice_2> XY_spins(Lx, Ly, T, h);
  //XY_spins.ini_rand(myran);
  XY_spins.ini_ordered();


  double vm, v_angle;
  XY_spins.cal_order_para(vm, v_angle);

  std::cout << "0\t" << vm << "\t" << v_angle << std::endl;


  for (int i = 1; i <= n_step; i++) {
    XY_spins.update(myran);
    if (i % snap_dt == 0) {
      XY_spins.cal_order_para(vm, v_angle);
      std::cout << i * h << "\t" << vm << "\t" << v_angle << std::endl;
    }
  }
}
#pragma once
#include "kuramoto_2.h"

template <typename TLattice>
class NR_Kuramoto_2 : public Kuramoto_2<TLattice> {
public:
  NR_Kuramoto_2(int Lx, int Ly, double T,
    double sigma, double h, double omega0=0.,
    double J0 = 1.0, double J1 = 0.);

  template <typename TRan>
  void update(TRan& myran);

protected:
  double J1_;
};

template<typename TLattice>
NR_Kuramoto_2<TLattice>::NR_Kuramoto_2(int Lx, int Ly, double T,
                                       double sigma, double h,
                                       double omega0, double J0, double J1)
  : Kuramoto_2<TLattice>(Lx, Ly, T, sigma, h, omega0, J0), J1_(J1) {}


template<typename TLattice>
template<typename TRan>
void NR_Kuramoto_2<TLattice>::update(TRan& myran) {
  // align with nearest neighbors
  auto NR_force = [this](size_t i, size_t j, double angle_ij) {
    XY_2<TLattice>::align(i, j, angle_ij, XY_2<TLattice>::J0_, J1_);
    };
  XY_2<TLattice>::lat_.for_each_pair_NR(NR_force);

  // update theta
  for (size_t i = 0; i < XY_2<TLattice>::N_; i++) {
    XY_2<TLattice>::update_theta(i, myran, Kuramoto_2<TLattice>::omega_[i]);
    XY_2<TLattice>::tangle_theta(i);
    XY_2<TLattice>::tau_[i] = 0.;
  }
}

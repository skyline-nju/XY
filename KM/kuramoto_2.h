#pragma once
#include "xy_2.h"

template <typename TLattice>
class Kuramoto_2 : public XY_2<TLattice> {
public:
  Kuramoto_2(int Lx, int Ly, double T, double sigma, double h, double J = 1.0);

  ~Kuramoto_2();

  template <typename TSnap, typename TRan>
  void ini(const std::string& ini_mode, const TSnap& snap, TRan& myran, double theta0=0.);

  template <typename TRan>
  void ini_omega(TRan& myran, const std::string& mode = "bimodal");

  template <typename TRan>
  void update(TRan& myran);

  template <typename TRan>
  void update_theta(size_t i, TRan& myran);

  template <typename TFloat>
  void get_pos(TFloat* pos) const;

protected:
  double sigma_;
  double h_;
  double* omega_;
};


template<typename TLattice>
Kuramoto_2<TLattice>::Kuramoto_2(int Lx, int Ly,
                                 double T, double sigma,
                                 double h, double J)
  : XY_2<TLattice>(Lx, Ly, T, h, J), sigma_(sigma), h_(h) {
  omega_ = new double[XY_2<TLattice>::N_] {};
}

template<typename TLattice>
Kuramoto_2<TLattice>::~Kuramoto_2() {
  delete[] omega_;
}

template<typename TLattice>
template<typename TSnap, typename TRan>
void Kuramoto_2<TLattice>::ini(const std::string& ini_mode, const TSnap& snap,
                               TRan& myran, double theta0){
  if (ini_mode == "rand") {
    for (size_t i = 0; i < XY_2<TLattice>::N_; i++) {
      XY_2<TLattice>::theta_[i] = (myran.doub() - 0.5) * 2 * PI;
    }
    ini_omega(myran);
  } else if (ini_mode == "ordered") {
    for (size_t i = 0; i < XY_2<TLattice>::N_; i++) {
      XY_2<TLattice>::theta_[i] = theta0;
    }
    ini_omega(myran);
  } else if (ini_mode == "resume") {
    double* pos = new double[XY_2<TLattice>::N_ * 3];
    snap.read_last_frame(pos, XY_2<TLattice>::theta_, XY_2<TLattice>::N_);
    for (size_t i = 0; i < XY_2<TLattice>::N_; i++) {
      omega_[i] = pos[i * 3 + 2];
    }
    delete[] pos;
  } else {
    std::cout << "Error, ini_mode must be one of rand, ordered and resume" << std::endl;
    exit(1);
  }
}


template<typename TLattice>
template<typename TRan>
void Kuramoto_2<TLattice>::ini_omega(TRan& myran, const std::string& mode) {
  if (mode == "bimodal") {
    size_t n_half = XY_2<TLattice>::N_ / 2;
    for (size_t i = 0; i < n_half; i++) {
      omega_[i] = sigma_;
    }
    for (size_t i = n_half; i < XY_2<TLattice>::N_; i++) {
      omega_[i] = -sigma_;
    }
    shuffle(omega_, XY_2<TLattice>::N_, myran);
  }
}


template<typename TLattice>
template<typename TRan>
void Kuramoto_2<TLattice>::update(TRan& myran) {
  // align with nearest neighbors
  auto pair_force = [this](size_t i, size_t j) {
    XY_2<TLattice>::align(i, j);
    };
  XY_2<TLattice>::lat_.for_each_pair(pair_force);

  // update theta
  for (size_t i = 0; i < XY_2<TLattice>::N_; i++) {
    update_theta(i, myran);
    XY_2<TLattice>::tangle_theta(i);
    XY_2<TLattice>::tau_[i] = 0.;
  }
}

template<typename TLattice>
template<typename TRan>
void Kuramoto_2<TLattice>::update_theta(size_t i, TRan& myran) {
  double dtheta = omega_[i] * h_
    + XY_2<TLattice>::J_h_ * XY_2<TLattice>::tau_[i]
    + XY_2<TLattice>::sqrt_24_T_h_ * (myran.doub() - 0.5);
  XY_2<TLattice>::theta_[i] += dtheta;
}

template<typename TLattice>
template<typename TFloat>
void Kuramoto_2<TLattice>::get_pos(TFloat* pos) const {
  int Lx = XY_2<TLattice>::lat_.get_Lx();
  int Ly = XY_2<TLattice>::lat_.get_Ly();

  for (int row = 0; row < Ly; row++) {
    TFloat y = row + 0.5 - Ly / 2.;
    for (int col = 0; col < Lx; col++) {
      size_t idx = col + size_t(row) * Lx;
      TFloat x = col + 0.5 - Lx / 2.;
      pos[3 * idx] = x;
      pos[3 * idx + 1] = y;
      pos[3 * idx + 2] = omega_[idx];
    }
  }
}

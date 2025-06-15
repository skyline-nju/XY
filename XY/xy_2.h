#pragma once
#include "comn.h"

template <typename TLattice>
class XY_2 {
public:
  XY_2(int Lx, int Ly, double T, double h, double J=1.0);
  ~XY_2();

  template <typename TRan>
  void ini_rand(TRan& myran);

  void ini_ordered(double theta0 = 0.);

  void align(size_t i, size_t j);

  template <typename TRan>
  void update(TRan& myran);

  void cal_order_para(double& v_module, double& v_angle) const;

private:
  TLattice lat_;
  size_t N_;
  double J_h_;           // J * h
  double sqrt_24_T_h_;   // sqrt(24 * T * h)

  double* theta_;
  double* tau_;
};


template<typename TLattice>
XY_2<TLattice>::XY_2(int Lx, int Ly, double T, double h, double J): lat_(Lx, Ly), N_(size_t(Lx) * size_t(Ly)) {
  J_h_ = h * J;
  sqrt_24_T_h_ = std::sqrt(24. * T * h);
  theta_ = new double[N_] {};
  tau_ = new double[N_] {};
}

template<typename TLattice>
XY_2<TLattice>::~XY_2() {
  delete[] theta_;
  delete[] tau_;
}

template<typename TLattice>
template<typename TRan>
void XY_2<TLattice>::ini_rand(TRan& myran) {
  for (size_t i = 0; i < N_; i++) {
    theta_[i] = myran.doub() * 2 * PI;
  }
}

template<typename TLattice>
void XY_2<TLattice>::ini_ordered(double theta0) {
  for (size_t i = 0; i < N_; i++) {
    theta_[i] = theta0;
  }
}


template<typename TLattice>
void XY_2<TLattice>::align(size_t i, size_t j) {
  double tmp = std::sin(theta_[j] - theta_[i]);
  tau_[i] += tmp;
  tau_[j] -= tmp;
}

template<typename TLattice>
template<typename TRan>
void XY_2<TLattice>::update(TRan& myran) {
  // align with nearest neighbors
  auto pair_force = [this](size_t i, size_t j) {
    align(i, j);
  };
  lat_.for_each_pair(pair_force);

  // update theta
  for (size_t i = 0; i < N_; i++) {
    theta_[i] += (J_h_ * tau_[i] + sqrt_24_T_h_ * (myran.doub() - 0.5));
    tau_[i] = 0.;
  }
}

template<typename TLattice>
void XY_2<TLattice>::cal_order_para(double& v_module, double& v_angle) const {
  double vx = 0;
  double vy = 0;
  for (size_t i = 0; i < N_; i++) {
    vx += std::cos(theta_[i]);
    vy += std::sin(theta_[i]);
  }

  vx /= N_;
  vy /= N_;
  v_module = std::sqrt(vx * vx + vy * vy);
  v_angle = std::atan2(vy, vx);
}
#pragma once
#include "comn.h"

template <typename TLattice>
class XY_2 {
public:
  XY_2(int Lx, int Ly, double T, double h, double J0=1.0);
  ~XY_2();

  template <typename TSnap, typename TRan>
  void ini(const std::string& ini_mode, const TSnap& snap, TRan& myran, double theta0=0.);

  void align(size_t i, size_t j);

  void align(size_t i, size_t j, double psi, double J0, double J1);

  template <typename TRan>
  void update(TRan& myran);

  void cal_order_para(double& v_module, double& v_angle) const;

  template <typename TFloat>
  void get_theta(TFloat* theta) const;

  template <typename TFloat>
  void get_pos(TFloat* pos) const;

  size_t get_n() const { return N_; }

  template <typename TRan>
  void update_theta(size_t i, TRan& myran, double intrinsic_freq = 0.);

  void tangle_theta(size_t i);

protected:
  TLattice lat_;
  size_t N_;
  double h_;            
  double J0_;           // J0
  double sqrt_24_T_h_;   // sqrt(24 * T * h)

  double* theta_;
  double* tau_;
};


template<typename TLattice>
XY_2<TLattice>::XY_2(int Lx, int Ly, double T, double h, double J0): lat_(Lx, Ly), N_(size_t(Lx) * size_t(Ly)), h_(h), J0_(J0) {
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
template<typename TSnap, typename TRan>
void XY_2<TLattice>::ini(const std::string& ini_mode, const TSnap& snap, TRan& myran, double theta0) {
  if (ini_mode == "rand") {
    for (size_t i = 0; i < N_; i++) {
      theta_[i] = (myran.doub() - 0.5) * 2 * PI;
    }
  } else if (ini_mode == "ordered") {
    for (size_t i = 0; i < N_; i++) {
      theta_[i] = theta0;
    }
  } else if (ini_mode == "resume") {
    double* pos = nullptr;
    snap.read_last_frame(pos, theta_, N_);
  } else {
    std::cout << "Error, ini_mode must be one of rand, ordered and resume" << std::endl;
    exit(1);
  }
}

template<typename TLattice>
void XY_2<TLattice>::align(size_t i, size_t j) {
  double tmp = J0_ * std::sin(theta_[j] - theta_[i]);
  tau_[i] += tmp;
  tau_[j] -= tmp;
}

// calculate the torque felt by i resulted by j
// tau_i = J * sin(theta_j - theta_i)
// J = J0 + J1 * cos(psi)
// psi = theta_i - arg(r_j - r_i)
template<typename TLattice>
void XY_2<TLattice>::align(size_t i, size_t j, double angle_ij, double J0, double J1){
  double psi_ij = theta_[i] - angle_ij;
  double psi_ji = theta_[j] - (angle_ij + PI);
  double Ji = J0 + J1 * std::cos(psi_ij);
  double Jj = J0 + J1 * std::cos(psi_ji);
  double sin_dtheta = std::sin(theta_[j] - theta_[i]);
  tau_[i] += Ji * sin_dtheta;
  tau_[j] -= Jj * sin_dtheta;
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
    update_theta(i, myran);
    tangle_theta(i);
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

template<typename TLattice>
template<typename TFloat>
void XY_2<TLattice>::get_theta(TFloat* theta) const {
  for (size_t i = 0; i < N_; i++) {
    theta[i] = theta_[i];
  }
}

template<typename TLattice>
template<typename TFloat>
void XY_2<TLattice>::get_pos(TFloat* pos) const {
  int Lx = lat_.get_Lx();
  int Ly = lat_.get_Ly();

  for (int row = 0; row < Ly; row++) {
    TFloat y = row + 0.5 - Ly / 2.;
    for (int col = 0; col < Lx; col++) {
      size_t idx = col + size_t(row) * Lx;
      TFloat x = col + 0.5 - Lx / 2.;
      pos[3 * idx] = x;
      pos[3 * idx + 1] = y;
      pos[3 * idx + 2] = 0.;
    }
  }
}

template<typename TLattice>
template<typename TRan>
void XY_2<TLattice>::update_theta(size_t i, TRan& myran, double intrinsic_freq) {
  theta_[i] += (h_ * (tau_[i] + intrinsic_freq) + sqrt_24_T_h_ * (myran.doub() - 0.5));
}

template<typename TLattice>
void XY_2<TLattice>::tangle_theta(size_t i) {
  if (theta_[i] >= PI) {
    theta_[i] -= 2 * PI;
  } else if (theta_[i] < -PI) {
    theta_[i] += 2 * PI;
  }
}
#pragma once
#include "xy_2.h"

template <typename TLattice>
class Kuramoto_2 : public XY_2<TLattice> {
public:
  Kuramoto_2(int gl_Lx, int gl_Ly, double sigma,
             int proc_nx, int proc_ny, MPI_Comm group_comm,
             double T, double h, double J = 1.0)
    : XY_2<TLattice>(gl_Lx, gl_Ly, proc_nx, proc_ny, group_comm, T, h, J),
      sigma_(sigma), omega_(new double[this->get_n()] {}) {};

  ~Kuramoto_2();

  template <typename TSnap, typename TRan>
  void ini(const std::string& ini_mode,
           const TSnap& snap, TRan& myran,
           double theta0 = 0.);

  template <typename TRan>
  void update(TRan& myran) {
    XY_2<TLattice>::update_tau();
    XY_2<TLattice>::update_theta(myran, omega_);
  }

  template <typename TFloat>
  void get_pos(TFloat* pos_gl) const {
    XY_2<TLattice>::get_pos(pos_gl, omega_);
  }

protected:
  double sigma_;
  double* omega_;
};


template <typename TLattice>
Kuramoto_2<TLattice>::~Kuramoto_2() {
  delete[] omega_;
}

template<typename TLattice>
template<typename TSnap, typename TRan>
void Kuramoto_2<TLattice>::ini(const std::string& ini_mode,
                               const TSnap& snap,
                               TRan& myran, double theta0) {
  double* gl_theta = nullptr;
  double* gl_omega = nullptr;
  size_t n_gl = XY_2<TLattice>::get_gl_n();
  std::cout << "N=" << n_gl << std::endl;

  if (XY_2<TLattice>::is_root()) {
    gl_theta = new double[n_gl] {};
    gl_omega = new double[n_gl] {};

    if (ini_mode == "resume") {
      snap.read_last_frame(n_gl, gl_theta, gl_omega);
    } else if (ini_mode == "rand") {
      get_rand_fields(gl_theta, n_gl, myran, -PI, PI);
      get_bimodal_fields(gl_omega, n_gl, myran, -sigma_, sigma_);
    } else if (ini_mode == "ordered") {
      get_uniform_fields(gl_theta, n_gl, theta0);
      get_bimodal_fields(gl_omega, n_gl, myran, -sigma_, sigma_);
    } else {
      std::cout << "Error, ini_mode must be one of rand, ordered and resume" << std::endl;
      exit(1);
    }
  }
  XY_2<TLattice>::lat_.scatter_fields(gl_omega, omega_);
  XY_2<TLattice>::lat_.scatter_fields(gl_theta, this->theta_);

  delete[] gl_theta;
  delete[] gl_omega;
}

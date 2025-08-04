#pragma once
#include "comn.h"
#include "lattice_2.h"
#include <mpi.h>


template <typename TLattice>
class XY_2 {
public:
  XY_2(int gl_Lx, int gl_Ly,
       int proc_nx, int proc_ny,
       MPI_Comm group_comm,
       double T, double h, double J = 1.0);
  ~XY_2();

  template <typename TSnap, typename TRan>
  void ini(const std::string& ini_mode, const TSnap& snap,
    TRan& myran, double theta0 = 0.);

  void align(size_t i, size_t j);

  void update_tau();

  template<typename TRan>
  void update_theta(TRan& myran, const double* omega = nullptr);

  template <typename TRan>
  void update(TRan& myran) {
    update_tau();
    update_theta(myran);
  }

  void cal_order_para(double& v_module, double& v_angle) const;

  void get_theta(double* theta_gl) const;

  void get_theta(float* theta_gl) const;

  template <typename TFloat>
  void get_pos(TFloat* pos_gl, const double* omega) const;

  template <typename TFloat>
  void get_pos(TFloat* pos_gl) const {
    get_pos(pos_gl, nullptr);
  }

  size_t get_n() const { return N_; }
  size_t get_gl_n() const { return gl_N_; }

  bool is_root() const { return lat_.is_root(); }

protected:
  TLattice lat_;
  size_t N_;
  size_t gl_N_;
  double h_;
  double J_h_;           // J * h
  double sqrt_24_T_h_;   // sqrt(24 * T * h)

  double* theta_;
  double* tau_;
};


template<typename TLattice>
XY_2<TLattice>::XY_2(int gl_Lx, int gl_Ly,
                     int proc_nx, int proc_ny,
                     MPI_Comm group_comm,
                     double T, double h, double J)
  : lat_(gl_Lx, gl_Ly, proc_nx, proc_ny, group_comm), N_(lat_.get_n()), h_(h) {
  gl_N_ = size_t(gl_Lx) * size_t(gl_Ly);
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
template<typename TSnap, typename TRan>
void XY_2<TLattice>::ini(const std::string& ini_mode,
                         const TSnap& snap,
                         TRan& myran, double theta0) {
  double* gl_theta = nullptr;

  if (is_root()) {
    gl_theta = new double[gl_N_];

    if (ini_mode == "resume") {
      snap.read_last_frame(gl_N_, gl_theta);
    } else if (ini_mode == "rand") {
      get_rand_fields(gl_theta, gl_N_, myran, -PI, PI);
    } else if (ini_mode == "ordered") {
      get_uniform_fields(gl_theta, gl_N_, theta0);
    } else {
      std::cout << "Error, ini_mode must be one of rand, ordered and resume" << std::endl;
      exit(1);
    }
  }

  lat_.scatter_fields(gl_theta, theta_);
  delete[] gl_theta;
}

template<typename TLattice>
void XY_2<TLattice>::align(size_t i, size_t j) {
  double tmp = std::sin(theta_[j] - theta_[i]);
  tau_[i] += tmp;
  tau_[j] -= tmp;
}

template<typename TLattice>
void XY_2<TLattice>::update_tau() {
  // copy theta in padded regions to nearest proc
  lat_.comm_before_cal_force(theta_);

  // align with nearest neighbors
  auto pair_force = [this](size_t i, size_t j) {
    align(i, j);
    };
  lat_.for_each_pair(pair_force);
}

template<typename TLattice>
template<typename TRan>
void XY_2<TLattice>::update_theta(TRan& myran, const double* omega) {
  auto updating = [this, &myran, omega](size_t i) {
    this->theta_[i] += this->J_h_ * this->tau_[i] + this->sqrt_24_T_h_ * (myran.doub() - 0.5);
    if (omega) {
      this->theta_[i] += omega[i] * this->h_;
    }
    tangle_theta(this->theta_[i]);
    this->tau_[i] = 0.;
  };
  lat_.for_each_real_site(updating);

  auto set_zero = [this](size_t i) {
    this->tau_[i] = 0.;
  };
  lat_.for_each_padded_site(set_zero);
}


template<typename TLattice>
void XY_2<TLattice>::cal_order_para(double& v_module, double& v_angle) const {
  double my_v[2]{};
  int count = 0;

  auto sum_v = [this, &my_v, &count](size_t idx) {
    my_v[0] += std::cos(theta_[idx]);
    my_v[1] += std::sin(theta_[idx]);
    count++;
    };

  lat_.for_each_real_site(sum_v);
  if (count != lat_.get_real_n()) {
    std::cout << "Error, found " << count << " sites" << std::endl;
    exit(1);
  }

  int root = 0;
  double v_gl[2]{};
  int n_gl = 0;
  MPI_Comm comm = lat_.get_communicator();
  MPI_Reduce(my_v, v_gl, 2, MPI_DOUBLE, MPI_SUM, root, comm);
  MPI_Reduce(&count, &n_gl, 1, MPI_INT, MPI_SUM, root, comm);

  if (lat_.is_root()) {
    v_gl[0] /= n_gl;
    v_gl[1] /= n_gl;

    v_module = std::sqrt(v_gl[0] * v_gl[0] + v_gl[1] * v_gl[1]);
    v_angle = std::atan2(v_gl[1], v_gl[0]);
  }
}

template<typename TLattice>
void XY_2<TLattice>::get_theta(double* theta_gl) const {
  lat_.gather_fields(theta_gl, theta_);
}

template<typename TLattice>
void XY_2<TLattice>::get_theta(float* theta_gl) const {
  double* theta2 = nullptr;
  if (is_root()) {
    theta2 = new double[gl_N_] {};
  }
  get_theta(theta2);

  if (is_root()) {
    for (size_t i = 0; i < gl_N_; i++) {
      theta_gl[i] = theta2[i];
    }
    delete[] theta2;
  }
}

template<typename TLattice>
template<typename TFloat>
void XY_2<TLattice>::get_pos(TFloat* pos_gl, const double* omega) const {
  double* omega_gl = nullptr;
  if (is_root()) {
    omega_gl = new double[gl_N_] {};
  }
  if (omega) {
    lat_.gather_fields(omega_gl, omega);
  }

  if (is_root()) {
    int nrows = lat_.get_gl_Ly();
    int ncols = lat_.get_gl_Lx();
    for (int row = 0; row < nrows; row++) {
      TFloat y = row + 0.5 - nrows / 2.;
      for (int col = 0; col < ncols; col++) {
        size_t idx = col + size_t(row) * ncols;
        TFloat x = col + 0.5 - ncols / 2.;
        pos_gl[3 * idx] = x;
        pos_gl[3 * idx + 1] = y;
        pos_gl[3 * idx + 2] = omega_gl[idx];
      }
    }
    delete[] omega_gl;
  }
}
#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <mpi.h>
#include "comn.h"
#include "gsd.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

namespace io {

  /**
    * @brief Basic class for exporting data.
    *
    * Define the timming to dump data.
    */
  class ExporterBase {
  public:
    ExporterBase(int n_step, int sep, int start, MPI_Comm group_comm);

    bool need_export(const int i_step);

    bool need_export_log_scale(const int i_step);

    void set_log_scale_frames(double h, double log_sep = 0.1);

    bool is_root() const { return my_rank_ == 0; }

  protected:
    int n_step_;    // total steps to run
    int sep_;
    int start_ = 0; // The first step 
    int tot_proc_ = 1;
    std::vector<int> frames_;
    int cur_frame_ = 0;
    double log_sep_ = -1;
    MPI_Comm comm_;
    int my_rank_;
  };

  /**
    * @brief Exporter to output log
    *
    * Output the parameters after the initialization.
    * Output the beginning and endding time of the simulation.
    * Record time every certain time steps.
    */
  class LogExporter : public ExporterBase {
  public:
    LogExporter(const std::string& outfile,
      int n_step, int sep, int start,
      int np_gl, MPI_Comm group_comm);

    ~LogExporter();

    void record(int i_step);

    std::ofstream fout;
  private:
    std::chrono::time_point<std::chrono::system_clock> t_start_;
    int n_par_;
    int step_count_ = 0;
  };


  class OrderParaExporter : public ExporterBase {
  public:
    OrderParaExporter(const std::string& outfile, int start, int n_step, int sep,
      int flush_sep, MPI_Comm group_comm);

    ~OrderParaExporter() { if (my_rank_ == 0) fout_.close(); };

    template <typename TSpins>
    void dump(int i_step, const TSpins& spins);

  private:
    std::ofstream fout_;
    int flush_sep_;
  };


  template <typename TSpins>
  void OrderParaExporter::dump(int i_step, const TSpins& spins) {
    if (need_export(i_step)) {
      double p, theta;
      spins.cal_order_para(p, theta);
      if (is_root()) {
        fout_ << i_step << "\t" << std::setprecision(8) << p << "\t" << theta;
        if (i_step % flush_sep_ == 0) {
          fout_ << std::endl;
        } else {
          fout_ << "\n";
        }
      }
    }
  }

  class Snap_GSD_2 : public ExporterBase {
  public:
    Snap_GSD_2(const std::string& filename,
               int n_step, int sep, int& start,
               double h, double log_sep,
               int Lx, int Ly,
               const std::string& open_flag,
               MPI_Comm group_comm);

    ~Snap_GSD_2();

    uint64_t get_time_step();

    int reset_start_time_step();

    template <typename TSpins>
    void dump(int i_step, const TSpins& spins, bool flag_out_pos=false);

    template <typename TFloat>
    void read(int i_frame, int n, TFloat* theta_out, TFloat* omega_out=nullptr) const;

    template <typename TFloat>
    void read_last_frame(int n, TFloat* theta_out, TFloat* omega_out=nullptr) const;

  private:
    gsd_handle* handle_ = nullptr;
  };


  template<typename TSpins>
  void Snap_GSD_2::dump(int i_step, const TSpins& spins, bool flag_out_pos) {
    if (need_export(i_step)) {
      int nframes;
      if (is_root()) {
        nframes= gsd_get_nframes(handle_);
      }
      MPI_Bcast(&nframes, 1, MPI_INT, 0, comm_);

      size_t n_gl = spins.get_gl_n();
      float* theta = nullptr;
      if (is_root()) {
        theta = new float[n_gl];
      }
      spins.get_theta(theta);

      float* pos = nullptr;
      if (nframes == 0 || flag_out_pos) {
        if (is_root()) {
          pos = new float[n_gl * 3];
        }
        spins.get_pos(pos);
      }

      uint64_t step = start_ + i_step;
      if (is_root()) {
        gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
        gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n_gl);
        gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n_gl, 1, 0, theta);
        if (nframes == 0 || flag_out_pos) {
          gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n_gl, 3, 0, pos);
        }
        gsd_end_frame(handle_);
        std::cout << "dump frame " << nframes << " at time step " << step << std::endl;
        delete[] theta;
        delete[] pos;
      }
      MPI_Barrier(comm_);

    }
  }


  template<typename TFloat>
  void io::Snap_GSD_2::read(int i_frame, int n, TFloat* theta_out, TFloat* omega) const {
    if (is_root()) {
      uint32_t n_par;
      const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
      gsd_read_chunk(handle_, &n_par, chunk);
      if (n == n_par) {
        std::cout << "frame " << i_frame << ": find " << n_par << " spins" << std::endl;
      } else {
        std::cout << "Error, " << n_par << " spins were found, while " << n << " spins are expected" << std::endl;
        exit(1);
      }

      float* theta = new float[n_par];
      chunk = gsd_find_chunk(handle_, i_frame, "particles/charge");
      gsd_read_chunk(handle_, theta, chunk);

      for (size_t i = 0; i < n_par; i++) {
        theta_out[i] = theta[i];
      }
      delete[] theta;

      if (omega) {
        float* pos = new float[n_par * 3];
        chunk = gsd_find_chunk(handle_, 0, "particles/position");
        gsd_read_chunk(handle_, pos, chunk);

        for (size_t i = 0; i < n_par; i++) {
          omega[i] = pos[3 * i + 2];
        }
        delete[] pos;
      }
    }
  }


  template<typename TFloat>
  void io::Snap_GSD_2::read_last_frame(int n, TFloat* theta_out, TFloat* omega) const {
    if (is_root()) {
      int nframes = gsd_get_nframes(handle_);
      if (nframes < 1) {
        std::cout << "Error, nframes=" << nframes << std::endl;
        exit(1);
      } else {
        read(nframes - 1, n, theta_out, omega);
      }
    }
  }

}

#pragma once
#include <vector>
#include <chrono>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
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
  ExporterBase(int n_step, int sep, int start);

  bool need_export(const int i_step) {
    if (log_sep_ < 0) {
      return i_step % sep_ == 0;
    } else {
      return need_export_log_scale(i_step);
    }
  }

  bool need_export_log_scale(const int i_step) {
    bool res = false;
    if (i_step + start_ == frames_[cur_frame_]) {
      cur_frame_++;
      res = true;
    }
    return res;
  }


  void set_log_scale_frames(double h, double log_sep = 0.1);

protected:
  int n_step_;    // total steps to run
  int sep_;
  int start_ = 0; // The first step 
  int my_rank_ = 0;
  int tot_proc_ = 1;
  std::vector<int> frames_;
  int cur_frame_ = 0;
  double log_sep_ = -1;
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
    int np_gl);

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
  OrderParaExporter(const std::string& outfile, int start, int n_step, int sep, int flush_sep);

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
    fout_ << i_step << "\t" << std::setprecision(8) << p << "\t" << theta;
    if (i_step % flush_sep_ == 0) {
      fout_ << std::endl;
    } else {
      fout_ << "\n";
    }
  }
}

class Snap_GSD_2 : public ExporterBase {
public:
  Snap_GSD_2(const std::string& filename,
    int n_step, int sep, int& start,
    double h, double log_sep,
    int Lx, int Ly,
    const std::string& open_flag);

  ~Snap_GSD_2() { gsd_close(handle_); delete handle_; };


  uint64_t get_time_step();

  int reset_start_time_step();

  template <typename TSpins>
  void dump(int i_step, const TSpins& spins);

  template <typename TFloat>
  void read(int i_frame, TFloat* pos_out, TFloat* theta_out, int n) const;

  template <typename TFloat>
  void read_last_frame(TFloat* pos_out, TFloat* theta_out, int n) const;

private:
  gsd_handle* handle_ = nullptr;
};


template<typename TSpins>
void Snap_GSD_2::dump(int i_step, const TSpins& spins) {
  if (need_export(i_step)) {
    int nframes = gsd_get_nframes(handle_);

    size_t n = spins.get_n();
    float* theta = new float[n];
    spins.get_theta(theta);

    uint64_t step = start_ + i_step;
    gsd_write_chunk(handle_, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &step);
    gsd_write_chunk(handle_, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n);
    gsd_write_chunk(handle_, "particles/charge", GSD_TYPE_FLOAT, n, 1, 0, theta);

    float* pos = nullptr;
    if (nframes == 0) {
      pos = new float[n * 3];
      spins.get_pos(pos);
      gsd_write_chunk(handle_, "particles/position", GSD_TYPE_FLOAT, n, 3, 0, pos);
    }

    gsd_end_frame(handle_);

    std::cout << "dump frame " << nframes << " at time step " << step << std::endl;
    delete[] theta;
    delete[] pos;
  }
}


template<typename TFloat>
void io::Snap_GSD_2::read(int i_frame, TFloat* pos_out, TFloat* theta_out, int n) const {
  uint32_t n_par;
  const gsd_index_entry* chunk = gsd_find_chunk(handle_, i_frame, "particles/N");
  gsd_read_chunk(handle_, &n_par, chunk);
  if (n == n_par) {
    std::cout << "frame " << i_frame << ": find " << n_par << " spins" << std::endl;
  } else {
    std::cout << "Error, " << n_par << " spins were found, while " << n << " spins are expected" << std::endl;
    exit(1);
  }
  float* pos = new float[n_par * 3];
  chunk = gsd_find_chunk(handle_, 0, "particles/position");
  gsd_read_chunk(handle_, pos, chunk);

  float* theta = new float[n_par];
  chunk = gsd_find_chunk(handle_, i_frame, "particles/charge");
  gsd_read_chunk(handle_, theta, chunk);

  if (pos_out != nullptr) {
    for (size_t i = 0; i < n_par; i++) {
      pos_out[i * 3] = pos[i * 3];
      pos_out[i * 3 + 1] = pos[i * 3 + 1];
      pos_out[i * 3 + 2] = pos[i * 3 + 2];
    }
  }

  for (size_t i = 0; i < n_par; i++) {
    theta_out[i] = theta[i];
  }
  delete[] pos;
  delete[] theta;
}


template<typename TFloat>
void io::Snap_GSD_2::read_last_frame(TFloat* pos_out, TFloat* theta_out, int n) const {
  int nframes = gsd_get_nframes(handle_);
  if (nframes < 1) {
    std::cout << "Error, nframes=" << nframes << std::endl;
    exit(1);
  } else {
    read(nframes - 1, pos_out, theta_out, n);
  }
}

}

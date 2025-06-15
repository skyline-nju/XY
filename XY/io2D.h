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
void io::OrderParaExporter::dump(int i_step, const TSpins& spins) {
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

}

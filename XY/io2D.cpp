#include "io2D.h"

io::ExporterBase::ExporterBase(int n_step, int sep, int start)
  : n_step_(n_step), sep_(sep), start_(start) {}

void io::ExporterBase::set_log_scale_frames(double h, double log_sep) {
  log_sep_ = log_sep;
  int frame_size = int(20 / log_sep);
  frames_.reserve(frame_size);
  std::vector<double> log10_t;
  log10_t.reserve(frame_size);
  log10_t.push_back(1);
  for (int i = 1; i < frame_size; i++) {
    log10_t.push_back(log10_t[i-1]+log_sep);
  }
  for (int i = 0; i < frame_size; i++) {
    double t = pow(10, log10_t[i]);
    int i_step = round(t/h);
    frames_.push_back(i_step);
  }
  int cur_step = start_ + 1;
  if (cur_step <= frames_[0]) {
    cur_frame_ = 0;
  } else {
    for (int i = 1; i < frames_.size(); i++) {
      if (cur_step <= frames_[i] && cur_step > frames_[i-1]) {
        cur_frame_ = i;
        break;
      }
    }
  }
}

io::LogExporter::LogExporter(const std::string& outfile, 
                             int n_step, int sep, int start,
                             int np_gl)
  : ExporterBase(n_step, sep, start), n_par_(np_gl){
  if (my_rank_ == 0) {
    fout.open(outfile);
    t_start_ = std::chrono::system_clock::now();
    auto start_time = std::chrono::system_clock::to_time_t(t_start_);
    char str[100];
    tm now_time;
#ifdef _MSC_VER
    localtime_s(&now_time, &start_time);
#else
    localtime_r(&start_time, &now_time);
#endif
    std::strftime(str, 100, "%c", &now_time);
    fout << "Started simulation at " << str << "\n";
    std::cout << "log file: " << outfile << std::endl;
  }
}


io::LogExporter::~LogExporter() {
  if (my_rank_ == 0) {
    const auto t_now = std::chrono::system_clock::now();
    auto end_time = std::chrono::system_clock::to_time_t(t_now);
    char str[100];
    tm now_time;
#ifdef _MSC_VER
    localtime_s(&now_time, &end_time);
#else
    localtime_r(&end_time, &now_time);
#endif
    std::strftime(str, 100, "%c", &now_time);
    fout << "Finished simulation at " << str << "\n";
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    fout << "speed=" << std::scientific << step_count_ * double(n_par_) / elapsed_seconds.count()
      << " particle time step per second per core\n";
    fout.close();
  }
}

void io::LogExporter::record(int i_step) {
  if (need_export(i_step + start_) && my_rank_ == 0) {
    const auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = t_now - t_start_;
    const auto dt = elapsed_seconds.count();
    const auto hour = floor(dt / 3600);
    const auto min = floor((dt - hour * 3600) / 60);
    const int sec = dt - hour * 3600 - min * 60;
    fout << i_step + start_ << "\t" << hour << ":" << min << ":" << sec << std::endl;
  }
  step_count_++;
}


io::OrderParaExporter::OrderParaExporter(const std::string& outfile, int start, int n_step, int sep, int flush_sep)
  : ExporterBase(n_step, sep, start), flush_sep_(flush_sep){
  if (my_rank_ == 0) {
    fout_.open(outfile);
    std::cout << "dump order para to " << outfile << std::endl;
  }
}



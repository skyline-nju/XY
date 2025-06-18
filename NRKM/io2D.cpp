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


io::Snap_GSD_2::Snap_GSD_2(const std::string& filename, int n_step, int sep, int& start,
                           double h, double log_sep, int Lx, int Ly, const std::string& open_flag)
  : ExporterBase(n_step, sep, start) {
  unsigned int version = gsd_make_version(1, 4);
  handle_ = new gsd_handle;
  if (open_flag != "resume") {
    int flag = gsd_create(filename.c_str(), "cpp", "hoomd", version);
    if (flag != 0) {
      std::cout << "Error when create " << filename << "; state=" << flag << std::endl;
      exit(1);
    }
    flag = gsd_open(handle_, filename.c_str(), GSD_OPEN_READWRITE);
    if (flag != 0) {
      std::cout << "Error when open " << filename << "; state=" << flag << std::endl;
      exit(1);
    }

    float box[6] = { Lx, Ly, 1, 0, 0, 0 };
    gsd_write_chunk(handle_, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, box);

  } else {
    int flag = gsd_open(handle_, filename.c_str(), GSD_OPEN_READWRITE);
    if (flag != 0) {
      std::cout << "Error when open " << filename << "; state=" << flag << std::endl;
      exit(1);
    } else {
      std::cout << "open " << filename << std::endl;
    }
  }

  start = reset_start_time_step();
  if (log_sep > 0) {
    set_log_scale_frames(h, log_sep);
  }
}


uint64_t io::Snap_GSD_2::get_time_step() {
  uint64_t step;
  size_t n_frame = gsd_get_nframes(handle_);
  if (n_frame == 0) {
    step = sep_;
  } else {
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, n_frame - 1, "configuration/step");
    if (chunk) {
      gsd_read_chunk(handle_, &step, chunk);
      step += sep_;
    } else {
      step = sep_;
    }
  }
  return step;
}

int io::Snap_GSD_2::reset_start_time_step() {
  size_t n_frame = gsd_get_nframes(handle_);
  if (n_frame == 0) {
    start_ = 0;
  } else {
    const gsd_index_entry* chunk = gsd_find_chunk(handle_, n_frame - 1, "configuration/step");
    if (chunk) {
      uint64_t last_step;
      gsd_read_chunk(handle_, &last_step, chunk);
      start_ = last_step;
    } else {
      std::cout << "Warning, failed to read the time step of the last frame" << std::endl;
      start_ = 0;
    }
  }
  return start_;
}

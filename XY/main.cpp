#include "lattice_2.h"
#include "xy_2.h"
#include "rand.h"
#include "io2D.h"


int main(int argc, char* argv[]) {
  // set parameters
  int Lx = 128;
  int Ly = 128;

#ifdef _MSC_VER
  double T = 0.2;
  double h = 0.1;
  unsigned long long seed = 1234;
  int snap_dt = 1000;
  int n_step = 30000;
  std::string ini_mode = "resume";
#else
  double T = atof(argv[1]);
  double h = atof(argv[2]);
  unsigned long long seed = atoi(argv[3]);
  int snap_dt = atoi(argv[4]);
  int n_step = atoi(argv[5]);
  std::string ini_mode = argv[6];
#endif
  
  // set output
  char basename[255];
  char log_file[255];
  char op_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "D:\\code\\XY\\data\\";
#else
  char folder[255];
  snprintf(folder, 255, "/home/ps/data/LatticeModel/XY/L%d/", Lx);
#endif
  
  char log_folder[255];
  char op_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  snprintf(op_folder, 255, "%sop%s", folder, delimiter.c_str());
  mkdir(log_folder);
  mkdir(op_folder);
  
  int start = 0;
  int log_dt = 10000;
  int op_dt = 100;
  int op_flush_dt = snap_dt;
  
  snprintf(basename, 255, "L%d_%d_T%g_h%g_S%llu", Lx, Ly, T, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, h, -1, Lx, Ly, ini_mode);
  
  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  io::LogExporter log(log_file, n_step, log_dt, start, Lx * Ly);
  
  snprintf(op_file, 255, "%s%s_t%d.dat", op_folder, basename, start);
  io::OrderParaExporter op(op_file, start, n_step, op_dt, op_flush_dt);
  
  
  // initialize spins
  Ran myran(seed);
  XY_2<SquareLattice_2> XY_spins(Lx, Ly, T, h);
  XY_spins.ini(ini_mode, gsd, myran);

  // run
  for (int i = 1; i <= n_step; i++) {
    XY_spins.update(myran);
    log.record(i);
    op.dump(i, XY_spins);
    gsd.dump(i, XY_spins);
  }
}
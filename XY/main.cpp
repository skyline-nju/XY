#include "lattice_2.h"
#include "xy_2.h"
#include "rand.h"
#include "io2D.h"


int main(int argc, char* argv[]) {
  int Lx = 64;
  int Ly = 64;

  double T = 0.9;
  double h = 0.05;

  unsigned long long seed = 1234;
  Ran myran(seed);

  int snap_dt = 1000;
  int n_step = 100000;

  XY_2<SquareLattice_2> XY_spins(Lx, Ly, T, h);
  //XY_spins.ini_rand(myran);
  XY_spins.ini_ordered();


  // set output
  char basename[255];
  char log_file[255];
  char op_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "D:\\code\\XY\\data\\";
#else
  char folder[255];
  snprintf(folder, 255, "/home/ps/data/NRXY/L%g/", gl_l.x);
#endif

  snprintf(basename, 255, "L%d_%d_T%g_h%g_S%llu", Lx, Ly, T, h, seed);

  char log_folder[255];
  char op_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  snprintf(op_folder, 255, "%sop%s", folder, delimiter.c_str());
  mkdir(log_folder);
  mkdir(op_folder);

  int start = 0;
  int log_dt = 10000;
  int op_dt = 100;
  int op_flush_dt = 10000;

  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  io::LogExporter log(log_file, n_step, log_dt, start, Lx * Ly);

  snprintf(op_file, 255, "%s%s_t%d.dat", op_folder, basename, start);
  io::OrderParaExporter op(op_file, start, n_step, op_dt, op_flush_dt);


  for (int i = 1; i <= n_step; i++) {
    XY_spins.update(myran);
    log.record(i);
    op.dump(i, XY_spins);
  }
}
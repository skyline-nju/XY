#include "lattice_2.h"
#include "NR_kuramoto_2.h"
#include "rand.h"
#include "io2D.h"


int main(int argc, char* argv[]) {
  // set parameters
  int Lx = 400;
  int Ly = 400;

#ifdef _MSC_VER
  double T = 0.3;
  double sigma = 0.;
  double omega0 = 0.15;
  double J0 = 1;
  double J1 = 1;
  double h = 0.1;
  unsigned long long seed = 1000;
  int snap_dt = 1000;
  int n_step = 80000;
  std::string ini_mode = "resume";
#else
  double T = atof(argv[1]);
  double sigma = atof(argv[2]);
  double omega0 = atof(argv[3]);
  double J0 = 1.;
  double J1 = atof(argv[4]);
  double h = atof(argv[5]);
  unsigned long long seed = atoi(argv[6]);
  int snap_dt = atoi(argv[7]);
  int n_step = atoi(argv[8]);
  std::string ini_mode = argv[9];
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
  snprintf(folder, 255, "/home/ps/data/LatticeModel/NRKM/L%d/", Lx);
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
  
  snprintf(basename, 255, "L%d_%d_T%g_J%g_%g_o%g_s%g_h%g_S%llu", Lx, Ly, T, J0, J1, omega0, sigma, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, h, -1, Lx, Ly, ini_mode);
  
  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  io::LogExporter log(log_file, n_step, log_dt, start, Lx * Ly);
  
  snprintf(op_file, 255, "%s%s_t%d.dat", op_folder, basename, start);
  io::OrderParaExporter op(op_file, start, n_step, op_dt, op_flush_dt);
  
  
  // initialize spins
  Ran myran(seed);
  NR_Kuramoto_2<SquareLattice_2> spins(Lx, Ly, T, sigma, h, omega0, J0, J1);
  spins.ini(ini_mode, gsd, myran);

  // run
  for (int i = 1; i <= n_step; i++) {
    spins.update(myran);
    log.record(i);
    op.dump(i, spins);
    gsd.dump(i, spins);
  }
}
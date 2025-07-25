#include <mpi.h>
#include "lattice_2.h"
#include "kuramoto_2.h"
#include "rand.h"
#include "io2D.h"


int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int my_rank, tot_proc;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &tot_proc);

  ////! set parameters
  int Lx = 256;
  int Ly = 256;

#ifdef _MSC_VER
  double T = 0.3;
  double h = 0.1;
  double sigma = 0.2;
  unsigned long long seed = 1001;
  int snap_dt = 100;
  int n_step = 10000;
  std::string ini_mode = "resume";
  int proc_nx = 2;
  int proc_ny = 2;
#else
  double T = atof(argv[1]);
  double sigma = atof(argv[2]);
  double h = atof(argv[3]);
  unsigned long long seed = atoi(argv[4]);
  int snap_dt = atoi(argv[5]);
  int n_step = atoi(argv[6]);
  std::string ini_mode = argv[7];
  int proc_nx = 2;
  int proc_ny = 2;
#endif
  if (tot_proc != proc_nx * proc_ny) {
    if (my_rank == 0) {
      std::cout << "Error, total processor number should be " << proc_nx * proc_ny 
        << " instead of " << tot_proc <<  std::endl;
    }
    exit(1);
  }

  // set output
  char basename[255];
  char log_file[255];
  char op_file[255];
  char gsd_file[255];
#ifdef _MSC_VER
  char folder[] = "D:\\code\\XY\\KM_MPI\\data\\";
#else
  char folder[255];
  snprintf(folder, 255, "/home/ps/data/LatticeModel/XY/L%d/", Lx);
#endif

  char log_folder[255];
  char op_folder[255];
  snprintf(log_folder, 255, "%slog%s", folder, delimiter.c_str());
  snprintf(op_folder, 255, "%sop%s", folder, delimiter.c_str());
  if (my_rank == 0) {
    mkdir(log_folder);
    mkdir(op_folder);
  }

  int start = 0;
  int log_dt = 10000;
  int op_dt = 100;
  int op_flush_dt = snap_dt;

  snprintf(basename, 255, "L%d_%d_T%g_s%g_h%g_S%llu", Lx, Ly, T, sigma, h, seed);
  snprintf(gsd_file, 255, "%s%s.gsd", folder, basename);
  io::Snap_GSD_2 gsd(gsd_file, n_step, snap_dt, start, h, -1, Lx, Ly, ini_mode, MPI_COMM_WORLD);

  snprintf(log_file, 255, "%s%s_t%d.dat", log_folder, basename, start);
  io::LogExporter log(log_file, n_step, log_dt, start, Lx * Ly, MPI_COMM_WORLD);

  snprintf(op_file, 255, "%s%s_t%d.dat", op_folder, basename, start);
  io::OrderParaExporter op(op_file, start, n_step, op_dt, op_flush_dt, MPI_COMM_WORLD);

  {
    // initialize spins
    Ran myran(seed + my_rank);
    Kuramoto_2 <SquareLatticePadded_2> spins(Lx, Ly, sigma,
                                             proc_nx, proc_ny,
                                             MPI_COMM_WORLD, T, h);
    spins.ini(ini_mode, gsd, myran);

    // run
    for (int i = 1; i <= n_step; i++) {
      spins.update(myran);
      log.record(i);
      op.dump(i, spins);
      gsd.dump(i, spins);
    }
  }

  MPI_Finalize();
}

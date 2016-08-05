#include <mpi.h>
#include <adios.h>

struct adiosamrinfo {
  char *name;
  MPI_Comm comm;
  int rank;
  int nprocs;
  int tsteps;
  
  int numxvars;
  int maxxvars;
  char **xvarnames;

  int bufallocsize;

  int64_t gid;
};

void adiosamr_init(struct adiosamrinfo *nfo, char *method, char *name, MPI_Comm comm, int rank, int nprocs, int tsteps);

void adiosamr_addxvar(struct adiosamrinfo *nfo, char *varname);

void adiosamr_write(struct adiosamrinfo *nfo, int tstep, uint64_t cnpoints, float *points, float **xvals);
  
void adiosamr_finalize(struct adiosamrinfo *nfo);

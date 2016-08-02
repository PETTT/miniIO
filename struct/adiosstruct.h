#include <mpi.h>
#include <adios.h>

struct adiosstructinfo {
  char *name;
  MPI_Comm comm;
  int rank;
  int nprocs;
  int tsteps;
  
  int nvars;
  int maxvars;
  char **varnames;
  float **datas;

  int bufallocsize;

  int64_t gid;
};

void adiosstruct_init(struct adiosstructinfo *nfo, char *method, char *name, MPI_Comm comm, int rank, int nprocs, int tsteps);

void adiosstruct_addxvar(struct adiosstructinfo *nfo, char *varname, float *data);

void adiosstruct_write(struct adiosstructinfo *nfo, int tstep, uint64_t cnpoints, uint64_t npoints, uint64_t cstart, uint64_t *mask);
  
void adiosstruct_finalize(struct adiosstructinfo *nfo);


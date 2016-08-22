#include <mpi.h>
#include <adios.h>

struct adiosstructinfo {
  char *name;
  MPI_Comm comm;
  int rank;
  int nprocs;
  int tsteps;
  int ni, nj, nk;
  int is, cni, js, cnj, ks, cnk;
  float deltax, deltay, deltaz;
  
  int nvars;
  int maxvars;
  char **varnames;
  float **datas;

  int bufallocsize;

  int64_t gid;
};

void adiosstruct_init(struct adiosstructinfo *nfo, char *method,
               char *name, MPI_Comm comm, int rank, int nprocs,
               int tsteps, int ni, int nj, int nk, int is, int cni, int js, int cnj,
               int ks, int cnk, float deltax, float deltay, float deltaz);

void adiosstruct_addxvar(struct adiosstructinfo *nfo, char *varname, float *data);

void adiosstruct_write(struct adiosstructinfo *nfo, int tstep, int *mask);
  
void adiosstruct_finalize(struct adiosstructinfo *nfo);


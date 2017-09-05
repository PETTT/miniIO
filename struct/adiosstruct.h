/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/


#include <mpi.h>
#include <adios.h>

struct adiosstructinfo {
  char *transform;
  char *name;
  MPI_Comm comm;
  int rank;
  int nprocs;
  int tsteps;
  int ni, nj, nk;
  int is, cni, js, cnj, ks, cnk;
  float deltax, deltay, deltaz, fillvalue;

  int maxvars;

  int nrealvars;
  char **realvarnames;
  float **realdatas;
  
  int nintvars;
  char **intvarnames;
  int **intdatas;

  int bufallocsize;

  int64_t gid;
};

void adiosstruct_init(struct adiosstructinfo *nfo, char *method, char *transform,
               char *name, MPI_Comm comm, int rank, int nprocs,
               int tsteps, int ni, int nj, int nk, int is, int cni, int js, int cnj,
	       int ks, int cnk, float deltax, float deltay, float deltaz, float fillvalue);

void adiosstruct_addrealxvar(struct adiosstructinfo *nfo, char *varname, float *data);

void adiosstruct_addintxvar(struct adiosstructinfo *nfo, char *varname, int *data);

void adiosstruct_write(struct adiosstructinfo *nfo, int tstep);
  
void adiosstruct_finalize(struct adiosstructinfo *nfo);


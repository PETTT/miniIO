/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <mpi.h>

struct hdf5amrinfo {
  char *name;
  MPI_Comm comm;
  int rank;
  int nprocs;
  int tsteps;
  
  int numxvars;
  int maxxvars;
  char **xvarnames;

};

void hdf5_addxvar(struct hdf5amrinfo *nfo, char *varname);

void hdf5_write(struct hdf5amrinfo *nfo, int tstep, uint64_t cnpoints, float *points, float **xvals);

void hdf5_finalize(struct hdf5amrinfo *nfo);

void hdf5_init(struct hdf5amrinfo *nfo, char *name,
	       MPI_Comm comm, int rank, int nprocs, int tsteps);

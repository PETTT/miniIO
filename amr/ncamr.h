/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#ifndef NCAMR_H
#define NCAMR_H
#include <stdint.h>
#include <mpi.h>

struct ncamrinfo {
  char *name;
  MPI_Comm comm;
  int rank;
  int nprocs;
  int tsteps;

  int numxvars;
  int maxxvars;
  char **xvarnames;
};

void nc_addxvar(struct ncamrinfo *nfo, char *varname);

void nc_write(struct ncamrinfo *nfo, int tstep, uint64_t cnpoints, float *points, float **xvals);

void nc_finalize(struct ncamrinfo *nfo);

void nc_init(struct ncamrinfo *nfo, char *name,
	       MPI_Comm comm, int rank, int nprocs, int tsteps);

#define MAX(a,b) ((a) > (b) ? (a) : (b))

#endif

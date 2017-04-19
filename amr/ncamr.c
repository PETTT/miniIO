/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include "netcdf.h"
#include "ncamr.h"

static const int fnstrmax = 4095;

void nc_init(struct ncamrinfo *nfo, char *name,
		   MPI_Comm comm, int rank, int nprocs, int tsteps) {
  /* Set up struct, haven't decided if using all of these yet */
  nfo->name = name;
  nfo->comm = comm;
  nfo->rank = rank;
  nfo->nprocs = nprocs;
  nfo->tsteps = tsteps;

  nfo->numxvars = 0;
  nfo->maxxvars = 1000;
  nfo->xvarnames = (char **) malloc(nfo->maxxvars * sizeof(char *));

}

void nc_addxvar(struct ncamrinfo *nfo, char *varname) {
  if(nfo->numxvars > nfo->maxxvars)
    return;   /* Just ignore too many variables, for now */
  nfo->xvarnames[nfo->numxvars] = varname;
  nfo->numxvars++;
}

void nc_write(struct ncamrinfo *nfo, int tstep, uint64_t cnpoints, float *points, float **xvals) {

}

void nc_finalize(struct ncamrinfo *nfo) {
  free(nfo->xvarnames);
}

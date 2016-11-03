/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include "adiosamr.h"

static const int fnstrmax = 4095;

void adiosamr_init(struct adiosamrinfo *nfo, char *method, char *name,
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
  nfo->bufallocsize = 0;

  /* Set up ADIOS */
  adios_init_noxml(comm); 
  adios_declare_group(&nfo->gid, nfo->name, "", adios_flag_no);
  adios_select_method(nfo->gid, method, "", "");

  /* Define output variables */
  adios_define_var(nfo->gid, "rank", "", adios_integer, "", "", "");
  adios_define_var(nfo->gid, "tstep", "", adios_integer, "", "", "");
  adios_define_var(nfo->gid, "cnpoints", "", adios_unsigned_long, "", "", "");
  adios_define_var(nfo->gid, "npoints", "", adios_unsigned_long, "", "", "");
  adios_define_var(nfo->gid, "cstart", "", adios_unsigned_long, "", "", "");
}    

void adiosamr_addxvar(struct adiosamrinfo *nfo, char *varname) {
  if(nfo->numxvars > nfo->maxxvars)
    return;   /* Just ignore too many variables, for now */
  nfo->xvarnames[nfo->numxvars] = varname;
  nfo->numxvars++;
  adios_define_var(nfo->gid, varname, "", adios_real, "cnpoints", "npoints", "cstart");
}

void adiosamr_write(struct adiosamrinfo *nfo, int tstep, uint64_t cnpoints, float *points, float **xvals) {
  char fname[fnstrmax+1];
  int timedigits = 4;
  uint64_t groupsize, totalsize;
  int64_t handle;
  int ret;
  int bufneeded;
  uint64_t i, totalpoints,  cstart, *npts_all;
  

  /* numCoords = npoints * 8; */
  groupsize = sizeof(int) * 2 /*rank-numtask*/ +
    sizeof(uint64_t)*3 /*cnpoints-cstart*/ +
    sizeof(float)*cnpoints*nfo->numxvars; /*xvars*/ 

  /* Allocate buffer large enough for all data to write, if not done already */
  bufneeded = (int)(groupsize/(1024*1024));
  bufneeded += bufneeded/10 + 5;   /* Add an extra 10% & 5MB to be sure */
  if(nfo->bufallocsize < bufneeded) {
    adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, bufneeded);
    nfo->bufallocsize = bufneeded;
  }

  /* Determine global sizes */
  npts_all = (uint64_t *) malloc(nfo->nprocs * sizeof(uint64_t));
  MPI_Allgather(&cnpoints, 1, MPI_UNSIGNED_LONG_LONG,
                npts_all, 1, MPI_UNSIGNED_LONG_LONG, nfo->comm);
  for(totalpoints = 0, i = 0; i < nfo->nprocs; ++i)
      totalpoints += npts_all[i];
  for(cstart = 0, i = 0; i < nfo->rank; ++i)
      cstart += npts_all[i];
  free(npts_all);

  /* Set filename */
  snprintf(fname, fnstrmax, "%s.out.%0*d.bp", nfo->name, timedigits, tstep);

  /* Open & Write */
  ret = adios_open(&handle, nfo->name, fname, "w", nfo->comm);
  if(ret) {
    fprintf(stderr, "Error opening ADIOS file: %s\n", fname);
    return;
  }
  adios_group_size(handle, groupsize, &totalsize);

  adios_write(handle, "rank", &nfo->rank);
  adios_write(handle, "tstep", &tstep);
  adios_write(handle, "cnpoints", &cnpoints);
  adios_write(handle, "npoints", &totalpoints);
  adios_write(handle, "cstart", &cstart);
  for(i = 0; i < nfo->numxvars; i++){
    adios_write(handle, nfo->xvarnames[i], xvals[i]);
  }
  adios_close(handle);
  
}

void adiosamr_finalize(struct adiosamrinfo *nfo) {
  free(nfo->xvarnames);
  adios_finalize(nfo->rank);
}

/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>

#include "ncamr.h"

#include "netcdf.h"
#include "netcdf_par.h"

#define NCERR {if(err != NC_NOERR) {printf("(rank %d) Error at line %d: %s\n",rank,__LINE__,nc_strerror(err)); fflush(stdout); MPI_Abort(nfo->comm,1);}}

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

  char fname[fnstrmax+1];
  char dimname[fnstrmax+1];
  int timedigits = 4;
  uint64_t i, totalpoints, *npts_all;

  int ncid;
  size_t start[1] = {0};
  size_t count[1] = {0};
  size_t dims[1] = {0};
  int err = 0;
  MPI_Info info = MPI_INFO_NULL;
  uint64_t attr_data[2] = {0,0};
  int *dimids;
  int *varids;
  int varid;
  int dimid;
  int rank = nfo->rank;

  npts_all = (uint64_t *) malloc(nfo->nprocs * sizeof(uint64_t));

  MPI_Allgather(&cnpoints, 1, MPI_UNSIGNED_LONG_LONG,
                npts_all, 1, MPI_UNSIGNED_LONG_LONG, nfo->comm);

  /* Set filename */
  snprintf(fname, fnstrmax, "%s.%0*d.nc", nfo->name, timedigits, tstep);

  //if(nfo->rank == 0) {

  for(totalpoints = 0, i = 0; i < nfo->nprocs; ++i)
    totalpoints += npts_all[i];

  err = nc_create_par(fname,NC_NETCDF4|NC_MPIIO,nfo->comm,info,&ncid); NCERR;

  dims[0] = totalpoints;
  count[0] = nfo->nprocs;

  err = nc_def_dim(ncid,"phony_dim_0",nfo->nprocs,&dimid); NCERR;
  err = nc_def_var(ncid,"cnpoints",NC_UINT64,1,&dimid,&varid); NCERR;

  /* A little bit of memory allocation. */
  dimids = (int*)malloc(sizeof(int)*nfo->numxvars);
  varids = (int*)malloc(sizeof(int)*nfo->numxvars);

  /* Create the dataset with default properties, one dimension per variable.*/
  for(i = 0; i < nfo->numxvars; i++){

    /* Define phony dim. */
    snprintf(dimname,fnstrmax,"phony_dim_%jd",i+1);
    err = nc_def_dim(ncid,dimname,totalpoints,&dimids[i]); NCERR;

    /* Define variable.*/

    err = nc_def_var(ncid,nfo->xvarnames[i],NC_FLOAT,1,&dimids[i],&varids[i]); NCERR;
    err = nc_var_par_access(ncid,varids[i],NC_COLLECTIVE); NCERR;
  }

  err = nc_put_vara_long(ncid,varid,start,count,npts_all); NCERR;

  attr_data[0] = tstep;
  attr_data[1] = totalpoints;

  /* Create the data space for the attribute. */
  dims[0] = 2;
  err = nc_put_att(ncid,NC_GLOBAL,"tstep, totalpoints",NC_UINT64,2,attr_data);

  err = nc_enddef(ncid); NCERR;

  MPI_Barrier(nfo->comm);

  MPI_Info_create(&info);
  MPI_Info_set(info, "striping_factor", "1");

  MPI_Barrier(nfo->comm);

  for(start[0] = 0, i = 0; i < nfo->rank; ++i) {
    start[0] += (size_t)npts_all[i];
  }

  count[0] = (size_t)(cnpoints);

  /* Create property list for collective dataset write. */
  for(i = 0; i < nfo->numxvars; i++){
    err = nc_put_vara_float(ncid,varids[i],start,count,xvals[i]); NCERR;
  }

  MPI_Barrier(nfo->comm);
  nc_close(ncid); NCERR;

  free(npts_all);
  free(dimids);
  free(varids);

}

void nc_final(struct ncamrinfo *nfo) {
  free(nfo->xvarnames);
}

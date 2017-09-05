/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include "nccartiso.h"
#include "netcdf.h"
#include "netcdf_par.h"

static const int fnstrmax = 4095;

#define NCERR {if(err != NC_NOERR) {printf("(rank %d) Error at line %d: %s\n",rank,__LINE__,nc_strerror(err)); fflush(stdout); MPI_Abort(comm,1);}}


/*! Write nci

 */
void writenci(char *name, char *varname, MPI_Comm comm, int rank, int nprocs,
              int tstep, int ni, int nj, int nk, int is, int ie, int js, int je,
              int ks, int ke, float deltax, float deltay, float deltaz, int nci,
              int ncj, int nck, float *data)
{
  char fname[fnstrmax+1];
  int timedigits = 4;
  MPI_Info info = MPI_INFO_NULL;
  int dimid = 0;
  int ncid = 0;
  int varid = 0;
  int dimlen = 0;
  size_t start[3] = {0,0,0};
  size_t count[3] = {0,0,0};
  int dims[3] = {0,0,0};
  int dimids[3] = {0,0,0};
  int err = 0;
  int i;

  snprintf(fname, fnstrmax, "cart.%s_t%0*d.nc", varname, timedigits, tstep);
  dimlen = MAX(ni,MAX(nj,nk));
  /* Set up MPI info */
  MPI_Info_create(&info);

  err = nc_create_par(fname,NC_NETCDF4|NC_MPIIO,comm,info,&ncid); NCERR;

  /* Create dataset */
  dims[0] = (size_t)nk;
  dims[1] = (size_t)nj;
  dims[2] = (size_t)ni;

  /* Create the dimension, variable. */

  err = nc_def_dim(ncid,"phony_dim_0",dimlen,&dimid); NCERR;

  for(i = 0; i < 3; i++) {
    dimids[i] = dimid;
  }

  err = nc_def_var(ncid,varname,NC_FLOAT,3,&dimids[0],&varid); NCERR;

  err = nc_enddef(ncid); NCERR;
  MPI_Barrier(comm);

  /* Set up MPI info */
  MPI_Info_create(&info);
  //MPI_Info_set(info, "striping_factor", "1");

  /* Set up parallel access for the variable. */
  err = nc_var_par_access(ncid,varid,NC_COLLECTIVE); NCERR;

  start[0] = ks;
  start[1] = js;
  start[2] = is;

  count[0] = (size_t)nck;
  count[1] = (size_t)ncj;
  count[2] = (size_t)nci;

  err = nc_put_vara_float(ncid,varid,start,count,data); NCERR;

  MPI_Barrier(comm);

  err = nc_close(ncid); NCERR;
}

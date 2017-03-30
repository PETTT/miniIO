/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include "netcdf.h"
#include "netcdf_par.h"
#include "netcdfstruct.h"

#define NC_MAX_NAME_LENGTH 4095

#define ERR {if(err != NC_NOERR)printf("Error at line %d: %s\n",__LINE__,nc_strerror(err)); MPI_Abort(comm,1);}

void writenc(const int num_varnames, char **varnames, MPI_Comm comm, int rank,
             int nprocs, int tstep,
             int is, int js, int ks,
             int ni, int nj, int nk, int cni, int cnj, int cnk) {


  char fname[NC_MAX_NAME_LENGTH +1];
  int timedigits = 4;
  MPI_Info info = MPI_INFO_NULL;
  int err=0;
  size_t start[3] = {0,0,0};
  size_t count[3];
  size_t dimsf[3];
  size_t dimsm[3];
  int ncid;

  snprintf(fname, NC_MAX_NAME_LENGTH, "struct_t%0*d.nc", timedigits, tstep);

  dimsf[0] = nk;
  dimsf[1] = nj;
  dimsf[2] = ni;
  dimsm[0] = cnk;
  dimsm[1] = cnj;
  dimsm[2] = cni;

  if(rank == 0) {
    err = nc_create_par(fname, NC_NETCDF4|NC_MPIIO, comm, info, &ncid); ERR;
  }

  err = nc_close(ncid);

  printf("Finished writenc()\n");

}

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
#define NDIMS 3

#define ERR {if(err != NC_NOERR) {printf("Error at line %d: %s\n",__LINE__,nc_strerror(err)); fflush(stdout); MPI_Abort(comm,1);}}

/*! Write structured grid via netCDF4 API.
 *
 * Writing a structured grid for the struct mini app.
 *
 * @param[in] num_varnames Number of variable names.
 * @param[in] varnames Variable names.
 * @param[in] comm MPI_Comm variable.
 * @param[in] rank MPI_Comm_rank variable.
 * @param[in] nprocs MPI_Comm_size variable.
 * @param[in] tstep
 * @param[in] is global index starting point.
 * @param[in] js global index starting point.
 * @param[in] ks global index starting point.
 * @param[in] ni global grid size.
 * @param[in] nj global grid size.
 * @param[in] nk global grid size.
 * @param[in] cni points in this task.
 * @param[in] cnj points in this task.
 * @param[in] cnk points in this task.
 */
void writenc(const int num_varnames, char **varnames, MPI_Comm comm, int rank,
               int nprocs, int tstep,
               int is, int js, int ks,
               int ni, int nj, int nk, int cni, int cnj, int cnk) {


  char fname[NC_MAX_NAME_LENGTH +1];
  int i = 0;
  int j = 0;
  int timedigits = 4;
  MPI_Info info = MPI_INFO_NULL;
  int err=0;
  size_t start[3] = {0,0,0};
  size_t count[3];
  size_t dimsf[3];
  size_t dimsm[3];
  int ncid;
  int dimids[NDIMS];
  int *varids;
  snprintf(fname, NC_MAX_NAME_LENGTH, "struct_t%0*d.nc", timedigits, tstep);

  dimsf[0] = nk;
  dimsf[1] = nj;
  dimsf[2] = ni;
  dimsm[0] = cni;
  dimsm[1] = cnj;
  dimsm[2] = cnk;

  /* Rank = 0, create a new dataset. */
  if(rank == 0) {
    varids = (int*)malloc(sizeof(int) * num_varnames);
    /* Create file, serial fashion. */
    err = nc_create(fname,NC_NETCDF4,&ncid);
    /*create dimensions. */
    err = nc_def_dim(ncid,"k", nk, &dimids[0]); ERR;
    err = nc_def_dim(ncid,"j", nj, &dimids[1]); ERR;
    err = nc_def_dim(ncid,"i", ni, &dimids[2]); ERR;

    for(i = 0; i < num_varnames; i++) {
      	if(strcmp(varnames[i],"data") == 0 || strcmp(varnames[i],"height") == 0) {
          nc_def_var(ncid,varnames[i], NC_FLOAT, 3, &dimids[0], &varids[i]);
        } else {
          nc_def_var(ncid,varnames[i], NC_INT, 3, &dimids[0], &varids[i]);
        }
    }

    /* End definition mode & close file.*/
    err = nc_enddef(ncid); ERR;
    err = nc_close(ncid); ERR;
  }
  /* Wait for everybody to get to this point,
   * reopen file and do that thing. */
  MPI_Barrier(comm);

}

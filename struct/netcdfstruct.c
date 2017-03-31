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

#define ERR {if(err != NC_NOERR)printf("Error at line %d: %s\n",__LINE__,nc_strerror(err)); MPI_Abort(comm,1);}

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

  snprintf(fname, NC_MAX_NAME_LENGTH, "struct_t%0*d.nc", timedigits, tstep);

  dimsf[0] = nk;
  dimsf[1] = nj;
  dimsf[2] = ni;
  dimsm[0] = cnk;
  dimsm[1] = cnj;
  dimsm[2] = cni;

  /* Rank = 0, create a new dataset. */
  if(rank == 0) {
    /* Create file */
    printf("Rank == 0: Defining file.\n"); fflush(stdout);
    err = nc_create_par(fname, NC_NETCDF4|NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid); ERR;

    /*create dimensions. */
    err = nc_def_dim(ncid,"x", NC_UNLIMITED, dimids); ERR;
    err = nc_def_dim(ncid,"y", NC_UNLIMITED, &dimids[1]); ERR;
    err = nc_def_dim(ncid,"z", NC_UNLIMITED, &dimids[2]); ERR;

    /* End definition mode & close file.*/
    err = nc_enddef(ncid); ERR;
    err = nc_close(ncid); ERR;
  }

  MPI_Barrier(comm);
  MPI_Info_create(&info);
  MPI_Info_set(info, "striping_factor", "1");

  printf("Finished writenc()\n");

}

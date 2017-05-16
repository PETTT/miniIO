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
#include <inttypes.h>

#include "netcdf.h"
#include "netcdf_par.h"
static const int fnstrmax = 4095;

#define NCERR {if(err != NC_NOERR) {printf("(rank %d) Error at line %d: %s\n",rank,__LINE__,nc_strerror(err)); fflush(stdout); MPI_Abort(comm,1);}}

/*! Write netCDF P data/info

 */
void writencp(char *name, char *varname, MPI_Comm comm, int rank, int nprocs,
              int tstep, uint64_t ntris, float *points, float *norms,
              float *xvals, char *xname)
{
  char fname[fnstrmax+1];
  char fname_xdmf[fnstrmax+1];
  int timedigits = 4;
  MPI_Info info = MPI_INFO_NULL;
  uint64_t *rntris;   /* All triangle counts from each task */
  uint64_t tot_tris;
  int grp_id_grid_pts;
  int var_id_xyz;
  int var_id_normals;
  int var_id_conn;
  int var_id_xname;
  size_t start[1];
  size_t count[1];

  int var_dim_phony_0;
  int var_dim_phony_1;

  int ncid = 0;

  int j;
  uint64_t *temparr;
  int err = 0;

  snprintf(fname, fnstrmax, "cartiso_t%0*d.nc", timedigits, tstep);

  /* Gather tri counts, in case some are zero, to leave them out */
  rntris = (uint64_t *) malloc(nprocs*sizeof(uint64_t));

  MPI_Allgather(&ntris, 1, MPI_UNSIGNED_LONG_LONG, rntris, 1, MPI_LONG_LONG, comm); NCERR;

  uint64_t i;

  tot_tris=0;
  for (i=0; i<nprocs; i++) {
    tot_tris = tot_tris + rntris[i];
  }

  //if(rank == 0) {
    /* Create file with parallel I/O. */
    err = nc_create_par(fname,NC_NETCDF4|NC_MPIIO,comm,info,&ncid); NCERR;

    /* Create Grid Group */
    err = nc_def_grp(ncid,"grid points",&grp_id_grid_pts); NCERR;

    /*
     * Create top level dimension, variables.
     */
    err = nc_def_dim(ncid,"Phony_Dimension_1",NC_UNLIMITED,&var_dim_phony_1); NCERR;
    err = nc_def_var(ncid,"conn",NC_UINT64,1,&var_dim_phony_1,&var_id_conn); NCERR;

    if(xvals) {
      err = nc_def_var(ncid,xname,NC_FLOAT,1,&var_dim_phony_1,&var_id_xname); NCERR;
    }

    /*
     * Create group dimension, variables.
     */

    /* Create phony dimension. */
    err = nc_def_dim(grp_id_grid_pts,"Phony_Dimension_0",NC_UNLIMITED,&var_dim_phony_0); NCERR;

    /* Create group variables xyz and normals */
    err = nc_def_var(grp_id_grid_pts,"xyz",NC_FLOAT,1,&var_dim_phony_0,&var_id_xyz);
    err = nc_def_var(grp_id_grid_pts,"Normals",NC_FLOAT,1,&var_dim_phony_0,&var_id_normals);

    err = nc_enddef(ncid); NCERR;

    //}

  MPI_Barrier(comm);

  if(tot_tris == 0)
    return;

  /* Set up MPI info */
  MPI_Info_create(&info);
  MPI_Info_set(info, "striping_factor", "1");

  /* Set up variables for parallel I/O Access. */
  err = nc_var_par_access(ncid,var_id_conn,NC_COLLECTIVE); NCERR;
  err = nc_var_par_access(ncid,var_id_xname,NC_COLLECTIVE); NCERR;
  err = nc_var_par_access(grp_id_grid_pts,var_id_xyz,NC_COLLECTIVE); NCERR;
  err = nc_var_par_access(grp_id_grid_pts,var_id_normals,NC_COLLECTIVE); NCERR;

  /*
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  /* XYZ */
  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 9*rntris[j];
  }

  count[0] = (size_t)ntris*9;
  err = nc_put_vara_float(grp_id_grid_pts,var_id_xyz,start,count,points); NCERR;

  /* End XYZ */

  /*
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */
  /* Normals */
  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 3*rntris[j];
  }

  count[0] = (size_t)ntris*3;

  err = nc_put_vara_float(grp_id_grid_pts,var_id_normals,start,count,norms); NCERR;

  /* End Normals */

  /* Conn */
   /*
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  count[0] = (size_t)ntris*3;

  /* Select hyperslab in the file.*/
  temparr = (uint64_t *) malloc(ntris*3*sizeof(uint64_t));

  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 3*rntris[j];
  }

  for(j = 0; j < ntris*3; j++)
    temparr[j] = start[0] + j;

  err = nc_put_vara_long(ncid,var_id_conn,start,count,temparr); NCERR;

  /* End Conn */


  /* xvals */
  if(xvals) {

    start[0] = 0;
    count[0] = (size_t)ntris*3;

    err = nc_put_vara_float(ncid,var_id_xname,start,count,xvals);
  }
  /* End xvals */

  err = nc_close(ncid); NCERR;

  free(temparr);
  free(rntris);

}

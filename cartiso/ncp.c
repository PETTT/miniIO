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

  int var_dim_phony_0;
  int var_dim_phony_1;

  int ncid = 0;

  int j;
  uint64_t *temparr;
  int err = 0;

  snprintf(fname, fnstrmax, "cartiso_t%0*d.h5", timedigits, tstep);

  /* Gather tri counts, in case some are zero, to leave them out */
  rntris = (uint64_t *) malloc(nprocs*sizeof(uint64_t));

  MPI_Allgather(&ntris, 1, MPI_UNSIGNED_LONG_LONG, rntris, 1, MPI_LONG_LONG, comm); NCERR;

  uint64_t i;

  tot_tris=0;
  for (i=0; i<nprocs; i++) {
    tot_tris = tot_tris + rntris[i];
  }

  if(rank == 0) {
    /* Create file with parallel I/O. */
    err = nc_create_par(fname,NC_NETCDF4|NC_MPIIO,comm,info,&ncid); NCERR;

    /* Create Grid Group */
    err = nc_def_grp(ncid,"grid points",&grp_id_grid_pts); NCERR;

    /*
     * Create top level dimension, variables.
     */
    err = nc_def_dim(ncid,"Phony Dimension 1",&var_dim_phony_1); NCERR;
    err = nc_def_var(ncid,"conn",NC_UINT64,&var_dim_phony_1,&var_id_conn); NCERR;

    if(xvals) {
      err = nc_def_var(ncid,xname,NC_FLOAT,&var_dim_phony_1,&var_id_xname); NCERR;
    }

    /*
     * Create group dimension, variables.
     */

    /* Create phony dimension. */
    err = nc_def_dim(grp_id_grid_pts,"Phony Dimension 0",&var_dim_phony_0); NCERR;

    /* Create group variables xyz and normals */
    err = nc_def_var(grp_id_grid_pts,"xyz",NC_FLOAT,&var_dim_phony_0,&var_id_xyz);
    err = nc_def_var(grp_id_grid_pts,"Normals",NC_FLOAT,&var_dim_phony_0,&var_id_normals);

    err = nc_enddef(ncid); NCERR;

  }

  MPI_Barrier(comm);

  if(tot_tris == 0)
    return;

  /* Set up MPI info */
  MPI_Info_create(&info);
  MPI_Info_set(info, "striping_factor", "1");

  /* Set up variables for parallel I/O Access. */
  err = nc_var_par_access(ncid,var_id_conn,NC_COLLECTIVE); NCERR;
  err = nc_var_par_access(ncid,var_id_xname,NC_COLLECTIVE); NCERR;
  err = nc_var_par_access(grp_id_grid_points,var_id_xyz,NC_COLLECTIVE); NCERR;
  err = nc_var_par_access(grp_id_grid_points,var_id_normals,NC_COLLECTIVE); NCERR;

  /*
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 9*rntris[j];
  }

  count[0] = (hsize_t)ntris*9;


  /* Select hyperslab in the file.*/

  err = nc_put_vara_float(var_id_xyz,start,count,points); NCERR;
  // err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, points);
  //  if( err < 0) {
  //    fprintf(stderr, "writencp error: could not write datset %s \n", "xyz");
  //    MPI_Abort(comm, 1);
  //  }

  //  err = H5Dclose(did);



  //  did = H5Dopen(group_id, "Normals",H5P_DEFAULT);
  /*
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  /* XYZ */
  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 9*rntris[j];
  }

  count[0] = (hsize_t)ntris*9;

  //memspace = H5Screate_simple(1, count, NULL);

  /* Select hyperslab in the file.*/
  //filespace = H5Dget_space(did);
  //H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );
  err = nc_put_vara_float(var_id_xyz,start,count,points); NCERR;

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

  count[0] = (hsize_t)ntris*3;

  err = nc_put_vara_float(var_id_normals,start,count,norms); NCERR;

  /* End Normals */

  /* Conn */
  did = H5Dopen(file_id, "conn",H5P_DEFAULT);

  /*
   * Each process defines dataset in memory and writes it to the hyperslab
   * in the file.
   */

  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 3*rntris[j];
  }

  count[0] = (hsize_t)ntris*3;

  memspace = H5Screate_simple(1, count, NULL);

  /* Select hyperslab in the file.*/
  filespace = H5Dget_space(did);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );
  temparr = (uint64_t *) malloc(ntris*3*sizeof(uint64_t));



  start[0] = 0;
  for (j=0; j<rank; j++) {
    start[0] = start[0] + 3*rntris[j];
  }

  for(j = 0; j < ntris*3; j++)
    temparr[j] = start[0] + j;

  err = H5Dwrite(did, H5T_NATIVE_ULLONG, memspace, filespace, plist_id, temparr);

  err = H5Dclose(did);

  if(xvals) {
    did = H5Dopen(file_id, xname, H5P_DEFAULT);
    err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xvals);
    H5Dclose(did);
  }

  err = H5Sclose(filespace);
  err = H5Sclose(memspace);

  if(H5Pclose(plist_id) < 0)
    printf("writehdf5p error: Could not close property list \n");

  free(temparr);
  free(rntris);
  if(H5Fclose(file_id) != 0)
    printf("writehdf5p error: Could not close HDF5 file \n");

}

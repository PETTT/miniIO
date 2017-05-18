/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>

#include "ncamr.h"

#include "netcdf.h"
#include "netcdf_par.h"


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

  npts_all = (uint64_t *) malloc(nfo->nprocs * sizeof(uint64_t));

  MPI_Allgather(&cnpoints, 1, MPI_UNSIGNED_LONG_LONG,
                npts_all, 1, MPI_UNSIGNED_LONG_LONG, nfo->comm);

  /* Set filename */
  snprintf(fname, fnstrmax, "%s.%0*d.nc", nfo->name, timedigits, tstep);

  //if(nfo->rank == 0) {

  for(totalpoints = 0, i = 0; i < nfo->nprocs; ++i)
    totalpoints += npts_all[i];

  err = nc_create_par(fname,NC_NETCDF4|NC_MPIIO,comm,info,&ncid); NCERR;

  dims[0] = totalpoints;
  count[0] = nfo->nprocs;

  err = nc_def_dim(ncid,"phony_dim_0",nfo->nprocs,&dimid); NCERR;
  err = nc_def_var(ncid,"cnpoints",NC_LONG,1,&dimid,&varid); NCERR;

  /* A little bit of memory allocation. */
  dimids = (int*)malloc(sizeof(int)*nfo->numxvars);
  varids = (int*)malloc(sizeof(int)*nfo->numxvars);

  /* Create the dataset with default properties, one dimension per variable.*/
  for(i = 0; i < nfo->numxvars; i++){

    /* Define phony dim. */
    snprintf(dimname,fnstrmax,"phony_dim_%d",i+1);
    err = nc_def_dim(ncid,dimname,NC_UNLIMITED,&dimids[i]); NCERR;

    /* Define variable.*/

    err = nc_def_var(ncid,nfo->xvarnames[i],NC_FLOAT,1,&dimids[i],&varids[i]); NCERR;
    err = nc_var_par_access(ncid,varids[i],NC_COLLECTIVE); NCERR;
  }


  //dims[0] = nfo->nprocs;
  //filespace = H5Screate_simple(1, dims, NULL);
  err = nc_enddef(ncid); NCERR;

  MPI_Barrier(comm);

  //did = H5Dcreate(file_id, "cnpoints", H5T_NATIVE_ULLONG, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  err = nc_put_vara_long(ncid,varid,start,count,npts_all); NCERR;
  //err = H5Dwrite(did, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, npts_all);

  /** WRITE ATTRIBUTE INFORMATION **/

  attr_data[0] = tstep;
  attr_data[1] = totalpoints;

  /* Create the data space for the attribute. */
  dims[0] = 2;
  filespace = H5Screate_simple(1, dims, NULL);

  /* Create a dataset attribute. */
  attr_id = H5Acreate2 (file_id, "tstep, totalpoints", H5T_NATIVE_ULLONG, filespace,
                        H5P_DEFAULT, H5P_DEFAULT);

  /* Write the attribute data. */
  err = H5Awrite(attr_id, H5T_NATIVE_ULLONG, attr_data);

  /* Close the attribute. */
  err = H5Aclose(attr_id);

  /* Close the dataspace. */
  err = H5Sclose(filespace);

  /* Close the file */
  H5Fclose(file_id);

  //} /*rank == 0 */

  /* Set up MPI info */
  MPI_Info_create(&info);
  MPI_Info_set(info, "striping_factor", "1");

  /* Set up file access property list with parallel I/O access */
  if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
    printf("hdf5 error: Could not create property list \n");
    MPI_Abort(nfo->comm, 1);
  }

  H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

  if(H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info) < 0) {
    printf("hdf5 error: Could not create property list \n");
    MPI_Abort(nfo->comm, 1);
  }

  MPI_Barrier(nfo->comm);
  if( (file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id)) < 0) {
    fprintf(stderr, "writehdf5p error: could not open %s \n", fname);
    MPI_Abort(nfo->comm, 1);
  }

  if(H5Pclose(plist_id) < 0) {
    printf("hdf5 error: Could not close property list \n");
    MPI_Abort(nfo->comm, 1);
  }
  for(start[0] = 0, i = 0; i < nfo->rank; ++i) {
    start[0] += (hsize_t)npts_all[i];
  }

  count[0] = (hsize_t)(cnpoints);

  memspace = H5Screate_simple(1, count, NULL);

  /* Create property list for collective dataset write. */
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

  for(i = 0; i < nfo->numxvars; i++){

    did = H5Dopen(file_id, nfo->xvarnames[i], H5P_DEFAULT);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */

    /* Select hyperslab in the file.*/
    filespace = H5Dget_space(did);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );

    err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xvals[i]);

    if( err < 0) {
      fprintf(stderr, "hdf5 error: could not write datset %s \n", nfo->xvarnames[i]);
      MPI_Abort(nfo->comm, 1);
    }

    err = H5Dclose(did);
  }

  if(H5Sclose(memspace)  != 0)
    printf("hdf5 error: Could not close memory space \n");

  if(H5Pclose(plist_id) < 0)
    printf("hdf5 error: Could not close property list \n");

  if(H5Fclose(file_id) != 0)
    printf("hdf5 error: Could not close HDF5 file \n");

  free(npts_all);
  free(dimids);
  free(varids);

}

void nc_finalize(struct ncamrinfo *nfo) {
  free(nfo->xvarnames);
}

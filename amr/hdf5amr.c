/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include "hdf5.h"
#include "hdf5amr.h"

static const int fnstrmax = 4095;

void hdf5_init(struct hdf5amrinfo *nfo, char *name,
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

void hdf5_addxvar(struct hdf5amrinfo *nfo, char *varname) {
  if(nfo->numxvars > nfo->maxxvars)
    return;   /* Just ignore too many variables, for now */
  nfo->xvarnames[nfo->numxvars] = varname;
  nfo->numxvars++;
}

void hdf5_write(struct hdf5amrinfo *nfo, int tstep, uint64_t cnpoints, float *points, float **xvals) {
  char fname[fnstrmax+1];
  int timedigits = 4;
  uint64_t i, totalpoints, *npts_all;
  
  hid_t file_id;
  hid_t plist_id;
  hid_t memspace;
  hid_t filespace;
  hid_t attr_id;
  hid_t did;
  hsize_t start[1], count[1];
  hsize_t dims[1];
  herr_t err;
  MPI_Info info = MPI_INFO_NULL;
  uint64_t attr_data[2];

  npts_all = (uint64_t *) malloc(nfo->nprocs * sizeof(uint64_t));

  MPI_Allgather(&cnpoints, 1, MPI_UNSIGNED_LONG_LONG,
                npts_all, 1, MPI_UNSIGNED_LONG_LONG, nfo->comm);

  /* Set filename */
  snprintf(fname, fnstrmax, "%s.%0*d.h5", nfo->name, timedigits, tstep);

  if(nfo->rank == 0) {

    for(totalpoints = 0, i = 0; i < nfo->nprocs; ++i)
      totalpoints += npts_all[i];

    if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	fprintf(stderr, "writehdf5 error: could not create %s \n", fname);
	MPI_Abort(nfo->comm, 1);
      }

    dims[0] = totalpoints; 
    filespace = H5Screate_simple(1, dims, NULL);

    /* Create the dataset with default properties */
      
    for(i = 0; i < nfo->numxvars; i++){

      did = H5Dcreate(file_id, nfo->xvarnames[i], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dclose(did); 

    }
    /* Close filespace. */
    H5Sclose(filespace);
    
    dims[0] = nfo->nprocs; 
    filespace = H5Screate_simple(1, dims, NULL);

    did = H5Dcreate(file_id, "cnpoints", H5T_NATIVE_ULLONG, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    err = H5Dwrite(did, H5T_NATIVE_ULLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, npts_all);

    H5Dclose(did);
    H5Sclose(filespace);

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

  }

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
  
  H5Pset_all_coll_metadata_ops(plist_id, 1 );
  H5Pset_coll_metadata_write(plist_id, 1);

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

}

void hdf5_finalize(struct hdf5amrinfo *nfo) {
  free(nfo->xvarnames);
}

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

#include "hdf5.h"

static const int fnstrmax = 4095;

void
write_xdmf_xml(char *fname, char *fname_xdmf, uint64_t npoints, char *xname);

void writehdf5p(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, uint64_t ntris, float *points, float *norms, 
		float *xvals, char *xname, hsize_t *h5_chunk, int hdf5_compress)
{
    char fname[fnstrmax+1];
    char fname_xdmf[fnstrmax+1];
    int timedigits = 4;
    MPI_Info info = MPI_INFO_NULL;
    uint64_t *rntris;   /* All triangle counts from each task */
    uint64_t tot_tris;

    hid_t file_id;
    hid_t plist_id;
    hid_t group_id;
    hid_t memspace;
    hid_t filespace;
    hid_t did;
    hsize_t start[1], count[1];
    hsize_t dims[1];
    int j;
    uint64_t *temparr;
    herr_t err;
    hid_t chunk_pid;
    

    snprintf(fname, fnstrmax, "cartiso_t%0*d.h5", timedigits, tstep);
    snprintf(fname_xdmf, fnstrmax, "cartiso_t%0*d.xmf", timedigits, tstep);

    /* Gather tri counts, in case some are zero, to leave them out */
    rntris = (uint64_t *) malloc(nprocs*sizeof(uint64_t));

    MPI_Allgather(&ntris, 1, MPI_UNSIGNED_LONG_LONG, rntris, 1, MPI_LONG_LONG, comm);

    uint64_t i;
   
    tot_tris=0;
    for (i=0; i<nprocs; i++) {
      tot_tris = tot_tris + rntris[i];
    }

    if(rank == 0) {

      if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	fprintf(stderr, "writehdf5p error: could not create %s \n", fname);
	MPI_Abort(comm, 1);
      }

      dims[0] = 9*(hsize_t)tot_tris;

      filespace = H5Screate_simple(1, dims, NULL);

      /* Create Grid Group */
      group_id = H5Gcreate(file_id, "grid points", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);

      if(h5_chunk && tot_tris != 0) {

	H5Pset_layout(chunk_pid, H5D_CHUNKED);

	if( H5Pset_fill_time(chunk_pid, H5D_FILL_TIME_NEVER) < 0 ) {
	  printf("writehdf5 error: Could not set fill time\n");
	  MPI_Abort(comm, 1);
	}
	hsize_t chunk = dims[0]*(h5_chunk[0]/100.);
	H5Pset_chunk(chunk_pid, 1, &chunk);

	if(hdf5_compress == 1) {

	  /* Set ZLIB / DEFLATE Compression using compression level 6. */
	  H5Pset_deflate (chunk_pid, 6);
	  
	  /* Uncomment these lines to set SZIP Compression
	     szip_options_mask = H5_SZIP_NN_OPTION_MASK;
	     szip_pixels_per_block = 16;
	     status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
	  */
	}
      }
      /* Create the dataset with default properties */
      did = H5Dcreate(group_id, "xyz", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Dclose(did);

      /* Create the dataset with default properties and close filespace. */
      did = H5Dcreate(group_id, "Normals", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Dclose(did);

      if(h5_chunk)
	H5Pclose(chunk_pid);

      H5Sclose(filespace);
      H5Gclose(group_id);

      /* Create connectivity dataset */
      dims[0] = 3*(hsize_t)tot_tris;
      filespace = H5Screate_simple(1, dims, NULL);

      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);

      if(h5_chunk && tot_tris != 0) {

	H5Pset_layout(chunk_pid, H5D_CHUNKED);

	if( H5Pset_fill_time(chunk_pid, H5D_FILL_TIME_NEVER) < 0 ) {
	  printf("writehdf5 error: Could not set fill time\n");
	  MPI_Abort(comm, 1);
	}
	hsize_t chunk = dims[0]*(h5_chunk[0]/100.);
	H5Pset_chunk(chunk_pid, 1, &chunk);

	if(hdf5_compress == 1) {

	  /* Set ZLIB / DEFLATE Compression using compression level 6. */
	  H5Pset_deflate (chunk_pid, 6);
	  
	  /* Uncomment these lines to set SZIP Compression
	     szip_options_mask = H5_SZIP_NN_OPTION_MASK;
	     szip_pixels_per_block = 16;
	     status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
	  */
	}
      }

      /* Create the dataset with default properties and close filespace. */
      did = H5Dcreate(file_id, "conn", H5T_NATIVE_ULLONG, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Dclose(did);

      if(xvals) {
	did = H5Dcreate(file_id, xname, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
	H5Dclose(did);
      }   
      if(h5_chunk)
	H5Pclose(chunk_pid);
      H5Sclose(filespace);
      H5Fclose(file_id);

      /* Create xdmf file for timestep */
      write_xdmf_xml(fname, fname_xdmf, tot_tris, xname);
    }
    
    MPI_Barrier(comm);

    if(tot_tris == 0)
      return;
    
    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");

    /* Set up file access property list with parallel I/O access */
    if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

    if(H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    if( (file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id)) < 0) {
      fprintf(stderr, "writehdf5p error: could not open %s \n", fname);
      MPI_Abort(comm, 1);
    }

    if(H5Pclose(plist_id) < 0) {
      printf("writehdf5p error: Could not close property list \n");
      MPI_Abort(comm, 1);
    }

    group_id = H5Gopen(file_id, "grid points", H5P_DEFAULT);
    did = H5Dopen(group_id, "xyz",H5P_DEFAULT);
    
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */

    start[0] = 0;
    for (j=0; j<rank; j++) {
      start[0] = start[0] + 9*rntris[j];
    }

    count[0] = (hsize_t)ntris*9;

    memspace = H5Screate_simple(1, count, NULL);
    if(count[0] == 0)
      H5Sselect_none(memspace);

    /* Select hyperslab in the file.*/
    filespace = H5Dget_space(did);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, points);
    if( err < 0) {
      fprintf(stderr, "writehdf5p error: could not write datset %s \n", "xyz");
      MPI_Abort(comm, 1);
    }
    err = H5Dclose(did);
    
    did = H5Dopen(group_id, "Normals",H5P_DEFAULT);
    /* 
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */

    start[0] = 0;
    for (j=0; j<rank; j++) {
      start[0] = start[0] + 9*rntris[j];
    }


    count[0] = (hsize_t)ntris*9;
      
    memspace = H5Screate_simple(1, count, NULL);
      
    /* Select hyperslab in the file.*/
    filespace = H5Dget_space(did);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );

    err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, norms);
    err = H5Dclose(did);

    err = H5Sclose(filespace);
    err = H5Sclose(memspace);
    err = H5Gclose(group_id);

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

void
write_xdmf_xml(char *fname, char *fname_xdmf, uint64_t npoints, char *xname)
{
    FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh.
     */
 
    xmf = fopen(fname_xdmf, "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "  <Grid Name=\"Unstructured Mesh\">\n");
    fprintf(xmf, "    <Topology TopologyType=\"Triangle\" NumberOfElements=\"%" PRIu64"\">\n", npoints);
    fprintf(xmf, "      <DataItem Dimensions=\"%" PRIu64"\" Format=\"HDF\">\n", npoints*3);
    fprintf(xmf, "      %s:/conn\n",fname);
    fprintf(xmf, "      </DataItem>\n");
    fprintf(xmf, "    </Topology>\n");
    fprintf(xmf, "    <Geometry GeometryType=\"XYZ\">\n");
    fprintf(xmf, "       <DataItem Name=\"XYZ\" Dimensions=\"%" PRIu64"\" Format=\"HDF\">\n", 9*npoints);
    fprintf(xmf, "        %s:/grid points/xyz\n",fname);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "    </Geometry>\n");
    fprintf(xmf, "    <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", xname);
    fprintf(xmf, "       <DataItem Dimensions=\"%" PRIu64"\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", npoints*3);
    fprintf(xmf, "        %s:/%s\n", fname, xname);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "    </Attribute>\n");
    fprintf(xmf, "    <Attribute Name=\"Normals\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%" PRIu64"\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", npoints*9);
    fprintf(xmf, "        %s:/grid points/Normals\n", fname);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "    </Attribute>\n");
    fprintf(xmf, "  </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

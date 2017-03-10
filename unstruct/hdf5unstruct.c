/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <stdint.h>
#include <stdio.h>
#include <math.h>

#include <pdirs.h>
#include "hdf5unstruct.h"

uint64_t nelems_in[2];
uint64_t nelems_out[2];

void
write_xdmf_xml(char *fname, char *fname_xdmf, uint64_t npoints);

static const int fnstrmax = 4095;

void writehdf5(char *name, MPI_Comm comm, int tstep, uint64_t npoints, uint64_t nptstask, 
               float *xpts, float *ypts, float *zpts, uint64_t nelems3, uint64_t *conns3,
               uint64_t nelems2, uint64_t *conns2, char *varname, float *data, hsize_t *h5_chunk, int hdf5_compress)
{
    char dirname[fnstrmax+1];
    char fname[fnstrmax+1];
    char rel_fname[fnstrmax+1];
    char fname_xdmf[fnstrmax+1];
    int rank, nprocs;
    int timedigits = 4;
    MPI_Info info = MPI_INFO_NULL;

    hid_t file_id;
    hid_t plist_id;
    hid_t group_id;
    hid_t memspace;
    hid_t filespace;
    
    hid_t did[3];
    hsize_t start[1], count[1];
    hsize_t block, *pblock=NULL;
    hsize_t dims[1];
    herr_t err;
    hid_t chunk_pid;
    hsize_t chunk;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* Make dir for all output and subdir for timestep */
    snprintf(dirname, fnstrmax, "%s.hdf5.d", name);
    mkdir1task(dirname, comm);
    snprintf(dirname, fnstrmax, "%s.hdf5.d/t%0*d.d", name, timedigits, tstep);
    mkdir1task(dirname, comm);

    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");    

    chkdir1task(dirname, comm);

    snprintf(fname, fnstrmax, "unstruct.hdf5.d/t%0*d.d/r.h5", timedigits, tstep);
    snprintf(rel_fname, fnstrmax, "t%0*d.d/r.h5", timedigits, tstep);
    snprintf(fname_xdmf, fnstrmax, "unstruct.hdf5.d/t%0*d.d.xmf", timedigits, tstep);

    /* Set up file access property list with parallel I/O access */
    if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
      printf("writehdf5 error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

    if(H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info) < 0) {
      printf("writehdf5 error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }
    
    if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id)) < 0) {
      fprintf(stderr, "writehdf5 error: could not open %s \n", fname);
      MPI_Abort(comm, 1);
    }
    
    if(H5Pclose(plist_id) < 0) {
      printf("writehdf5 error: Could not close property list \n");
      MPI_Abort(comm, 1);
    }

    /* Optional grid points */
    if(xpts && ypts && zpts) {
      /* Create the dataspace for the dataset. */
      dims[0] = (hsize_t)npoints;
      filespace = H5Screate_simple(1, dims, NULL);

      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);
      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	chunk = h5_chunk[0];
	H5Pset_chunk(chunk_pid, 1, &chunk);
      } 

      if(hdf5_compress == 1) {

	/* Set ZLIB / DEFLATE Compression using compression level 6. */
	H5Pset_deflate (chunk_pid, 6);

	/* Uncomment these lines to set SZIP Compression
	   szip_options_mask = H5_SZIP_NN_OPTION_MASK;
	   szip_pixels_per_block = 16;
	   status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
	*/
      }

      /* Create Grid Group */
      group_id = H5Gcreate(file_id, "grid points", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      /* Create the dataset with default properties and close filespace. */
      did[0] = H5Dcreate(group_id, "x", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      did[1] = H5Dcreate(group_id, "y", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      did[2] = H5Dcreate(group_id, "z", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Sclose(filespace);

      /* 
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */    
      start[0] =(hsize_t)(nptstask*rank);
      count[0] =(hsize_t)nptstask;
      if(h5_chunk) {
	block = 1;
	pblock = &block;
      }
      memspace = H5Screate_simple(1, count, NULL);
      
      /* Select hyperslab in the file.*/
      filespace = H5Dget_space(did[0]);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, pblock );

      /* Create property list for collective dataset write. */
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

      err = H5Dwrite (did[0], H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xpts);
      err = H5Dwrite (did[1], H5T_NATIVE_FLOAT, memspace, filespace, plist_id, ypts);
      err = H5Dwrite (did[2], H5T_NATIVE_FLOAT, memspace, filespace, plist_id, zpts);
      
      err = H5Dclose(did[0]);
      err = H5Dclose(did[1]);
      err = H5Dclose(did[2]);
      
      err = H5Sclose(filespace);
      err = H5Sclose(memspace);
      err = H5Gclose(group_id);
    
      if(h5_chunk)
	H5Pclose(chunk_pid);
    }
    if(H5Pclose(plist_id) < 0)
      printf("writehdf5 error: Could not close property list \n");

    nelems_in[0] = nelems3 ;
    nelems_in[1] = nelems2 ;

    MPI_Allreduce( nelems_in, nelems_out, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
   
    //MSB is it possible that some processors have 0?

    /* Optional grid connections, writes a 64-bit 0 if no connections */
    if(conns3 && nelems3) {

    /* Create the dataspace for the dataset. */
      dims[0] = (hsize_t)nelems_out[0]*6;
      filespace = H5Screate_simple(1, dims, NULL);


      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);
      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	if(dims[0]%h5_chunk[1] == 0) {
	  chunk = dims[0]/h5_chunk[1];
	} else {
	  printf("writehdf5 error: conns3 not evenly divisible by chunk size [1] \n");
	  MPI_Abort(comm, 1);
	}
	H5Pset_chunk(chunk_pid, 1, &chunk);
      } 
      
      /* Create the dataset with default properties and close filespace. */
      did[0] = H5Dcreate(file_id, "conns3", H5T_NATIVE_ULLONG, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Sclose(filespace);

      /* 
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */    
      start[0] =(hsize_t)(nelems3*6*rank);
      count[0] =(hsize_t)nelems3*6;
      if(h5_chunk) {
	block = 1;
	pblock = &block;
      }
      
      memspace = H5Screate_simple(1, count, NULL);
      
      /* Select hyperslab in the file.*/
      filespace = H5Dget_space(did[0]);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );

      /* Create property list for collective dataset write. */
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      err = H5Dwrite (did[0], H5T_NATIVE_ULLONG, memspace, filespace, plist_id, conns3);
      

      err = H5Dclose(did[0]);
      
      err = H5Sclose(filespace);
      err = H5Sclose(memspace);
      if(h5_chunk)
	H5Pclose(chunk_pid);

    }

/*     Optional 2D surface triangle connections, writes a 64-bit 0 if none */
    if(conns2 && nelems2) {

      /* Create the dataspace for the dataset. */
      dims[0] = (hsize_t)nelems_out[1]*3;
      filespace = H5Screate_simple(1, dims, NULL);

      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);
      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	if(dims[0]%h5_chunk[2] == 0) {
	  chunk = dims[0]/h5_chunk[2];
	} else {
	  printf("writehdf5 error: conns2 not evenly divisible by chunk size [2] \n");
	  MPI_Abort(comm, 1);
	}

	H5Pset_chunk(chunk_pid, 1, &chunk);
      }

      /* Create the dataset with default properties and close filespace. */
      did[0] = H5Dcreate(file_id, "conns2", H5T_NATIVE_ULLONG, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Sclose(filespace);

      /* 
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */  
      start[0] =(hsize_t)(nelems2*3*rank);
      count[0] =(hsize_t)nelems2*3;
      if(h5_chunk) {
	block = 1;
	pblock = &block;
      }
      
      
      memspace = H5Screate_simple(1, count, NULL);
      
      /* Select hyperslab in the file.*/
      filespace = H5Dget_space(did[0]);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );

      /* Create property list for collective dataset write. */
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      err = H5Dwrite (did[0], H5T_NATIVE_ULLONG, memspace, filespace, plist_id, conns3);
      

      err = H5Dclose(did[0]);
      if(h5_chunk)
	H5Pclose(chunk_pid);
      err = H5Sclose(filespace);
      err = H5Sclose(memspace);

    } 

    /* Optional variable data, starting with number of variables */
    if(data && varname) {
      /* Create the dataspace for the dataset. */
      dims[0] = (hsize_t)npoints;
      filespace = H5Screate_simple(1, dims, NULL);

      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);
      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	chunk = dims[0];
	H5Pset_chunk(chunk_pid, 1, &chunk);
      }

      /* Create the dataset with default properties and close filespace. */
      did[0] = H5Dcreate(file_id, "vars", H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Sclose(filespace);

      /* 
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */
      start[0] =(hsize_t)(nptstask*rank);
      count[0] =(hsize_t)nptstask;
      if(h5_chunk) {
	block = 1;
	pblock = &block;
      }
      
      memspace = H5Screate_simple(1, count, NULL);
      
      /* Select hyperslab in the file.*/
      filespace = H5Dget_space(did[0]);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );

      /* Create property list for collective dataset write. */
      plist_id = H5Pcreate(H5P_DATASET_XFER);
      H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      if(H5Dwrite (did[0], H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data) < 0) {
	printf("writehdf5 error: Could not write HDF5 file \n");
	MPI_Abort(comm, 1);
      }
      
      if(H5Dclose(did[0]) ){
	printf("writehdf5 error: Could not close HDF5 data space \n");
	MPI_Abort(comm, 1);
      }
      
      if(H5Sclose(filespace)) {
	printf("writehdf5 error: Could not close HDF5 file space \n");
	MPI_Abort(comm, 1);
      }
      if(H5Sclose(memspace)) {
	printf("writehdf5 error: Could not close HDF5 memory space \n");
	MPI_Abort(comm, 1);
      }
      if(h5_chunk)
	H5Pclose(chunk_pid);

    }

    if(H5Fclose(file_id) != 0)
      printf("writehdf5 error: Could not close HDF5 file \n");

    /* Create xdmf file for timestep */
    if(rank == 0)
      write_xdmf_xml(rel_fname, fname_xdmf, npoints);

}

void
write_xdmf_xml(char *fname, char *fname_xdmf, uint64_t npoints)
{
    FILE *xmf = 0;
 
    /*
     * Open the file and write the XML description of the mesh.
     */
 
    xmf = fopen(fname_xdmf, "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, "<Domain>\n");
    fprintf(xmf, "<Grid Name=\"Unstructured Mesh\">\n");
    fprintf(xmf, "<Topology TopologyType=\"Wedge\" NumberOfElements=\"%d\">\n", nelems_out[0]);
    fprintf(xmf, "<DataItem Dimensions=\"%d\" Format=\"HDF\">\n", nelems_out[0]*6);
    fprintf(xmf, "%s:/conns3\n",fname);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Topology>\n");
    fprintf(xmf, "<Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "<DataItem Name=\"X\" Dimensions=\"%d\" Format=\"HDF\">\n", npoints);
    fprintf(xmf, "%s:/grid points/x\n",fname);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Name=\"Y\" Dimensions=\"%d\" Format=\"HDF\">\n", npoints);
    fprintf(xmf, "%s:/grid points/y\n",fname);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "<DataItem Name=\"Z\" Dimensions=\"%d\" Format=\"HDF\">\n", npoints);
    fprintf(xmf, "%s:/grid points/z\n",fname);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Geometry>\n");
    fprintf(xmf, "<Attribute Name=\"Scalar\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "<DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", npoints);
    fprintf(xmf, "%s:/vars\n", fname);
    fprintf(xmf, "</DataItem>\n");
    fprintf(xmf, "</Attribute>\n");
    fprintf(xmf, "</Grid>\n");
    fprintf(xmf, "</Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

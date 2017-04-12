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

#include "hdf5.h"

static const int fnstrmax = 4095;

void
write_xdmf_xml_value(char *fname, char *fname_xdmf, char *varname, float deltax, float deltay, float deltaz, int ni, int nj, int nk);

void writehdf5i(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, int ni, int nj, int nk, int is, int ie, int js, int je,
		int ks, int ke, float deltax, float deltay, float deltaz, int nci, int ncj, int nck, float *data, hsize_t *h5_chunk, int hdf5_compress)
{
    char fname[fnstrmax+1];
    char fname_xdmf[fnstrmax+1];
    int timedigits = 4;
    MPI_Info info = MPI_INFO_NULL;

    hid_t file_id;
    hid_t plist_id;
    hid_t memspace;
    hid_t filespace;
    hid_t did;
    hsize_t start[3], count[3];
    hsize_t dims[3];
    herr_t err;
    hid_t chunk_pid;

    snprintf(fname, fnstrmax, "cart.%s_t%0*d.h5", varname, timedigits, tstep);
    snprintf(fname_xdmf, fnstrmax, "cart.%s_t%0*d.xmf", varname, timedigits, tstep);

    /* Create pvti file */
    if(rank == 0) {

      if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	fprintf(stderr, "writehdf5i error: could not create %s \n", fname);
	MPI_Abort(comm, 1);
      }

      /* Create dataset */
      dims[0] = (hsize_t)nk;
      dims[1] = (hsize_t)nj;
      dims[2] = (hsize_t)ni;
      filespace = H5Screate_simple(3, dims, NULL);
      
      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);

      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	if( H5Pset_fill_time(chunk_pid, H5D_FILL_TIME_NEVER) < 0 ) {
	  printf("writehdf5i error: Could not set fill time\n");
	  MPI_Abort(comm, 1);
	}
	H5Pset_chunk(chunk_pid, 3, h5_chunk);

	if(hdf5_compress == 1) {
	  /* Set ZLIB / DEFLATE Compression using compression level 6. */
	  if( H5Pset_deflate (chunk_pid, 6) < 0 ) {
	    printf("writehdf5i error: Could not set compression\n");
	    MPI_Abort(comm, 1);
	  }
	  
	  /* Uncomment these lines to set SZIP Compression
	     szip_options_mask = H5_SZIP_NN_OPTION_MASK;
	     szip_pixels_per_block = 16;
	     status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
	  */
	}
      }

      /* Create the dataset with default properties and close filespace. */
      did = H5Dcreate(file_id, varname, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
      H5Dclose(did);

      if(h5_chunk)
	H5Pclose(chunk_pid);

      H5Sclose(filespace);

      H5Fclose(file_id);

    } /* end rank==0 */

    MPI_Barrier(comm);

    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");

    /* Set up file access property list with parallel I/O access */
    if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
      printf("writehdf5i error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }
    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    if(H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info) < 0) {
      printf("writehdf5i error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    if( (file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id)) < 0) {
      fprintf(stderr, "writehdf5i error: could not open %s \n", fname);
      MPI_Abort(comm, 1);
    }

    if(H5Pclose(plist_id) < 0) {
      printf("writehdf5i error: Could not close property list \n");
      MPI_Abort(comm, 1);
    }

    start[0] = ks;
    start[1] = js;
    start[2] = is;

    count[0] = (hsize_t)nck;
    count[1] = (hsize_t)ncj;
    count[2] = (hsize_t)nci;

    if( (memspace = H5Screate_simple(3, count, NULL)) < 0) {
      printf("writehdf5i error: Could not create memory space \n");
      MPI_Abort(comm, 1);
    };

    if( (did = H5Dopen(file_id, varname, H5P_DEFAULT)) < 0) {
      printf("writehdf5i error: Could not open data space \n");
      MPI_Abort(comm, 1);
    };
      
    /* Select hyperslab in the file.*/
    filespace = H5Dget_space(did);
    if( H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL ) < 0) {
      printf("writehdf5i error: Could not select hyperslab \n");
      MPI_Abort(comm, 1);
    };

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
    if( err < 0) {
      fprintf(stderr, "writehdf5i error: could not write dataset %s \n", varname);
      MPI_Abort(comm, 1);
    }


    err = H5Dclose(did);

    err = H5Sclose(filespace);
    err = H5Sclose(memspace);
    
    if(H5Pclose(plist_id) < 0)
      printf("writehdf5i error: Could not close property list \n");

    if(H5Fclose(file_id) != 0)
      printf("writehdf5i error: Could not close HDF5 file \n");

    /* Create xdmf file for timestep */
    if(rank == 0) {
      write_xdmf_xml_value(fname, fname_xdmf, varname, deltax, deltay, deltaz, ni, nj, nk);
    }

}


void
write_xdmf_xml_value(char *fname, char *fname_xdmf, char *varname, float deltax, float deltay, float deltaz, int ni, int nj, int nk)
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
    fprintf(xmf, "  <Grid Name=\"Structured Grid\" GridType=\"Uniform\">\n");
    fprintf(xmf, "    <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n", ni, nj, nk );
    fprintf(xmf, "      <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n");
    fprintf(xmf, "        <DataItem Dimensions=\"3 \" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
    fprintf(xmf, "           0 0 0 \n");
    fprintf(xmf, "        </DataItem> \n");
    fprintf(xmf, "       <DataItem Dimensions=\"3 \" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">\n");
    fprintf(xmf, "            %f %f %f \n",deltax, deltay, deltaz);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "    </Geometry>\n");
    fprintf(xmf, "    <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", varname);
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", ni, nj, nk);
    fprintf(xmf, "        %s:/%s\n", fname, varname);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "    </Attribute>\n");
    fprintf(xmf, "  </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}


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

#include "hdf5.h"
#include "hdf5struct.h"

static const int fnstrmax = 4095;

void
write_xdmf_xml(char *fname, char *fname_xdmf, int num_xname, char **xname,
	       int ni, int nj, int nk,
	       float deltax, float deltay, float deltaz);

void writehdf5(const int num_varnames, char **varnames, MPI_Comm comm, int rank, int nprocs, int tstep, 
	       int is, int js, int ks,
               int ni, int nj, int nk, int cni, int cnj, int cnk, 
               float deltax, float deltay, float deltaz, 
               float *data, float *height, int *ola_mask, int *ol_mask)
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
    int j;
    herr_t err;
    
    snprintf(fname, fnstrmax, "struct_t%0*d.h5", timedigits, tstep);
    snprintf(fname_xdmf, fnstrmax, "struct_t%0*d.xmf", timedigits, tstep);

    if(rank == 0) {

      if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	fprintf(stderr, "writehdf5 error: could not create %s \n", fname);
	MPI_Abort(comm, 1);
      }

      dims[0] = nk;
      dims[1] = nj;
      dims[2] = ni;
      filespace = H5Screate_simple(3, dims, NULL);

      /* Create the dataset with default properties */
      
      for (j=0; j<num_varnames; j++) {

	if(strcmp(varnames[j],"data") == 0 || strcmp(varnames[j],"height") == 0) {
	  did = H5Dcreate(file_id, varnames[j], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	} else {
	  did = H5Dcreate(file_id, varnames[j], H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}

	H5Dclose(did);
      }

      /* Close filespace. */
      H5Sclose(filespace);
      
      /* Close the file */
      H5Fclose(file_id);
    }
    
    MPI_Barrier(comm);

    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");

    /* Set up file access property list with parallel I/O access */
    if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

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

    start[0] = (hsize_t)(ks);
    start[1] = (hsize_t)(js);
    start[2] = (hsize_t)(is);
    count[0] = (hsize_t)(cnk);
    count[1] = (hsize_t)(cnj);
    count[2] = (hsize_t)(cni);

    memspace = H5Screate_simple(3, count, NULL);

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    for (j=0; j<num_varnames; j++) {

      did = H5Dopen(file_id, varnames[j],H5P_DEFAULT);

      /* 
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */

      /* Select hyperslab in the file.*/
      filespace = H5Dget_space(did);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL );
      
      err = 0;
      if(strcmp(varnames[j],"data") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
      } else if(strcmp(varnames[j],"height") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, height);
      } else if(strcmp(varnames[j],"ola_mask") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_INT, memspace, filespace, plist_id, ola_mask);
      } else if(strcmp(varnames[j],"ol_mask") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_INT, memspace, filespace, plist_id, ol_mask);
      } else {
	printf("writehdf5 error: Unknown how to handle variable %s \n", varnames[j]);
	MPI_Abort(comm, 1);
      }

      if( err < 0) {
	fprintf(stderr, "writehdf5 error: could not write datset %s \n", varnames[j]);
	MPI_Abort(comm, 1);
      }

      err = H5Dclose(did);
    }

    if(H5Sclose(memspace)  != 0)
      printf("writehdf5 error: Could not close memory space \n");

    if(H5Pclose(plist_id) < 0)
      printf("writehdf5 error: Could not close property list \n");
      
    if(H5Fclose(file_id) != 0)
      printf("writehdf5 error: Could not close HDF5 file \n");

      /* Create xdmf file for timestep */
    if(rank == 0) {
      write_xdmf_xml(fname, fname_xdmf, num_varnames, varnames, ni, nj, nk, deltax, deltay, deltaz);
    }
}

void
write_xdmf_xml(char *fname, char *fname_xdmf, int num_xname, char **varnames, 
	       int ni, int nj, int nk,
	       float deltax, float deltay, float deltaz)
{
    FILE *xmf = 0;
    int j;
 
    /*
     * Open the file and write the XML description of the mesh.
     */
 
    xmf = fopen(fname_xdmf, "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"3.0\">\n");
    fprintf(xmf, " <Domain>\n\n");
    fprintf(xmf, "   <Grid Name =\"grid\" GridType=\"Uniform\">\n");
    fprintf(xmf, "    <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\">\n", nk, nj, ni);
    fprintf(xmf, "    </Topology>\n\n");
    fprintf(xmf, "    <Geometry Type=\"ORIGIN_DXDYDZ\">\n");
    fprintf(xmf, "        <!-- Origin -->\n");
    fprintf(xmf, "        <DataItem Format=\"XML\" Dimensions=\"3\">\n");
    fprintf(xmf, "                    0.0 0.0 0.0 \n");
    fprintf(xmf, "        </DataItem>\n");
    fprintf(xmf, "        <!-- DxDyDz -->\n");
    fprintf(xmf, "        <DataItem Format=\"XML\" Dimensions=\"3\">\n");
    fprintf(xmf, "                  %.6f %.6f %.6f \n", deltax, deltay, deltaz);
    fprintf(xmf, "        </DataItem>\n");
    fprintf(xmf, "    </Geometry>\n");
    for (j=0; j<num_xname; j++) {
      fprintf(xmf, "    <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n", varnames[j]);
      if(strcmp(varnames[j],"data") == 0 || strcmp(varnames[j],"height") == 0) {
	fprintf(xmf, "       <DataItem Dimensions=\"%d \" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", nk*nj*ni);
      } else {
	fprintf(xmf, "       <DataItem Dimensions=\"%d \" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">\n", nk*nj*ni);
      }
      fprintf(xmf, "        %s:/%s\n", fname, varnames[j]);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "    </Attribute>\n");
    }
    fprintf(xmf, "  </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

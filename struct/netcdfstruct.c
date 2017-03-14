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
#include "netcdfstruct.h"



/*! Save netcdf structured data summary to XML.
 *
 * Save metadata description of a mesh to XML.
 */
void write_xdmf_xml(char *fname, char *fname_xdmf, int num_xname, char **varnames,
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


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include "pdirs.h"
#include "pvti.h"

static const int fnstrmax = 4095;

void writepvti(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, int ni, int nj, int nk, int is, int ie, int js, int je,
               int ks, int ke, float deltax, float deltay, float deltaz, float *data)
{
    char dirname[fnstrmax+1];
    char fname[fnstrmax+1];
    char line[fnstrmax+1];
    int rankdigits = nprocs > 1 ? (int)(log10(nprocs-1)+1.5) : 1;
    int timedigits = 4;
    char typestr[] = "Float32";
    char endianstr[] = "LittleEndian";
    FILE *f;
    MPI_File mf;
    MPI_Status mstat;
    MPI_Info info = MPI_INFO_NULL;
    int r;
    int subset[6] = { is, ie, js, je, ks, ke };  /* Local point indices */
    int *rsubsets = NULL;   /* All point indices from each task */
    int ret;
    uint64_t ijkelems;
    int32_t ijkbytes32;

    /* Make directory for timestep */
    snprintf(dirname, fnstrmax, "%s.%s.%0*d.d", name, varname, timedigits, tstep);
    mkdir1task(dirname, comm);

    /* Gather subset sizes from all processes of communicator to rank 0 */
    if(rank == 0) {
        rsubsets = (int *) malloc(nprocs*6*sizeof(int));
    }
    MPI_Gather(subset, 6, MPI_INT, rsubsets, 6, MPI_INT, 0, comm);

    /* Create pvti file */
    if(rank == 0) {
        snprintf(fname, fnstrmax, "%s.%s.%0*d.pvti", name, varname, timedigits, tstep);
        if( ! (f = fopen(fname, "w")) ) {
            fprintf(stderr, "writepvti error: Could not create .pvti file.\n");
            MPI_Abort(comm, 1);
        }
        fprintf(f, "<?xml version=\"1.0\"?>\n"
                   "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"%s\">\n",
                   endianstr);
        fprintf(f, "  <PImageData WholeExtent=\"0 %d 0 %d 0 %d\" GhostLevel=\"0\" "
                   "Origin=\"0 0 0\" Spacing=\"%f %f %f\">\n", ni-1, nj-1, nk-1,
                   deltax, deltay, deltaz);
        fprintf(f, "    <PPointData Scalars=\"%s\">\n"
                   "      <PDataArray type=\"%s\" Name=\"%s\"/>\n"
                   "    </PPointData>\n", varname, typestr, varname);
        /* Write info for each block */
        for(r = 0; r < nprocs; ++r) {
            int *rsub = rsubsets + r*6;    /* Remote subset indices */
            snprintf(fname, fnstrmax, "%s.%s.%0*d.d/%0*d.vti", name, varname, timedigits,
                     tstep, rankdigits, r);
            fprintf(f, "    <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\"/>\n",
                    rsub[0], rsub[1], rsub[2], rsub[3], rsub[4], rsub[5], fname);
        }
        /* Finish the .pvti file */
        fprintf(f, "  </PImageData>\n</VTKFile>\n");
        fclose(f);
        free(rsubsets);
    } /* end rank==0 */

    chkdir1task(dirname, comm);

    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");

    /* Write .vti files */
    snprintf(fname, fnstrmax, "%s.%s.%0*d.d/%0*d.vti", name, varname, timedigits,
             tstep, rankdigits, rank);
    ret = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                        info, &mf);
    if(ret) {
        fprintf(stderr, "writepvti error: could not open %s\n", fname);
        MPI_Abort(comm, 1);
    }
    snprintf(line, fnstrmax, "<?xml version=\"1.0\"?>\n"
             "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"%s\">\n",
             endianstr);
    MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
    snprintf(line, fnstrmax, "  <ImageData WholeExtent=\"%d %d %d %d %d %d\" "
             "Origin=\"0 0 0\" Spacing=\"%f %f %f\">\n    "
             "<Piece Extent=\"%d %d %d %d %d %d\">\n", is, ie, js, je, ks, ke,
             deltax, deltay, deltaz, is, ie, js, je, ks, ke);
    MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
    snprintf(line, fnstrmax, "      <PointData Scalars=\"%s\">\n", varname);
    MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
    snprintf(line, fnstrmax, "        <DataArray type=\"%s\" Name=\"%s\" "
             "format=\"appended\" offset=\"0\"/>\n", typestr, varname);
    MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
    snprintf(line, fnstrmax, "      </PointData>\n"
             "      <CellData></CellData>\n"
             "    </Piece>\n"
             "  </ImageData>\n"
             "  <AppendedData encoding=\"raw\">\n"
             "   _");
    MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
    ijkelems = (ie-is+1) * (je-js+1) * (ke-ks+1);
    ijkbytes32 = (int32_t)(ijkelems * sizeof(float));
    MPI_File_write(mf, &ijkbytes32, 1, MPI_INT, &mstat);
    MPI_File_write(mf, data, ijkelems, MPI_FLOAT, &mstat);
    snprintf(line, fnstrmax, "\n  </AppendedData>\n</VTKFile>\n"); 
    MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
    MPI_File_close(&mf);
}




#include <stdio.h>
#include <math.h>

#include <pdirs.h>
#include "przm.h"

static const int fnstrmax = 4095;

void writeprzm(char *name, MPI_Comm comm, int tstep, uint64_t npoints,
               float *xpts, float *ypts, float *zpts, uint64_t nelems3, uint64_t *conns3,
               char *varname, float *data)
{
    char dirname[fnstrmax+1];
    char fname[fnstrmax+1];
    char line[fnstrmax+1];
    int rank, nprocs;
    int rankdigits;
    int timedigits = 4;
    FILE *f;
    MPI_File mf;
    MPI_Status mstat;
    MPI_Info info = MPI_INFO_NULL;
    int r;
    uint64_t *rntris;   /* All triangle counts from each task */
    int ret;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    rankdigits = nprocs > 1 ? (int)(log10(nprocs-1)+1.5) : 1;

    /* Make dir for all output and subdir for timestep */
    snprintf(dirname, fnstrmax, "%s.przm", name);
    mkdir1task(dirname, comm);
    snprintf(dirname, fnstrmax, "%s.przm/t%0*d.d", name, timedigits, tstep);
    mkdir1task(dirname, comm);

    /* Placeholder to create a metadata file, if we decide we need one for reading */

    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");    

    chkdir1task(dirname, comm);

    /* Write files per rank */

    snprintf(fname, fnstrmax, "%s/r%0*d.dat", dirname, rankdigits, rank);
    ret = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                        info, &mf);
    if(ret) {
        fprintf(stderr, "writeprzm error: could not open %s\n", fname);
        MPI_Abort(comm, 1);
    }

    MPI_File_write(mf, &npoints, 1, MPI_UNSIGNED_LONG_LONG, &mstat);

    /* Optional grid points */
    if(xpts && ypts && zpts) {
        uint32_t hasgrid = 1;
        MPI_File_write(mf, &hasgrid, 1, MPI_UNSIGNED, &mstat);
        MPI_File_write(mf, xpts, npoints, MPI_FLOAT, &mstat);
        MPI_File_write(mf, ypts, npoints, MPI_FLOAT, &mstat);
        MPI_File_write(mf, zpts, npoints, MPI_FLOAT, &mstat);
    } else {
        uint32_t hasgrid = 0;
        MPI_File_write(mf, &hasgrid, 1, MPI_UNSIGNED, &mstat);
    }

    /* Optional grid connections, writes a 64-bit 0 if no new connections */
    if(conns3 && nelems3) {
        MPI_File_write(mf, &nelems3, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
        MPI_File_write(mf, conns3, nelems3*6, MPI_UNSIGNED_LONG_LONG, &mstat);
    } else {
        uint64_t hasconn = 0;
        MPI_File_write(mf, &hasconn, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
    }

    /* Optional variable data, starting with number of variables */
    if(data && varname) {
        uint32_t nvars = 1;
        MPI_File_write(mf, &nvars, 1, MPI_UNSIGNED, &mstat);
        MPI_File_write(mf, data, npoints, MPI_FLOAT, &mstat);
    } else {
        uint32_t nvars = 0;
        MPI_File_write(mf, &nvars, 1, MPI_UNSIGNED, &mstat);
    }

    MPI_File_close(&mf);
}


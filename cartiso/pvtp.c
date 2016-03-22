
#include <stdio.h>
#include <errno.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "pvtp.h"

static const int fnstrmax = 4095;

static void mkdir1task(char *dirname, MPI_Comm comm, int rank)
{
    int ret;

    if(rank == 0) {
        ret = mkdir(dirname, 0777);
        if( ret && errno != EEXIST ) {  /* An error for anything but pre-existing */
            fprintf(stderr, "pvtk error: Could not create directory.\n");
            MPI_Abort(comm, 1);
        }
    }
}

static void chkdir1task(char *dirname, MPI_Comm comm, int rank, int nprocs)
{
    int ret;
    struct stat fsbuf;

    /* Have the last rank ensure that the directory is available */
    if(rank == nprocs-1) { 
        do {
            int retry = 0;
            usleep(200000);
            ret = stat(dirname, &fsbuf);
            if(ret && retry > 50) {   /* 50 retries * 200000 us = 10 s timeout */
                fprintf(stderr, "pvtk error: timeout waiting for directory %s\n",
                        dirname);
                MPI_Abort(comm, 1);
            }
        } while(ret);
    }
    MPI_Barrier(comm);   /* Everyone is safe to use directory after barrier */
}


void writepvtp(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, uint64_t ntris, float *points, float *norms, 
               float *xvals, char *xname)
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
    int r;
    uint64_t *rntris;   /* All triangle counts from each task */
    int ret;
    
    /* Make directory for timestep */
    snprintf(dirname, fnstrmax, "%s.%s.%0*d.d", name, varname, timedigits, tstep);
    mkdir1task(dirname, comm, rank);

    /* Gather tri counts, in case some are zero, to leave them out */
    if(rank == 0)
        rntris = (uint64_t *) malloc(nprocs*sizeof(uint64_t));
    MPI_Gather(&ntris, 1, MPI_UNSIGNED_LONG_LONG, rntris, 1, MPI_LONG_LONG, 0, comm);

    /* Create pvtp file */
    if(rank == 0) {
        snprintf(fname, fnstrmax, "%s.%s.%0*d.pvtp", name, varname, timedigits, tstep);
        if( ! (f = fopen(fname, "w")) ) {
            fprintf(stderr, "writepvtp error: Could not create .pvtp file.\n");
            MPI_Abort(comm, 1);
        }
        fprintf(f, "<?xml version=\"1.0\"?>\n"
                   "<VTKFile type=\"PPolyData\" version=\"1.0\" byte_order=\"%s\" "
                   "header_type=\"UInt64\">\n",
                   endianstr);
        fprintf(f, "  <PPolyData GhostLevel=\"0\">\n");
        fprintf(f, "    <PPointData Normals=\"Normals\">\n"
                   "      <PDataArray type=\"%s\" Name=\"Normals\" NumberOfComponents=\"3\"/>\n",
                   typestr);
        if(xvals) {
            fprintf(f, "      <PDataArray type=\"%s\" Name=\"%s\" NumberOfComponents=\"1\"/>\n",
                    typestr, xname);
        }
        fprintf(f, "    </PPointData>\n");
        fprintf(f, "    <PPoints>\n"
                   "      <PDataArray type=\"%s\" Name=\"Points\" NumberOfComponents=\"3\"/>\n"
                   "    </PPoints>\n", typestr);
        /* Write info for each block */
        for(r = 0; r < nprocs; ++r) {
            if(rntris[r] > 0) {   /* But only for blocks that have triangles */
                snprintf(fname, fnstrmax, "%s.%s.%0*d.d/%0*d.vtp", name, varname, timedigits,
                         tstep, rankdigits, r);
                fprintf(f, "    <Piece Source=\"%s\"/>\n", fname);
            }
        }
        /* Finish the .pvtp file */
        fprintf(f, "  </PPolyData>\n</VTKFile>\n");
        fclose(f);
        free(rntris);
    } /* end rank==0 */

    chkdir1task(dirname, comm, rank, nprocs);

    /* Write .vtp files, if we have polys */
    if(ntris > 0) {
        uint64_t i, offsets = 0;   /* Offsets into the binary portion for each field */
        uint64_t *temparr;
        snprintf(fname, fnstrmax, "%s.%s.%0*d.d/%0*d.vtp", name, varname, timedigits,
                 tstep, rankdigits, rank);
        ret = MPI_File_open(MPI_COMM_SELF, fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,
                            MPI_INFO_NULL, &mf);
        if(ret) {
            fprintf(stderr, "writepvtp error: could not open %s\n", fname);
            MPI_Abort(comm, 1);
        }
        snprintf(line, fnstrmax, "<?xml version=\"1.0\"?>\n"
                 "<VTKFile type=\"PolyData\" version=\"2.0\" byte_order=\"%s\" "
                 "header_type=\"UInt64\">\n"
                 "  <PolyData>\n", endianstr);
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        snprintf(line, fnstrmax, "    <Piece NumberOfPoints=\"%llu\" NumberOfVerts=\"0\" "
                 "NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"%llu\">\n",
                 ntris*3, ntris);
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        snprintf(line, fnstrmax, "      <PointData Normals=\"Normals\">\n"
                 "        <DataArray type=\"%s\" Name=\"Normals\" NumberOfComponents=\"3\" "
                 "format=\"appended\" offset=\"0\"/>\n", typestr);
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        offsets += ntris*9*sizeof(float)+sizeof(uint64_t);   /* Add size of normals for points */
        if(xvals) {
            snprintf(line, fnstrmax, "        <DataArray type=\"%s\" Name=\"%s\" "
                     "NumberOfComponents=\"1\" "
                     "format=\"appended\" offset=\"%llu\">\n", typestr, xname, offsets);
            MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
            offsets += ntris*3*sizeof(float)+sizeof(uint64_t);   /* Add size of xdata array */
        }
        snprintf(line, fnstrmax, "      </PointData>\n"
                 "      <Points>\n"
                 "        <DataArray type=\"%s\" Name=\"Points\" NumberOfComponents=\"3\" "
                 "format=\"appended\" offset=\"%llu\"/>\n"
                 "      </Points>\n", typestr, offsets);
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        offsets += ntris*9*sizeof(float)+sizeof(uint64_t);   /* Add size of points */
        snprintf(line, fnstrmax, "      <Polys>\n"
                 "        <DataArray type=\"UInt64\" Name=\"connectivity\" "
                 "format=\"appended\" offset=\"%llu\"/>\n", offsets);
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        offsets += ntris*3*sizeof(uint64_t)+sizeof(uint64_t);   /* Add size of connections */
        snprintf(line, fnstrmax, "        <DataArray type=\"UInt64\" Name=\"offsets\" "
                 "format=\"appended\" offset=\"%llu\"/>\n"
                 "      </Polys>\n    </Piece>\n  </PolyData>\n"
                 "  <AppendedData encoding=\"raw\">\n"
                 "   _", offsets);
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        
        /* Normals header & array */
        offsets = ntris*9*sizeof(float);
        MPI_File_write(mf, &offsets, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
        MPI_File_write(mf, norms, ntris*9, MPI_FLOAT, &mstat);
        /* xvals header & array */
        if(xvals) {
            offsets = ntris*3*sizeof(float);
            MPI_File_write(mf, &offsets, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
            MPI_File_write(mf, xvals, ntris*3, MPI_FLOAT, &mstat);
        }
        /* Points header & array */
        offsets = ntris*9*sizeof(float);
        MPI_File_write(mf, &offsets, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
        MPI_File_write(mf, points, ntris*9, MPI_FLOAT, &mstat);
        /* Create connections & write */
        temparr = (uint64_t *) malloc(ntris*3*sizeof(uint64_t));
        for(i = 0; i < ntris*3; i++)
            temparr[i] = i;
        offsets = ntris*3*sizeof(uint64_t);
        MPI_File_write(mf, &offsets, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
        MPI_File_write(mf, temparr, ntris*3, MPI_UNSIGNED_LONG_LONG, &mstat);
        /* Create offsets & write */
        for(i = 1; i <= ntris; i++)
            temparr[i] = i*3;
        offsets = ntris*sizeof(uint64_t);
        MPI_File_write(mf, &offsets, 1, MPI_UNSIGNED_LONG_LONG, &mstat);
        MPI_File_write(mf, temparr, ntris, MPI_UNSIGNED_LONG_LONG, &mstat);
        /* Finish and close */
        free(temparr);
        snprintf(line, fnstrmax, "\n  </AppendedData>\n</VTKFile>\n"); 
        MPI_File_write(mf, line, strlen(line), MPI_CHAR, &mstat);
        MPI_File_close(&mf);
    }

}


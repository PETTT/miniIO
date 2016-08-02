
#include <stdio.h>
#include <stdlib.h>

#include "adiosstruct.h"
#include "pdirs.h"

static const int fnstrmax = 4095;

void adiosstruct_init(struct adiosstructinfo *nfo, char *method, char *name, MPI_Comm comm, int rank, int nprocs, int tsteps) {
  nfo->name = name;
  nfo->comm = comm;
  nfo->rank = rank;
  nfo->nprocs = nprocs;
  nfo->tsteps = tsteps;
    
  nfo->nvars = 0;
  nfo->maxvars = 1000;
  nfo->varnames = (char **) malloc(nfo->maxvars * sizeof(char *));
  nfo->datas = (float **) malloc(nfo->maxvars * sizeof(float *));
  nfo->bufallocsize = 0;

  /* Set up ADIOS */
  adios_init_noxml(comm);
  adios_declare_group(&nfo->gid, nfo->name, "", adios_flag_no);
  adios_select_method(nfo->gid, method, "", "");

  /* Define output variables */
  adios_define_var(nfo->gid, "rank", "", adios_integer, "", "", "");
  adios_define_var(nfo->gid, "tstep", "", adios_integer, "", "", "");
  adios_define_var(nfo->gid, "cnpoints", "", adios_unsigned_long, "", "", "");
  adios_define_var(nfo->gid, "npoints", "", adios_unsigned_long, "", "", "");
  adios_define_var(nfo->gid, "cstart", "", adios_unsigned_long, "", "", "");

  adios_define_var(nfo->gid, "ola_mask", "", adios_unsigned_long,  "cnpoints", "npoints", "cstart");
  
  
}

void adiosstruct_addxvar(struct adiosstructinfo *nfo, char *varname, float *data) {

    if(nfo->nvars >= nfo->maxvars)
        return;   /* Just ignore too many variables, for now */
    nfo->varnames[nfo->nvars] = varname;
    nfo->datas[nfo->nvars] = data;
    nfo->nvars++;
    adios_define_var(nfo->gid, varname, "", adios_real, "cnpoints", "npoints", "cstart");
}

void adiosstruct_write(struct adiosstructinfo *nfo, int tstep, uint64_t cnpoints, uint64_t npoints, uint64_t cstart, uint64_t *mask) {
    char fname[fnstrmax+1];
    int timedigits = 4;
    uint64_t groupsize, totalsize;
    int64_t handle;
    int ret, i;
    int bufneeded;

    /* Set sizes */
    groupsize = sizeof(int) * 2 /*rank-tstep*/ +
                sizeof(uint64_t) * 3 /*npoints-cstart*/ +
                sizeof(uint64_t) * cnpoints /*ola_mask*/ +
                sizeof(float) * cnpoints * nfo->nvars ;
    
    /* adios_allocate_buffer (ADIOS_BUFFER_ALLOC_NOW, 10); */
    /* Allocate buffer large enough for all data to write, if not done already */
    bufneeded = (int)(groupsize/(1024*1024));
    bufneeded += bufneeded/10 + 10;   /* Add an extra 10% & 10MB to be sure */
    if(nfo->bufallocsize < bufneeded) {
        adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, bufneeded);
        nfo->bufallocsize = bufneeded;
    }

    /* Set filename */
    snprintf(fname, fnstrmax, "%s.%0*d.bp", nfo->name, timedigits, tstep);

    /* Open & Write */
    ret = adios_open(&handle, nfo->name, fname, "w", nfo->comm);
    if(ret) {
        fprintf(stderr, "Error opening ADIOS file: %s\n", fname);
        return;
    }
    adios_group_size(handle, groupsize, &totalsize);

    adios_write(handle, "rank", &nfo->rank);
    adios_write(handle, "tstep", &tstep);
    adios_write(handle, "cnpoints", &cnpoints);
    adios_write(handle, "npoints", &npoints);
    adios_write(handle, "cstart", &cstart);
    adios_write(handle, "ola_mask", mask);
    for(i = 0; i < nfo->nvars; i++) {
        adios_write(handle, nfo->varnames[i], nfo->datas[i]);
    }
    adios_close(handle);
}

void adiosstruct_finalize(struct adiosstructinfo *nfo) {
    free(nfo->varnames);
    free(nfo->datas);
    adios_finalize(nfo->rank);
}


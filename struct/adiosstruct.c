
#include <stdio.h>
#include <stdlib.h>

#include "adiosstruct.h"
#include "pdirs.h"

static const int fnstrmax = 4095;

void adiosstruct_init(struct adiosstructinfo *nfo, char *method,
               char *name, MPI_Comm comm, int rank, int nprocs,
               int tsteps, int ni, int nj, int nk, int is, int cni, int js, int cnj,
               int ks, int cnk, float deltax, float deltay, float deltaz) {
  
  nfo->name = name;
  nfo->comm = comm;
  nfo->rank = rank;
  nfo->nprocs = nprocs;
  nfo->tsteps = tsteps;
  nfo->ni = ni;
  nfo->nj = nj;
  nfo->nk = nk;
  nfo->is = is;
  nfo->cni = cni;
  nfo->js = js;
  nfo->cnj = cnj;
  nfo->ks = ks;
  nfo->cnk = cnk;
  nfo->deltax = deltax;
  nfo->deltay = deltay;
  nfo->deltaz = deltaz;
  
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
  /* Define output variables */
    adios_define_var(nfo->gid, "rank", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "tstep", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "ni", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "nj", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "nk", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "is", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "cni", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "js", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "cnj", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "ks", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "cnk", "", adios_integer, "", "", "");
    adios_define_var(nfo->gid, "deltax", "", adios_real, "", "", "");
    adios_define_var(nfo->gid, "deltay", "", adios_real, "", "", "");
    adios_define_var(nfo->gid, "deltaz", "", adios_real, "", "", "");

  adios_define_var(nfo->gid, "ola_mask", "", adios_integer,   "cnk,cnj,cni",
                     "nk,nj,ni", "ks,js,is");
  
  
}

void adiosstruct_addxvar(struct adiosstructinfo *nfo, char *varname, float *data) {

    if(nfo->nvars >= nfo->maxvars)
        return;   /* Just ignore too many variables, for now */
    nfo->varnames[nfo->nvars] = varname;
    nfo->datas[nfo->nvars] = data;
    nfo->nvars++;
    adios_define_var(nfo->gid, varname, "", adios_real, "cnk,cnj,cni",
                     "nk,nj,ni", "ks,js,is");
}

void adiosstruct_write(struct adiosstructinfo *nfo, int tstep, int *mask) {
    char fname[fnstrmax+1];
    int timedigits = 4;
    uint64_t ijkelems, groupsize, totalsize;
    int64_t handle;
    int ret, i;
    int bufneeded;

    ijkelems = (uint64_t) nfo->cni * nfo->cnj * nfo->cnk;
    groupsize = sizeof(int) /*rank*/ + sizeof(int) /*tstep*/ + 
                sizeof(int)*3 /*ni,nj,nk*/ + sizeof(int)*6 /*is-cnk*/ +
                sizeof(float)*3 /*deltax,y,z*/ +
                sizeof(int) * ijkelems /* mask */ +
                sizeof(float) * ijkelems * nfo->nvars;
    
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
    adios_write(handle, "ni", &nfo->ni);
    adios_write(handle, "nj", &nfo->nj);
    adios_write(handle, "nk", &nfo->nk);
    adios_write(handle, "is", &nfo->is);
    adios_write(handle, "cni", &nfo->cni);
    adios_write(handle, "js", &nfo->js);
    adios_write(handle, "cnj", &nfo->cnj);
    adios_write(handle, "ks", &nfo->ks);
    adios_write(handle, "cnk", &nfo->cnk);
    adios_write(handle, "deltax", &nfo->deltax);
    adios_write(handle, "deltay", &nfo->deltay);
    adios_write(handle, "deltaz", &nfo->deltaz);
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


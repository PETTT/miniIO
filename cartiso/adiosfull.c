/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */
 
#include <stdio.h>
#include <stdlib.h>

#include "adiosfull.h"
#include "pdirs.h"

static const int fnstrmax = 4095;

void adiosfull_init(struct adiosfullinfo *nfo, char *method,
               char *name, MPI_Comm comm, int rank, int nprocs,
               int tsteps, int ni, int nj, int nk, int is, int cni, int js, int cnj,
               int ks, int cnk, float deltax, float deltay, float deltaz, char *adiosopts)
{
    char dirname[fnstrmax+1];
    static char emptystr[] = "";

    /* Make directory for output collection, all timesteps */
    /*snprintf(dirname, fnstrmax, "%s_adios.d/", name);
    mkdir1task(dirname, comm);*/

    /* Set up struct */
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

    if(!adiosopts)  adiosopts = emptystr;
    
    /* Set up ADIOS */ 
    adios_init_noxml(comm);
    adios_declare_group(&nfo->gid, nfo->name, "", adios_flag_no);
    /*chkdir1task(dirname, comm);*/
    adios_select_method(nfo->gid, method, adiosopts, /*dirname*/"");

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
}

void adiosfull_addvar(struct adiosfullinfo *nfo, char *varname, float *data)
{
    if(nfo->nvars >= nfo->maxvars)
        return;   /* Just ignore too many variables, for now */
    nfo->varnames[nfo->nvars] = varname;
    nfo->datas[nfo->nvars] = data;
    nfo->nvars++;
    adios_define_var(nfo->gid, varname, "", adios_real, "cni,cnj,cnk",
                     "ni,nj,nk", "is,js,ks");
}

void adiosfull_write(struct adiosfullinfo *nfo, int tstep)
{
    char fname[fnstrmax+1];
    int timedigits = 4;
    uint64_t ijkelems, groupsize, totalsize;
    int64_t handle;
    int ret, i;
    int bufneeded;

    /* Set sizes */
    ijkelems = (uint64_t)nfo->cni * nfo->cnj * nfo->cnk;
    groupsize = sizeof(int) /*rank*/ + sizeof(int) /*tstep*/ + 
                sizeof(int)*3 /*ni,nj,nk*/ + sizeof(int)*6 /*is-cnk*/ +
                sizeof(float)*3 /*deltax,y,z*/ +
                sizeof(float) * ijkelems * nfo->nvars + 1024;

    /* Allocate buffer large enough for all data to write, if not done already */
    bufneeded = (int)(groupsize/(1024*1024));
    bufneeded += bufneeded/8 + 2;   /* Add an extra 12.5% & 2MB to be sure */
    if(nfo->bufallocsize < bufneeded) {
#       if ADIOS_VERSION_GE(1,10,0)
        adios_set_max_buffer_size(bufneeded);
#       else
        adios_allocate_buffer(ADIOS_BUFFER_ALLOC_NOW, bufneeded);
#       endif
        nfo->bufallocsize = bufneeded;
    }

    /* Set filename */
    snprintf(fname, fnstrmax, "%s_%0*d", nfo->name, timedigits, tstep);

    /* Open & Write */
    ret = adios_open(&handle, nfo->name, fname, "w", nfo->comm);
    if(ret) {
        fprintf(stderr, "Error opening ADIOS file: %s\n", fname);
        return;
    }
#   if ADIOS_VERSION_LE(1,9,0)
    adios_group_size(handle, groupsize, &totalsize);
#   endif

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
    for(i = 0; i < nfo->nvars; i++)
        adios_write(handle, nfo->varnames[i], nfo->datas[i]);

    adios_close(handle);
}

void adiosfull_finalize(struct adiosfullinfo *nfo)
{
    free(nfo->varnames);
    free(nfo->datas);
    adios_finalize(nfo->rank);
}


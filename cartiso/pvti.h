/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#define __STDC_FORMAT_MACROS
#include <inttypes.h>

void writepvti(char *name, char *varname, MPI_Comm comm, int rank, int nprocs, 
               int tstep, int ni, int nj, int nk, int is, int ie, int js, int je,
               int ks, int ke, float deltax, float deltay, float deltaz, float *data);

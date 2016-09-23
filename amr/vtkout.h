/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */


void writevtk(char *name, MPI_Comm comm, int rank, int nprocs, int tstep,
              uint64_t npoints, uint64_t ncubes, float *points, float *xvals, char *xname, int debug);

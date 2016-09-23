/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include "vtkout.h"

static const int fnstrmax = 4095;
static const int cubeVertices = 8;
static const int vertexDims= 3;
static const int cubeVertexDims = cubeVertices * vertexDims;

void writevtk(char *name, MPI_Comm comm, int rank, int nprocs, int tstep,
	      uint64_t npoints, uint64_t ncubes, float *points, float *xvals, char *xname, int debug) {

  char fname[fnstrmax+1];
  int timedigits = 4;
  FILE *f;
  int  i;
  uint64_t allcubes;        /* total cube count: sum over all ranks */
  uint64_t *rncubes=NULL;   /* cube count for currrent rank */
  
  int nallpoints=0;
  int  nallxvals=0;
  float *allpoints=NULL;
  float *allxvals=NULL;

  int cpointcnts;
  int cpointoffsets;
  int *pointcnts=NULL;
  int *pointoffsets=NULL;
  
  int cxvalcnts;
  int cxvaloffsets;
  int *xvalcnts=NULL;
  int *xvaloffsets=NULL;

  if (debug) {
    printf("Debug 0.0 (rank=%d): Entered writevtk routine\n", rank);
  }
  
  pointcnts = (int *) malloc(nprocs*sizeof(int));
  pointoffsets = (int *) malloc(nprocs*sizeof(int));
  xvalcnts = (int *) malloc(nprocs*sizeof(int));
  xvaloffsets = (int *) malloc(nprocs*sizeof(int));

  if(rank == 0) {
    rncubes = (uint64_t *) malloc(nprocs*sizeof(uint64_t));
  }

  MPI_Gather(&ncubes, 1, MPI_UNSIGNED_LONG_LONG, rncubes, 1, MPI_LONG_LONG, 0, comm);
  
  cpointcnts = ncubes * cubeVertexDims;
  MPI_Allgather(&cpointcnts, 1, MPI_INT, pointcnts, 1, MPI_INT, comm);

  cxvalcnts = ncubes * cubeVertices;
  MPI_Allgather(&cxvalcnts, 1, MPI_INT, xvalcnts, 1, MPI_INT, comm);

  cpointoffsets = 0;
  cxvaloffsets = 0;
  for(i = 1; i < rank+1; i++) {
    cpointoffsets += pointcnts[i-1];
    cxvaloffsets += xvalcnts[i-1];
  }
 
  MPI_Allgather(&cpointoffsets, 1, MPI_INT, pointoffsets, 1, MPI_INT, comm);
  MPI_Allgather(&cxvaloffsets, 1, MPI_INT, xvaloffsets, 1, MPI_INT, comm);
  
  if(rank == 0) {
    allcubes=0;
    for(i = 0; i < nprocs; i++) {
      allcubes += rncubes[i];
    }
    nallpoints = allcubes*cubeVertexDims;
    nallxvals  = allcubes*cubeVertices;
    
    allpoints = (float *) malloc(nallpoints*sizeof(float));
    allxvals = (float *) malloc(nallxvals*sizeof(float));
    
    for(i = 0; i <allcubes*cubeVertexDims; i++) {
      allpoints[i] = 0;
    }
    
    for(i = 0; i <allcubes*cubeVertices; i++) {
      allxvals[i] = 0;
    }

    if (debug) {
      printf("Debug 1.0 (rank=%d): nallpoints=%d  nallxvals=%d  allcubes:%llu\n", rank, nallpoints, nallxvals, allcubes);
      for(i = 0; i < nprocs; i++) {
	printf("Debug 2.0 (rank=%d): rncubes:%llu timestep=%d\n", i, rncubes[i], tstep);
      }
    }
  }

  if (debug) {
    for(i = 0; i < nprocs; i++) {
      printf("Debug 3.0 (rank=%d): i=%d  pointcnts=%d pointoffsets=%d\n", rank, i, pointcnts[i], pointoffsets[i]);
    }
  }
    
  MPI_Gatherv(points, pointcnts[rank], MPI_FLOAT, allpoints, pointcnts, pointoffsets, MPI_FLOAT, 0, comm);
  MPI_Gatherv(xvals, xvalcnts[rank], MPI_FLOAT, allxvals, xvalcnts, xvaloffsets, MPI_FLOAT, 0, comm);

  /* Create VTK file */
  if(rank == 0) {
  
    snprintf(fname, fnstrmax, "%s.%s.%0*d.vtk", name, xname, timedigits, tstep);

    if( ! (f = fopen(fname, "w")) ) {
      fprintf(stderr, "writepvtp error: Could not create .vtk file.\n");
      MPI_Abort(comm, 1);
    }
    
    fprintf(f, "# vtk DataFile Version 2.0\n");
    fprintf(f, "vtk output\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
    fprintf(f, "POINTS %llu float\n", allcubes*cubeVertices);
    for(i = 0; i <allcubes*cubeVertexDims ; i+=vertexDims) {
      fprintf(f, "%f %f %f\n", allpoints[i], allpoints[i+1], allpoints[i+2]);
    }
    
    fprintf(f, "\nCELLS %llu %llu\n", allcubes, (allcubes*cubeVertices) + allcubes);
    for(i = 0; i <allcubes*cubeVertices ; i+=cubeVertices) {
      fprintf(f, "8 %d %d %d %d %d %d %d %d\n", i, i+1, i+2, i+3, i+4, i+5, i+6, i+7);
    }

    fprintf(f, "\nCELL_TYPES %llu\n", allcubes);
    for(i = 0; i <allcubes ; i++) {
      fprintf(f, "11\n");
    }

    fprintf(f, "\nPOINT_DATA %llu\n", allcubes*cubeVertices);
    fprintf(f, "SCALARS %s float 1\n", xname);
    fprintf(f, "LOOKUP_TABLE default\n");
    for(i = 0; i <allcubes*cubeVertices ; i+=cubeVertices) {
      fprintf(f, "%f %f %f %f %f %f %f %f\n", allxvals[i], allxvals[i+1], allxvals[i+2], allxvals[i+3], allxvals[i+4], allxvals[i+5], allxvals[i+6], allxvals[i+7]);
    }

    fclose(f);
    free(rncubes);
    free(allpoints);
    free(allxvals);
     
  } /* end rank==0 */
  
  free(pointcnts);
  free(pointoffsets);
  free(xvalcnts);
  free(xvaloffsets);
  
  if (debug) {
    printf("Debug 11 (rank=%d): Exiting  writevtk routine\n", rank);
  }
  
}

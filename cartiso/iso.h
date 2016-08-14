/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <stdint.h>

struct isoinfo {
    float x0, y0, z0;    /* Starting coordinates */
    float xd, yd, zd;    /* Coordinate deltas */
    int xdim, ydim, zdim;   /* Dimensions of local data */
    uint64_t ntris;        /* Number of triangles */
    float *points;     /* Points of triangles */
    float *norms;     /* Normals of triangles */
    float *xvals;     /* Extra values on triangle points */
};

/* Initialize info for isosurface
 *   Must be done before the first call to isosurf.
 *   Does not need to be repeated for subsequent calls to isosurf,
 *   unless the coordinate, deltas, or dimensions change, but in
 *   that case, isofree must be called first */ 

void isoinit(struct isoinfo *nfo, float x0, float y0, float z0,
        float xd, float yd, float zd, int xdim, int ydim, int zdim,
        int numxarrays);

/* Free up resources from isosurface */

void isofree(struct isoinfo *nfo);

/* Generate an isosurface */

void isosurf(struct isoinfo *nfo,   /* Isosurface info */
             float thresh,          /* Isosurface threshold */
             const float *data,    /* Pointer to data, i,j,k order */
             const float *xdata);    /* Extra data */


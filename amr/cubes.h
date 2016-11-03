/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

typedef struct cubeInfo {
  int debug;
  uint64_t ncubes;       /* Number of cubes */
  uint64_t npoints;      /* Number of cube points */
  float *points;         /* Points of cube */
  float *data;           /* Values on points */
}cubeInfo;


/* holds an entry in the stack (in this case cube information) */
typedef struct cubeItem {
  int r_id;
  int level;
  float x;
  float y;
  float z;
  float dx;
  float dy;
  float dz;
} cubeItem;

/* A stack that holds and array of stack items (in this case an array of cube information) */
typedef struct {
  cubeItem *cubes;
  int top;
  int size;
  int maxsize;
  int debug;
} stack;


void cubesinit(cubeInfo *nfo, int task, int levels, int debug);

void cubesfree(cubeInfo *nfo);

void refine(cubeInfo *nfo, int t, int rpId, float thres, int level_start, float x_start, float y_start,
	    float z_start, float dx_start, float dy_start, float dz_start, struct osn_context *osn, int maxLevel,
	    double noisespacefreq, double noisetimefreq);

void cubeprint(cubeInfo *nfo);

void stack_new(stack *astack, int maxsize, int debug);

int stack_isempty(stack *s);

int stack_isfull(stack *s);

void stack_delete(stack *s);

void stack_push(stack *s,  float xval, float yval, float zval, float deltax, float deltay, float deltaz, int level, int r_id);

void stack_pop(stack *s, float *xval, float *yval, float *zval, float *deltax, float *deltay, float *deltaz, int *level, int * r_id);


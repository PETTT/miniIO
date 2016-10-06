/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.
 * See LICENSE file for details.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "open-simplex-noise.h"
#include "cubes.h"

void cubesinit(cubeInfo *nfo, int task, int levels, int debug) {
  int maxpoints=0;

  maxpoints = pow(8, levels+1);
  nfo->ncubes = 0;
  nfo->npoints = 0;
  nfo->debug = debug;

  if (nfo->debug) {
    printf("Init cubes maxpoints=%d task=%d +++++++++++++++\n", maxpoints, task);
  }
    
  /* Allocate list for all possible cubes and points */
  nfo->points = (float *) malloc((size_t)task*maxpoints*3*sizeof(float) );
  nfo->data = (float *) malloc((size_t)task*maxpoints*sizeof(float) );

}

void cubesfree(cubeInfo *nfo) {
  nfo->debug =0;
  nfo->ncubes = 0;
  nfo->npoints = 0;
  free(nfo->points);
  free(nfo->data);
  nfo->points = NULL;
  nfo->data = NULL;
}

void refine(cubeInfo *nfo, int t, int rpId, float thres, int level_start, float x_start, float y_start, float z_start, float dx_start, float dy_start, float dz_start, struct osn_context *osn, int maxLevel) {
  double noisespacefreq = 10;    /* Spatial frequency of noise */
  double noisetimefreq = 0.25;    /* Temporal frequency of noise */
  float x_center, y_center, z_center, center_val;
  float xpts[8], ypts[8], zpts[8];
  float curdata[8];
  int i;
  int split;
  int refinePathID = 0;
  int new_refinePathID = 0;
  int inside;
  stack octStack;
  float x, y, z, dx, dy, dz;
  int level;
  int stacksize;


  if (nfo->debug)
      printf("Start start initial refinement on Block = %d\n", rpId);
  
  stacksize = 8 * (maxLevel+1);  
  stack_new(&octStack, stacksize, nfo->debug);
  
  stack_push(&octStack, x_start, y_start, z_start, dx_start, dy_start, dz_start, level_start,  rpId);
 
  while (!stack_isempty(&octStack)){
    split = 0;
    inside = 0;
    
    /* get next cube to test */
    stack_pop(&octStack, &x, &y, &z, &dx, &dy, &dz, &level,  &refinePathID);

    /* refinePathID = rpId; */

    if (nfo->debug)
      printf("Block refinePathID = %d\n", refinePathID);
  
    
    /* given a unit cube (8 points) */
    /* - calulate poistion x,y,z for all points */
    xpts[0] = x;
    ypts[0] = y;
    zpts[0] = z;
    xpts[1] = x + dx;
    ypts[1] = y;
    zpts[1] = z;
    xpts[2] = x;
    ypts[2] = y + dy;
    zpts[2] = z;
    xpts[3] = x + dx;
    ypts[3] = y + dy;
    zpts[3] = z;
    xpts[4] = x;
    ypts[4] = y;
    zpts[4] = z + dz;
    xpts[5] = x + dx;
    ypts[5] = y;
    zpts[5] = z + dz;
    xpts[6] = x;
    ypts[6] = y + dy;
    zpts[6] = z + dz;
    xpts[7] = x + dx;
    ypts[7] = y + dy;
    zpts[7] = z + dz;

    /* Get block center value */
    x_center = x + dx/2;
    y_center = y + dy/2;
    z_center = z + dz/2;
    center_val = (float)open_simplex_noise4(osn, x_center*noisespacefreq,
					    y_center*noisespacefreq, z_center*noisespacefreq, t*noisetimefreq);

    if (center_val > thres)
      inside = 1;
      
    if (inside && (level < maxLevel))
      split = 1;

    if (nfo->debug)
      printf("Block_Center: (%f %f, %f)=%f thres=%f  level=%d  inside=%d  split=%d\n", x_center, y_center, z_center, center_val, thres, level, inside, split);
    
    /*   Refine if criteria satisfied */
    if (split) {
      xpts[0] = x;
      ypts[0] = y;
      zpts[0] = z;
      xpts[1] = x + dx/2;
      ypts[1] = y;
      zpts[1] = z;
      xpts[2] = x;
      ypts[2] = y + dy/2;
      zpts[2] = z;
      xpts[3] = x + dx/2;
      ypts[3] = y + dy/2;
      zpts[3] = z;
      xpts[4] = x;
      ypts[4] = y;
      zpts[4] = z + dz/2;
      xpts[5] = x + dx/2;
      ypts[5] = y;
      zpts[5] = z + dz/2;
      xpts[6] = x;
      ypts[6] = y + dy/2;
      zpts[6] = z + dz/2;
      xpts[7] = x + dx/2;
      ypts[7] = y + dy/2;
      zpts[7] = z + dz/2;

      level++;
      for (i=0; i<8; i++) {
	new_refinePathID = (refinePathID *10)+ (i+1);
	  
	if (nfo->debug)
	  printf("Refine on cell %d: (%f %f, %f) ... refinePathID=%d\n", i+1, xpts[i], ypts[i], zpts[i], new_refinePathID);

	stack_push(&octStack, xpts[7-i], ypts[7-i], zpts[7-i], dx/2.0, dy/2.0, dz/2.0, level, new_refinePathID);
      }
    }
    else {
      nfo->ncubes++;
      
      /* - calulate value using open_simplex_noise4 */
      for (i=0; i<8; i++) {
	uint64_t npoints3 = nfo->npoints * 3;
	
	nfo->data[nfo->npoints] =  curdata[i] = (float)open_simplex_noise4(osn, xpts[i]*noisespacefreq,
									   ypts[i]*noisespacefreq, zpts[i]*noisespacefreq, t*noisetimefreq);
	
	nfo->points[npoints3]   = xpts[i];
	nfo->points[npoints3+1] = ypts[i];
	nfo->points[npoints3+2] = zpts[i];
	
	if (nfo->debug)
	  printf("Cube_Point_%d: (%f %f, %f) = %f \n", i,nfo->points[npoints3],nfo->points[npoints3+1] , nfo->points[npoints3+2] , nfo->data[nfo->npoints] );

	nfo->npoints++;
      }

      if (nfo->debug)
	printf("CubeID %llu, refinePathID %d connot be Refined ... refinement level=%d  npoints=%llu +++++++++++++++\n", nfo->ncubes, refinePathID, level, nfo->npoints);
    
    }
  }
  
  if (nfo->debug)
    printf("Ended refinement on initial Block = %d\n", rpId);

  stack_delete(&octStack);


}

void cubeprint(cubeInfo *nfo) {
  int i;

  printf("NumBlock=%llu NumPoints=%llu +++++++++++++++\n", nfo->ncubes, nfo->npoints);
  /* for (i=0; i< nfo->npoints; i++) { */
  /*   printf("%f \n", nfo->data[i]); */
  /* } */
  
}

/* create a new stack */
void stack_new(stack *s, int maxsize, int debug) {
  cubeItem *newcubes;

  newcubes = (cubeItem *) malloc((size_t)maxsize*sizeof(cubeItem) );
 
  s->cubes = newcubes;
  s->top = -1;
  s->maxsize = maxsize;
  s->size = 0;
  s->debug = debug;

  if (s->debug) {
    printf("Created stack\n");
  }  
}

/* check if stack is empty */
int stack_isempty(stack *s) {
  int empty=0;

  if (s->top < 0)
    empty=1;

  return empty;
}

/* check if stack is full */
int stack_isfull(stack *s) {
  int full=0;

  if (s->top >= s->maxsize)
    full=1;

  return full;
}

/* delete a stack */
void stack_delete(stack *s) {

  if (s->debug) {
    printf("Start to delete stack\n");
  }  
  free(s->cubes);
  s->cubes = NULL;
  s->top = -1;
  s->maxsize = 0;
  s->size = 0;
  s->debug = 0;
  
  if (s->debug) {
    printf("Deleted stack\n");
  }  
}

/* pushes ne element onto stack */
void stack_push(stack *s,  float xval, float yval, float zval, float deltax, float deltay, float deltaz, int level, int r_id) {

  cubeItem newCube;
 
  if (stack_isfull(s)) {
    printf("ERROR: Stack FULL .... Could not push new cube data onto stack ...\n");
    exit(1);
  }

  newCube.r_id = r_id;
  newCube.level = level;
  newCube.x = xval;
  newCube.y = yval;
  newCube.z = zval;
  newCube.dx = deltax;
  newCube.dy = deltay;
  newCube.dz = deltaz;

  s->top++;
  s->size++;
  s->cubes[s->top] = newCube;

  if (s->debug) {
    printf("Pushed Cube LRcorner position  (%f %f, %f) stack_size=%d\n", s->cubes[s->top].x, s->cubes[s->top].y, s->cubes[s->top].z, s->size);
    printf("Pushed Cube delta values (%f %f, %f)\n", s->cubes[s->top].dx, s->cubes[s->top].dy, s->cubes[s->top].dz);
  }
}

/* removes top element from stack,  */
void stack_pop(stack *s, float *xval, float *yval, float *zval, float *deltax, float *deltay, float *deltaz, int *level, int *r_id) {
  
  if (stack_isempty(s))  {
    printf("ERROR: Stack EMPTY ... Could not pop cube data off stack ...\n");
    exit(1);
  }

  
  /* update data */
  *r_id =  s->cubes[s->top].r_id;
  *level = s->cubes[s->top].level;
  *xval = s->cubes[s->top].x;
  *yval = s->cubes[s->top].y;
  *zval = s->cubes[s->top].z;
  *deltax = s->cubes[s->top].dx;
  *deltay = s->cubes[s->top].dy;
  *deltaz = s->cubes[s->top].dz;
  
  if (s->debug) {
    printf("Popped Cube LRcorner position (%f %f, %f) stack_size=%d\n",  s->cubes[s->top].x, s->cubes[s->top].y, s->cubes[s->top].z, s->size);
    printf("Popped Cube delta values (%f %f, %f)\n", s->cubes[s->top].dx, s->cubes[s->top].dy, s->cubes[s->top].dz);
  }

   s->top--;
   s->size--;
}

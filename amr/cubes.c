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

  if (nfo->debug)
    printf("Init cubes maxpoints=%d +++++++++++++++\n", maxpoints);
    
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

void refine(cubeInfo *nfo, int t, int rpId, float thres, float value, int level_start, float x_start, float y_start, float z_start, float dx_start, float dy_start, float dz_start, struct osn_context *osn, int maxLevel) {
  double noisespacefreq = 10;    /* Spatial frequency of noise */
  double noisetimefreq = 0.25;    /* Temporal frequency of noise */
  float x_center, y_center, z_center, center_val;
  float xpts[8], ypts[8], zpts[8];
  float curdata[8];
  float conns[12];
  uint64_t nocts=0;        /* Number of cubes */
  int i;
  int below = 0;
  int above = 0;
  int split = 0;
  int exact = 0;
  int refinePathID = 0;
  int new_refinePathID = 0;
  int inside = 0;
  stack * octStack;
  float x, y, z, dx, dy, dz;
  int level;


  octStack =  new_stack(nfo->debug);

  stack_push(octStack, x_start, y_start, z_start, dx_start, dy_start, dz_start, level_start);
  nfo->npoints = 0;
  nfo->ncubes = 0;

  while (!stack_isempty(octStack)){
    split = 0;

    /* get next cube to test */
    stack_pop(octStack, &x, &y, &z, &dx, &dy, &dz, &level);

    refinePathID = rpId;

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

    if (nfo->debug)
      printf("Block_Center: (%f %f, %f)=%f  level=%d\n", x_center, y_center, z_center, center_val, level);
    
    if (center_val > thres)
      inside = 1;
      
    if (inside && (level < maxLevel))
      split = 1;
  
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

	stack_push(octStack, xpts[7-i], ypts[7-i], zpts[7-i], dx/2.0, dy/2.0, dz/2.0, level);
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
	  printf("Point_%d: (%f %f, %f) = %f \n", i,nfo->points[npoints3],nfo->points[npoints3+1] , nfo->points[npoints3+2] , nfo->data[nfo->npoints] );

	nfo->npoints++;
      }

      if (nfo->debug)
	printf("Block %d connot be Refined ... refinement level=%d +++++++++++++++\n", refinePathID, level);
    
    }
  }
  stack_delete(octStack);
}

void cubeprint(cubeInfo *nfo) {
  int i;

  printf("NumBlock=%llu NumPoints=%llu +++++++++++++++\n", nfo->ncubes, nfo->npoints);
  /* for (i=0; i< nfo->npoints; i++) { */

  /*   printf("%f \n", nfo->data[i]); */
  /* } */
  
}

/* create a new stack */
stack *new_stack(int debug) {

  stack *s = malloc(sizeof(stack));
  s->top = NULL;
  s->size = 0;
  s->debug = debug;

  return s;
}

/* check if stack is empty */
int stack_isempty(stack *s) {

  return s->top == NULL ? 1 : 0;
  
}

void stack_clean(stack *s)
{
  float *xval, *yval, *zval, *deltax, *deltay, *deltaz;
  int *level;
  
  while(!stack_isempty(s))
    stack_pop(s, xval, yval, zval, deltax, deltay, deltaz, level);
}

/* completly delete a stack */
void stack_delete(stack *s) {

   stack_clean(s);
   free(s);
}

/* pushes a value onto the stack */
/* return 1 if element was pushed on successfully, 0 otherwise */
void stack_push(stack *s,  float xval, float yval, float zval, float deltax, float deltay, float deltaz, int level) {

 cubeItem *newNode = malloc(sizeof(cubeItem));
 
  if (newNode == NULL) {
    printf("ERROR: Could not push new cube data onto stack\n");
    exit(0);
  }
  
  newNode->level = level;
  newNode->x = xval;
  newNode->y = yval;
  newNode->z = zval;
  newNode->dx = deltax;
  newNode->dy = deltay;
  newNode->dz = deltaz;
      
  newNode->next = s->top;
  s->top = newNode;
  s->size += 1;

  if (s->debug) {
    printf("Push Cube (%f %f, %f) stackszie=%d\n", s->top->x, s->top->y, s->top->z, s->size);
    printf("Deltas (%f %f, %f)\n", s->top->dx, s->top->dy, s->top->dz);
  }
}

/* removes top element from stack,  */

void stack_pop(stack *s, float *xval, float *yval, float *zval, float *deltax, float *deltay, float *deltaz, int *level) {

  cubeItem *oldTop;
  
  if (s == NULL || s->top == NULL )  {
    printf("Stack size=%d\n", s->size);
    printf("ERROR: Could not pop cube data off stack\n");
    exit(0);
  }

  /* update data */
  *level = s->top->level;
  *xval = s->top->x;
  *yval = s->top->y;
  *zval = s->top->z;
  *deltax = s->top->dx;
  *deltay = s->top->dy;
  *deltaz = s->top->dz;

  if (s->debug) {
    printf("Pop Cube (%f %f, %f)\n", s->top->x, s->top->y, s->top->z);
    printf("Deltas (%f %f, %f)\n", s->top->dx, s->top->dy, s->top->dz);
  }
  
  oldTop = s->top;
  s->top = oldTop->next;
  s->size -= 1;
  free(oldTop);
  oldTop = NULL;

  if (s->debug) {
    printf("Stack size=%d\n", s->size);
  }
}

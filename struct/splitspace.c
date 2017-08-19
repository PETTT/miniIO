/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

#include <stdlib.h>  
#include <stdio.h>
#include <math.h>
#include "splitspace.h"


void init_primeFactors(struct factorstruct *factors, int n)
{
  int max_size;

  max_size = (n+1)/2;
  
  factors->n = n;
  factors->size = 0;
  factors->fac = (int *) malloc(max_size * sizeof(int));
}


void free_primeFactors(struct factorstruct *factors)
{
 
  factors->n = 0;
  factors->size = 0;
  free(factors->fac);
  factors->fac = NULL;

}


void  get_primeFactors(struct factorstruct *factors)
{
  int d, n;
  
  n = factors->n;
  d = 2;
  while (n > 1) {
      while (n % d == 0) {
	factors->fac[factors->size] = d;
	factors->size++;
	n = n/d;
      }
      d = d+1;
      if (d*d > n) {
	if (n > 1) {
	  factors->fac[factors->size] = n;
	  factors->size++;
	}
	break;
      }
  }

}

void  print_primeFactors(struct factorstruct *factors)
{
  int i;

  if (factors->n !=1) {
    printf("List of Prime Factors of %d = (", factors->n);
    for (i=0; i<factors->size; i++) {
      if (i == factors->size-1) {
	printf("%d)\n", factors->fac[i]);
      }
      else {
	printf("%d,", factors->fac[i]);
      }
    }
  }
}


void init_rects(recList *l, struct factorstruct *factors)
{
  l->top = 0;
  l->size = 0;
  l->max_size = factors->n;
  
  l->recs = (recItem *) malloc(l->max_size * sizeof(recItem));
}

void empty_rects(recList *l)
{
  l->top = 0;
  l->size = 0;
  l->max_size = 0;
}

void free_rects(recList *l)
{
  empty_rects(l);
  free(l->recs);
  l->recs = NULL;
}

void print_rects(recList *l)
{
  int i;

  printf("\nList of rects (%d):\n", l->size);
  for (i=0; i<l->size; i++) {
    printf("rect %d = (%d,%d;%d,%d)\n", i, l->recs[i].x0,  l->recs[i].y0,  l->recs[i].x1,  l->recs[i].y1);
  
  }
}

void add_rect(recList *l, int x0, int y0, int x1, int y1) {
  int debug=0;

  l->recs[l->size].x0 = x0;
  l->recs[l->size].y0 = y0;
  l->recs[l->size].x1 = x1;
  l->recs[l->size].y1 = y1;
  l->size++;

  if (debug) {
    printf("Added rect %d = (%d,%d;%d,%d)\n", l->size, l->recs[l->size-1].x0,  l->recs[l->size-1].y0,  l->recs[l->size-1].x1,  l->recs[l->size-1].y1);
  }
}

void copy_rectList(recList *l, recList *l2) {
  int i;

  /* printf ("Copy rectList\n"); */
  empty_rects(l);
  for (i=0; i<l2->size; i++) {
    add_rect(l, l2->recs[i].x0,  l2->recs[i].y0,  l2->recs[i].x1,  l2->recs[i].y1);
  }
}

void append_rectList(recList *l, recList *l2) {
  int i;
  
  for (i=0; i<l2->size; i++) {
    add_rect(l, l2->recs[i].x0,  l2->recs[i].y0,  l2->recs[i].x1,  l2->recs[i].y1);
  }
}

void init_splitspace(struct factorstruct *factors, recList *rects_list, int nprocs) {

  init_primeFactors(factors, nprocs);
  get_primeFactors(factors);
  init_rects(rects_list, factors);
}

void splitspace(recList *l, int x_dims, int y_dims, struct factorstruct *factors, int *data)
{
  int i, j;
  recList splitList, tmpList;
  recItem rec;
  
  init_rects(&splitList, factors);
  init_rects(&tmpList, factors);
  add_rect(l, 0, 0, x_dims-1, y_dims-1);

  /* loop through every factor */
  for (i=0; i<factors->size; i++) {
    empty_rects(&splitList);
    /* loop through every rect */
    for (j=0; j< l->size; j++) {
      rec =  l->recs[j];
      subdivide_rects(&tmpList, rec, data, factors->fac[i]);
      append_rectList(&splitList, &tmpList);
      empty_rects(&tmpList);
    }
    copy_rectList(l, &splitList);
    /* print_rects(l); */
  }
 
  free_rects(&splitList);
  free_rects(&tmpList);
}


void free_splitspace(struct factorstruct *factors, recList *rects_list) {
  free_primeFactors(factors);
  free_rects(rects_list);
}

void subdivideArray(int *cvals, subArrayItem *subarray, int size, int n) {
  int i, j, ii, subarray_size, move_size, cur_id;
  int *movedlt, *movedrt;
  int chunksize, remainder;
  int done;
  int sumA, sumA_size, sumA_last_index, subA_last_val;
  int sumB, sumB_size, sumB_first_index, subB_first_val;
  int debug = 0;
  
  move_size = n-1;
  movedlt = (int *) malloc( move_size * sizeof(int));
  movedrt = (int *) malloc( move_size * sizeof(int));

  chunksize = size /n;
  remainder = size%n;
  
  if (debug) {
    printf("size=%d n_split=%d chunksize=%d  remainder=%d\n", size, n, chunksize, remainder);
  }
  
  
  subarray[0].start_index =j;
  subarray[0].size = 0;
  cur_id = 0;
  for(j = 0; j < n; j++) {
    subarray_size = chunksize;
    if (remainder > 0) {
      subarray_size++;
      remainder--;
    }
    subarray[j].start_index = cur_id;
    subarray[j].size = 0;
    for(i = 0; i < subarray_size; i++) {
      cur_id++;
      subarray[j].size++;
    }
  }

  if (debug) {
    for(i = 0; i < n ; i++) {
      printf("Initial subarray[%d]: start_index=%d  size=%d\n", i, subarray[i].start_index, subarray[i].size);
    }
    printf("\n ");
  }
   
  for (i=0; i<move_size; i++) {
    movedlt[i] = 0;
    movedrt[i] = 0;
  }
  
  done = 0;
  while ( done == 0) {
    done = 1;
    for (i=0; i<move_size; i++) {

      sumA = subarray_sum(cvals, subarray[i]);
      sumA_size = subarray[i].size;
      sumA_last_index = subarray[i].start_index + subarray[i].size-1;
      subA_last_val = cvals[ sumA_last_index];
      sumB = subarray_sum(cvals, subarray[i+1]);
      sumB_size = subarray[i+1].size;
      sumB_first_index = subarray[i+1].start_index;
      subB_first_val = cvals[sumB_first_index];
      
      if( ((sumB - sumA) >= subB_first_val) && (sumB_size > 1) && (movedrt[i]==0) ) {

	subarray[i].size++;
	subarray[i+1].size--;
	subarray[i+1].start_index++;
	movedlt[i] = 1;
	done = 0;

	if (debug) {
	  printf("** Shift LEFT:\n ");
	  for(ii = 0; ii < n ; ii++) {
	    printf("** subarray[%d]: start_index=%d  size=%d\n", ii, subarray[ii].start_index, subarray[ii].size);
	  }
	}
      }
      else if ( ((sumA - sumB) >= subA_last_val) && (sumA_size > 1) && (movedlt[i]==0) ) {

	if (debug) {
	  printf("** Initial:\n"); 
	  for(ii = 0; ii < n ; ii++) {
	    printf("subarray[%d]: start_index=%d  size=%d\n", ii, subarray[ii].start_index, subarray[ii].size);
	  }
	}
	
	subarray[i+1].size++;
	subarray[i].size--;
	subarray[i+1].start_index--;
	movedrt[i] = 1;
	done = 0;
	
	if (debug) {
	  printf("** Shift RIGHT:\n ");
	  for(ii = 0; ii < n ; ii++) {
	    printf("** subarray[%d]: start_index=%d  size=%d\n", ii, subarray[ii].start_index, subarray[ii].size);
	  }
	}
      }
    }
    
  }

  if (debug) {
    for(i = 0; i < n ; i++) {
      printf("Ending subarray[%d]: start_index=%d  size=%d\n", i, subarray[i].start_index, subarray[i].size);
    }
    printf("\n ");
  }
  
  free(movedlt);
  free(movedrt);
  
}

int subarray_sum(int *cvals, subArrayItem subarray) {
  int i, sum;
  int start_index;
  int end_index;

  sum = 0;
  start_index = subarray.start_index;
  end_index = start_index + subarray.size;
  
  for (i=start_index; i < end_index; i++) {
    sum += cvals[i];
  }

  return sum;
}

float subarrays_var(int *cvals, subArrayItem *subarray, int n) {
  int i;
  float total;
  float mean;
  float var;

  var = 0.0;
  mean = 0.0;
  total = 0.0;
  for (i=0; i<n; i++) {
    total += subarray_sum(cvals, subarray[i]);
  }
  mean = total/ n;
  
  total = 0;
  for (i=0; i<n; i++) {
    total += pow( (subarray_sum(cvals, subarray[i]) - mean), 2);
  }
  var = total/n;

  return var;
 
}

/* Split a Rect into n sub-Rects that have roughly the same data value sums. */
void subdivide_rects(recList *l, recItem rec, int *data, int n) {
  int i, j, x0, y0, point_id, cur_id;
  int x_dims, y_dims, xy_dims, x_start, x_end, y_start, y_end;
  int *cx, *cy;
  float cxvar, cyvar;
  subArrayItem *cx_subarray, *cy_subarray;
  int splitx=0;
  int debug=0;

  x_start = rec.x0;
  y_start = rec.y0;
  x_end = rec.x1+1;
  y_end = rec.y1+1;
  x_dims = x_end - x_start;
  y_dims = y_end - y_start;
  xy_dims = x_dims * y_dims;
  
  cx = (int *) malloc(x_dims * sizeof(int));
  cx_subarray = (subArrayItem *) malloc( n * sizeof(subArrayItem));
  cy_subarray = (subArrayItem *) malloc( n * sizeof(subArrayItem));
  cy = (int *) malloc(y_dims * sizeof(int));

  if (debug) {
    printf("\nx_dims=%d y_dims=%d xy_dims=%d\n", x_dims, y_dims, xy_dims);
    for(i = x_start; i < x_end; i++) {
      for(j = y_start; j < y_end; j++) {
	/* calculate index */
	point_id = (j * x_dims) + i;
	if( j == 0) {
	  printf("[ %d ", data[point_id]);
	}
	else  {
	  printf("%d ",data[point_id]);
	}
      }
      printf("]\n");
    }
  }

  /* Init array */ 
  for(i = 0; i < x_dims; i++) {
    cx[i] = 0;
  }

  /* sum  across x */
  cur_id =0;
  for(i = x_start; i < x_end; i++) {
    for(j = y_start; j < y_end; j++) {
      /* calculate point index */
      point_id = (j * x_dims) + i;
      cx[cur_id] += data[point_id];
    }
    cur_id++;
  }

  if (debug) {
    printf("cx =[ ");
    for(i = 0; i < x_dims; i++) {
      printf("%d ", cx[i]);
    }
    printf(" ]\n");
  }

  subdivideArray(cx, cx_subarray, x_dims, n);
  
  cxvar = subarrays_var(cx, cx_subarray, n);
  if (debug) { printf("cxvar %f\n", cxvar);}


  /* Init array */ 
  for(i = 0; i < y_dims; i++) {
    cy[i] = 0;
  }

   /* sum  across y */
  cur_id =0;
  for(j = y_start; j < y_end; j++) {
    for(i = x_start; i < x_end; i++) {
      /* calculate point index */
      point_id = (j * x_dims) + i;
      cy[cur_id] += data[point_id];
    }
    cur_id++;
  }

  if (debug) {
    printf("cy =[ ");
    for(i = 0; i < y_dims; i++) {
      printf("%d ", cy[i]);
    }
    printf(" ]\n");
  }

  subdivideArray(cy, cy_subarray, y_dims, n);
  cyvar = subarrays_var(cy, cy_subarray, n);
  if (debug) { printf("cyvar %f\n", cyvar); }

  /* Pick the axis that gives the best load balance */
  if (cxvar < cyvar) {
    splitx = 1;
  }
  else if (cxvar > cyvar) {
    splitx = 0;
  }
  else if ( (rec.x1 - rec.x0) >= (rec.y1 - rec.y0) ) {   /* if equal balance, use largest axis */
    splitx = 1;
  }
  else{
      splitx = 0;
  }
  
  if (debug) {
    printf("splitx = %d\n", splitx);
    printf("rect = (%d,%d;%d,%d)\n", rec.x0,  rec.y0,  rec.x1, rec.y1);
  }

  
  /* Convert array splits to list of rect's */
  if (splitx) {
    add_rect(l, rec.x0, rec.y0, rec.x0+cx_subarray[0].size-1, rec.y1);
  }
  else {
    add_rect(l, rec.x0, rec.y0, rec.x1, rec.y0+cy_subarray[0].size-1);
  }

  for (i=1; i<n; i++) {
    if (splitx) {
      x0 = l->recs[i-1].x1+1;
      add_rect(l, x0, rec.y0, x0+cx_subarray[i].size-1, rec.y1);
    }
    else {
      y0 = l->recs[i-1].y1 + 1;
      add_rect(l,rec.x0, y0, rec.x1, y0+cy_subarray[i].size-1);
    }
  }
   
  free (cx);
  free (cx_subarray);
  free (cy);
  free (cy_subarray); 
}

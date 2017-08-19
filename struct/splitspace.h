/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

struct factorstruct {
  int *fac;
  int size;
  int n;
};

typedef struct subArrayItem {
 
  int start_index;
  int size;

} subArrayItem;

/* holds an entry in the stack (in this case cube information) */
typedef struct recItem {
 
  int x0;
  int y0;
  int x1;
  int y1;
} recItem;

/* A stack that holds and array of stack items (in this case an array of cube information) */
typedef struct {
  recItem *recs;
  int top;
  int size;
  int max_size;
} recList;

void init_primeFactors(struct factorstruct *factors, int n);
void free_primeFactors(struct factorstruct *factors);
void get_primeFactors(struct factorstruct *factors);
void print_primeFactors(struct factorstruct *factors);
void init_rects(recList *l, struct factorstruct *factors);
void add_rect(recList *l, int x0, int y0, int x1, int y1);
void free_rects(recList *l);
void empty_rects(recList *l);
void print_rects(recList *l);
void count_rects(recList *l);
void copy_rectList(recList *l, recList *l2);
void append_rectList(recList *l, recList *l2);
void subdivideArray(int *cvals, subArrayItem *subarray, int size, int n);
int subarray_sum(int *cvals, subArrayItem subarray);
float subarrays_var(int *cvals, subArrayItem *subarray, int n);
void subdivide_rects(recList *l, recItem rec, int *data, int n) ;
void init_splitspace(struct factorstruct *factors, recList *rects_list, int nprocs);
void splitspace(recList *l, int x_dims, int y_dims, struct factorstruct *factors, int *data);
void free_splitspace(struct factorstruct *factors, recList *rects_list);

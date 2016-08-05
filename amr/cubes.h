typedef struct cubeInfo {
  int debug;
  uint64_t ncubes;       /* Number of cubes */
  uint64_t npoints;      /* Number of cube points */
  float *points;         /* Points of cube */
  float *data;           /* Values on points */
}cubeInfo;


// holds an entry in the stack and point to the next element if there is any
typedef struct cubeItem {
  int level;
  float x;
  float y;
  float z;
  float dx;
  float dy;
  float dz;
  struct cubeItem * next;
} cubeItem;

// the actual stack, holds the reference to the top and other meta information
// all operation will only require a reference to this
typedef struct {
    cubeItem *top;
    int size;
    int debug;
} stack;


void cubesinit(cubeInfo *nfo, int task, int levels, int debug);

void cubesfree(cubeInfo *nfo);

void refine(cubeInfo *nfo, int t, int rpId, float thres, float value, int level_start, float x_start, float y_start,
	    float z_start, float dx_start, float dy_start, float dz_start, struct osn_context *osn, int maxLevel);

void cubeprint(cubeInfo *nfo);

stack * new_stack(int debug);

int stack_isempty(stack *s);

void stack_clean(stack *s);

void stack_delete(stack *s);

void stack_push(stack *s,  float xval, float yval, float zval, float deltax, float deltay, float deltaz, int level);

void stack_pop(stack *s, float *xval, float *yval, float *zval, float *deltax, float *deltay, float *deltaz, int *level);


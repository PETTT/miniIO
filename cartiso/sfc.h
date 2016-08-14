/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

struct sfc3_ctx {
    int ni, nj, nk;   /* Total size of space */
    int reverse;      /* Whether the curve advances in reverse from the end */
    int ndx;          /* Current index */
    int valid;        /* Are the current indices valid */
    int i, j, k;      /* Current coordinates in space */
    int idir;         /* Current direction along i, +1 or -1 */
    int jdir;         /* Current direction along j, +1 or -1 */
    int kdir;         /* Current direction along k, +1 or -1 */
};

void sfc3_serpentine_init(struct sfc3_ctx *ctx, int ni, int nj, int nk, int reverse);
void sfc3_serpentine_next(struct sfc3_ctx *ctx);


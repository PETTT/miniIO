/*
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include "sfc.h"

void sfc3_serpentine_init(struct sfc3_ctx *ctx, int ni, int nj, int nk, int reverse)
{
    ctx->ni = ni;
    ctx->nj = nj;
    ctx->nk = nk;
    ctx->reverse = reverse;
    ctx->ndx = 0;
    ctx->valid = 1;
    ctx->i = 0;
    ctx->j = 0;
    ctx->k = 0;
    ctx->idir = 1;
    ctx->jdir = 1;
    ctx->kdir = 1;
    if(reverse) {
        ctx->k = nk - 1;
        ctx->kdir = -1;
        if(nk % 2) {      /* nk odd => j ends opposite */
            ctx->j = nj - 1;
            ctx->jdir = -1;
            if(nj % 2) {       /* nj & nk odd => i ends opposite */
                ctx->i = ni - 1;
                ctx->idir = -1;
            }
        }
    }
}

void sfc3_serpentine_next(struct sfc3_ctx *ctx)
{
    ctx->i += ctx->idir;
    if(ctx->i == ctx->ni) {
        ctx->i -= 1;
        ctx->j += ctx->jdir;
        ctx->idir = -1;
    } else if(ctx->i == -1) {
        ctx->i += 1;
        ctx->j += ctx->jdir;
        ctx->idir = 1;
    }
    if(ctx->j == ctx->nj) {
        ctx->j -= 1;
        ctx->k += ctx->kdir;
        ctx->jdir = -1;
    } else if(ctx->j == -1) {
        ctx->j += 1;
        ctx->k += ctx->kdir;
        ctx->jdir = 1;
    }
    ctx->ndx += 1;
    if(ctx->k >= ctx->nk || ctx->k < 0)    /* In either case, we're done */
        ctx->valid = 0;
}

#ifdef SFC_SERPENTINE_TEST

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
    struct sfc3_ctx sfc;
    int ni, nj, nk, reverse;

    if(argc != 5) {
        printf("Invalid arguments.\nUsage: %s ni nj nk reverse\n", argv[0]);
        printf("   reverse = 0 or 1     1 runs the curve in reverse from the end");
        return 1;
    }

    ni = atoi(argv[1]);
    nj = atoi(argv[2]);
    nk = atoi(argv[3]);
    reverse = atoi(argv[4]);

    if(ni < 1 || nj < 1 || nk < 1) {
        printf("Invalid number or ni, nj, or nk.\n");
        return 1;
    }

    sfc3_serpentine_init(&sfc, ni, nj, nk, reverse);

    do {
        printf("%d %d %d\n", sfc.i, sfc.j, sfc.k);
        sfc3_serpentine_next(&sfc);
    } while(sfc.valid);

    return 0;
}

#endif

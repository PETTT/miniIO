/* 
 * Functions to assist with argument handling
 *
 * Copyright (c) DoD HPCMP PETTT.  All rights reserved.  
 * See LICENSE file for details.
 */

#include <stdlib.h>
#include <string.h>

/* Convert a string with/out trailing 'x' to an integer and bool (with/out the x)
 *    s: pointer to input string
 *    i: pointer to an output integer
 *    hasx: pointer to output boolean: true of string ends in x
 */

void atoix(const char *s, int *i, int *hasx)
{
    char *ss = strdup(s);
    int len = strlen(ss);
    if(ss[len-1] == 'x' || ss[len-1] == 'X') {
        *hasx = 1; 
        ss[len-1] = 0;
        len--;
    } else
        *hasx = 0;
    *i = atoi(ss);
    free(ss);
}

/* Convert a string with/out trailing 'x' to a double and bool (with/out the x)
 *    s: pointer to input string
 *    d: pointer to an output integer
 *    hasx: pointer to output boolean: true of string ends in x
 */

void atodx(const char *s, double *d, int *hasx)
{
    char *ss = strdup(s);
    int len = strlen(ss);
    if(ss[len-1] == 'x' || ss[len-1] == 'X') {
        *hasx = 1; 
        ss[len-1] = 0;
        len--;
    } else
        *hasx = 0;
    *d = strtod(ss, NULL);
    free(ss);
}


#include <stdio.h>

/*  ---- globals required by bozorth3.c --------------------------- */
FILE *errorfp = NULL;          /* NBIS prints warnings here; NULL disables them */

int  m1_xyt  = 0;              /* 0 = CW angle math (matches NBIS default)      */
int  min_computable_minutiae = 15;  /* same hard-coded default NBIS uses       */
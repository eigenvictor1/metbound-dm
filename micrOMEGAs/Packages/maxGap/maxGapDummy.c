#include <stdio.h>
#include <math.h>

extern float upperlim_(float *CL,int *If,int *N,float*FC,float *muB,float*FB,int*Iflag);

float upperlim_(float *CL,int *If,int *N,float*FC,float *muB,float*FB,int*Iflag)
{  static int nWorn=0;
   if(nWorn==0)  printf("The CRESST_III exclusion cannot be determined because the Optimal Interval  Method needs Fortran.\n");
   nWorn=1;
   return NAN;
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>

//typedef struct
//{
//   double dat[2];
//} gsl_complex;

int main()
{
   double t;
   gsl_complex z;   

   z=gsl_complex_rect(1,2);
   
   return 0;
}

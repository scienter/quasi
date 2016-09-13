#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void loadLaser(Domain *D,LaserList *L,double t)
{
  void loadLaser1D(Domain *D,LaserList *L,double t);

  switch(D->dimension)  {
  case 1 :
    loadLaser1D(D,L,t);
    break;
  default :
    printf("In loadLaser, what is field_type? and what is dimension?\n");
  }
}


void loadLaser1D(Domain *D,LaserList *L,double t)
{
   double rU,rD,flat,x,dx,offset,a0;
   int i,j=0,k=0,istart,iend,ibound;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;
   iend=D->iend;

   rU=L->rU;
   rD=L->rD;
   flat=L->flat;
   dx=D->dx;
   a0=L->amplitude;

   offset=3;
   ibound=D->nx-(int)(offset*rU*D->division);
   for(i=istart; i<ibound; i++)  {
     x=(D->domainMinX+i-istart)*dx;
     D->aNow[i][j][k].real=a0*exp(-(x-D->domainMaxX+offset*rU)*(x-D->domainMaxX+offset*rU)/rD/rD);
     D->aNow[i][j][k].img=0.0;
//     D->aNow[i][j][k].img=-1.0/L->k0*D->aNow[i][j][k].real;
     D->aOld[i][j][k].real=D->aNow[i][j][k].real;
     D->aOld[i][j][k].img=D->aNow[i][j][k].img;
     D->aOld2[i][j][k].real=D->aNow[i][j][k].real;
     D->aOld2[i][j][k].img=D->aNow[i][j][k].img;
     D->aOld3[i][j][k].real=D->aNow[i][j][k].real;
     D->aOld3[i][j][k].img=D->aNow[i][j][k].img;
   }

   for(i=ibound; i<iend; i++)  {
     x=(D->domainMinX+i-istart)*dx;
     D->aNow[i][j][k].real=a0*exp(-(x-D->domainMaxX+offset*rU)*(x-D->domainMaxX+offset*rU)/rU/rU);
     D->aNow[i][j][k].img=0.0;
//     D->aNow[i][j][k].img=-1.0/L->k0*D->aNow[i][j][k].real;
     D->aOld[i][j][k].real=D->aNow[i][j][k].real;
     D->aOld[i][j][k].img=D->aNow[i][j][k].img;
     D->aOld2[i][j][k].real=D->aNow[i][j][k].real;
     D->aOld2[i][j][k].img=D->aNow[i][j][k].img;
     D->aOld3[i][j][k].real=D->aNow[i][j][k].real;
     D->aOld3[i][j][k].img=D->aNow[i][j][k].img;
   }
/*
   for(i=istart; i<iend; i++)  {
     x=(D->domainMinX+i-istart)*dx;
     printf("%g %g\n",x,D->aNow[i][j][k].real);
   }
*/   
}


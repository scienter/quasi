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
     D->aNow[i][j][k].real=a0*exp(-(x-D->domainMaxX+offset*rU)*(x-D->domainMaxX+offset*rU)/rD/rD)*cos(L->k0*(i-iend)*dx);
     D->aNow[i][j][k].img=a0*exp(-(x-D->domainMaxX+offset*rU)*(x-D->domainMaxX+offset*rU)/rD/rD)*sin(L->k0*(i-iend)*dx);
   }
   for(i=ibound; i<iend; i++)  {
     x=(D->domainMinX+i-istart)*dx;
     D->aNow[i][j][k].real=a0*exp(-(x-D->domainMaxX+offset*rU)*(x-D->domainMaxX+offset*rU)/rU/rU)*cos(L->k0*(i-iend)*dx);
     D->aNow[i][j][k].img=a0*exp(-(x-D->domainMaxX+offset*rU)*(x-D->domainMaxX+offset*rU)/rU/rU)*sin(L->k0*(i-iend)*dx);
   }
/*
   for(i=istart; i<iend; i++)  {
     x=(D->domainMinX+i-istart)*dx;
     printf("%g %g\n",x,D->aNow[i][j][k].real);
   }
*/   
}


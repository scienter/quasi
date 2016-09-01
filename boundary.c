#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>


void boundary(Domain *D)
{
  FILE *out;
  int nxSub1D,nySub2D,nzSub3D;
  int myrank, nTasks,a;
  double ***memoryAsign();
  Laser ***memoryAsignLaser();
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

  D->nxSub=D->nx;
  D->istart=2;
  D->iend=D->nxSub+2;

  nxSub1D=D->nxSub+5;
  nySub2D=1;
  nzSub3D=1;

  D->maxXSub=D->minXSub+D->nx;

  // laser setting
  D->aNext=memoryAsignLaser(nxSub1D,nySub2D,nzSub3D);
  D->aNow=memoryAsignLaser(nxSub1D,nySub2D,nzSub3D);
  D->aOld=memoryAsignLaser(nxSub1D,nySub2D,nzSub3D);

}

double ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   double ***field;

   field = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)
   {
     field[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       field[i][j] = (double *)malloc((nz)*sizeof(double ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k]=0.0;
       }

   return field;
}

Laser ***memoryAsignLaser(int nx, int ny, int nz)
{
   int i,j,k;
   Laser ***field;

   field = (Laser ***)malloc((nx)*sizeof(Laser **));
   for(i=0; i<nx; i++)
   {
     field[i] = (Laser **)malloc((ny)*sizeof(Laser *));
     for(j=0; j<ny; j++)
       field[i][j] = (Laser *)malloc((nz)*sizeof(Laser ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++){
         field[i][j][k].real=0.0;
         field[i][j][k].img=0.0;
       }

   return field;
}

#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>


void boundary(Domain *D)
{
  FILE *out;
  int i,j,k,s,nxSub1D,nySub2D,nzSub3D;
  int startj,startk;
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
   D->aOld2=memoryAsignLaser(nxSub1D,1,1);
   D->aOld3=memoryAsignLaser(nxSub1D,1,1);

   // current setting
   D->Jx=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->Jy=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->Jz=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);

   // field setting
   D->Ex=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->Ey=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->Ez=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->Bx=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->By=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->Bz=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->BxNow=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->ByNow=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);
   D->BzNow=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);

   // gamma,density setting
   D->gamma=memoryAsign(nxSub1D,nySub2D,nzSub3D,1.0);
   D->Rho=memoryAsign(nxSub1D,nySub2D,nzSub3D,0.0);


   // Particle setting
   nxSub1D=D->nxSub+3;
   nySub2D=1;
   nzSub3D=1;
   startj=0;
   startk=0;

   D->particle = (Particle ***)malloc((nxSub1D)*sizeof(Particle **));
   for(i=0; i<nxSub1D; i++) {
     D->particle[i] = (Particle **)malloc((nySub2D)*sizeof(Particle *));
     for(j=0; j<nySub2D; j++)
       D->particle[i][j] = (Particle *)malloc((nzSub3D)*sizeof(Particle ));
   }

   // setting up particle's pointer
   if(D->dimension==1)
   {
     j=k=0;
     for(i=0; i<D->iend+1; i++)       //i starts at 0 because of boost frame
     {
       D->particle[i][j][k].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
       for(s=0; s<D->nSpecies; s++)       {
         D->particle[i][j][k].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
         D->particle[i][j][k].head[s]->pt = NULL;
       }
     }
   }




}

double ***memoryAsign(int nx, int ny, int nz,double value)
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
         field[i][j][k]=value;
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

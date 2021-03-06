#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"

void fieldSolve(Domain *D)
{
   int myrank, nTasks,rank,rankM,rankN;
   void Bsolve1D(Domain *D);
   void Esolve1D(Domain *D);

   MPI_Status status;
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   switch(D->dimension) {
   //1D field
   case 1:
     Bsolve1D(D);
     Esolve1D(D);
     break;

  }

}

void Bsolve1D(Domain *D)
{
    int i,j,k,istart,iend,nxSub;
    double dx,dy,dz,dt,oldBx,oldBy,oldBz;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;

    istart=D->istart;
    iend=D->iend;
    j=k=0;

    for(i=istart; i<iend; i++)
    {
      oldBx=D->Bx[i][j][k];
      oldBy=D->By[i][j][k];
      oldBz=D->Bz[i][j][k];
      D->Bx[i][j][k]+=0.5*dt/dx*(-3.0*D->Bx[i][j][k]+4.0*D->Bx[i+1][j][k]-D->Bx[i+2][j][k]);
      D->By[i][j][k]+=0.5*dt/dx*(-3.0*D->By[i][j][k]+4.0*D->By[i+1][j][k]-D->By[i+2][j][k])+dt/dx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
      D->Bz[i][j][k]+=0.5*dt/dx*(-3.0*D->Bz[i][j][k]+4.0*D->Bz[i+1][j][k]-D->Bz[i+2][j][k])-dt/dx*(D->Ey[i+1][j][k]-D->Ey[i][j][k]);
      D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
      D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
      D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
    }
}

void Esolve1D(Domain *D)
{
    int i,j,k,istart,iend,nxSub;
    double dx,dy,dz,dt;

    dx=D->dx;
    dt=D->dt;
    nxSub=D->nxSub;

    istart=D->istart;
    iend=D->iend;

    j=k=0;
    for(i=istart; i<iend; i++)
    {
      D->Ex[i][j][k]+=0.5*dt/dx*(-3.0*D->Ex[i][j][k]+4.0*D->Ex[i+1][j][k]-D->Ex[i+2][j][k])-dt*D->Jx[i][j][k];
      D->Ey[i][j][k]+=0.5*dt/dx*(-3.0*D->Ey[i][j][k]+4.0*D->Ey[i+1][j][k]-D->Ey[i+2][j][k])-dt/dx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-dt*D->Jy[i][j][k];
      D->Ez[i][j][k]+=0.5*dt/dx*(-3.0*D->Ez[i][j][k]+4.0*D->Ez[i+1][j][k]-D->Ez[i+2][j][k])+dt/dx*(D->By[i][j][k]-D->By[i-1][j][k])-dt*D->Jz[i][j][k];
    }
}



void solveLaser(Domain *D)
{
  int i,j,k,n,istart,iend;
  double dx,dt,a,b,c,d,e,alpha1,alpha2,alpha3,beta1,beta2,beta3,kp;
  double densityRatio;
  LaserList *L;

  dx=D->dx;
  dt=D->dt;
  kp=D->kp;
  istart=D->istart;
  iend=D->iend;

  j=k=0;

  L=D->laserList;
  while(L->next)  {
    for(i=istart; i<iend; i++)  {
      D->aOld2[i][0][0].real=D->aOld[i][j][k].real;
      D->aOld2[i][0][0].img=D->aOld[i][j][k].img;
      D->aOld[i][j][k].real=D->aNow[i][j][k].real;
      D->aOld[i][j][k].img=D->aNow[i][j][k].img;
      D->aNow[i][j][k].real=D->aNext[i][j][k].real;
      D->aNow[i][j][k].img=D->aNext[i][j][k].img;
    }
      
    for(i=istart; i<iend; i++)  {
      densityRatio=-dt*dt*0.5*D->Rho[i-1][j][k]/D->gamma[i-1][j][k];
      c=-dt/dx;
      a=1.0+c+densityRatio;
      d=-1.0+c-densityRatio;
      b=L->k0*dt;       
      e=a*a+b*b;

      alpha1=a*c/e;
      beta1=-b*c/e;
      alpha2=2.0*alpha1;
      beta2=-2.0*b/e;
      alpha3=(a*d+b*b)/e;
      beta3=(b*a-b*d)/e;

      D->aNext[i][j][k].real=(D->aNow[i][j][k].real-D->aOld2[i][j][k].real)*alpha1-(D->aNow[i][j][k].img-D->aOld2[i][j][k].img)*beta1+D->aNow[i-1][j][k].real*alpha2-D->aNow[i-1][j][k].img*beta2+D->aOld[i-1][j][k].real*alpha3-D->aOld[i-1][j][k].img*beta3;
      D->aNext[i][j][k].img=(D->aNow[i][j][k].real-D->aOld2[i][j][k].real)*beta1+(D->aNow[i][j][k].img-D->aOld2[i][j][k].img)*alpha1+D->aNow[i-1][j][k].real*beta2+D->aNow[i-1][j][k].img*alpha2+D->aOld[i-1][j][k].real*beta3+D->aOld[i-1][j][k].img*alpha3;
    }

    for(n=0; n<1; n++)  {
      for(i=istart; i<iend; i++)  {
        densityRatio=-dt*dt*0.5*D->Rho[i-1][j][k]/D->gamma[i-1][j][k];
        c=-dt/dx;
        a=1.0+c+densityRatio;
        d=-1.0+c-densityRatio;
        b=L->k0*dt;       
        e=a*a+b*b;

        alpha1=a*c/e;
        beta1=-b*c/e;
        alpha2=2.0*alpha1;
        beta2=-2.0*b/e;
        alpha3=(a*d+b*b)/e;
        beta3=(b*a-b*d)/e;

        D->aNext[i-1][j][k].real=(D->aNext[i][j][k].real-D->aOld[i][j][k].real)*alpha1-(D->aNext[i][j][k].img-D->aOld[i][j][k].img)*beta1+D->aNow[i-1][j][k].real*alpha2-D->aNow[i-1][j][k].img*beta2+D->aOld[i-1][j][k].real*alpha3-D->aOld[i-1][j][k].img*beta3;
        D->aNext[i-1][j][k].img=(D->aNext[i][j][k].real-D->aOld[i][j][k].real)*beta1+(D->aNext[i][j][k].img-D->aOld[i][j][k].img)*alpha1+D->aNow[i-1][j][k].real*beta2+D->aNow[i-1][j][k].img*alpha2+D->aOld[i-1][j][k].real*beta3+D->aOld[i-1][j][k].img*alpha3;
      }
    }
    L=L->next;
  }

}


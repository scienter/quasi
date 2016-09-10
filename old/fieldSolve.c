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
  double dx,dt,a,b,c,d,e,f,g,h,o;
  double alpha1,alpha2,alpha3,alpha4,beta1,beta2,beta3,beta4,kp;
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
    for(i=0; i<iend+3; i++)  {
      D->aOld3[i][0][0].real=D->aOld2[i][j][k].real;
      D->aOld3[i][0][0].img=D->aOld2[i][j][k].img;
      D->aOld2[i][0][0].real=D->aOld[i][j][k].real;
      D->aOld2[i][0][0].img=D->aOld[i][j][k].img;
      D->aOld[i][j][k].real=D->aNow[i][j][k].real;
      D->aOld[i][j][k].img=D->aNow[i][j][k].img;
      D->aNow[i][j][k].real=D->aNext[i][j][k].real;
      D->aNow[i][j][k].img=D->aNext[i][j][k].img;
    }
     
    for(i=istart; i<iend; i++)  {
//      i=iend;
      densityRatio=dt*dt*0.5*D->Rho[i+1][j][k]/D->gamma[i+1][j][k];
      c=dt/dx;
      a=0.5+c;
      b=-0.5*L->k0*dt;       
      d=-b;
      e=-densityRatio;
      f=1.0+e;
      g=c-0.5;
      h=b;
      o=a*a+b*b;

      alpha1=(a*c+b*d)/o;
      beta1=(a*d-b*c)/o;
      alpha2=a*e/o;
      beta2=-e*b/o;
      alpha3=a*f/o;
      beta3=-f*b/o;
      alpha4=(a*g+b*h)/o;
      beta4=(a*h-b*g)/o;

      D->aNext[i][j][k].real=(D->aOld[i+1][j][k].real-D->aOld3[i+1][j][k].real)*alpha1-(D->aOld[i+1][j][k].img-D->aOld3[i+1][j][k].img)*beta1
        +D->aNow[i+1][j][k].real*alpha2-D->aNow[i+1][j][k].img*beta2
        +D->aNow[i][j][k].real*alpha3-D->aNow[i][j][k].img*beta3
        +D->aOld[i][j][k].real*alpha4-D->aOld[i][j][k].img*beta4;
      D->aNext[i][j][k].img=(D->aOld[i+1][j][k].img-D->aOld3[i+1][j][k].img)*alpha1+(D->aOld[i+1][j][k].real-D->aOld3[i+1][j][k].real)*beta1
        +D->aNow[i+1][j][k].img*alpha2+D->aNow[i+1][j][k].real*beta2
        +D->aNow[i][j][k].img*alpha3+D->aNow[i][j][k].real*beta3
        +D->aOld[i][j][k].img*alpha4+D->aOld[i][j][k].real*beta4;
    }

    for(n=0; n<9; n++)  {
      for(i=istart; i<iend; i++)  {
//      for(i=iend; i>=istart; i--)  {
        densityRatio=dt*dt*0.5*D->Rho[i][j][k]/D->gamma[i][j][k];
        c=dt/dx;
        a=0.5+c;
        b=-0.5*L->k0*dt;       
        d=-b;
        e=-densityRatio;
        f=1.0+e;
        g=c-0.5;
        h=b;
        o=a*a+b*b;

        alpha1=(a*c+b*d)/o;
        beta1=(a*d-b*c)/o;
        alpha2=a*e/o;
        beta2=-e*b/o;
        alpha3=a*f/o;
        beta3=-f*b/o;
        alpha4=(a*g+b*h)/o;
        beta4=(a*h-b*g)/o;

        D->aNext[i][j][k].real=(D->aNext[i+1][j][k].real-D->aOld[i+1][j][k].real)*alpha1-(D->aNext[i+1][j][k].img-D->aOld[i+1][j][k].img)*beta1
          +D->aNow[i+1][j][k].real*alpha2-D->aNow[i+1][j][k].img*beta2
          +D->aNow[i][j][k].real*alpha3-D->aNow[i][j][k].img*beta3
          +D->aOld[i][j][k].real*alpha4-D->aOld[i][j][k].img*beta4;
        D->aNext[i][j][k].img=(D->aNext[i+1][j][k].img-D->aOld[i+1][j][k].img)*alpha1+(D->aNext[i+1][j][k].real-D->aOld[i+1][j][k].real)*beta1
          +D->aNow[i+1][j][k].img*alpha2+D->aNow[i+1][j][k].real*beta2
          +D->aNow[i][j][k].img*alpha3+D->aNow[i][j][k].real*beta3
          +D->aOld[i][j][k].img*alpha4+D->aOld[i][j][k].real*beta4;
      }
    }
/*
      for(i=istart; i<iend; i++)  {
        densityRatio=dt*dt*0.5*D->Rho[i-1][j][k]/D->gamma[i-1][j][k];
        c=dt/dx;
        a=0.5+c;
        b=-0.5*L->k0*dt;       
        d=-b;
        e=-densityRatio;
        f=1.0+e;
        g=c-0.5;
        h=b;
        o=a*a+b*b;

        alpha1=(a*c+b*d)/o;
        beta1=(a*d-b*c)/o;
        alpha2=a*e/o;
        beta2=-e*b/o;
        alpha3=a*f/o;
        beta3=-f*b/o;
        alpha4=(a*g+b*h)/o;
        beta4=(a*h-b*g)/o;

        D->aNext[i-1][j][k].real=(D->aNext[i][j][k].real-D->aOld[i][j][k].real)*alpha1-(D->aNext[i][j][k].img-D->aOld[i][j][k].img)*beta1
          +D->aNow[i][j][k].real*alpha2-D->aNow[i][j][k].img*beta2
          +D->aNow[i-1][j][k].real*alpha3-D->aNow[i-1][j][k].img*beta3
          +D->aOld[i-1][j][k].real*alpha4-D->aOld[i-1][j][k].img*beta4;
        D->aNext[i-1][j][k].img=(D->aNext[i][j][k].img-D->aOld[i][j][k].img)*alpha1+(D->aNext[i][j][k].real-D->aOld[i][j][k].real)*beta1
          +D->aNow[i][j][k].img*alpha2+D->aNow[i][j][k].real*beta2
          +D->aNow[i-1][j][k].img*alpha3+D->aNow[i-1][j][k].real*beta3
          +D->aOld[i-1][j][k].img*alpha4+D->aOld[i-1][j][k].real*beta4;
      }
*/
    L=L->next;
  }

}


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"


void solveLaser(Domain *D)
{
  int i,j,k,n,istart,iend;
  double dx,dt,a,b,c,d,e,alpha1,alpha2,alpha3,beta1,beta2,beta3,kp;
  LaserList *L;

  dx=D->dx;
  dt=D->dt;
  kp=D->kp;
  istart=D->istart;
  iend=D->iend;

  j=k=0;

  L=D->laserList;
//  while(L->next)  {
    for(i=istart; i<iend; i++)  {
      D->aOld2[i][0][0].real=D->aOld[i][j][k].real;
      D->aOld2[i][0][0].img=D->aOld[i][j][k].img;
      D->aOld[i][j][k].real=D->aNow[i][j][k].real;
      D->aOld[i][j][k].img=D->aNow[i][j][k].img;
      D->aNow[i][j][k].real=D->aNext[i][j][k].real;
      D->aNow[i][j][k].img=D->aNext[i][j][k].img;
    }
      
    c=-dt/dx;
    a=1.0+c;
    d=-1.0+c;
    b=L->k0*dt;       
    e=a*a+b*b;

    alpha1=a*c/e;
    beta1=-b*c/e;
    alpha2=2.0*alpha1;
    beta2=-2.0*b/e;
    alpha3=(a*d+b*b)/e;
    beta3=(b*a-b*d)/e;

    for(i=istart; i<iend; i++)  {
      D->aNext[i][j][k].real=(D->aNow[i][j][k].real-D->aOld2[i][j][k].real)*alpha1-(D->aNow[i][j][k].img-D->aOld2[i][j][k].img)*beta1+D->aNow[i-1][j][k].real*alpha2-D->aNow[i-1][j][k].img*beta2+D->aOld[i-1][j][k].real*alpha3-D->aOld[i-1][j][k].img*beta3;
      D->aNext[i][j][k].img=(D->aNow[i][j][k].real-D->aOld2[i][j][k].real)*beta1+(D->aNow[i][j][k].img-D->aOld2[i][j][k].img)*alpha1+D->aNow[i-1][j][k].real*beta2+D->aNow[i-1][j][k].img*alpha2+D->aOld[i-1][j][k].real*beta3+D->aOld[i-1][j][k].img*alpha3;
    }

    for(n=0; n<1; n++)  {
      for(i=istart; i<iend; i++)  {
        D->aNext[i-1][j][k].real=(D->aNext[i][j][k].real-D->aOld[i][j][k].real)*alpha1-(D->aNext[i][j][k].img-D->aOld[i][j][k].img)*beta1+D->aNow[i-1][j][k].real*alpha2-D->aNow[i-1][j][k].img*beta2+D->aOld[i-1][j][k].real*alpha3-D->aOld[i-1][j][k].img*beta3;
        D->aNext[i-1][j][k].img=(D->aNext[i][j][k].real-D->aOld[i][j][k].real)*beta1+(D->aNext[i][j][k].img-D->aOld[i][j][k].img)*alpha1+D->aNow[i-1][j][k].real*beta2+D->aNow[i-1][j][k].img*alpha2+D->aOld[i-1][j][k].real*beta3+D->aOld[i-1][j][k].img*alpha3;
      }
    }
//    L=L->next;
//  }

}


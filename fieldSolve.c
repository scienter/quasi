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

  L=D->laserList;
  while(L->next)  {
    for(i=istart; i<iend; i++)  {
      c=-dt/dx;
      a=1.0+c;
      d=-1.0+c;
      b=L->k0/kp*dt;       
      e=a*a+b*b;

      alpha1=a*c/e;
      beta1=-b*c/e;
      alpha2=2.0*alpha1;
      beta2=-2.0*b/e;
      alpha3=(a*d+b*b)/e;
      beta3=(b*a-b*d)/e;

      D->aNext[i-1][j][k].real=(D->aNow[i][j][k].real-D->aOld[i][j][k].real)*alpha1-(D->aNow[i][j][k].img-D->aOld[i][j][k].img)*beta1+D->aNow[i-1][j][k].real*alpha2-D->aNow[i-1][j][k].img*beta2+D->aOld[i-1][j][k].real*alpha3-D->aOld[i-1][j][k].img*beta3;
      D->aNext[i-1][j][k].img=(D->aNow[i][j][k].real-D->aOld[i][j][k].real)*beta1+(D->aNow[i][j][k].img-D->aOld[i][j][k].img)*alpha1+D->aNow[i-1][j][k].real*beta2+D->aNow[i-1][j][k].img*alpha2+D->aOld[i-1][j][k].real*beta3+D->aOld[i-1][j][k].img*alpha3;
      for(n=0; n<9; n++)  {
        D->aNext[i-1][j][k].real=(D->aNow[i][j][k].real-D->aOld[i][j][k].real)*alpha1-(D->aNow[i][j][k].img-D->aOld[i][j][k].img)*beta1+D->aNow[i-1][j][k].real*alpha2-D->aNow[i-1][j][k].img*beta2+D->aOld[i-1][j][k].real*alpha3-D->aOld[i-1][j][k].img*beta3;
        D->aNext[i-1][j][k].img=(D->aNow[i][j][k].real-D->aOld[i][j][k].real)*beta1+(D->aNow[i][j][k].img-D->aOld[i][j][k].img)*alpha1+D->aNow[i-1][j][k].real*beta2+D->aNow[i-1][j][k].img*alpha2+D->aOld[i-1][j][k].real*beta3+D->aOld[i-1][j][k].img*alpha3;
      }
    }
    L=L->next;
  }

}

/*
void fieldSolve(Domain *D)
{
  void solveField1D_DSX();
  void solveField2DC_DSX();
  void solveField2D_DSX();
  void solveField3DC_DSX();
  void solveField3D_DSX();
  void Bsolve2D_Yee();
  void Esolve2D_Yee();
  void solveField2DC_boostIon_DSX();
  void solveField2D_boostIon_DSX();

  int myrank, nTasks,rank,rankM,rankN;
  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch((D->fieldType-1)*3+D->dimension) {
  //1D field
  case (Split-1)*3+1:
    solveField1D_DSX(D);
    break;
  //split 2D
  case (Split-1)*3+2:
    MPI_TransferF_DSX_Yminus(D,D->Sr,D->Sl,D->Ex,D->Pr,D->Pl,D->Bx,D->nx+5,1,3);
    MPI_TransferF_DSX_Yplus(D,D->Pr,D->Pl,D->Bx,D->Sr,D->Sl,D->Ex,D->nx+5,1,3);

    if(D->boostIon==OFF && D->boostOn==ON)
      solveField2DC_boostIon_DSX(D);
    else 
      solveField2DC_DSX(D);
    if(D->pmlOn==ON)
    {
      if(myrank==nTasks-1)
        absorb2DC(D,UP);
      else	;
      if(myrank==0)
        absorb2DC(D,DOWN);
      else	;
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);
    //D->nx+5, 1 is number z, 3 is share.
    MPI_TransferF_DSX_YminusC(D,D->SrC,D->SlC,D->ExC,D->nx+5,1,3);
    MPI_TransferF_DSX_YplusC(D,D->PrC,D->PlC,D->BxC,D->nx+5,1,3);


    if(D->boostIon==OFF && D->boostOn==ON)
      solveField2D_boostIon_DSX(D);
    else 
      solveField2D_DSX(D);
    if(D->pmlOn==ON)
    {
      if(myrank==nTasks-1)
        absorb2D(D,UP);
      else	;
      if(myrank==0)
        absorb2D(D,DOWN);
      else	;
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);

    break;

  //3D
  case (Split-1)*3+3:
//    if(D->boostIon==OFF && D->boostOn==ON)
//      solveField2DC_boostIon_DSX(D);
//    else 
    if(D->M>1)
    {
      MPI_TransferF_DSX_Yminus(D,D->Sr,D->Sl,D->Ex,D->Pr,D->Pl,D->Bx,D->nx+5,D->nzSub+5,3);
      MPI_TransferF_DSX_Yplus(D,D->Pr,D->Pl,D->Bx,D->Sr,D->Sl,D->Ex,D->nx+5,D->nzSub+5,3);
    }
    else	;
    if(D->N>1)
    {
      MPI_TransferF_DSX_Zminus(D,D->Sr,D->Sl,D->Ex,D->Pr,D->Pl,D->Bx,D->nx+5,D->nySub+5,3);
      MPI_TransferF_DSX_Zplus(D,D->Pr,D->Pl,D->Bx,D->Sr,D->Sl,D->Ex,D->nx+5,D->nySub+5,3);
    }
    else	;
      solveField3DC_DSX(D);

    if(D->pmlOn==ON)
    {
      rankM=myrank%D->M;
      rankN=(int)(myrank/D->M);
      if(rankM==D->M-1)
        absorb3DC(D,UP);
      else	;
//      if(rankM==D->M-1 && rankN==D->N-1)
//        absorb3DC(D,UPFRONT);
//      if(rankM==D->M-1 && rankN==0)
//        absorb3DC(D,UPBACK);
      if(rankM==0)
        absorb3DC(D,DOWN);
      else	;
//      if(rankM==0 && rankN==D->N-1)
//        absorb3DC(D,DOWNFRONT);
//      if(rankM==0 && rankN==0)
//        absorb3DC(D,DOWNBACK);
      if(rankN==D->N-1)
        absorb3DC(D,FRONT);
      else	;
      if(rankN==0)
        absorb3DC(D,BACK);
      else	;
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);

    //3 is share. But it works on '..C'.
    if(D->M>1)
    {
      MPI_TransferF_DSX_YminusC(D,D->SrC,D->SlC,D->ExC,D->nx+5,D->nzSub+5,3);
      MPI_TransferF_DSX_YplusC(D,D->PrC,D->PlC,D->BxC,D->nx+5,D->nzSub+5,3);
    }
    else	;
    if(D->N>1)
    {
      MPI_TransferF_DSX_ZminusC(D,D->PrC,D->PlC,D->ExC,D->nx+5,D->nySub+5,3);
      MPI_TransferF_DSX_ZplusC(D,D->SrC,D->SlC,D->BxC,D->nx+5,D->nySub+5,3);
    }
    else	;


//    if(D->boostIon==OFF && D->boostOn==ON)
//      solveField2D_boostIon_DSX(D);
//    else 
      solveField3D_DSX(D);
    
    if(D->pmlOn==ON)
    {
      rankM=myrank%D->M;
      rankN=(int)(myrank/D->M);
      if(rankM==D->M-1)
        absorb3D(D,UP);
      else	;
//      if(rankM==D->M-1 && rankN==D->N-1)
//        absorb3D(D,UPFRONT);
//      if(rankM==D->M-1 && rankN==0)
//        absorb3D(D,UPBACK);
      if(rankM==0)
        absorb3D(D,DOWN);
      else	;
//      if(rankM==0 && rankN==D->N-1)
//        absorb3D(D,DOWNFRONT);
//      if(rankM==0 && rankN==0)
//        absorb3D(D,DOWNBACK);
      if(rankN==D->N-1)
        absorb3D(D,FRONT);
      else	;
      if(rankN==0)
        absorb3D(D,BACK);
      else	;
    }
    else	;
    MPI_Barrier(MPI_COMM_WORLD);

    break;

  //Yee 2D
  case (Yee-1)*3+2:
    Bsolve2D_Yee(D);

    Esolve2D_Yee(D);

    break;
  default:
    printf("what fieldType? and what dimension?\n");
  }
}


void solveField2DC_boostIon_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,preB1C;
    double nowPr,nowSr,prevPr,prevSr;
    double nowPrC,nowSrC,prevPrC,prevSrC;

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    // PrC,PlC,E1C,SrC,SlC,B1C
    k=0;
      for(j=jstart; j<jend; j++)
      {
        nowPrC=D->PrC[istart-1][j][k];
        nowSrC=D->SrC[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->ExC[i][j][k]+=dt/dy*(D->Pr[i][j][k]-D->Pr[i][j-1][k]-D->Pl[i][j][k]+D->Pl[i][j-1][k])-pi*dt*(D->JxOld[i][j][k]+D->Jx[i][j][k]-D->JxBoost[i][j][k]);
          D->BxC[i][j][k]+=-dt/dy*(D->Sr[i][j+1][k]-D->Sr[i][j][k]+D->Sl[i][j+1][k]-D->Sl[i][j][k]);
          prevPrC=nowPrC;
          nowPrC=D->PrC[i][j][k];
          D->PrC[i][j][k]=prevPrC+0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]-D->JyBoost[i][j][k]);
          D->PlC[i-1][j][k]=D->PlC[i][j][k]-0.25*dt/dy*(D->Ex[i][j+1][k]+D->Ex[i-1][j+1][k]-D->Ex[i][j][k]-D->Ex[i-1][j][k])-0.5*pi*dt*(D->JyOld[i][j][k]+D->Jy[i][j][k]-D->JyBoost[i][j][k]);
          prevSrC=nowSrC;
          nowSrC=D->SrC[i][j][k];
          D->SrC[i][j][k]=prevSrC-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]-D->JzBoost[i][j][k]);
          D->SlC[i-1][j][k]=D->SlC[i][j][k]-0.25*dt/dy*(D->Bx[i][j][k]+D->Bx[i-1][j][k]-D->Bx[i][j-1][k]-D->Bx[i-1][j-1][k])-0.5*pi*dt*(D->JzOld[i][j][k]+D->Jz[i][j][k]-D->JzBoost[i][j][k]);
        }	//End of i
      }		//End of j
}



void solveField2D_boostIon_DSX(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,preB1C;
    double nowPr,nowSr,prevPr,prevSr;
    double nowPrC,nowSrC,prevPrC,prevSrC;

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    // Pr,Pl,E1,Sr,Sl,B1
    k=0;
      for(j=jstart; j<jend; j++)
      {
        nowPr=D->Pr[istart-1][j][k];
        nowSr=D->Sr[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]+=dt/dy*(D->PrC[i][j][k]-D->PrC[i][j-1][k]-D->PlC[i][j][k]+D->PlC[i][j-1][k])-2*pi*dt*(D->Jx[i][j][k]-D->JxBoost[i][j][k]);
          D->Bx[i][j][k]+=-dt/dy*(D->SrC[i][j+1][k]-D->SrC[i][j][k]+D->SlC[i][j+1][k]-D->SlC[i][j][k]);
          prevPr=nowPr;
          nowPr=D->Pr[i][j][k];
          D->Pr[i][j][k]=prevPr+0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*(D->Jy[i][j][k]-D->JyBoost[i][j][k]);
          D->Pl[i-1][j][k]=D->Pl[i][j][k]-0.25*dt/dy*(D->ExC[i][j+1][k]+D->ExC[i-1][j+1][k]-D->ExC[i][j][k]-D->ExC[i-1][j][k])-pi*dt*(D->Jy[i][j][k]-D->JyBoost[i][j][k]);
          prevSr=nowSr;
          nowSr=D->Sr[i][j][k];
          D->Sr[i][j][k]=prevSr-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])-pi*dt*(D->Jz[i][j][k]-D->JzBoost[i][j][k]);
          D->Sl[i-1][j][k]=D->Sl[i][j][k]-0.25*dt/dy*(D->BxC[i][j][k]+D->BxC[i-1][j][k]-D->BxC[i][j-1][k]-D->BxC[i-1][j-1][k])-pi*dt*(D->Jz[i][j][k]-D->JzBoost[i][j][k]);
        }	//End of i
      }		//End of j
}

void Bsolve2D_Yee(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt,oldBx,oldBy,oldBz;

    dx=D->dx;
    dy=D->dy;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    k=0;
    //Solving B field
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        oldBx=D->Bx[i][j][k];
        oldBy=D->By[i][j][k];
        oldBz=D->Bz[i][j][k];
        D->Bx[i][j][k]+=-dt/dy*(D->Ez[i][j+1][k]-D->Ez[i][j][k]);
        D->By[i][j][k]+=dt/dx*(D->Ez[i+1][j][k]-D->Ez[i][j][k]);
        D->Bz[i][j][k]+=-dt/dx*(D->Ey[i+1][j][k]-D->Ey[i][j][k])+dt/dy*(D->Ex[i][j+1][k]-D->Ex[i][j][k]);
        D->BxNow[i][j][k]=0.5*(D->Bx[i][j][k]+oldBx);
        D->ByNow[i][j][k]=0.5*(D->By[i][j][k]+oldBy);
        D->BzNow[i][j][k]=0.5*(D->Bz[i][j][k]+oldBz);
      }
}

void Esolve2D_Yee(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,nxSub,nySub,nzSub;  
    double dx,dy,dz,dt;

    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dt=D->dt;
    nxSub=D->nxSub;
    nySub=D->nySub;
    nzSub=D->nzSub;
    
    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    k=0;
    //Solving E field
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        D->Ex[i][j][k]+=dt/dy*(D->Bz[i][j][k]-D->Bz[i][j-1][k])-2*pi*dt*D->Jx[i][j][k];
        D->Ey[i][j][k]+=-dt/dx*(D->Bz[i][j][k]-D->Bz[i-1][j][k])-2*pi*dt*D->Jy[i][j][k];
        D->Ez[i][j][k]+=dt/dx*(D->By[i][j][k]-D->By[i-1][j][k])-dt/dy*(D->Bx[i][j][k]-D->Bx[i][j-1][k])-2*pi*dt*D->Jz[i][j][k];
      }
}

void solveField1D_DSX(Domain *D)
{
    int i,j,k,istart,iend,nxSub;  
    double dx,dt;
    double nowPr,nowSr,prevPr,prevSr;
    int nTasks,myrank;
    MPI_Status status;          
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);    

    dx=D->dx;
    dt=D->dt;
    nxSub=D->nxSub;
    
    istart=D->istart;
    iend=D->iend;

    j=k=0;
        nowPr=D->Pr[istart-1][j][k];
        nowSr=D->Sr[istart-1][j][k];
        for(i=istart; i<iend; i++)
        {
          D->Ex[i][j][k]+=-2*pi*dt*D->Jx[i][j][k];
          prevPr=nowPr;
          nowPr=D->Pr[i][j][k];
          D->Pr[i][j][k]=prevPr-pi*dt*D->Jy[i][j][k];
          D->Pl[i-1][j][k]=D->Pl[i][j][k]-pi*dt*D->Jy[i][j][k];
          prevSr=nowSr;
          nowSr=D->Sr[i][j][k];
          D->Sr[i][j][k]=prevSr-pi*dt*D->Jz[i][j][k];
          D->Sl[i-1][j][k]=D->Sl[i][j][k]-pi*dt*D->Jz[i][j][k];
        }	//End of i
}
*/

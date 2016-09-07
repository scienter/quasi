#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void interpolation(Domain *D)
{
  void interpolation1D(Domain *D);

  switch(D->dimension)  {
  case 1 :
    interpolation1D(D);
    break;
  default :
    ;
  }
}

void interpolation1D(Domain *D)
{
   int i,j,k,i1,istart,iend,s;
   double E1,E2,E3,B1,B2,B3,aReal,aImg,a,x,x1,aNow,aNext;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;

   Particle ***particle;
   particle=D->particle;

   j=k=0;
   for(i=istart; i<iend; i++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x; // y=p->y; // z=p->z;
           i1=(int)(i+x+0.5);
           x1=x+0.5-((int)(x+0.5));

           E1=(1-x1)*D->Ex[i1-1][j][k]
             +    x1*D->Ex[i1][j][k];
           E2=(1-x)*D->Ey[i][j][k]
             +    x*D->Ey[i+1][j][k];
           E3=(1-x)*D->Ez[i][j][k]
             +    x*D->Ez[i+1][j][k];
           B1=0;
           B2=(1-x1)*D->By[i1-1][j][k]
             +    x1*D->By[i1][j][k];
           B3=(1-x1)*D->Bz[i1-1][j][k]
             +    x1*D->Bz[i1][j][k];
         //  aReal=0.5*(D->aNow[i][j][k].real+D->aOld[i][j][k].real);
//           aImg=0.5*(D->aNow[i][j][k].img+D->aOld[i][j][k].img);
           aReal=D->aNow[i][j][k].real;
           aImg=D->aNow[i][j][k].img;
           aNow=D->aNow[i][j][k].real*D->aNow[i][j][k].real+D->aNow[i][j][k].img*D->aNow[i][j][k].img;
           aNext=D->aNow[i+1][j][k].real*D->aNow[i+1][j][k].real+D->aNow[i+1][j][k].img*D->aNow[i+1][j][k].img;
           a=(1-x)*aNow+x*aNext;

           p->E1=E1; p->E2=E2; p->E3=E3;
           p->B1=B1; p->B2=B2; p->B3=B3;
           p->aa=a;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}



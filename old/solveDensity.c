#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"


void solveDensity(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s,ii,jj,kk,i1,j1,k1;
    int l,m,n;
    double x,y,z,Wx[3],Wy[3],Wz[3],dt,dx;
    char name[100];
    Particle ***particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    dx=D->dx;
    dt=D->dt;

    double rho0[D->nSpecies];
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       rho0[s]=LL->charge*LL->density/D->density/((double)LL->numberInCell);
       LL=LL->next;
       s++;
    }

    switch (D->dimension)  {
    case 1 :
      j=k=0;
      for(i=0; i<iend+3; i++)
        D->Rho[i][j][k]=0.0;
      
      //initializing density
      for(i=istart; i<iend; i++)
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            x=p->x;
            i1=(int)(i+x+0.5);
            x=i+x-i1;
            Wx[0]=0.5*(0.5-x)*(0.5-x);
            Wx[1]=0.75-x*x;
            Wx[2]=0.5*(x+0.5)*(x+0.5);

            for(ii=0; ii<3; ii++)
            {
              l=i1-1+ii;
              if(istart<=l && l<iend)
                D->Rho[l][j][k]+=Wx[ii]*rho0[s]*dt/dx;
            }
            p=p->next;
          }
        }
      break;
/*
    case 2 :
      k=0;
      //initializing density
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
          D->Rho[i][j][k]=0.0;

      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(s=0; s<D->nSpecies; s++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              x=p->x; y=p->y; z=p->z;
              i1=(int)(i+x+0.5);
              j1=(int)(j+y+0.5);
              x=i+x-i1;
              y=j+y-j1;
              Wx[0]=0.5*(0.5-x)*(0.5-x);
              Wx[1]=0.75-x*x;
              Wx[2]=0.5*(x+0.5)*(x+0.5);
              Wy[0]=0.5*(0.5-y)*(0.5-y);
              Wy[1]=0.75-y*y;
              Wy[2]=0.5*(y+0.5)*(y+0.5);

              for(jj=0; jj<3; jj++)
                for(ii=0; ii<3; ii++)
                {
                  l=i1-1+ii;
                  m=j1-1+jj;
                  if(istart<=l && l<iend && jstart<=m && m<jend)
                    D->Rho[l][m][k]+=Wx[ii]*Wy[jj]*rho0[s];
                }
              p=p->next;
            }
          }

      sprintf(name,"rho%d_%d",iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart-1; i<=iend; i++)
      {
        for(j=jstart-1; j<=jend; j++)
        {
          x=(i-istart+D->minXSub)*D->dx*D->lambda;
          y=(j-jstart+D->minYSub)*D->dy*D->lambda;
          fprintf(out,"%g %g %g\n",x,y,D->Rho[i][j][k]);    
        }           
        fprintf(out,"\n");    
      }
      fclose(out);

      break;
    case 3 :
      //initializing density
      for(i=0; i<=iend; i++)
        for(j=jstart-1; j<=jend; j++)
          for(k=kstart-1; k<=kend; k++)
            D->Rho[i][j][k]=0.0;

      for(i=istart-1; i<iend-1; i++)
        for(j=jstart-1; j<jend-1; j++)
          for(k=kstart-1; k<=kend; k++)
            for(s=0; s<D->nSpecies; s++)
            {
              p=particle[i][j][k].head[s]->pt;
              while(p)
              {
                x=p->x; y=p->y; z=p->z;
                i1=(int)(i+x+0.5);
                j1=(int)(j+y+0.5);
                k1=(int)(k+z+0.5);
                x=i+x-i1;
                y=j+y-j1;
                z=k+z-k1;
                Wx[0]=0.5*(0.5-x)*(0.5-x);
                Wx[1]=0.75-x*x;
                Wx[2]=0.5*(x+0.5)*(x+0.5);
                Wy[0]=0.5*(0.5-y)*(0.5-y);
                Wy[1]=0.75-y*y;
                Wy[2]=0.5*(y+0.5)*(y+0.5);
                Wz[0]=0.5*(0.5-z)*(0.5-z);
                Wz[1]=0.75-z*z;
                Wz[2]=0.5*(z+0.5)*(z+0.5);

                for(ii=0; ii<3; ii++)
                  for(jj=0; jj<3; jj++)
                    for(kk=0; kk<3; kk++)
                      D->Rho[i1-1+ii][j1-1+jj][k1-1+kk]
                            +=Wx[ii]*Wy[jj]*Wz[kk]*rho0[s];
                p=p->next;
              }
            }

      sprintf(name,"rho%d_%d",iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart-1; i<=iend; i++)
      {
        for(j=jstart-1; j<=jend; j++)
        {
          for(k=kstart-1; k<=kend; k++)
          {
            x=(i-istart+D->minXSub)*D->dx*D->lambda;
            y=(j-jstart+D->minYSub)*D->dy*D->lambda;
            z=(k-kstart+D->minZSub)*D->dz*D->lambda;
            fprintf(out,"%g %g %g %g\n",x,y,z,D->Rho[i][j][k]);    
          }           
          fprintf(out,"\n");    
        }           
        fprintf(out,"\n");    
      }
      fclose(out);
      break;
*/
    }

}

#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void saveFile(Domain *D,int iteration)
{
  int nTasks,myrank;
  void saveLaser(Domain *D,int iteration);

  MPI_Comm_size(MPI_COMM_WORLD, &nTasks); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 

  if(iteration%D->saveStep==0 && iteration>=D->saveStart)
  {
    if(D->fieldSave==ON) { 
      if(D->saveMode==TXT)         saveLaser(D,iteration);
      else        ;
      if(myrank==0)
        printf("laser%d is made.\n",iteration);
    }
    else  ;

/*
          if(D.particleSave==ON) { 
            saveParticle(&D,iteration);
//            if(D.saveMode==TXT)               saveParticle(&D,iteration);
//            else if(D.saveMode==HDF)  saveParticleHDF(D,iteration,D.nSpecies);
//            else ;
            if(myrank==0)
              printf("particle%d is made.\n",iteration); 
          }
          else  ;

          if(D.rhoSave==ON) { 
            if(D.saveMode==TXT)         saveRho(&D,iteration);
            else if(D.saveMode==HDF)    saveDensityHDF(D,iteration,D.nSpecies);
            else                        ;
            if(myrank==0)
              printf("rho%d is made.\n",iteration);    
          }
          else  ;


          if(D.dumpSave==ON && iteration>=D.dumpStart) { 
              saveDump(D,iteration);
//            if(D.saveMode==TXT)       saveDump(D,iteration);
//            else if(D.saveMode==HDF)  saveDumpHDF(D,iteration);
//            else      ;
            if(myrank==0)
              printf("dump%d is made.\n",iteration);              
          }

          if(D.probeNum>0) { 
//            saveProbe(&D,iteration);
//            if(myrank==0)
//              printf("probe%d is made.\n",iteration);              
          }
*/




  }

}

/*
void saveProbe(Domain *D,int iteration)
{
    int i,j,n;
    char name[100];
    double t,Ex,Ey,Ez,Bx,By,Bz,Pr,Pl,Sr,Sl,x,y;
    double omega,frequency,dt;
    FILE *out;
    int myrank, nprocs;    
    Probe **probe;
    probe=D->probe;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    omega=2*pi*velocityC/D->lambda;
    frequency=omega/2.0/pi;   
    dt=1.0/frequency*D->dt;

    for(n=0; n<D->probeNum; n++)
    {
      if(D->probeX[n]>=D->minXSub && D->probeX[n]<D->maxXSub && 
         D->probeY[n]>=D->minYSub && D->probeY[n]<D->maxYSub)
      {
        x=D->probeX[n]*D->dx*D->lambda;
        y=D->probeY[n]*D->dy*D->lambda;
        sprintf(name,"probeRaman%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Pr=probe[n][i].Pr;
          Pl=probe[n][i].Pl;
          Sr=probe[n][i].Sr;
          Sl=probe[n][i].Sl;
          fprintf(out,"%g %g %g %g %g %g %g\n",t,Pr,Pl,Sr,Sl,x,y);
        }
        fclose(out);

        sprintf(name,"probe%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Ex=probe[n][i].E1;
          Bx=probe[n][i].B1;
          Ey=probe[n][i].Pr+probe[n][i].Pl;
          Ez=probe[n][i].Sr+probe[n][i].Sl;
          By=probe[n][i].Sl-probe[n][i].Sr;
          Bz=probe[n][i].Pr-probe[n][i].Pl;
          fprintf(out,"%g %g %g %g %g %g %g %g %g\n",t,Ex,Ey,Ez,Bx,By,Bz,x,y);
        }             
        fclose(out);
      }
    }
}
*/

/*
void saveRho(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,s,ii,jj,kk,i1,j1,k1;
    int l,m,n;
    double x,y,z,Wx[3],Wy[3],Wz[3];
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

    double rho0[D->nSpecies];
    s=0;
    LL=D->loadList;
    while(LL->next)
    {
       rho0[s]=LL->charge*LL->density/((double)LL->numberInCell);
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
                    D->Rho[l][j][k]+=Wx[ii]*rho0[s];
                }
              p=p->next;
            }
          }

      sprintf(name,"rho%d_%d",iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart-1; i<=iend; i++)
      {
          x=(i-istart+D->minXSub)*D->dx*D->lambda;
          fprintf(out,"%g %g\n",x,D->Rho[i][j][k]);    
      }
      fclose(out);

      break;
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
    }

}
*/
/*
void boostSaveField(Domain *D,int labSaveStep)
{
    int i,j,istart,iend,jstart,jend,show;
    char name[100];
    double x,y,e1,pr,pl,b1,sr,sl;
    double factor,dx;
    FILE *out;
    int myrank, nprocs;    
    Boost **boost;
    boost=D->boost;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"bField%d_%d",labSaveStep,myrank);
    out = fopen(name,"w");

    for(i=istart; i<iend; i++)
    {
      for(j=jstart; j<jend; j++)
      {
        x=boost[i][j].x;
        y=boost[i][j].y;
        e1=boost[i][j].E1;
        pr=boost[i][j].Pr;
        pl=boost[i][j].Pl;
        b1=boost[i][j].B1;
        sr=boost[i][j].Sr;
        sl=boost[i][j].Sl;
        if(x>0) {
          show=1;
          fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,e1,pr,pl,b1,sr,sl);
        }
        else show=0;
      }  
      if(show==1)           
        fprintf(out,"\n");
    }
    fclose(out);
    
    if(myrank==0)
      printf("bField%d is saved.\n",labSaveStep);
}
*/

/*
void saveField(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    char name[100];
    double x,y,z,Ex,Ey,Ez,Bx,By,Bz,factor;
    FILE *out,*out1;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    sprintf(name,"field%d_%d",iteration,myrank);
    out = fopen(name,"w");
    factor=D->gamma*(1+D->beta);

    switch((D->fieldType-1)*3+D->dimension) {
    case (Split-1)*3+1:
      j=k=0;
      for(i=istart; i<iend; i++)
      {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          Ex=D->Ex[i][j][k];    
          Ey=D->Pr[i][j][k]+D->Pl[i][j][k];
          Ez=D->Sr[i][j][k]+D->Sl[i][j][k];
          Bx=0.0;    
          By=D->Sl[i][j][k]-D->Sr[i][j][k];
          Bz=D->Pr[i][j][k]-D->Pl[i][j][k];
          fprintf(out,"%g %g %g %g %g %g %g\n",x,Ex,Ey,Ez,Bx,By,Bz);
      }
      fclose(out);
      break;
    case (Split-1)*3+2:
      k=0;
      for(i=istart; i<iend; i++)
      {
        for(j=jstart; j<jend; j++)
        {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          y=(j-2+D->minYSub)*D->dy*D->lambda;
          Ex=D->Ex[i][j][k];    
          Ey=D->Pr[i][j][k]+D->Pl[i][j][k];
          Ez=D->Sr[i][j][k]+D->Sl[i][j][k];
          Bx=D->Bx[i][j][k];    
          By=D->Sl[i][j][k]-D->Sr[i][j][k];
          Bz=D->Pr[i][j][k]-D->Pl[i][j][k];
          fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,Ex,Ey,Ez,Bx,By,Bz);
        }
        fprintf(out,"\n");                 
      }
      fclose(out);
      break;
    case (Split-1)*3+3:
//      k=(int)((kstart+kend)*0.5);
//      j=(int)((jstart+jend)*0.5);
      for(i=istart; i<iend; i++)
      {
        for(j=jstart; j<jend; j++)
        {
          for(k=kstart; k<kend; k++)
          {
            x=(i-2+D->minXSub)*D->dx*D->lambda;
            y=(j-2+D->minYSub)*D->dy*D->lambda;
            z=(k-2+D->minZSub)*D->dz*D->lambda;
            Ex=D->Ex[i][j][k];    
            Ey=D->Pr[i][j][k]+D->Pl[i][j][k];
            Ez=D->Sr[i][j][k]+D->Sl[i][j][k];
            Bx=D->Bx[i][j][k];    
            By=D->Sl[i][j][k]-D->Sr[i][j][k];
            Bz=D->Pr[i][j][k]-D->Pl[i][j][k];
            fprintf(out,"%g %g %g %g %g %g %g %g %g\n",x,y,z,Ex,Ey,Ez,Bx,By,Bz);
          }
          fprintf(out,"\n");                 
        }
        fprintf(out,"\n");                 
      }
      fclose(out);
      break;
    case (Yee-1)*3+2 :
      k=0;
      i=(int)(0.5*(istart+iend-1));
//      for(i=istart; i<iend; i++)
//      {
        for(j=jstart; j<jend; j++)
        {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          y=(j-2+D->minYSub)*D->dy*D->lambda;
          Ex=D->Ex[i][j][k];    
          Ey=D->Ey[i][j][k];
          Ez=D->Ez[i][j][k];
          Bx=D->Bx[i][j][k];    
          By=D->By[i][j][k];    
          Bz=D->Bz[i][j][k];    
          fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,Ex,Ey,Ez,Bx,By,Bz);
        }
        fprintf(out,"\n");                 
//      }
      fclose(out);
      break;
    default :
      printf("what field_type? and what dimension?\n");
    }
}
*/

void saveLaser(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    char name[100];
    double x,y,z,real,img;
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;
    iend=D->iend;

    sprintf(name,"laser%d_%d",iteration,myrank);
    out = fopen(name,"w");

    switch(D->dimension) {
    case 1 :
      j=k=0;
      for(i=istart; i<iend; i++)
      {
        x=(i-2+D->minXSub)*D->dx/D->kp;
        real=D->aNow[i][j][k].real;
        img=D->aNow[i][j][k].img;
        fprintf(out,"%g %g %g\n",x,real,img);
      }
      fclose(out);
      break;
    default :
      ;
    }
}

/*
void saveParticle(Domain *D,int iteration)
{
  int i,j,k,istart,iend,jstart,jend,kstart,kend,s,core,index;
  char name[100];
  double x,y,z,p1,p2,p3,gamma,mc;
  double Pr,Pl,E1,Sr,Sl,B1;
  double minPx[D->nSpecies];
  Particle ***particle;
  particle=D->particle;
  ptclList *p;
  LoadList *LL;
  FILE *out;
  int myrank, nprocs;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;

  LL=D->loadList;
  s=0;
  while(LL->next)
  {
    minPx[s]=LL->givenMinPx;
    LL=LL->next;
    s++;
  }

  switch (D->dimension)
  {
  case 1:
    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
      out = fopen(name,"w");    
      j=k=0;
      for(i=istart; i<iend; i++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            x=((i-istart+D->minXSub)+p->x)*D->dx*D->lambda; 
            p1=p->p1;
            p2=p->p2;    
            p3=p->p3;
            gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
            index=p->index;
            fprintf(out,"%g %g %g %g %g %d\n",x,p1,p2,p3,gamma,index);               
//              fprintf(out,"%g %g %g %g %g\n",x,y,p->E1,p->E2,p->E3);               
            p=p->next;
          }	//End of while(p)
        }
    }				//End of for(s)
    break;
  case 2:
    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
      out = fopen(name,"w");    
      k=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
        {
          p=particle[i][j][k].head[s]->pt;
          while(p)
          {
            x=((i-istart+D->minXSub)+p->x)*D->dx*D->lambda; 
            y=((j-jstart+D->minYSub)+p->y)*D->dy*D->lambda; 
            p1=p->p1;
            p2=p->p2;    
            p3=p->p3;
            gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
            index=p->index;
            core=p->core;
            if(p1>=minPx[s])
              fprintf(out,"%g %g %g %g %g %g %9d %d\n",x,y,p1,p2,p3,gamma,index,core);               
//              fprintf(out,"%g %g %g %g %g\n",x,y,p->E1,p->E2,p->E3);               
            p=p->next;
          }	//End of while(p)
        }	//End of for(i,j)
      fclose(out);
    }				//End of for(s)
    break;
  case 3:
    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
          for(k=kstart; k<kend; k++)
          {
            p=particle[i][j][k].head[s]->pt;
            while(p)
            {
              x=((i-istart+D->minXSub)+p->x)*D->dx*D->lambda; 
              y=((j-jstart+D->minYSub)+p->y)*D->dy*D->lambda; 
              z=((k-kstart+D->minZSub)+p->z)*D->dz*D->lambda; 
              p1=p->p1;
              p2=p->p2;    
              p3=p->p3;
              gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
              index=p->index;
              core=p->core;
              if(p1>=minPx[s])
                fprintf(out,"%g %g %g %g %g %g %g %9d %d\n",x,y,z,p1,p2,p3,gamma,index,core);               
//              fprintf(out,"%g %g %g %g %g %g\n",x,y,z,p->E1,p->E2,p->E3);               
              p=p->next;
            }	//End of while(p)
          }	//End of for(i,j,k)
      fclose(out);
    }				//End of for(s)
    break;
  default:
    ;
  }
}
*/

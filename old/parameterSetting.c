#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>


double randomV()
{
   double r;
   int intRand, randRange=1000, rangeDev;

   intRand = rand() % randRange;
   r = ((double)intRand)/randRange;

   return r;
}

void parameterSetting(Domain *D,char *input)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   FILE *out;
   int FindParameters();
//   int findLoadParameters();
//   int findLaserParameters();
   int whatONOFF();
   int whatSaveMode();
   int whatFieldType();
   double minX,maxX,minY,maxY,minZ,maxZ;
   double positionX,factor,pMinX,pMaxX,pPosition;
   double normalB,normalE,Ex,Ey,Ez,Bx,By,Bz,dyoverdx,dzoverdx;
   char str[100],name[100];
   int rank,minT,maxT,tmpInt,fail=0,cnt;
   int i,j,k,n,numProbeX,numProbeY,numProbeZ,probeType;
   double lambda,tmpDouble,probeDx,probeDy,probeDz,maxProbeX,minProbeX,maxProbeY,minProbeY,maxProbeZ,minProbeZ,tmp;
   LoadList *LL, *New;
   LaserList *L, *LNew;


   //initially
   if(FindParameters("Domain",1,"dimension",input,str)) D->dimension=atoi(str);
   else  {
      printf("in [Domain], dimension=?  (1:1D, 2:2D, 3:3D)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_mode",input,str))
     D->saveMode=whatSaveMode(str);
   else
     D->saveMode=TXT;
   if(FindParameters("Domain",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Domain",1,"max_step",input,str)) D->maxStep=atoi(str);
   else  {
      printf("In [Domain], maxStep=? \n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_step",input,str)) D->saveStep=atoi(str);
   else  {
      printf("In [Domain], save_step=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  {
      printf("In [Domain], save_start=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dump_save",input,str))
     D->dumpSave=whatONOFF(str);
   else
     D->dumpSave=OFF;
   if(FindParameters("Domain",1,"dump_start",input,str))
     D->dumpStart=atoi(str);
   else
     D->dumpStart=D->saveStart;
   D->dumpStep=0;
   if(FindParameters("Domain",1,"field_save",input,str))
     D->fieldSave=whatONOFF(str);
   else
     D->fieldSave=ON;
   if(FindParameters("Domain",1,"particle_save",input,str))
     D->particleSave=whatONOFF(str);
   else
     D->particleSave=ON;
   if(FindParameters("Domain",1,"density_save",input,str))
     D->rhoSave=whatONOFF(str);
   else
     D->rhoSave=ON;
   if(FindParameters("Domain",1,"current_save",input,str))
     D->currentSave=whatONOFF(str);
   else
     D->currentSave=ON;

   //domain
   if(FindParameters("Domain",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Domain",1,"max_step",input,str)) D->maxStep=atoi(str);
   else  {
      printf("In [Domain], maxStep=? \n");
      fail=1;
   }
   if(FindParameters("Domain",1,"minX",input,str))
      D->domainMinX=atof(str);
   else  {
      printf("In [Domain], minX=? [m].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxX",input,str))
      D->domainMaxX=atof(str);
   else  {
      printf("In [Domain], maxX=? [m].\n");
      fail=1;
   }

   if(FindParameters("Domain",1,"density",input,str))
      D->density=atof(str);
   else  {
      printf("In [Domain], density=? [m^-3]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"division",input,str))
      D->division=atof(str);
   else  {
      printf("In [Domain], division=? [number of devision of lambda_p]\n");
      fail=1;
   }

   //additional domain parameters
   D->wp=sqrt(D->density*eCharge*eCharge/eps0/eMass);
   D->kp=D->wp/velocityC;
   D->dx=1.0/D->division;
   D->domainMaxX*=D->kp;
   D->domainMinX*=D->kp;
   D->domainX=D->domainMaxX-D->domainMinX;
   D->nx=((int)D->domainX*D->division);
   D->dt=D->dx;

   while(D->nx*D->dx<D->domainX)
     D->nx+=1;

   D->minXSub=D->domainMinX*D->division;
printf("nx=%d, nx*dx=%g\n",D->nx,D->nx*D->dx/D->kp);

   //Laser parameter setting
   D->laserList = (LaserList *)malloc(sizeof(LaserList));
   D->laserList->next = NULL;
   L = D->laserList;
   rank = 1;
   while(findLaserParameters(rank,L,D,input))
   {
      LNew = (LaserList *)malloc(sizeof(LaserList));
      LNew->next = NULL;
      L->next=LNew;
      L=L->next;
      rank ++;
   }
   D->nLaser = rank-1;

   //Plasma parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
   LL = D->loadList;
   rank = 1;
   while(findLoadParameters(rank, LL, D,input))
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   if(fail==1)
     exit(0);
   else	;
}


int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   int FindParameters();
   char name[100], str[100];
   int fail=0,mode=0;
   double lambda;

   if(FindParameters("Laser",rank,"mode",input,str)) mode=atoi(str);

   if(mode)
   {
     if(FindParameters("Laser",rank,"wavelength",input,str)) 
        lambda=atof(str);
     else  {
        printf("in [Laser], lambda=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"a0",input,str)) 
        L->amplitude=atof(str);
     else  {
        printf("in [Laser], a0=??\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
     else  {
        printf("in [Laser], rU=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
     else  {
        printf("in [Laser], rD=? [# of basic wavelength]\n");
        fail=1;
     }
     if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
     else  L->flat=0.0;

     //additional laser parameters
     L->mode=mode;
     L->k0=2*pi/lambda/D->kp;
//     L->rayleighLength=pi/(lambda/D->gamma/(1.0+D->beta))*L->beamWaist*L->beamWaist/lambda;
//     L->beamWaist=L->beamWaist/D->lambda;
//     L->focus=L->focus/D->lambda;
     L->rU*=lambda*D->kp;
     L->rD*=lambda*D->kp;
     if(fail==1)
        exit(0);
   }
   return mode;
}


int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   int FindParameters();
   LoadList *New;
   int whatSpecies();
   int whatPlasmaType();
   double whatMass(int species);
   double pointPosition,wp,pDt;
   int whatCharge();
   int whatFunctionMode();
   int whatDefineMode();
   char name[100], str[100];
   int i,n,cnt,species,fail=0;
   double tmp,max,min;
   double *shareDouble;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Plasma",rank,"type",input,name)) 
   {
     LL->type = whatPlasmaType(name);
//     if(D->boostOn==ON)
//       LL->type = BoostFrame;
   }
   else LL->type=0;

   if(LL->type>0)
   {
      if(FindParameters("Plasma",rank,"density",input,str)) 
         LL->density=atof(str);
      else  {
         printf("in [Plasma], density=? [m-3]\n");
         fail=1;
      }

      if(FindParameters("Plasma",rank,"species",input,name)) 
         species = whatSpecies(name);
      else  species = 0;
      LL->species=species;
      if(FindParameters("Plasma",rank,"numberInCell",input,str)) 
         LL->numberInCell=atoi(str);
      else  {
         printf("in [Plasma], numberInCell=? \n");
         fail=1;
      }
      if(FindParameters("Plasma",rank,"startIndex",input,str)) 
         LL->index=atoi(str);
      else  
         LL->index=0;
      if(FindParameters("Plasma",rank,"temperature",input,str))  
         LL->temperature=atof(str);
      else   LL->temperature=0.0;	
      LL->mass=whatMass(species);
      LL->charge=whatCharge(species);
//      LL->superP=LL->density*D->lambda*D->dx*D->lambda*D->dy*D->lambda*D->dz/LL->numberInCell;

      srand(1*(myrank+1));
      switch (LL->type)  {
      case Polygon :
        if(FindParameters("Plasma",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
        else  {
          printf("in [Plasma], Xnodes=?\n");
          printf("Each nodes indicates the point of plasma density changing.\n");
          fail=1;
        }
        if(LL->xnodes>0)
        {
          LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
          LL->xn = (double *)malloc(LL->xnodes*sizeof(double));   
          for(i=0; i<LL->xnodes; i++)
          {
            sprintf(name,"X%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xpoint[i] = atof(str)*D->kp*D->division;
            else 
            { printf("X%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Xn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xn[i] = atof(str);
            else 
            { printf("Xn%d should be defined.\n",i);  fail=1; } 
          }
        }
/*
        if(D->dimension>1)
        {
          if(FindParameters("Plasma",rank,"Ynodes",input,str)) LL->ynodes=atoi(str);
          else  {
            printf("in [Plasma], Ynodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          if(LL->ynodes>0)
          {
            LL->ypoint = (double *)malloc(LL->ynodes*sizeof(double));
            LL->yn = (double *)malloc(LL->ynodes*sizeof(double));   
            for(i=0; i<LL->ynodes; i++)
            {
              sprintf(name,"Y%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) {
                LL->ypoint[i] = atof(str)/D->lambda/D->dy;
              }
              else 
              { printf("Y%d should be defined.\n",i);  fail=1; }
 
              sprintf(name,"Yn%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) 
                LL->yn[i] = atof(str);
              else 
              { printf("Yn%d should be defined.\n",i);  fail=1; } 
            }
          }
        }
        if(D->dimension>2)
        {
          if(FindParameters("Plasma",rank,"Znodes",input,str)) LL->znodes=atoi(str);
          else  {
            printf("in [Plasma], Znodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          LL->zpoint = (double *)malloc(LL->znodes*sizeof(double));
          LL->zn = (double *)malloc(LL->znodes*sizeof(double));   
          for(i=0; i<LL->znodes; i++)
          {
            sprintf(name,"Z%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) {
              LL->zpoint[i] = atof(str)/D->lambda/D->dz;
            }
            else 
            { printf("Z%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Zn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->zn[i] = atof(str);
            else 
            { printf("Zn%d should be defined.\n",i);  fail=1; } 
          }
        }

        if(FindParameters("Plasma",rank,"centerX",input,str))  
          LL->centerX=atof(str)/D->lambda/D->dx;
        else   LL->centerX=0.0;	
        if(FindParameters("Plasma",rank,"centerZ",input,str))  
          LL->centerZ=atof(str)/D->lambda/D->dz;
        else   LL->centerZ=0.0;	
        if(FindParameters("Plasma",rank,"gaussCoefX",input,str))  
          LL->gaussCoefX=atof(str)/D->lambda/D->dx;
        else   LL->gaussCoefX=1.0;
        if(FindParameters("Plasma",rank,"polyCoefX2",input,str))  
          LL->polyCoefX2=atof(str)/D->lambda/D->dx;
        else   LL->polyCoefX2=0.0;	
        if(FindParameters("Plasma",rank,"functionModeX",input,str))  
          LL->modeX=whatFunctionMode(str);
        else   LL->modeX=0;	
        if(FindParameters("Plasma",rank,"functionModeYZ",input,str))  
          LL->modeYZ=whatFunctionMode(str);
        else   LL->modeYZ=0;	
        if(D->dimension>1)
        {	
          if(FindParameters("Plasma",rank,"centerY",input,str))  
            LL->centerY=atof(str)/D->lambda/D->dy;
          else   LL->centerY=0.0;	
          if(FindParameters("Plasma",rank,"gaussCoefYZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/D->dy;
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"polyCoefYZ2",input,str))  
            LL->polyCoefYZ2=atof(str)/D->lambda/D->dy;
          else   LL->polyCoefYZ2=0.0;	
        }
        else if(D->dimension>2)
        {	
          if(FindParameters("Plasma",rank,"centerY",input,str))  
            LL->centerZ=atof(str)/D->lambda/D->dz;
          else   LL->centerZ=0.0;	
          if(FindParameters("Plasma",rank,"gaussCoefYZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/sqrt(D->dy*D->dy+D->dz*D->dz);
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"polyCoefYZ2",input,str))  
            LL->polyCoefYZ2=atof(str)/D->lambda/sqrt(D->dy*D->dy+D->dz*D->dz);
          else   LL->polyCoefYZ2=0.0;	
        }
        else 	;
*/
        break;
/*
      case Defined :
//        srand(time(NULL)*(myrank+1));
        if(FindParameters("Plasma",rank,"define_mode",input,str))  
          LL->defineMode=whatDefineMode(str);
        else   LL->defineMode=byNumber;	
        if(FindParameters("Plasma",rank,"number_defined",input,str))  
          LL->numDefined=atoi(str);
        else   LL->numDefined=0;	
        if(LL->defineMode==byDensity)
        {
          if(FindParameters("Plasma",rank,"minX",input,str))  
            LL->minX=atof(str);
          else   {   printf("minX=?  [m]\n");  exit(0);   }	
          if(FindParameters("Plasma",rank,"maxX",input,str))  
            LL->maxX=atof(str);
          else   {   printf("maxX=?  [m]\n");  exit(0);   }	
          if(D->dimension>1)
          {
            if(FindParameters("Plasma",rank,"minY",input,str))  
              LL->minY=atof(str);
            else   {   printf("minY=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxY",input,str))  
              LL->maxY=atof(str);
            else   {   printf("maxY=?  [m]\n");  exit(0);   }	
          }
          if(D->dimension>2)
          {
            if(FindParameters("Plasma",rank,"minZ",input,str))  
              LL->minZ=atof(str);
            else   {   printf("minZ=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxZ",input,str))  
              LL->maxZ=atof(str);
            else   {   printf("maxZ=?  [m]\n");  exit(0);   }	
          }
        }
        else 	;

        if(FindParameters("Plasma",rank,"xlength_particle",input,str))  
          LL->xLengDef=atof(str)/D->lambda;
        else   LL->xLengDef=0.0;
        if(LL->numDefined>0)	
        {
          LL->xPosition=(double *)malloc(LL->numDefined*sizeof(double));
          shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
          for(i=0; i<LL->numDefined; i++)
          {
            if(LL->defineMode==byNumber)
            {
              sprintf(name,"xPosition%d",i);
              if(FindParameters("Plasma",rank,name,input,str))  
                LL->xPosition[i]=atof(str)/D->lambda;
              else
              { printf("xPosition%d should be defined.\n",i); fail=1;}
            }
            else if(LL->defineMode==byDensity)
            {
              if(myrank==0)
                shareDouble[i]=(LL->minX+randomV()*(LL->maxX-LL->minX))/D->lambda;
              else 	;
            }
            else 	;
          }
          if(LL->defineMode==byDensity)
          {
            MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            for(i=0; i<LL->numDefined; i++)
              LL->xPosition[i]=shareDouble[i];
          }
          else	;
          free(shareDouble);
        }
        if(D->dimension>1)
        {
          if(FindParameters("Plasma",rank,"ylength_particle",input,str))  
            LL->yLengDef=atof(str)/D->lambda;
          else   LL->yLengDef=0.0;	
          if(LL->numDefined>0)
          {
            LL->yPosition=(double *)malloc(LL->numDefined*sizeof(double));
            shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
            for(i=0; i<LL->numDefined; i++)
            {
              if(LL->defineMode==byNumber)
              {
                sprintf(name,"yPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str))  
                  LL->yPosition[i]=atof(str)/D->lambda;
                else
                { printf("yPosition%d should be defined.\n",i); fail=1;}
              }
              else if(LL->defineMode==byDensity)
              {
                if(myrank==0)
                  shareDouble[i]=(LL->minY+randomV()*(LL->maxY-LL->minY))/D->lambda;
                else 	;
              }
              else 	;
            }
            if(LL->defineMode==byDensity)
            {
              MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
              for(i=0; i<LL->numDefined; i++)
                LL->yPosition[i]=shareDouble[i];
            }
            else	;
            free(shareDouble);
          }
        }		//End of demension>1
        if(D->dimension>2)
        {
          if(FindParameters("Plasma",rank,"zlength_particle",input,str))  
            LL->zLengDef=atof(str)/D->lambda;
          else   LL->zLengDef=0.0;	
          if(LL->numDefined>0)
          {
            LL->zPosition=(double *)malloc(LL->numDefined*sizeof(double));
            shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
            for(i=0; i<LL->numDefined; i++)
            {
              if(LL->defineMode==byNumber)
              {
                sprintf(name,"zPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str))  
                  LL->zPosition[i]=atof(str)/D->lambda;
                else
                { printf("zPosition%d should be defined.\n",i); fail=1;}
              }
              else if(LL->defineMode==byDensity)
              {
                if(myrank==0)
                  shareDouble[i]=(LL->minZ+randomV()*(LL->maxZ-LL->minZ))/D->lambda;
                else 	;
              }
              else 	;
            }
            if(LL->defineMode==byDensity)
            {
              MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
              for(i=0; i<LL->numDefined; i++)
                LL->zPosition[i]=shareDouble[i];
            }
            else	;
            free(shareDouble);  
          }
        }
        if(D->dimension==2)
        {
          LL->numDefPtcls=(int)(LL->xLengDef*LL->yLengDef/D->dx/D->dy*LL->numberInCell);
          cnt=LL->numDefPtcls;
          if(LL->numDefined>0)
          {
            LL->define=(double **)malloc(LL->numDefined*sizeof(double *));
            for(i=0; i<LL->numDefined; i++)
              LL->define[i]=(double *)malloc((cnt*2)*sizeof(double ));
 
            max=0;
            min=(double)D->maxStep;
            for(i=0; i<LL->numDefined; i++)
              for(n=0; n<cnt; n++)
              {
                tmp=(double)(randomV());
                LL->define[i][n]=LL->xPosition[i]-0.5*LL->xLengDef+tmp*LL->xLengDef;
                if(LL->define[i][n]>max) max=LL->define[i][n];
                if(LL->define[i][n]<min) min=LL->define[i][n];
                tmp=(double)(randomV());
                LL->define[i][n+cnt]=LL->yPosition[i]-0.5*LL->yLengDef+tmp*LL->yLengDef;
              }
          }
        }  
        else if(D->dimension==3)
        {
          LL->numDefPtcls=(int)(LL->xLengDef*LL->yLengDef*LL->zLengDef/D->dx/D->dy/D->dz*LL->numberInCell);
          cnt=LL->numDefPtcls;
          if(LL->numDefined>0)
          {
            LL->define=(double **)malloc(LL->numDefined*sizeof(double *));
            for(i=0; i<LL->numDefined; i++)
              LL->define[i]=(double *)malloc((LL->numDefPtcls*3)*sizeof(double ));
            max=0;
            min=(double)D->maxStep;
            for(i=0; i<LL->numDefined; i++)
              for(n=0; n<cnt; n++)
              {
                tmp=(double)(randomV());
                LL->define[i][n]=LL->xPosition[i]-0.5*LL->xLengDef+tmp*LL->xLengDef;
                if(LL->define[i][n]>max) max=LL->define[i][n];
                if(LL->define[i][n]<min) min=LL->define[i][n];
                tmp=(double)(randomV());
                LL->define[i][n+cnt]=LL->yPosition[i]-0.5*LL->yLengDef+tmp*LL->yLengDef;
                tmp=(double)(randomV());
                LL->define[i][n+2*cnt]=LL->zPosition[i]-0.5*LL->zLengDef+tmp*LL->zLengDef;
              }
          }
        } 	//End fo dimension=3 : defined Plasma 
        LL->maxLoadTime=(int)((max+1)*D->divisionLambda);
        LL->minLoadTime=(int)((min-1)*D->divisionLambda);
        break;
*/
      }
   
   }	//end of if(species)

   if(fail==1)
      exit(0);

   return LL->type;
}


int whatDefineMode(char *str)
{
   if(strstr(str,"by_number")) 		return byNumber;
   else if(strstr(str,"by_density"))   	return byDensity;
   else 				return byNumber;
}

int whatONOFF(char *str)
{
   if(strstr(str,"ON")) 		return ON;
   else if(strstr(str,"OFF"))   	return OFF;
   else 				return OFF;
}

int whatSaveMode(char *str)
{
   if(strstr(str,"TXT")) 		return TXT;
   else if(strstr(str,"HDF"))   	return HDF;
   else 				return TXT;
}

int whatFieldType(char *str)
{
   if(strstr(str,"Split")) 		return Split;
   else if(strstr(str,"Yee"))   	return Yee;
   else 				return 0;
}

int whatSpecies(char *str)
{
   if(strstr(str,"Electron")) 		return Electron;
   else if(strstr(str,"HPlus0"))   	return HPlus0;
   else if(strstr(str,"HPlus1"))   	return HPlus1;
   else if(strstr(str,"HePlus0"))   	return HePlus0;
   else if(strstr(str,"HePlus1"))   	return HePlus1;
   else if(strstr(str,"HePlus2"))   	return HePlus1;
   else if(strstr(str,"CPlus0"))   	return CPlus0;
   else if(strstr(str,"CPlus1"))   	return CPlus1;
   else if(strstr(str,"CPlus2"))   	return CPlus2;
   else if(strstr(str,"CPlus3"))   	return CPlus3;
   else if(strstr(str,"CPlus4"))   	return CPlus4;
   else if(strstr(str,"CPlus5"))   	return CPlus5;
   else if(strstr(str,"CPlus6"))   	return CPlus6;
   else if(strstr(str,"AlPlus11"))   	return AlPlus11;
   else return 0;
}

int whatPlasmaType(char *str)
{
   if(strstr(str,"Polygon"))         return Polygon;
   else if(strstr(str,"Defined"))   	return Defined;
   else
   {
     printf("No Plasma type!\n"); 
     exit(0);
   }
   return 0;
}

double whatMass(int species)
{
   if(species == Electron) 		return 1;
   else if(species == HPlus0)  		return 1.00794/eMassU;
   else if(species == HPlus1)  		return (1.00794-1*eMassU)/eMassU;
   else if(species == HePlus0)  	return (4.00260-0*eMassU)/eMassU;
   else if(species == HePlus1)  	return (4.00260-1*eMassU)/eMassU;
   else if(species == HePlus2)  	return (4.00260-2*eMassU)/eMassU;
   else if(species == CPlus0)           return (12.0111-0*eMassU)/eMassU;
   else if(species == CPlus1)           return (12.0111-1*eMassU)/eMassU;
   else if(species == CPlus2)           return (12.0111-2*eMassU)/eMassU;
   else if(species == CPlus3)           return (12.0111-3*eMassU)/eMassU;
   else if(species == CPlus4)           return (12.0111-4*eMassU)/eMassU;
   else if(species == CPlus5)           return (12.0111-5*eMassU)/eMassU;
   else if(species == CPlus6)           return (12.0111-6*eMassU)/eMassU;
   else if(species == AlPlus11)         return (26.9815-11*eMassU)/eMassU;
   else {  printf("Species' mass not defined\n");  exit(0);  }
}

int whatCharge(int species)
{
   int fail;

   if(species == Electron) 		return -1;
   else if(species == HPlus0)  		return 0;
   else if(species == HPlus1)  		return 1;
   else if(species == HePlus0)  	return 0;
   else if(species == HePlus1)  	return 1;
   else if(species == HePlus2)  	return 2;
   else if(species == CPlus0)           return 0;
   else if(species == CPlus1)           return 1;
   else if(species == CPlus2)           return 2;
   else if(species == CPlus3)           return 3;
   else if(species == CPlus4)           return 4;
   else if(species == CPlus5)           return 5;
   else if(species == CPlus6)           return 6;
   else if(species == AlPlus11)           return 11;
   else {  printf("Species' charge not defined\n");  exit(0);  }
}

int whatFunctionMode(char *str)
{
   if(strstr(str,"Constant")) 		return Constant;
   else if(strstr(str,"Gaussian"))   	return Gaussian;
   else if(strstr(str,"Polynomial"))   	return Polynomial;
   else return 0;
}


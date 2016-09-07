#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include <time.h>

int main(int argc, char *argv[])
{
    int i,j,k,n,s,iteration=0,filter,boost,filterStep,labSaveStep;
    int rnk,suddenDump=OFF;
    double factor,time_spent;
    clock_t begin,end;
    double t;
    char name[100];
    FILE *out;
    Domain D;  
    LaserList *L;
    LoadList *LL;
    External Ext;
    UPML UPml;
    DPML DPml;
    int myrank, nTasks;
    MPI_Status status; 

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    begin=clock();

    if(argc < 2) 
    {  
      printf("mpirun -np N show [inputFile] [dumpNum]\n"); 
      exit(0); 
    }
    if(FindParameters("Domain",1,"filter",argv[1],name)) filter=atoi(name);
    else  filter=0;
    if(FindParameters("Domain",1,"filterStep",argv[1],name)) filterStep=atoi(name);
    else  filterStep=10;

    //parameter setting
    parameterSetting(&D,argv[1]);
    MPI_Barrier(MPI_COMM_WORLD);

    //create mesh
    boundary(&D);
    MPI_Barrier(MPI_COMM_WORLD);

    //load laser
    L=D.laserList;
    while(L->next)  {
      loadLaser(&D,L,t); 
//         if(L->direction==1)     loadLaser2D(&D,L,t); 
//         else if(L->direction==-1)     loadLaserOpp2D(&D,L,t); 
      L=L->next;
    }
       
/*
    //load plasma or load dump file
    if(argc >= 3)
    {   
      
      D.dumpStep = atoi(argv[2]);
      iteration=D.dumpStep;
      MPI_Barrier(MPI_COMM_WORLD);
      restoreData(&D,iteration);
//      if(D.saveMode==TXT)
//        restoreData(&D,iteration);
//      else if(D.saveMode==HDF)
//        restoreDumpHDF(&D,iteration);
//      else 	;

      t=D.dt*iteration;
    }
    else 
    {
      LL=D.loadList;
      s=0;
      while(LL->next)
      {
        loadPlasma(&D,LL,s,iteration);
        LL=LL->next;
        s++;
      }
      t=0;
    }
*/
    //rooping time 
    while(iteration<=D.maxStep)
    {
       //save File      
       saveFile(&D,iteration);

//       fieldSolve(&D);
       interpolation(&D);
       particlePush(&D);
//       updateCurrent(&D);


       LL=D.loadList;
       s=0;
       while(LL->next)   {
//         loadPlasma(&D,LL,s,iteration);
         LL=LL->next;
         s++;
       }
       rearrangeParticles(&D);
       removeEdge(&D);

       solveDensity(&D);
       solveLaser(&D);

       //time update
       if(iteration%10==0 && myrank==0)  
          printf("iteration = %d\n",iteration);           
       iteration+=1;
       t=D.dt*iteration;  
       D.minXSub++;
    }     //end of time roop                  

//    if(D.tracking==ON)
//      saveTracking(&D);

    end=clock();
    time_spent=(end-begin)/CLOCKS_PER_SEC;

    //make 'report' file
    if(myrank==0)
    {
      sprintf(name,"report");
      out = fopen(name,"w");
      fprintf(out,"nx=%d\n",D.nx);
      fprintf(out,"ny=%d\n",D.ny);
      fprintf(out,"nz=%d\n",D.nz);
      fprintf(out,"cores=%d\n",nTasks);
      fprintf(out,"running time=%gm\n",time_spent/60.0);
      fprintf(out,"dx=%g,kp=%g\n",D.dx/D.kp,D.kp);
      fclose(out);
    }
    else	;

//    cleanMemory(&D);
    
    MPI_Finalize();

    return 0;
}

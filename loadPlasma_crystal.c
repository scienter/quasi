#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>

double maxwellianVelocity(double temterature); 
void loadMovingPlasma_crystal(Domain *D,LoadList *LL,int s)
{
   int i,j,k,l,intNum,cnt,np,istart,iend;
   double space,position,positionX,x,ne,minX;
   double leftIndex,rightIndex,nc;
   double v1,v2,v3,gamma,mass;
   Particle ***particle;
   particle=D->particle;
   ptclList *New,*p;         

   istart=D->istart;
   iend=D->iend;

   i=iend-1;
   j=k=0;
   while(LL->next)
   {  
     for(l=0; l<LL->xnodes-1; l++)
     {
       position=i+D->minXSub;
 
       if(position>=LL->xpoint[l] && position<LL->xpoint[l+1])
       {
         nc=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])
            *(position-LL->xpoint[l])+LL->xn[l]);
         nc*=LL->numberInCell;	//it is the double number of superparticles.
         space=1.0/nc;
         p=particle[iend-1][j][k].head[s]->pt;
         x=-2;
         while(p)
         {
           if(p->x>x)  x=p->x;
           p=p->next;
         }
   
         while(x<1)
         {
           x+=space;
           positionX=x;

           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;

           New->x = positionX;   New->oldX=i+positionX;
           New->E1=New->E2=New->E3=New->B1=New->B2=New->B3=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
           New->p1=v1;
           New->p2=v2;
           New->p3=v3;
           LL->index++;
           New->index=LL->index;            
         }
       }
     }
     LL=LL->next;
     s++;
   }		//End of while(LL)  
}

double maxwellianVelocity(double temperature)
{
   double vth,r,prob,v,random;
   int intRand,randRange=1e5;

   vth=sqrt(2.0*eCharge*temperature/eMass);

   r=1.0;
   prob=0.0;
   while (r>prob)  {
      intRand = rand() % randRange;
      r = ((double)intRand)/randRange;
      intRand = rand() % randRange;
      random = ((double)intRand)/randRange;
      v = 6.0*(random-0.5);
      prob=exp(-v*v);
   }
   return vth*v;
}


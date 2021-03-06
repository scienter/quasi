#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>

void particlePush(Domain *D)
{
  void particlePush1D();
  void particlePush2D();
  void particlePush3D();

  switch(D->dimension)  {
  case 1 :
    particlePush1D(D);
    break;
  case 2 :
//    particlePush2D(D);
    break;
  case 3 :
//    particlePush3D(D);
    break;
  default :
    printf("In particlePush, what dimension(%d)?\n",D->dimension);
  }

}


void particlePush1D(Domain *D)
{
    int i,j,k,istart,iend,l,m,s,shift,cnt;
    double x,shiftX,dt,dx,gamma,gamma0,oldGamma,fluidGamma,coef,k1,k2;
    double A1,A2,A3,a,b;
    double Ex,Ey,Ez,Bx,By,Bz,p1,p2,p3,aa,nextA,presA,devAX;
    Particle ***particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p, *New, *tmp, *prev;

    istart=D->istart;
    iend=D->iend;

    dt=D->dt;
    dx=D->dx;

    double mass[D->nSpecies];
    int charge[D->nSpecies];
    LL=D->loadList;
    s=0;
    while(LL->next)
    {
       mass[s]=LL->mass;
       charge[s]=LL->charge;
       s++;
       LL=LL->next;
    }

    shiftX=0;

    j=k=0;
    for(i=istart; i<iend; i++) 
    { 
      A1=D->aNow[i][j][k].real*D->aNow[i][j][k].real+D->aNow[i][j][k].img*D->aNow[i][j][k].img;
      A2=D->aNow[i+1][j][k].real*D->aNow[i+1][j][k].real+D->aNow[i+1][j][k].img*D->aNow[i+1][j][k].img;
      A3=D->aNow[i+2][j][k].real*D->aNow[i+2][j][k].real+D->aNow[i+2][j][k].img*D->aNow[i+2][j][k].img;
      a=(A1-A2)+0.5*(A3-A1);
      b=2.0*(A2-A1)+0.5*(A1-A3);

      cnt=0; 
      fluidGamma=0.0; 
        for(s=0; s<D->nSpecies; s++)
        {
          coef=charge[s]/mass[s]*dt;
          p=particle[i][j][k].head[s]->pt;    
          while(p)
          {   
            oldGamma=p->gamma;
            devAX=(2.0*a*p->x+b)*0.25; 
            Ex=p->E1;
            Ey=p->E2;
            Ez=p->E3;
            Bx=p->B1;
            By=p->B2;
            Bz=p->B3;
            p1=p->p1;
            p2=p->p2;
            p3=p->p3;
            aa=p->aa;
            gamma0=sqrt(1.0+p1*p1+p2*p2+p3*p3+0.5*aa);
            
            //Calculate vector px
            k1=coef*(Ex+p2*Bz-p3*By+devAX)/gamma0;
            gamma=sqrt(1.0+(p1+0.5*k1)*(p1+0.5*k1)+p2*p2+p3*p3+0.5*aa);
            k2=coef*(Ex+p2*Bz-p3*By+devAX)/gamma;
            p->p1+=k2;

            //Calculate vector py
            k1=coef*(Ey-p1*Bz+p3*Bx)/gamma0;
            gamma=sqrt(1.0+p1*p1+(p2+0.5*k1)*(p2+0.5*k1)+p3*p3+0.5*aa);
            k2=coef*(Ey-p1*Bz+p3*Bx)/gamma;
            p->p2+=k2;

            //Calculate vector pz
            k1=coef*(Ez+p1*By-p2*Bx)/gamma0;
            gamma=sqrt(1.0+p1*p1+p2*p2+(p3+0.5*k1)*(p3+0.5*k1)+0.5*aa);
            k2=coef*(Ez+p1*By-p2*Bx)/gamma;
            p->p3+=k2;

            p->gamma=gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3+0.5*aa);
            //Translation
            shiftX=p->p1/gamma;    //dt is ignored because of dx=dt=1 in cell.
             //dt is ignored because of dx=dt=1 in cell.
            if(shiftX<=-1.0 || shiftX>=1.0)  {
              printf("particle's movement exceeds C velocity\n");
              printf("i=%d,shiftX=%g\n",i,shiftX);
              exit(0);
            } 
            p->oldX=i+p->x;
            p->x+=(shiftX-1.0)*dt/dx;

            cnt++;
            fluidGamma+=0.5*(oldGamma+gamma);
            p=p->next;
          }		//End of while(p)
        }		//End of for(s)

        if(cnt==0)  
          D->gamma[i][j][k]=1.0;
        else 
          D->gamma[i][j][k]=fluidGamma/((double)cnt);
    }
}

/*
void particlePush2D(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,l,m,s,shift,cnt;
    double x,dt,dx,dy,dz,sqrT,coef,dxOverdy,dxOverdz;
    double shiftX,shiftY,shiftZ,gamma;
    double E1,Pr,Pl,B1,Sr,Sl;
    double pMinus[3],T[3],S[3],operate[3][3],pPlus[3];
    Particle ***particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p, *New, *tmp, *prev;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    for(i=0; i<3; i++)   {
       pMinus[i]=0.0;
       T[i]=0.0;
       S[i]=0.0;
       pPlus[i]=0.0;
    }

    dt=D->dt;
    dx=D->dx;
    dy=D->dy;
    dxOverdy=D->dx/D->dy;

    double mass[D->nSpecies];
    int charge[D->nSpecies];
    LL=D->loadList;
    s=0;

    while(LL->next)
    {
       mass[s]=LL->mass;
       charge[s]=LL->charge;
       s++;
       LL=LL->next;
    }

    shiftX=shiftY=shiftZ=0.0;

    k=0;
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(s=0; s<D->nSpecies; s++)
        {
          coef=pi*charge[s]/mass[s]*dt;
          p=particle[i][j][k].head[s]->pt;     
          while(p)
          {    
            //Calculate vector P- 
            pMinus[0]=p->p1+coef*(p->E1);
            pMinus[1]=p->p2+coef*(p->E2);    
            pMinus[2]=p->p3+coef*(p->E3);

            //Calculate vector T 
            gamma=sqrt(1.0+pMinus[0]*pMinus[0]+pMinus[1]*pMinus[1]+pMinus[2]*pMinus[2]);           
            T[0]=coef/gamma*(p->B1);   
            T[1]=coef/gamma*(p->B2);
            T[2]=coef/gamma*(p->B3);

            //Calculate vector S
            sqrT=1.0+T[0]*T[0]+T[1]*T[1]+T[2]*T[2];
            for(l=0; l<3; l++)  
              S[l]=2.0*T[l]/sqrT;
  
            //Calculate operator A from P+=A.P-
            operate[0][0]=1.0-S[2]*T[2]-S[1]*T[1];      
            operate[0][1]=S[1]*T[0]+S[2];    
            operate[0][2]=S[2]*T[0]-S[1];     
            operate[1][0]=S[0]*T[1]-S[2];        
            operate[1][1]=1.0-S[0]*T[0]-S[2]*T[2];          
            operate[1][2]=S[2]*T[1]+S[0];         
            operate[2][0]=S[0]*T[2]+S[1];    
            operate[2][1]=S[1]*T[2]-S[0];    
            operate[2][2]=1-S[0]*T[0]-S[1]*T[1]; 
            //Calculate vector P+
            for(l=0; l<3; l++)  {
              pPlus[l]=0.0;
              for(m=0; m<3; m++)   
                pPlus[l]+=operate[l][m]*pMinus[m];
            }
            //Updated momentum              
            p->p1=pPlus[0]+coef*(p->E1); 
            p->p2=pPlus[1]+coef*(p->E2);    
            p->p3=pPlus[2]+coef*(p->E3); 
    
            //Translation
            gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
            shiftX=p->p1/gamma;    //dt is ignored because of dx=dt=1 in cell.
             //dt is ignored because of dx=dt=1 in cell.
            shiftY=p->p2/gamma*dxOverdy;    
            if(shiftX>=1 || shiftY>=1*dxOverdy)  {
              printf("particle's movement exceeds C velocity\n");
              printf("i=%d,j=%d,k=%d,shiftX=%g,shiftY=%g,shiftZ=%g\n",i,j,k,shiftX,shiftY,shiftZ);
              exit(0);
            } 
            p->oldX=i+p->x;
            p->x+=shiftX;
            p->oldY=j+p->y;
            p->y+=shiftY;
            p=p->next;
          }		//End of while(p)
        }		//End of for(s)
      }      	//End of for(i,j)

}



void particlePush3D(Domain *D)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend,l,m,s,shift,cnt;
    double x,shiftX,shiftY,shiftZ,dt,dx,dy,dz,gamma,sqrT,coef,dxOverdy,dxOverdz;
    double E1,Pr,Pl,B1,Sr,Sl;
    double pMinus[3],T[3],S[3],operate[3][3],pPlus[3];
    Particle ***particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p, *New, *tmp, *prev;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;
    kstart=D->kstart;
    kend=D->kend;

    for(i=0; i<3; i++)   {
       pMinus[i]=0.0;
       T[i]=0.0;
       S[i]=0.0;
       pPlus[i]=0.0;
    }

    dt=D->dt;
    dx=D->dx;
    dy=D->dy;
    dz=D->dz;
    dxOverdy=D->dx/D->dy;
    dxOverdz=D->dx/D->dz;

    double mass[D->nSpecies];
    int charge[D->nSpecies];
    LL=D->loadList;
    s=0;

    while(LL->next)
    {
       mass[s]=LL->mass;
       charge[s]=LL->charge;
       s++;
       LL=LL->next;
    }

    shiftX=shiftY=shiftZ=0.0;

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        for(k=kstart; k<kend; k++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            coef=pi*charge[s]/mass[s]*dt;
            p=particle[i][j][k].head[s]->pt;     
            while(p)
            {    
              //Calculate vector P- 
              pMinus[0]=p->p1+coef*(p->E1);
              pMinus[1]=p->p2+coef*(p->E2);    
              pMinus[2]=p->p3+coef*(p->E3);
 
             //Calculate vector T 
             gamma=sqrt(1.0+pMinus[0]*pMinus[0]+pMinus[1]*pMinus[1]+pMinus[2]*pMinus[2]);           
             T[0]=coef/gamma*(p->B1);   
             T[1]=coef/gamma*(p->B2);
             T[2]=coef/gamma*(p->B3);

             //Calculate vector S
             sqrT=1.0+T[0]*T[0]+T[1]*T[1]+T[2]*T[2];
             for(l=0; l<3; l++)  
                S[l]=2.0*T[l]/sqrT;
  
             //Calculate operator A from P+=A.P-
             operate[0][0]=1.0-S[2]*T[2]-S[1]*T[1];      
             operate[0][1]=S[1]*T[0]+S[2];    
             operate[0][2]=S[2]*T[0]-S[1];     
             operate[1][0]=S[0]*T[1]-S[2];        
             operate[1][1]=1.0-S[0]*T[0]-S[2]*T[2];          
             operate[1][2]=S[2]*T[1]+S[0];         
             operate[2][0]=S[0]*T[2]+S[1];    
             operate[2][1]=S[1]*T[2]-S[0];    
             operate[2][2]=1-S[0]*T[0]-S[1]*T[1]; 
             //Calculate vector P+
             for(l=0; l<3; l++)  {
                pPlus[l]=0.0;
                for(m=0; m<3; m++)   
                   pPlus[l]+=operate[l][m]*pMinus[m];
                }
             //Updated momentum              
             p->p1=pPlus[0]+coef*(p->E1); 
             p->p2=pPlus[1]+coef*(p->E2);    
             p->p3=pPlus[2]+coef*(p->E3); 
    
             //Translation
             gamma=sqrt(1.0+p->p1*p->p1+p->p2*p->p2+p->p3*p->p3);
             shiftX=p->p1/gamma;    //dt is ignored because of dx=dt=1 in cell.
             //dt is ignored because of dx=dt=1 in cell.
             shiftY=p->p2/gamma*dxOverdy;    
             shiftZ=p->p3/gamma*dxOverdz;  
             if(shiftX>=1 || shiftY>=1*dxOverdy || shiftZ>=1*dxOverdz)  {
                printf("particle's movement exceeds C velocity\n");
                printf("i=%d,j=%d,k=%d,shiftX=%g,shiftY=%g,shiftZ=%g\n",i,j,k,shiftX,shiftY,shiftZ);
                exit(0);
             } 
             p->oldX=i+p->x;
             p->x+=shiftX;
             p->oldY=j+p->y;
             p->y+=shiftY;
             p->oldZ=k+p->z;
             p->z+=shiftZ;
             p=p->next;
          }		//End of while(p)
        }		//End of for(s)
      }      	//End of for(i,j)

}
*/
            
double ponderoPushX(ptclList *p,double coef,double devAX)
{
   double a,px,py,pz,temp,k1,k2;
   
   a=p->aa;
   px=p->p1;
   py=p->p2;
   pz=p->p3;
   temp=1.0+a*0.5+px*px+py*py+pz*pz;
   k1=1.0/sqrt(temp)*coef*devAX;
   temp=1.0+a*0.5+(px+0.5*k1)*(px+0.5*k1)+py*py+pz*pz;
   k2=1.0/sqrt(temp)*coef*devAX;
   return k2;
} 

double ponderoPushY(ptclList *p,double coef,double devAY)
{
   double a,px,py,pz,temp,k1,k2;
   
   a=p->aa;
   px=p->p1;
   py=p->p2;
   pz=p->p3;
   temp=1.0+a*0.5+px*px+py*py+pz*pz;
   k1=1.0/sqrt(temp)*coef*devAY;
   temp=1.0+a*0.5+px*px+(py+0.5*k1)*(py+0.5*k1)+pz*pz;
   k2=1.0/sqrt(temp)*coef*devAY;
   return k2;
} 

double ponderoPushZ(ptclList *p,double coef,double devAZ)
{
   double a,px,py,pz,temp,k1,k2;
   
   a=p->aa;
   px=p->p1;
   py=p->p2;
   pz=p->p3;
   temp=1.0+a*0.5+px*px+py*py+pz*pz;
   k1=1.0/sqrt(temp)*coef*devAZ;
   temp=1.0+a*0.5+px*px+py*py+(pz+0.5*k1)*(pz+0.5*k1);
   k2=1.0/sqrt(temp)*coef*devAZ;
   return k2;
} 

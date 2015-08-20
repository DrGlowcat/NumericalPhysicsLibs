//library for Ordinary Differential Equations

#include "prototypes.h"

void PoorEulerODE(double t,double *Y,void (*dYdt)(double,double *,double *),double dt,const int N)
{
  double R[N];
  dYdt(t,Y,R);
  for(int i=0;i<N;i++)Y[i]+=(dt*R[i]);
}

void RungeKutta4(double t,double *Y,void (*dYdt)(double,double *,double *),double h,int N)
{
  double R[N];
  double Y1[N],Y2[N],Y3[N];
  double k1[N],k2[N],k3[N],k4[N];
  
  dYdt(t,Y,k1);
  for(int i=0;i<N;i++)
    {
      Y1[i]=Y[i]+0.5*h*k1[i];
    }

  dYdt(t+h*0.5,Y1,k2);
  for(int i=0;i<N;i++)
    {
      Y2[i]=Y[i]+0.5*h*k2[i];
    }

      dYdt(t+0.5*h,Y2,k3);
   for(int i=0;i<N;i++)
     {
      Y3[i]=Y[i]+0.5*h*k3[i];
     }
      
   dYdt(t+0.5*h,Y3,k4);
  for(int i=0;i<N;i++)
    {
        Y[i]+=(h/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
    }
}

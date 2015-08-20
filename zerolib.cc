//libreria di algoritmi che trovano gli zeri di una funzione
#include "prototypes.h"

int Bracket(double (*funz)(double x),double &inf,double &sup) //& indica che la funz crea un alias per inf e sup --> li modifica davvero
{
   double N=100.; //numero step
   double finf,fsup;
   double delta=(sup-inf)/N;
   finf=funz(inf);
   for(int i=0;i<=N;i++)
    {
      fsup=funz(inf+delta);
      if(finf*fsup<0)return 0;
      else
	{
	  finf=fsup;
	  inf+=delta;
	  sup=inf+delta;
	}
    }//end for
   return 1; //essendoci return 0 nel ciclo questo dovrebbe attivarsi solo se la funzione arriva alla fine del ciclo senza trovare lo zero
}

//algoritmo di bisezione
double Bisect(double (*fz)(double y),double a,double b,double precis,int &counter)
{
  double fa=fz(a);
  if(fa==0.)return a;
  double fb=fz(b);
  if(fb==0.)return b;
  double c,fc;
do
  {
    c=0.5*(a+b);
    fc=fz(c);
    if(fa*fc<0)
      {
	b=c;
	fb=fc;
      }
    if(fb*fc<0)
      {
	a=c;
	fa=fc;
      }
    counter++;
  }
 while(fabs(b-a)>1.e-9&&fabs(fc)>precis);

  return c;
}
//algoritmo della false position
double Falsepos(double (*fz)(double y),double a,double b,double precis,int &counter)
{
  double fa=fz(a);
  if(fa==0.)return a;
  double fb=fz(b);
  if(fb==0.)return b;
  double c=(a*fb-b*fa)/(fb-fa);
  double fc=fz(c);
do
  {
    c=(a*fb-b*fa)/(fb-fa);
    fc=fz(c);
    if(fa*fc<0)
      {
	b=c;
	fb=fc;
      }
    if(fc*fb<0) 
      {
	a=c;
	fa=fc;
      }
    counter++;
  }
 while(fabs(b-a)>1.e-9&&fabs(fc)>precis);

 return c;
}
//algoritmo di newton
double Newton(double (*fz)(double y),double (*dfz)(double),double a,double b,double precis,int &counter)
{
  double fa=fz(a);
  if(fa==0.)return a;
  double fb=fz(b);
  if(fb==0.)return b;
  double x0=0.5*(a+b);
  double fx0=fz(x0);
  double dfx0=dfz(x0);
  double xk=x0;
  for(int i=0;i<15;i++)
  {
    x0=xk;
    fx0=fz(x0);
    dfx0=dfz(x0);
    xk=x0-(fx0/dfx0);
    counter++;
    if((fabs(b-a)<1.e-9||fabs(fx0)<precis))return xk;
  }
}

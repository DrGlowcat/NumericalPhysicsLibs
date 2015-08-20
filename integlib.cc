//libreria degli algoritmi per l'integrazione
#include "prototypes.h"

double PoorEuler(double (*fz)(double x),double inf,double sup,int &N,int &count)
{
  double ak,bk,Ik=0,Iprev,fzak;
  double dx=(sup-inf)/N;

do
  {
    ak=inf;
    bk=ak+dx;
    Iprev=Ik;
    Ik=0.;
    for(int i=0;i<=N;i++)
      {
	fzak=fz(ak);
	Ik+=(bk-ak)*fzak;
	ak=bk;
	bk+=dx;
      }
  } while(fabs(Ik-Iprev)>1.e-9);
  return Ik;
}

double MidPoint(double (*fz)(double x),double inf,double sup,int &N,int &count)
{
  double ak,bk,Ik=0,Iprev,ck,fzck;
  double dx=(sup-inf)/N;

  do
    {
      ak=inf;
      bk=ak+dx;
      Iprev=Ik;
      Ik=0.;
      for(int i=0;i<=N;i++)
	{
	  ck=0.5*(ak+bk);
	  fzck=fz(ck);
	  Ik+=(bk-ak)*fzck;
	  ak=bk;
	  bk+=dx;
	}
    } while(fabs(Ik-Iprev)>1.e-9);
  return Ik;
}

double Trapezoidal(double (*fz)(double x),double inf,double sup,int &N,int &count)
{
  double ak,bk,fzak,fzbk,Ik=0,Iprev;
  double dx=(sup-inf)/N;

 do
   {
     ak=inf;
     bk=ak+dx;
     fzak=fz(ak);
     fzbk=fz(bk);
     Iprev=Ik;
     Ik=0.;
     for(int i=0;i<=N;i++)
       {
	 Ik+=0.5*(bk-ak)*(fzak+fzbk);
	 ak=bk;
	 bk+=dx;
	 fzak=fzbk;
	 fzbk=fz(bk);
       }//end for
  N*=2;
  count++;
   } while(fabs(Ik-Iprev)>1.e-9);
  return Ik;
}

double Simpson(double (*fz)(double x),double inf,double sup,int &N,int &count)
{
  double ak,bk,fzak,fzbk,ck,fzck,Ik=0,Iprev;
  do{
  double dx=(sup-inf)/N;

  ak=inf;
  bk=ak+dx;
  Iprev=Ik;
  Ik=0.;
  /* for(int i=1;i<=N;i++)
    {
      fzak=fz(ak);
      fzbk=fz(bk);
      ck=0.5*(ak+bk);
      fzck=fz(ck);
      Ik+=((bk-ak)/6)*(fzak+4*fzck+fzbk);
      ak=bk;
      bk+=dx;
      }*/
  //metodo più efficiente:
  int coeff=4;
  Ik=fz(inf)+fz(sup);
  for(int i=1;i<=N-1;i++)
    {
      Ik+=fz(inf+i*dx)*coeff;
      coeff=6-coeff;
    }
  Ik*=(bk-ak)/3.;
  N*=2;
  count++;
  }while(fabs(Ik-Iprev)>1.e-9);
  return Ik;
}

double GaussLegendre(double (*fz)(double x),double inf,double sup,int &N,int &count)
{
  static const int n=8;//grado del polinomio
  double xi[n]={-0.1834346424956498,0.1834346424956498,-0.5255324099163290,0.5255324099163290,-0.7966664774136267,0.7966664774136267,-0.9602898564975363,0.9602898564975363};
  double wi[n]={0.3626837833783620,0.3626837833783620,0.3137066458778873,0.3137066458778873,0.2223810344533745,0.2223810344533745,0.1012285362903763,0.1012285362903763};
  double I=0.,Ik,Iprev,ak,bk,ck;
  do{
    double delta=(sup-inf)/(double)N;
    ak=inf;
    bk=inf+delta;
    Iprev=I;
    I=0.;
      for(int ii=1;ii<=N;ii++)
	{
	  Ik=0.;
	  for(int i=0;i<8;i++)
	    {
	      ck=0.5*(bk-ak)*xi[i]+0.5*(ak+bk);
	      Ik+=wi[i]*fz(ck);
	    }//end i for
	  Ik*=0.5*(bk-ak);
	  I+=Ik;
	  ak=bk;
	  bk+=delta;
	}//end ii for
      N*=2;
      count++;
  }while(fabs(I-Iprev)>1.e-9);
  return I;
}

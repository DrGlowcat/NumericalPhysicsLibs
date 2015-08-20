#include <stdio.h>
#include <omp.h>
#include "prototypes.h"
using namespace std;

void dYdt(double,double*,double*);
void AnaliticSolution(double &,double);

int G_n; //indice politropico
#pragma omp_set_dynamics(0) //disattiva creazione dinamica di thread
#pragma omp threadprivate(G_n) //ogni thread possiede una copia di G_n

int main(){
	
#pragma omp parallel num_threads(8) //creazione di 8 thread
{
  G_n=omp_get_thread_num(); //ogni thread assegna a G_n il valore del suo thread number

  const int Neq=2;//numero incognite
  double start=0.,end=10.,step=1000.; //intervallo su cui varia csi
  double csi=0.,deltacsi=(end-start)/step;
  double Y[Neq],R[Neq]; //funzioni numeriche
  double theta=0.; //funzioni analitiche
  double err; //discrepanza soluzioni analitica e numerica
  char filename[32];
  ofstream f;
  sprintf(filename,"n=%d.dat",G_n);
  f.precision(5);
  f.open(filename,ios::out);
  
    //condizioni iniziali:
  Y[0]=1.; //theta
  Y[1]=0.; //phi
  AnaliticSolution(theta,csi);
  
  err=fabs(Y[0]-theta);
  
  f<<csi<<"   "<<Y[0]<<"   "<<Y[1]<<"   "<<theta<<"   "<<err<<endl;

  do
    {
      theta=0.;
      err=0.;
      RungeKutta4(csi,Y,dYdt,deltacsi,Neq);
      csi+=deltacsi;
      AnaliticSolution(theta,csi);
      err=fabs(Y[0]-theta);
      f<<csi<<"   "<<Y[0]<<"   "<<Y[1]<<"   "<<theta<<"   "<<err<<endl;
    }while(csi<end);
  f.close();
}//fine parte parallela

return 0;
}

void dYdt(double csi,double *Y,double *R)  //funzione derivata
{
  if(fabs(csi)<1.e-16)
    {
      R[0]=0.;
      R[1]=0.;
    }
  else
   {
      R[1]= pow(Y[0],G_n)*csi*csi;
      R[0]=-Y[1]/(csi*csi); //d(theta)/d(csi)
   }
   //R[0]=-csi/3.;
}

void AnaliticSolution(double &sol,double x)  //soluzioni esatte quando possibile
{
  switch(G_n)
    {
  case 0:
    sol=1.-1./6*x*x;
    break;
  case 1:
      if(x==0)sol=1.;
      else    sol=sin(x)/x;
    break;
  case 5:
    sol=1./sqrt(1.+1./3*x*x);
    break;
  default:
  sol=0;
    break;
  }//fine switch
}

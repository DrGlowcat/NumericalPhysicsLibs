//prototipi delle funzioni del corso di Algoritmi Numerici
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include <limits>

//funzioni per la ricerca di zeri:
int Bracket(double (*fz)(double),double &,double &); //trova l'intervallo che contiene lo zero
double Bisect(double (*fz)(double),double,double,double,int&);
double Falsepos(double (*fz)(double),double,double,double,int&);
double Newton(double (*fz)(double),double (*dfz)(double),double,double,double,int&);

//funzioni per l'integrazione:
double PoorEuler(double (*fz)(double),double,double,int&,int&);//metodo di Eulero
double MidPoint(double (*fz)(double),double,double,int&,int&);//metodo mid point
double Trapezoidal(double (*fz)(double),double,double,int&,int&);
double Simpson(double (*fz)(double),double,double,int&,int&);
double GaussLegendre(double (*fz)(double),double,double,int&,int&);

//funzioni per le ODE
void PoorEulerODE(double /*variabile*/,double *  /*array soluzione*/,void(*derivefz)(double,double *,double *),double/*delta*/,int/*#step*/);
void RungeKutta4(double,double *,void (*dYdt)(double,double *,double *),double,int);//Runge-Kutta 4 ordine
//se c'è tempo prova a implementare leapfrog e verlet

//funzioni per Boudary Value Problems


//elenco di costanti
#define CONST_me 0.510999; //massa elettrone MeV
#define CONST_mp 938.272; //massa protone in MeV
#define CONST_ht 4.1357e-15; //h tagliato in eV*s
#define CONST_c 3.e8; //velocità luce in m/s
#define CONST_g 9.81 //acc gravità in m/s^2
#define CONST_pi 3.1415 //pigreco
#define CONST_G 6.67e-8 //G in cgs
#define CONST_R0 6.96e10 //raggio solare in cm

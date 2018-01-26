#ifndef MC_H
#define MC_H


#include <stdlib.h>
#include <math.h>

double square(double x);

double cube(double x);

double dot(double a[], double b[]);

/* Generate a random double in (0, 1]*/
double ranf();

/* Generate normally distributed double (mean 0, variance 1)*/
double rang();

double nni_delta(double i, double j, double L);

//double pair_energy(double i[], double j[], double L);

double r_six(double i[], double j[], double L);

int mcmove(double old[], double& energy, double& vir, double& rhoN, double k[], int N, double DELTA, double BETA, double L, int& count_accept, int& count_times);


#endif // MC_H

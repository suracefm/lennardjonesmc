#ifndef MC_H
#define MC_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sstream>

double square(double x);

double cube(double x);

double dot(double a[], double b[]);

/* Generate a random double in (0, 1]*/
double ranf();

/* Generate normally distributed double (mean 0, variance 1)*/
double rang();

double nni_delta(double i, double j, double L);

double r_six(double i[], double j[], double L);

int thermalizationmove(double old[], double& energy, double& vir, double& rhok, double k[], int N, double DELTA, double BETA, double L, int& count_accept, int& count_times);

int mcmove(double old[], double& energy, double& vir, int N, double DELTA, double BETA, double L, int& count_accept, int& count_times);

void strtoin(std::string line, int &N1, double &L, double &DELTA, double &BETA, int &MCTIME, int &BTIME, int &BINMIN);

#endif // MC_H

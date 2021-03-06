#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "mc.h"
#include <fstream>
#include <sstream>

#define PI 3.14159265

double square(double x){
    return x*x;
}

double cube(double x){
    return x*x*x;
}

double dot(double a[], double b[]){
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

/* Generate a random double in (0, 1]*/
double ranf(){
    int i = rand();
    if (i!=0) return (double)i/(double)(RAND_MAX);
    else return ranf();
}

/* Generate normally distributed double (mean 0, variance 1)*/
double rang(){
    double u1 = ranf();
    double u2 = ranf();
    return sqrt(-2*log(u1))*cos(2*PI*u2);
}

double nni_delta(double i, double j, double L){
    return i-j-L*rint((i-j)/L);
}


double r_six(double i[], double j[], double L){
    return 1/cube(square(nni_delta(i[0],j[0], L))+square(nni_delta(i[1],j[1], L))+square(nni_delta(i[2], j[2], L)));
}

int mcmove(double old[], double& energy, double& vir, int N, double DELTA, double BETA, double L, int& count_accept, int& count_times){
    int i = rand()%N; //not safe!
    double new_pos[3];
    new_pos[0] = old[3*i]+DELTA*rang();
    new_pos[1] = old[3*i+1]+DELTA*rang();
    new_pos[2] = old[3*i+2]+DELTA*rang();
    double delta_energy = 0, delta_vir = 0, r_six_new, r_six_old, cutoff=cube(square(2/L));
    for (int j=0; j<N; j++){
        if (i!=j){
	    r_six_new=r_six(new_pos, old+3*j, L);
	    r_six_old=r_six(old+3*i, old+3*j, L);
	    if (r_six_old>cutoff){
		delta_energy+=4*(-square(r_six_old)+r_six_old);
		delta_vir+=8*(-2*square(r_six_old)+r_six_old);
	    }
	    if (r_six_new>cutoff){
 		delta_energy+=4*(square(r_six_new)-r_six_new);
		delta_vir+=8*(2*square(r_six_new)-r_six_new);
	    }
	}
    }
    if (delta_energy>0){
        count_times++;
        if (ranf()>exp(-BETA*delta_energy)) return 0;
        else count_accept++;
    }
    old[3*i] = new_pos[0];
    old[3*i+1] = new_pos[1];
    old[3*i+2] = new_pos[2];
    energy+=delta_energy/N;
    vir+=delta_vir/N;   
    return 0;
}

int thermalizationmove(double old[], double& energy, double& vir, double& rhok, double k[], int N, double DELTA, double BETA, double L, int& count_accept, int& count_times){
    int i = rand()%N; //not safe!
    double new_pos[3];
    new_pos[0] = old[3*i]+DELTA*rang();
    new_pos[1] = old[3*i+1]+DELTA*rang();
    new_pos[2] = old[3*i+2]+DELTA*rang();
    double delta_energy = 0;
    double delta_vir = 0;
    double r_six_new, r_six_old;
    double cutoff=cube(square(double(2)/L));
    for (int j=0; j<N; j++){
        if (i!=j){
	    r_six_new=r_six(new_pos, old+3*j, L);
	    r_six_old=r_six(old+3*i, old+3*j, L);
	    if (r_six_old>cutoff){
		delta_energy+=4*(-square(r_six_old)+r_six_old);
		delta_vir+=8*(-2*square(r_six_old)+r_six_old);
	    }
	    if (r_six_new>cutoff){
 		delta_energy+=4*(square(r_six_new)-r_six_new);
		delta_vir+=8*(2*square(r_six_new)-r_six_new);
	    }
	}
    }
    if (delta_energy>0){
        count_times++;
        if (ranf()>exp(-BETA*delta_energy)) return 0;
        else count_accept++;
    }
    rhok+=(cos(dot(k, new_pos))-cos(dot(k, old+3*i)))/N;
    old[3*i] = new_pos[0];
    old[3*i+1] = new_pos[1];
    old[3*i+2] = new_pos[2];
    energy+=delta_energy/N;
    vir+=delta_vir/N;   
    return 0;
}


void strtoin(std::string line, int &N1, double &L, double &DELTA, double &BETA, int &MCTIME, int &BTIME, int &THBIN)
{
    if (line.length()>8){
     	std::istringstream iss(line.substr(8));
        if (line.substr(0,8)=="N1      ") {iss >> N1;}
        else if (line.substr(0,8)=="L       ") {iss >> L;}
        else if (line.substr(0,8)=="DELTAIN ") {iss >> DELTA;}
        else if (line.substr(0,8)=="BETA    ") {iss >> BETA;}
	else if (line.substr(0,8)=="MCTIME  ") {iss >> MCTIME;}
	else if (line.substr(0,8)=="BTIME   ") {iss >> BTIME;}
	else if (line.substr(0,8)=="THBIN   ") {iss >> THBIN;}
    }
}

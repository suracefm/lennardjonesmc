#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mc.h"

#define PI 3.14159265

using namespace std;


void strtoin(string line, int &N1, double &L, double &DELTA, double &BETA, int &NSTEPS, int &EQTIME){
        if (line.length()>8){
               	istringstream iss(line.substr(8));
                if (line.substr(0,8)=="N1      ") {iss >> N1;}
                else if (line.substr(0,8)=="L       ") {iss >> L;}
                else if (line.substr(0,8)=="DELTA   ") {iss >> DELTA;}
                else if (line.substr(0,8)=="BETA    ") {iss >> BETA;}
		else if (line.substr(0,8)=="NSTEPS  ") {iss >> NSTEPS;}
		else if (line.substr(0,8)=="EQTIME  ") {iss >> EQTIME;}
        }
}


int main(int argc, char* argv[]){
 if(argc != 3) {
        std::cerr << "Usage: " << argv[0]
                << " file.in file.out " << std::endl;
        return 1;
    }

int N1, NSTEPS, EQTIME;
double L, DELTA, BETA;

  ifstream input_file (argv[1]);
  if (input_file.is_open())
  {
   	string line;
        while (std::getline(input_file, line)){strtoin(line, N1, L, DELTA, BETA, NSTEPS, EQTIME);}
    input_file.close();
  }
  else cout<<"Unable to open file.\n";

    srand(time(NULL));
    int N = N1*N1*N1;
    double n1=double(N1)/L;
    double pos[3*N];
    double k[]={2*PI*n1, 0,0};

    // Initial positions and initial energy
    for (int i=0; i<N; i++){
        pos[3*i] = (i/(N1*N1))/n1;
        pos[3*i+1] = ((i%(N1*N1))/N1)/n1;
        pos[3*i+2] = (i%N1)/n1;
    }
    double energy_init = 0;
    double vir_init = 0;
    double rhoN_init = 0;
    double rmsix;
    for (int i=0; i<N; i++){
        for (int j=i+1; j<N; j++){
		rmsix=r_six(pos+3*i, pos+3*j, L);
		energy_init+=4*(square(rmsix)-rmsix);
		vir_init+=24*(2*square(rmsix)-rmsix);
        }
	rhoN_init+=cos(dot(k,pos+3*i));
    }
    energy_init/=N;
    vir_init/=N;
    rhoN_init/=N;

    int count_accept, count_times, counter;
    double entot, ensqtot, energybin, energy, vir, rhoN;


    ofstream outfile;
    outfile.open (argv[2]);
    outfile<<"# N1 "<<N1<<endl;
    outfile<<"# L "<<L<<endl;
    outfile<<"# DELTA "<<DELTA<<endl;
    outfile<<"# BETA "<<BETA<<endl;
    outfile<<"# NSTEPS "<<NSTEPS<<endl;
    outfile<<"# EQTIME "<<EQTIME<<endl;
    outfile<<"#\n#BINSTEP\tvar\tk\tacceptance"<<endl;

    int binsteps[] = {10,20,50,100,200,500,1000, 1250, 1600, 2000, 2500, 3125, 4000, 5000, 6250, 8000, 10000};

    for(const int &BINSTEP : binsteps){
	count_accept=0;
    	count_times=0;
    	energy = energy_init;
	vir = vir_init;
	rhoN =rhoN_init;
    	entot=0;
    	ensqtot=0;
    	energybin =0;
    	counter=0;
	    for(int i=-EQTIME; i<NSTEPS; i++){
		mcmove(pos, energy, vir, rhoN, k, N, DELTA, BETA, L, count_accept, count_times);
		if (i>=0){
			energybin+=energy;
			if (i%BINSTEP == BINSTEP-1){
			    energybin/=BINSTEP;
			    entot+=energybin;
			    ensqtot+=square(energybin);
			    counter++;
			    energybin=0;
			}
		}
	    }
	double var = (ensqtot/counter-square(entot/counter))/(counter-1);
    	outfile<<BINSTEP<<"\t"<<sqrt(var)<<"\t"<<counter<<"\t"<<double(count_accept)/double(count_times)<<endl;
	for (int i=0; i<N; i++){
        	pos[3*i] = (i/(N1*N1))/n1;
        	pos[3*i+1] = ((i%(N1*N1))/N1)/n1;
        	pos[3*i+2] = (i%N1)/n1;
    	}
    }
    
    outfile.close();

    return 0;
}


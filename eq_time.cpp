/*
Compute energy, rho as a function of time in order to get equilibration time
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mc.h"

#define PI 3.14159265


using namespace std;

// get input from string
void strtoin(string line, int &N1, double &L, double &DELTA, double &BETA, int &NSTEPS, int &BINSTEP){
        if (line.length()>8){
               	istringstream iss(line.substr(8));
                if (line.substr(0,8)=="N1      ") {iss >> N1;}
                else if (line.substr(0,8)=="L       ") {iss >> L;}
                else if (line.substr(0,8)=="DELTA   ") {iss >> DELTA;}
                else if (line.substr(0,8)=="BETA    ") {iss >> BETA;}
		else if (line.substr(0,8)=="NSTEPS  ") {iss >> NSTEPS;}
		else if (line.substr(0,8)=="BINSTEP ") {iss >> BINSTEP;}
        }
}


int main(int argc, char* argv[]){
// check syntax
 if(argc != 3) {
        std::cerr << "Usage: " << argv[0]
                << " file.in file.out " << std::endl;
        return 1;
    }

int N1, NSTEPS, BINSTEP;
double L, DELTA, BETA;

// read input
  ifstream input_file (argv[1]);
  if (input_file.is_open())
  {
   	string line;
        while (std::getline(input_file, line)){strtoin(line, N1, L, DELTA, BETA, NSTEPS, BINSTEP);}
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
    double energy = 0;
    double vir = 0;
    double rhoN = 0;
    double rmsix;
    for (int i=0; i<N; i++){
        for (int j=i+1; j<N; j++){
		rmsix=r_six(pos+3*i, pos+3*j, L);
		energy+=4*(square(rmsix)-rmsix);
		vir+=24*(2*square(rmsix)-rmsix);
        }
	rhoN+=cos(dot(k,pos+3*i));
    }
    energy/=N;
    vir/=N;
    rhoN/=N;
    //std::cout<<energy/N<<std::endl;

    int count_accept=0;
    int count_times=0;
    double entot=0;
    double ensqtot=0;
    //double energybin =0;
    //double ensqbin =0;
    int counter=0;

// Print input in output file
    ofstream energyfile;
    energyfile.open (argv[2]);
    energyfile<<"# N1 "<<N1<<endl;
    energyfile<<"# L "<<L<<endl;
    energyfile<<"# DELTA "<<DELTA<<endl;
    energyfile<<"# BETA "<<BETA<<endl;
    energyfile<<"# NSTEPS "<<NSTEPS<<endl;
    energyfile<<"# BINSTEP "<<BINSTEP<<endl;
    energyfile<<"#\n# ENERGY\tVIRIAL\tRHO"<<endl;

// MC cycle
    for(int i=0; i<NSTEPS; i++){
        mcmove(pos, energy, vir, rhoN, k, N, DELTA, BETA, L, count_accept, count_times);
        //energybin+=energy;
        if (i%BINSTEP == BINSTEP-1){
            entot+=energy;
            ensqtot+=square(energy);
            counter++;
            energyfile << energy<<"\t"<<vir<<"\t"<<rhoN<<endl;
            //energybin=0;
            //cout<<energy/N<<endl;
        }
    }
    
    energyfile<<"# acceptance "<<double(count_accept)/double(count_times)<<endl;
    double var = (ensqtot/counter-square(entot/counter))/(counter-1);
    energyfile<<"# avg energy per particle "<<entot/counter/N<<endl;
    energyfile<<"# error "<<sqrt(var)/N<<endl;
    energyfile.close();

    return 0;
}


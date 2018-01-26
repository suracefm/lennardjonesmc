#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <map>
#include "mc.h"

#define PI 3.14159265

using namespace std;


void strtoin(string line, int &N1, double &L, int &BINSTEP, int &NSTEPS, int &EQTIME, map<double,double> &bdtable){
	double BETA, DELTA;
        if (line.length()>8){
               	istringstream iss(line.substr(8));
                if (line.substr(0,8)=="N1      ") {iss >> N1;}
                else if (line.substr(0,8)=="L       ") {iss >> L;}
                else if (line.substr(0,8)=="BINSTEP ") {iss >> BINSTEP;}
		else if (line.substr(0,8)=="NSTEPS  ") {iss >> NSTEPS;}
		else if (line.substr(0,8)=="EQTIME  ") {iss >> EQTIME;}
		else if (line.substr(0,8)=="BD PAIR ") {
			iss>>BETA>>DELTA;
			bdtable[BETA]=DELTA;
		}
        }
}


int main(int argc, char* argv[]){
 if(argc != 3) {
        std::cerr << "Usage: " << argv[0]
                << " file.in file.out " << std::endl;
        return 1;
    }

int N1, NSTEPS, EQTIME, BINSTEP;
double L;
map<double,double> bdtable;

  ifstream input_file (argv[1]);
  if (input_file.is_open())
  {
   	string line;
        while (std::getline(input_file, line)){strtoin(line, N1, L, BINSTEP, NSTEPS, EQTIME, bdtable);}
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
    double rhok_init = 0;
    double rmsix;
    for (int i=0; i<N; i++){
        for (int j=i+1; j<N; j++){
		rmsix=r_six(pos+3*i, pos+3*j, L);
		energy_init+=4*(square(rmsix)-rmsix);
		vir_init+=24*(2*square(rmsix)-rmsix);
        }
	rhok_init+=cos(dot(k,pos+3*i));
    }
    energy_init/=N;
    vir_init/=N;
    rhok_init/=N;

    int count_accept, count_times, counter;
    double energy1, energy2, energy3, energy4, energy_step, vir1, vir2, vir_step, rhok, rhok_step, BETA, DELTA;


    ofstream outfile;
    outfile.open (argv[2]);
    outfile<<"# N1 "<<N1<<endl;
    outfile<<"# L "<<L<<endl;
    outfile<<"# BINSTEP "<<BINSTEP<<endl;
    outfile<<"# NSTEPS "<<NSTEPS<<endl;
    outfile<<"# EQTIME "<<EQTIME<<endl;
    outfile<<"#\n#BETA\tDELTA\tENERGY1\t\tENERGY2\t\tENERGY3\t\tENERGY4\t\tVIR\t\tRHON\t\tacceptance"<<endl;


    for(const pair<double,double> bd : bdtable){
	count_accept=0;
    	count_times=0;
    	energy_step = energy_init;
	vir_step = vir_init;
	rhok_step = rhok_init;
	BETA=bd.first;
	DELTA=bd.second;

	energy1=0;
	energy2=0;
	energy3=0;
	energy4=0;
	vir1=0;
	vir2=0;
	rhok=0;
    	counter=0;
	    for(int i=-EQTIME; i<NSTEPS; i++){
		mcmove(pos, energy_step, vir_step, rhok_step, k, N, DELTA, BETA, L, count_accept, count_times);
		if ((i>=0)&&(i%BINSTEP == BINSTEP-1)){
		    energy1+=energy_step;
		    energy2+=square(energy_step);
		    energy3+=square(energy_step)*energy_step;
		    energy4+=square(square(energy_step));
		    vir1+=vir_step;
		    vir2+=square(vir_step);
		    rhok+=rhok_step;
		    counter++;
		}
	    }
    	outfile<<BETA<<"\t"<<DELTA<<"\t"<<energy1/counter<<"\t\t"<<energy2/counter<<"\t\t"<<energy3/counter<<"\t\t"<<energy4/counter<<"\t\t"<<vir1/counter<<"\t\t"<<vir2/counter<<"\t\t"<<rhok/counter<<"\t\t"<<double(count_accept)/double(count_times)<<endl;
    for (int i=0; i<N; i++){
        pos[3*i] = (i/(N1*N1))/n1;
        pos[3*i+1] = ((i%(N1*N1))/N1)/n1;
        pos[3*i+2] = (i%N1)/n1;
    }
    }
    
    outfile.close();

    return 0;
}


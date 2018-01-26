#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mc.h"

#define PI 3.14159265

using namespace std;


void strtoin(string line, int &N1, double &L, double &DELTA, double &BETA, int &MCTIME, int &BTIME, int &BINMIN){
        if (line.length()>8){
               	istringstream iss(line.substr(8));
                if (line.substr(0,8)=="N1      ") {iss >> N1;}
                else if (line.substr(0,8)=="L       ") {iss >> L;}
                else if (line.substr(0,8)=="DELTA   ") {iss >> DELTA;}
                else if (line.substr(0,8)=="BETA    ") {iss >> BETA;}
		else if (line.substr(0,8)=="MCTIME  ") {iss >> MCTIME;}
		else if (line.substr(0,8)=="BTIME   ") {iss >> BTIME;}
		else if (line.substr(0,8)=="BINMIN  ") {iss >> BINMIN;}
        }
}


int main(int argc, char* argv[]){
 if(argc != 3) {
        std::cerr << "Usage: " << argv[0]
                << " file.in file.out " << std::endl;
        return 1;
    }

  int N1, MCTIME, BTIME, BINMIN;
  double L, DELTA, BETA;

  ifstream input_file (argv[1]);
  if (input_file.is_open())
  {
   	string line;
        while (std::getline(input_file, line)){strtoin(line, N1, L, DELTA, BETA, MCTIME, BTIME, BINMIN);}
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
    double rhok = 0;
    double rmsix;
    for (int i=0; i<N; i++){
        for (int j=i+1; j<N; j++){
		rmsix=r_six(pos+3*i, pos+3*j, L);
		energy+=4*(square(rmsix)-rmsix);
		vir+=24*(2*square(rmsix)-rmsix);
        }
	rhok+=cos(dot(k,pos+3*i));
    }
    energy/=N;
    vir/=N;
    rhok/=N;

    int count_accept, count_times, counter;
    double energy1, energy2, energy3, energy4, vir1, vir2, rhok, L, energybin, virbin, rhokbin;


    ofstream outfile;
    outfile.open (argv[2]);
    outfile<<"# N1 "<<N1<<endl;
    outfile<<"# L "<<L<<endl;
    outfile<<"# DELTA "<<DELTA<<endl;
    outfile<<"# BETA "<<BETA<<endl;
    outfile<<"# MCTIME "<<MCTIME<<endl;
    outfile<<"# BTIME "<<BTIME<<endl;
    outfile<<"# BINMIN "<<BTIME<<endl;
    //outfile<<"#\n#BINSTEP\tvar\tk\tacceptance"<<endl;

    count_accept=0;
    count_times=0;
    entot=0;
    ensqtot=0;
    energybin=0;
    virbin=0;
    rhokbin=0;    
    counter=0;



    // Thermalization

    outfile<<"#\n#\n# *****THERMALIZATION*****\n#\n# ENERGY\tVIRIAL\tRHO"<<endl;
    int thermtime=0;    
    while(rhokbin>0){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int i=0; i<BTIME; i++){
	    mcmove(pos, energy, vir, rhok, k, N, DELTA, BETA, L, count_accept, count_times);
	    thermtime++;
	    energybin+=energy;
	    virbin+=vir;
	    rhokbin+=rhok;
	 }
	outfile << energybin/BTIME<<"\t"<<virbin/BTIME<<"\t"<<rhokbin/BTIME<<endl;
	}
	thermtime*=2;
    for(int i=0;i<thermtime; i++){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int i=0; i<BTIME; i++){
	    mcmove(pos, energy, vir, rhok, k, N, DELTA, BETA, L, count_accept, count_times);
	    thermtime++;
	    energybin+=energy;
	    virbin+=vir;
	    rhokbin+=rhok;
	 }
	outfile << energybin/BTIME<<"\t"<<virbin/BTIME<<"\t"<<rhokbin/BTIME<<endl;
	}
    outfile <<"#\n#THERMALIZATION TIME "<<thermtime<<endl;




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

   outfile<<"# ACCEPTANCE"<<double(count_accept)/double(count_times)<<endl;
    
    outfile.close();

    return 0;
}


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mc.h"

#define PI 3.14159265

using namespace std;

int main(int argc, char* argv[])
{
    // Read input
    if(argc != 3) {
        std::cerr << "Usage: " << argv[0]<< " file.in file.out " << std::endl;
        return 1;
    }
    int N1, MCTIME, BTIME, BINMIN;
    double L, DELTA, BETA;
    ifstream input_file (argv[1]);
    if (input_file.is_open()){
   	string line;
        while (std::getline(input_file, line)){
	    strtoin(line, N1, L, DELTA, BETA, MCTIME, BTIME, BINMIN);}
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
    double energy = 0, vir=0, rhok=0, rmsix;
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

    // Writing input data to output
    ofstream outfile;
    outfile.open (argv[2]);
    outfile<<"# N1 "<<N1<<endl;
    outfile<<"# L "<<L<<endl;
    outfile<<"# DELTAIN "<<DELTA<<endl;
    outfile<<"# BETA "<<BETA<<endl;
    outfile<<"# MCTIME "<<MCTIME<<endl;
    outfile<<"# BTIME "<<BTIME<<endl;
    outfile<<"# BINMIN "<<BINMIN<<endl;


    // THERMALIZATION
    outfile<<"#\n#\n# *****THERMALIZATION*****\n#\n# ENERGY\tVIRIAL\tRHO"<<endl;
    int count_accept=0, count_times=0, counter=0, thermtime =0;
    double energybin, virbin, rhokbin=1;   
    while(rhokbin>0){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int i=0; i<BINMIN; i++){
	    thermalizationmove(pos, energy, vir, rhok, k, N, DELTA, BETA, L, count_accept, count_times);
	    thermtime++;
	    energybin+=energy;
	    virbin+=vir;
	    rhokbin+=rhok;
	}
	outfile << energybin/BINMIN<<"\t"<<virbin/BINMIN<<"\t"<<rhokbin/BINMIN<<endl;
	count_accept=0;
	count_times=0;
    }
    for(int i=0;i<thermtime; i+=BINMIN){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int i=0; i<BINMIN; i++){
	    thermalizationmove(pos, energy, vir, rhok, k, N, DELTA, BETA, L, count_accept, count_times);
	    energybin+=energy;
	    virbin+=vir;
	    rhokbin+=rhok;
	 }
	outfile << energybin/BINMIN<<"\t"<<virbin/BINMIN<<"\t"<<rhokbin/BINMIN<<endl;
    }
    outfile <<"#\n# THERMTIME "<<thermtime<<endl;
    outfile <<"# DELTA "<<DELTA<<endl;


    // BINNING TECHNIQUE: error bar estimate
    outfile<<"#\n#\n# ***** ERROR BAR ESTIMATE *****\n#\n# ENERGY\tVIRIAL"<<endl;
    double entot=0, ensqtot=0, virtot=0;
    for(int i=0; i<BTIME; i+=BINMIN){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int j=0; j<BINMIN; j++){
	    mcmove(pos, energy, vir, N, DELTA, BETA, L, count_accept, count_times);
	    energybin+=energy;
	    virbin+=vir;
	    entot+=energy;
	    ensqtot+=square(energy);
	    virtot+=vir;
	    counter++;
	}
	outfile << energybin/BINMIN<<"\t"<<virbin/BINMIN<<"\t"<<endl;
    }


    // MC SIMULATION
    for(int i=BTIME; i<MCTIME; i++){
	mcmove(pos, energy, vir, N, DELTA, BETA, L, count_accept, count_times);
	entot+=energy;
	ensqtot+=square(energy);
	virtot+=vir;
	counter++;
    }

    outfile<<"# ACCEPTANCE "<<double(count_accept)/double(count_times)<<endl;
    outfile<<"# ENERGY "<<entot/counter<<endl;
    outfile<<"# ENERGYSQUARE "<<ensqtot/counter<<endl;
    outfile<<"# VIRIAL1 "<<vir/counter<<endl;
    outfile<<"# MCTIME "<<counter<<endl;
    
    outfile.close();

    return 0;
}


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
    int N1, MCTIME, BTIME, THBIN;
    double L, DELTA, BETA;
    ifstream input_file (argv[1]);
    if (input_file.is_open()){
   	string line;
        while (std::getline(input_file, line)){
	    strtoin(line, N1, L, DELTA, BETA, MCTIME, BTIME, THBIN);}
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
		vir+=8*(2*square(rmsix)-rmsix);
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
    outfile<<"# THBIN "<<THBIN<<endl;


    // THERMALIZATION
    outfile<<"#\n#\n# *****THERMALIZATION*****\n#\n# ENERGY\tVIRIAL\tRHO"<<endl;
    int count_accept=0, count_times=0, counter=0, thermtime =0;
    double energybin, virbin, rhokbin=1, ensqbin;   
    while(rhokbin>0){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int i=0; i<THBIN; i++){
	    thermalizationmove(pos, energy, vir, rhok, k, N, DELTA, BETA, L, count_accept, count_times);
	    thermtime++;
	    energybin+=energy;
	    virbin+=vir;
	    rhokbin+=rhok;
	}
	outfile << energybin/THBIN<<"\t"<<virbin/THBIN<<"\t"<<rhokbin/THBIN<<endl;
	count_accept=0;
	count_times=0;
    }
    for(int i=0;i<thermtime; i+=THBIN){
	energybin=0;
	virbin=0;
	rhokbin=0;
	for(int i=0; i<THBIN; i++){
	    thermalizationmove(pos, energy, vir, rhok, k, N, DELTA, BETA, L, count_accept, count_times);
	    energybin+=energy;
	    virbin+=vir;
	    rhokbin+=rhok;
	 }
	outfile << energybin/THBIN<<"\t"<<virbin/THBIN<<"\t"<<rhokbin/THBIN<<endl;
    }
    outfile <<"#\n# THERMTIME "<<thermtime<<endl;
    outfile <<"# DELTA "<<DELTA<<endl;


    // MC SIMLULATION
    time_t begin, end;
    time(&begin);
    outfile<<"#\n#\n# ***** ERROR BAR ESTIMATE *****\n#\n# ENERGY\t ENERGYSQUARE\tVIRIAL"<<endl;
    for(int i=0; i<MCTIME; i+=BTIME){
	energybin=0;
	ensqbin=0;
	virbin=0;
	for(int j=0; j<BTIME; j++){
	    mcmove(pos, energy, vir, N, DELTA, BETA, L, count_accept, count_times);
	    energybin+=energy;
	    virbin+=vir;
	    ensqbin+=square(energy);
	    counter++;
	}
	outfile << energybin/BTIME<<"\t"<< ensqbin/BTIME<<"\t"<<virbin/BTIME<<endl;
    }
    time(&end);

    outfile<<"# ACCEPTANCE "<<double(count_accept)/double(count_times)<<endl;
    outfile<<"# MCTIME "<<counter<<endl;
    outfile<<"# SIMTIME "<<difftime(end,begin)<<endl; //in seconds, thermalization time not included
    
    outfile.close();

    return 0;
}


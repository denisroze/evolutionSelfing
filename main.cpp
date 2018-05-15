// main() function: reads parameter values in input file,
// broadcast to all processors, initiate sprng random number generator
// and runs the simulation.

#include "fisher.h"
#include "MersenneTwister.h"
#include <string>
#include <iostream>
using namespace std;

// input and output files:

MTRand rnd;
FILE * fichierS;

int main(int argc, char* argv[])
    {
        // definitions of variables:
        
        int Nt, n, m, nbS, NbGen, NbGenPreli, pas, smp, nbneut, nbself;
        double s, sbar, omega, Q, U, L, Uself, sigself, kappa, Ve, uneut, sel, h, Udel;
        
        // opens output files:
        
        ouvrirFichierS();
        
        
        // Input parameters
        
        Nt = atoi(argv[1]);
        s = atof(argv[2]) ;
        n = atoi(argv[3]);
        m = atoi(argv[4]);
        sbar = atof(argv[5]);
        omega = atof(argv[6]);
        Q = atof(argv[7]);
        U = atof(argv[8]);
        nbS = atoi(argv[9]);
        Udel = atof(argv[10]);
        sel = atof(argv[11]);
        h = atof(argv[12]);
        uneut = atof(argv[13]);
        nbneut = atoi(argv[14]);
        Uself = atof(argv[15]) ;
        sigself = atof(argv[16]);
        nbself = atoi(argv[17]);
        Ve = atof(argv[18]);
        kappa = atof(argv[19]) ;
        L = atof(argv[20]) ;
        NbGen = atoi(argv[21]) ;
        NbGenPreli = atoi(argv[22]) ;
        pas = atoi(argv[23]) ;
        smp = atoi(argv[24]) ;
        

    ecrireParametres(Nt, s, n, m, sbar, omega, Q, U, nbS, Udel, sel, h, uneut, nbneut, Uself, sigself, nbself, Ve, kappa, L, NbGen, NbGenPreli, pas, smp);
    recursion(Nt, s, n, m, sbar, omega, Q, U, nbS, Udel, sel, h, uneut, nbneut, Uself, sigself, nbself, Ve, kappa, L, NbGen, NbGenPreli, pas, smp);
        cout<<rec<<endl;
	
	// closes files:

	fclose(fichierS);
	
	return 0;
}

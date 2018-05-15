// Functions to open input and output files,
// read parameter values from input file and 
// write them in output file.

#include "fisher.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierS;


// opens output file:

void ouvrirFichierS()   
{
	fichierS = fopen(fichierEcriture,"a");
}



// writes parameter values in output file:

void ecrireParametres(int Nv, double sv, int nv, int mv, double sbarv, double omegav, double Qv, double Uv, int nbSv, double Udelv, double selv, double hv, double uneutv, int nbneutv, double Uselfv, double sigselfv, int nbseflv, double Vev, double kappav, double Lv, int NbGenv,int NbGenPreliv, int pasv, int smpv)
{
	fprintf(fichierS,"\n_________________________________________\n");
	fprintf(fichierS,"\nN = %d", Nv);
	fprintf(fichierS,", s = %g", sv);
	fprintf(fichierS,", n = %d", nv);
	fprintf(fichierS,", m = %d", mv);
	fprintf(fichierS,", sbar = %g", sbarv);
	fprintf(fichierS,", omega = %g", omegav);
	fprintf(fichierS,", Q = %g", Qv);
	fprintf(fichierS,", U = %g", Uv);
	fprintf(fichierS,", nbS = %d", nbSv);
	fprintf(fichierS,", L = %g", Lv);
	fprintf(fichierS,"\n Uself = %g", Uselfv);
	fprintf(fichierS,", generations = %d", NbGenv);
	fprintf(fichierS,", pas = %d", pasv);
	fprintf(fichierS,", smp = %d", smpv);
}

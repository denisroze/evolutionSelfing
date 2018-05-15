// Header file: definitions of global variables, function prototypes

#ifndef FISHER_H
#define FISHER_H

#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>	/* for dynamic_bitset objects (www.boost.org)*/
#include "MersenneTwister.h"
using namespace std;


// Global variables:

#define fichierLecture "parametres.txt"     // names of input
#define fichierEcriture "resultats.txt"		// and output files

// "chr": represents a genome segment

struct chr
{
    boost::dynamic_bitset<> sel; // selected loci (chain of 0 and 1)
    vector<double> incond; // vector of positions of unconditionally deleterious mutations with fixed effect
    int nbchr;   // number of copies of this genome segment in the population
    float * selfing; // allelic values at loci controlling the selfing rate
    float * neut; //allelic values at neutral loci
};

// "Nall" is used to count allele frequencies at the neutral loci:

struct Nall
{
	double all;
	double freq;
};

// Function prototypes:

void ouvrirFichierE();
void ouvrirFichierS();
void ecrireParametres(int Nv, double sv, int nv, int mv, double sbarv, double omegav, double Qv, double Uv, int nbSv, double Udelv, double selv, double hv, double uneutv, int nbneutv, double Uselfv, double sigselfv, int nbselfv,double Vev, double kappav, double Lv, int NbGenv, int NbGenPreliv, int pasv, int smpv);

void recursion(int Nv, double sv, int nv, int mv, double sbarv, double omegav, double Qv, double Uv, int nbSv, double Udelv, double selv, double hv, double uneutv, int nbneutv, double Uselfv, double sigselfv, int nbselfv, double Vev, double kappav, double Lv,int NbGenv, int NbGenPreliv, int pasv, int smpv);
double gammln(const double xx);
double poisdev(const double xm);
double gasdev();
double binldev(const double pp, const int n);
void rec(chr &res, chr &c1, chr &c2, int nbCo, int nS, int nself, int nneut, float * posU, float * posneut);
void cntl_c_handler(int bidon);
double fitness(chr &c1, chr &c2, double wHe, double wHo);

#endif

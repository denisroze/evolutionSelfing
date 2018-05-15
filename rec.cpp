#include "fisher.h"
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <cmath>
#include <algorithm>
using namespace std;

extern MTRand rnd;

//fitness function calculates the multiplicative fitness due to uncondtionally deleterious mutations

double fitness(chr &c1, chr &c2, double wHe, double wHo)
{
    int s1 = c1.incond.size() - 1;
    int s2 = c2.incond.size() - 1;
    
    double w = 1.0;
    
    while ((s1 > -1) || (s2 > -1))
    {
        if (s1 == -1)
        {
            w *= pow(wHe, s2 + 1);
            break;         }
        
        if (s2 == -1)
        {
            w *= pow(wHe, s1 + 1);
            break;
        }
        
        // heterozygous mutation on c1:
        
        if (c1.sel[s1] > c2.sel[s2])
        {
            w *= wHe;
            s1--;
            continue;
        }
        
        // heterozygous mutation on c2:
        
        if (c1.sel[s1] < c2.sel[s2])
        {
            w *= wHe;
            s2--;
            continue;
        }
        
        // homozygous mutation:
        
        w *= wHo;
        s1--; 
        s2--;
    }
    
    return w;
}

// rec function: generates recombinant genome segment "res"
// from parental genome segments "c1" and "c2"
// nbCo is the number of cross-overs in the genome segment,
// nS the number of selected loci

void rec(chr &res, chr &c1, chr &c2, int nbCo, int nS, int nself, int nneut, float * posU, float * posneut)
{
    
	vector<int> pos;
	int j, locus;
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
    int nbsel1 = c1.incond.size()-1;
    int nbsel2 = c2.incond.size()-1;

	res.sel.clear();
    res.incond.clear();
    
    
	// vector "pos" holds the positions of cross-overs:
	
	for (j = 0; j < nbCo; j++)
		pos.push_back(int(nS * rnd.rand()));
	sort(pos.begin(), pos.end());
    

	// creates recombination mask:

	for (j = 0; j < nbCo; j++)
		rec.resize(pos[j], (j % 2) == 0 ? 1 : 0);
	    rec.resize(nS, (nbCo % 2) == 0 ? 1 : 0);

	off1 = (c1.sel & rec);
	rec.flip();
	off2 = (c2.sel & rec);
	res.sel = (off1 | off2);
	
    // loci with inconditionally deleterious effects:
   int locus1 = 0;
    int locus2 = 0;
    for (j = 0; j < nbCo; j++)
    {

        if (j % 2 == 0){
            while ((locus1 < nbsel1) && (c1.incond[locus1] < pos[j]))
            {
                res.incond.push_back(c1.incond[locus1]);
                locus1++;
            }
            while ((locus2 < nbsel2) && (c2.incond[locus2] < pos[j]))
            {
                locus2++;
            }
            
        }
        else{
            while ((locus2 < nbsel2) && (c2.incond[locus2] < pos[j]))
            {
                res.incond.push_back(c2.incond[locus2]);
                locus2++;
            }
            while ((locus1 < nbsel1) && (c1.incond[locus1] < pos[j]))
            {
                locus1++;
            }
        }
    }
    
    if (nbCo % 2 == 0)
        while (locus1 < nbsel1)
        {
            res.incond.push_back(c1.incond[locus1]);
            locus1++;
        }
    else
        while (locus2 < nbsel2)
        {
            res.incond.push_back(c2.incond[locus2]);
            locus2++;
        }

    sort(res.incond.begin(), res.incond.end());

    
    
    // loci affecting the selfing rate:
    locus = 0;
    for (j = 0; j < nbCo; j++)
    {
        
        if (j % 2 == 0)
            while ((locus < nself) && (posU[locus] < pos[j]))
            {
                res.selfing[locus] = c1.selfing[locus];
                locus++;
            }
        else
            while ((locus < nself) && (posU[locus] < pos[j]))
            {
                res.selfing[locus] = c2.selfing[locus];
                locus++;
            }
       
    }
    
    if (nbCo % 2 == 0)
        while (locus < nself)
        {
            res.selfing[locus] = c1.selfing[locus];
            locus++;
        }
    else
        while (locus < nself)
        {
            res.selfing[locus] = c2.selfing[locus];
            locus++;
        }
    
    
    // neutral loci:
    
    locus = 0;
    for (j = 0; j < nbCo; j++)
    {
        if (j % 2 == 0)
            while ((locus < nneut) && (posneut[locus] < pos[j]))
            {
                res.neut[locus] = c1.neut[locus];
                locus++;
            }
        else
            while ((locus < nneut) && (posneut[locus] < pos[j]))
            {
                res.neut[locus] = c2.neut[locus];
                locus++;
            }
    }
    
    if (nbCo % 2 == 0)
        while (locus < nneut)
        {
            res.neut[locus] = c1.neut[locus];
            locus++;
        }
    else
        while (locus < nneut)
        {
            res.neut[locus] = c2.neut[locus];
            locus++;
        }
}

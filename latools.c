/* lamap - fast two-way and three-way admixture likelihood functions
   Copyright (C) 2013 Giulio Genovese

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.   */

/* Written by Giulio Genovese <giulio.genovese@gmail.com> */

// lkl is a vector with likelihoods computed from ancestry proportions
// n is the number of samples
// frqA is the frequency of the alternate allele for ancestry A
// frqB is the frequency of the alternate allele for ancestry B
// frqC is the frequency of the alternate allele for ancestry C
// ancAA is the sample proportion of ancestry AA
// ancAB is the sample proportion of ancestry AB
// ancBB is the sample proportion of ancestry BB
// ancAC is the sample proportion of ancestry AC
// ancBC is the sample proportion of ancestry BC
// ancCC is the sample proportion of ancestry CC
// la is a local ancestry assignment
// gl0 is the likelihood of homozygous reference
// gl1 is the likelihood of heterozygous alternate
// gl2 is the likelihood of homozygous alternate

#include <stdio.h>
#include <stdlib.h>

// structure with emission frequencies for two ancestries
typedef struct {
  double aa0;
  double aa1;
  double aa2;
  double ab0;
  double ab1;
  double ab2;
  double bb0;
  double bb1;
  double bb2;
} snp_emis2_t;

// structure with emission frequencies for three ancestries
typedef struct {
  double aa0;
  double aa1;
  double aa2;
  double ab0;
  double ab1;
  double ab2;
  double bb0;
  double bb1;
  double bb2;
  double ac0;
  double ac1;
  double ac2;
  double bc0;
  double bc1;
  double bc2;
  double cc0;
  double cc1;
  double cc2;
} snp_emis3_t;

snp_emis2_t emis2(double frqA, double frqB) {
    snp_emis2_t f;
    f.aa0 = (1-frqA)*(1-frqA);
    f.aa1 = 2*frqA*(1-frqA);
    f.aa2 = frqA*frqA;
    f.ab0 = (1-frqA)*(1-frqB);
    f.ab1 = frqA*(1-frqB)+(1-frqA)*frqB;
    f.ab2 = frqA*frqB;
    f.bb0 = (1-frqB)*(1-frqB);
    f.bb1 = 2*frqB*(1-frqB);
    f.bb2 = frqB*frqB;
    return f;
}

snp_emis3_t emis3(double frqA, double frqB, double frqC) {
    snp_emis3_t f;
    f.aa0 = (1-frqA)*(1-frqA);
    f.aa1 = 2*frqA*(1-frqA);
    f.aa2 = frqA*frqA;
    f.ab0 = (1-frqA)*(1-frqB);
    f.ab1 = frqA*(1-frqB)+(1-frqA)*frqB;
    f.ab2 = frqA*frqB;
    f.bb0 = (1-frqB)*(1-frqB);
    f.bb1 = 2*frqB*(1-frqB);
    f.bb2 = frqB*frqB;
    f.ac0 = (1-frqA)*(1-frqC);
    f.ac1 = frqA*(1-frqC)+(1-frqA)*frqC;
    f.ac2 = frqA*frqC;
    f.bc0 = (1-frqB)*(1-frqC);
    f.bc1 = frqB*(1-frqC)+(1-frqB)*frqC;
    f.bc2 = frqB*frqC;
    f.cc0 = (1-frqC)*(1-frqC);
    f.cc1 = 2*frqC*(1-frqC);
    f.cc2 = frqC*frqC;
    return f;
}

// compute the likelihood for the genotypes from ancestry proportions
void lkl2(double *lkl, long n, double frqA, double frqB, double *ancAA, double *ancAB, double *ancBB, double *gl0, double *gl1, double *gl2) {
    snp_emis2_t f = emis2(frqA, frqB);
    long i;
    for (i=0; i<n; i++) {
	  lkl[i] = ( (f.aa0 * gl0[i] + f.aa1 * gl1[i] + f.aa2 * gl2[i]) * ancAA[i] +
                     (f.ab0 * gl0[i] + f.ab1 * gl1[i] + f.ab2 * gl2[i]) * ancAB[i] +
                     (f.bb0 * gl0[i] + f.bb1 * gl1[i] + f.bb2 * gl2[i]) * ancBB[i] ) /
                   (ancAA[i] + ancAB[i] + ancBB[i]);
    }
}

// compute the odds ratio for the genotypes from local ancestry assignments and ancestry proportions
double odd2(long n, double frqA, double frqB, char *la, double *gl0, double *gl1, double *gl2, double *lkl) {
    snp_emis2_t f = emis2(frqA, frqB);
    long i;
    double lod = 1;
    for (i=0; i<n; i++) {
        switch( la[i] ) {
        case 1: // AA
            lod *= ( f.aa0 * gl0[i] + f.aa1 * gl1[i] + f.aa2 * gl2[i] ) / lkl[i];
            break;
        case 2: // AB
            lod *= ( f.ab0 * gl0[i] + f.ab1 * gl1[i] + f.ab2 * gl2[i] ) / lkl[i];
            break;
        case 3: // BB
            lod *= ( f.bb0 * gl0[i] + f.bb1 * gl1[i] + f.bb2 * gl2[i] ) / lkl[i];
            break;
      }
    }
    return lod;
}

// compute the likelihood for the genotypes from ancestry proportions
void lkl3(double *lkl, long n, double frqA, double frqB, double frqC, double *ancAA, double *ancAB, double *ancBB, double *ancAC, double *ancBC, double *ancCC, double *gl0, double *gl1, double *gl2) {
    snp_emis3_t f = emis3(frqA, frqB, frqC);
    long i;
    for (i=0; i<n; i++) {
	  lkl[i] = ( (f.aa0 * gl0[i] + f.aa1 * gl1[i] + f.aa2 * gl2[i]) * ancAA[i] +
                     (f.ab0 * gl0[i] + f.ab1 * gl1[i] + f.ab2 * gl2[i]) * ancAB[i] +
                     (f.bb0 * gl0[i] + f.bb1 * gl1[i] + f.bb2 * gl2[i]) * ancBB[i] +
                     (f.ac0 * gl0[i] + f.ac1 * gl1[i] + f.ac2 * gl2[i]) * ancAC[i] +
                     (f.bc0 * gl0[i] + f.bc1 * gl1[i] + f.bc2 * gl2[i]) * ancBC[i] +
                     (f.cc0 * gl0[i] + f.cc1 * gl1[i] + f.cc2 * gl2[i]) * ancCC[i] ) /
                   (ancAA[i] + ancAB[i] + ancBB[i] + ancAC[i] + ancBC[i] + ancCC[i]);
    }
}

// compute the odds ratio for the genotypes from local ancestry assignments and ancestry proportions
double odd3(long n, double frqA, double frqB, double frqC, char *la, double *gl0, double *gl1, double *gl2, double *lkl) {
    snp_emis3_t f = emis3(frqA, frqB, frqC);
    long i;
    double lod = 1;
    for (i=0; i<n; i++) {
        switch( la[i] ) {
        case 1: // AA
            lod *= ( f.aa0 * gl0[i] + f.aa1 * gl1[i] + f.aa2 * gl2[i] ) / lkl[i];
            break;
        case 2: // AB
            lod *= ( f.ab0 * gl0[i] + f.ab1 * gl1[i] + f.ab2 * gl2[i] ) / lkl[i];
            break;
        case 3: // BB
            lod *= ( f.bb0 * gl0[i] + f.bb1 * gl1[i] + f.bb2 * gl2[i] ) / lkl[i];
            break;
        case 4: // AC
            lod *= ( f.ac0 * gl0[i] + f.ac1 * gl1[i] + f.ac2 * gl2[i] ) / lkl[i];
            break;
        case 5: // BC
            lod *= ( f.bc0 * gl0[i] + f.bc1 * gl1[i] + f.bc2 * gl2[i] ) / lkl[i];
            break;
        case 6: // CC
            lod *= ( f.cc0 * gl0[i] + f.cc1 * gl1[i] + f.cc2 * gl2[i] ) / lkl[i];
            break;
      }
    }
    return lod;
}

// fast code (4x improvement over load numpy.loadtxt) to load the local ancestry inference matrix
void getla(long *chrom, long *start, long *end, char* la, long n, long m, char** str) {
    long i, j, offset;
    for (i=0; i<n; i++) {
        if(sscanf(*str,"%ld%ld%ld%ln",chrom,start,end,&offset)!=3) {
            fprintf(stderr, "Wrong format at line %ld in the bedgraph file ...\n",i+2);
            exit(EXIT_FAILURE);
        }
        *str += offset;
        for (j=0; j<m; j++) {
            if(sscanf(*str,"%hhd%ln",la,&offset)!=1) {
                fprintf(stderr, "Wrong format at line %ld, column %ld in the bedgraph file ...\n",i+2,j+4);
                exit(EXIT_FAILURE);
            }
            *str += offset;
            la++;
        }
        chrom++; start++; end++; str++;
    }
}

void getsig(long *chrom, long *pos, char* sig, long T, long n, char** str) {
    long i, j, offset;
    char buffer[256];
    for (i=0; i<T; i++) {
        if(sscanf(*str,"%ld%ld%s%ln",chrom,pos,buffer,&offset)!=3) {
            fprintf(stderr, "Wrong format at line %ld in the signature file ...\n",i+1);
            exit(EXIT_FAILURE);
        }
        *str += offset;
        for (j=0; j<n; j++) {
            if(sscanf(*str,"%hhd%ln",sig,&offset)!=1) {
                fprintf(stderr, "Wrong format at line %ld, column %ld in the signature file ...\n",i+1,j+4);
                exit(EXIT_FAILURE);
            }
            *str += offset;
            sig++;
        }
        chrom++; pos++; str++;
    }
}

// compute the Viterbi path for the local ancestry deconvolution
void getpath(char *path, long T, long n, double *trans, char *phased, char *matobs, char *patobs) {
    long t, i, j, n2 = n*n;
    double prob[n2];
    double newprob[n2];
    char ptr[(T-1)*n2];
    double tmp;

    for (i=0; i<n2; i++)
        prob[i] = matobs[i/n] + patobs[i%n];
    for (t=1; t<T; t++) {
        for (i=0; i<n2; i++)
            for (j=0; j<n2; j++) {
                if (phased[t] || matobs[t*n+i/n] + patobs[t*n+i%n] < matobs[t*n+i%n] + patobs[t*n+i/n])
                    tmp = prob[j] + trans[j*n2+i] + (double)(matobs[t*n+i/n] + patobs[t*n+i%n]);
                else
                    tmp = prob[j] + trans[j*n2+i] + (double)(matobs[t*n+i%n] + patobs[t*n+i/n]);
                if (j==0 || tmp<newprob[i]) {
                    newprob[i] = tmp;
                    ptr[(t-1)*n2+i] = j;
                }
            }
        for (i=0; i<n2; i++)
            prob[i] = newprob[i];
    }
    for (i=0; i<n2; i++)
        if (i==0 || prob[i]<prob[path[T-1]])
            path[T-1] = i;

    for (t=T-2; t>=0; t--)
        path[t] = ptr[t*n2+path[t+1]];
}

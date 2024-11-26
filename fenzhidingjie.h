#ifndef _FENZHIDINGJIE_H_
#define _FENZHIDINGJIE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define R 10       // Define the number of rounds in the algorithm
#define N 16      // Define the number of s-boxes
#define s_in_L 4  // Define the number of input bits for the s-box
#define s_out_L 4 // Define the number of output bits for the s-box
#define L (1 << (s_in_L + s_out_L))
#define L1 (1 << s_in_L)
#define L2 (1 << s_out_L)
#define base (L2 - 1)

extern int Rr;
extern long double DiffP[R][N]; // Differential probability for each round
extern int Alpha[R][N];    // Input differential for each round
extern int Beta[R][N];     // Output differential for each round
extern int maxAlpha[R][N]; // Optimal input differential for each round
extern int maxBeta[R][N];  // Optimal output differential for each round
extern long double Pestim;
extern int flag;
extern long double all_sbdiff_p[N][L]; // Differential probability for all s-boxes
extern long double pn_sb[N];           // Optimal probability for n s-boxes
extern long double Pbest[R];           // Optimal probability for each round
extern int actIndex[R][N];        // Index of active s-boxes
extern int actIndexLength[R];     // Number of active s-boxes
extern long double sboxdiffprobability[L];
extern int maxbacktable[N][L2];
extern int maxforwardtable[N][L1];
extern int *combian;
extern int *combianSize;
extern int combianIndex;
extern int data[N];
extern long double *arr; 
extern int decreasort[L1][L2];
extern int pbox[64];

void Round(int r);           
void SubRound(int r, int n); 
void FirstRound(int M[], int n);
void init(); 
void calculateDifferential();
void all_sbdiff_pfun();
void Decreasesort();
void maxback();
void maxforward();
void pn_sbfun();
void calculateCombianAndCombianSize();
void generateCombinations(int n, int start, int *buffer, int bufferIndex);

int sumcombian(int n);
void activesbox(int begin, int end);
void update(int begin, int end);
long double product(int begin, int end);
long double fmaxfun(long double data[], int begin, int end);
void sortDoubleArray(long double array[], int size);
void sort_indices(long double *array, int *indices, int n);
int compare1(const void *a, const void *b);
int compare2(const void *a, const void *b);
void pboxfun(int *data);
long double LastRound();
int OptTrailSearch();
long double OptTrail();

#endif
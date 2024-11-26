#define printf __mingw_printf
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include "fenzhidingjie.h"

int main()
{
    
    init();
    calculateDifferential();
    all_sbdiff_pfun();
    Decreasesort();
    maxback();
    maxforward();
    calculateCombianAndCombianSize();
    pn_sbfun();
    Pbest[0] = pn_sb[0];

    clock_t start, end;
    double cpu_time_used = 0;

    for (Rr = 2; Rr < R + 1; Rr++)
    {
        start = clock();
        Pbest[Rr - 1] = OptTrail();
        printf("Pestime= %Le\n", Pestim);
        end = clock();
        cpu_time_used = cpu_time_used + ((double)(end - start)) / (CLOCKS_PER_SEC / 1000);

        printf("%d round Optimal Differential Characteristics:\n", Rr);
        for (int i = 0; i < Rr; i++)
        {
            printf("0x");
            for (int j = 0; j < N; j++)
            {
                printf("%X", maxAlpha[i][j]);
            }
            printf("\n");
            printf("0x");
            for (int j = 0; j < N; j++)
            {
                printf("%X", maxBeta[i][j]);
            }
            printf("\n");
        }
    }
    // print Pbest
    printf("Optimal Differential Transition Probability:\n");
    for (int i = 1; i < R; i++) {
        printf("%d: %.*Le\n",i+1, LDBL_DIG, Pbest[i]);
    }
    printf("%d round running time: %g ms\n", R, cpu_time_used);
    free(combian);
    free(combianSize);
    // printf("Size of long double: %zu bytes\n", sizeof(long double));
    // printf("Precision of long double: %d decimal digits\n", LDBL_DIG);
    // printf("Minimum positive value of long double: %Le\n", LDBL_MIN);
    // printf("Maximum positive value of long double: %Le\n", LDBL_MAX);
    // printf("Epsilon value of long double: %Le\n", LDBL_EPSILON);
    return 0;
}
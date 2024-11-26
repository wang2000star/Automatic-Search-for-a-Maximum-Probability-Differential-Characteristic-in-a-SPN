#include "fenzhidingjie.h"

int Rr;
long double DiffP[R][N]; // Differential probability for each round
int Alpha[R][N];         // Input differential for each round
int Beta[R][N];          // Output differential for each round
int maxAlpha[R][N];      // Optimal input differential for each round
int maxBeta[R][N];       // Optimal output differential for each round
long double Pestim;
int flag;
long double all_sbdiff_p[N][L]; // Differential probability for all s-boxes
long double pn_sb[N];           // Optimal probability for n s-boxes
long double Pbest[R];           // Optimal probability for each round
int actIndex[R][N];             // Index of active s-boxes
int actIndexLength[R];          // Number of active s-boxes
long double sboxdiffprobability[L];
int maxbacktable[N][L2];
int maxforwardtable[N][L1];
int *combian;
int *combianSize;
int combianIndex = 0;
int data[N];
long double *arr;
int decreasort[L1][L2];

// Comparison function for qsort
int compare1(const void *a, const void *b)
{
    int index1 = *(int *)a;
    int index2 = *(int *)b;
    if (arr[index1] < arr[index2])
        return 1;
    if (arr[index1] > arr[index2])
        return -1;
    return 0;
}
// Comparison function
int compare2(const void *a, const void *b)
{
    long double arg1 = *(const long double *)a;
    long double arg2 = *(const long double *)b;
    if (arg1 < arg2)
        return -1;
    if (arg1 > arg2)
        return 1;
    return 0;
}
// Sorting function
void sortDoubleArray(long double array[], int size)
{
    qsort(array, size, sizeof(long double), compare2);
}
long double fmaxfun(long double data[], int begin, int end)
{
    long double max = data[begin];
    for (int i = begin; i <= end; i++)
    {
        if (data[i] > max)
        {
            max = data[i];
        }
    }
    return max;
}
void sort_indices(long double *array, int *indices, int n)
{
    arr = array; // Set global variable
    for (int i = 0; i < n; i++)
    {
        indices[i] = i;
    }
    qsort(indices, n, sizeof(int), compare1);
}
void generateCombinations(int n, int start, int *buffer, int bufferIndex)
{
    if (bufferIndex == n)
    {
        for (int i = 0; i < n; i++)
        {
            combian[combianIndex++] = buffer[i];
        }
        return;
    }
    for (int i = start; i < N; i++)
    {
        buffer[bufferIndex] = i;
        generateCombinations(n, i + 1, buffer, bufferIndex + 1);
    }
}
void calculateCombianAndCombianSize()
{
    int totalSize = N * (1 << (N - 1));        // Calculate the total size of combian
    combian = malloc(totalSize * sizeof(int)); // Allocate size for combian
    combianSize = malloc(N * sizeof(int));     // Store the number of combinations for each n
    int prevIndex = 0;                         // Used to calculate the number of combinations for each n
    for (int n = 1; n <= N; n++)
    {
        int *buffer = malloc(n * sizeof(int));
        generateCombinations(n, 0, buffer, 0);
        combianSize[n - 1] = (combianIndex - prevIndex) / n; // Store the number of combinations for n
        prevIndex = combianIndex;                            // Update prevIndex
        free(buffer);
    }
}
int sumcombian(int n)
{
    int sum = 0;
    if (n == 1)
    {
        return 0;
    }
    for (int i = 0; i < n - 1; i++)
    {
        sum = sum + (i + 1) * combianSize[i];
    }
    return sum;
}
void init()
{
    for (int i = 0; i < R; i++)
    {
        for (int j = 0; j < N; j++)
        {
            Alpha[i][j] = 0;
            Beta[i][j] = 0;
            DiffP[i][j] = 1.0;
            actIndex[i][j] = j;
        }
        actIndexLength[i] = N;
    }
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            all_sbdiff_p[i][j] = 0.0;
        }
    }
}
//Assuming that the s-box is m*n, then sboxdiffprobability is 2^m*2^n, where p(i->j)=sboxdiffprobability[2^m*i+j]
int sbox[L1] = {0xC, 0x5, 0x6, 0xB, 0x9, 0x0, 0xA, 0xD, 0x3, 0xE, 0xF, 0x8, 0x4, 0x7, 0x1, 0x2};
void sboxfun(int *data)
{
    // data size is 2, each occupies two bits, a total of m bits, replace each one with sbox
    for (int i = 0; i < (1 << L1); i++)
    {
        data[i] = sbox[data[i]];
    }
}
long double calculateDifferentialProbability(int inputDiff, int outputDiff)
{
    int count = 0;
    for (int x = 0; x < L1; x++)
    {
        int y = x ^ inputDiff;
        if (outputDiff == (sbox[x] ^ sbox[y]))
            count++;
    }
    long double probability = (long double)count / L1;
    return probability;
}
void calculateDifferential()
{
    for (int i = 0; i < L1; i++)
    {
        for (int j = 0; j < L2; j++)
        {
            sboxdiffprobability[i * L1 + j] = calculateDifferentialProbability(i, j);
        }
    }
}
void all_sbdiff_pfun()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            all_sbdiff_p[i][j] = sboxdiffprobability[j];
        }
    } 
 // note: the s-boxes are the same, for simplicity, only one can be used, but for generality, we use N s-boxes
}
long double product(int begin, int end)
{
    long double p = 1.0;
    for (int i = begin - 1; i < end; i++)
    {
        for (int j = 0; j < actIndexLength[i]; j++)
        {
            p = p * DiffP[i][actIndex[i][j]];
        }
    }
    return p;
}
void update(int begin, int end)
{
    int index1, index2;
    for (int r = begin - 1; r < end; r++)
    {
        for (int i = 0; i < actIndexLength[r]; i++)
        {
            index1 = actIndex[r][i];
            index2 = Beta[r][index1] + Alpha[r][index1] * L2;
            DiffP[r][index1] = all_sbdiff_p[index1][index2];
        }
    }
}
void activesbox(int begin, int end)
{
    for (int r = begin - 1; r < end; r++)
    {
        int count = 0;
        for (int i = 0; i < N; i++)
        {
            if (Alpha[r][i] != 0)
            {
                actIndex[r][count] = i;
                count++;
            }
        }
        actIndexLength[r] = count;
    }
}
void Decreasesort()
{
    long double array[L2];
    int *indices = malloc(L2 * sizeof(int)); // index of sorted array
    for (int i = 0; i < L1; i++)
    {
        for (int j = 0; j < L2; j++)
        {
            array[j] = sboxdiffprobability[i * L1 + j];
        }
        sort_indices(array, indices, L2);
        for (int k = 0; k < L2; k++)
        {
            decreasort[i][k] = indices[k];
        }
    }
    free(indices);
}
void pn_sbfun()
{
    long double temp[N];
    for (int i = 0; i < N; i++)
    {
        temp[i] = fmaxfun(sboxdiffprobability, (1 << s_out_L), L - 1);
    }
    sortDoubleArray(temp, N);
    pn_sb[0] = temp[N - 1];
    for (int i = 1; i < N; i++)
    {
        pn_sb[i] = pn_sb[i - 1] * temp[N - i - 1];
    }
}
void maxback()
{
    for (int sbox_index = 0; sbox_index < N; sbox_index++)
    {
        for (int beta = 0; beta < L2; beta++)
        {
            int alpha = 0;
            long double pmax = 0;
            for (int i = 0; i < L1; i++)
            {
                if (all_sbdiff_p[sbox_index][i * L2 + beta] >= pmax)
                {
                    alpha = i;
                    pmax = all_sbdiff_p[sbox_index][beta + i * L2];
                }
            }
            maxbacktable[sbox_index][beta] = alpha;
        }
    }
    return;
}
void maxforward()
{
    for (int sbox_index = 0; sbox_index < N; sbox_index++)
    {
        for (int alpha = 0; alpha < L1; alpha++)
        {
            int beta = 0;
            long double pmax = 0;
            for (int i = 0; i < (1 << s_out_L); i++)
            {
                if (all_sbdiff_p[sbox_index][i + alpha * L1] >= pmax)
                {
                    beta = i;
                    pmax = all_sbdiff_p[sbox_index][i + alpha * L1];
                }
            }
            maxforwardtable[sbox_index][alpha] = beta;
        }
    }
    return;
}
int pbox[64] = {0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51, 4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55, 8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59, 12, 28, 44, 60, 13, 29, 45, 61, 14, 30, 46, 62, 15, 31, 47, 63};
void pboxfun(int *data)
{
    int bits[64];
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            bits[i * 4 + j] = (data[i] >> (3 - j)) & 1;
        }
    }
    int temp[64];
    for (int i = 0; i < 64; i++)
    {
        temp[i] = bits[i];
    }
    for (int i = 0; i < 64; i++)
    {
        bits[pbox[i]] = temp[i];
    }
    for (int i = 0; i < 16; i++)
    {
        data[i] = 0;
        for (int j = 0; j < 4; j++)
        {
            data[i] |= bits[i * 4 + j] << (3 - j);
        }
    }
}
long double LastRound()
{
    if (actIndexLength[Rr - 1] > 0)
    {
        int index1 = 0, index2 = 0;
        for (int i = 0; i < N; i++)
        {
            Beta[Rr - 1][i] = 0;
        }
        for (int i = 0; i < actIndexLength[Rr - 1]; i++)
        {
            index1 = actIndex[Rr - 1][i];
            Beta[Rr - 1][index1] = maxforwardtable[index1][Alpha[Rr - 1][index1]];
        }
        long double p = 1.0;
        int m;
        for (int i = 0; i < Rr; i++)
        {
            for (int j = 0; j < actIndexLength[i]; j++)
            {
                index1 = actIndex[i][j];
                index2 = Beta[i][index1] + Alpha[i][index1] * L2;
                p = p * all_sbdiff_p[index1][index2];
            }
        }
        if (p >= Pestim)
        {
            printf("p*:%.10Lf\n", p);
            printf("Pestim*:%.10Lf\n", Pestim);
            flag = 0;
            Pestim = p;
            activesbox(1, Rr);
            for (int i = 0; i < Rr; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    maxAlpha[i][j] = Alpha[i][j];
                    maxBeta[i][j] = Beta[i][j];
                }
            }
        }
    }
    else
    {
        printf("all zero error\n");
    }
}
void Round(int r)
{
    for (int i = 0; i < N; i++)
    {
        Beta[r - 1][i] = 0;
    }
    activesbox(1, r);
    SubRound(r, 1);
}
void SubRound(int r, int n)
{
    if (n > actIndexLength[r - 1])
    {
        int data[N];
        for (int i = 0; i < N; i++)
        {
            data[i] = Beta[r - 1][i];
        }
        pboxfun(data); // linear transformation
        for (int i = 0; i < N; i++)
        {
            Alpha[r][i] = data[i];
        }
        activesbox(1, r + 1);
        update(1, r + 1);
        if (r + 1 < Rr)
        {
            Round(r + 1);
        }
        else
        {
            LastRound();
        }
    }
    else
    {
        for (int i = 0; i < L2; i++)
        {
            Beta[r - 1][actIndex[r - 1][n - 1]] = decreasort[Alpha[r - 1][actIndex[r - 1][n - 1]]][i];
            long double p = 1;
            for (int j = 0; j < n; j++)
            {
                p = p * all_sbdiff_p[actIndex[r - 1][j]][Beta[r - 1][actIndex[r - 1][j]] + Alpha[r - 1][actIndex[r - 1][j]] * L2];
            }
            for (int i = 0; i < r - 1; i++)
            {
                for (int j = 0; j < actIndexLength[i]; j++)
                {
                    p = p * all_sbdiff_p[actIndex[i][j]][Beta[i][actIndex[i][j]] + Alpha[i][actIndex[i][j]] * L2];
                }
            }
            p = p * Pbest[Rr - r - 1];
            if (actIndexLength[r - 1] > n)
            {
                p = p * pn_sb[actIndexLength[r - 1] - n - 1];
            }
            if (p < Pestim)
            {
                break;
            }
            int data[N];
            for (int i = 0; i < N; i++)
            {
                data[i] = Beta[r - 1][i];
            }
            pboxfun(data); // linear transformation
            int sum = 0;
            for (int i = 0; i < 16; i++)
            {
                if (data[i] != 0)
                {
                    sum++;
                }
            }
            p = p * pn_sb[sum - 1];
            if (Rr > r + 1)
            {
                p = p * Pbest[Rr - r - 2];
            }
            if (p < Pbest[Rr - r - 1] * Pestim)
            {
                break;
            }
            else
            {
                update(r, r);
                activesbox(1, r);
                SubRound(r, n + 1);
            }
        }
    }
    return;
}
void FirstRound(int M[], int n)
{
    int sbox_index;
    for (int i = 0; i < n; i++)
    {
        sbox_index = M[i];
        Alpha[0][sbox_index] = maxbacktable[sbox_index][Beta[0][sbox_index]];
    }
    for (int i = 0; i < N; i++)
    {
        data[i] = Beta[0][i];
    }
    pboxfun(data); // linear transformation
    for (int i = 0; i < N; i++)
    {
        Alpha[1][i] = data[i];
    }
    if (Rr > 2)
    {
        activesbox(1, 2);
        update(1, 1);
        Round(2);
    }
    else
    {
        for (int i = 0; i < N; i++)
        {
            Beta[Rr - 1][i] = 0;
        }
        activesbox(1, Rr);
        LastRound();
    }
    return;
}
int OptTrailSearch() // r>=2
{
    int begin;
    long limit;
    long temp;
    for (int n = 1; n < N + 1; n++)
    {
        if (pn_sb[n - 1] * Pbest[Rr - 2] < Pestim)
        {
            break;
        }
        else
        {
            for (int i = 0; i < N; i++)
            {
                Beta[0][i] = 0;
                Alpha[0][i] = 0;
            }
            begin = sumcombian(n);
            int M[n];
            for (int i = 0; i < combianSize[n - 1]; i++)
            {
                for (int j = 0; j < n; j++)
                    M[j] = combian[begin + n * i + j];
            }
            limit = powl(base, n);
            for (long k = 0; k < limit; k++)
            {
                temp = k;
                for (int index = 0; index < n; index++)
                {
                    Beta[0][M[index]] = temp % base + 1;
                    temp = temp / base;
                }
                FirstRound(M, n);
            }
        }
        printf("n=%d\n", n);
    }
    return flag;
}
long double OptTrail()
{
    Pestim = Pbest[Rr - 2] * pn_sb[0];
    flag = -1;
    flag = OptTrailSearch();
    printf("flag=%d\n", flag);
    while (1)
    {
        if (~flag)
        {
            break;
        }
        else
        {
            printf("Pestim*=%Le\n", Pestim);
            Pestim = Pestim / 2;
            printf("Pestim*=%Le\n", Pestim);
            flag = OptTrailSearch();
            printf("flag=%d\n", flag);
        }
    }
    return Pestim;
}
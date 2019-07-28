#include <time.h>
#include <iostream>
#include <mpich/mpi.h>
#include <fstream>

#include "gf2misc.h"
#include "mqmisc.h"

#define SEARCH_SPACE 57

FILE *fr = fopen("challenge-1-74-0-resident", "rb");
poly fullpoly[M][N + 1];
__int128 partialDerivative[M][N] = {0};

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

int main() {
    ffile2poly(fr, fullpoly, EQUATION_NUM);
    loadPD(fullpoly, partialDerivative);
    fclose(fr);

    __int128 value[EQUATION_NUM] = {0};
    int32_t matrix[EQUATION_NUM] = {0};
    int64_t z;
    int64_t pre = 0;
    int64_t los = 0;
    int64_t guessMask = 0x1FFFFFFFFFFFFFF;

    // MPI initialization
    int rank, size;
    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int64_t final_key = 0x200000000000000;

    // key initialization
    int64_t init_key = final_key * rank / 10000;
    int64_t end_key = final_key * (rank + 1) / 10000;
    for (int i = 0; i < EQUATION_NUM; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < fullpoly[i][j].length; k++) {
                fullpoly[i][j].p[k].data[j / 32] ^= (1 << (j % 32));
            }
    for (int i = 0; i < EQUATION_NUM; i++) {
        for (int j = 0; j < SEARCH_SPACE; j++) {
            if ((init_key >> j) & 1)
                for (int k = 0; k < fullpoly[i][j].length; k++)
                    value[i] ^= (fullpoly[i][j].p[k].data[0] ^ ((__int128) fullpoly[i][j].p[k].data[1] << 32)
                        ^ ((__int128) fullpoly[i][j].p[k].data[2] << 64));
        }
    }

    for (int64_t key = init_key; key < end_key; key++) {
        memset(matrix, 0, EQUATION_NUM * 4);
        uint32_t clist[11] = {0};
        z = (key >> 1) ^ key ^ pre;
        los = __builtin_ffsl(z);
        for (int j = 0; j < EQUATION_NUM; j++) {
            value[j] ^= partialDerivative[j][los] & pre;
            matrix[j] = __builtin_parityl(value[j] & guessMask) ^ ((__int128) value[j] >> 67);
            matrix[j] |= (value[j] >> 56) & 0x7FE;
        }
        for (int i = 0; i < VARIABLE_NUM + 1; i++) {
            for (int j = 0; j < EQUATION_NUM; j++) {
                clist[i] |= ((matrix[j] >> (10 - i)) & 1) << j;
            }
        }
        uint32_t mask = 0x3FFFFF;
        checkConsist_22x10(clist, mask);
        if (!(mask & clist[10])) { // not singular
            std::ofstream file;
            file.open("solution"+currentDateTime()+".txt");
            file << key << std::endl;
            for (auto col: clist)
                file << col << std::endl;
            file.close();
        }
    }
    MPI_Finalize();

    return EXIT_SUCCESS;
}
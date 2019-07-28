#include <time.h>
#include <iostream>

#include "gf2misc.h"
#include "mqmisc.h"

FILE *fr = fopen("challenge-1-74-0-resident", "rb");
poly fullpoly[M][N + 1];
int64_t partialDerivative[M][N] = {0};

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

    __int128 value[EQUATION_NUM] = {0};
    int32_t matrix[EQUATION_NUM] = {0};
    int64_t z;
    int64_t pre = 0;
    int64_t los = 0;

    int64_t init_key;
    int64_t end_key;
    for (int i = 0; i < EQUATION_NUM; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < fullpoly[i][j].length; k++) {
                fullpoly[i][j].p[k].data[j / 32] ^= (1 << (j % 32));
            }
    for (int i = 0; i < EQUATION_NUM; i++) {
        for (int j = 0; j < 58; j++) {
            if ((init_key >> j) & 1)
                for (int k = 0; k < fullpoly[i][j].length; k++)
                    value[i] ^= (fullpoly[i][j].p[k].data[0] ^ ((__int128) fullpoly[i][j].p[k].data[1] << 32)
                        ^ ((__int128) fullpoly[i][j].p[k].data[2] << 64));
        }
    }

    for (int64_t key = init_key; key < end_key; key++) {
        memset(matrix, 0, EQUATION_NUM * 4);
        z = (key >> 1) ^ key ^ pre;
        los = __builtin_ffsl(z);
        for (int j = 0; j < EQUATION_NUM; j++)
            value[j] ^= partialDerivative[j][los] & pre;

    }

    return EXIT_SUCCESS;
}
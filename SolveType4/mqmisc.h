#ifndef PROJECT_MQMISC_H
#define PROJECT_MQMISC_H

#include <string.h>
#include "mqparam.h"

// MQ polynomial arithmetic
void addPoly(poly &dstPoly, poly &poly);
void addTerm(poly &dstPoly, term &term);
void mulPoly(poly &dstPoly, poly &poly1, poly &poly2);
void file2poly(FILE *fr, poly destPoly[M][N + 1]);
void poly2file(FILE *fw, poly poly);
void repeatPoly(poly &spoly);
void ffile2poly(FILE *fr, poly poly[][N + 1], int m);
void printpoly(FILE *fw, poly poly);
void simplifyPoly(poly poly[M][N + 1]);
bool verifyPoly(int64_t guess[N], poly poly[M][N + 1]);
int getLos(int64_t x);

inline void loadPD(poly fullpoly[M][N + 1], __int128 partialDerivative[M][N]) {
    for (int i = 0; i < EQUATION_NUM; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < fullpoly[i][j].length; k++) {
                int los = j;
                for (int l = j + 1; l < N; l++) {
                    int word = l / 32;
                    int bit = l % 32;
                    if ((fullpoly[i][j].p[k].data[word] >> bit) & 1) {
                        los = l;
                        break;
                    }
                }
                if (los > j) {
                    partialDerivative[i][j] ^= (int64_t(1) << los);
                    partialDerivative[i][los] ^= (int64_t(1) << j);
                } else
                    partialDerivative[i][j] ^= (int64_t(1) << N);
            }
        }
    }
}

#endif //PROJECT_MQMISC_H

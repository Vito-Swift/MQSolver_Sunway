#include <stdint.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "../ttmath-0.9.3/ttmath/ttmath.h"
#include "mpi.h"

typedef ttmath::Int<1> int32;
typedef ttmath::Int<2> int64;
typedef ttmath::Int<4> int128;

// structure to store one term
struct term {
  uint32_t data[3];
};

// structure to store a MQ polynomial
struct poly {
  // term array
  term *p;
  int length;
};

#define M 67
#define N 67
#define LEN 10001
#define VARIABLE_NUM 10
#define EQUATION_NUM 22
#define SEARCH_SPACE 57

// mq arithmetic
void repeatPoly(poly &spoly);
void ffile2poly(FILE *fr, poly poly[][N + 1], int m);
void addPoly(poly &dstPoly, poly &poly);
void addTerm(poly &dstPoly, term &term);
void loadPD(poly fullpoly[M][N + 1], int128 partialDerivative[M][N]);

// gf2 arithmetic
void checkConsist_22x10(uint32_t clist[11], uint32_t &mask);
const std::string currentDateTime();
const int parityCheck(const int128 v);

// Main loop
int main() {
    FILE *fr = fopen("mq-resident4-67-1.txt", "rb");
    poly fullpoly[M][N + 1];
    int128 partialDerivative[M][N] = {0};
    ffile2poly(fr, fullpoly, EQUATION_NUM);
    loadPD(fullpoly, partialDerivative);
    fclose(fr);

    int128 value[EQUATION_NUM] = {0};
    int matrix[EQUATION_NUM] = {0};
    int64 z;
    int64 pre = 0;
    int los = 0;
    int64 guessMask = "144115188075855871"; // 0x1FFFFFFFFFFFFFF

    // MPI Init
    int rank, size;
    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Key Boundary Init
    int64 final_key = "144115188075855872"; // 0x200000000000000
    int64 init_key = (final_key / 10000) * rank;
    int64 end_key = (final_key / 10000) * (rank + 1);

    // Equation Init
    for (int i = 0; i < EQUATION_NUM; i++)
        for (int j = 0; j < N; j++)
            for (int k = 0; k < fullpoly[i][j].length; k++)
                fullpoly[i][j].p[k].data[j / 32] ^= (1 << (j % 32));
    for (int i = 0; i < EQUATION_NUM; i++)
        for (int j = 0; j < SEARCH_SPACE; j++)
            if (((init_key >> j) & 1) == 1)
                for (int k = 0; k < fullpoly[i][j].length; k++)
                    value[i] ^= ((int32) fullpoly[i][j].p[k].data[0] ^ ((int128) fullpoly[i][j].p[k].data[1] << 32)
                        ^ ((int128) fullpoly[i][j].p[k].data[2] << 64));

    if (init_key != 0)
        pre = (init_key - 1) ^ ((init_key - 1) >> 1);

    for (int64 key = init_key; key < end_key; key++) {
        memset(matrix, 0, EQUATION_NUM * 4);
        uint32_t clist[11] = {0};
        z = (key >> 1) ^ key ^ pre;
        ttmath::uint tableid, index;
        z.FindLowestBit(tableid, index);
        los = tableid * 32 + index;
        for (int j = 0; j < EQUATION_NUM; j++) {
            value[j] ^= partialDerivative[j][los] & pre;
            int128 c = value[j] >> 67;
            int128 v = value[j] >> 56;
            uint32_t cuint, vuint;
            c.ToUInt(cuint);
            v.ToUInt(vuint);
            matrix[j] = parityCheck(value[j] & guessMask) ^ cuint;
            matrix[j] |= vuint & 0x7FE;
        }
        for (int i = 0; i < VARIABLE_NUM + 1; i++)
            for (int j = 0; j < EQUATION_NUM; j++)
                clist[i] |= ((matrix[j] >> (10 - i)) & 1) << j;

        uint32_t mask = 0x3FFFFF;
        checkConsist_22x10(clist, mask);
        if (!(mask & clist[10])) {
            std::ofstream file;
            std::string filename =
                "/home/export/base/suntnt/Vito/MQSolver_Sunway/solution" + std::to_string(rank) + currentDateTime()
                    + ".txt";
            file.open(filename.c_str());
            file << key << std::endl;
            for (int i = 0; i < 11; i++)
                file << clist[i] << std::endl;
            file.close();
        }
        pre = key ^ (key >> 1);
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
}

// Implementations of essential functions
const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}
void loadPD(poly fullpoly[M][N + 1], int128 partialDerivative[M][N]) {
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
                    partialDerivative[i][j] ^= (int64(1) << los);
                    partialDerivative[i][los] ^= (int64(1) << j);
                } else
                    partialDerivative[i][j] ^= (int64(1) << N);
            }
        }
    }
}
void addPoly(poly &poly1, poly &poly2) {
    if (poly2.length > 0) {
        int len = poly1.length + poly2.length;
        if ((poly1.p = (term *) realloc(poly1.p, len * sizeof(term))) == NULL) {
            printf("Out of Memory!");
            exit(1);
        }
        for (int i = 0; i < (poly2.length); i++) {
            for (int k = 0; k < 3; k++)
                poly1.p[i + poly1.length].data[k] = poly2.p[i].data[k];
        }
        poly1.length = len;
    }
}
void addTerm(poly &poly1, term &term1) {
    int len = poly1.length + 1;
    if ((poly1.p = (term *) realloc(poly1.p, len * sizeof(term))) == NULL) {
        printf("Out of Memory!");
        exit(1);
    }
    for (int k = 0; k < 3; k++)
        poly1.p[poly1.length].data[k] = term1.data[k];
    poly1.length = len;
}
void repeatPoly(poly &spoly) {
    int llen = 0;
    int len = 0;
    poly Chain[LEN];
    for (int i = 0; i < LEN; i++) {
        Chain[i].p = NULL;
        Chain[i].length = 0;
    }
    for (int s = 0; s < spoly.length; s++) {
        poly p1;
        p1.length = 1;
        p1.p = (term *) malloc(p1.length * sizeof(term));
        for (int i = 0; i < 3; i++)
            p1.p[0].data[i] = spoly.p[s].data[i];
        uint32_t hash = 0;
        for (int i = 0; i < 3; i++)
            hash = (hash + (p1.p[0].data[i] % LEN)) % LEN;
        if (Chain[hash].length == 0) {
            addPoly(Chain[hash], p1);
        } else {
            bool flag = 1;
            for (int j = 0; j < Chain[hash].length; j++) {
                bool tag = 1;
                for (int k = 0; k < 3; k++) {
                    if (Chain[hash].p[j].data[k] != p1.p[0].data[k]) {
                        tag = 0;
                        break;
                    }
                }
                if (tag) {
                    for (int k = j; k < Chain[hash].length - 1; k++) {
                        for (int l = 0; l < 3; l++)
                            Chain[hash].p[k].data[l] = Chain[hash].p[k + 1].data[l];
                    }
                    Chain[hash].length = Chain[hash].length - 1;
                    flag = 0;
                    break;
                }
            }
            if (flag)
                addPoly(Chain[hash], p1);
        }
        p1.length = 0;
        free(p1.p);
        p1.p = NULL;
    }
    for (int i = 0; i < LEN; i++) {
        if (Chain[i].length > 0) {
            for (int j = 0; j < Chain[i].length; j++) {
                for (int k = 0; k < 3; k++)
                    spoly.p[len].data[k] = Chain[i].p[j].data[k];
                len++;
            }
        }
        Chain[i].length = 0;
        free(Chain[i].p);
    }
    spoly.length = len;
}
void checkConsist_22x10(uint32_t clist[11], uint32_t &mask) {
    for (int i = 0; i < 10; i++) {
        uint32_t ci = clist[i] & mask;
        if (ci == 0) continue;
        uint32_t x = __builtin_clz(ci);
        uint32_t xp = (0x80000000 >> x);
        uint32_t xormask = ci ^xp;
        for (int j = 0; j < 11; j++)
            clist[j] = (clist[j] & xp) ? (clist[j] ^ xormask) : clist[j];
        mask ^= xp;
    }
}
const int parityCheck(const int128 v) {
    uint32_t lo = 0;
    uint32_t hi = 0;
    int128 hiv = v >> 64;
    v.ToUInt(lo);
    hiv.ToUInt(hi);
    return __builtin_parity(hi) ^ __builtin_parity(lo);
}

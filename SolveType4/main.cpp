#include <stdint.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <bitset>
#include "../ttmath-0.9.3/ttmath/ttmath.h"

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
    int128 partialDerivative[EQUATION_NUM][N] = {0};
    ffile2poly(fr, fullpoly, EQUATION_NUM);
    loadPD(fullpoly, partialDerivative);
    fclose(fr);

    int128 value[EQUATION_NUM] = {0};
    int matrix[EQUATION_NUM] = {0};
    int128 z;
    int128 pre = ((int128) 0x7FF << 57);
    int los = 0;
    int128 guessMask = "144115188075855871"; // 0x1FFFFFFFFFFFFFF

    // MPI Init
    int rank = 0;

    // Key Boundary Init
    int64 final_key = "144115188075855872"; // 0x200000000000000
    int64 init_key = (final_key / 10000) * rank;
    int64 end_key = (final_key / 10000) * (rank + 1);

    // Equation Init
    for (int i = 0; i < EQUATION_NUM; i++) {
        for (int j = 0; j < N; j++)
            for (int k = 0; k < fullpoly[i][j].length; k++) {
                int los = k;
                for (int ll = k + 1; ll < N; ll++) {
                    if ((fullpoly[i][j].p[k].data[ll / 32] >> (ll % 32)) & 1) {
                        los = ll;
                        break;
                    }
                }
                if (los > k)
                    value[i] ^= ((int128) 1 << los);
                else
                    value[i] ^= ((int128) 1 << N);
            }
        value[i] ^= ((int128) fullpoly[i][N].length << N);
    }

    if (init_key != 0)
        pre = ((int128) 0x7FF << 57) | (((int128) init_key - 1) ^ (((int128) init_key - 1) >> 1));

    for (int64 key = init_key; key < 10; key++) {
        memset(matrix, 0, EQUATION_NUM * 4);
        uint32_t clist[11] = {0};
        z = (key >> 1) ^ key ^ pre;

        // equiv to los = _bitScanForward(z);
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
        for (auto& c: clist)
            std::cout << std::bitset<22>(c) << std::endl;
        std::cout << std::endl;
        if (!(mask & clist[10])) {
//            std::ofstream file;
//            std::string filename =
//                "/home/export/base/suntnt/Vito/MQSolver_Sunway/solution" + std::to_string(rank) + currentDateTime()
//                    + ".txt";
//            file.open(filename.c_str());
//            file << key << std::endl;
//            for (int i = 0; i < 11; i++)
//                file << clist[i] << std::endl;
//            file.close();
        }
        pre = ((int128) key ^ ((int128) key >> 1)) | ((int128) 0x7FF << 57);
    }
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
void ffile2poly(FILE *fr, poly spoly[][N + 1], int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < N + 1; j++) {
            fscanf(fr, "%u", &(spoly[i][j].length));
            spoly[i][j].p = (term *) malloc(spoly[i][j].length * sizeof(term));
            for (int k = 0; k < spoly[i][j].length; k++)
                for (int l = 0; l < 3; l++)
                    fscanf(fr, "%u", &(spoly[i][j].p[k].data[l]));
        }
    }
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
                    partialDerivative[i][j] ^= (int128(1) << los);
                    partialDerivative[i][los] ^= (int128(1) << j);
                } else
                    partialDerivative[i][j] ^= (int128(1) << N);
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
    uint32_t v1 = 0;
    uint32_t v2 = 0;
    uint32_t v3 = 0;
    uint32_t v4 = 0;
    int128 hv1 = v >> 32;
    int128 hv2 = v >> 64;
    int128 hv3 = v >> 96;
    v.ToUInt(v4);
    hv1.ToUInt(v3);
    hv2.ToUInt(v2);
    hv3.ToUInt(v1);
    return __builtin_parity(v1) ^ __builtin_parity(v2) ^ __builtin_parity(v3) ^ __builtin_parity(v4);
}

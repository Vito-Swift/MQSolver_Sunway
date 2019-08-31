#ifndef PROJECT_MAIN_H
#define PROJECT_MAIN_H

#include <iostream>
#include <stdint.h>
#include <time.h>
#include <fstream>
#include <string>
#include <bitset>
#include <cstring>

#include <emmintrin.h>
#include <immintrin.h>

#include "../ttmath-0.9.3/ttmath/ttmath.h"

typedef ttmath::Int<1> ttint64_t;
typedef ttmath::Int<2> ttint128_t;

#define M 72
#define N 72
#define VARIABLE_NUM 16
#define EQUATION_NUM 22
#define RESIDENT_NUM 56

struct poly {
  ttint128_t *p;
  uint32_t length;
};

static int rank;
static int size;
const std::string ResidentFilename = "mq-resident1-72-1.txt";
const std::string ChallengeFilename = "challenge-1-72-1";
static poly ResidentPoly[M][N + 1];
static poly VerifyPoly[M][N + 1];
static uint32_t ColPartialDrt[N + 1][RESIDENT_NUM];
uint32_t ColVal[N + 1];
static uint64_t startKey;
static uint64_t endKey;
uint64_t binKey;
uint64_t grayCodeKey;

inline uint32_t ttmath_lzcntll(ttint128_t v) {
    ttint128_t hiv = v >> 64;
    uint64_t hi, lo;
    hiv.ToUInt(hi);
    v.ToUInt(lo);
    lo = (hi == 0) ? lo : -1ULL;
    return _lzcnt_u64(hi) + _lzcnt_u64(lo);
}

inline uint32_t ttmath_parityll(ttint128_t v) {
    ttint128_t hi = v >> 64;
    uint64_t h, l;
    v.ToUInt(l);
    hi.ToUInt(h);
    return __builtin_parityl(l) ^ __builtin_parityl(h);
}

inline uint32_t ttmath_ffsll(ttint128_t v) {
    uint64_t table, index;
    v.FindLowestBit(table, index);
    return (uint32_t) (table * 64 + index);
}

inline void addTerm(poly &dstPoly, ttint128_t term) {
    uint32_t len = dstPoly.length + 1;
    if ((dstPoly.p = (ttint128_t *) realloc(dstPoly.p, len * sizeof(ttint128_t))) == nullptr) {
        std::cerr << "Out of memory!" << std::endl;
        exit(1);
    }
    dstPoly.p[dstPoly.length] = term;
    dstPoly.length = len;
}

inline void file2poly(FILE *fr, poly dstPoly[M][N + 1]) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N + 1; j++) {
            dstPoly[i][j].length = 0;
            dstPoly[i][j].p = nullptr;
        }

        uint32_t temp;
        for (int j = 0; j < N; j++) {
            for (int k = 0; k <= j; k++) {
                fscanf(fr, "%u", &temp);
                if (temp) {
                    ttint128_t term = ttint128_t(1) << j;
                    addTerm(dstPoly[i][k], term);
                }
            }
        }
        for (int j = 0; j < N; j++) {
            fscanf(fr, "%u", &temp);
            if (temp) {
                ttint128_t term = ttint128_t(1) << j;
                addTerm(dstPoly[i][j], term);
            }
        }
        fscanf(fr, "%u", &temp);
        if (temp) {
            ttint128_t term = ttint128_t(1) << N;
            addTerm(dstPoly[i][N], term);
        }
    }
}

inline void ffile2poly(FILE *fr, poly dstPoly[M][N + 1], int m) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < N + 1; j++) {
            fscanf(fr, "%u", &(dstPoly[i][j].length));
            dstPoly[i][j].p = (ttint128_t *) malloc(dstPoly[i][j].length * sizeof(ttint128_t));
            for (int k = 0; k < dstPoly[i][j].length; k++) {
                uint32_t temp[3];
                for (int l = 0; l < 3; l++)
                    fscanf(fr, "%u", &(temp[l]));
                dstPoly[i][j].p[k] =
                    (ttint128_t(temp[2]) << 64) | (ttint128_t(temp[1]) << 32) | (ttint128_t(temp[0]));
            }
        }
    }
}

inline void loadPD(poly fullpoly[M][N + 1],
                   uint32_t colPD[N + 1][RESIDENT_NUM]) {
    // store the partial derivative in row form
    ttint128_t tmpRowPD[M][N];
    for (int i = 0; i < EQUATION_NUM; i++) {
        for (int j = 0; j < 56; j++) {
            for (int k = 0; k < fullpoly[i][j].length; k++) {
                uint64_t los;
                ttint128_t t = fullpoly[i][j].p[k] ^(ttint128_t(1) << j);
                if (t == ttint128_t(0))
                    los = j;
                else
                    los = ttmath_ffsll(t);
                if (los > j) {
                    tmpRowPD[i][j] ^= (ttint128_t(1) << los);
                    tmpRowPD[i][los] ^= (ttint128_t(1) << j);
                } else
                    tmpRowPD[i][j] ^= (ttint128_t(1) << N);
            }
        }
    }

    // transpose row-form into column-form
    for (int valIndex = 0; valIndex < N + 1; valIndex++) {
        for (int bitIndex = 0; bitIndex < RESIDENT_NUM; bitIndex++) {
            for (int eqIndex = 0; eqIndex < EQUATION_NUM; eqIndex++) {
                uint32_t rec;
                ttint128_t prod = tmpRowPD[eqIndex][bitIndex] >> valIndex;
                prod.ToUInt(rec);
                colPD[valIndex][bitIndex] |= (rec & 1) << eqIndex;
            }
        }
    }
}

inline void checkConsist(uint32_t clist[17], uint32_t &mask) {
    uint32_t _mask = 0x3fffff;
    for (int i = 0; i < 16; i++) {
        uint32_t ci = clist[i] & mask;
        if (ci == 0) continue;
        uint32_t x = __builtin_clz(ci);
        uint32_t xp = (0x80000000 >> x);
        uint32_t xormask = (clist[i] & _mask) ^xp;
        for (int j = 0; j < 17; j++)
            clist[j] = (clist[j] & xp) ? (clist[j] ^ xormask) : clist[j];
        mask ^= xp;
    }
}

inline void extractSolution(const uint32_t clist[17], uint32_t sol[16]) {
    for (int i = 0; i < 16; i++) {
        if (clist[i] == 0) continue;
        uint32_t xp = 0x80000000 >> __builtin_clz(clist[i]);
        sol[i] = (bool) (clist[16] & xp);
    }
}

inline bool verifyPoly(const uint32_t guess[N], poly spoly[M][N + 1]) {
    for (int i = 0; i < M; i++) {
        int8_t res = 0;
        for (int j = 0; j < N; j++)
            for (int k = 0; k < spoly[i][j].length; k++)
                res ^= guess[127 - ttmath_lzcntll(spoly[i][j].p[k])] & guess[j];
        if (spoly[i][N].length)
            res ^= 1;
        if (res)
            return false;
    }
    return true;
}

inline void mqInit(uint64_t start, uint64_t end) {
    FILE *ResidentFile = fopen(ResidentFilename.c_str(), "rb");
    ffile2poly(ResidentFile, ResidentPoly, EQUATION_NUM);
    loadPD(ResidentPoly, ColPartialDrt);
    fclose(ResidentFile);
    std::cout << "[INFO] partial derivative loaded from " << ResidentFilename << std::endl;

    FILE *ChallengeFile = fopen(ChallengeFilename.c_str(), "rb");
    file2poly(ChallengeFile, VerifyPoly);
    fclose(ChallengeFile);
    std::cout << "[INFO] challenge loaded from " << ChallengeFilename << std::endl;

    startKey = start;
    endKey = end;
    binKey = start;
    if (binKey != 0)
        grayCodeKey = (binKey - 1) ^ ((binKey - 1) >> 1);
    else
        grayCodeKey = 0;

    for (int eqIndex = 0; eqIndex < EQUATION_NUM; eqIndex++) {
        for (int varIndex = 0; varIndex < RESIDENT_NUM; varIndex++) {
            if ((grayCodeKey >> varIndex) & 1) {
                for (int termIndex = 0; termIndex < ResidentPoly[eqIndex][varIndex].length; termIndex++) {
                    uint32_t los = 127 - ttmath_lzcntll(ResidentPoly[eqIndex][varIndex].p[termIndex]);
                    if (los == varIndex)
                        ColVal[N] ^= 1 << eqIndex;
                    else if (los < RESIDENT_NUM)
                        ColVal[N] ^= ((grayCodeKey >> los) & 1) << eqIndex;
                    else
                        ColVal[los] ^= 1 << eqIndex;
                }
            }
        }
        for (int varIndex = RESIDENT_NUM; varIndex < N; varIndex++)
            if (ResidentPoly[eqIndex][varIndex].length)
                ColVal[varIndex] ^= 1 << eqIndex;
        ColVal[N] ^= ResidentPoly[eqIndex][N].length << eqIndex;
    }
}

inline void mqLoop() {
    uint64_t loopCount = 0;
    while (binKey < endKey) {
        uint32_t Val[VARIABLE_NUM + 1];
        std::copy(ColVal + RESIDENT_NUM, ColVal + N + 1, Val);

        for (int varIndex = 0; varIndex < RESIDENT_NUM; varIndex++)
            Val[VARIABLE_NUM] ^= ColVal[varIndex];

        uint32_t mask = 0x3fffff;
        checkConsist(Val, mask);

        if (!(mask & Val[VARIABLE_NUM])) {
            uint32_t sol[9] = {0};
            extractSolution(Val, sol);
            uint32_t guess[N] = {0};
            for (int varIndex = 0; varIndex < RESIDENT_NUM; varIndex++)
                guess[varIndex] = (uint32_t) ((grayCodeKey >> varIndex) & 1);
            for (int varIndex = 0; varIndex < VARIABLE_NUM; varIndex++)
                guess[RESIDENT_NUM + varIndex] = sol[varIndex];
            if (verifyPoly(guess, VerifyPoly)) {
                std::ofstream file("solution.txt");
                for (int g = 0; g < N; g++)
                    file << guess[g] << std::endl;
                file.close();
            }
        }

        uint64_t z = ((uint64_t) binKey >> 1) ^  binKey ^ grayCodeKey;
        uint32_t los = __builtin_ffsl(z) - 1;
        for (int varIndex = 0; varIndex < RESIDENT_NUM; varIndex++)
            if ((grayCodeKey >> varIndex) & 1)
                ColVal[varIndex] ^= ColPartialDrt[varIndex][los];
        for (int varIndex = RESIDENT_NUM; varIndex < N + 1; varIndex++)
            ColVal[varIndex] ^= ColPartialDrt[varIndex][los];
        grayCodeKey = (binKey ^ ((uint64_t) binKey >> 1));
        binKey++;

        if ((loopCount++ % 0x10000000) == 0)
            std::cout << "rank: " << rank << "\toffset: " << loopCount << std::endl;
    }
}

#endif //PROJECT_MAIN_H

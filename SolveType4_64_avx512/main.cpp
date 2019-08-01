#include <stdint.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <immintrin.h>
#include <emmintrin.h>
#include <bitset>
#include <cstring>

#define M 64
#define N 64
#define LEN 10001
#define VARIABLE_NUM 10
#define EQUATION_NUM 19
#define SEARCH_SPACE 54

struct term {
  uint32_t data[3];
};
struct poly {
  term *p;
  int length;
};

void repeatPoly(poly &spoly);
void ffile2poly(FILE *fr, poly poly[][N + 1], int m);
void addPoly(poly &dstPoly, poly &poly);
void addTerm(poly &dstPoly, term &term);
void loadPD(poly fullpoly[M][N + 1], __int128 partialDerivative[M][N]);
bool verifyPoly(uint32_t guess[N], poly spoly[M][N + 1]);
void file2poly(FILE *fr, poly spoly[M][N + 1]);

void checkConsist_19x10(__m512i clist[11], __m512i &masks);
void extractSolution_19x10(const uint32_t clist[11], uint32_t sol[10]);
const std::string currentDateTime();

int main() {
    FILE *fr = fopen("mq-resident4-64-0.txt", "rb");
    FILE *frr = fopen("mq4-64-0-f.txt", "rb");
    __int128 partialDerivative[EQUATION_NUM][N] = {0};
    poly fullpoly[M][N + 1];
    ffile2poly(fr, fullpoly, EQUATION_NUM);
    loadPD(fullpoly, partialDerivative);
    poly verifypoly[M][N + 1];
    ffile2poly(frr, verifypoly, M);
    fclose(fr);
    fclose(frr);

    __int128 value[16][EQUATION_NUM] = {0};
    int64_t z;
    __int128 pre[16] = {(__int128(0x7FF) << SEARCH_SPACE)};
    int los = 0;
    __int128 guessMask = 0x3fffffffffffff;

    int64_t keys[16];
    int64_t init_key = 0;
    int64_t final_key;
    int64_t seg_size = (final_key - init_key) / 16;
    for (int i = 0; i < 16; i++)
        keys[i] = init_key + seg_size * i;
    for (int i = 0; i < 16; i++) {
        if (keys[i] != 0)
            pre[i] = ((__int128) (0x7FF) << SEARCH_SPACE)
                | ((__int128) (keys[i] - 1) ^ (__int128) (keys[i] - 1) >> 1);
        for (int j = 0; j < EQUATION_NUM; j++) {
            for (int k = 0; k < SEARCH_SPACE; k++) {
                if (((pre[i] >> k) & 1) == 1) {
                    for (int l = 0; l < fullpoly[j][k].length; l++) {
                        int los = k;
                        for (int ll = k + 1; ll < N; ll++) {
                            if ((fullpoly[j][k].p[l].data[ll / 32] >> (ll % 32)) & 1) {
                                los = ll;
                                break;
                            }
                        }
                        if (los == k)
                            value[i][j] ^= ((__int128) 1 << N);
                        else if (los < SEARCH_SPACE)
                            value[i][j] ^= (((pre[i] >> los) & 1) << N);
                        else
                            value[i][j] ^= (((pre[i] >> los) & 1) << los);
                    }
                }
            }
            for (int k = SEARCH_SPACE; k < N; k++)
                if (fullpoly[j][k].length)
                    value[i][j] ^= (__int128) 1 << j;
            value[i][j] ^= (__int128) fullpoly[j][N].length << N;
        }
    }

    uint32_t matrix[16][EQUATION_NUM] = {0};
    for (int64_t key = init_key; key < init_key + seg_size; key++) {
        // set matrices
        uint32_t clist[16][11] = {0};
        for (int i = 0; i < 16; i++) {
            memset(matrix[i], 0, EQUATION_NUM * 2);
            for (int j = 0; j < EQUATION_NUM; j++) {
                matrix[i][j] = __builtin_parityl(value[i][j] & guessMask) ^ (value[i][j] >> N);
                matrix[i][j] |= (value[i][j] >> 53) & 0x7FE;
            }
            for (int j = 0; j < VARIABLE_NUM + 1; j++)
                for (int k = 0; k < EQUATION_NUM; k++)
                    clist[i][j] |= ((matrix[i][k] >> (10 - j)) & 1) << k;
        }

        // gaussian elimination
        __m512i m512clist[11];
        for (int j = 0; j < 11; j++) {
            m512clist[j] = _mm512_set_epi32(clist[0][j],
                                            clist[1][j],
                                            clist[2][j],
                                            clist[3][j],
                                            clist[4][j],
                                            clist[5][j],
                                            clist[6][j],
                                            clist[7][j],
                                            clist[8][j],
                                            clist[9][j],
                                            clist[10][j],
                                            clist[11][j],
                                            clist[12][j],
                                            clist[13][j],
                                            clist[14][j],
                                            clist[15][j]);
        }
        __m512i masks = _mm512_set1_epi32(0X7FFFF);
        checkConsist_19x10(m512clist, masks);

        // verify result
        for (int i = 0; i < 16; i++) {
            if (!(masks[i] & m512clist[10][i])) {
                uint32_t sol[10] = {0};
                uint32_t ci[11] = {0};
                for (int j = 0; j < 11; j++)
                    ci[j] = m512clist[i][j];
                extractSolution_19x10(ci, sol);
                uint32_t guess[N] = {0};
                for (int j = 0; j < SEARCH_SPACE; j++)
                    guess[j] = (uint32_t) (pre[i] >> j) & 1;
                for (int j = 0; j < 10; j++)
                    guess[SEARCH_SPACE + j] = sol[VARIABLE_NUM - j - 1];
                if (verifyPoly(guess, verifypoly)) {
                    std::ofstream file("solution.txt");
                    for (auto &g: guess)
                        file << g << std::endl;
                    file.close();
                }
            }
        }

        // update value
        for (int i = 0; i < 16; i++) {
            z = ((uint64_t) keys[i] >> 1) ^ keys[i] ^ pre[i];
            los = __builtin_ffsl(z) - 1;
            for (int j = 0; j < EQUATION_NUM; j++)
                value[i][j] ^= (partialDerivative[j][los] & pre[i]);
            pre[i] = ((__int128) key ^ ((__int128) key >> 1)) | ((__int128) 0x7FF << SEARCH_SPACE);
        }
    }

    return 0;
}

const std::string currentDateTime() {
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}
bool verifyPoly(uint32_t guess[N], poly spoly[M][N + 1]) {
    for (int i = 0; i < M; i++) {
        int8_t res = 0;
        for (int j = 0; j < N; j++) {
            if (guess[j]) {
                for (int k = 0; k < spoly[i][j].length; k++) {
                    int8_t r = 1;
                    for (int l = j + 1; l < N; l++) {
                        int word = l / 32;
                        int bit = l % 32;
                        if ((spoly[i][j].p[k].data[word] >> bit) & 1)
                            r &= guess[l];
                    }
                    res ^= r;
                }
            }
        }

        if (spoly[i][N].length)
            res ^= 1;
        if (res)
            return false;
    }
    return true;
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
void file2poly(FILE *fr, poly spoly[M][N + 1]) {
    uint32_t temp;
    term tterm;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N + 1; j++) {
            spoly[i][j].length = 0;
            spoly[i][j].p = NULL;
        }
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k <= j; k++) {
                fscanf(fr, "%u", &temp);
                if (temp) {
                    memset(tterm.data, 0, 12);
                    int word = j / 32;
                    int bit = j % 32;
                    tterm.data[word] |= (1 << bit);
                    word = k / 32;
                    bit = k % 32;
                    tterm.data[word] |= (1 << bit);
                    addTerm(spoly[i][k], tterm);
                }
            }
        }
        for (int j = 0; j < N; j++) {
            fscanf(fr, "%u", &temp);
            if (temp) {
                memset(tterm.data, 0, 12);
                int word = j / 32;
                int bit = j % 32;
                tterm.data[word] |= (1 << bit);
                addTerm(spoly[i][j], tterm);
            }
        }
        fscanf(fr, "%u", &temp);
        if (temp) {
            memset(tterm.data, 0, 12);
            int word = N / 32;
            int bit = N % 32;
            tterm.data[word] |= (1 << bit);
            addTerm(spoly[i][N], tterm);
        }
    }
}
void loadPD(poly fullpoly[M][N + 1], __int128 partialDerivative[M][N]) {
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
                    partialDerivative[i][j] ^= (__int128(1) << los);
                    partialDerivative[i][los] ^= (__int128(1) << j);
                } else
                    partialDerivative[i][j] ^= (__int128(1) << N);
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
void checkConsist_19x10(__m512i clist[11], __m512i &masks) {
    for (int i = 0; i < 10; i++) {
        __m512i ci = _mm512_and_epi32(clist[i], masks);

    }
}
void extractSolution_19x10(const uint32_t clist[11], uint32_t sol[10]) {
    for (int i = 0; i < 10; i++) {
        if (clist[i] == 0) continue;
        uint32_t xp = (0x80000000 >> (__builtin_clz(clist[i])));
        sol[i] = (bool) (clist[10] & xp);
    }
}

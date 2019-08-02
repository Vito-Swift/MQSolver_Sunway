#include <stdint.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <bitset>
#include <cstring>

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

#define M 55
#define N 55
#define LEN 10001
#define VARIABLE_NUM 9
#define EQUATION_NUM 19
#define SEARCH_SPACE 46

// mq arithmetic
void repeatPoly(poly &spoly);

void ffile2poly(FILE *fr, poly poly[][N + 1], int m);

void addPoly(poly &dstPoly, poly &poly);

void addTerm(poly &dstPoly, term &term);

void loadPD(poly fullpoly[M][N + 1], int64_t partialDerivative[M][N]);

bool verifyPoly(uint32_t guess[N], poly spoly[M][N + 1]);

void file2poly(FILE *fr, poly spoly[M][N + 1]);

// gf2 arithmetic
void checkConsist_19x9(uint32_t clist[10], uint32_t &mask);

void extractSolution_19x9(const uint32_t clist[10], uint32_t sol[9]);

const std::string currentDateTime();

int main() {
    // MPI Init
    FILE *fr = fopen("mq-resident4-55-0.txt", "rb");
    FILE *frr = fopen("mq4-55-0-f.txt", "rb");
    poly fullpoly[M][N + 1];
    int64_t partialDerivative[EQUATION_NUM][N] = {0};
    ffile2poly(fr, fullpoly, EQUATION_NUM);
    loadPD(fullpoly, partialDerivative);
    fclose(fr);

    poly verifypoly[M][N + 1];
    ffile2poly(frr, verifypoly, M);
    fclose(frr);

    uint64_t value[EQUATION_NUM] = {0};
    uint64_t z;
    uint64_t pre = ((uint64_t) 0x3FF << SEARCH_SPACE);
    int los = 0;
    uint64_t guessMask = 0x3fffffffffff;

    uint64_t init_key = 0x1812fc4e0000;
    uint64_t end_key = 0x1812fc4e4950;

    // Equation Init
    if (init_key != 0)
        pre = ((uint64_t) 0x3FF << SEARCH_SPACE) | (((init_key - 1)) ^ (((init_key - 1)) >> 1));
    std::cout << pre << std::endl;
    for (int i = 0; i < EQUATION_NUM; i++) {
        for (int j = 0; j < SEARCH_SPACE; j++) {
            if (((pre >> j) & 1) == 1) {
                for (int k = 0; k < fullpoly[i][j].length; k++) {
                    int los = j;
                    for (int ll = j + 1; ll < N; ll++) {
                        if ((fullpoly[i][j].p[k].data[ll / 32] >> (ll % 32)) & 1) {
                            los = ll;
                            break;
                        }
                    }
                    if (los == j)
                        value[i] ^= ((uint64_t) 1 << N);
                    else if (los < SEARCH_SPACE)
                        value[i] ^= (((pre >> los) & 1) << N);
                    else
                        value[i] ^= (((pre >> los) & 1) << los);
                }
            }

        }
        for (int j = SEARCH_SPACE; j < N; j++) {
            if (fullpoly[i][j].length)
                value[i] ^= ((uint64_t) 1 << j);
        }
        value[i] ^= ((uint64_t) fullpoly[i][N].length << N);
    }
    std::cout << init_key << std::endl;
    uint32_t matrix[EQUATION_NUM] = {0};
    for (int64_t key = init_key; key < end_key; key++) { // Exhaustive Search
        // Step 1: set matrix
        memset(matrix, 0, EQUATION_NUM * 4);
        uint32_t clist[10] = {0};
        for (int j = 0; j < EQUATION_NUM; j++) {
            matrix[j] = __builtin_parityl(value[j] & guessMask) ^ (value[j] >> N);
            matrix[j] |= (value[j] >> 45) & 0x3FE;
        }
        for (int i = 0; i < VARIABLE_NUM + 1; i++)
            for (int j = 0; j < EQUATION_NUM; j++)
                clist[i] |= ((matrix[j] >> (9 - i)) & 1) << j;

        // Step 2: Gaussian Elimination
        uint32_t mask = 0x7FFFF;
        checkConsist_19x9(clist, mask);
        if (!(mask & clist[9])) {
            for (auto& c: clist)
                std::cout << std::bitset<19> (c) << std::endl;
            std::cout << std::endl;
            uint32_t sol[9] = {0};
            extractSolution_19x9(clist, sol);
            uint32_t guess[N] = {0};
            for (int i = 0; i < SEARCH_SPACE; i++)
                guess[i] = (pre >> i) & 1;
            for (int i = 0; i < 9; i++)
                guess[SEARCH_SPACE + i] = sol[VARIABLE_NUM - i - 1];
            if (verifyPoly(guess, verifypoly)) {
                std::ofstream file("solution.txt");
                for (auto &g: guess)
                    file << g << std::endl;
                file.close();
            }
        }

        // Step 3: Update Polynomial
        z = ((uint64_t) key >> 1) ^ key ^ pre;
        los = __builtin_ffsl(z) - 1;
        for (int j = 0; j < EQUATION_NUM; j++)
            value[j] ^= partialDerivative[j][los] & pre;
        pre = (key ^ ((uint64_t) key >> 1)) | ((uint64_t) 0x3FF << SEARCH_SPACE);
    }
    return EXIT_SUCCESS;
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

void loadPD(poly fullpoly[M][N + 1], int64_t partialDerivative[M][N]) {
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

void checkConsist_19x9(uint32_t clist[10], uint32_t &mask) {
    uint32_t _mask = 0x7FFFF;
    for (int i = 0; i < 9; i++) {
        uint32_t ci = clist[i] & mask;
        if (ci == 0) continue;
        uint32_t x = __builtin_clz(ci);
        uint32_t xp = (0x80000000 >> x);
        uint32_t xormask = (clist[i] & _mask) ^xp;
        for (int j = 0; j < 10; j++)
            clist[j] = (clist[j] & xp) ? (clist[j] ^ xormask) : clist[j];
        mask ^= xp;
    }
}

void extractSolution_19x9(const uint32_t clist[10], uint32_t sol[9]) {
    for (int i = 0; i < 9; i++) {
        if (clist[i] == 0) continue;
        uint32_t xp = (0x80000000 >> (__builtin_clz(clist[i])));
        sol[i] = (bool) (clist[9] & xp);
    }
}

/**
 * @filename mqfuncs.cpp
 * @headerfile mqfuncs.h
 * @author Chenhao Wu
 *         (email: chenhaowu[at]link.cuhk.edu.cn)
 *
 *   Various methods for MQ polynomial problem
 *
 *   See COPYRIGHT file at the top of the source tree.
 *
 *   This program is a free software: you can redistribute
 * it and/or modify it under the terms of GNU General Public
 * License as published by the Free Software Foundation.
 */

#include "mqmisc.h"

using namespace std;

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

void poly2file(FILE *fw, poly spoly) {
    fprintf(fw, "%u ", spoly.length);
    for (int k = 0; k < spoly.length; k++) {
        for (int l = 0; l < 3; l++)
            fprintf(fw, "%u ", spoly.p[k].data[l]);
    }
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

void mulPoly(poly &poly3, poly &poly1, poly &poly2) {
    poly3.length = 0;
    if ((poly3.p = (term *) malloc(poly1.length * poly2.length * sizeof(term))) == NULL) {
        printf("Out of Memory!\n");
        exit(1);
    }
    for (int i = 0; i < poly1.length; i++) {
        for (int j = 0; j < poly2.length; j++) {
            for (int l = 0; l < 3; l++) {
                poly3.p[poly3.length].data[l] = (poly1.p[i].data[l]) | (poly2.p[j].data[l]);
            }
            int word = N / 32;
            int bit = N % 32;
            uint32_t tt = (((poly1.p[i].data[word]) & (poly2.p[j].data[word]) & (1 << bit)) ^ (1 << bit) ^ 0xffffffff);
            poly3.p[poly3.length].data[word] &= tt;
            (poly3.length)++;
        }
    }
}

void printPoly(FILE *fw, poly spoly) {
    for (int i = 0; i < spoly.length; i++) {
        fprintf(fw, "+");
        for (int j = 0; j < N + 1; j++) {
            int word = j / 32;
            int bit = j % 32;
            if ((spoly.p[i].data[word] >> bit) & 1)
                fprintf(fw, "x%d", j);
        }
    }
    fprintf(fw, "\n");
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

void simplifyPoly(poly spoly[M][N + 1]) {
    for (int i = 0; i < M; i++) {
        for (int j = 1; j < N; j++) {
            int len = 0;
            for (int k = 0; k < spoly[i][j].length; k++) {
                int tag = 1;
                for (int los = 0; los < j; los++) {
                    int word = los / 32;
                    int bit = los % 32;
                    if ((spoly[i][j].p[k].data[word] >> bit) & 1) {
                        addTerm(spoly[i][los], spoly[i][j].p[k]);
                        tag = 0;
                        break;
                    }
                }
                if (tag) {
                    for (int l = 0; l < 3; l++)
                        spoly[i][j].p[len].data[l] = spoly[i][j].p[k].data[l];
                    len++;
                }
            }
            spoly[i][j].length = len;
        }
        int len = 0;
        for (int k = 0; k < spoly[i][N].length; k++) {
            int tag = 1;
            for (int los = 0; los < N; los++) {
                int word = los / 32;
                int bit = los % 32;
                if ((spoly[i][N].p[k].data[word] >> bit) & 1) {
                    addTerm(spoly[i][los], spoly[i][N].p[k]);
                    tag = 0;
                    break;
                }
            }
            if (tag) {
                for (int l = 0; l < 3; l++)
                    spoly[i][N].p[len].data[l] = spoly[i][N].p[k].data[l];
                len++;
            }
        }
        spoly[i][N].length = len;
    }
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N + 1; j++) {
            repeatPoly(spoly[i][j]);
        }
        if (spoly[i][N].length == 1) {
            spoly[i][N].p[0].data[N / 32] |= (1 << (N % 32));
        }
    }
}

bool verifyPoly(int64_t guess[N], poly spoly[M][N + 1]) {
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

int getLos(int64_t x) {
    switch (x) {
        case 0x1: return 0;
        case 0x2: return 1;
        case 0x4: return 2;
        case 0x8: return 3;
        case 0x10: return 4;
        case 0x20: return 5;
        case 0x40: return 6;
        case 0x80: return 7;
        case 0x100: return 8;
        case 0x200: return 9;
        case 0x400: return 10;
        case 0x800: return 11;
        case 0x1000: return 12;
        case 0x2000: return 13;
        case 0x4000: return 14;
        case 0x8000: return 15;
        case 0x10000: return 16;
        case 0x20000: return 17;
        case 0x40000: return 18;
        case 0x80000: return 19;
        case 0x100000: return 20;
        case 0x200000: return 21;
        case 0x400000: return 22;
        case 0x800000: return 23;
        case 0x1000000: return 24;
        case 0x2000000: return 25;
        case 0x4000000: return 26;
        case 0x8000000: return 27;
        case 0x10000000: return 28;
        case 0x20000000: return 29;
        case 0x40000000: return 30;
        case 0x80000000: return 31;
        case 0x100000000: return 32;
        case 0x200000000: return 33;
        case 0x400000000: return 34;
        case 0x800000000: return 35;
        default: return 36;
    }
}
#include "gf2misc.h"

void checkConsist_28x16(uint32_t clist[17], uint32_t& mask) {
    for (int64_t i = 0; i < 16; i++) {
        uint32_t ci = clist[i] & mask;
        if (ci == 0) continue;
        uint32_t x = __lzcnt32(ci);
        uint32_t xp = (0x80000000 >> x);
        uint32_t xormask = ci ^ xp;
        for (int64_t j = 0; j < 17; j++)
            clist[j] = (clist[j] & xp) ? (clist[j] ^ xormask) : clist[j];
        mask ^= xp;
    }
}

void checkConsist_22x10(uint32_t clist[11], uint32_t& mask) {
    for (int64_t i = 0; i < 10; i++) {
        uint32_t ci = clist[i] & mask;
        if (ci == 0) continue;
        uint32_t x = __lzcnt32(ci);
        uint32_t xp = (0x80000000 >> x);
        uint32_t xormask = ci ^ xp;
        for (int64_t j = 0; j < 11; j++)
            clist[j] = (clist[j] & xp) ? (clist[j] ^ xormask) : clist[j];
        mask ^= xp;
    }
}

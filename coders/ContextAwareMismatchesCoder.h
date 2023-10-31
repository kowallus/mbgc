#ifndef MBGC_CONTEXTAWAREMISMATCHESCODER_H
#define MBGC_CONTEXTAWAREMISMATCHESCODER_H

#include "../utils/helper.h"

class ContextAwareMismatchesCoder {
private:

    static uint8_t sym2val[256];
    static constexpr char val2sym[128] = "ACGTNacgtnUYRKMBDHVSWuyrkmbdhvsw";
    static constexpr int val2symLen = 32; // length or val2sym

    static constexpr char mis2code[5][5] =  { { 'A', 2, 0, 1, 3 },
                             { 1, 'C', 2, 0, 3 },
                             { 0, 2, 'G', 1, 3 },
                             { 1, 0, 2, 'T', 3 },
                             { 1, 2, 3, 0, 'N' }};
    static char code2mis[5][5];

    static struct StaticInitializer { StaticInitializer(); } initializer;

    inline uint8_t symbol2value(char symbol) {
        return sym2val[symbol];
    };

    inline char value2symbol(uint8_t value) {
        return val2sym[value];
    }

public:

    ContextAwareMismatchesCoder();

    uint8_t mismatch2code(char actual, char mismatch, char context = 'N');
    char code2mismatch(char actual, uint8_t code);
    char encodeMismatchWithExclusion(char actual, char mismatch);
    char decodeMismatchWithExclusion(char actual, char mismatch);

    static ContextAwareMismatchesCoder defaultInstance;

};

#endif //MBGC_CONTEXTAWAREMISMATCHESCODER_H

#include "ContextAwareMismatchesCoder.h"

uint8_t ContextAwareMismatchesCoder::sym2val[256];
char ContextAwareMismatchesCoder::code2mis[5][5];
ContextAwareMismatchesCoder ContextAwareMismatchesCoder::defaultInstance;
ContextAwareMismatchesCoder::StaticInitializer ContextAwareMismatchesCoder::initializer;

ContextAwareMismatchesCoder::StaticInitializer::StaticInitializer() {
    memset(sym2val, UINT8_MAX, 256);
    for(int i = 0; i < val2symLen; i++)
        sym2val[val2sym[i]] = i;
    for(char i = 0; i < 5; i++)
        for(char j = 0; j < 5; j++) {
            if (mis2code[i][j] < 5)
                code2mis[i][mis2code[i][j]] = val2sym[j];
        }
}

ContextAwareMismatchesCoder::ContextAwareMismatchesCoder() {
}

uint8_t ContextAwareMismatchesCoder::mismatch2code(char actual, char mismatch, char context) {
    if (actual == mismatch || sym2val[actual] > 4 || sym2val[mismatch] > 4)
        return (uint8_t) mismatch;
    uint8_t mismatchCode = mis2code[sym2val[actual]][sym2val[mismatch]];
    return mismatchCode;
}

char ContextAwareMismatchesCoder::code2mismatch(char actual, uint8_t code) {
    if (code >= val2symLen)
        return (char) code;
    char mismatch = code2mis[sym2val[actual]][code];
    return mismatch;
}

char ContextAwareMismatchesCoder::encodeMismatchWithExclusion(char actual, char mismatch) {
    return value2symbol(mismatch2code(actual, mismatch));
}

char ContextAwareMismatchesCoder::decodeMismatchWithExclusion(char actual, char encoded) {
    return code2mismatch(actual, symbol2value(encoded));
}

#ifndef PGTOOLS_PPMDCODER_H
#define PGTOOLS_PPMDCODER_H

#include "LzmaCoder.h"

class PpmdCoderProps: public CoderProps {
public:
    uint32_t memSize;
    const uint8_t order;

    PpmdCoderProps(uint32_t memSize, uint8_t order)
            : CoderProps(PPMD7_CODER, 1), memSize(memSize), order(order) { }

    string log() {
        return " ppmd (mem = " + to_string(memSize >> 20) + "MB; ord = " + to_string(order) + ")";
    }

    virtual ~PpmdCoderProps() { };
};

MY_STDAPI Ppmd7Compress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                        PpmdCoderProps* coder_props, double estimated_compression);

MY_STDAPI PpmdUncompress(unsigned char *dest, size_t *destLen, istream &src, size_t *srcLen,
        ostream* logout = PgHelpers::logout);

#endif //PGTOOLS_PPMDCODER_H

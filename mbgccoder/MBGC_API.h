#ifndef MBGC_MBGC_API_H
#define MBGC_MBGC_API_H

#include <iostream>
#include <stdint.h>

using namespace std;

class MBGC_Decoder_API {

public:
    MBGC_Decoder_API() = default;
    virtual ~MBGC_Decoder_API() = default;

    virtual void decode() = 0;

    virtual bool iterateNext(const char* name, uint8_t* &res, size_t &size) = 0;
};


#endif //MBGC_MBGC_API_H

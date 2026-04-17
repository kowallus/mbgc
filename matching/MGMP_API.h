#ifndef MGMP_MBGC_API_H
#define MGMP_MBGC_API_H

#include <iostream>
#include <stdint.h>

using namespace std;

class MGMP_Decoder_API {

public:
    MGMP_Decoder_API() = default;
    virtual ~MGMP_Decoder_API() = default;

    virtual void decode() = 0;

    virtual bool iterateNext(const char* name, uint8_t* &res, size_t &size) = 0;
};


#endif //MGMP_MBGC_API_H

#ifndef MGMP_INPUT_WITH_LIBDEFLATE_WRAPPER_H
#define MGMP_INPUT_WITH_LIBDEFLATE_WRAPPER_H

#include "../utils/helper.h"

#include "../coders/libdeflate/libdeflate.h"
#include "MGMP_API.h"

struct mgmpInFile {
    uint8_t* fileOut;
    size_t fileSize;

    // current view
    uint8_t* out;
    size_t size;
    size_t pos;
};

#define GZIP_ID1		0x1F
#define GZIP_ID2		0x8B

extern MGMP_Decoder_API* iterableDecoder;

static uint32_t
load_u32_gzip(const uint8_t *p)
{
    return ((uint32_t)p[0] << 0) | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

mgmpInFile mgmpInOpen(const char *filename);

size_t mgmpInRead(mgmpInFile &gzf, void* buf, size_t nbyte);

mgmpInFile& mgmpInReset(mgmpInFile &gzf);

int mgmpInClose(mgmpInFile &gzf);

string gzcompress(string& src, int clevel);

// for inner-file iteration purposes
mgmpInFile& mgmpInSplit_init(mgmpInFile &gzf);
mgmpInFile mgmpInSplit_next(mgmpInFile &itergzf, size_t minSplitSize, char splitChar);

#endif //MGMP_INPUT_WITH_LIBDEFLATE_WRAPPER_H

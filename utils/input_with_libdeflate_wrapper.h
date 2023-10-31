#ifndef MBGC_INPUT_WITH_LIBDEFLATE_WRAPPER_H
#define MBGC_INPUT_WITH_LIBDEFLATE_WRAPPER_H

#include "helper.h"

#include "../coders/libdeflate/libdeflate.h"
#include "../mbgccoder/MBGC_API.h"

struct mbgcInFile {
    uint8_t* fileOut;
    size_t fileSize;

    // current view
    uint8_t* out;
    size_t size;
    size_t pos;
};

#define GZIP_ID1		0x1F
#define GZIP_ID2		0x8B

extern MBGC_Decoder_API* iterableDecoder;

static uint32_t
load_u32_gzip(const uint8_t *p)
{
    return ((uint32_t)p[0] << 0) | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

mbgcInFile mbgcInOpen(const char *filename);

size_t mbgcInRead(mbgcInFile &gzf, void* buf, size_t nbyte);

mbgcInFile& mbgcInReset(mbgcInFile &gzf);

int mbgcInClose(mbgcInFile &gzf);

string gzcompress(string& src, int clevel);

// for inner-file iteration purposes
mbgcInFile& mbgcInSplit_iter(mbgcInFile &gzf, size_t minSplitSize, char splitChar);
mbgcInFile mbgcInSplit_next(mbgcInFile &itergzf, size_t minSplitSize, char splitChar);

#endif //MBGC_INPUT_WITH_LIBDEFLATE_WRAPPER_H

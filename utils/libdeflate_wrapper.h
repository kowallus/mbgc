#ifndef MBGC_LIBDEFLATE_WRAPPER_H
#define MBGC_LIBDEFLATE_WRAPPER_H

#include "helper.h"
#include <libdeflate.h>

struct gzFile {
    uint8_t* out;
    size_t size;
    size_t pos;
};

#define GZIP_ID1		0x1F
#define GZIP_ID2		0x8B

static uint32_t
load_u32_gzip(const uint8_t *p)
{
    return ((uint32_t)p[0] << 0) | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

gzFile gzopen(const char *filename, const char *mode);

size_t gzread(gzFile &gzf, void* buf, size_t nbyte);

int gzclose(gzFile &gzf);


#endif //MBGC_LIBDEFLATE_WRAPPER_H

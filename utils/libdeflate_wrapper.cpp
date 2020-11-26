#include "libdeflate_wrapper.h"

gzFile gzopen(const char *filename, const char *mode) {
    gzFile gzf;
    struct libdeflate_decompressor* d;
    d = libdeflate_alloc_decompressor();
    if (d == NULL) {
        fprintf(stderr, "Cannot allocate decompressor.\n");
        exit(EXIT_FAILURE);
    }

    FILE* f = fopen(filename, mode);
    if (f == NULL) {
        string tmpname = filename + string(".gz");
        f = fopen(tmpname.c_str(), mode);
        if (f == NULL) {
            fprintf(stderr, "Cannot open file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
    }
    fseek(f, 0, SEEK_END);
    size_t size = ftell(f);
    fseek(f, 0, SEEK_SET);
    uint8_t* in = (uint8_t*) malloc(size);
    size_t read_bytes = fread(in, 1, size, f);
    if (read_bytes != size) {
        fprintf(stderr, "Problem reading from file: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    fclose(f);

    if (in[0] == GZIP_ID1 && in[1] == GZIP_ID2) {
        size_t uncompressed_size = load_u32_gzip(&in[size - 4]);
        if (uncompressed_size == 0)
            uncompressed_size = 1;
        gzf.out = (uint8_t*) malloc(uncompressed_size);
        size_t in_size;
        enum libdeflate_result res = libdeflate_gzip_decompress_ex(d, in, size, gzf.out, uncompressed_size, &in_size, &gzf.size);
        if (res != LIBDEFLATE_SUCCESS) {
            fprintf(stderr, "Error decompressing gz file: %d.\n", res);
            exit(EXIT_FAILURE);
        }
        free(in);
    } else { //input is not compressed?
        gzf.out = in;
        gzf.size = size;
    }

    libdeflate_free_decompressor(d);
    gzf.pos = 0;
    return gzf;
}

size_t gzread(gzFile &gzf, void *buf, size_t nbyte) {
    if (gzf.size - gzf.pos < nbyte)
        nbyte = gzf.size - gzf.pos;
    memcpy(buf, gzf.out + gzf.pos, nbyte);
    gzf.pos += nbyte;
    return nbyte;
}

int gzclose(gzFile &gzf) {
    free(gzf.out);
    return 0;
}

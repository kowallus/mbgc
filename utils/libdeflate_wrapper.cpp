#include "libdeflate_wrapper.h"
#include "../mbgccoder/MBGC_Params.h"

gzFile gzopen(const char *filename) {
    gzFile gzf;
    uint8_t* in;
    size_t size;

    if (strcmp(filename, MBGC_Params::STANDARD_IO_POSIX_ALIAS) == 0) {
        int c;
        size_t capacity = 4096, i = 0;
        void *newPtr = NULL;
        in = (uint8_t*) malloc(capacity * sizeof (char));
        while (in != NULL && (c = getchar()) != EOF)
        {
            if (i == capacity * sizeof (char))
            {
                capacity *= 2;
                if ((in = (uint8_t*) realloc(in, capacity * sizeof (char))) == NULL)
                    break;
            }
            in[i++] = c;
        }
        if (in == NULL) {
            fprintf(stderr, "Cannot read from stdin\n");
            exit(EXIT_FAILURE);
        }
        size = i;
    } else {
        FILE *f = fopen(filename, "rb");
        if (f == NULL) {
#ifndef NO_GZ_SUPPORT
            string tmpname = filename + string(".gz");
            f = fopen(tmpname.c_str(), "rb");
            if (f == NULL)
#endif
            {
                fprintf(stderr, "Cannot open file: %s\n", filename);
                exit(EXIT_FAILURE);
            }
        }
        fseek(f, 0, SEEK_END);
        size = ftell(f);
        fseek(f, 0, SEEK_SET);
        in = (uint8_t *) malloc(size);
        size_t read_bytes = fread(in, 1, size, f);
        if (read_bytes != size) {
            fprintf(stderr, "Problem reading from file: %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fclose(f);
    }
#ifndef NO_GZ_SUPPORT
    if (in[0] == GZIP_ID1 && in[1] == GZIP_ID2) {
        size_t uncompressed_size = load_u32_gzip(&in[size - 4]);
        if (uncompressed_size == 0)
            uncompressed_size = 1;
        gzf.out = (uint8_t*) malloc(uncompressed_size);
        size_t in_size;
        struct libdeflate_decompressor* d;
        d = libdeflate_alloc_decompressor();
        if (d == NULL) {
            fprintf(stderr, "Cannot allocate decompressor.\n");
            exit(EXIT_FAILURE);
        }
        enum libdeflate_result res = libdeflate_gzip_decompress_ex(d, in, size, gzf.out, uncompressed_size, &in_size, &gzf.size);
        if (res != LIBDEFLATE_SUCCESS) {
            fprintf(stderr, "Error decompressing gz file: %d.\n", res);
            exit(EXIT_FAILURE);
        }
        free(in);
        libdeflate_free_decompressor(d);
    } else
#endif
    { //input is not compressed?
        gzf.out = in;
        gzf.size = size;
    }
    gzf.pos = 0;
    gzf.fileOut = gzf.out;
    gzf.fileSize = gzf.size;
    return gzf;
}

size_t gzread(gzFile &gzf, void *buf, size_t nbyte) {
    if (gzf.size - gzf.pos < nbyte)
        nbyte = gzf.size - gzf.pos;
    memcpy(buf, gzf.out + gzf.pos, nbyte);
    gzf.pos += nbyte;
    return nbyte;
}

gzFile& gzreset(gzFile &gzf) {
    gzf.pos = 0;
    return gzf;
};


int gzclose(gzFile &gzf) {
    free(gzf.fileOut);
    return 0;
}

gzFile &gzsplit_iter(gzFile &gzf, size_t minSplitSize, char splitChar) {
    gzf.fileOut = gzf.out;
    gzf.fileSize = gzf.size;
    gzf.size = 0;
    gzf.pos = 0;
    gzsplit_next(gzf, minSplitSize, splitChar);
    return gzf;
}

gzFile gzsplit_next(gzFile &itergzf, size_t minSplitSize, char splitChar) {
    gzFile gzf = itergzf;
    gzf.out += gzf.size;
    gzf.pos = 0;
    size_t pos = gzf.out - gzf.fileOut;
    if (gzf.out - gzf.fileOut == gzf.fileSize)
        gzf.size = 0;
    else {
        size_t endPos = pos + minSplitSize;
        if (endPos >= gzf.fileSize)
            gzf.size = gzf.fileSize - pos;
        else {
            void *split = memchr(gzf.out + minSplitSize, splitChar, gzf.fileSize - endPos);
            if (split == nullptr)
                gzf.size = gzf.fileSize - pos;
            else
                gzf.size = (uint8_t*) split - gzf.out;
        }
    }
    itergzf = gzf;
    return gzf;
}


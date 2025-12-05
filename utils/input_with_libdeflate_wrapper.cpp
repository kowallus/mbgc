#include "input_with_libdeflate_wrapper.h"
#include "../mbgccoder/MBGC_Params.h"

MBGC_Decoder_API* iterableDecoder = nullptr;

mbgcInFile mbgcInOpen(const char *filename) {
    mbgcInFile gzf;
    uint8_t* in;
    size_t size;

    if (iterableDecoder != nullptr) {
        iterableDecoder->iterateNext(filename, in, size);
    } else if (strcmp(filename, MBGC_Params::STANDARD_IO_POSIX_ALIAS) == 0) {
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
            string tmpname = filename + string(".gz");
            f = fopen(tmpname.c_str(), "rb");
            if (f == NULL) {
                fprintf(stderr, "ERROR: Cannot open file: %s\n", filename);
                switch (filename[0]) {
                case ' ': case '\t': case '\r': case '\n':
                    fprintf(stderr, "Invalid sequences list file? File name cannot start from '%c' whitespace.\n",
                        filename[0]);
                    exit(EXIT_FAILURE);
                default: ;
                }
                for (int i = 0; i < strlen(filename); i++) {
                    switch (filename[i]) {
                    case '<': case '>': case ':': case '"': case '|': case '?': case '*':
                        fprintf(stderr, "Invalid sequences list file? '%s' is not a valid file name.\n", filename);
                        exit(EXIT_FAILURE);
                    default: ;
                    }
                }
                exit(EXIT_FAILURE);
            }
        }
#ifdef __MINGW32__
        _fseeki64(f, 0, SEEK_END);
        size = _ftelli64(f);
        _fseeki64(f, 0, SEEK_SET);
#else
        fseeko(f, 0, SEEK_END);
        size = ftello(f);
        fseeko(f, 0, SEEK_SET);
#endif
        in = (uint8_t *) malloc(size);
        size_t read_bytes = fread(in, 1, size, f);
        if (read_bytes != size) {
            fprintf(stderr, "Problem reading from file: %s (read_bytes: %zu < size: %zu; ferror code: %d; feof code: %d)\n", filename, read_bytes, size, ferror( f ), feof(f));
            exit(EXIT_FAILURE);
        }
        fclose(f);
    }
    if (in[0] == GZIP_ID1 && in[1] == GZIP_ID2) {
        size_t uncompressed_size = load_u32_gzip(&in[size - 4]);
        if (uncompressed_size == 0)
            uncompressed_size = size * 4;
        if ((gzf.out = (uint8_t*) malloc(uncompressed_size)) == NULL) {
            fprintf(stderr, "Error: out of memory during %s decompression.\n", filename);
            exit(EXIT_FAILURE);
        };
        gzf.size = 0;
        size_t inPos = 0;
        size_t in_size, ret_size;
        struct libdeflate_decompressor* d;
        d = libdeflate_alloc_decompressor();
        do {
            if (d == NULL) {
                fprintf(stderr, "Cannot allocate decompressor.\n");
                exit(EXIT_FAILURE);
            }
            enum libdeflate_result res = libdeflate_gzip_decompress_ex(d, in + inPos, size - inPos, gzf.out + gzf.size,
                uncompressed_size - gzf.size, &in_size, &ret_size);
            if (res == LIBDEFLATE_INSUFFICIENT_SPACE) {
                uncompressed_size *= 2;
                if ((gzf.out = (uint8_t*) realloc(gzf.out, uncompressed_size)) == NULL) {
                    fprintf(stderr, "Error: out of memory during %s decompression.\n", filename);
                    exit(EXIT_FAILURE);
                };
                continue;
            }
            if (res != LIBDEFLATE_SUCCESS) {
                fprintf(stderr, "Error decompressing gz file: %d.\n", res);
                exit(EXIT_FAILURE);
            }
            inPos += in_size;
            gzf.size += ret_size;
        } while (inPos < size);
        free(in);
        libdeflate_free_decompressor(d);
    } else
    { //input is not compressed?
        gzf.out = in;
        gzf.size = size;
    }
    gzf.pos = 0;
    gzf.fileOut = gzf.out;
    gzf.fileSize = gzf.size;
    return gzf;
}

size_t mbgcInRead(mbgcInFile &gzf, void *buf, size_t nbyte) {
    if (gzf.size - gzf.pos < nbyte)
        nbyte = gzf.size - gzf.pos;
    memcpy(buf, gzf.out + gzf.pos, nbyte);
    gzf.pos += nbyte;
    return nbyte;
}

mbgcInFile& mbgcInReset(mbgcInFile &gzf) {
    gzf.pos = 0;
    return gzf;
};


int mbgcInClose(mbgcInFile &gzf) {
    free(gzf.fileOut);
    return 0;
}

mbgcInFile &mbgcInSplit_init(mbgcInFile &gzf) {
    gzf.fileOut = gzf.out;
    gzf.fileSize = gzf.size;
    gzf.size = 0;
    gzf.pos = 0;
    return gzf;
}

mbgcInFile mbgcInSplit_next(mbgcInFile &itergzf, size_t minSplitSize, char splitChar) {
    mbgcInFile gzf = itergzf;
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

string gzcompress(string &src, int clevel) {
    struct libdeflate_compressor* c;
    c = libdeflate_alloc_compressor(clevel);
    if (c == NULL) {
        fprintf(stderr, "Cannot allocate compressor.\n");
        exit(EXIT_FAILURE);
    }
    string dest;
    dest.resize(src.size() / 3);
    size_t res = 0;
    while (res == 0) {
        res = libdeflate_gzip_compress(c, src.data(), src.size(), dest.data(), dest.size());
        if (res == 0)
            dest.resize(dest.size() * 2);
    }
    dest.resize(res);
    libdeflate_free_compressor(c);
    return dest;
}


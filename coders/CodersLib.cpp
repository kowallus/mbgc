#include <assert.h>
#include "CodersLib.h"
#include "LzmaCoder.h"
#include "PpmdCoder.h"
#include "VarLenDNACoder.h"
#include <omp.h>

#ifdef DEVELOPER_BUILD
bool dump_after_decompression = false;
int dump_after_decompression_counter = 1;
string dump_after_decompression_prefix;
#endif

/*
Compress
------------

outPropsSize -
In:  the pointer to the size of outProps buffer; *outPropsSize = LZMA_PROPS_SIZE = 5.
Out: the pointer to the size of written properties in outProps buffer; *outPropsSize = LZMA_PROPS_SIZE = 5.

Out:
destLen  - processed output size
        Returns:
SZ_OK               - OK
SZ_ERROR_MEM        - Memory allocation error
SZ_ERROR_PARAM      - Incorrect paramater
SZ_ERROR_OUTPUT_EOF - output buffer overflow
SZ_ERROR_THREAD     - errors in multithreading functions (only for Mt version)
*/

/*
Uncompress
--------------
In:
  dest     - output data
  destLen  - output data size
  src      - input data
  srcLen   - input data size
Out:
  destLen  - processed output size
  srcLen   - processed input size
Returns:
  SZ_OK                - OK
  SZ_ERROR_DATA        - Data error
  SZ_ERROR_MEM         - Memory allocation arror
  SZ_ERROR_UNSUPPORTED - Unsupported properties
  SZ_ERROR_INPUT_EOF   - it needs more bytes in input buffer (src)
*/

using namespace PgHelpers;

unsigned char* Compress(size_t &destLen, const unsigned char *src, size_t srcLen, CoderProps* props,
        double estimated_compression, std::ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    unsigned char* dest = 0;
    unsigned char* compSeq = 0;
    int res = 0;
    switch (props->getCoderType()) {
        case LZMA_CODER:
            res = LzmaCompress(dest, destLen, src, srcLen,
                               ((LzmaCoderProps*) props)->getProps(), estimated_compression);
            break;
        case PPMD7_CODER:
            res = Ppmd7Compress(dest, destLen, src, srcLen, (PpmdCoderProps*) props,
                    estimated_compression);
            break;
        case VARLEN_DNA_CODER:
            res = VarLenDNACoder::Compress(dest, destLen, src, srcLen,
                    (VarLenDNACoderProps*) props);
            estimated_compression = VarLenDNACoder::COMPRESSION_ESTIMATION;
            break;
        case PARALLEL_BLOCKS_CODER_TYPE:
            res = parallelBlocksCompress(dest, destLen, src, srcLen,
                    (ParallelBlocksCoderProps*) props, estimated_compression, logout);
            break;
        case COMPOUND_CODER_TYPE:
            *logout << "\n\t\t";
            compSeq = Compress(destLen, (const unsigned char*)  src, srcLen,
                               ((CompoundCoderProps*) props)->primaryCoder, estimated_compression, logout);
            if (destLen >= srcLen) {
                destLen = srcLen;
                memcpy(compSeq, src, srcLen);
            }
            ((CompoundCoderProps*) props)->compLen = destLen;
            *logout << "\t\t";
            dest = Compress(destLen, (const unsigned char*)  compSeq, destLen,
                            ((CompoundCoderProps*) props)->secondaryCoder, estimated_compression, logout);
            if (destLen >= ((CompoundCoderProps*) props)->compLen) {
                destLen = ((CompoundCoderProps*) props)->compLen;
                swap(dest, compSeq);
            }
            *logout << "\t\t";
            delete[] compSeq;
            res = SZ_OK;
            break;
        case LZMA2_CODER:
        default:
            fprintf(stderr, "Unsupported coder type: %d.\n", props->getCoderType());
            exit(EXIT_FAILURE);
    }

    if (res != SZ_OK) {
        fprintf(stderr, "Error during compression (code: %d).\n", res);
        fprintf(stderr,"%s\tsrcLen: %zu\tdestLen: %zu\tcoder: %d\n",
                props->log().c_str(), srcLen, destLen, (int) props->getCoderType());
        exit(EXIT_FAILURE);
    }
    const double ratio = ((double) destLen) / srcLen;
    if (destLen < srcLen) {
        *logout << props->log() << " ... ";
        *logout << "compressed " << srcLen << " bytes to " << destLen << " bytes (ratio "
                << PgHelpers::toString(ratio, 3) << " vs estimated "
                << PgHelpers::toString(estimated_compression, 3) << ") in "
                << PgHelpers::time_millis(start_t) << " msec." << endl;
        if (ratio > estimated_compression)
            *logout << "WARNING: compression ratio " << PgHelpers::toString(ratio / estimated_compression, 5)
                    << " times greater than estimation." << endl;
    } else {
        *logout << " ignored" << props->log() << " due ratio >=1 ... " << srcLen << " bytes disregarding compression in "
                << PgHelpers::time_millis(start_t) << " msec." << endl;
    }
    return dest;
}

MY_STDAPI NoCoderUncompress(unsigned char *dest, size_t *destLen, unsigned char *src, size_t *srcLen, ostream* logout) {
    *destLen = *srcLen;
    *logout << "... raw ... ";
    memcpy(dest, (const void *) src, *destLen);
    return SZ_OK;
}

void Uncompress(unsigned char* dest, size_t destLen, istream &src, size_t srcLen, uint8_t coder_type,
        std::ostream* logout) {
    switch (coder_type) {
    default:
        string srcString;
        srcString.resize(srcLen);
        PgHelpers:readArray(src, (void*) srcString.data(), srcLen);
        Uncompress(dest, destLen, (unsigned char*) srcString.data(), srcLen, coder_type, logout);
    }
}

void Uncompress(unsigned char* dest, size_t destLen, unsigned char* src, size_t srcLen, uint8_t coder_type,
                std::ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int res = SZ_OK;
    size_t outLen = destLen;
    switch (coder_type) {
        case VARLEN_DNA_CODER:
            res = PgHelpers::VarLenDNACoder::Uncompress(dest, &outLen, src, &srcLen);
            break;
        case NO_CODER:
            res = NoCoderUncompress(dest, &outLen, src, &srcLen, logout);
            break;
        case LZMA_CODER:
            res = LzmaUncompress(dest, &outLen, src, &srcLen, logout);
            break;
        case PPMD7_CODER:
            res = PpmdUncompress(dest, &outLen, src, &srcLen, logout);
            break;
        case PARALLEL_BLOCKS_CODER_TYPE:
            res = parallelBlocksDecompress(dest, &outLen, src, &srcLen, logout);
            break;
        case LZMA2_CODER:
        default:
            fprintf(stderr, "Unsupported coder type: %d.\n", coder_type);
            exit(EXIT_FAILURE);
    }
    assert(outLen == destLen);

    if (res != SZ_OK) {
        fprintf(stderr, "Error during decompression (code: %d).\n", res);
        fprintf(stderr, "srcLen: %zu\tdestLen: %zu\tcoder: %d\n", srcLen, destLen, (int) coder_type);
        exit(EXIT_FAILURE);
    }
    *logout << "uncompressed " << srcLen << " bytes to " << destLen << " bytes in "
                         << PgHelpers::time_millis(start_t) << " msec." << endl;
}

void writeCompressed(ostream &dest, const char *src, size_t srcLen, CoderProps* props, double estimated_compression) {
    if (srcLen == 0) {
        PgHelpers::writeValue<uint64_t>(dest, 0);
        *PgHelpers::devout << "skipped compression (0 bytes)." << endl;
        return;
    }
    size_t compLen = 0;
    unsigned char* compSeq = Compress(compLen, (const unsigned char*) src, srcLen, props, estimated_compression);
    writeHeader(dest, srcLen, compLen, props);
    if (compLen < srcLen) {
        PgHelpers::writeArray(dest, (void*) compSeq, compLen);
    } else {
        PgHelpers::writeArray(dest, (void *) src, srcLen);
    }
    delete[] compSeq;
}

void writeCompressed(ostream &dest, const string& srcStr, CoderProps* props, double estimated_compression) {
    writeCompressed(dest, srcStr.data(), srcStr.length(), props, estimated_compression);
}

void writeHeader(ostream &dest, size_t srcLen, size_t destLen, CoderProps *props) {
    PgHelpers::writeValue<uint64_t>(dest, srcLen);
    if (props->getCoderType() == COMPOUND_CODER_TYPE) {
        PgHelpers::writeValue<uint64_t>(dest,destLen + ((CompoundCoderProps*) props)->primaryCoder->getHeaderLen());
        PgHelpers::writeValue<uint8_t>(dest, COMPOUND_CODER_TYPE);
        PgHelpers::writeValue<uint64_t>(dest, ((CompoundCoderProps*) props)->compLen);
        PgHelpers::writeValue<uint8_t>(dest, ((CompoundCoderProps*) props)->compLen == srcLen ?
            NO_CODER : ((CompoundCoderProps*) props)->primaryCoder->getCoderType());
        writeHeader(dest, ((CompoundCoderProps*) props)->compLen, destLen, ((CompoundCoderProps*) props)->secondaryCoder);
    } else if (destLen >= srcLen) {
        PgHelpers::writeValue<uint64_t>(dest, srcLen);
        PgHelpers::writeValue<uint8_t>(dest, NO_CODER);
    } else {
        PgHelpers::writeValue<uint64_t>(dest, destLen);
        PgHelpers::writeValue<uint8_t>(dest, props->getCoderType());
    }
}

void readCompressed(istream &src, string& dest, ostream* logout) {
    uint64_t destLen = 0;
    uint64_t srcLen = 0;
    uint8_t coder_type = 0;
    PgHelpers::readValue<uint64_t>(src, destLen);
    dest.resize(destLen);
    if (destLen == 0)
        return;
    PgHelpers::readValue<uint64_t>(src, srcLen);
    PgHelpers::readValue<uint8_t>(src, coder_type);
    if (coder_type == COMPOUND_CODER_TYPE) {
        uint64_t compLen = 0;
        uint8_t primary_coder_type = 0;
        PgHelpers::readValue<uint64_t>(src, compLen);
        PgHelpers::readValue<uint8_t>(src, primary_coder_type);
        string component;
        *logout << "\t\t";
        readCompressed(src, component, logout);
        assert(compLen == component.length());
        *logout << "\t\t";
        Uncompress((unsigned char *) dest.data(), destLen, (unsigned char*) component.data(), compLen, primary_coder_type, logout);
    } else {
        Uncompress((unsigned char *) dest.data(), destLen, src, srcLen, coder_type, logout);
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10?"0":"");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                      (void*) dest.data(), destLen);
    }
#endif
}

void readCompressed(unsigned char *src, string& dest, ostream* logout) {
    uint64_t destLen = 0;
    uint64_t srcLen = 0;
    uint8_t coder_type = 0;
    PgHelpers::readValue<uint64_t>(src, destLen);
    dest.resize(destLen);
    if (destLen == 0)
        return;
    PgHelpers::readValue<uint64_t>(src, srcLen);
    PgHelpers::readValue<uint8_t>(src, coder_type);
    if (coder_type == COMPOUND_CODER_TYPE) {
        uint64_t compLen = 0;
        uint8_t primary_coder_type = 0;
        PgHelpers::readValue<uint64_t>(src, compLen);
        PgHelpers::readValue<uint8_t>(src, primary_coder_type);
        string component;
        *logout << "\t\t";
        readCompressed(src, component, logout);
        assert(compLen == component.length());
        *logout << "\t\t";
        Uncompress((unsigned char *) dest.data(), destLen, (unsigned char*) component.data(), compLen, primary_coder_type, logout);
    } else {
        Uncompress((unsigned char *) dest.data(), destLen, src, srcLen, coder_type, logout);
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10?"0":"");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                      (void*) dest.data(), destLen);
    }
#endif
}

double simpleUintCompressionEstimate(uint64_t dataMaxValue, uint64_t typeMaxValue) {
    const double dataBits = 64 - ((double) __builtin_clzl(dataMaxValue));
    const double typeBits = 64 - __builtin_clzl(typeMaxValue);
    return dataBits / typeBits;
}

int parallelBlocksCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
        ParallelBlocksCoderProps *props, double estimated_compression, ostream* logout) {
    props->prepare(srcLen);
    size_t blockSize = ((srcLen / props->numOfBlocks) / props->blockAlignment) * props->blockAlignment;
    stringstream destOut;
    PgHelpers::writeValue(destOut, props->numOfBlocks);
    vector<CompressionJob> cJobs;
    size_t offset = 0;
    for (int i = 0; i < props->numOfBlocks - 1; i++) {
        cJobs.push_back(CompressionJob("block " + toString(i + 1) + "... ", src + offset, blockSize,
                props->blocksCoder));
        offset += blockSize;
    }
    cJobs.push_back(CompressionJob("block " + toString(props->numOfBlocks) + "... ", src + offset, srcLen - offset,
                    props->blocksCoder));
    CompressionJob::writeCompressedCollectiveParallel(destOut, cJobs, &null_stream);
    destLen = destOut.tellp();
    dest = new unsigned char[destLen];
    destOut.seekg(0);
    destOut.read((char*) dest, destLen);

    return SZ_OK;
}

int parallelBlocksDecompress(unsigned char *dest, size_t *destLen, unsigned char* src, size_t *srcLen, ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int numOfBlocks;
    PgHelpers::readValue<int>(src, numOfBlocks);
    *logout << "... parallel_blocks (no = " << numOfBlocks << ") of ";
    size_t offset = 0;
#pragma omp parallel
    {
#pragma omp single
        {
            for (int i = 0; i < numOfBlocks; i++) {
                unsigned char* destPtr = dest + offset;
                uint64_t blockLen = 0;
                uint64_t srcLen = 0;
                uint8_t coder_type = 0;
                unsigned char* srcString;
                PgHelpers::readValue<uint64_t>(src, blockLen);
                offset += blockLen;
                if (blockLen > 0) {
                    PgHelpers::readValue<uint64_t>(src, srcLen);
                    PgHelpers::readValue<uint8_t>(src, coder_type);
                    srcString = src;
                    src += srcLen;
                }
                if (blockLen == 0)
                    continue;
#pragma omp task
                {
                    ostringstream tmpout;
                    ostream* currentOut = i == 0?&tmpout:&null_stream;
                    Uncompress(destPtr, blockLen, srcString, srcLen, coder_type, currentOut);
                    if (i == 0) {
                        string log = tmpout.str();
                        *logout << log.substr(4, log.find("...", 4));
                    }
                }
            }
        }
    }
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10 ? "0" : "");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                    (void *) dest, *destLen);
    }
#endif
    return SZ_OK;
}

CompressionJob::CompressionJob(string label, const unsigned char *src, size_t srcLen,
        CoderProps *props, double estimated_compression) : log(label), src(src), srcLen(srcLen),
        props(props), estimated_compression(estimated_compression) {}

CompressionJob::CompressionJob(string label, const string &src,
        CoderProps *props, double estimated_compression) : log(label), src((const unsigned char*) src.data()),
        srcLen(src.length()), props(props), estimated_compression(estimated_compression) { }

void CompressionJob::writeCompressedCollectiveParallel(ostream &dest, vector<CompressionJob> &cJobs, ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    *logout << "collective compression of streams..." << endl;
    vector<unsigned char*> compSeqs;
    vector<size_t> compLens;
    compSeqs.resize(cJobs.size());
    compLens.resize(cJobs.size(), 0);
#ifdef __APPLE__
    omp_set_max_active_levels(4);
#else
    omp_set_nested(true);
#endif
#pragma omp parallel for
    for(int i = 0; i < cJobs.size(); i++) {
        ostringstream localLogOut;
        compSeqs[i] = Compress(compLens[i], cJobs[i].src, cJobs[i].srcLen, cJobs[i].props,
                cJobs[i].estimated_compression, &localLogOut);
        cJobs[i].log.append(localLogOut.str());
    }
    for(int i = 0; i < cJobs.size(); i++) {
        if (cJobs[i].srcLen == 0) {
            PgHelpers::writeValue<uint64_t>(dest, 0);
            *PgHelpers::devout << "skipped compression (0 bytes)." << endl;
            continue;
        }
        writeHeader(dest, cJobs[i].srcLen, compLens[i], cJobs[i].props);
        if (compLens[i] < cJobs[i].srcLen) {
            *logout << "\t" << cJobs[i].log;
            PgHelpers::writeArray(dest, (void *) compSeqs[i], compLens[i]);
        } else {
            *logout << "\t" << "skipped compression (" << cJobs[i].srcLen << " bytes)." << endl;
            PgHelpers::writeArray(dest, (void *) cJobs[i].src, cJobs[i].srcLen);
        }

        delete[] compSeqs[i];
    }
    *logout << "collective compression finished in " << PgHelpers::time_millis(start_t) << " msec." << endl;
    if (logout != &null_stream)
        *PgHelpers::logout << PgHelpers::time_millis(start_t) << "       \t";
}

void readCompressedCollectiveParallel(istream &src, vector<string*>& destStrings) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    *PgHelpers::devout << "collective decompression of streams..." << endl;
    vector<ostringstream> logOuts;
    logOuts.resize(destStrings.size());
#ifdef __APPLE__
    omp_set_max_active_levels(4);
#else
    omp_set_nested(true);
#endif
#pragma omp parallel
    {
#pragma omp single
        {
            for (int i = 0; i < destStrings.size(); i++) {
                uint64_t destLen = 0;
                uint64_t srcLen = 0;
                uint64_t compLen = 0;
                uint8_t coder_type = 0;
                uint8_t primary_coder_type = 0;
                string srcString;
                PgHelpers::readValue<uint64_t>(src, destLen);
                destStrings[i]->resize(destLen);
                if (destLen == 0)
                    continue;
                PgHelpers::readValue<uint64_t>(src, srcLen);
                PgHelpers::readValue<uint8_t>(src, coder_type);
                srcString.resize(srcLen);
                if (coder_type == COMPOUND_CODER_TYPE) {
                    PgHelpers::readValue<uint64_t>(src, compLen);
                    PgHelpers::readValue<uint8_t>(src, primary_coder_type);
                }
                PgHelpers::readArray(src, (void *) srcString.data(), srcLen);
#pragma omp task
                {
                    if (coder_type == COMPOUND_CODER_TYPE) {
                        string component;
                        logOuts[i] << "\t";
                        readCompressed((unsigned char*) srcString.data(), component, &logOuts[i]);
                        assert(compLen == component.length());
                        logOuts[i] << "\t\t";
                        Uncompress((unsigned char *) destStrings[i]->data(), destLen, (unsigned char*) component.data(), compLen,
                                   primary_coder_type, &logOuts[i]);
                    } else
                        Uncompress((unsigned char *) destStrings[i]->data(), destLen, (unsigned char*) srcString.data(), srcLen, coder_type,
                                   &logOuts[i]);
                }
            }
        }
    }
    for (int i = 0; i < destStrings.size(); i++) {
        *PgHelpers::devout << "\t" << logOuts[i].str();
#ifdef DEVELOPER_BUILD
        if (dump_after_decompression) {
            string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10 ? "0" : "");
            PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                                        (void *) destStrings[i]->data(), destStrings[i]->size());
        }
#endif
    }
    *PgHelpers::devout << "collective decompression finished in " << PgHelpers::time_millis(start_t) << " msec."
                       << endl;
    *PgHelpers::logout << PgHelpers::time_millis(start_t) << "       \t";
}


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
            ((CompoundCoderProps*) props)->compLen = destLen;
            *logout << "\t\t";
            dest = Compress(destLen, (const unsigned char*)  compSeq, destLen,
                            ((CompoundCoderProps*) props)->secondaryCoder, estimated_compression, logout);
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
        exit(EXIT_FAILURE);
    }
    *logout << props->log() << " ... ";
    const double ratio = ((double) destLen) / srcLen;
    *logout << "compressed " << srcLen << " bytes to " << destLen << " bytes (ratio "
         << PgHelpers::toString(ratio, 3) << " vs estimated "
         << PgHelpers::toString(estimated_compression, 3) << ") in "
         << PgHelpers::time_millis(start_t) << " msec." << endl;
    if (ratio > estimated_compression)
        *logout << "WARNING: compression ratio " << PgHelpers::toString(ratio / estimated_compression, 5)
        << " times greater than estimation." << endl;

    return dest;
}

void Uncompress(unsigned char* dest, size_t destLen, istream &src, size_t srcLen, uint8_t coder_type,
        std::ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int res = 0;
    size_t outLen = destLen;
    switch (coder_type) {
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
        exit(EXIT_FAILURE);
    }
    *logout << "uncompressed " << srcLen << " bytes to " << destLen << " bytes in "
         << PgHelpers::time_millis(start_t) << " msec." << endl;
}

void Uncompress(unsigned char* dest, size_t destLen, string& src, size_t srcLen, uint8_t coder_type) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int res = SZ_OK;
    size_t outLen = destLen;
    switch (coder_type) {
        case VARLEN_DNA_CODER:
            res = PgHelpers::VarLenDNACoder::Uncompress(dest, &outLen,
                                                        (unsigned char*) src.data(), &srcLen);
            break;
        default:
            istringstream srcStr(std::move(src));
            Uncompress(dest, outLen, srcStr, srcLen, coder_type);
    }
    assert(outLen == destLen);

    if (res != SZ_OK) {
        fprintf(stderr, "Error during decompression (code: %d).\n", res);
        exit(EXIT_FAILURE);
    }
    *PgHelpers::devout << "uncompressed " << srcLen << " bytes to " << destLen << " bytes in "
                         << PgHelpers::time_millis(start_t) << " msec." << endl;
}

void writeCompressed(ostream &dest, const char *src, size_t srcLen, CoderProps* props, double estimated_compression) {
    if (srcLen == 0) {
        PgHelpers::writeValue<uint64_t>(dest, 0, false);
        *PgHelpers::devout << "skipped compression (0 bytes)." << endl;
        return;
    }
    size_t compLen = 0;
    unsigned char* compSeq = Compress(compLen, (const unsigned char*) src, srcLen, props, estimated_compression);
    writeHeader(dest, srcLen, compLen, props);
    PgHelpers::writeArray(dest, (void*) compSeq, compLen);
    delete[] compSeq;
}

void writeCompressed(ostream &dest, const string& srcStr, CoderProps* props, double estimated_compression) {
    writeCompressed(dest, srcStr.data(), srcStr.length(), props, estimated_compression);
}

void writeHeader(ostream &dest, size_t srcLen, size_t destLen, CoderProps *props) {
    PgHelpers::writeValue<uint64_t>(dest, srcLen, false);
    if (props->getCoderType() == COMPOUND_CODER_TYPE) {
        PgHelpers::writeValue<uint64_t>(dest,destLen + ((CompoundCoderProps*) props)->primaryCoder->getHeaderLen(), false);
        PgHelpers::writeValue<uint8_t>(dest, COMPOUND_CODER_TYPE, false);
        PgHelpers::writeValue<uint64_t>(dest, ((CompoundCoderProps*) props)->compLen, false);
        PgHelpers::writeValue<uint8_t>(dest,
                                       ((CompoundCoderProps*) props)->primaryCoder->getCoderType(), false);
        writeHeader(dest, ((CompoundCoderProps*) props)->compLen, destLen, ((CompoundCoderProps*) props)->secondaryCoder);
    } else {
        PgHelpers::writeValue<uint64_t>(dest, destLen, false);
        PgHelpers::writeValue<uint8_t>(dest, props->getCoderType(), false);
    }
}

void readCompressed(istream &src, string& dest) {
    size_t destLen = 0;
    size_t srcLen = 0;
    uint8_t coder_type = 0;
    PgHelpers::readValue<uint64_t>(src, destLen, false);
    dest.resize(destLen);
    if (destLen == 0)
        return;
    PgHelpers::readValue<uint64_t>(src, srcLen, false);
    PgHelpers::readValue<uint8_t>(src, coder_type, false);
    if (coder_type == COMPOUND_CODER_TYPE) {
        size_t compLen = 0;
        uint8_t primary_coder_type = 0;
        PgHelpers::readValue<uint64_t>(src, compLen, false);
        PgHelpers::readValue<uint8_t>(src, primary_coder_type, false);
        string component;
        readCompressed(src, component);
        assert(compLen == component.length());
        Uncompress((unsigned char *) dest.data(), destLen, component, compLen, primary_coder_type);
    } else {
        Uncompress((unsigned char *) dest.data(), destLen, src, srcLen, coder_type);
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

int parallelBlocksDecompress(unsigned char *dest, size_t *destLen, istream &src, size_t *srcLen, ostream* logout) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int numOfBlocks;
    PgHelpers::readValue(src, numOfBlocks, false);
    *logout << "... parallel_blocks (no = " << numOfBlocks << ") of ";
    size_t offset = 0;
#pragma omp parallel
    {
#pragma omp single
        {
            for (int i = 0; i < numOfBlocks; i++) {
                unsigned char* destPtr = dest + offset;
                size_t blockLen = 0;
                size_t srcLen = 0;
                uint8_t coder_type = 0;
                string srcString;
                PgHelpers::readValue<uint64_t>(src, blockLen, false);
                offset += blockLen;
                if (blockLen > 0) {
                    PgHelpers::readValue<uint64_t>(src, srcLen, false);
                    PgHelpers::readValue<uint8_t>(src, coder_type, false);
                    srcString.resize(srcLen);
                    PgHelpers:readArray(src, (void*) srcString.data(), srcLen);
                }
                if (blockLen == 0)
                    continue;
#pragma omp task
                {
                    istringstream srcStream(std::move(srcString));
                    ostringstream tmpout;
                    ostream* currentOut = i == 0?&tmpout:&null_stream;
                    Uncompress(destPtr, blockLen, srcStream, srcLen, coder_type, currentOut);
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
    omp_set_nested(true);
#pragma omp parallel for
    for(int i = 0; i < cJobs.size(); i++) {
        ostringstream localLogOut;
        compSeqs[i] = Compress(compLens[i], cJobs[i].src, cJobs[i].srcLen, cJobs[i].props,
                cJobs[i].estimated_compression, &localLogOut);
        cJobs[i].log.append(localLogOut.str());
    }
    for(int i = 0; i < cJobs.size(); i++) {
        if (cJobs[i].srcLen == 0) {
            PgHelpers::writeValue<uint64_t>(dest, 0, false);
            *PgHelpers::devout << "skipped compression (0 bytes)." << endl;
            continue;
        }
        *logout << "\t" << cJobs[i].log;
        writeHeader(dest, cJobs[i].srcLen, compLens[i], cJobs[i].props);
        PgHelpers::writeArray(dest, (void *) compSeqs[i], compLens[i]);
        delete[] compSeqs[i];
    }
    *logout << "collective compression finished in " << PgHelpers::time_millis(start_t) << " msec." << endl;
}

void readCompressedCollectiveParallel(istream &src, vector<string*>& destStrings) {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    *PgHelpers::devout << "collective decompression of streams..." << endl;
    vector<ostringstream> logOuts;
    logOuts.resize(destStrings.size());
    omp_set_nested(true);
#pragma omp parallel
    {
#pragma omp single
        {
            for (int i = 0; i < destStrings.size(); i++) {
                size_t destLen = 0;
                size_t srcLen = 0;
                size_t compLen = 0;
                uint8_t coder_type = 0;
                uint8_t primary_coder_type = 0;
                string srcString;
                PgHelpers::readValue<uint64_t>(src, destLen, false);
                destStrings[i]->resize(destLen);
                if (destLen == 0)
                    continue;
                PgHelpers::readValue<uint64_t>(src, srcLen, false);
                PgHelpers::readValue<uint8_t>(src, coder_type, false);
                srcString.resize(srcLen);
                if (coder_type == COMPOUND_CODER_TYPE) {
                    PgHelpers::readValue<uint64_t>(src, compLen, false);
                    PgHelpers::readValue<uint8_t>(src, primary_coder_type, false);
                }
                PgHelpers::readArray(src, (void *) srcString.data(), srcLen);
#pragma omp task
                {
                    istringstream srcStream(std::move(srcString));
                    if (coder_type == COMPOUND_CODER_TYPE) {
                        string component;
                        readCompressed(srcStream, component);
                        assert(compLen == component.length());
                        Uncompress((unsigned char *) destStrings[i]->data(), destLen, component, compLen,
                                   primary_coder_type);
                    } else
                        Uncompress((unsigned char *) destStrings[i]->data(), destLen, srcStream, srcLen, coder_type,
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
}


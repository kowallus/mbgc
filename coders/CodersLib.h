#ifndef CODERS_LIB_H
#define CODERS_LIB_H

#include "../utils/helper.h"
#include <vector>

using namespace std;

#ifdef DEVELOPER_BUILD
extern bool dump_after_decompression;
extern int dump_after_decompression_counter;
extern string dump_after_decompression_prefix;
#endif

const static uint8_t NO_CODER = 0;
const static uint8_t LZMA_CODER = 1;
const static uint8_t LZMA2_CODER = 2;
const static uint8_t PPMD7_CODER = 3;
const static uint8_t VARLEN_DNA_CODER = 11;
const static uint8_t COMPOUND_CODER_TYPE = 77;
const static uint8_t PARALLEL_BLOCKS_CODER_TYPE = 88;

const static int LZMA_DATAPERIODCODE_8_t = 0;
const static int LZMA_DATAPERIODCODE_16_t = 1;
const static int LZMA_DATAPERIODCODE_32_t = 2;
const static int LZMA_DATAPERIODCODE_64_t = 3;
const static int LZMA_DATAPERIODCODE_128_t = 4;

const static uint8_t CODER_LEVEL_FAST = 1;
const static uint8_t CODER_LEVEL_NORMAL = 2;
const static uint8_t CODER_LEVEL_MAX = 3;

const static double COMPRESSION_ESTIMATION_UINT8_BITMAP = 0.125;
const static double COMPRESSION_ESTIMATION_BASIC_DNA = 0.250;
const static double COMPRESSION_ESTIMATION_VAR_LEN_DNA = 0.83;
const static double COMPRESSION_ESTIMATION_MIS_CNT = 0.5;
const static double COMPRESSION_ESTIMATION_MIS_SYM = 0.250;

const static int MINIMAL_PARALLEL_BLOCK_LENGTH = 1 << 20;
const static int PARALLEL_BLOCKS_ALIGNMENT = 1 << 4;

class CoderProps {
private:
    uint8_t coder_type;
protected:
    int numThreads;

public:
    CoderProps(uint8_t coder_type, int numThreads = 1) : coder_type(coder_type), numThreads(numThreads) {}

    virtual ~CoderProps() {}

    uint8_t getCoderType() const { return coder_type; }
    virtual size_t getHeaderLen() { return 2 * sizeof(uint64_t) + sizeof(uint8_t); }
    virtual int getNumThreads() const { return numThreads; }

    virtual string log() = 0;
};

double simpleUintCompressionEstimate(uint64_t dataMaxValue, uint64_t typeMaxValue);

unsigned char* Compress(size_t &destLen, const unsigned char *src, size_t srcLen, CoderProps* props,
        double estimated_compression = 1, ostream* logout = PgHelpers::devout);
void writeHeader(ostream &dest, size_t srcLen, size_t destLen, CoderProps* props);
void writeCompressed(ostream &dest, const char *src, size_t srcLen, CoderProps* props, double estimated_compression = 1);
void writeCompressed(ostream &dest, const string& srcStr, CoderProps* props, double estimated_compression = 1);

void Uncompress(unsigned char* dest, size_t destLen, istream &src, size_t srcLen, uint8_t coder_type,
        ostream* logout = PgHelpers::devout);
void Uncompress(unsigned char* dest, size_t destLen, unsigned char* src, size_t srcLen, uint8_t coder_type,
        ostream* logout = PgHelpers::devout);
void readCompressed(istream &src, string& dest, ostream* logout = PgHelpers::devout);
void readCompressed(unsigned char *src, string& dest, ostream* logout = PgHelpers::devout);

class CompoundCoderProps: public CoderProps {
public:

    CoderProps* primaryCoder;
    CoderProps* secondaryCoder;
    uint64_t compLen = 0;

    CompoundCoderProps(CoderProps* primaryCoderProps, CoderProps* secondaryCoderProps)
            : CoderProps(COMPOUND_CODER_TYPE, -1), primaryCoder(primaryCoderProps),
              secondaryCoder(secondaryCoderProps) {
        if (primaryCoderProps->getCoderType() == COMPOUND_CODER_TYPE) {
            fprintf(stderr, "Primary compound coder cannot be of compound coder type.\n");
            exit(EXIT_FAILURE);
        }
    }

    string log() {
        return " compound_coders (" + secondaryCoder->log() + " over " + primaryCoder->log() + ")";
    }

    size_t getHeaderLen() {
        return 3 * sizeof(uint64_t) + 2 * sizeof(uint8_t) + primaryCoder->getHeaderLen();
    }

    virtual ~CompoundCoderProps() { };

};

template<typename T>
void readCompressed(istream &src, vector<T>& dest) {
    uint64_t destLen = 0;
    uint64_t srcLen = 0;
    uint8_t coder_type = 0;
    PgHelpers::readValue<uint64_t>(src, destLen);
    if (destLen % sizeof(T)) {
        fprintf(stderr, "Invalid output size %zu for decompressing to the vector of %zu-byte elements",
                destLen, sizeof(T));
    }
    dest.resize(destLen / sizeof(T));
    if (destLen == 0)
        return;
    PgHelpers::readValue<uint64_t>(src, srcLen);
    PgHelpers::readValue<uint8_t>(src, coder_type);
    Uncompress((unsigned char*) dest.data(), destLen, src, srcLen, coder_type);
#ifdef DEVELOPER_BUILD
    if (dump_after_decompression) {
        string dumpFileName = dump_after_decompression_prefix + (dump_after_decompression_counter < 10?"0":"");
        PgHelpers::writeArrayToFile(dumpFileName + PgHelpers::toString(dump_after_decompression_counter++),
                dest.data(), destLen);
    }
#endif
}

class ParallelBlocksCoderProps: public CoderProps {
public:

    int numOfBlocks;
    CoderProps* blocksCoder;
    const uint32_t minBlockLength;
    const uint32_t blockAlignment;

    ParallelBlocksCoderProps(int numOfBlocks, CoderProps* blocksCoder,
            uint32_t minBlockLength = MINIMAL_PARALLEL_BLOCK_LENGTH,
            uint32_t blockAlignment = PARALLEL_BLOCKS_ALIGNMENT)
            : CoderProps(PARALLEL_BLOCKS_CODER_TYPE, -1), blocksCoder(blocksCoder),
              numOfBlocks(numOfBlocks), minBlockLength(minBlockLength),
              blockAlignment(blockAlignment) { }

    string log() {
        if (numOfBlocks > 1)
            return " parallel_blocks (no = " + to_string(numOfBlocks) + ") of" + blocksCoder->log();
        else
            return " (skipped blocks)" + blocksCoder->log();
    }

    virtual ~ParallelBlocksCoderProps() { };

    void prepare(size_t srcLen) {
        int maxBlocks = srcLen / minBlockLength;
        if (numOfBlocks > maxBlocks)
            numOfBlocks = maxBlocks;
        else if (numOfBlocks > maxBlocks)
            numOfBlocks = maxBlocks;
        if (numOfBlocks == 0)
            numOfBlocks = 1;
    }
};

int parallelBlocksCompress(unsigned char *&dest, size_t &destLen, const unsigned char *src, size_t srcLen,
                             ParallelBlocksCoderProps* props, double estimated_compression = 1,
                             ostream* logout = PgHelpers::devout);
int parallelBlocksDecompress(unsigned char *dest, size_t *destLen, unsigned char *src, size_t *srcLen,
        ostream* logout = PgHelpers::devout);

class CompressionJob {
private:
    string log;
    const unsigned char *src;
    size_t srcLen;
    CoderProps* props;
    double estimated_compression;

public:
    CompressionJob(string label, const unsigned char *src, size_t srcLen, CoderProps *props,
            double estimated_compression = 1);
    CompressionJob(string label, const string& src, CoderProps *props, double estimated_compression = 1);

    void appendLog(string txt) { log.append(txt); };

    static void writeCompressedCollectiveParallel(ostream &dest, vector<CompressionJob> &cJobs,
            ostream* logout = PgHelpers::devout);
};

void readCompressedCollectiveParallel(istream &src, vector<string*> &destStrings);

#endif
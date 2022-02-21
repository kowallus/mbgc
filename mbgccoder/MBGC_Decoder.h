#ifndef MBGC_MBGC_DECODER_H
#define MBGC_MBGC_DECODER_H

#include <iostream>

#include "../utils/helper.h"
#include "MBGC_Params.h"
#include <vector>
#include <atomic>

class MBGC_Decoder {
private:
    MBGC_Params *params;

    string outBuffer;

    char mbgcVersionMajor;
    char mbgcVersionMinor;
    char mbgcVersionRevision;

    string refStr;
    size_t refPos = 0;
    int rcStart;

    uint64_t refTotalLength;
    uint32_t filesCount;
    size_t totalFilesLength;
    uint32_t largestContigLength;
    uint32_t largestFileLength;

    vector<uint8_t> unmatchedFractionFactors;
    string literalStr;
    size_t literalPos = 0;
    string headersStr;
    size_t headersPos = 0;
    string headersTemplates;
    size_t hTemplatesPos = 0;
    string namesStr;
    size_t namesPos = 0;

    int fileIdx = 0;

    istringstream mapOffSrc, mapOff5thByteSrc, mapLenSrc, seqsCountSrc, matchingLocksPosSrc;

    bool isStdinMode();
    bool isStdoutMode();

    void writeDNA(const char *sequence, int64_t length);
    void decodeHeader(string& headerTemplate);
    bool moveToFile(const string& filepath, string& src, int thread_no, bool append = false);
    void initReference(const string &name);
    void decodeReference(const string &name);
    uint32_t decodeSequenceAndReturnUnmatchedChars(string &dest);
    void decodeFile(uint8_t unmatchedFractionFactor);
    void loadRef(const char *seqText, size_t seqLength, size_t i1);
    void extractFiles();

    void writeFilesParallelTask(int thread_no);
    void extractFilesParallel();
    static const int DEFAULT_WRITING_BUFFER_SIZE = 16;
    int64_t writingBufferSize = DEFAULT_WRITING_BUFFER_SIZE;
    vector<vector<string>> contentsBuf;
    vector<vector<string>> namesBuf;
    vector<uint32_t> in, out;
    vector<uint32_t> extractedFilesCount;
    bool isDecoding;

    void readParamsAndStats(istream &in);

public:

    explicit MBGC_Decoder(MBGC_Params *mbgcParams);

    void decode();

};


#endif //MBGC_MBGC_DECODER_H

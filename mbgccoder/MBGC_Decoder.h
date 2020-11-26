#ifndef MBGC_MBGC_DECODER_H
#define MBGC_MBGC_DECODER_H

#include <iostream>
#include <fstream>

#include "../utils/helper.h"
#include "MBGC_Params.h"
#include <vector>

class MBGC_Decoder {
private:
    MBGC_Params *params;

    string outBuffer;

    char mbgcVersionMajor;
    char mbgcVersionMinor;
    char mbgcVersionRevision;

    string refStr;
    int rcStart;

    uint32_t refTotalLength;
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

    istringstream mapOffSrc, mapLenSrc, seqsCountSrc;

    void writeDNA(const char *sequence, int64_t length);
    void decodeHeader(string& headerTemplate);
    bool moveToFile(const string& filepath, string& src);
    void decodeReference(const string &name);
    uint32_t decodeSequenceAndReturnUnmatchedChars(string &dest);
    void decodeFile(uint8_t unmatchedFractionFactor);
    void extractFiles();

    void writeFilesParallelTask(const int thread_no);
    void extractFilesParallel();
    static const int WRITING_BUFFER_SIZE = 32;
    vector<vector<string>> contentsBuf;
    vector<vector<string>> namesBuf;
    vector<uint8_t> in, out;
    bool isDecoding;

    void readParamsAndStats(fstream &fin);

public:

    MBGC_Decoder(MBGC_Params *mbgcParams);

    void decode();

};


#endif //MBGC_MBGC_DECODER_H

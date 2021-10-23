#ifndef MBGC_MBGC_ENCODER_H
#define MBGC_MBGC_ENCODER_H

#include <iostream>

#include "../utils/helper.h"
#include "../matching/SlidingWindowSparseEMMatcher.h"
#include "MBGC_Params.h"

class MBGC_Encoder {
private:
    MBGC_Params *params;

    SlidingWindowSparseEMMatcher* matcher;

    ostringstream seqsCountDest;

    int rcStart;
#ifdef DEVELOPER_BUILD
    int32_t currentRefExtGoal;
    uint32_t refExtLengthPerFile;
#endif
    uint64_t refFinalTotalLength;
    uint32_t filesCount = 0;
    uint32_t targetsCount = 0;
    size_t totalFilesLength = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize = 0;
    uint32_t largestFileLength = 0;

    string headersStr;
    string headersTemplates;

    string namesStr;

    size_t resCount = 0;
    string locksPosStream;
    string mapOffStream;
    string mapOff5thByteStream;
    string mapLenStream;

    vector<string> fileNames;
    vector<string> fileHeadersTemplates;
    vector<string> fileHeaders;

    vector<uint8_t> unmatchedFractionFactors;
    vector<string> targetRefExtensions;
    vector<string> targetLiterals;
    vector<ostringstream> targetMapOffDests;
    vector<ostringstream> targetMapOff5thByteDests;
    vector<ostringstream> targetMapLenDests;
    vector<uint32_t> targetSeqsCounts;
    vector<size_t> matchingLocksPos;

    bool compressToStdout();

    void processFileName(string &fileName);
    void updateHeadersTemplate(uint32_t fileIndex, const string& header);
    void processHeader(uint32_t fileIndex, const string &header);

    void loadRef(string& refName);
    void appendRef(string& refExtRes, const char* extPtr, size_t length);

    void processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes);
    size_t processMatches(vector<PgTools::TextMatch>& textMatches, char *destPtr, size_t destLen, int i);

    void applyTemplatesToHeaders();

    size_t prepareAndCompressStreams();
    void writeParamsAndStats(ostream &out) const;

    int claimedTargetsCount;
    int64_t processedTargetsCount;

    uint32_t masterTargetsStats = 0;
    uint32_t taskTargetsStats = 0;
    uint32_t masterRefExtensionsStats = 0;
    uint32_t taskRefExtensionsStats = 0;
    int readingThreadsCount;
    vector<uint32_t> out, in;
    void readFilesParallelTask(int thread_no);
    static const int READING_BUFFER_SIZE = 32;

    void interleaveOrderOfFiles();
    void loadFileNames();
    void initParallelEncoding();
    void finalizeParallelEncodingInSingleFastaFileMode();
    inline void encodeTargetSequence(int i);
    inline int finalizeParallelEncodingOfTarget();
    void tryClaimAndProcessTarget(bool calledByMaster);
    void encodeTargetsWithParallelIO();
    void encodeTargetsParallel();
    void encodeTargetsBruteParallel();

public:

    explicit MBGC_Encoder(MBGC_Params *mbgcParams);

    void encode();

};


#endif //MBGC_MBGC_ENCODER_H


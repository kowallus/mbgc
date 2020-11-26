#ifndef MBGC_MBGC_ENCODER_H
#define MBGC_MBGC_ENCODER_H

#include <iostream>
#include <fstream>

#include "../utils/helper.h"
#include "../matching/SparseEMMatcher.h"
#include "MBGC_Params.h"

class MBGC_Encoder {
private:
    MBGC_Params *params;

    SparseEMMatcher* matcher;

    ostringstream seqsCountDest;

    string refStr;
    int rcStart;

    int32_t currentRefExtGoal;
    uint32_t refExtLengthPerFile;

    uint32_t refTotalLength;
    uint32_t filesCount = 0;
    uint32_t targetsCount = 0;
    size_t totalFilesLength = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize = 0;
    uint32_t largestFileLength = 0;

    string literalStr;
    string headersStr;
    string headersTemplates;

    string namesStr;

    size_t resCount = 0;
    string mapOff;
    string mapLen;

    vector<string> fileNames;
    vector<string> fileHeadersTemplates;
    vector<string> fileHeaders;

    vector<uint8_t> unmatchedFractionFactors;
    vector<string> targetRefExtensions;
    vector<string> targetLiterals;
    vector<ostringstream> targetMapOffDests;
    vector<ostringstream> targetMapLenDests;
    vector<uint32_t> targetSeqsCounts;

    void processFileName(string &fileName);
    void updateHeadersTemplate(uint32_t fileIndex, const string& header);
    void processHeader(uint32_t fileIndex, const string &header);

    void loadRef(string& refName);
    void appendRef(string& refExtRes, const char* extPtr, size_t length);

    void processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes);
    size_t processMatches(vector<PgTools::TextMatch>& textMatches, char *destPtr, size_t destLen,
                          string& literalRes, ostringstream& mapOffDest, ostringstream& mapLenDest, string& refExtRes);

    void buildHeadersTemplates();
    void applyTemplatesToHeaders();

    size_t prepareAndCompressStreams();
    void writeParamsAndStats(fstream &out) const;

    void encodeTargets(ifstream& listSrc);

    int claimedTargetsCount;
    int processedTargetsCount;

    uint32_t masterTargetsStats = 0;
    uint32_t taskTargetsStats = 0;
    uint32_t masterRefExtensionsStats = 0;
    uint32_t taskRefExtensionsStats = 0;
    int readingThreadsCount;
    vector<uint32_t> out, in;
    void readFilesParallelTask(const int thread_no);
    static const int READING_BUFFER_SIZE = 32;

    void interleaveOrderOfFiles();
    void loadFileNames(ifstream &listSrc);
    void initParallelEncoding();
    void finalizeParallelEncoding();
    inline void finalizeRefExtensionsOfTarget(int i);
    inline void encodeTargetSequence(int i);
    inline int finalizeParallelEncodingOfTarget();
    void tryClaimAndProcessTarget(bool calledByMaster);
    void encodeTargetsWithParallelIO();
    void encodeTargetsParallel();
    void encodeTargetsBruteParallel();

public:

    MBGC_Encoder(MBGC_Params *mbgcParams);

    void encode();

};


#endif //MBGC_MBGC_ENCODER_H


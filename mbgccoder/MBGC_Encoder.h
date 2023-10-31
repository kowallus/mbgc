#ifndef MBGC_MBGC_ENCODER_H
#define MBGC_MBGC_ENCODER_H

#include <iostream>

#include "../utils/helper.h"
#include "../matching/SlidingWindowSparseEMMatcher.h"
#include "MBGC_Params.h"
#include "../coders/ContextAwareMismatchesCoder.h"
#include "MBGC_Decoder.h"
#include <unordered_set>

class MBGC_Encoder {
private:
    MBGC_Params *params;
    static constexpr timespec SLEEP_TIME[] {{0, 100L}};

    SlidingWindowSparseEMMatcher* matcher;

    ContextAwareMismatchesCoder* mismatchesCoder = &ContextAwareMismatchesCoder::defaultInstance;

    ostringstream seqsCountsDest, dnaLineLengthsDest;

    uint64_t refG0InitPos;
#ifdef DEVELOPER_BUILD
    int64_t currentRefExtGoal;
    uint64_t refExtLengthPerFile;
#endif
    uint64_t refFinalTotalLength;
    uint8_t refBuffersCount = 1;
    uint32_t filesCount = 0;
    bool singleFastaFileMode;
    uint32_t targetsCount = 0;
    int targetsShift = 1;
    uint64_t totalFilesLength = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize = 0;
    uint64_t largestFileLength = 0;

    uint32_t appendedFilesCount = 0;
    uint32_t appendedTargetsCount = 0;

    string headersStr;
    string headersTemplates;

    string namesStr;

    size_t resCount = 0;
    string locksPosStream;
    string mapOffStream;
    string mapLenStream;
    string singularGapFlagsStream;
    string gapFlagsStream;
    string gapMismatchesFlagsStream;

    vector<string> fileNames;
    unordered_set<string> fileNamesSet;
    vector<string> fileHeadersTemplates;
    vector<string> fileHeaders;

    vector<uint8_t> unmatchedFractionFactors;
    vector<string> targetRefExtensions;
    vector<string> targetLiterals;
    vector<ostringstream> targetMapOffDests;
    vector<string> targetMapOff5thByte;
    vector<ostringstream> targetMapLenDests;
    vector<string> targetGapDeltas;
    vector<string> targetGapMismatchesFlags;
    vector<uint32_t> seqsCounts;
    vector<uint64_t> dnaLineLength;
    vector<size_t> matchingLocksPos;
    ostringstream refExtSizeDest;
    vector<size_t> refExtLoadedPosArr;

    vector<vector<string>*> usedByteStreams;
    void enrollByteStream(vector<string>& stream) {
        usedByteStreams.push_back(&stream);
    }
    vector<pair<string&, vector<ostringstream>&>> usedStreamsWithDests;
    void enrollStream(string& stream, vector<ostringstream>& dests) {
        usedStreamsWithDests.emplace_back(pair<string& , vector<ostringstream>&>(stream, dests));
    }

    void processFileName(string &fileName);
    void updateHeadersTemplate(uint32_t fileIndex, const string& header);
    void processHeader(uint32_t fileIndex, const string &header);

    void loadG0Ref(string& refName);
    void initMatcher(const char* refStrPtr, const size_t refStrSize, size_t basicRefLength);

    size_t getMatchLoadedPos(size_t pos);

    void processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes);
    size_t processMatches(vector<PgTools::TextMatch>& textMatches, char *destStart, size_t destLen, int targetIdx,
        size_t matchingLockPos = SIZE_MAX);
    uint64_t extendMatchLeft(const char *destStart, uint64_t length, const TextMatch &match,
                         int targetIdx, size_t matchingLockPos,
                         uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches);

    uint64_t extendMatchRight(const char *gapStartPtr, const TextMatch &coreMatch, const TextMatch &match,
                              uint64_t &length, int targetIdx, bool isGap, bool gapStart, bool gapMiddle, bool gapEnd,
                              uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches);

    void performMatching();

    void applyTemplateToHeaders(uint32_t fileIdx);
    void prepareHeadersStreams();
    void prepareAndCompressStreams();
    void writeStats(ostream &out) const;

    int claimedTargetsCount = 0;
    int64_t processedTargetsCount = 0;
    int allowedTargetsOutrun = INT32_MAX;

    uint32_t masterTargetsStats = 0;
    uint32_t taskTargetsStats = 0;
    uint32_t masterRefExtensionsStats = 0;
    uint32_t taskRefExtensionsStats = 0;
    int readingThreadsCount;
    vector<uint32_t> out, in;
    void readFilesParallelTask(int thread_no);
    static const int DEFAULT_READING_BUFFER_SIZE = 16;
    int readingBufferSize = DEFAULT_READING_BUFFER_SIZE;

    void interleaveOrderOfFiles();
    void loadFileNames();
    void initParallelEncoding();
    void finalizeParallelEncodingInSingleFastaFileMode();
    inline void encodeTargetSequence(int i, bool calledByMaster);
    inline int finalizeParallelEncodingOfTarget(bool calledByMaster);
    void tryClaimAndProcessTarget(bool calledByMaster);
    void encodeTargetsWithParallelIO();
    void encodeTargetsParallel();
    void encodeTargetsBruteParallel();

    template<bool lazyMode>
    void appendTemplate(MBGC_Decoder<lazyMode>* decoder, MBGC_Params& newParams);

    template<bool lazyMode>
    void repackTemplate(MBGC_Decoder<lazyMode>* decoder);

public:

    explicit MBGC_Encoder(MBGC_Params *mbgcParams);

    void encode();

    void append(MBGC_Params& newParams);

    void repack(MBGC_Params& inParams);

};


#endif //MBGC_MBGC_ENCODER_H


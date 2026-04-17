#ifndef MULTIPLEGENOMEMATCHINGPROCESSOR_H
#define MULTIPLEGENOMEMATCHINGPROCESSOR_H

#include "../utils/helper.h"
#include <unordered_set>

#include "MGMP_Params.h"

#include "SlidingWindowSparseEMMatcher.h"

class MultipleGenomeMatchingProcessor {
private:
    MGMP_Params *params;
    static constexpr timespec SLEEP_TIME[] {{0, 100L}};

    int claimedTargetsCount = 0;
    int allowedTargetsOutrun = INT32_MAX;

#ifdef DEVELOPER_BUILD
    int64_t currentRefExtGoal;
    uint64_t refExtLengthPerFile;
#endif
    uint32_t masterTargetsStats = 0;
    uint32_t taskTargetsStats = 0;
    uint32_t masterRefExtensionsStats = 0;
    uint32_t taskRefExtensionsStats = 0;
    int readingThreadsCount;
    vector<uint32_t> out, in;
    void readFilesParallelTask(int thread_no);

    void interleaveOrderOfFiles();
    inline int finalizeParallelProcessingOfTarget(bool calledByMaster);
    inline void processTarget(int i, bool calledByMaster);
    void tryClaimAndProcessTarget(bool calledByMaster);
    void processTargetsWithParallelIO();
    void processTargetsParallel();
    void processTargetsBruteParallel();

protected:
    SlidingWindowSparseEMMatcher* matcher;

    vector<string> fileNames;
    unordered_set<string> fileNamesSet;
    uint32_t filesCount = 0;
    bool singleFastaFileMode;
    uint32_t targetsCount = 0;
    uint32_t targetsCountShift = 0;
    int targetsShift = 1;

    static const int DEFAULT_READING_BUFFER_SIZE = 16;
    int readingBufferSize = DEFAULT_READING_BUFFER_SIZE;

    int64_t processedTargetsCount = 0;
    size_t openFirstFile(const char* filename, bool isTarget = false);
    void switchToSplitIterMode();

    void validate_kseq_status(string& fileName, const char* seq_name, int kseq_status);
    virtual void print_invalid_kseq_status_message(int kseq_status) = 0;

    uint64_t refG0InitPos;
    uint64_t refFinalTotalLength;
    uint8_t refBuffersCount = 1;
    uint64_t totalFilesLength = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize = 0;
    uint64_t largestFileLength = 0;

    size_t resCount = 0;

    vector<uint32_t> seqsCounts;
    vector<uint64_t> dnaLineLength;
    vector<size_t> matchingLocksPos;
    vector<uint8_t> unmatchedFractionFactors;
    vector<string> targetRefExtensions;

    void initMatcher(const char* refStrPtr, const size_t refStrSize, size_t basicRefLength,
        bool useOfExistingReference = false);
    void loadG0Ref(string& refName);
    virtual void setProteinsProfile() = 0;
    virtual void initStreamsForG0Ref() = 0;
    virtual void processG0RefContig(const char *seq, size_t len) = 0;

    size_t unmatchedCharsAll = 0;
    size_t totalMatchedAll = 0;
    size_t totalDestOverlapAll = 0;
    size_t totalDestLenAll = 0;

    void performMatching();

    virtual void printAdditionalMatchingStats() = 0;

    virtual void processFileName(string &fileName) = 0;
    virtual void processTargetMeta(uint32_t seqCount, uint64_t dnaLineLength) = 0;
    virtual void processHeader(uint32_t fileIndex, const string &header) = 0;
    virtual void processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes) = 0;

    const size_t PROCESSING_MATCHES_SKIPPED_DUE_TO_CONTIG_DISSIMILARITY = SIZE_MAX;

    virtual size_t processMatches(vector<PgTools::TextMatch>& textMatches, char *destStart, size_t destLen, int targetIdx,
        size_t matchingLockPos = SIZE_MAX) = 0;

    virtual void initProcessTargetsWithParallelIO() = 0;
    virtual void finalizeProcessTargetsWithParallelIO() = 0;
    virtual void processAfterTargetWithParallelIO(size_t matcherLoaderStartPos) = 0;
    virtual void processAfterTarget(uint32_t targetIdx) = 0;
    virtual void processAfterSequence(uint32_t targetIdx) = 0;

    virtual void initProcessTarget(uint32_t targetIdx) = 0;
    virtual void initParallelProcessing();
    virtual void finalizeParallelProcessingInSingleFastaFileMode();
    virtual void finalizeParallelProcessingOfTarget(uint32_t targetIdx, size_t matcherLoaderStartPos) = 0;

    virtual void applyTemplateToHeaders(uint32_t fileIdx) = 0;

    virtual void loadFileNames();


public:

    explicit MultipleGenomeMatchingProcessor(MGMP_Params *mgmpParams);

};

#endif //MULTIPLEGENOMEMATCHINGPROCESSOR_H

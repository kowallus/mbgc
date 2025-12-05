#ifndef MBGC_MBGC_DECODER_H
#define MBGC_MBGC_DECODER_H

#include "MBGC_API.h"

#include "../utils/helper.h"
#include "MBGC_Params.h"
#include <vector>
#include <atomic>
#include "../coders/ContextAwareMismatchesCoder.h"

static const int NO_TARGET_TO_SCHEDULE = -1;

template<bool lazyMode>
class MBGC_Decoder: public MBGC_Decoder_API {
private:
    static const int MASTER_THREAD_ID = 0;
    static constexpr timespec SLEEP_TIME[] {{0, 100L}};
    static constexpr timespec WRITE_THREAD_SLEEP_TIME[] {{0, 1000L}};

    MBGC_Params *params;
    istream* inStream;
    bool sequentialDecoding = false;

    ContextAwareMismatchesCoder* mismatchesCoder = &ContextAwareMismatchesCoder::defaultInstance;

    string g0name;
    vector<string *> destStrings;
    vector<uint32_t> mapLenArr;

    vector<string> threadOutBuffer;
    string refStr;
    size_t refPos = 0;
    const size_t REF_SHIFT = 1;
    uint64_t refG0InitPos;
    uint8_t refBuffersCount = 1;
    int reachedRefLengthCount = 0;

    size_t getLoadedRefLength() { return reachedRefLengthCount * (refTotalLength - REF_SHIFT) + (refPos - REF_SHIFT); };

    uint64_t refTotalLength;
    uint32_t filesCount;
    uint64_t totalFilesLength;
    uint32_t largestContigLength;
    uint64_t largestFileLength;
    uint64_t refContributingFilesCount = 0;

#ifdef DEVELOPER_BUILD
    size_t matchesPerTargetEstimate;
    vector<vector<uint64_t>> threadLitRefExtPos;
    vector<vector<uint64_t>> threadLitRefExtLen;
    void processLiteral(size_t pos, size_t literalLen);
    void prepareLiteralRefExtension();
#endif
    string literalStr;
    vector<int64_t> tIdLiteralPos;
    uint64_t literalsLength = SIZE_MAX;
    string headersStr;
    vector<size_t> tIdHeadersPos;
    string headersTemplates;
    vector<size_t> tIdHTemplatesPos;
    string namesStr;
    vector<size_t> tIdNamesPos;

    vector<string> threadSeqStr, threadExtStr, threadGzOut;

    string seqsCount, dnaLineLengthsStr, unmatchedFractionFactorsStr;
    string gapDeltaStr, gapMismatchesFlagsStr;
    string mapOffStream, mapOff5thByteStream, mapLenStream, refExtSizeStream;
    string matchingLocksPosStream;
    uint8_t* unmatchedFractionFactorsArray;
    size_t targetsCount;
    vector<uint8_t*> tIdMapOff5thBytePtr;
    uint8_t* mapOffRefIDArr;
    vector<WrapperStrStream> tIdMapOffSrc;
    vector<uint32_t*> tIdMapLenPtr;
    uint32_t* seqsCountArr;
    size_t* dnaLineLengthsArr;
    size_t* matchingLocksPosArr;
    uint8_t* gapDeltaGuard;
    uint8_t* gapMismatchesFlagsGuard;
    vector<uint8_t*> tIdGapDeltaPtr;
    vector<uint8_t*> tIdGapMismatchesFlagsPtr;
    vector<size_t> refExtLoadedPosArr;
    vector<size_t> refExtPosArr;

    void writeDNA(const char *sequence, const int64_t length, const int64_t dnaLineLength, int threadId);
    void decodeHeader(string& headerTemplate, int tId);
    void moveToFile(const string& filepath, string& src, int thread_no, bool append = false);
    void initReference(const string &name);
    void decodeReference(const string &name);


    int64_t pairedGapOffsetDeltaInit[MBGC_Params::MAX_GAP_DEPTH];
    uint32_t decodeSequenceAndReturnUnmatchedChars(string &dest, size_t refLockPos, const int tId);

    vector<string> threadMatchExtension;
    size_t extendMatchLeft(char *extStrEndPtr, uint64_t& matchSrcPos, bool skipOffset,
                           size_t refLockPos, size_t markPos, const int tId);
    size_t extendMatchRight(string& dest, int64_t& offsetDelta,
                            bool isGap, bool gapStart, bool gapMiddle, bool gapEnd,
                            size_t guardLitPos, const int tId);

    void finalizeLazyTargetDecoding(int targetIdx, bool targetActuallyDecoded);
    void decodeTarget(size_t targetIdx, int tId, bool movePartsToFile = false);
    void loadRef(const char *seqText, size_t seqLength, size_t refLockPos, size_t& refPos1, bool loadRCRef);

    size_t namePos = 0;
    string currentName;
    void extractFilesSequentially();

    void extractFilesParallel();

    void workerParallelTask(int thread_no);

    int masterThreadTargetIdx;
    int lazyDecompressionTargetsGuard = 0;
    int scheduledTargetsGuard = 0;
    int scheduledTargetsCount = 0;
    vector<bool> isTargetScheduled;
    vector<bool> isTargetClaimed;
    int targetsDecodedBeforeScheduling;
    int lowestTargetForProcessing = 0;
    vector<bool> isTargetProcessed;
    vector<bool> isTargetForParallelTaskThreadDecoding;
    vector<int> targetsForParallelTaskThreadDecoding;
    int64_t lowestNotClaimedParallelIdx = INT64_MAX;
    vector<int> highestTargetMapping;
    vector<int> tIdProcessesTarget;
    int64_t lowestContributingNotDecodedTargetBeforeScheduling;
    int64_t lowestContributingNotDecodedTarget = 0;
    int64_t decodedTargetsCount = 0;
    bool tryClaimingTarget(int tIdx);
    int tryClaimingTarget();

    uint32_t masterTargetsStats = 0;
    static const int DEFAULT_WORKER_THREADS_COUNT = 4;
    static const int DEFAULT_WRITING_THREADS_COUNT = 2;
    int8_t workerThreadsCount = lazyMode ? DEFAULT_WORKER_THREADS_COUNT : DEFAULT_WRITING_THREADS_COUNT;
    static const int DEFAULT_WRITING_BUFFER_SIZE = lazyMode ? 16 : 32;
    int64_t writingBufferSize = DEFAULT_WRITING_BUFFER_SIZE;
    vector<vector<string>> contentsBuf;
    vector<vector<string>> namesBuf;
    vector<uint32_t> in, out;
    vector<uint32_t> extractedFilesCount;
    bool isDecoding;

    void readStats(istream &in);

    void copyStreamsPositions(const int srcId, const int destId);
    void initLazyDecompression();
    void fastForwardTargetStreams(int targetIdx, int tId);
    int findStreamsPositions(int prevTargetSrcId, int startTargetIdx, int endTargetIdx, size_t namesEndGuard = SIZE_MAX);
    int prepareStreamsForNextFileToDecompress(const int streamsSrcId, int streamsTargetIdx);
    void scheduleDependencyChain(int64_t tId);
    void scheduleParallelLazyDecompression();
    void printSchedulingStats();

    void decodeMapLenStream(string &mapLenStream, vector<uint32_t> &mapLenArr);
    void decodeRefExtSizeStream(string &refExtSizeStream);

    bool isFilterListActive = false;
    void processFilterPatterns();
    vector<uint8_t> isFileSelectedFlag;
    uint8_t* isTargetSelectedPtr = nullptr;
    size_t namesEndGuard = std::string::npos;
    static const int FIRST_FILE_TARGET_IDX = -1;
    bool isTargetSelected(const string &filename, int targetIdx);
    bool isFileSelected(const string &filename, int fileIdx);
    void moveToFileIfSelected(const string& filepath, string& src, int thread_no,  int targetIdx = FIRST_FILE_TARGET_IDX, bool append = false);
    void getSelectedFilesNames(vector<string> &filenames);

    void list();

    void decodeInit();
    string iterCurFilename;
    string iterCurFileContent;
    void iterateInit();

public:

    explicit MBGC_Decoder(MBGC_Params *mbgcParams, istream* inStream);
    ~MBGC_Decoder();

    static MBGC_Decoder_API* getInstance(MBGC_Params* params, string optionalWarning = "");;

    virtual void decode();

    virtual bool iterateNext(const char* name, uint8_t* &res, size_t &size);

    friend class MBGC_Encoder;

    void extractNextSequentially(int tId);
};

#endif //MBGC_MBGC_DECODER_H

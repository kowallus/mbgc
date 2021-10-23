#ifndef PGTOOLS_SWSMEMMATCHER_H
#define PGTOOLS_SWSMEMMATCHER_H

#include "TextMatchers.h"
#include <deque>

using namespace PgTools;

enum verbosity { v0, v1, v2 };
enum reverseMode { no, yes, both };

static const int OVERLAP_MATCH_MAX_LENGTH = 1 << 13;
static const int HASH_SIZE_MIN_ORDER = 24;
static const int HASH_SIZE_MAX_ORDER = 31;
typedef std::pair<std::string, size_t> SequenceItem;
typedef std::vector<SequenceItem> SequenceVector;

template<class MyUINT1, class MyUINT2>
using HashBuffer =  std::pair<MyUINT1*, MyUINT2* >;

class SlidingWindowSparseEMMatcher {
private:
    char* start1;
    int64_t pos1;
    const size_t maxRefLength;
    int reachedRefLengthCount = 0;
    const size_t REF_SHIFT = 1;
    size_t samplingPos = REF_SHIFT;
    int bigRef;
    const int L;
    int K, k1, k2, skipMargin;
    std::uint32_t(*hashFunc32)(const char*);
    std::uint32_t(*hashFuncMatrix[64][6])(const char*);
    const int H = 1;
    std::uint32_t hash_size;
    std::uint32_t hash_size_minus_one;

    int LK2 = (L - K) / 2;
    int LK2_MINUS_4 = LK2 - 4;
    int K_PLUS_LK24 = K + LK2_MINUS_4;

    // sliding window properties
    size_t swSize, swEnd;
    static const size_t SW_END_ERASED_FLAG = SIZE_MAX;
    static const uint8_t SW_WIDTH_FACTOR = 16;
    std::deque<size_t> workersSwEndPositions;

    void initHashFuncMatrix();
    void initParams(uint32_t minMatchLength, int _k1 = -1, int _k2 = -1);
    void calcCoprimes();
    void displayParams();

    inline std::uint32_t hashFunc(const char* str) { return hashFunc32(str) & hash_size_minus_one; };

    std::pair<std::uint64_t*, std::uint64_t*> buffer2;
    std::pair<std::uint64_t*, std::uint32_t*> buffer1;
    std::pair<std::uint32_t*, std::uint32_t*> buffer0;

    template<typename MyUINT1, typename MyUINT2>
    void processIgnoreCollisionsRef(HashBuffer<MyUINT1, MyUINT2>& buffer);

    template <class MyUINT1, class MyUINT2>
    void deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf);

    template<typename MyUINT1, typename MyUINT2, bool checkOverlaps>
    void processExactMatchQueryIgnoreCollisionsTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches,
                             const char* start2, size_t N2,
                             bool destIsRef, bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos);

    template<typename MyUINT1, typename MyUINT2>
    void processExactMatchQueryIgnoreCollisionsTightTemplate(HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches,
                                                     const char* start2, size_t N2,
                                                     bool destIsRef, bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos);

public:
    SlidingWindowSparseEMMatcher(const size_t refLengthLimit, const uint32_t targetMatchLength,
                                 int _k1 = -1, int _k2 = -1, int skipMargin = 0, uint32_t minMatchLength = UINT32_MAX);

    void disableSlidingWindow() { swSize = 0; swEnd = 0; };

    void setSlidingWindowSize(uint8_t factor) { swSize = this->maxRefLength / factor; };

    void loadRef(const char *refText, size_t refLength);

    size_t getRefLength() { return reachedRefLengthCount ? maxRefLength : pos1; };
    size_t getLoadedRefLength() { return reachedRefLengthCount * maxRefLength + pos1; };

    size_t acquireWorkerMatchingLockPos();
    void releaseWorkerMatchingLockPos(size_t lockValue);

    virtual ~SlidingWindowSparseEMMatcher();

    void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsRef, bool revComplMatching,
                    uint32_t minMatchLength, size_t matchingLockPos = SW_END_ERASED_FLAG);

    void matchTexts(vector<TextMatch> &resMatches, const char* destText, size_t destLen, bool destIsRef,
            bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos = SW_END_ERASED_FLAG);

};

#endif //PGTOOLS_SWSMEMMATCHER_H

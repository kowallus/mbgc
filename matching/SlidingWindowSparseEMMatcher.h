#ifndef PGTOOLS_SWSMEMMATCHER_H
#define PGTOOLS_SWSMEMMATCHER_H

#include "TextMatchers.h"
#include <deque>

using namespace PgTools;

typedef std::pair<std::string, size_t> SequenceItem;
typedef std::vector<SequenceItem> SequenceVector;

class SlidingWindowSparseEMMatcher {
public:
    const size_t REF_SHIFT = 1;

private:

    static const int OVERLAP_MATCH_MAX_LENGTH = 1 << 13;
    static const int HASH_SIZE_MIN_ORDER = 24;
    static const int HASH_SIZE_MAX_ORDER = 31;
    static const int TOTAL_RAM_REF_LIMIT_PERCENT = 60;

    static const int MIN_MATCH_LENGTH = 16;

    enum verbosity { v0, v1, v2 };
    enum reverseMode { no, yes, both };

protected:
    char* start1;
    int64_t pos1;
    size_t maxRefLength;
    int reachedRefLengthCount = 0;
    bool ht64bitFlag = false;
    const int L;
    int K, k1, k2, skipMargin;
    std::uint32_t(*hashFunc32)(const char*);
    std::uint32_t(*hashFuncMatrix[64][6])(const char*);
    const int H = 1;
    std::uint32_t hash_size;
    std::uint32_t hash_size_minus_one;

    int LK2, LK2_MINUS_4, K_PLUS_LK24;

    // sliding window properties
    size_t swSize, swEnd;
    static const size_t SW_END_ERASED_FLAG = SIZE_MAX;
    static const uint8_t SW_WIDTH_FACTOR = 16;
    std::deque<size_t> workersSwEndPositions;
    bool circularBuffer = true;

    void initHashFuncMatrix();
    void initParams(uint32_t minMatchLength, int _k1 = -1, int _k2 = -1);
    void calcCoprimes();
    void displayParams();

    inline std::uint32_t hashFunc(const char* str) { return hashFunc32(str) & hash_size_minus_one; };

    std::uint64_t* ht64bit;
    std::uint32_t* ht32bit;

    template<typename uint_ht_t>
    void processIgnoreCollisionsRef(uint_ht_t* ht, size_t& refSamplingPos, int64_t& refPos1);

    template<typename uint_ht_t, bool checkOverlaps>
    void processExactMatchQueryIgnoreCollisionsTight(uint_ht_t* ht, vector<TextMatch> &resMatches,
                            const char* start2, size_t N2,
                            bool destIsRef, bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos);

    template<typename uint_ht_t>
    void processExactMatchQueryIgnoreCollisionsTightTemplate(uint_ht_t* ht, vector<TextMatch> &resMatches,
                             const char* start2, size_t N2,
                             bool destIsRef, bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos);

    virtual inline uint64_t htRePos(uint64_t pos) const { return pos; };
    virtual inline uint64_t htEncodePos(uint64_t pos) const { return pos; };
    virtual inline uint64_t htDecodePos(uint64_t code) const { return code; };

    size_t samplingPos = REF_SHIFT;

    SlidingWindowSparseEMMatcher(const size_t refLengthLimit, const uint32_t targetMatchLength,
                                 int _k1, int _k2, int skipMargin, uint32_t minMatchLength, bool skipHtInit);

    size_t loadRef(const char *refText, size_t refLength, size_t& refSamplingPos, int64_t& refPos1,
                 bool rcRef, bool addRegionSeparators, char regionSeparator);

public:

    SlidingWindowSparseEMMatcher(const size_t refLengthLimit, const uint32_t targetMatchLength,
                                 int _k1 = -1, int _k2 = -1, int skipMargin = 0, uint32_t minMatchLength = UINT32_MAX):
            SlidingWindowSparseEMMatcher(refLengthLimit, targetMatchLength, _k1, _k2, skipMargin, minMatchLength, false)
    {};

    void disableSlidingWindow() { swSize = 0; swEnd = circularBuffer ? 0 : maxRefLength; };

    void disableCircularBuffer() { circularBuffer = false; swEnd = maxRefLength; };

    void setSlidingWindowSize(uint8_t factor) { swSize = this->maxRefLength / factor; };

    virtual void loadRef(const char *refText, size_t refLength, bool loadRCRef,
                         bool addRegionSeparators, char regionSeparator);

    void loadSeparator(char regionSeparator);

    const char* getRef() { return start1; };
    size_t getMaxRefLength() { return maxRefLength; };
    size_t getRefLength() { return reachedRefLengthCount ? maxRefLength : pos1; };
    size_t getLoadingPosition() { return pos1; };
    size_t getLoadedRefLength() { return reachedRefLengthCount * (maxRefLength - REF_SHIFT) + (pos1 - REF_SHIFT); };

    void setPosition(size_t refPos, int reachedRefLengthCount) {
        this->pos1 = refPos;
        this->reachedRefLengthCount = reachedRefLengthCount;
    }

    size_t acquireWorkerMatchingLockPos();
    void releaseWorkerMatchingLockPos(size_t lockValue);

    virtual ~SlidingWindowSparseEMMatcher();

    void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsRef, bool revComplMatching,
                    uint32_t minMatchLength, size_t matchingLockPos = SW_END_ERASED_FLAG);

    void matchTexts(vector<TextMatch> &resMatches, const char* destText, size_t destLen, bool destIsRef,
            bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos = SW_END_ERASED_FLAG);

};


class SlidingWindowExpSparseEMMatcher: public SlidingWindowSparseEMMatcher {
protected:
    const int k1ord;
    virtual inline uint64_t htEncodePos(uint64_t pos) const override { return pos >> k1ord; };
    virtual inline uint64_t htDecodePos(uint64_t code) const override { return code << k1ord; };

public:
    SlidingWindowExpSparseEMMatcher(const size_t refLengthLimit, const uint32_t targetMatchLength, int _k1,
                                    int _k2, int skipMargin, uint32_t minMatchLength = UINT32_MAX);
};

class SlidingWindowExpRandSparseEMMatcher: public SlidingWindowExpSparseEMMatcher {
protected:

    const uint32_t POPCOUNT_MASK = 0x00002CB3;

    inline uint64_t htRePos(uint64_t pos) const override {
        return htDecodePos(htEncodePos(pos));
    };

    inline uint64_t htDecodePos(uint64_t code) const override {
        size_t randDelta = __builtin_popcount((uint32_t) code & POPCOUNT_MASK);
        return (code << k1ord) - randDelta;
    };

public:
    SlidingWindowExpRandSparseEMMatcher(const size_t refLengthLimit, const uint32_t targetMatchLength, int _k1,
                                    int _k2, int skipMargin, uint32_t minMatchLength = UINT32_MAX);
};

#endif //PGTOOLS_SWSMEMMATCHER_H

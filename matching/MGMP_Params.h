#ifndef MGMP_PARAMS_H
#define MGMP_PARAMS_H

#include <string>
#include <omp.h>

class MGMP_Params {
public:

    static constexpr char *const STANDARD_IO_POSIX_ALIAS = (char*) "-";

    bool refRegionSeparators = true;

    static const char REF_REGION_SEPARATOR = 0;

    static const int SPEED_MODE_MATCHER_NO_OF_THREADS = 12;
    static const int DEFAULT_MATCHER_NO_OF_THREADS = 8;
    static const int PARALLELIO_MODE_THREADS_LIMIT = 3;
    static const int SINGLEFILE_PARALLEL_MIN_TARGETS = 4;

    bool matcherWorkingThreadsFixed = false;
    int matcherWorkingThreads = DEFAULT_MATCHER_NO_OF_THREADS;

    static const int ADJUSTED_REFERENCE_FACTOR_FLAG = -1;

    static const int DEFAULT_KMER_LENGTH = 32;
    static const int PROTEINS_PROFILE_KMER_LENGTH = 16;

    static const int DEFAULT_MODE_REF_SAMPLING = 16;
    static const int ALT_REF_SAMPLING = 7;

    static const int MIN_MATCH_LENGTH = 16;
    static const int DEFAULT_SKIP_MARGIN = 16;
    static const int MAX_MODE_SKIP_MARGIN = 24;

    static const int DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR = 16;

    static const int DEFAULT_BIG_REFERENCE_COMPRESSOR_RATIO = 16;
    static const int MAX_BIG_REFERENCE_COMPRESSOR_RATIO = 4;

    static const int DEFAULT_UNMATCHED_LENGTH_FACTOR = 128;
    static const int DEFAULT_UNMATCHED_LENGTH_RC_FACTOR = 8;

    static const uint64_t DEFAULT_REF_LITERAL_BEFORE_AFTER_EXT = 32;
    static const uint64_t DEFAULT_REF_LITERAL_MINIMAL_LENGTH_EXT = SIZE_MAX;

    static const uint64_t MIN_REF_INIT_SIZE = 1 << 16;
    static const uint64_t MIN_BASIC_BLOCK_SIZE = 1 << 21;
    static const uint64_t AVG_REF_INIT_SIZE = 1 << 23;
    static const uint64_t SEQ_BLOCK_SIZE = 1 << 17;

    static const uint8_t FRACTION_REF_EXTENSION_STRATEGY = 0;
    static const uint8_t NO_RC_REF_EXTENSION_STRATEGY_MASK = 2;

    static const uint64_t REFERENCE_LENGTH_LIMIT = (uint64_t) UINT32_MAX << 8;

    static const uint8_t DEFAULT_ALLOWED_TARGETS_OUTRUN_FACTOR = 8;
    static const uint8_t REPO_MODE_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS = 0;
    static const uint8_t REPO_MODE_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS = 2;
    static const uint8_t DEFAULT_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS = 1;
    static const uint8_t DEFAULT_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS = 16;
    static const uint8_t SPEED_MODE_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS = 4;
    static const uint8_t SPEED_MODE_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS = 32;
    static const uint64_t DEFAULT_MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS = 1024;

    uint8_t refExtensionStrategy = FRACTION_REF_EXTENSION_STRATEGY;
    uint64_t refLiteralBeforeAfterExt = DEFAULT_REF_LITERAL_BEFORE_AFTER_EXT;
    uint64_t refLiteralMinimalLengthExt = DEFAULT_REF_LITERAL_MINIMAL_LENGTH_EXT;

    int dynamicUnmatchedFractionFactorLimit = DEFAULT_UNMATCHED_LENGTH_FACTOR;
    int currentUnmatchedFractionFactor = DEFAULT_UNMATCHED_LENGTH_FACTOR;
    bool unmatchedFractionFactorFixed = false;
    int unmatchedFractionRCFactor = DEFAULT_UNMATCHED_LENGTH_RC_FACTOR;
    bool unmatchedFractionRCFactorFixed = false;
    bool contigsIndivduallyReversed = true;

    uint8_t allowedTargetsOutrunFactor = DEFAULT_ALLOWED_TARGETS_OUTRUN_FACTOR;
    uint8_t allowedTargetsOutrunForDissimilarContigs = DEFAULT_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS;
    uint64_t minimalLengthForDissimilarContigs = DEFAULT_MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS;
    uint8_t unmatchedFractionFactorTweakForDissimilarContigs = DEFAULT_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS;

    inline bool isRCinReferenceDisabled() const {
        return refExtensionStrategy & NO_RC_REF_EXTENSION_STRATEGY_MASK;
    }

    const static int MIN_PROBE_LEN = 2 << 7;
    const static int MAX_PROBE_LEN = 2 << 15;
    const static int MAX_NON_STANDARD_SYMBOLS_IN_PERCENT = 10;

    int probe_remaining = MAX_PROBE_LEN;
    int probe_non_std_count = 0;

    bool probeProteinsProfile(char* seq, uint64_t len)
    {
        if (probe_remaining == 0)
            return false;
        len = len > probe_remaining ? probe_remaining : len;
        probe_remaining -= len;

        for (int i = 0; i < len; i++)
        {
            switch (seq[i])
            {
            case 'a':
            case 'c':
            case 'g':
            case 't':
            case 'u':
            case 'A':
            case 'C':
            case 'G':
            case 'T':
            case 'U':
            case 'N': break;
            default: probe_non_std_count++;
            }
        }
        int probe_len = MAX_PROBE_LEN - probe_remaining;
        int non_std_symbols_in_percents = probe_non_std_count * 100 / (int) len;
        if (probe_len >= MIN_PROBE_LEN && k != PROTEINS_PROFILE_KMER_LENGTH
                && non_std_symbols_in_percents > MAX_NON_STANDARD_SYMBOLS_IN_PERCENT)
        {
            probe_remaining = 0;
            return true;
        }
        return false;
    }

#ifdef DEVELOPER_BUILD
    static const uint8_t COMBINED_REF_EXTENSION_STRATEGY_MASK = 1;
    static const uint8_t ONLY_LITERAL_REF_EXTENSION_STRATEGY_MASK = 4;
    static const int INITIAL_REF_EXT_GOAL_FACTOR = 8;
    static const int FINAL_REF_EXT_END_MARGIN_FILES = 16;
    static const int MINIMAL_UNMATCHED_LENGTH_FACTOR = 128;
    static const int RELAX_FRACTION_FACTOR_FACTOR = 4;
    static const int TIGHTEN_FRACTION_FACTOR_REDUCTION = 4;
    static const uint8_t BLOCK_REF_EXTENSION_STRATEGY_MASK = 8;
    static const uint8_t DYNAMIC_REF_EXT_FACTOR_MASK = 16;

    inline bool usesCombinedRefExtensionStrategy() const {
        return refExtensionStrategy & COMBINED_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool useOnlyLiteralsinReference() const {
        return refExtensionStrategy & ONLY_LITERAL_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool splitContigsIntoBlocks() const {
        return refExtensionStrategy & BLOCK_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool dynamicRefExtStrategy() const {
        return refExtensionStrategy & DYNAMIC_REF_EXT_FACTOR_MASK;
    }

    void relaxUnmatchedFractionFactor() {
        int tmp = currentUnmatchedFractionFactor;
        if (dynamicRefExtStrategy())
            tmp *= RELAX_FRACTION_FACTOR_FACTOR;
        if (tmp > dynamicUnmatchedFractionFactorLimit)
            tmp = dynamicUnmatchedFractionFactorLimit;
        currentUnmatchedFractionFactor = tmp;
    }

    void tightenUnmatchedFractionFactor() {
        int tmp = currentUnmatchedFractionFactor;
        if (dynamicRefExtStrategy())
            tmp -= TIGHTEN_FRACTION_FACTOR_REDUCTION;
        if (tmp < 1)
            tmp = 1;
        currentUnmatchedFractionFactor = tmp;
    }
#endif

    inline bool isLiteralProperForRefExtension(uint64_t litLength) {
        return litLength > refLiteralMinimalLengthExt;
    }

    inline bool isContigProperForRefExtension(uint64_t seqLength, uint64_t unmatchedChars, int unmatchedFractionFactor) {
        return
#ifdef DEVELOPER_BUILD
                !useOnlyLiteralsinReference() && (!usesCombinedRefExtensionStrategy() ||
                                                  unmatchedChars > MINIMAL_UNMATCHED_LENGTH_FACTOR) &&
                #endif
                unmatchedChars * unmatchedFractionFactor > seqLength;
    }

    inline bool isContigProperForRefRCExtension(uint64_t seqLength, uint64_t unmatchedChars, int unmatchedFractionRCFactor) {
        return unmatchedChars * unmatchedFractionRCFactor > seqLength;
    }


    inline bool isContigDissimilar(uint64_t seqLength, int64_t unmatchedChars, int unmatchedFractionFactor) {
        return seqLength > minimalLengthForDissimilarContigs &&
               unmatchedChars * (unmatchedFractionFactor / unmatchedFractionFactorTweakForDissimilarContigs) > seqLength;
    }

    // INPUT PARAMETERS
    int k = DEFAULT_KMER_LENGTH;
    int k1 = DEFAULT_MODE_REF_SAMPLING;
    bool referenceSamplingStepFixed = false;
    int k2 = 1;
    uint8_t skipMargin = DEFAULT_SKIP_MARGIN;
    bool skipMarginFixed = false;
    int referenceFactor = ADJUSTED_REFERENCE_FACTOR_FLAG;
    int referenceSlidingWindowFactor = DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR;
    bool referenceSlidingWindowFactorOption = false;
    bool enable40bitReference = true;
    bool circularReference = true;
    uint8_t bigReferenceCompressorRatio = DEFAULT_BIG_REFERENCE_COMPRESSOR_RATIO;
    bool separateRCBuffer = false;

    bool allowLossyParsing = false;
    bool enableDNALineLengthDetection = true;

    bool uppercaseDNA = false;

    bool sequentialMatching = false;
    bool g0IsTarget = false;
    bool interleaveFileOrder = false;
    bool bruteParallel = false;

    void setSequentialMatchingMode() {
        sequentialMatching = true;
    }


};

#endif //MGMP_PARAMS_H

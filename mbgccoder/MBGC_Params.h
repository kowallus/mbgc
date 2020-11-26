#ifndef MBGC_MBGCPARAMS_H
#define MBGC_MBGCPARAMS_H

#include <string>
#include "../coders/CodersLib.h"

using namespace std;

class MBGC_Params {
public:
    static const char MBGC_VERSION_MODE = '#';
    static const char MBGC_VERSION_MAJOR = 1;
    static const char MBGC_VERSION_MINOR = 0;
    static const char MBGC_VERSION_REVISION = 2;

    static constexpr char *const MBGC_HEADER = (char*) "MBGC";
    static constexpr char *const TEMPORARY_FILE_SUFFIX = (char*) ".temp";

    static const int DEFAULT_NO_OF_THREADS = 8;

    // MATCHING PARAMS

    static const char MATCH_MARK = '%' + 128;
    static const char SEQ_SEPARATOR_MARK = '"' + 128;
    static const char FILE_SEPARATOR_MARK = ';' + 128;

    static const int ADJUSTED_REFERENCE_FACTOR_FLAG = -1;
    static const int BOOSTED_ADJUSTED_REFERENCE_FACTOR_FLAG = -2;

    static const int DEFAULT_KMER_LENGTH = 32;

    static const int DEFAULT_MODE_REF_SAMPLING = 16;
    static const int MAX_MODE_REF_SAMPLING = 7;

    static const int INITIAL_REF_EXT_GOAL_FACTOR = 8;
    static const int FINAL_REF_EXT_END_MARGIN_FILES = 16;
    static const int RELAX_FRACTION_FACTOR_FACTOR = 4;
    static const int TIGHTEN_FRACTION_FACTOR_REDUCTION = 4;
    static const int DEFAULT_UNMATCHED_LENGTH_FACTOR = 192;
    static const int MINIMAL_UNMATCHED_LENGTH_FACTOR = 128;
    static const int REF_LITERAL_BEFORE_AFTER_EXT = 0;
    static const int REF_LITERAL_MINIMAL_LENGTH_EXT = 64;
    static const size_t SEQ_BLOCK_SIZE = 1 << 17;

    static const uint8_t FRACTION_REF_EXTENSION_STRATEGY = 0;
    static const uint8_t COMBINED_REF_EXTENSION_STRATEGY_MASK = 1;
    static const uint8_t NO_RC_REF_EXTENSION_STRATEGY_MASK = 2;
    static const uint8_t LITERAL_REF_EXTENSION_STRATEGY_MASK = 4;
    static const uint8_t BLOCK_REF_EXTENSION_STRATEGY_MASK = 8;
    static const uint8_t DYNAMIC_REF_EXT_FACTOR_MASK = 16;

    uint8_t refExtensionStrategy = FRACTION_REF_EXTENSION_STRATEGY;

    inline bool usesCombinedRefExtensionStrategy() const {
        return refExtensionStrategy & COMBINED_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool dontUseRCinReference() const {
        return refExtensionStrategy & NO_RC_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool useLiteralsinReference() const {
        return refExtensionStrategy & LITERAL_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool splitContigsIntoBlocks() const {
        return refExtensionStrategy & BLOCK_REF_EXTENSION_STRATEGY_MASK;
    }

    inline bool dynamicRefExtStrategy() const {
        return refExtensionStrategy & DYNAMIC_REF_EXT_FACTOR_MASK;
    }

    uint8_t dynamicUnmatchedFractionFactorLimit = DEFAULT_UNMATCHED_LENGTH_FACTOR;
    uint8_t currentUnmatchedFractionFactor = DEFAULT_UNMATCHED_LENGTH_FACTOR;

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

    inline bool isContigProperForRefExtension(size_t seqLength, size_t unmatchedChars, int unmatchedFractionFactor) {
        return
#ifdef DEVELOPER_BUILD
            !useLiteralsinReference() && (!usesCombinedRefExtensionStrategy() ||
            unmatchedChars > MINIMAL_UNMATCHED_LENGTH_FACTOR) &&
#endif
            unmatchedChars > seqLength / unmatchedFractionFactor;
    }

    inline bool isLiteralProperForRefExtension(size_t seqLength, size_t litLength) {
        return useLiteralsinReference() && litLength > REF_LITERAL_MINIMAL_LENGTH_EXT &&
                (!usesCombinedRefExtensionStrategy() ||
                        litLength > seqLength / MBGC_Params::dynamicUnmatchedFractionFactorLimit);
    }

    // HEADERS TEMPLATES PARAMS

    static const int MINIMAL_PATTERN_LENGTH = 10;
    static const int MAX_PATTERN_SHIFT = 3;

    // INPUT PARAMETERS
    int k = DEFAULT_KMER_LENGTH;
    int k1 = DEFAULT_MODE_REF_SAMPLING;
    bool referenceSamplingStepFixed = false;
    int k2 = 1;
    int referenceFactor = ADJUSTED_REFERENCE_FACTOR_FLAG;

    uint8_t coderLevel = CODER_LEVEL_NORMAL;

    string seqListFileName;
    string archiveFileName;
    string outputPath;

    bool sequentialMatching = false;
    bool interleaveFileOrder = false;
    bool bruteParallel = false;
#ifdef DEVELOPER_BUILD
    bool validationMode = false;
    uint32_t validFilesCount = 0;
    uint32_t invalidFilesCount = 0;

    bool disableDNAformatting = false;
#endif

    void setBoostedReferenceFactorFlag() {
        if (referenceFactor != ADJUSTED_REFERENCE_FACTOR_FLAG) {
            fprintf(stderr, "Boosting reference factor flag doesn't work with o flag "
                            "(i.e., reference factor binary order).\n\n");
            exit(EXIT_FAILURE);
        }
        referenceFactor = BOOSTED_ADJUSTED_REFERENCE_FACTOR_FLAG;
    }

    void setMixedCollectionMode() {
        interleaveFileOrder = true;
    }

    void setKmerLength(int k) {
        if (k < 24 || k > 40) {
            fprintf(stderr, "k - matching kmer length - should be an integer between 24 and 40.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::k = k;
    }

    void setReferenceSamplingStep(int k1) {
        if (k1 <= 0) {
            fprintf(stderr, "s - reference sampling step - should be a positive integer.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::k1 = k1;
        referenceSamplingStepFixed = true;
    }

    void setReferenceFactorBinaryOrder(int o) {
        if (o < 0 || o > 12) {
            fprintf(stderr, "o - reference factor binary order - should be an integer between 0 and 12.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::referenceFactor = ((size_t) 1) << o;
    }

    void setUnmatchedFractionFactor(int u) {
        if (u < 1 || u > 255) {
            fprintf(stderr, "u - unmatched fraction factor - should be an integer between 1 and 255.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::dynamicUnmatchedFractionFactorLimit = u;
        MBGC_Params::currentUnmatchedFractionFactor = u;
    }

    void setArchiveFileName(const string &archiveFileName) {
        MBGC_Params::archiveFileName = archiveFileName;
    }

    void setSeqListFileName(const string &seqListFileName) {
        MBGC_Params::seqListFileName = seqListFileName;
    }

    void setOutputPath(string outputPath) {
        if (outputPath.size() && outputPath.back() != '/')
            outputPath.push_back('/');
        MBGC_Params::outputPath = outputPath;
    }

    void setSequentialMatchingMode() {
        sequentialMatching = true;
    }

    void setBruteParallelMode() {
        bruteParallel = true;
    }

    void setCompressionLevel(int coderLevel) {
        if (coderLevel > CODER_LEVEL_MAX || coderLevel < CODER_LEVEL_FAST) {
            fprintf(stderr, "Generate quality coefficient should be between %d and %d.\n",
                    CODER_LEVEL_FAST, CODER_LEVEL_MAX);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::coderLevel = coderLevel;
        if (coderLevel == CODER_LEVEL_MAX) {
            setSequentialMatchingMode();
            setBoostedReferenceFactorFlag();
            if (!referenceSamplingStepFixed) {
                MBGC_Params::k1 = MAX_MODE_REF_SAMPLING;
            }
        }
    }

};

#endif //MBGC_MBGCPARAMS_H

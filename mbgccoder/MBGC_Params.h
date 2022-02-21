#ifndef MBGC_MBGCPARAMS_H
#define MBGC_MBGCPARAMS_H

#include <string>
#include "../coders/CodersLib.h"

using namespace std;

class MBGC_Params {
public:
    static const char MBGC_VERSION_MODE = '#';
    static const char MBGC_VERSION_MAJOR = 1;
    static const char MBGC_VERSION_MINOR = 2;
    static const char MBGC_VERSION_REVISION = 2;

    static constexpr char *const MBGC_HEADER = (char*) "MBGC";
    static constexpr char *const TEMPORARY_FILE_SUFFIX = (char*) ".temp";

    static constexpr char *const STANDARD_IO_POSIX_ALIAS = (char*) "-";

    static const int FAST_MODE_NO_OF_THREADS = 12;
    static const int DEFAULT_NO_OF_THREADS = 8;
    static const int PARALLELIO_MODE_THREADS_LIMIT = 3;
    bool noOfThreadsFixed = false;

    // MATCHING PARAMS

    static const char MATCH_MARK = '%' + 128;
    static const char SEQ_SEPARATOR_MARK = '"' + 128;
    static const char FILE_SEPARATOR_MARK = ';' + 128;

    static const int ADJUSTED_REFERENCE_FACTOR_FLAG = -1;

    static const int DEFAULT_KMER_LENGTH = 32;

    static const int DEFAULT_MODE_REF_SAMPLING = 16;
    static const int ALT_REF_SAMPLING = 7;

    static const int DEFAULT_SKIP_MARGIN = 16;
    static const int MAX_MODE_SKIP_MARGIN = 24;

    static const int DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR = 16;

    static const int DEFAULT_BIG_REFERENCE_COMPRESSOR_RATIO = 16;
    static const int MAX_BIG_REFERENCE_COMPRESSOR_RATIO = 4;

    static const int DEFAULT_UNMATCHED_LENGTH_FACTOR = 192;
    static const int MINIMAL_UNMATCHED_LENGTH_FACTOR = 128;

    static const size_t MIN_REF_INIT_SIZE = 1 << 21;
    static const size_t AVG_REF_INIT_SIZE = 1 << 23;
    static const size_t SEQ_BLOCK_SIZE = 1 << 17;

    static const uint8_t FRACTION_REF_EXTENSION_STRATEGY = 0;
    static const uint8_t NO_RC_REF_EXTENSION_STRATEGY_MASK = 2;

    static const size_t REFERENCE_LENGTH_LIMIT = (size_t) UINT32_MAX << 8;

    uint8_t refExtensionStrategy = FRACTION_REF_EXTENSION_STRATEGY;

    uint8_t dynamicUnmatchedFractionFactorLimit = DEFAULT_UNMATCHED_LENGTH_FACTOR;
    uint8_t currentUnmatchedFractionFactor = DEFAULT_UNMATCHED_LENGTH_FACTOR;

    inline bool dontUseRCinReference() const {
        return refExtensionStrategy & NO_RC_REF_EXTENSION_STRATEGY_MASK;
    }

#ifdef DEVELOPER_BUILD
    static const uint8_t COMBINED_REF_EXTENSION_STRATEGY_MASK = 1;
    static const uint8_t LITERAL_REF_EXTENSION_STRATEGY_MASK = 4;
    static const int INITIAL_REF_EXT_GOAL_FACTOR = 8;
    static const int FINAL_REF_EXT_END_MARGIN_FILES = 16;
    static const int RELAX_FRACTION_FACTOR_FACTOR = 4;
    static const int TIGHTEN_FRACTION_FACTOR_REDUCTION = 4;
    static const int REF_LITERAL_BEFORE_AFTER_EXT = 0;
    static const int REF_LITERAL_MINIMAL_LENGTH_EXT = 64;
    static const uint8_t BLOCK_REF_EXTENSION_STRATEGY_MASK = 8;
    static const uint8_t DYNAMIC_REF_EXT_FACTOR_MASK = 16;

    inline bool usesCombinedRefExtensionStrategy() const {
        return refExtensionStrategy & COMBINED_REF_EXTENSION_STRATEGY_MASK;
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

    inline bool isLiteralProperForRefExtension(size_t seqLength, size_t litLength) {
        return useLiteralsinReference() && litLength > REF_LITERAL_MINIMAL_LENGTH_EXT &&
                (!usesCombinedRefExtensionStrategy() ||
                        litLength > seqLength / MBGC_Params::dynamicUnmatchedFractionFactorLimit);

    }
#endif

    inline bool isContigProperForRefExtension(size_t seqLength, size_t unmatchedChars, int unmatchedFractionFactor) {
        return
#ifdef DEVELOPER_BUILD
            !useLiteralsinReference() && (!usesCombinedRefExtensionStrategy() ||
            unmatchedChars > MINIMAL_UNMATCHED_LENGTH_FACTOR) &&
#endif
            unmatchedChars > seqLength / unmatchedFractionFactor;
    }

    // HEADERS TEMPLATES PARAMS

    static const int MINIMAL_PATTERN_LENGTH = 10;
    static const int MAX_PATTERN_SHIFT = 3;

    // INPUT PARAMETERS
    int k = DEFAULT_KMER_LENGTH;
    int k1 = DEFAULT_MODE_REF_SAMPLING;
    bool referenceSamplingStepFixed = false;
    int k2 = 1;
    uint8_t skipMargin = DEFAULT_SKIP_MARGIN;
    bool skipMarginFixed = false;
    int referenceFactor = ADJUSTED_REFERENCE_FACTOR_FLAG;
    int referenceSlidingWindowFactor = DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR;
    bool enable40bitReference = true;
    uint8_t bigReferenceCompressorRatio = DEFAULT_BIG_REFERENCE_COMPRESSOR_RATIO;

    uint8_t coderLevel = CODER_LEVEL_NORMAL;
    bool fastDecoder = false;
    bool ultraStreamsCompression = false;

    int64_t dnaLineLength = -1;
    bool enableDNAformatting = false;

    string inputFastaFileName;
    string seqListFileName;
    bool singleFastaFileMode = false;
    string archiveFileName;
    string outputPath;
    string filterPattern;

    bool sequentialMatching = false;
    bool interleaveFileOrder = false;
    bool bruteParallel = false;
#ifdef DEVELOPER_BUILD
    bool validationMode = false;
    uint32_t validFilesCount = 0;
    uint32_t invalidFilesCount = 0;
    bool concatHeadersAndSequencesMode = false;
#endif

    void setInterleaveFileOrder() {
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

    void setSkipMargin(int skipMargin) {
        if (skipMargin < 0 || skipMargin > UINT8_MAX) {
            fprintf(stderr, "m - skip margin - should be an integer between 0 and 255.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::skipMargin = skipMargin;
        skipMarginFixed = true;
    }

    void setReferenceFactorBinaryOrder(int o) {
        if (o < 0 || o > 12) {
            fprintf(stderr, "o - reference factor binary order - should be an integer between 0 and 12.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::referenceFactor = ((size_t) 1) << o;
        MBGC_Params::bigReferenceCompressorRatio = false;
    }

    void setUnmatchedFractionFactor(int u) {
        if (u < 1 || u > 255) {
            fprintf(stderr, "u - unmatched fraction factor - should be an integer between 1 and 255.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::dynamicUnmatchedFractionFactorLimit = u;
        MBGC_Params::currentUnmatchedFractionFactor = u;
    }

    void setReferenceSlidingWindowFactor(int w) {
        if (w < 2 || w > 255) {
            fprintf(stderr, "w - reference sliding window factor - should be an integer between 1 and 255.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::referenceSlidingWindowFactor = w;
    }

    void setArchiveFileName(const string &archiveFileName) {
        MBGC_Params::archiveFileName = archiveFileName;
    }

    void setSeqListFileName(const string &seqListFileName) {
        if (singleFastaFileMode) {
            fprintf(stderr, "ERROR: sequences list file used in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::seqListFileName = seqListFileName;
    }

    void setInputFastaFileName(const string &inputFastaFileName) {
        if (!seqListFileName.empty()) {
            fprintf(stderr, "ERROR: sequences list file used in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        if (bruteParallel) {
            fprintf(stderr, "Brute parallel mode not supported in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        singleFastaFileMode = true;
        MBGC_Params::inputFastaFileName= inputFastaFileName;
    }

    void setOutputPath(string outputPath) {
        if (outputPath.size() && outputPath != STANDARD_IO_POSIX_ALIAS && outputPath.back() != '/')
            outputPath.push_back('/');
        MBGC_Params::outputPath = outputPath;
    }

    void setFilterPattern(string filterPattern) {
        MBGC_Params::filterPattern = filterPattern;
    }

    void setDNALineLength(int dnaLineLength) {
        if (dnaLineLength <= 0) {
            fprintf(stderr, "l - dna line length - should be a positive integer.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::dnaLineLength = dnaLineLength;
        enableDNAformatting = true;
    }

    void setSequentialMatchingMode() {
        sequentialMatching = true;
    }

    void setBruteParallelMode() {
        if (singleFastaFileMode) {
            fprintf(stderr, "Brute parallel mode not supported in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        bruteParallel = true;
    }

    void limit32bitReference() {
        enable40bitReference = false;
    }

    void setCompressionMode(int coderMode) {
        if (coderMode > CODER_LEVEL_MAX || coderMode < 0) {
            fprintf(stderr, "Compression mode should be between %d and %d.\n",
                    0, CODER_LEVEL_MAX);
            exit(EXIT_FAILURE);
        }
        coderLevel = (coderMode >= CODER_LEVEL_MAX ? CODER_LEVEL_MAX : coderMode + 1);
        fastDecoder = !(coderMode & 1);
        if (coderLevel == CODER_LEVEL_MAX) {
            setSequentialMatchingMode();
            MBGC_Params::bigReferenceCompressorRatio = MAX_BIG_REFERENCE_COMPRESSOR_RATIO;
            if (!skipMarginFixed) {
                MBGC_Params::skipMargin = MAX_MODE_SKIP_MARGIN;
            }
        }
        if (coderLevel == CODER_LEVEL_FAST) {
            if (!noOfThreadsFixed) {
                PgHelpers::numberOfThreads = FAST_MODE_NO_OF_THREADS;
            }
        }
    }

    void setUltraStreamsCompression() {
        ultraStreamsCompression = true;
    }
};

#endif //MBGC_MBGCPARAMS_H

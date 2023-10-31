#ifndef MBGC_MBGCPARAMS_H
#define MBGC_MBGCPARAMS_H

#include <string>
#include "../coders/CodersLib.h"
#include <omp.h>

using namespace std;

class MBGC_Params {
public:
    static const char MBGC_VERSION_MODE = '#';
    static const char MBGC_VERSION_MAJOR = 2;
    static const char MBGC_VERSION_MINOR = 0;
    static const char MBGC_VERSION_REVISION = 0;

    char mbgcVersionMajor = MBGC_VERSION_MAJOR;
    char mbgcVersionMinor = MBGC_VERSION_MINOR;
    char mbgcVersionRevision = MBGC_VERSION_REVISION;

    static constexpr char *const MBGC_HEADER = (char*) "MBGC";
    static constexpr char *const TEMPORARY_FILE_SUFFIX = (char*) ".temp";

    static constexpr char *const STANDARD_IO_POSIX_ALIAS = (char*) "-";

    const static uint8_t SPEED_CODER_MODE = 0;
    const static uint8_t DEFAULT_CODER_MODE = 1;
    const static uint8_t REPO_CODER_MODE = 2;
    const static uint8_t MAX_CODER_MODE = 3;

    static const int SPEED_MODE_MATCHER_NO_OF_THREADS = 12;
    static const int DEFAULT_MATCHER_NO_OF_THREADS = 8;
    static const int PARALLELIO_MODE_THREADS_LIMIT = 3;
    bool noOfThreadsLimited = false;
    int coderThreads = -1;
    int backendThreads = -1;
    bool matcherWorkingThreadsFixed = false;
    int matcherWorkingThreads = DEFAULT_MATCHER_NO_OF_THREADS;

    bool appendCommand = false;
    bool infoCommand = false;
    bool repackCommand = false;

    bool storeFileSeparatorMarksInHeadersStream = false;
    bool lazyDecompressionSupport = true;
    bool disableLazyDecompression = false;

    // MATCHING PARAMS

    static const char MATCH_MARK = (char) ('%' + 128);
    static const char SEQ_SEPARATOR_MARK = (char) ('"' + 128);
    static const char FILE_SEPARATOR_MARK = (char) (';' + 128);
    static const char RC_MATCH_MARK = (char) ('$' + 128);
    static const char REF_REGION_SEPARATOR = 0;

    static const int ADJUSTED_REFERENCE_FACTOR_FLAG = -1;

    static const int DEFAULT_KMER_LENGTH = 32;
    static const int PROTEINS_PROFILE_KMER_LENGTH = 16;

    static const int DEFAULT_MODE_REF_SAMPLING = 16;
    static const int ALT_REF_SAMPLING = 7;

    static const int MIN_MATCH_LENGTH = 16;
    static const int DEFAULT_SKIP_MARGIN = 16;
    static const int MAX_MODE_SKIP_MARGIN = 24;

    static const int MAX_GAP_DEPTH = 128;
    static const int DEFAULT_GAP_DEPTH_OFFSET_ENCODING = 64;
    static const int DEFAULT_GAP_DEPTH_MISMATCHES_ENCODING = 2;
    static const int DEFAULT_GAP_BREAKING_MATCH_MIN_LENGTH = 256;

    static const int MAX_EXTEND_MATCH_LEFT_LENGTH = 1 << 24;

    static const int DEFAULT_MISMATCHES_LIMIT_IN_PERCENT = 50;
    static const int DEFAULT_MAX_CONSECUTIVE_MISMATCHES = 10;
    static const int DEFAULT_MISMATCHES_INITIAL_SCORE_IN_PERCENT = 25;

    static const int DEFAULT_RC_MATCH_MINIMUM_LENGTH = 55;

    static const int DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR = 16;

    static const int DEFAULT_BIG_REFERENCE_COMPRESSOR_RATIO = 16;
    static const int MAX_BIG_REFERENCE_COMPRESSOR_RATIO = 4;

    static const int DEFAULT_UNMATCHED_LENGTH_FACTOR = 128;
    static const int DEFAULT_UNMATCHED_LENGTH_RC_FACTOR = 8;
    static const int FAST_LAZY_DECODING_UNMATCHED_LENGTH_FACTOR = 16;
    static const int FAST_LAZY_DECODING_UNMATCHED_LENGTH_RC_FACTOR = 2;
    static const uint64_t DEFAULT_REF_LITERAL_BEFORE_AFTER_EXT = 32;
    static const uint64_t DEFAULT_REF_LITERAL_MINIMAL_LENGTH_EXT = SIZE_MAX;

    static const uint64_t MIN_REF_INIT_SIZE = 1 << 21;
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
    bool referenceSlidingWindowFactorOption = false;
    bool enable40bitReference = true;
    bool circularReference = true;
    uint8_t bigReferenceCompressorRatio = DEFAULT_BIG_REFERENCE_COMPRESSOR_RATIO;
    bool separateRCBuffer = false;
    bool frugal64bitLenEncoding = true;

    uint8_t coderMode = DEFAULT_CODER_MODE;
    bool ultraStreamsCompression = false;
    bool enableExtensionsWithMismatches = true;
    bool mismatchesWithExclusion = true;

    uint8_t gapDepthOffsetEncoding = DEFAULT_GAP_DEPTH_OFFSET_ENCODING;
    uint8_t gapDepthMismatchesEncoding = DEFAULT_GAP_DEPTH_MISMATCHES_ENCODING;
    uint64_t gapBreakingMatchMinLength = DEFAULT_GAP_BREAKING_MATCH_MIN_LENGTH;

    bool isMismatchesScoringParamPresent = false;
    uint8_t mismatchesLimitInPercent = DEFAULT_MISMATCHES_LIMIT_IN_PERCENT;
    uint8_t maxConsecutiveMismatches = DEFAULT_MAX_CONSECUTIVE_MISMATCHES;
    bool maxConsecutiveMismatchesFixed = false;
    uint8_t mismatchesInitialScoreInPercent = DEFAULT_MISMATCHES_INITIAL_SCORE_IN_PERCENT;

    int mmsMatchBonus, mmsMismatchPenalty, mmsMismatchesScoreThreshold, mmsMismatchesInitialScore;

    void initMismatchesMatchingScoreParams() {
        mmsMatchBonus = mismatchesLimitInPercent;
        mmsMismatchPenalty = 100 - mmsMatchBonus;
        mmsMismatchesScoreThreshold = (int) maxConsecutiveMismatches * mmsMismatchPenalty;
        mmsMismatchesInitialScore = mmsMismatchesScoreThreshold * (int) mismatchesInitialScoreInPercent / 100;
    }

    uint64_t rcMatchMinLength = 0;
    bool rcMatchMinLengthFixed = false;
    bool rcRedundancyRemoval = rcMatchMinLength != 0;

    static const bool checkIfDNAisWellFormed = false;
    bool enableDNALineLengthEncoding = true;

    int64_t dnaLineLength = -1;
    bool enableCustomDNAformatting = false;
    bool uppercaseDNA = false;
    bool ignoreFastaFilesPath = false;
    static const int RECOMMENDED_GZ_COMPRESSION_LEVEL = 3;
    int decompressionToGzCoderLevel = 0;
    bool listHeadersMode = false;

    string inputFileName;
    string seqListFileName;
    string inArchiveFileName;
    string outArchiveFileName;
    string outputPath;
    string masterFilterPattern;

    bool sequentialMatching = false;
    bool g0IsTarget = false;
    bool interleaveFileOrder = false;
    bool bruteParallel = false;
#ifdef DEVELOPER_BUILD
    bool validationMode = false;
    static const int VALIDATION_LOG_LIMIT = 100;
    static const int VALIDATION_DUMP_LIMIT = 3;
    static constexpr char *const VALIDATION_DUMP_PATH = (char*) "!validation_dump/";
    bool skipActualValidation = false;
    uint32_t dumpedFilesCount = 0;
    uint32_t validFilesCount = 0;
    uint32_t invalidFilesCount = 0;
    bool concatHeadersAndSequencesMode = false;

    bool logFileEnabled = false;
    bool logVersionCorrect = true;

    virtual ~MBGC_Params() {
        if (logFileEnabled) {
            paramsToLog();
            *PgHelpers::logout << (logVersionCorrect ? "" : "WARNING: wrong log version!") << endl;
            PgHelpers::closeLogFile();
        }
    }

    void enableLogFile(char* filename) {
        ifstream in(filename);
        if (in) {
            int major, minor, revision;
            in >> major >> minor >> revision;
            if (major != MBGC_VERSION_MAJOR || minor != MBGC_VERSION_MINOR || revision != MBGC_VERSION_REVISION) {
                logVersionCorrect = false;
                fprintf(stderr, "WARNING: Log file %s version %d.%d.%d incorrect (expected %d.%d.%d).\n",
                        filename, major, minor, revision, MBGC_VERSION_MAJOR, MBGC_VERSION_MINOR, MBGC_VERSION_REVISION);
            }
        }
        PgHelpers::openLogFile(filename);
        logFileEnabled = true;
        if (!in) {
            generateLogHeader();
        }
    }
#endif

    constexpr static const char * const CLI_COMMANDS[] {
            "compress FASTA file(s) into archive",
            "decompress FASTA file(s) from archive",
            "info about archive contents (FASTA file names & headers)",
            "append FASTA file(s) to the given archive",
            "repack selected FASTA files from existing to new archive"
#ifdef DEVELOPER_BUILD
            , "validate contents of archive"
#endif
    };
    constexpr static const int cliCommandsCount = sizeof(CLI_COMMANDS) / sizeof(CLI_COMMANDS[0]);

    static const char COMPRESS_CMD = CLI_COMMANDS[0][0];
    static const char DECOMPRESS_CMD = CLI_COMMANDS[1][0];
    static const char INFO_CMD = CLI_COMMANDS[2][0];
    static const char APPEND_CMD = CLI_COMMANDS[3][0];
    static const char REPACK_CMD = CLI_COMMANDS[4][0];
#ifdef DEVELOPER_BUILD
    static const char VALIDATE_CMD = CLI_COMMANDS[5][0];
#endif

private:
    constexpr static const char *const COMMON_OPTS = "t:Ivh";
    constexpr static const char *const COMMON_DEV_BUILD_OPTS = "PO:";

    constexpr static const char *const COMPRESS_OPTS = "m:Qk:u:r:s:o:Cg:b:x:XR:";
    constexpr static const char *const COMPRESS_DEV_BUILD_OPTS = "BSM:Dw:G:y:Y:14E:N:8WA:j:J:Vp";

    constexpr static const char *const COMPRESS_APPEND_COMMON_OPTS = "i:";
    constexpr static const char *const ENCODING_OPTS = "T:U";

    constexpr static const char *const COMMON_FILTERING_OPTS = "f:F:";

    constexpr static const char *const DECOMPRESS_OPTS = "z:";
    constexpr static const char *const DECOMPRESS_DEV_BUILD_OPTS = "c";
    constexpr static const char *const DECODING_COMMON_OPTS = "l:L";

    constexpr static const char *const INFO_OPTS = "H";

    constexpr static const char *const REPACK_OPTS = "";

    constexpr static const char *const VALIDATION_DEV_BUILD_OPTS = "DVU";

public:
#ifdef DEVELOPER_BUILD
    inline static const string COMMON_SHORT_OPTS = string(COMMON_OPTS) + string(COMMON_DEV_BUILD_OPTS);

    inline static const string COMPRESS_SHORT_OPTS = string(COMPRESS_APPEND_COMMON_OPTS) + string(COMPRESS_OPTS)
                                                     + string(COMPRESS_DEV_BUILD_OPTS) + string(ENCODING_OPTS) + COMMON_SHORT_OPTS;

    inline static const string DECOMPRESS_SHORT_OPTS = string(DECOMPRESS_OPTS) + string(DECODING_COMMON_OPTS)
                                                       + string(DECOMPRESS_DEV_BUILD_OPTS) + string(COMMON_FILTERING_OPTS) + COMMON_SHORT_OPTS;
    inline static const string VALIDATION_SHORT_OPTS = string(VALIDATION_DEV_BUILD_OPTS) + string(DECODING_COMMON_OPTS)
                                                       + string(COMMON_FILTERING_OPTS) + COMMON_SHORT_OPTS;
#else
    inline static const string COMMON_SHORT_OPTS = string(COMMON_OPTS);

    inline static const string COMPRESS_SHORT_OPTS = string(COMPRESS_APPEND_COMMON_OPTS) + string(MBGC_Params::COMPRESS_OPTS)
                                                    + string(ENCODING_OPTS) + MBGC_Params::COMMON_SHORT_OPTS;
    inline static const string DECOMPRESS_SHORT_OPTS = string(DECOMPRESS_OPTS) + string(DECODING_COMMON_OPTS)
                                            + string(COMMON_FILTERING_OPTS) + COMMON_SHORT_OPTS;
#endif
    inline static const string INFO_SHORT_OPTS = string(INFO_OPTS) + string(COMMON_FILTERING_OPTS)
                                                 + COMMON_SHORT_OPTS;
    inline static const string APPEND_SHORT_OPTS = string(COMPRESS_APPEND_COMMON_OPTS)
                                                   + string(ENCODING_OPTS) + COMMON_SHORT_OPTS;
    inline static const string REPACK_SHORT_OPTS = string(REPACK_OPTS) + string(DECODING_COMMON_OPTS)
                                                   + COMPRESS_SHORT_OPTS + string(COMMON_FILTERING_OPTS) + COMMON_SHORT_OPTS;

    static const char NO_OF_THREADS_OPT = COMMON_OPTS[0];
    static const char IGNORE_FASTA_FILES_PATHS_OPT = COMMON_OPTS[2];
    static const char VERSION_OPT = COMMON_OPTS[3];
    static const char HELP_OPT = COMMON_OPTS[4];

    static const char PRINTOUT_PARAMS_OPT = COMMON_DEV_BUILD_OPTS[0];
    static const char PARAMS_TO_LOGFILE_OPT = COMMON_DEV_BUILD_OPTS[1];

    static const char INPUT_FILE_OPT = COMPRESS_APPEND_COMMON_OPTS[0];
    static const char COMPRESSION_MODE_OPT = COMPRESS_OPTS[0];
    static const char SEQ_MATCHING_MODE_OPT = COMPRESS_OPTS[2];
    static const char KMER_LENGTH_OPT = COMPRESS_OPTS[3];
    static const char UNMATCHED_FRACTION_FACTOR_OPT = COMPRESS_OPTS[5];
    static const char UNMATCHED_FRACTION_RC_FACTOR_OPT = COMPRESS_OPTS[7];
    static const char REF_SAMPLING_STEP_OPT = COMPRESS_OPTS[9];
    static const char REF_FACTOR_BIN_ORDER_OPT = COMPRESS_OPTS[11];
    static const char CIRCULAR_REFERENCE_BUFFER_OPT = COMPRESS_OPTS[13];
    static const char GAP_DEPTH_OFFSET_ENCODING_OPT = COMPRESS_OPTS[14];
    static const char GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT = COMPRESS_OPTS[16];
    static const char MAX_CONSECUTIVE_MISMATCHES_OPT = COMPRESS_OPTS[18];
    static const char DISABLE_MISMATCHES_WITH_EXCLUSION_OPT = COMPRESS_OPTS[20];
    static const char RC_MATCH_MINIMAL_LENGTH_OPT = COMPRESS_OPTS[21];

    static const char MATCHER_NO_OF_THREADS_OPT = ENCODING_OPTS[0];
    static const char UPPERCASE_DNA_OPT = ENCODING_OPTS[2];

    static const char BRUTE_PARALLEL_MODE_OPT = COMPRESS_DEV_BUILD_OPTS[0];
    static const char ULTRA_STREAMS_COMPRESSION_OPT = COMPRESS_DEV_BUILD_OPTS[1];
    static const char SKIP_MARGIN_OPT = COMPRESS_DEV_BUILD_OPTS[2];
    static const char LIMIT_32_BIT_REF_OPT = COMPRESS_DEV_BUILD_OPTS[4];
    static const char REF_SLIDING_WINDOW_FACTOR_OPT = COMPRESS_DEV_BUILD_OPTS[5];
    static const char GAP_DEPTH_MISMATCHES_ENCODING_OPT = COMPRESS_DEV_BUILD_OPTS[7];
    static const char MISMATCHES_LIMIT_IN_PERCENT_OPT = COMPRESS_DEV_BUILD_OPTS[9];
    static const char MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT = COMPRESS_DEV_BUILD_OPTS[11];

    static const char COMBINED_REF_EXTENSION_STRATEGY_OPT = COMPRESS_DEV_BUILD_OPTS[13];
    static const char ONLY_LITERAL_REF_EXTENSION_STRATEGY_OPT = COMPRESS_DEV_BUILD_OPTS[14];
    static const char REF_LITERAL_BEFORE_AFTER_EXT_OPT = COMPRESS_DEV_BUILD_OPTS[15];
    static const char REF_LITERAL_MINIMAL_LENGTH_OPT = COMPRESS_DEV_BUILD_OPTS[17];
    static const char BLOCK_REF_EXTENSION_STRATEGY_OPT = COMPRESS_DEV_BUILD_OPTS[19];
    static const char DYNAMIC_REF_EXTENSION_FACTOR_OPT = COMPRESS_DEV_BUILD_OPTS[20];

    static const char ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT = COMPRESS_DEV_BUILD_OPTS[21];
    static const char MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT = COMPRESS_DEV_BUILD_OPTS[23];
    static const char UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT = COMPRESS_DEV_BUILD_OPTS[25];
    static const char INTERLEAVE_FILE_ORDER_OPT = COMPRESS_DEV_BUILD_OPTS[27];
    static const char PROTEINS_PROFILE_OPT = COMPRESS_DEV_BUILD_OPTS[28];

    static const char FILTER_PATTERN_OPT = COMMON_FILTERING_OPTS[0];
    static const char PATTERNS_FILE_OPT = COMMON_FILTERING_OPTS[2];
    static const char GZ_DECOMPRESSION_OPT = DECOMPRESS_OPTS[0];
    static const char DNA_LINE_LENGTH_OPT = DECODING_COMMON_OPTS[0];
    static const char DISABLE_LAZY_DECODING_OPT = DECODING_COMMON_OPTS[2];

    static const char CONCAT_HEADERS_AND_SEQS_OPT = DECOMPRESS_DEV_BUILD_OPTS[0];
    static const char DUMP_STREAMS_OPT = VALIDATION_DEV_BUILD_OPTS[0];
    static const char SKIP_ACTUAL_VALIDATION_OPT = VALIDATION_DEV_BUILD_OPTS[1];

    static const char LIST_HEADERS_OPT = INFO_OPTS[0];

    void engageAppendCommand() {
        if (!isVersionAtLeast(MBGC_VERSION_MAJOR, MBGC_VERSION_MINOR)) {
            fprintf(stderr, "%s archive version too old (%d.%d < %d.%d). Please, first use 'repack' command.\n",
                    inArchiveFileName.c_str(), (int) mbgcVersionMajor, (int) mbgcVersionMinor,
                    (int) MBGC_VERSION_MAJOR, (int) MBGC_VERSION_MINOR);
            exit(EXIT_FAILURE);
        }
        appendCommand = true;
    }

    bool isStdinMode(bool decoding) {
        return (decoding && inArchiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }

    bool isStdoutMode(bool decoding) {
        return  (decoding && outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS) ||
                (!decoding && outArchiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }

    void printout() {
        fprintf(stderr, "\n* printing internal running parameters");
        fprintf(stderr, "\n  * general");
        fprintf(stderr, "\nversion\t%c\t%c_coder\t%c_back\t%c\tg0t\t%c\t%c\n%d.%d.%d\t%d\t%d\t%d\t%d\t%d\t%ld\t%d",
                NO_OF_THREADS_OPT, MATCHER_NO_OF_THREADS_OPT, NO_OF_THREADS_OPT, COMPRESSION_MODE_OPT, DNA_LINE_LENGTH_OPT, GZ_DECOMPRESSION_OPT,
                mbgcVersionMajor, mbgcVersionMinor, mbgcVersionRevision, PgHelpers::numberOfThreads, matcherWorkingThreads, backendThreads, coderMode, g0IsTarget, dnaLineLength, decompressionToGzCoderLevel);
        fprintf(stderr, "\n  * matcher");
        fprintf(stderr, "\n%c\t%c\t%c\t%c_fixed\t2^%c\t!%c\t%c\t%s\t%s\t%s\n%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d",
                KMER_LENGTH_OPT, REF_SAMPLING_STEP_OPT, SKIP_MARGIN_OPT, SKIP_MARGIN_OPT, REF_FACTOR_BIN_ORDER_OPT, LIMIT_32_BIT_REF_OPT, REF_SLIDING_WINDOW_FACTOR_OPT, "k2", "BRCR", "SepaRC",
                k, k1, skipMargin, skipMarginFixed, referenceFactor, enable40bitReference, referenceSlidingWindowFactor, k2, bigReferenceCompressorRatio, separateRCBuffer);
        fprintf(stderr, "\n  * matching management (+ dev_build reference extension strategies)");
        fprintf(stderr, "\n%c\t%c\t%c\t%c\t%c\t%c\t%c\t%s\t%c\t%c\t%c_limit\n%d\t%d\t%d\t%d\t%d\t%ld\t%d\t%d\t%ld\t%ld\t%d",
                UNMATCHED_FRACTION_FACTOR_OPT, UNMATCHED_FRACTION_RC_FACTOR_OPT, SEQ_MATCHING_MODE_OPT, BRUTE_PARALLEL_MODE_OPT, ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT, MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT, UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT, "ReExSt", REF_LITERAL_MINIMAL_LENGTH_OPT, REF_LITERAL_BEFORE_AFTER_EXT_OPT, DYNAMIC_REF_EXTENSION_FACTOR_OPT,
                currentUnmatchedFractionFactor, unmatchedFractionRCFactor, sequentialMatching, bruteParallel, allowedTargetsOutrunForDissimilarContigs, minimalLengthForDissimilarContigs, unmatchedFractionFactorTweakForDissimilarContigs, refExtensionStrategy, refLiteralMinimalLengthExt, refLiteralBeforeAfterExt, dynamicUnmatchedFractionFactorLimit);
        fprintf(stderr, "\n  * extending matches, mismatches & gaps encoding");
        fprintf(stderr, "\n!%c\t%c\t%c\t%c\t%c\t%c\t%c\n%d\t%d\t%d\t%d\t%d\t%d\t%ld",
                DISABLE_MISMATCHES_WITH_EXCLUSION_OPT, GAP_DEPTH_OFFSET_ENCODING_OPT, GAP_DEPTH_MISMATCHES_ENCODING_OPT, MISMATCHES_LIMIT_IN_PERCENT_OPT, MAX_CONSECUTIVE_MISMATCHES_OPT, MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT, GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT,
                mismatchesWithExclusion, gapDepthOffsetEncoding, gapDepthMismatchesEncoding, mismatchesLimitInPercent, maxConsecutiveMismatches, mismatchesInitialScoreInPercent, gapBreakingMatchMinLength);
        fprintf(stderr, "\n  * miscellaneous");
        fprintf(stderr, "\nrcRR\tFr64LE\t%c\tfsmh\tlazy\n%d\t%d\t%ld\t%d\t%d", RC_MATCH_MINIMAL_LENGTH_OPT,
                rcRedundancyRemoval, frugal64bitLenEncoding, rcMatchMinLength, storeFileSeparatorMarksInHeadersStream,
                lazyDecompressionSupport);
#ifdef DEVELOPER_BUILD
        fprintf(stderr, "\n  * developer build");
        fprintf(stderr, "\n%c\t%c\t%c\n%d\t%d\t%d",
                ULTRA_STREAMS_COMPRESSION_OPT, INTERLEAVE_FILE_ORDER_OPT, CONCAT_HEADERS_AND_SEQS_OPT,
                ultraStreamsCompression, interleaveFileOrder, concatHeadersAndSequencesMode);
#endif
    }

    void generateLogHeader() {
        *PgHelpers::logout << (int) MBGC_VERSION_MAJOR << "\t" << (int) MBGC_VERSION_MINOR << "\t" << (int) MBGC_VERSION_REVISION << endl;
        *PgHelpers::logout << "t_match [ms]\tt_back [ms]\tt_total [ms]\tratio\tinput size [B]\tarchive [B]\t";
        *PgHelpers::logout << "version\t" << NO_OF_THREADS_OPT << "\t" << MATCHER_NO_OF_THREADS_OPT << "_coder\t" <<
                           NO_OF_THREADS_OPT << "_back\t" << COMPRESSION_MODE_OPT << "\t" << "g0t" << "\t" << DNA_LINE_LENGTH_OPT << "\t" << GZ_DECOMPRESSION_OPT << "\t";
        *PgHelpers::logout << KMER_LENGTH_OPT << "\t" << REF_SAMPLING_STEP_OPT << "\t" << SKIP_MARGIN_OPT << "\t" <<
                           SKIP_MARGIN_OPT << "_fixed\t2^" << REF_FACTOR_BIN_ORDER_OPT << "\t!" << LIMIT_32_BIT_REF_OPT << "\t" << CIRCULAR_REFERENCE_BUFFER_OPT << "\t" <<
                           REF_SLIDING_WINDOW_FACTOR_OPT << "\t" << "k2" << "\t" << "BRCR" << "\t" <<  "SepaRC" << "\t";
        *PgHelpers::logout << UNMATCHED_FRACTION_FACTOR_OPT << "\t" << UNMATCHED_FRACTION_RC_FACTOR_OPT << "\t" << SEQ_MATCHING_MODE_OPT << "\t" <<
                           BRUTE_PARALLEL_MODE_OPT << "\t" << ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT << "\t" <<
                           MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT << "\t" <<
                           UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT << "\t" << "ReExSt" << "\t" <<
                           REF_LITERAL_MINIMAL_LENGTH_OPT << "\t" << REF_LITERAL_BEFORE_AFTER_EXT_OPT << "\t" <<
                           DYNAMIC_REF_EXTENSION_FACTOR_OPT << "_limit\t";
        *PgHelpers::logout << "!" << DISABLE_MISMATCHES_WITH_EXCLUSION_OPT << "\t" << GAP_DEPTH_OFFSET_ENCODING_OPT <<
                           "\t" << GAP_DEPTH_MISMATCHES_ENCODING_OPT << "\t" << MISMATCHES_LIMIT_IN_PERCENT_OPT << "\t" <<
                           MAX_CONSECUTIVE_MISMATCHES_OPT << "\t" << MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT << "\t" <<
                           GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT << "\t";
        *PgHelpers::logout << "rcRR\tFr64LE\t" << RC_MATCH_MINIMAL_LENGTH_OPT << "\tfsmh\tlazy\t";
        *PgHelpers::logout << ULTRA_STREAMS_COMPRESSION_OPT << "\t" << INTERLEAVE_FILE_ORDER_OPT << "\t" <<
                           CONCAT_HEADERS_AND_SEQS_OPT << "\t";
        *PgHelpers::logout << endl;
    }

    void paramsToLog() {
        *PgHelpers::logout << (int) mbgcVersionMajor << "." << (int) mbgcVersionMinor << "." <<
                           (int) mbgcVersionRevision << "\t" << PgHelpers::numberOfThreads << "\t" << matcherWorkingThreads << "\t" <<
                           backendThreads << "\t" << (int) coderMode << "\t" << g0IsTarget << "\t" << dnaLineLength << "\t" << decompressionToGzCoderLevel << "\t";
        *PgHelpers::logout << k << "\t" << k1 << "\t" << (int) skipMargin << "\t" << skipMarginFixed << "\t" <<
                           referenceFactor << "\t" << enable40bitReference << "\t" << circularReference << "\t" << referenceSlidingWindowFactor << "\t" <<
                           k2 << "\t" << (int) bigReferenceCompressorRatio << "\t" << separateRCBuffer << "\t";
        *PgHelpers::logout << (int) currentUnmatchedFractionFactor << "\t" << (int) unmatchedFractionRCFactor << "\t" << sequentialMatching << "\t" <<
                           bruteParallel << "\t" << (int) allowedTargetsOutrunForDissimilarContigs << "\t" <<
                           minimalLengthForDissimilarContigs << "\t" << (int) unmatchedFractionFactorTweakForDissimilarContigs <<
                           "\t" <<  (int) refExtensionStrategy << "\t" << (int) refLiteralMinimalLengthExt << "\t" <<
                           refLiteralBeforeAfterExt << "\t" << (int) dynamicUnmatchedFractionFactorLimit << "\t";
        *PgHelpers::logout << mismatchesWithExclusion << "\t" << (int) gapDepthOffsetEncoding << "\t" <<
                           (int) gapDepthMismatchesEncoding << "\t" << (int) mismatchesLimitInPercent << "\t" <<
                           (int) maxConsecutiveMismatches << "\t" << (int) mismatchesInitialScoreInPercent << "\t" <<
                           gapBreakingMatchMinLength << "\t";
        *PgHelpers::logout << rcRedundancyRemoval << "\t" << frugal64bitLenEncoding << "\t" << rcMatchMinLength << "\t" <<
                           storeFileSeparatorMarksInHeadersStream << "\t" << lazyDecompressionSupport << "\t";
#ifdef DEVELOPER_BUILD
        *PgHelpers::logout << ultraStreamsCompression << "\t" << interleaveFileOrder << "\t" <<
                           concatHeadersAndSequencesMode << "\t";
#endif
    }

    bool isVersionAtLeast(char major, char minor) {
        return (mbgcVersionMajor > major ||
                (mbgcVersionMajor == major && mbgcVersionMinor >= minor));
    }

    bool isLazyDecompressionEnabled() {
        return lazyDecompressionSupport && !disableLazyDecompression &&
                (PgHelpers::numberOfThreads > 1 || masterFilterPattern.empty());
    }

    void write(ostream &out) const {
        out.write(MBGC_Params::MBGC_HEADER, strlen(MBGC_Params::MBGC_HEADER));
        out.put(MBGC_Params::MBGC_VERSION_MODE);
        out.put(MBGC_Params::MBGC_VERSION_MAJOR);
        out.put(MBGC_Params::MBGC_VERSION_MINOR);
        out.put(MBGC_Params::MBGC_VERSION_REVISION);
        out.put(coderMode);
        out.put(matcherWorkingThreads);
        out.put(backendThreads);
        out.put(appendCommand ? g0IsTarget : sequentialMatching);
        out.put(k);
        PgHelpers::writeValue(out, k1);
        PgHelpers::writeValue(out, k2);
        PgHelpers::writeValue(out, skipMargin);
        PgHelpers::writeValue(out, referenceFactor);
        out.put(refExtensionStrategy);
        PgHelpers::writeValue(out, refLiteralMinimalLengthExt);
        if (refLiteralMinimalLengthExt != SIZE_MAX)
            PgHelpers::writeValue(out, refLiteralBeforeAfterExt);
        out.put(frugal64bitLenEncoding);
#ifdef DEVELOPER_BUILD
        if (usesCombinedRefExtensionStrategy())
            PgHelpers::writeValue(out, useOnlyLiteralsinReference() ?
                                       dynamicUnmatchedFractionFactorLimit : MBGC_Params::MINIMAL_UNMATCHED_LENGTH_FACTOR);
        if (splitContigsIntoBlocks())
            PgHelpers::writeValue(out, MBGC_Params::SEQ_BLOCK_SIZE);
#endif
        PgHelpers::writeValue(out, bigReferenceCompressorRatio);
        PgHelpers::writeValue(out, allowedTargetsOutrunForDissimilarContigs);
        PgHelpers::writeValue(out, minimalLengthForDissimilarContigs);
        PgHelpers::writeValue(out, unmatchedFractionFactorTweakForDissimilarContigs);
        PgHelpers::writeValue(out, gapDepthOffsetEncoding);
        out.put(enableExtensionsWithMismatches);
        if (enableExtensionsWithMismatches) {
            PgHelpers::writeValue(out, gapDepthMismatchesEncoding);
            PgHelpers::writeValue(out, mismatchesLimitInPercent);
            PgHelpers::writeValue(out, maxConsecutiveMismatches);
            PgHelpers::writeValue(out, mismatchesInitialScoreInPercent);
            PgHelpers::writeValue(out, gapBreakingMatchMinLength);
            out.put(mismatchesWithExclusion);
        }
        PgHelpers::writeValue(out, rcRedundancyRemoval ? rcMatchMinLength : (uint64_t) 0);
        out.put(storeFileSeparatorMarksInHeadersStream);
        out.put(lazyDecompressionSupport);
    }

    void read(istream &in) {
        for (int i = 0; i < strlen(MBGC_Params::MBGC_HEADER); i++) {
            if (MBGC_Params::MBGC_HEADER[i] != in.get()) {
                fprintf(stderr, "Error processing header.\n");
                exit(EXIT_FAILURE);
            }
        }
        char ch = in.get();
        if (ch != MBGC_Params::MBGC_VERSION_MODE) {
            fprintf(stderr, "Error processing header.\n");
            exit(EXIT_FAILURE);
        }
        mbgcVersionMajor = in.get();
        mbgcVersionMinor = in.get();
        mbgcVersionRevision = in.get();
        if (isVersionAtLeast(MBGC_Params::MBGC_VERSION_MAJOR, MBGC_Params::MBGC_VERSION_MINOR + 1)) {
            fprintf(stderr, "Exiting. Unsupported newer version of MBGC archive (%d.%d.%d)\n",
                    (int) mbgcVersionMajor, (int) mbgcVersionMinor,
                    (int) mbgcVersionRevision);
            exit(EXIT_FAILURE);
        }
        coderMode = in.get();
        matcherWorkingThreads = in.get();
        if (isVersionAtLeast(2, 0))
            backendThreads = in.get();
        g0IsTarget = (bool) in.get();
        if (!isVersionAtLeast(1, 2)) {
            sequentialMatching = g0IsTarget;
            g0IsTarget = false;
        }
        k = in.get();
        PgHelpers::readValue<int>(in, k1);
        PgHelpers::readValue<int>(in, k2);
        if (isVersionAtLeast(1, 1))
            PgHelpers::readValue<uint8_t>(in, skipMargin);
        PgHelpers::readValue<int>(in, referenceFactor);
        refExtensionStrategy = in.get();
#ifdef DEVELOPER_BUILD
        if (refExtensionStrategy & !MBGC_Params::NO_RC_REF_EXTENSION_STRATEGY_MASK) {
#else
            if (refExtensionStrategy) {
#endif
            fprintf(stderr, "Error: unsupported reference extension strategy (code: %d).\n",
                    (int) refExtensionStrategy);
            exit(EXIT_FAILURE);
        }
        if (isVersionAtLeast(2, 0)) {
            PgHelpers::readValue<uint64_t>(in, refLiteralMinimalLengthExt);
            if (refLiteralMinimalLengthExt != SIZE_MAX) {
#ifdef DEVELOPER_BUILD
                PgHelpers::readValue<uint64_t >(in, refLiteralBeforeAfterExt);
#else
                fprintf(stderr, "Using literals ref extension in decompression only enabled in developer mode.\n");
            exit(EXIT_FAILURE);
#endif
            }
        } else
            refLiteralMinimalLengthExt = SIZE_MAX;
        if (isVersionAtLeast(2, 0))
            frugal64bitLenEncoding = (bool) in.get();
        else
            frugal64bitLenEncoding = false;
#ifdef DEVELOPER_BUILD
        if (useOnlyLiteralsinReference()) {
            fprintf(stderr, "Literals ref extension strategy decompression not implemented.\n");
            exit(EXIT_FAILURE);
        }
        if (splitContigsIntoBlocks()) {
            fprintf(stderr, "Split contigs into blocks strategy decompression not implemented.\n");
            exit(EXIT_FAILURE);
        }
        int tmp;
        if (usesCombinedRefExtensionStrategy()) {
            PgHelpers::readValue<int>(in, tmp);
            if (tmp != MBGC_Params::MINIMAL_UNMATCHED_LENGTH_FACTOR) {
                fprintf(stderr, "Invalid minimal unmatched length factor: %d (expected %d).\n", tmp,
                        MBGC_Params::MINIMAL_UNMATCHED_LENGTH_FACTOR);
                exit(EXIT_FAILURE);
            }
        }
#endif
        if (isVersionAtLeast(2, 0)) {
            PgHelpers::readValue<uint8_t>(in, bigReferenceCompressorRatio);
            PgHelpers::readValue<uint8_t>(in, allowedTargetsOutrunForDissimilarContigs);
            PgHelpers::readValue<uint64_t>(in, minimalLengthForDissimilarContigs);
            PgHelpers::readValue<uint8_t>(in, unmatchedFractionFactorTweakForDissimilarContigs);
            PgHelpers::readValue<uint8_t>(in, gapDepthOffsetEncoding);
            enableExtensionsWithMismatches = (bool) in.get();
            if (enableExtensionsWithMismatches) {
                PgHelpers::readValue<uint8_t>(in, gapDepthMismatchesEncoding);
                PgHelpers::readValue<uint8_t>(in, mismatchesLimitInPercent);
                PgHelpers::readValue<uint8_t>(in, maxConsecutiveMismatches);
                PgHelpers::readValue<uint8_t>(in, mismatchesInitialScoreInPercent);
                initMismatchesMatchingScoreParams();
                PgHelpers::readValue<uint64_t>(in, gapBreakingMatchMinLength);
                mismatchesWithExclusion = (bool) in.get();
            }
            PgHelpers::readValue<uint64_t>(in, rcMatchMinLength);
            rcRedundancyRemoval = rcMatchMinLength != 0;
            storeFileSeparatorMarksInHeadersStream = (bool) in.get();
            lazyDecompressionSupport = (bool) in.get();
        } else {
            disableGapMismatches();
            disableRCRedundancyRemoval();
            storeFileSeparatorMarksInHeadersStream = true;
            lazyDecompressionSupport = false;
            PgHelpers::revertToMBGC121();
        }
    }

    void disableGapMismatches() {
        enableExtensionsWithMismatches = false;
        gapDepthOffsetEncoding = 0;
        gapDepthMismatchesEncoding = 0;
        gapBreakingMatchMinLength = 0;
        mismatchesInitialScoreInPercent = 0;
        mismatchesLimitInPercent = 0;
        maxConsecutiveMismatches = 0;
    }

    void disableRCRedundancyRemoval() {
        rcRedundancyRemoval = false;
        rcMatchMinLength = 0;
    }

    void setInterleaveFileOrder() {
        if (repackCommand) {
            fprintf(stderr, "Interleaving files order option %c not available in repack command.\n\n",
                    INTERLEAVE_FILE_ORDER_OPT);
            return;
        }
        interleaveFileOrder = true;
    }

    void setKmerLength(int k) {
        if (k < MIN_MATCH_LENGTH || k > 40) {
            fprintf(stderr, "%c - matching kmer length - should be an integer between %d and 40.\n\n",
                    KMER_LENGTH_OPT, MIN_MATCH_LENGTH);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::k = k;
    }

    void setReferenceSamplingStep(int k1) {
        if (k1 <= 0) {
            fprintf(stderr, "%c - reference sampling step - should be a positive integer.\n\n",
                    REF_SAMPLING_STEP_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::k1 = k1;
        referenceSamplingStepFixed = true;
    }

    void setSkipMargin(int skipMargin) {
        if (skipMargin < 0 || skipMargin > UINT8_MAX) {
            fprintf(stderr, "%c - skip margin - should be an integer between 0 and 255.\n\n",
                    SKIP_MARGIN_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::skipMargin = skipMargin;
        skipMarginFixed = true;
    }

    void setReferenceFactorBinaryOrder(int o) {
        if (o < 0 || o > 12) {
            fprintf(stderr, "%c - reference factor binary order - should be an integer between 0 and 12.\n\n",
                    REF_FACTOR_BIN_ORDER_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::referenceFactor = ((uint64_t) 1) << o;
        MBGC_Params::bigReferenceCompressorRatio = 1;
    }

    void setUnmatchedFractionFactor(int u) {
        if (u > 255) {
#ifdef DEVELOPER_BUILD
            fprintf(stderr, "WARNING: unmatched fraction factor above 255 (experimental, decompression not supported).\n");
#else
            fprintf(stderr, "u - unmatched fraction factor - should be an integer between 1 and 255.\n\n");
            exit(EXIT_FAILURE);
#endif
        }
        if (u < 1) {
            fprintf(stderr, "%c - unmatched fraction factor - should be an integer between 1 and 255.\n\n",
                    UNMATCHED_FRACTION_FACTOR_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::dynamicUnmatchedFractionFactorLimit = u;
        MBGC_Params::currentUnmatchedFractionFactor = u;
        MBGC_Params::unmatchedFractionFactorFixed = true;
    }

    void setUnmatchedFractionRCFactor(int r) {
        if (r < 0 || r > 255) {
            fprintf(stderr, "%c - unmatched fraction RC factor - should be an integer between 0 and 255.\n\n",
                    UNMATCHED_FRACTION_RC_FACTOR_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::unmatchedFractionRCFactor = r;
        MBGC_Params::unmatchedFractionRCFactorFixed = true;
        if (r == 0)
            refExtensionStrategy |= MBGC_Params::NO_RC_REF_EXTENSION_STRATEGY_MASK;
    }

    void setGapDepthOffsetEncoding(int g_depth) {
        if (g_depth < 0 || g_depth > MAX_GAP_DEPTH) {
            fprintf(stderr, "%c - gap depth offset encoding - should be an integer between 0 and %d.\n\n", GAP_DEPTH_OFFSET_ENCODING_OPT, MAX_GAP_DEPTH);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::gapDepthOffsetEncoding = g_depth;
    }

    void setGapDepthMismatchesEncoding(int G_depth) {
        if (G_depth < 0 || G_depth > 255) {
            fprintf(stderr, "%c - gap depth mismatches encoding - should be an integer between 0 and %d.\n\n", GAP_DEPTH_MISMATCHES_ENCODING_OPT, MAX_GAP_DEPTH);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::gapDepthMismatchesEncoding = G_depth;
    }

    void setGapBreakingMatchMinimalLength(int a) {
        if (a < 0) {
            fprintf(stderr, "%c - gap breaking match minimal length - should be a non-negative integer.\n\n", GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::gapBreakingMatchMinLength = a;
    }

    void setMismatchesLimitInPercent(int y) {
        if (MBGC_Params::maxConsecutiveMismatches < 2) {
            fprintf(stderr, "%c - mismatches limit in percent - is disabled if %c (max consecutive mismatches) is below 2.\n\n",
                    MISMATCHES_LIMIT_IN_PERCENT_OPT, MAX_CONSECUTIVE_MISMATCHES_OPT);
            exit(EXIT_FAILURE);
        }
        if (y <= 0 || y >= 100) {
            fprintf(stderr, "%c - mismatches limit in percent - should be an integer between 1 and 99.\n\n",
                    MISMATCHES_LIMIT_IN_PERCENT_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::mismatchesLimitInPercent = y;
        MBGC_Params::isMismatchesScoringParamPresent = true;
    }

    void setMaxConsecutiveMismatches(int x) {
        if (x < 0 || x > 255) {
            fprintf(stderr, "%c - max consecutive mismatches - should be an integer between 0 and 255.\n\n",
                    MAX_CONSECUTIVE_MISMATCHES_OPT);
            exit(EXIT_FAILURE);
        }
        if (!MBGC_Params::mismatchesWithExclusion && x == 0) {
            fprintf(stderr, "encoding mismatches is disabled - %c option cannot be used.\n\n",
                    DISABLE_MISMATCHES_WITH_EXCLUSION_OPT);
            exit(EXIT_FAILURE);
        }
        if (MBGC_Params::isMismatchesScoringParamPresent && x < 2) {
            fprintf(stderr, "%c - max consecutive mismatches - cannot be set below 2 is other mismatches scoring params are used.\n\n",
                    MAX_CONSECUTIVE_MISMATCHES_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::maxConsecutiveMismatches = x;
        MBGC_Params::enableExtensionsWithMismatches = x > 0;
        if (x == 1)
            MBGC_Params::mismatchesLimitInPercent = 100;
        MBGC_Params::maxConsecutiveMismatchesFixed = true;
    }

    void setMismatchesInitialScoreInPercent(int Y) {
        if (MBGC_Params::maxConsecutiveMismatches < 2) {
            fprintf(stderr, "%c - mismatches initial score in percent - is disabled if %c (max consecutive mismatches) is below 2.\n\n",
                    MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT, MAX_CONSECUTIVE_MISMATCHES_OPT);
            exit(EXIT_FAILURE);
        }
        if (Y < 0 || Y > 99) {
            fprintf(stderr, "%c - mismatches initial score in percent - should be an integer between 0 and 99.\n\n",
                    MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::mismatchesInitialScoreInPercent = Y;
        MBGC_Params::isMismatchesScoringParamPresent = true;
    }

    void disableMismatchesEncodingWithExclusion() {
        if (!enableExtensionsWithMismatches) {
            fprintf(stderr, "encoding mismatches is disabled - %c option cannot be used.\n\n",
                    DISABLE_MISMATCHES_WITH_EXCLUSION_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::mismatchesWithExclusion = false;
    }

    void disableCircularReference() {
        if (referenceSlidingWindowFactorOption) {
            fprintf(stderr, "%c - reference sliding window factor - "
                            "cannot be used with disabled circular reference buffer (option %c).\n\n",
                    REF_SLIDING_WINDOW_FACTOR_OPT, CIRCULAR_REFERENCE_BUFFER_OPT);
            exit(EXIT_FAILURE);
        }
        circularReference = false;
    }

    void setReferenceSlidingWindowFactor(int w) {
        if (w < 2 || w > 255) {
            fprintf(stderr, "%c - reference sliding window factor - should be an integer between 1 and 255.\n\n",
                    REF_SLIDING_WINDOW_FACTOR_OPT);
            exit(EXIT_FAILURE);
        }
        if (!circularReference) {
            fprintf(stderr, "%c - reference sliding window factor - "
                            "cannot be used with disabled circular reference buffer (option %c).\n\n",
                    REF_SLIDING_WINDOW_FACTOR_OPT, CIRCULAR_REFERENCE_BUFFER_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::referenceSlidingWindowFactor = w;
        referenceSlidingWindowFactorOption = true;
    }

    void setRCMatchMinimalLength(int param_r) {
        if (param_r < 24 && param_r != 0) {
            fprintf(stderr, "%c - rc match minimal length - should be an integer larger than 23 (or 0 - disable).\n\n",
                    RC_MATCH_MINIMAL_LENGTH_OPT);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::rcMatchMinLength = param_r;
        MBGC_Params::rcRedundancyRemoval = param_r != 0;
        MBGC_Params::rcMatchMinLengthFixed = true;
    }

    void setInArchiveFileName(const string &archiveFileName) {
        MBGC_Params::inArchiveFileName = archiveFileName;
    }

    void setOutArchiveFileName(const string &archiveFileName) {
        MBGC_Params::outArchiveFileName = archiveFileName;
    }

    void setSeqListFileName(const string &seqListFileName) {
        if (!inputFileName.empty()) {
            fprintf(stderr, "ERROR: sequences list file used in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::seqListFileName = seqListFileName;
    }

    void setInputFileName(const string &inputFileName) {
        if (!seqListFileName.empty()) {
            fprintf(stderr, "ERROR: sequences list file used in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        if (bruteParallel) {
            fprintf(stderr, "Brute parallel mode not supported in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::inputFileName= inputFileName;
    }

    void setOutputPath(string outputPath) {
        if (outputPath.size() && outputPath != STANDARD_IO_POSIX_ALIAS && outputPath.back() != '/')
            outputPath.push_back('/');
        MBGC_Params::outputPath = outputPath;
    }

    void setFilterPattern(string filterPattern) {
        MBGC_Params::masterFilterPattern = filterPattern;
    }

    void setDNALineLength(int64_t dnaLineLength) {
        if (dnaLineLength < 0) {
            fprintf(stderr, "l - dna line length - should be a non-negative integer.\n\n");
            exit(EXIT_FAILURE);
        }
        MBGC_Params::dnaLineLength = dnaLineLength;
        enableCustomDNAformatting = true;
    }

    void setSequentialMatchingMode() {
        sequentialMatching = true;
    }

    void setBruteParallelMode() {
        if (!inputFileName.empty()) {
            fprintf(stderr, "Brute parallel mode not supported in single fasta file mode.\n\n");
            exit(EXIT_FAILURE);
        }
        bruteParallel = true;
    }

    void limit32bitReference() {
        enable40bitReference = false;
    }

    void setThreadsLimit(int t) {
        if (t < 1) {
            fprintf(stderr, "%c - threads number - should be a positive integer.\n\n", NO_OF_THREADS_OPT);
            exit(EXIT_FAILURE);
        }
        if (matcherWorkingThreadsFixed && t < MBGC_Params::matcherWorkingThreads) {
            fprintf(stderr, "%c - threads number (%d) - cannot smaller than number of worker threads (%c = %d).\n\n",
                    NO_OF_THREADS_OPT, t, MATCHER_NO_OF_THREADS_OPT, matcherWorkingThreads);
            exit(EXIT_FAILURE);
        }
        noOfThreadsLimited = true;
        MBGC_Params::coderThreads = t;
        MBGC_Params::backendThreads = t;
        if (MBGC_Params::matcherWorkingThreads > t)
            MBGC_Params::matcherWorkingThreads = t;
    }

    void setMatcherThreads(int T) {
        if (T < 1) {
            fprintf(stderr, "%c - threads number - should be a positive integer.\n\n", NO_OF_THREADS_OPT);
            exit(EXIT_FAILURE);
        }
        if (noOfThreadsLimited && MBGC_Params::coderThreads < T) {
            fprintf(stderr, "%c - threads number (%d) - cannot smaller than number of worker threads (%c = %d).\n\n",
                    NO_OF_THREADS_OPT, coderThreads, MATCHER_NO_OF_THREADS_OPT, T);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::matcherWorkingThreads = T;
        MBGC_Params::matcherWorkingThreadsFixed = true;
    }

    void enableOmpThreads(int t) {
        if (noOfThreadsLimited && t != MBGC_Params::coderThreads) {
            fprintf(stderr, "ERROR: incorrect number threads enabled: %d (expected %d).\n\n", t, MBGC_Params::coderThreads);
            exit(EXIT_FAILURE);
        }
        omp_set_num_threads(t);
        PgHelpers::numberOfThreads = t;
#pragma omp parallel
#pragma omp single
        *PgHelpers::devout << "number of threads: " << omp_get_num_threads() << endl;
    }

    void setCompressionMode(int coderMode) {
        if (coderMode > MAX_CODER_MODE || coderMode < SPEED_CODER_MODE) {
            fprintf(stderr, "Compression mode should be between %d and %d.\n",
                    SPEED_CODER_MODE, MAX_CODER_MODE);
            exit(EXIT_FAILURE);
        }
        MBGC_Params::coderMode = coderMode;
        if (coderMode == SPEED_CODER_MODE) {
            if (!matcherWorkingThreadsFixed)
                matcherWorkingThreads = SPEED_MODE_MATCHER_NO_OF_THREADS;
            if (noOfThreadsLimited && matcherWorkingThreads > coderThreads)
                matcherWorkingThreads = coderThreads;
            allowedTargetsOutrunForDissimilarContigs = SPEED_MODE_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS;
            unmatchedFractionFactorTweakForDissimilarContigs = SPEED_MODE_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS;
            frugal64bitLenEncoding = false;
        }
        if (coderMode == REPO_CODER_MODE) {
            allowedTargetsOutrunForDissimilarContigs = REPO_MODE_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS;
            unmatchedFractionFactorTweakForDissimilarContigs = REPO_MODE_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS;
        }
        if (coderMode == REPO_CODER_MODE || coderMode == MAX_CODER_MODE) {
            MBGC_Params::bigReferenceCompressorRatio = MAX_BIG_REFERENCE_COMPRESSOR_RATIO;
            if (!skipMarginFixed)
                MBGC_Params::skipMargin = MAX_MODE_SKIP_MARGIN;
            if (!unmatchedFractionRCFactorFixed) {
                MBGC_Params::setUnmatchedFractionRCFactor(DEFAULT_UNMATCHED_LENGTH_FACTOR);
                unmatchedFractionRCFactorFixed = false;
            }
        }
        if (coderMode == MAX_CODER_MODE) {
            setSequentialMatchingMode();
            if (!rcMatchMinLengthFixed) {
                setRCMatchMinimalLength(DEFAULT_RC_MATCH_MINIMUM_LENGTH);
                rcMatchMinLengthFixed = false;
            }
        }
    }

    void setProteinsCompressionProfile() {
        setKmerLength(PROTEINS_PROFILE_KMER_LENGTH);
        disableMismatchesEncodingWithExclusion();
        if (!rcMatchMinLengthFixed) {
            setRCMatchMinimalLength(0);
            rcMatchMinLengthFixed = false;
        }
    }

    void setUltraStreamsCompression() {
        ultraStreamsCompression = true;
    }

    void setGzDecompressionMode(int gzCoderLevel) {
        if (gzCoderLevel > 12 || gzCoderLevel < 1) {
            fprintf(stderr, "Gz compression level should be between %d and %d.\n",
                    1, 12);
            exit(EXIT_FAILURE);
        }
        decompressionToGzCoderLevel = gzCoderLevel;
    }
};

#endif //MBGC_MBGCPARAMS_H

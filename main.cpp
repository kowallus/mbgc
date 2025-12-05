#include <cstdlib>
#include <unistd.h>

#include "mbgccoder/MBGC_Encoder.h"
#include "mbgccoder/MBGC_Decoder.h"
#include "mbgccoder/MBGC_Params.h"

#include <omp.h>

#ifdef __MINGW32__
#include <fcntl.h>
#endif

#define RELEASE_DATE "2025-12-02"

using namespace std;


void printVersion(bool details) {
    string date = RELEASE_DATE;
    if (!details)
        date.resize(4);
    fprintf(stderr, "Multiple Bacteria Genome Compressor (MBGC) v%d.%d.%d (c) Tomasz Kowalski, Szymon Grabowski, %s\n",
            (int) MBGC_Params::MBGC_VERSION_MAJOR, (int) MBGC_Params::MBGC_VERSION_MINOR,
            (int) MBGC_Params::MBGC_VERSION_REVISION, date.c_str());
}

void printCommonOpts() {
    fprintf(stderr,
            "\t[-%c <noOfThreads>] set limit of used threads\n"
            "\t[-%c] ignore FASTA file paths (use only filenames)\n",
            MBGC_Params::NO_OF_THREADS_OPT, MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT);
    fprintf(stderr,
            "\t[-%c] print full command help and exit\n"
            "\t[-%c] print version number and exit\n",
            MBGC_Params::HELP_OPT, MBGC_Params::VERSION_OPT);
}

#ifdef DEVELOPER_BUILD
void printCommonDevOpts() {
    fprintf(stderr, "[-%c] printout params\n", MBGC_Params::PRINTOUT_PARAMS_OPT);
    fprintf(stderr, "[-%c <logfilename>] printout params to tab-separated filelog\n", MBGC_Params::PARAMS_TO_LOGFILE_OPT);
}

void printCompressDeveloperOptions(bool repackCommand) {
    fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
    if (repackCommand)
        fprintf(stderr, "[-%c] disable lazy decompression\n", MBGC_Params::DISABLE_LAZY_DECODING_OPT);
    fprintf(stderr, "[-%c] brute parallel encoding mode (without producer-consumer approach)\n",
            MBGC_Params::BRUTE_PARALLEL_MODE_OPT);
    fprintf(stderr, "[-%c] stronger backend compression (slow)\n", MBGC_Params::ULTRA_STREAMS_COMPRESSION_OPT);
    fprintf(stderr, "[-%c] protein sequences compression profile\n", MBGC_Params::PROTEINS_PROFILE_OPT);
    fprintf(stderr, "[-%c <skipMargin>] (0 <= %c <= 255, default: %d)\n", MBGC_Params::SKIP_MARGIN_OPT, MBGC_Params::SKIP_MARGIN_OPT, MBGC_Params::DEFAULT_SKIP_MARGIN);
    fprintf(stderr, "[-%c] limit reference buffer size up to 2^32 bytes\n", MBGC_Params::LIMIT_32_BIT_REF_OPT);
    fprintf(stderr, "[-%c <referenceSlidingWindowFactor>] (2 <= %c <= 255, default: %d)\n",
            MBGC_Params::REF_SLIDING_WINDOW_FACTOR_OPT, MBGC_Params::REF_SLIDING_WINDOW_FACTOR_OPT, MBGC_Params::DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR);
    fprintf(stderr, "--- MISMATCHES ENCODING AND MATCHING SCORE OPTIONS ---\n");
    fprintf(stderr, "[-%c <gapDepthMismatchesEncoding>] (0 <= %c <= %d, default: %d, 0 - disable)\n"
                    "[-%c <mismatchesLimitInPercent>] (0 < %c < 100, default: %d)\n"
                    "[-%c <initialScoreInPercent>] (0 <= %c <= 100, default: %d)\n",
            MBGC_Params::GAP_DEPTH_MISMATCHES_ENCODING_OPT, MBGC_Params::GAP_DEPTH_MISMATCHES_ENCODING_OPT, MBGC_Params::MAX_GAP_DEPTH, MBGC_Params::DEFAULT_GAP_DEPTH_MISMATCHES_ENCODING,
            MBGC_Params::MISMATCHES_LIMIT_IN_PERCENT_OPT, MBGC_Params::MISMATCHES_LIMIT_IN_PERCENT_OPT, MBGC_Params::DEFAULT_MISMATCHES_LIMIT_IN_PERCENT,
            MBGC_Params::MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT, MBGC_Params::MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT, MBGC_Params::DEFAULT_MISMATCHES_INITIAL_SCORE_IN_PERCENT);
    fprintf(stderr, "--- EXPERIMENTAL REF EXTENSION STRATEGIES ---\n");
    fprintf(stderr, "[-%c] combined ref extension strategy\n", MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_OPT);
    fprintf(stderr, "[-%c] literals ref extension strategy (experimental - only for encoding)\n"
                    "[-%c <literalContextLength>] (0 <= %c, default: %ld)\n"
                    "[-%c <minimalLiteralLength>] (1 <= %c, default: DISABLED)\n",
            MBGC_Params::ONLY_LITERAL_REF_EXTENSION_STRATEGY_OPT,
            MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT_OPT, MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT_OPT, MBGC_Params::DEFAULT_REF_LITERAL_BEFORE_AFTER_EXT,
            MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_OPT, MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_OPT);
    fprintf(stderr, "[-%c] split contigs into blocks (experimental - only for encoding)\n", MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_OPT);
    fprintf(stderr, "[-%c] dynamic unmatched fraction ref extension strategy\n", MBGC_Params::DYNAMIC_REF_EXTENSION_FACTOR_OPT);
    fprintf(stderr, "--- DISSIMILAR CONTIGS TWEAK OPTIONS ---\n");
    fprintf(stderr, "[-%c <allowedTargetsOutrun>] (default: %d)\n",
            MBGC_Params::ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT, MBGC_Params::DEFAULT_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS);
    fprintf(stderr, "[-%c <minimalDissimilarContigLength>] (default: %ld)\n",
            MBGC_Params::MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT, MBGC_Params::DEFAULT_MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS);
    fprintf(stderr, "[-%c <dissimilarContigUnmatchedFractionTweakFactor>] (default: %d)\n",
            MBGC_Params::UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT, MBGC_Params::DEFAULT_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS);
    fprintf(stderr, "--- INPUT/OUTPUT OPTIONS ---\n");
    if (!repackCommand)
        fprintf(stderr, "[-%c] interleave input files order\n", MBGC_Params::INTERLEAVE_FILE_ORDER_OPT);
    printCommonDevOpts();
}
#endif

void printEncodeUsage(string base_cmd_name, bool details, bool repackMode) {
    printVersion(details);
    if (repackMode) {
        fprintf(stderr, "\nUsage for repacking part of archive:\n\t%s [-%c <compressionMode>]"
                        " [-%c <pattern>] [-%c <patternsFile>]\n\t\t<archiveFile> <outputArchiveFile>\n", base_cmd_name.c_str(),
                        MBGC_Params::COMPRESSION_MODE_OPT,
                        MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
        if (details) {
            fprintf(stderr,
                    "\n<archiveFile> mbgc archive filename for repacking\n"
                    "\tfor standard input, set <archiveFile> to %s\n"
                    "<patternsFile> name of text file with list of patterns (in separate lines)\n"
                    "\texcludes files not matching any pattern (does not invalidate -%c option)\n"
                    "<outputArchiveFile> name of repacked mbgc archive\n"
                    "\tfor standard output, set <outputArchiveFile> to %s\n",
                    MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::STANDARD_IO_POSIX_ALIAS);
        }
    } else {
        fprintf(stderr, "\nUsage for multiple file compression (list of files given as input):\n\t%s [-%c <compressionMode>]"
                        " [-%c] <fastaFileList> <archiveFile>\n", base_cmd_name.c_str(),
                MBGC_Params::COMPRESSION_MODE_OPT, MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT);
        fprintf(stderr, "Usage for single file compression:\n\t%s [-%c <compressionMode>]"
                        " [-%c] -%c <inputFastaFile> <archiveFile>\n", base_cmd_name.c_str(),
                MBGC_Params::COMPRESSION_MODE_OPT, MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT, MBGC_Params::INPUT_FILE_OPT);

        if (details) {
            fprintf(stderr,
                    "\n<fastaFileList> name of text file with a list of FASTA files (raw or gz)\n"
                    "\t(given in separate lines) for compression\n"
                    "<inputFastaFile> name of a FASTA file (raw or gz) for compression\n"
                    "\tfor standard input, set <inputFastaFile> to %s\n"
                    "<archiveFile> mbgc archive filename\n"
                    "\tfor standard output in compression, set <archiveFile> to %s\n",
                    MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS);
        }
    }
    fprintf(stderr, "\nBasic options:\n\t[-%c <compressionMode>] "
                    "(speed: 0; default: 1; repo: 2; max: 3)\n"
                    "\t[-%c] allow lossy compression\n"
                    "\t[-%c] overwrite an existing output file\n",
                    MBGC_Params::COMPRESSION_MODE_OPT, MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT,
                    MBGC_Params::FORCE_OVERWRITE_OPT);

    if (repackMode) {
        fprintf(stderr,
                "\t[-%c <basesPerRow>] custom format of repacked DNA (0 - unlimited)\n"
                "\t[-%c <pattern>] exclude files with names not containing pattern\n"
                "\t[-%c <patternsFile>] exclude files not matching any pattern\n",
                MBGC_Params::DNA_LINE_LENGTH_OPT, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
    }
    printCommonOpts();
    if (details) {
        fprintf(stderr, "\nCompression modes description:\n"
                        "\t(0) speed - for speed (fastest compression and decompression)\n"
                        "\t(1) default - regular mode (good ratio, fast)\n"
                        "\t(2) repo - for public repositories (better ratio, good speed)\n"
                        "\t(3) max - for long-term storage (best ratio, memory-frugal)\n");
    }
}

void printEncodeAdvancedDetails() {
    fprintf(stderr, "\n------------------ ADVANCED OPTIONS ----------------\n");
    fprintf(stderr, "[-%c <noOfWorkers>] worker threads used by matcher (default: %d)\n",
                    MBGC_Params::MATCHER_NO_OF_THREADS_OPT, MBGC_Params::DEFAULT_MATCHER_NO_OF_THREADS);
    fprintf(stderr, "[-%c] converts bases to uppercase\n",
            MBGC_Params::UPPERCASE_DNA_OPT);
    fprintf(stderr, "[-%c] disable parallel matching (does not apply to I/O and backend compression)\n"
                    "[-%c <matchingKmerLength>] (%d <= %c <= 40, default: %d)\n"
                    "[-%c <unmatchedFractionFactor>] (1 <= %c <= 255, default: %d)\n",
            MBGC_Params::SEQ_MATCHING_MODE_OPT, MBGC_Params::KMER_LENGTH_OPT,
            MBGC_Params::MIN_MATCH_LENGTH, MBGC_Params::KMER_LENGTH_OPT, MBGC_Params::DEFAULT_KMER_LENGTH,
            MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT, MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT, MBGC_Params::DEFAULT_UNMATCHED_LENGTH_FACTOR);
    fprintf(stderr, "[-%c <unmatchedFractionRCFactor>] (0 <= %c <= 255, default: %d,\n\t\t0 disables rc-matches in reference)\n",
            MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT, MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT, MBGC_Params::DEFAULT_UNMATCHED_LENGTH_RC_FACTOR);
    fprintf(stderr, "[-%c <referenceSamplingStep>] (%c > 0, default: %d)\n"
                    "[-%c <referenceFactorBinaryOrder>] (0 <= %c <= 12, auto-adjusted by default)\n"
                    "[-%c] disable circular reference buffer\n",
            MBGC_Params::REF_SAMPLING_STEP_OPT, MBGC_Params::REF_SAMPLING_STEP_OPT, MBGC_Params::DEFAULT_MODE_REF_SAMPLING,
            MBGC_Params::REF_FACTOR_BIN_ORDER_OPT, MBGC_Params::REF_FACTOR_BIN_ORDER_OPT, MBGC_Params::CIRCULAR_REFERENCE_BUFFER_OPT);
    fprintf(stderr,
            "[-%c <gapDepthOffsetEncoding>] (0 <= %c <= %d, default: %d, 0 - disable)\n"
            "[-%c <gapBreakingMatchMinimalLength>] (%c > 0, default: %d, 0 - disable)\n"
            "[-%c <maxConsecutiveMismatches>] (0 <= %c <= 255, default: %d; 0 - disable;\n\t\t1 - single mismatch mode)\n"
            "[-%c] disable encoding mismatches with exclusion\n"
            "[-%c <rcMatchMinimalLength>] (%c >= 24, default: %d, 0 - disable)\n",
            MBGC_Params::GAP_DEPTH_OFFSET_ENCODING_OPT, MBGC_Params::GAP_DEPTH_OFFSET_ENCODING_OPT, MBGC_Params::MAX_GAP_DEPTH, MBGC_Params::DEFAULT_GAP_DEPTH_OFFSET_ENCODING,
            MBGC_Params::GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT, MBGC_Params::GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT, MBGC_Params::DEFAULT_GAP_BREAKING_MATCH_MIN_LENGTH,
            MBGC_Params::MAX_CONSECUTIVE_MISMATCHES_OPT, MBGC_Params::MAX_CONSECUTIVE_MISMATCHES_OPT, MBGC_Params::DEFAULT_MAX_CONSECUTIVE_MISMATCHES,
            MBGC_Params::DISABLE_MISMATCHES_WITH_EXCLUSION_OPT,
            MBGC_Params::RC_MATCH_MINIMAL_LENGTH_OPT, MBGC_Params::RC_MATCH_MINIMAL_LENGTH_OPT, 0);
}

int encodeCommand(std::string base_cmd_name, int argc, char *argv[]) {
    base_cmd_name += string(" ") + string(argv[0]);
    char cmdId = argv[0][0];
    argv[0] = (char *) base_cmd_name.c_str();
    MBGC_Params params;
    params.repackCommand = cmdId == MBGC_Params::REPACK_CMD;
    const string SHORT_OPTS = params.repackCommand ? MBGC_Params::REPACK_SHORT_OPTS : MBGC_Params::COMPRESS_SHORT_OPTS;
    if (argc == 1) {
        printEncodeUsage(base_cmd_name, false, params.repackCommand);
        fprintf(stderr, "\n");
        return EXIT_SUCCESS;
    }
    int opt; // current option
    bool cliError = false;
    bool printParamsFlag = false;
    params.coderMode = MBGC_Params::DEFAULT_CODER_MODE;
    MBGC_Params repackInputParams;

    while ((opt = getopt(argc, argv, SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                params.forceOverwrite = true;
                break;
            case MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT:
                params.allowLossyCompression = true;
                break;
            case MBGC_Params::NO_OF_THREADS_OPT:
                params.setThreadsLimit(atoi(optarg));
                break;
            case MBGC_Params::MATCHER_NO_OF_THREADS_OPT:
                params.setMatcherThreads(atoi(optarg));
                break;
            case MBGC_Params::COMPRESSION_MODE_OPT:
                params.setCompressionMode(atoi(optarg));
                break;
            case MBGC_Params::INPUT_FILE_OPT:
                params.setInputFileName(string(optarg));
                break;
            case MBGC_Params::SEQ_MATCHING_MODE_OPT:
                params.setSequentialMatchingMode();
                break;
            case MBGC_Params::KMER_LENGTH_OPT:
                params.setKmerLength(atoi(optarg));
                break;
            case MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT:
                params.setUnmatchedFractionFactor(atoi(optarg));
                break;
            case MBGC_Params::REF_SAMPLING_STEP_OPT:
                params.setReferenceSamplingStep(atoi(optarg));
                break;
            case MBGC_Params::REF_FACTOR_BIN_ORDER_OPT:
                params.setReferenceFactorBinaryOrder(atoi(optarg));
                break;
            case MBGC_Params::CIRCULAR_REFERENCE_BUFFER_OPT:
                params.disableCircularReference();
                break;
            case MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT:
                params.setUnmatchedFractionRCFactor(atoi(optarg));
                break;
            case MBGC_Params::GAP_DEPTH_OFFSET_ENCODING_OPT:
                params.setGapDepthOffsetEncoding(atoi(optarg));
                break;
            case MBGC_Params::GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT:
                params.setGapBreakingMatchMinimalLength(atoi(optarg));
                break;
            case MBGC_Params::MAX_CONSECUTIVE_MISMATCHES_OPT:
                params.setMaxConsecutiveMismatches(atoi(optarg));
                break;
            case MBGC_Params::DISABLE_MISMATCHES_WITH_EXCLUSION_OPT:
                params.disableMismatchesEncodingWithExclusion();
                break;
            case MBGC_Params::RC_MATCH_MINIMAL_LENGTH_OPT:
                params.setRCMatchMinimalLength(atoi(optarg));
                break;
            case MBGC_Params::FILTER_PATTERN_OPT:
                repackInputParams.setFilterPattern(string(optarg));
                break;
            case MBGC_Params::PATTERNS_FILE_OPT:
                repackInputParams.setInputFileName(string(optarg));
                break;
            case MBGC_Params::DNA_LINE_LENGTH_OPT:
                repackInputParams.setDNALineLength(atoll(optarg));
                break;
            case MBGC_Params::UPPERCASE_DNA_OPT:
                params.uppercaseDNA = true;
                break;
#ifdef DEVELOPER_BUILD
            case MBGC_Params::DISABLE_LAZY_DECODING_OPT:
                repackInputParams.disableLazyDecompression = true;
                break;
            case MBGC_Params::BRUTE_PARALLEL_MODE_OPT:
                params.setBruteParallelMode();
                break;
            case MBGC_Params::ULTRA_STREAMS_COMPRESSION_OPT:
                params.setUltraStreamsCompression();
                break;
            case MBGC_Params::SKIP_MARGIN_OPT:
                params.setSkipMargin(atoi(optarg));
                break;
            case MBGC_Params::LIMIT_32_BIT_REF_OPT:
                params.limit32bitReference();
                break;
            case MBGC_Params::REF_SLIDING_WINDOW_FACTOR_OPT:
                params.setReferenceSlidingWindowFactor(atoi(optarg));
                break;
            case MBGC_Params::GAP_DEPTH_MISMATCHES_ENCODING_OPT:
                params.setGapDepthMismatchesEncoding(atoi(optarg));
                break;
            case MBGC_Params::MISMATCHES_LIMIT_IN_PERCENT_OPT:
                params.setMismatchesLimitInPercent(atoi(optarg));
                break;
            case MBGC_Params::MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT:
                params.setMismatchesInitialScoreInPercent(atoi(optarg));
                break;
            case MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_OPT:
                params.refExtensionStrategy |= MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_MASK;
                break;
            case MBGC_Params::ONLY_LITERAL_REF_EXTENSION_STRATEGY_OPT:
                params.refExtensionStrategy |= MBGC_Params::ONLY_LITERAL_REF_EXTENSION_STRATEGY_MASK;
                break;
            case MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT_OPT:
                params.refLiteralBeforeAfterExt = atoi(optarg);
                break;
            case MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_OPT:
                params.refLiteralMinimalLengthExt = atoi(optarg);
                break;
            case MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_OPT:
                params.refExtensionStrategy |= MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_MASK;
                break;
            case MBGC_Params::DYNAMIC_REF_EXTENSION_FACTOR_OPT:
                params.refExtensionStrategy |= MBGC_Params::DYNAMIC_REF_EXT_FACTOR_MASK;
                break;
            case MBGC_Params::INTERLEAVE_FILE_ORDER_OPT:
                params.setInterleaveFileOrder();
                break;
            case MBGC_Params::ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT:
                params.allowedTargetsOutrunForDissimilarContigs = atoi(optarg);
                break;
            case MBGC_Params::MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT:
                params.minimalLengthForDissimilarContigs = atoi(optarg);
                break;
            case MBGC_Params::UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT:
                params.unmatchedFractionFactorTweakForDissimilarContigs = atoi(optarg);
                break;
            case MBGC_Params::PROTEINS_PROFILE_OPT:
                params.setProteinsCompressionProfile();
                break;
            case MBGC_Params::CONCAT_HEADERS_AND_SEQS_OPT:
                params.concatHeadersAndSequencesMode = true;
                break;
            case MBGC_Params::PRINTOUT_PARAMS_OPT:
                printParamsFlag = true;
                break;
            case MBGC_Params::PARAMS_TO_LOGFILE_OPT:
                params.enableLogFile(optarg);
                break;
#endif
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printEncodeUsage(base_cmd_name, true, params.repackCommand);
                printEncodeAdvancedDetails();
#ifdef DEVELOPER_BUILD
                printCompressDeveloperOptions(params.repackCommand);
#endif
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n\n");
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
    if (params.repackCommand) {
        if (optind > (argc - 2) || optind < (argc - 2)) {
            fprintf(stderr, "%s: For 'repack' command expected 2 arguments after options "
                            "(found %d)\n", base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
    } else {
        if (params.inputFileName.empty() && (optind > (argc - 2) || optind < (argc - 2))) {
            fprintf(stderr, "%s: For 'compress' command using files list as input expected 2 arguments after options "
                            "(found %d)\n", base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
        if (!params.inputFileName.empty() && (optind > (argc - 1) || optind < (argc - 1))) {
            fprintf(stderr, "%s: For 'compress' command given single file as input expected 1 argument after options "
                            "(found %d)\n", base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

    if (!params.repackCommand && params.inputFileName.empty())
        params.setSeqListFileName(argv[optind++]);
    if (params.repackCommand)
        repackInputParams.setInArchiveFileName(argv[optind++]);
    params.setOutArchiveFileName(argv[optind++]);
    if (params.outArchiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
        PgHelpers::devout = &null_stream;
    }
    MBGC_Encoder encoder(&params);
    if (params.repackCommand)
        encoder.repack(repackInputParams);
    else
        encoder.encode();
    *PgHelpers::appout << "encoding time - " << PgHelpers::time_millis() << " [ms]" << endl;

    if (printParamsFlag)
        params.printout();
    return EXIT_SUCCESS;
}

#ifdef DEVELOPER_BUILD
void printValidateUsage(string base_cmd_name, bool details) {
    printVersion(details);
    fprintf(stderr, "\nUsage for validation:\n\t%s [-%c <noOfThreads>] [-%c <dnaLineLength>]"
                    " <archiveFile>\n", base_cmd_name.c_str(),
            MBGC_Params::NO_OF_THREADS_OPT, MBGC_Params::DNA_LINE_LENGTH_OPT);
    fprintf(stderr, "Usage for partial validation:\n\t%s [-%c <pattern>] [-%c <patternsFile>]"
                    " <archiveFile>\n", base_cmd_name.c_str(), MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::INPUT_FILE_OPT);

    if (details) {
        fprintf(stderr,
                "\n<archiveFile> mbgc archive filename\n"
                "\tfor standard input in validation, set <archiveFile> to %s\n"
                "<patternsFile> name of text file with a list of patterns (in separate lines)\n"
                "\tfiles matching one of the patterns are validated obligatory (does not invalidate -%c option)\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::FILTER_PATTERN_OPT);
    }
    fprintf(stderr, "\nBasic options:\n"
                    "\t[-%c <basesPerRow>] custom format of decoded sequences (0 - unlimited)\n"
                    "\t[-%c <pattern>] validate files with names containing pattern\n"
                    "\t[-%c <patternsFile>] validate files with names containing any pattern\n",
            MBGC_Params::DNA_LINE_LENGTH_OPT, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
    fprintf(stderr, "\t[-%c] converts bases to uppercase\n",
            MBGC_Params::UPPERCASE_DNA_OPT);
    printCommonOpts();
}

void printDecompressDeveloperOptions(bool isListCmd, bool isValidationCmd) {
    fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
    fprintf(stderr, "[-%c] disable lazy decompression\n", MBGC_Params::DISABLE_LAZY_DECODING_OPT);
    if (isValidationCmd) {
        fprintf(stderr, "[-%c] dump streams to files\n", MBGC_Params::DUMP_STREAMS_OPT);
        fprintf(stderr, "[-%c] skip actual validation (measure decoding only)\n", MBGC_Params::SKIP_ACTUAL_VALIDATION_OPT);
    }
    if (!isListCmd && !isValidationCmd)
        fprintf(stderr, "[-%c] concatenate sequences and headers in each file\n", MBGC_Params::CONCAT_HEADERS_AND_SEQS_OPT);
    printCommonDevOpts();
}
#endif

void printDecodingUsage(string base_cmd_name, bool details, bool isListCmd) {
#ifdef DEVELOPER_BUILD
    if (base_cmd_name.back() == MBGC_Params::VALIDATE_CMD) {
        printValidateUsage(base_cmd_name, false);
        return;
    }
#endif
    printVersion(details);
    if (isListCmd) {
        fprintf(stderr, "\nUsage for partial file listing:\n\t%s [-%c] [-%c <pattern>] [-%c <patternsFile>]"
                        " <archiveFile>\n", base_cmd_name.c_str(),
                        MBGC_Params::LIST_HEADERS_OPT, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
    } else {
        fprintf(stderr, "\nUsage for decompression:\n\t%s [-%c <gzLevel>]"
                        " <archiveFile> [<outputPath>]\n", base_cmd_name.c_str(), MBGC_Params::GZ_DECOMPRESSION_OPT);
        fprintf(stderr, "Usage for partial decompression (list of patterns given as input):\n\t%s [-%c <patternsFile>]"
                        " <archiveFile> [<outputPath>]\n", base_cmd_name.c_str(), MBGC_Params::PATTERNS_FILE_OPT);
    }

    if (details) {
        fprintf(stderr,
                "\n<archiveFile> mbgc archive filename\n"
                "\tfor standard input in decompression, set <archiveFile> to %s\n"
                "<patternsFile> name of text file with list of patterns (in separate lines)\n"
                "\texcludes files not matching any pattern (does not invalidate -%c option)\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::FILTER_PATTERN_OPT);
        if (!isListCmd)
            fprintf(stderr,
                    "<outputPath> extraction target path root (current directory by default)\n"
                    "\tfor standard output, set <outputPath> to %s (all files are concatenated)\n",
                    MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }
    fprintf(stderr, "\nBasic options:\n");
    if (isListCmd)
        fprintf(stderr,
                "\t[-%c <pattern>] exclude files with names not containing pattern\n"
                "\t[-%c <patternsFile>] exclude files not matching any pattern\n"
                "\t[-%c] list sequence headers (using convention: \">sequencename>filename\")\n",
                MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT, MBGC_Params::LIST_HEADERS_OPT);
    else
        fprintf(stderr,
                "\t[-%c] overwrite an existing output files\n"
                "\t[-%c <gzLevel>] extract FASTA files to gz archives\n\t\t(compression level: 1 <= %c <= 12, recommended: %d)\n"
                "\t[-%c <basesPerRow>] custom format of decoded sequences (0 - unlimited)\n"
                "\t[-%c <pattern>] exclude files with names not containing pattern\n"
                "\t[-%c <patternsFile>] exclude files not matching any pattern\n",
                MBGC_Params::FORCE_OVERWRITE_OPT, MBGC_Params::GZ_DECOMPRESSION_OPT, MBGC_Params::GZ_DECOMPRESSION_OPT,
                MBGC_Params::RECOMMENDED_GZ_COMPRESSION_LEVEL, MBGC_Params::DNA_LINE_LENGTH_OPT,
                MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
    printCommonOpts();
}

int decodeCommands(std::string base_cmd_name, int argc, char *argv[]) {
    string defaultCommandWarning;
    if (strlen(argv[0]) > 1)
        defaultCommandWarning = string("\nNote: selected default command '") + MBGC_Params::INFO_CMD + "'. Please use '" + base_cmd_name + "' to list all commands.\n";
    char cmdId = strlen(argv[0]) > 1 ? MBGC_Params::INFO_CMD : argv[0][0];
    base_cmd_name += string(" ") + cmdId;
    argv[0] = (char *) base_cmd_name.c_str();
    MBGC_Params params;
    params.infoCommand = cmdId == MBGC_Params::INFO_CMD;
#ifndef DEVELOPER_BUILD
    const string SHORT_OPTS = params.infoCommand ? MBGC_Params::INFO_SHORT_OPTS : MBGC_Params::DECOMPRESS_SHORT_OPTS;
#else
    if (params.infoCommand)
        PgHelpers::devout = &std::cerr;
    bool isValidateCmd = cmdId == MBGC_Params::VALIDATE_CMD;
    params.validationMode = isValidateCmd;
    const string SHORT_OPTS = isValidateCmd ? MBGC_Params::VALIDATION_SHORT_OPTS :
                              params.infoCommand ? MBGC_Params::INFO_SHORT_OPTS : MBGC_Params::DECOMPRESS_SHORT_OPTS;
#endif
    if (argc == 1) {
        printDecodingUsage(base_cmd_name, false, params.infoCommand);
        fprintf(stderr, "%s\n", defaultCommandWarning.c_str());
        return EXIT_SUCCESS;
    }
    int opt; // current option
    bool cliError = false;
    bool printParamsFlag = false;

    while ((opt = getopt(argc, argv, SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                params.forceOverwrite = true;
                break;
            case MBGC_Params::NO_OF_THREADS_OPT:
                params.setThreadsLimit(atoi(optarg));
                break;
            case MBGC_Params::FILTER_PATTERN_OPT:
                params.setFilterPattern(string(optarg));
                break;
            case MBGC_Params::PATTERNS_FILE_OPT:
                params.setInputFileName(string(optarg));
                break;
            case MBGC_Params::GZ_DECOMPRESSION_OPT:
                params.setGzDecompressionMode(atoi(optarg));
                break;
            case MBGC_Params::DNA_LINE_LENGTH_OPT:
                params.setDNALineLength(atoll(optarg));
                break;
            case MBGC_Params::LIST_HEADERS_OPT:
                params.listHeadersMode = true;
                break;
#ifdef DEVELOPER_BUILD
            case MBGC_Params::DISABLE_LAZY_DECODING_OPT:
                params.disableLazyDecompression = true;
                break;
            case MBGC_Params::UPPERCASE_DNA_OPT:
                params.uppercaseDNA = true;
                break;
            case MBGC_Params::CONCAT_HEADERS_AND_SEQS_OPT:
                params.concatHeadersAndSequencesMode = true;
                break;
            case MBGC_Params::DUMP_STREAMS_OPT:
                dump_after_decompression = true;
                break;
            case MBGC_Params::SKIP_ACTUAL_VALIDATION_OPT:
                params.skipActualValidation = true;
                break;
            case MBGC_Params::PRINTOUT_PARAMS_OPT:
                printParamsFlag = true;
                break;
            case MBGC_Params::PARAMS_TO_LOGFILE_OPT:
                params.enableLogFile(optarg);
                break;
#endif
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printDecodingUsage(base_cmd_name, true, params.infoCommand);
#ifdef DEVELOPER_BUILD
                printDecompressDeveloperOptions(params.infoCommand, isValidateCmd);
#endif
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n");
                fprintf(stderr, "%s\n", defaultCommandWarning.c_str());
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
#ifdef DEVELOPER_BUILD
    if (isValidateCmd && (optind < (argc - 1) || optind > (argc - 1))) {
        fprintf(stderr, "%s: For 'validation' command expected 1 argument after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    } else
#endif
    if (params.infoCommand && (optind < (argc - 1) || optind > (argc - 1))) {
        fprintf(stderr, "%s: For 'info' command expected 1 argument after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    } else if (optind < (argc - 2) || optind > (argc - 1)) {
        fprintf(stderr, "%s: For 'decompress' command expected 1 or 2 arguments after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information%s\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT,
                defaultCommandWarning.c_str());
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

#ifdef DEVELOPER_BUILD
    dump_after_decompression_prefix = string(argv[optind]) + "_dump_";
#endif
    params.setInArchiveFileName(argv[optind++]);
    params.setOutputPath(argc > optind ? argv[optind++] : "");
    if (params.outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
        PgHelpers::devout = &null_stream;
    }
    *PgHelpers::logout << "            \t";
    MBGC_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&params, defaultCommandWarning);
    decoder->decode();
    delete(decoder);
    if (params.infoCommand)
        *PgHelpers::devout << "listing" << " time - " << PgHelpers::time_millis() << " [ms]" << endl;
    else
        *PgHelpers::appout << "decoding" << " time - " << PgHelpers::time_millis() << " [ms]" << endl;

    if (printParamsFlag) {
        params.printout();
    }
    return EXIT_SUCCESS;
}

#ifdef DEVELOPER_BUILD
void printAppendDeveloperOptions() {
    fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
    fprintf(stderr, "[-%c <matcherNumberOfThreads>] (default: %d)\n",
            MBGC_Params::MATCHER_NO_OF_THREADS_OPT, MBGC_Params::DEFAULT_MATCHER_NO_OF_THREADS);
    printCommonDevOpts();
}
#endif

void printAppendUsage(string base_cmd_name, bool details) {
    printVersion(details);
    fprintf(stderr, "\nUsage for appending multiple FASTA archive (list of files given as input):\n\t%s"
                    " <fastaFileList> <archiveFile> [<outputArchiveFile>]\n", base_cmd_name.c_str());
    fprintf(stderr, "Usage for appending single FASTA file archive:\n\t%s"
                    " -%c <inputFastaFile> <archiveFile> [<outputArchiveFile>]\n", base_cmd_name.c_str(), MBGC_Params::INPUT_FILE_OPT);

    if (details) {
        fprintf(stderr,
                "\n<fastaFileList> name of text file with a list of FASTA files (raw or gz)\n"
                "\t(given in separate lines) for appending to the <archiveFile>\n"
                "<inputFastaFile> name of a FASTA file (raw or gz) appended to the <archiveFile>\n"
                "\tfor standard input, set <inputFastaFile> to %s\n"
                "<archiveFile> mbgc archive filename to be appended\n"
                "\tfor standard input (and output if <outputArchiveFile> is not defined)\n"
                "\tset <archiveFile> to %s\n"
                "<outputArchiveFile> if defined, <archiveFile> archive remains unchanged\n"
                "\tand new archive is created\n"
                "\tfor standard output, set <outputArchiveFile> to %s\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }
    fprintf(stderr, "\nBasic options:\n");
    printCommonOpts();
}

int appendCommand(std::string base_cmd_name, int argc, char *argv[]) {
    base_cmd_name += string(" ") + string(argv[0]);
    argv[0] = (char *) base_cmd_name.c_str();
    if (argc == 1) {
        printAppendUsage(base_cmd_name, false);
        fprintf(stderr, "\n");
        return EXIT_SUCCESS;
    }
    int opt; // current option
    bool cliError = false;
    bool printParamsFlag = false;
    MBGC_Params params, newParams;
    newParams.coderThreads = -1;
    newParams.matcherWorkingThreads = -1;
    params.coderMode = MBGC_Params::DEFAULT_CODER_MODE;

    while ((opt = getopt(argc, argv, MBGC_Params::APPEND_SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case MBGC_Params::NO_OF_THREADS_OPT:
                newParams.coderThreads = atoi(optarg);
                params.setThreadsLimit(newParams.coderThreads);
                break;
            case MBGC_Params::MATCHER_NO_OF_THREADS_OPT:
                newParams.matcherWorkingThreads = atoi(optarg);
                params.setMatcherThreads(newParams.matcherWorkingThreads);
                break;
            case MBGC_Params::INPUT_FILE_OPT:
                params.setInputFileName(string(optarg));
                break;
#ifdef DEVELOPER_BUILD
            case MBGC_Params::PRINTOUT_PARAMS_OPT:
                printParamsFlag = true;
                break;
            case MBGC_Params::PARAMS_TO_LOGFILE_OPT:
                params.enableLogFile(optarg);
                break;
#endif
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printAppendUsage(base_cmd_name, true);
#ifdef DEVELOPER_BUILD
                printAppendDeveloperOptions();
#endif
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n\n");
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
    if (params.inputFileName.empty() && (optind > (argc - 2) || optind < (argc - 3))) {
        fprintf(stderr, "%s: For 'append' command using files list as input expected 2 or 3 arguments after options "
                        "(found %d)\n", base_cmd_name.c_str(), argc - optind);
        cliError = true;
    }
    if (!params.inputFileName.empty() && (optind > (argc - 1) || optind < (argc - 2))) {
        fprintf(stderr, "%s: For 'append' command given single file as input expected 1 or 2 arguments after options "
                        "(found %d)\n", base_cmd_name.c_str(), argc - optind);
        cliError = true;
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

    if (params.inputFileName.empty())
        params.setSeqListFileName(argv[optind++]);
    params.setInArchiveFileName(argv[optind++]);
    params.setOutArchiveFileName(optind < argc ? argv[optind++] : params.inArchiveFileName);
    if (params.outArchiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
        PgHelpers::devout = &null_stream;
    }

    MBGC_Encoder encoder(&params);
    encoder.append(newParams);

    *PgHelpers::appout << "appending time - " << PgHelpers::time_millis() << " [ms]" << endl;

    if (printParamsFlag)
        params.printout();
    return EXIT_SUCCESS;
}

void printCommands(string base_toolname) {
    fprintf(stderr, "\nUsage:\t%s <command> [options]\n", base_toolname.c_str());
    fprintf(stderr, "\nAvailable commands (%c - default):\n", MBGC_Params::INFO_CMD);
    for(int i = 0; i < MBGC_Params::cliCommandsCount; i++) {
        fprintf(stderr, "\t%c\t%s\n",
                MBGC_Params::CLI_COMMANDS[i][0], MBGC_Params::CLI_COMMANDS[i]);
    }
}

int printCommandsUsage(string base_toolname, string unknown_command_arg = "") {
    bool isCommandWrong = !unknown_command_arg.empty();
    if (!isCommandWrong)
        printVersion(false);
    else
        fprintf(stderr, "%s: invalid command -- '%s'\n", base_toolname.c_str(), unknown_command_arg.c_str());
    printCommands(base_toolname);
    if (!isCommandWrong) {
        fprintf(stderr, "\nCommon options:\n");
        printCommonOpts();
#ifdef DEVELOPER_BUILD
        fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
        printCommonDevOpts();
#endif
    }
    fprintf(stderr, "\n");
    return isCommandWrong ? EXIT_FAILURE : EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    std::string path = argv[0];
    std::string base_cmd_name = path.substr(path.find_last_of("/\\") + 1);
    if (argc == 1)
        return printCommandsUsage(base_cmd_name);
#ifndef DEVELOPER_BUILD
    PgHelpers::devout = &null_stream;
#endif
    string command_arg = string(argv[1]);
    if (command_arg.size() > 1)
        return decodeCommands(base_cmd_name, argc, argv);

    switch (command_arg[0]) {
        case MBGC_Params::COMPRESS_CMD:
        case MBGC_Params::REPACK_CMD:
            return encodeCommand(base_cmd_name, argc - 1, argv + 1);
        case MBGC_Params::DECOMPRESS_CMD:
        case MBGC_Params::INFO_CMD:
            return decodeCommands(base_cmd_name, argc - 1, argv + 1);
        case MBGC_Params::APPEND_CMD:
            return appendCommand(base_cmd_name, argc - 1, argv + 1);
#ifdef DEVELOPER_BUILD
        case MBGC_Params::VALIDATE_CMD:
            return decodeCommands(base_cmd_name, argc - 1, argv + 1);
#endif
        default: /* '?' */
            printCommandsUsage(base_cmd_name, command_arg);
            return EXIT_FAILURE;
    }
}


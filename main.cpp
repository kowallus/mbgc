#include <cstdlib>
#include <unistd.h>

#include "mbgccoder/MBGC_Encoder.h"
#include "mbgccoder/MBGC_Decoder.h"

#include <omp.h>

#define RELEASE_DATE "22-02-2022"

using namespace std;


void printVersion() {
    fprintf(stderr, "Multiple Bacteria Genome Compressor (MBGC) v%d.%d.%dPE (c) Szymon Grabowski, Tomasz Kowalski, 20%c%c\n",
            (int) MBGC_Params::MBGC_VERSION_MAJOR, (int) MBGC_Params::MBGC_VERSION_MINOR,
            (int) MBGC_Params::MBGC_VERSION_REVISION, RELEASE_DATE[8], RELEASE_DATE[9]);
#ifdef NO_GZ_SUPPORT
    fprintf(stderr, "*the build doesn't support input from gz archives\n\n");
#endif
    fprintf(stderr, "\n");
}

void printUsage(string base_toolname, bool details) {
    printVersion();

    fprintf(stderr, "Usage for multiple file compression (list of files given as input):\n\t%s [-c compressionMode]"
                    " [-t noOfThreads] <sequencesListFile> <archiveFile>\n", base_toolname.c_str());
    fprintf(stderr, "Usage for single file compression:\n\t%s [-c compressionMode] [-t noOfThreads]"
                    " -i <inputFastaFile> <archiveFile>\n", base_toolname.c_str());
    fprintf(stderr, "Usage for decompression:\n\t%s -d [-t noOfThreads] [-f pattern] [-l dnaLineLength]"
                    " <archiveFile> [<outputPath>]\n\n", base_toolname.c_str());
    if (details) {
        fprintf(stderr,
                "<sequencesListFile> name of text file containing a list of FASTA files (raw or in gz archives)\n"
                "\t(given in separate lines) for compression\n"
#ifndef NO_GZ_SUPPORT
                "<inputFastaFile> name of a FASTA file (raw or in gz archive) for compression\n"
#else
                "<inputFastaFile> name of a FASTA file (only raw) for compression\n"
#endif
                "\tfor standard input set <inputFastaFile> to %s\n"
                "<archiveFile> mbgc archive filename\n"
                "\tfor standard input (resp. output) in compression (resp. decompression) set <archiveFile> to %s\n"
                "<outputPath> extraction target path root (if skipped the root path is the current directory)\n"
                "\tfor standard output set <outputPath> to %s (all files are concatenated)\n\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS,
                MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }
    fprintf(stderr, "Basic options:\n\t-c select compression mode (speed: 0; default: 1; repo: 2; max: 3)\n"
                    "\t-d decompression mode\n"
                    "\t-l format decompressed DNA (i.e., sets the number of bases per row)\n"
                    "\t-f decompress files with names containing the given pattern\n"
                    "\t-t number of threads used (default: %d)\n"
                    "\t-h print full help and exit\n"
                    "\t-v print version number and exit\n\n", MBGC_Params::DEFAULT_NO_OF_THREADS);
    if (details) {
        fprintf(stderr, "Compression modes description:\n"
                        "\t(0) speed  - for speed (fast compression and decompression)\n"
                        "\t(1) default - regular mode (fast compression and good ratio)\n"
                        "\t(2) repo - for public repositories (better ratio and fast decompression)\n"
                        "\t(3) max - for long-term storage (the best ratio)\n\n");
    }
}

void printAdvancedDetails() {
    fprintf(stderr, "------------------ ADVANCED OPTIONS ----------------\n");
    fprintf(stderr, "[-k matchingKmerLength] (24 <= k <= 40, default: 32)\n"
            "[-s referenceSamplingStep] (s > 0, default: 16)\n"
            "[-m skipMargin] (0 <= m <= 255, default: 16)\n"
            "[-o referenceFactorBinaryOrder] (0 <= o <= 12, auto-adjusted by default)\n"
            "[-u unmatchedFractionFactor] (1 <= u <= 255, default: 192)\n"
            "[-w referenceSlidingWindowFactor] (2 <= w <= 255, default: 16)\n"
            "[-r] limit reference buffer size up to 2^32 bytes\n"
            "[-p] disable matching parallelization (except for I/O and backend compression)\n\n");
}

#ifdef DEVELOPER_BUILD
void printDeveloperOptions() {
    fprintf(stderr, "------------------ DEVELOPER OPTIONS ----------------\n");
    fprintf(stderr, "[-b] brute parallel encoding mode (without producer-consumer approach)\n");
    fprintf(stderr, "[-R] no reverse-compliments in reference\n");
    fprintf(stderr, "[-S] stronger backend compression (slow)\n");
    fprintf(stderr, "[-L] literals ref extension strategy (experimental - only for encoding)\n");
    fprintf(stderr, "[-C] combined ref extension strategy\n");
    fprintf(stderr, "[-B] split contigs into blocks (experimental - only for encoding)\n");
    fprintf(stderr, "[-U] dynamic unmatched fraction ref extension strategy\n");
    fprintf(stderr, "[-I] interleave order of files \n\n");
    fprintf(stderr, "[-V] validation mode (disables file extraction)\n");
    fprintf(stderr, "[-O] concatenate sequences and headers in each file\n");
    fprintf(stderr, "[-D] dump streams to files during decompression\n\n");
}
#endif

int main(int argc, char *argv[]) {
    std::string path = argv[0];
    std::string base_toolname = path.substr(path.find_last_of("/\\") + 1);
    if (argc == 1) {
        printUsage(base_toolname, false);
        exit(EXIT_SUCCESS);
    }

    int opt; // current option
    bool decompressMode = false;
    bool decoderParamPresent = false;
    bool encoderParamPresent = false;
    bool cliError = false;
    bool versionFlag = false;
    MBGC_Params params;
    PgHelpers::numberOfThreads = MBGC_Params::DEFAULT_NO_OF_THREADS;
    params.coderLevel = CODER_LEVEL_NORMAL;

#ifndef DEVELOPER_BUILD
    PgHelpers::devout = &null_stream;
#endif

#ifdef DEVELOPER_BUILD
    while ((opt = getopt(argc, argv, "c:i:t:f:k:s:m:o:u:w:l:pdrvhbRLCBUSIDOV")) != -1) {
#else
    while ((opt = getopt(argc, argv, "c:i:t:f:k:s:m:o:u:w:l:pdrvh")) != -1) {
#endif
        switch (opt) {
            case 'c':
                encoderParamPresent = true;
                params.setCompressionMode(atoi(optarg));
                break;
            case 'd':
                decoderParamPresent = true;
                decompressMode = true;
                break;
            case 'i':
                encoderParamPresent = true;
                params.setInputFastaFileName(string(optarg));
                break;
            case 'f':
                decoderParamPresent = true;
                params.setFilterPattern(string(optarg));
                break;
            case 't':
                PgHelpers::numberOfThreads = atoi(optarg);
                params.noOfThreadsFixed = true;
                break;
            case 'k':
                encoderParamPresent = true;
                params.setKmerLength(atoi(optarg));
                break;
            case 's':
                encoderParamPresent = true;
                params.setReferenceSamplingStep(atoi(optarg));
                break;
            case 'm':
                encoderParamPresent = true;
                params.setSkipMargin(atoi(optarg));
                break;
            case 'o':
                encoderParamPresent = true;
                params.setReferenceFactorBinaryOrder(atoi(optarg));
                break;
            case 'p':
                encoderParamPresent = true;
                params.setSequentialMatchingMode();
                break;
            case 'r':
                encoderParamPresent = true;
                params.limit32bitReference();
                break;
            case 'u':
                encoderParamPresent = true;
                params.setUnmatchedFractionFactor(atoi(optarg));
                break;
            case 'w':
                encoderParamPresent = true;
                params.setReferenceSlidingWindowFactor(atoi(optarg));
                break;
            case 'l':
                decoderParamPresent = true;
                params.setDNALineLength(atoi(optarg));
                break;
#ifdef DEVELOPER_BUILD
            case 'b':
                encoderParamPresent = true;
                params.setBruteParallelMode();
                break;
            case 'R':
                encoderParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::NO_RC_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'L':
                encoderParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::LITERAL_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'C':
                encoderParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'B':
                encoderParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'U':
                encoderParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::DYNAMIC_REF_EXT_FACTOR_MASK;
                break;
            case 'S':
                encoderParamPresent = true;
                params.setUltraStreamsCompression();
                break;
            case 'I':
                encoderParamPresent = true;
                params.setInterleaveFileOrder();
                break;
            case 'V':
                params.validationMode = true;
                break;
            case 'D':
                dump_after_decompression = true;
                break;
            case 'O':
                params.concatHeadersAndSequencesMode = true;
                break;
#endif
            case 'v':
                versionFlag = true;
                break;
            case 'h':
                printUsage(base_toolname, true);
                printAdvancedDetails();
#ifdef DEVELOPER_BUILD
                printDeveloperOptions();
#endif
                fprintf(stderr, "The order of all selected options is arbitrary.\n\n");
                exit(EXIT_SUCCESS);
            default: /* '?' */
                cliError = true;
        }
    }
    if (versionFlag) {
        printVersion();
        exit(EXIT_SUCCESS);
    }
    if (decompressMode && (optind < (argc - 2) || optind > (argc - 1))) {
        fprintf(stderr, "%s: For decompression expected 1 or 2 arguments after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    }
    if (!decompressMode && !params.singleFastaFileMode && (optind > (argc - 2) || optind < (argc - 2))) {
        fprintf(stderr, "%s: For compression using files list as input expected 2 arguments after options "
                        "(found %d)\n", argv[0], argc - optind);
        cliError = true;
    }
    if (!decompressMode && params.singleFastaFileMode && (optind > (argc - 1) || optind < (argc - 1))) {
        fprintf(stderr, "%s: For compression using files list as input expected 1 argument after options "
                        "(found %d)\n", argv[0], argc - optind);
        cliError = true;
    }
    if (!params.filterPattern.empty() && encoderParamPresent) {
        fprintf(stderr, "Cannot use file filter pattern in compression mode.\n");
        cliError = true;
    }
    if (decoderParamPresent && encoderParamPresent) {
        fprintf(stderr, "Cannot use compression and decompression options together.\n");
        cliError = true;
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -h' for more information\n", base_toolname.c_str());
        exit(EXIT_FAILURE);
    }
    omp_set_num_threads(PgHelpers::numberOfThreads);
    PgHelpers::time_checkpoint();

    if (!decompressMode) {
        if (!params.singleFastaFileMode)
            params.setSeqListFileName(argv[optind++]);
        params.setArchiveFileName(argv[optind++]);
        if (params.archiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
            PgHelpers::appout = &null_stream;
            PgHelpers::devout = &null_stream;
        }
#pragma omp parallel
#pragma omp single
        *PgHelpers::devout << "number of threads: " << omp_get_num_threads() << endl;
        MBGC_Encoder encoder(&params);
        encoder.encode();
        *PgHelpers::appout << "encoding time - " << PgHelpers::time_millis() << " [ms]" << endl;
    }
#ifdef DEVELOPER_BUILD
    if (params.validationMode || decompressMode){
        dump_after_decompression_prefix = string(argv[2]) + "_dump_";
#else
    if (decompressMode){
#endif
        if (decompressMode)
            params.setArchiveFileName(argv[optind++]);
        params.setOutputPath(argc > optind ? argv[optind++] : "");
        if (params.outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
            PgHelpers::appout = &null_stream;
            PgHelpers::devout = &null_stream;
        }
#pragma omp parallel
#pragma omp single
        *PgHelpers::devout << "number of threads: " << omp_get_num_threads() << endl;
        MBGC_Decoder decoder(&params);
        decoder.decode();
        *PgHelpers::appout << (decompressMode?"decoding":"total") << " time - " << PgHelpers::time_millis() << " [ms]" << endl;
    }
    return EXIT_SUCCESS;
}
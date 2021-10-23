#include <cstdlib>
#include <unistd.h>

#include "mbgccoder/MBGC_Encoder.h"
#include "mbgccoder/MBGC_Decoder.h"

#include <omp.h>

#define RELEASE_DATE "2021-10-22"

using namespace std;

void printBasicUsage(string path) {
    std::string base_filename = path.substr(path.find_last_of("/\\") + 1);
    fprintf(stderr, "Usage for compression using files list as input: %s [-t noOfThreads] "
                    "<sequencesListFile> <archiveFile>\n", base_filename.c_str());
    fprintf(stderr, "Usage for single file compression: %s [-t noOfThreads] "
                    "-i <inputFastaFile> <archiveFile>\n\n", base_filename.c_str());
    fprintf(stderr, "Usage for decompression: %s -d [-t noOfThreads] [-f pattern] [-l dnaLineLength] "
                    "<archiveFile> [<outputPath>]\n\n", base_filename.c_str());

    fprintf(stderr, "-t number of threads used (default: %d)\n", MBGC_Params::DEFAULT_NO_OF_THREADS);
    fprintf(stderr, "-d decompression mode\n"
                    "-l format decompressed DNA (i.e., sets the number of bases per row)\n"
                    "-f decompress files with names containing the given pattern\n\n");
}

void printVersion() {
    fprintf(stderr, "MBGC %d.%d: Copyright (c) 20%c%c Szymon Grabowski, Tomasz Kowalski : %s\n\n",
            (int) MBGC_Params::MBGC_VERSION_MAJOR, (int) MBGC_Params::MBGC_VERSION_MINOR,
            RELEASE_DATE[2], RELEASE_DATE[3], RELEASE_DATE);
}

void printOptionDetails() {
    fprintf(stderr, "<sequencesListFile> name of text file containing a list of FASTA files (raw or in gz archives)\n"
                    "\t(given in separate lines) for compression\n"
                    "<inputFastaFile> name of a FASTA file (raw or in gz archive) for compression\n"
                    "\tfor standard input set <inputFastaFile> to %s\n"
                    "<archiveFile> mbgc archive filename\n"
                    "\tfor standard input (resp. output) in compression (resp. decompression) set <archiveFile> to %s\n"
                    "<outputPath> extraction target path root (if skipped the root path is the current directory)\n"
                    "\tfor standard output set <outputPath> to %s (all files are concatenated)\n\n",
            MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS,
            MBGC_Params::STANDARD_IO_POSIX_ALIAS);

    fprintf(stderr, "------------------ ADVANCED OPTIONS ----------------\n");
    fprintf(stderr, "[-c compression level] (fast: 1; default: 2; max: 3)\n"
            "[-v] print version number and exit\n"
            "[-k matchingKmerLength] (24 <= k <= 40, default: 32)\n"
            "[-s referenceSamplingStep] (s > 0, default: 16)\n"
            "[-m skipMargin] (0 <= m <= 255, default: 16)\n"
            "[-o referenceFactorBinaryOrder] (0 <= o <= 12, auto-adjusted by default)\n"
            "[-u unmatchedFractionFactor] (1 <= u <= 255, default: 192)\n"
            "[-w referenceSlidingWindowFactor] (2 <= w <= 255, default: 16)\n"
            "[-r] enable reference buffer size up to 2^40 bytes\n"
            "[-p] disable matching parallelization (except for I/O and backend compression)\n\n");
}

#ifdef DEVELOPER_BUILD
void printDeveloperOptions() {
    fprintf(stderr, "------------------ DEVELOPER OPTIONS ----------------\n");
    fprintf(stderr, "[-b] brute parallel encoding mode (without producer-consumer approach)\n");
    fprintf(stderr, "[-R] no reverse-compliments in reference\n");
    fprintf(stderr, "[-H] stronger backend compression of headers (slow)\n");
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
    if (argc == 1) {
        printBasicUsage(argv[0]);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
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
    while ((opt = getopt(argc, argv, "c:i:t:f:k:s:m:o:u:w:l:pdrvbRLCBUHIDOV?")) != -1) {
#else
    while ((opt = getopt(argc, argv, "c:i:t:f:k:s:m:o:u:w:l:pdrv?")) != -1) {
#endif
        switch (opt) {
            case 'c':
                encoderParamPresent = true;
                params.setCompressionLevel(atoi(optarg));
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
                params.enable40bitReference();
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
            case 'H':
                encoderParamPresent = true;
                params.setHeaderMaxCompression();
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
            case '?':
                printVersion();
                printBasicUsage(argv[0]);
                printOptionDetails();
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
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
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
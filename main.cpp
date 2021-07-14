#include <cstdlib>
#include <unistd.h>

#include "mbgccoder/MBGC_Encoder.h"
#include "mbgccoder/MBGC_Decoder.h"

#ifdef DEVELOPER_BUILD
#include "coders/CodersLib.h"
#endif

#include <omp.h>

#define RELEASE_DATE "2021-07-15"

using namespace std;

int main(int argc, char *argv[]) {

    int opt; // current option
    bool decompressMode = false;
    bool compressionParamPresent = false;
    MBGC_Params params;
    PgHelpers::numberOfThreads = MBGC_Params::DEFAULT_NO_OF_THREADS;
    params.coderLevel = CODER_LEVEL_NORMAL;

#ifndef DEVELOPER_BUILD
    NullBuffer null_buffer;
    std::ostream null_stream(&null_buffer);
    PgHelpers::logout = &null_stream;
#endif

#ifdef DEVELOPER_BUILD
    while ((opt = getopt(argc, argv, "c:t:f:k:s:m:o:u:w:pdbvRLCBUIFD?")) != -1) {
#else
    while ((opt = getopt(argc, argv, "c:t:f:k:s:m:o:u:w:pdb?")) != -1) {
#endif
        switch (opt) {
            case 'c':
                compressionParamPresent = true;
                params.setCompressionLevel(atoi(optarg));
                break;
            case 'd':
                decompressMode = true;
                break;
            case 'f':
                params.setFilterPattern(string(optarg));
                break;
            case 't':
                PgHelpers::numberOfThreads = atoi(optarg);
                params.noOfThreadsFixed = true;
                break;
            case 'k':
                compressionParamPresent = true;
                params.setKmerLength(atoi(optarg));
                break;
            case 's':
                compressionParamPresent = true;
                params.setReferenceSamplingStep(atoi(optarg));
                break;
            case 'm':
                compressionParamPresent = true;
                params.setSkipMargin(atoi(optarg));
                break;
            case 'o':
                compressionParamPresent = true;
                params.setReferenceFactorBinaryOrder(atoi(optarg));
                break;
            case 'b':
                compressionParamPresent = true;
                params.setBruteParallelMode();
                break;
            case 'p':
                compressionParamPresent = true;
                params.setSequentialMatchingMode();
                break;
            case 'u':
                compressionParamPresent = true;
                params.setUnmatchedFractionFactor(atoi(optarg));
                break;
            case 'w':
                compressionParamPresent = true;
                params.setReferenceSlidingWindowFactor(atoi(optarg));
                break;
#ifdef DEVELOPER_BUILD
            case 'R':
                compressionParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::NO_RC_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'L':
                compressionParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::LITERAL_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'C':
                compressionParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'B':
                compressionParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_MASK;
                break;
            case 'U':
                compressionParamPresent = true;
                params.refExtensionStrategy |= MBGC_Params::DYNAMIC_REF_EXT_FACTOR_MASK;
                break;
            case 'I':
                compressionParamPresent = true;
                params.setInterleaveFileOrder();
                break;
            case 'F':
                params.disableDNAformatting = true;
                break;
            case 'v':
                params.validationMode = true;
                break;
            case 'D':
                dump_after_decompression = true;
                break;
#endif
            case '?':
            default: /* '?' */
                fprintf(stderr, "MBGC %d.%d: Copyright (c) 2021 Szymon Grabowski, Tomasz Kowalski : %s\n\n",
                        (int) MBGC_Params::MBGC_VERSION_MAJOR, (int) MBGC_Params::MBGC_VERSION_MINOR, RELEASE_DATE);
                fprintf(stderr, "Usage for compression: %s [-t noOfThreads] "
                                "<sequencesListFile> <archiveFile>\n\n", argv[0]);
                fprintf(stderr, "Usage for decompression: %s -d [-t noOfThreads] [-f pattern] "
                                "<archiveFile> [<outputPath>]\n\n",
                        argv[0]);

                fprintf(stderr, "-t number of threads used (%d - default)\n", MBGC_Params::DEFAULT_NO_OF_THREADS);
                fprintf(stderr, "-d decompression mode\n"
                                "-f decompress files with names containing the given pattern\n\n");

                fprintf(stderr, "------------------ EXPERT OPTIONS ----------------\n");
                fprintf(stderr, "[-k matchingKmerLength] (24 <= k <= 40, 32 - default)\n"
                                "[-s referenceSamplingStep] (s > 0, 16 - default)\n"
                                "[-m skipMargin] (0 <= m <= 255, 16 - default)\n"
                                "[-o referenceFactorBinaryOrder] (0 <= o <= 12, adjusted by default)\n"
                                "[-c compression level] (1 - fast; 2 - default; 3 - max)\n"
                                "[-u unmatchedFractionFactor] (1 <= u <= 255, 192 - default)\n"
                                "[-w referenceSlidingWindowFactor] (2 <= w <= 255, 16 - default)\n"
                                "[-b] brute parallel encoding mode\n"
                                "[-p] disable matching parallelization (only I/O and backend compression)\n\n");
#ifdef DEVELOPER_BUILD
                fprintf(stderr, "------------------ DEVELOPER OPTIONS ----------------\n");
                fprintf(stderr, "[-R] no reverse-compliments in reference\n");
                fprintf(stderr, "[-L] literals ref extension strategy (experimental - only for encoding)\n");
                fprintf(stderr, "[-C] combined ref extension strategy\n");
                fprintf(stderr, "[-B] split contigs into blocks (experimental - only for encoding)\n");
                fprintf(stderr, "[-U] dynamic unmatched fraction ref extension strategy\n");
                fprintf(stderr, "[-I] interleave order of files \n\n");
                fprintf(stderr, "[-F] do not format decompressed DNA (i.e., 80 bases per row limit)\n");
                fprintf(stderr, "[-v] validation mode (disables decompression)\n");
                fprintf(stderr, "[-D] dump streams to files during decompression\n\n");
#endif
                fprintf(stderr, "The order of all selected options is arbitrary.\n\n");
                exit(EXIT_FAILURE);
        }
    }
    if (decompressMode && (optind < (argc - 2) || optind > (argc - 1))) {
        fprintf(stderr, "%s: For decompression expected 1 or 2 arguments after options (found %d)\n", argv[0],
                argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);

    }
    if (!decompressMode && (optind > (argc - 2) || optind < (argc - 2))) {
        fprintf(stderr, "%s: For compression expected 2 arguments after options (found %d)\n", argv[0], argc - optind);
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (decompressMode && compressionParamPresent) {
        fprintf(stderr, "Cannot use compression options in decompression mode.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (!params.filterPattern.empty() && compressionParamPresent) {
        fprintf(stderr, "Cannot use file filter pattern in compression mode.\n");
        fprintf(stderr, "try '%s -?' for more information\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    if (PgHelpers::numberOfThreads <= 0) {
        fprintf(stderr, "The number of threads must be positive.\n");
        exit(EXIT_FAILURE);
    }
    omp_set_num_threads(PgHelpers::numberOfThreads);
#pragma omp parallel
#pragma omp single
    *PgHelpers::logout << "number of threads: " << omp_get_num_threads() << endl;

    PgHelpers::time_checkpoint();

    if (!decompressMode) {
        params.setSeqListFileName(argv[optind++]);
        params.setArchiveFileName(argv[optind++]);
        MBGC_Encoder encoder(&params);
        encoder.encode();
        cout << "encoding time - " << PgHelpers::time_millis() << " [ms]" << endl;
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
        MBGC_Decoder decoder(&params);
        decoder.decode();
        cout << (decompressMode?"decoding":"total") << " time - " << PgHelpers::time_millis() << " [ms]" << endl;
    }
    return EXIT_SUCCESS;
}

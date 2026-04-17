#include "MultipleGenomeMatchingProcessor.h"

#include "../utils/kseq.h"
#include "input_with_libdeflate_wrapper.h"
#include <fstream>

KSEQ_INIT(mgmpInFile, mgmpInRead)

#define KSEQ_READ(seq) (params->enableDNALineLengthDetection ? \
    ( params->allowLossyParsing ? kseq_read_lossy(seq) : kseq_read_lossless_fasta(seq) ): kseq_read(seq))

#define KSEQ_DNA_LINE_LENGTH(seq) (params->enableDNALineLengthDetection ? \
    ( params->allowLossyParsing ? seq->maxLastDnaLineLen \
        : (seq->dnaLineLen == DNA_NOT_WELLFORMED ? 0 : seq->dnaLineLen )) : 0)


void MultipleGenomeMatchingProcessor::validate_kseq_status(string& fileName, const char* seq_name, int kseq_status)
{
    if (kseq_status <= -2) {
        fprintf(stderr, "Error parsing file %s", fileName.c_str());
        switch (kseq_status)
        {
        case -2:
        case -3: fprintf(stderr, " - expected FASTA format.\n");
            break;
        case -4: fprintf(stderr, "\nDetected inconsistent line length in sequences while processing:\n>%s\n",
            seq_name);
            break;
        default: fprintf(stderr, " - unknown error (code %d).\n", kseq_status);
        }
        print_invalid_kseq_status_message(kseq_status);
        exit(EXIT_FAILURE);
    }
}

#if !defined(__arm__) && !defined(__aarch64__) && !defined(__ARM_ARCH)
#include "../libs/asmlib.h"
#endif
#include <omp.h>

MultipleGenomeMatchingProcessor::MultipleGenomeMatchingProcessor(MGMP_Params *mgmpParams):
    params(mgmpParams) { }

string readHeader(const kseq_t *seq) {
    string headersRes(seq->name.s, seq->name.l);
    if (seq->comment.l) {
        headersRes.push_back(' ');
        headersRes.append(seq->comment.s, seq->comment.l);
    }
    return headersRes;
}

mgmpInFile firstFp;
bool firstFpIsFirstTarget = false;

size_t MultipleGenomeMatchingProcessor::openFirstFile(const char* filename, bool isTarget) {
    firstFp = mgmpInOpen(filename);
    firstFpIsFirstTarget |= isTarget;
    return firstFp.fileSize;
}

void MultipleGenomeMatchingProcessor::switchToSplitIterMode() {
    mgmpInSplit_init(firstFp);
}

void MultipleGenomeMatchingProcessor::loadG0Ref(string& refName) {
    string refStr;
    int seqCounter = 0;

    firstFp = mgmpInOpen(refName.c_str());
    if (!params->sequentialMatching) {
        largestFileLength = firstFp.size;
        if (singleFastaFileMode) {
            switchToSplitIterMode();
            mgmpInSplit_next(firstFp, params->MIN_REF_INIT_SIZE, '>');
        }
        totalFilesLength += firstFp.size;
    }
    kseq_t* seq = kseq_init(firstFp);
    initStreamsForG0Ref();
    int kseq_status = 0;
    while ((kseq_status = KSEQ_READ(seq)) >= 0) {
        if (params->uppercaseDNA)
            PgHelpers::upperSequence(seq->seq.s, seq->seq.l);
        if (params->probeProteinsProfile(seq->seq.s, seq->seq.l)) {
            this->setProteinsProfile();
        }
        largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
        refStr.append(seq->seq.s, seq->seq.l);
        processG0RefContig(seq->seq.s, seq->seq.l);
        if (params->sequentialMatching)
        {
            int probed_len = seq->seq.l;
            while (probed_len < MGMP_Params::MIN_PROBE_LEN && ((kseq_status = KSEQ_READ(seq)) >= 0))
            {
                params->probeProteinsProfile(seq->seq.s, seq->seq.l);
                probed_len += seq->seq.l;
            }
            validate_kseq_status(refName, seq->name.s, kseq_status);
            break;
        }
        string header = readHeader(seq);
        processHeader(0, header);
        seqCounter++;
    }
    validate_kseq_status(refName, seq->name.s, kseq_status);
    refG0InitPos = refStr.size();
    *PgHelpers::devout << "loaded reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    size_t basicRefLength = params->sequentialMatching && !singleFastaFileMode ? firstFp.size : refStr.size();
    basicRefLength = basicRefLength > params->MIN_BASIC_BLOCK_SIZE ? basicRefLength : params->MIN_BASIC_BLOCK_SIZE;
    int64_t elementsCount = filesCount;
    if (singleFastaFileMode) {
        int64_t remaining_size = firstFp.fileSize -
            (params->sequentialMatching ? params->MIN_BASIC_BLOCK_SIZE : firstFp.size);
        elementsCount = 1 + (remaining_size + params->MIN_BASIC_BLOCK_SIZE - 1) / params->MIN_BASIC_BLOCK_SIZE;
    }
    params->g0IsTarget = params->sequentialMatching;
    targetsCount = (params->sequentialMatching && singleFastaFileMode) ? 1 : elementsCount - (params->g0IsTarget ? 0 : 1);
    if (!params->g0IsTarget && filesCount == 1 && targetsCount < MGMP_Params::SINGLEFILE_PARALLEL_MIN_TARGETS) {
        fprintf(stderr, "Switching to sequential matching mode (input file too small).\n");
        largestFileLength = 0;
        totalFilesLength = 0;
        largestContigSize = 0;
        kseq_destroy(seq);
        mgmpInClose(firstFp);
        params->setSequentialMatchingMode();
        loadG0Ref(refName);
        return;
    }
    if (params->referenceFactor < 1) {
        int tmp = 15 - (__builtin_clz(elementsCount) / 3);
        tmp = tmp < 5 ? 5 : (tmp > 12 ? 12 : tmp);
        params->referenceFactor = (size_t) 1 << tmp;
    }
    initMatcher(refStr.data(), refStr.size(), basicRefLength);
    size_t dnaLineLength = KSEQ_DNA_LINE_LENGTH(seq);
    kseq_destroy(seq);
    if (!params->sequentialMatching) {
        processTargetMeta(seqCounter, dnaLineLength);
        applyTemplateToHeaders(0);
        if (!singleFastaFileMode)
            mgmpInClose(firstFp);
    } else
        firstFpIsFirstTarget = true;

    if (basicRefLength > params->AVG_REF_INIT_SIZE * 2)
        readingBufferSize = 1 + readingBufferSize * params->AVG_REF_INIT_SIZE * 2 / basicRefLength;

    *PgHelpers::appout << "processed reference dataset - " << PgHelpers::time_millis() << " [ms]" << endl;
}

void MultipleGenomeMatchingProcessor::initMatcher(const char* refStrPtr, const size_t refStrSize, size_t basicRefLength,
    bool useOfExistingReference) {
    size_t refLengthLimit = params->referenceFactor * basicRefLength;
    if (!useOfExistingReference && !__builtin_ctz((uint32_t) params->k1))
        params->separateRCBuffer = false;
    if (!useOfExistingReference && (!params->isRCinReferenceDisabled() || (!params->separateRCBuffer)))
        refLengthLimit *= 2;
    if (refLengthLimit <= UINT32_MAX)
        params->enable40bitReference = false;
    else if (!params->enable40bitReference)
        refLengthLimit = UINT32_MAX;
    else if (refLengthLimit > MGMP_Params::REFERENCE_LENGTH_LIMIT)
        refLengthLimit = MGMP_Params::REFERENCE_LENGTH_LIMIT;
    if (refLengthLimit > UINT32_MAX)
        refLengthLimit = UINT32_MAX + (refLengthLimit - UINT32_MAX) / params->bigReferenceCompressorRatio;
    if (refLengthLimit < refStrSize)
        refLengthLimit = refStrSize;

    if (__builtin_ctz((uint32_t) params->k1)) {
        matcher = new SlidingWindowExpSparseEMMatcher(refLengthLimit, params->k, params->k1, params->k2,
                                                          params->skipMargin);
    } else {
        matcher = new SlidingWindowSparseEMMatcher(refLengthLimit, params->k, params->k1, params->k2,
                                                   params->skipMargin);
    }
    if (params->sequentialMatching)
        matcher->disableSlidingWindow();
    else
        matcher->setSlidingWindowSize(params->referenceSlidingWindowFactor);
    if (!params->circularReference)
        matcher->disableCircularBuffer();
    matcher->loadRef(refStrPtr , refStrSize, !useOfExistingReference && !params->isRCinReferenceDisabled(),
                     !useOfExistingReference && params->refRegionSeparators, MGMP_Params::REF_REGION_SEPARATOR);
#ifdef DEVELOPER_BUILD
    currentRefExtGoal = refG0InitPos * MGMP_Params::INITIAL_REF_EXT_GOAL_FACTOR;
    refExtLengthPerFile = (refLengthLimit - refStrSize - currentRefExtGoal) /
                          (filesCount > MGMP_Params::FINAL_REF_EXT_END_MARGIN_FILES ?
                           (filesCount - MGMP_Params::FINAL_REF_EXT_END_MARGIN_FILES) : 1);
#endif
    *PgHelpers::devout << "initial sparseEM reference length: " << matcher->getRefLength() << endl;
}

void MultipleGenomeMatchingProcessor::loadFileNames() {
    filesCount = fileNames.size();
    if (!filesCount) {
        fprintf(stderr, "ERROR: filelist is empty.\n");
        exit(EXIT_FAILURE);
    }
    if (params->interleaveFileOrder)
        interleaveOrderOfFiles();
    fileNamesSet.reserve(fileNamesSet.size() + filesCount);
    for(string& fileName: fileNames)
        processFileName(fileName);
    fileNames.erase(std::remove_if(fileNames.begin(), fileNames.end(), [&](const auto &name) { return name.empty(); }),
                    fileNames.end());
    filesCount = fileNames.size();
    if (!filesCount) {
        fprintf(stderr, "ERROR: no files to add to the archive.\n");
        exit(EXIT_FAILURE);
    }
}

void MultipleGenomeMatchingProcessor::interleaveOrderOfFiles() {
    vector<string> interleaveNames;
    interleaveNames.reserve(filesCount);
    int stride = sqrt(filesCount);
    if (stride == 0)
        stride = 1;
    const int group = 1;
    for(int i = 0; i < stride * group; i += group) {
        int curr = i;
        while (curr < filesCount) {
            for (int g = 0; g < group && curr + g < filesCount; g++)
                interleaveNames.push_back(fileNames[curr + g]);
            curr += stride * group;
        }
    }
    fileNames = std::move(interleaveNames);
}

void MultipleGenomeMatchingProcessor::processTargetsWithParallelIO() {
    initProcessTargetsWithParallelIO();
    int l;
    mgmpInFile fp;
    kseq_t *seq;
    vector<PgTools::TextMatch> resMatches;
    targetRefExtensions.resize(1);
#ifdef DEVELOPER_BUILD
    targetRefExtensions[0].reserve(refG0InitPos / 16);
#endif
    int threadsLimit = PgHelpers::numberOfThreads < MGMP_Params::PARALLELIO_MODE_THREADS_LIMIT ?
            PgHelpers::numberOfThreads : MGMP_Params::PARALLELIO_MODE_THREADS_LIMIT;
#pragma omp parallel for ordered schedule(static, 1) private(fp, seq) num_threads(threadsLimit)
    for(int i = 0; i < filesCount; i++) {
        int seqCounter = 0;
        fp = (i || !firstFpIsFirstTarget) ? mgmpInOpen(fileNames[i].c_str()) : mgmpInReset(firstFp);
        seq = kseq_init(fp);
#pragma omp ordered
        {
            const size_t startPos = matcher->getLoadedRefLength();
            unmatchedFractionFactors.push_back(params->currentUnmatchedFractionFactor < 256 ? params->currentUnmatchedFractionFactor : 0);
            unmatchedFractionFactors.push_back(params->unmatchedFractionRCFactor);
            totalFilesLength += fp.size;
            largestFileLength = fp.size > largestFileLength ? fp.size : largestFileLength;
            int kseq_status = 0;
            while ((kseq_status = KSEQ_READ(seq)) >= 0) {
                if (params->uppercaseDNA)
                    PgHelpers::upperSequence(seq->seq.s, seq->seq.l);
                targetRefExtensions[0].clear();
                largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
                string header = readHeader(seq);
                processHeader(i, header);
                seqCounter++;
                char* seqPtr = seq->seq.s;
                size_t seqLeft = seq->seq.l;
#ifdef DEVELOPER_BUILD
                size_t bSize = params->splitContigsIntoBlocks() ? params->SEQ_BLOCK_SIZE : seqLeft;
#else
                size_t bSize = seqLeft;
#endif
                while (seqLeft > 0) {
                    bSize = seqLeft < bSize ? seqLeft : bSize;
                    matcher->matchTexts(resMatches, seqPtr, bSize, false, false, params->k);
                    resCount += resMatches.size();
                    size_t currentUnmatched = processMatches(resMatches, seqPtr, bSize, 0);
                    bool loadContigToRef = params->isContigProperForRefExtension(bSize, currentUnmatched,
                                                                                 params->currentUnmatchedFractionFactor);
                    if (loadContigToRef)
                        largestRefContigSize = bSize > largestRefContigSize ? bSize : largestRefContigSize;
                    char* extPtr = loadContigToRef ? seqPtr : (char*) targetRefExtensions[0].data();
                    size_t extSize = loadContigToRef ? bSize : targetRefExtensions[0].size();
                    bool loadContigRCToRef = !params->isRCinReferenceDisabled() &&
                            params->isContigProperForRefRCExtension(bSize, currentUnmatched,
                                                                    params->unmatchedFractionRCFactor);
                    matcher->loadRef(extPtr, extSize, loadContigRCToRef,
                                     params->refRegionSeparators, MGMP_Params::REF_REGION_SEPARATOR);
#ifdef DEVELOPER_BUILD
                    currentRefExtGoal -= extSize;
#endif
                    seqPtr += bSize;
                    seqLeft -= bSize;
                }
                processAfterSequence(0);
            }
            validate_kseq_status(fileNames[i], seq->name.s, kseq_status);

            processTargetMeta(seqCounter, KSEQ_DNA_LINE_LENGTH(seq));
#ifdef DEVELOPER_BUILD
            currentRefExtGoal += refExtLengthPerFile;
            if (currentRefExtGoal > 0)
                params->relaxUnmatchedFractionFactor();
            else
                params->tightenUnmatchedFractionFactor();
#endif
            processAfterTargetWithParallelIO(startPos);
        }
        kseq_destroy(seq);
        mgmpInClose(fp);
        applyTemplateToHeaders(i);
    }
    finalizeProcessTargetsWithParallelIO();
}

vector<mgmpInFile> fps;

void MultipleGenomeMatchingProcessor::initParallelProcessing() {
    if (singleFastaFileMode)
        fileNames.resize(targetsCount + targetsShift, fileNames[0] );
    unmatchedFractionFactors.resize(2 * (targetsCountShift + targetsCount));
    targetRefExtensions.resize(targetsCount);
    seqsCounts.resize(targetsCount, 0);
    dnaLineLength.resize(targetsCount, 0);
    matchingLocksPos.resize(targetsCount, SIZE_MAX);
    fps.resize(targetsCount);
    processedTargetsCount = 0;
}

void MultipleGenomeMatchingProcessor::finalizeParallelProcessingInSingleFastaFileMode() {
    fileNames.resize(targetsCount + targetsShift);
    unmatchedFractionFactors.resize(2 * (targetsCountShift + targetsCount));
    targetRefExtensions.resize(targetsCount);

    seqsCounts.resize(targetsCount);
    dnaLineLength.resize(targetsCount);
    matchingLocksPos.resize(targetsCount);
    fps.resize(targetsCount);
}

inline void MultipleGenomeMatchingProcessor::processTarget(int i, bool calledByMaster) {
    initProcessTarget(i);
    uint32_t resCount = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize= 0;
    mgmpInFile& fp = fps[i];
    int l;
    vector<TextMatch> resMatches;
    targetRefExtensions[i].reserve(refG0InitPos / 16);
    kseq_t* seq = kseq_init(fp);
    int unmatchedFractionFactor = params->currentUnmatchedFractionFactor;
    unmatchedFractionFactors[2 * (targetsCountShift + i)] = unmatchedFractionFactor < 256 ? unmatchedFractionFactor : 0;
    unmatchedFractionFactors[2 * (targetsCountShift + i) + 1] = params->unmatchedFractionRCFactor;
#pragma omp critical(acquireWorkerMatchingLockPos)
    {
        int64_t j = i;
        while (j >= 0 && matchingLocksPos[j] == SIZE_MAX)
            matchingLocksPos[j--] = matcher->acquireWorkerMatchingLockPos();
    }
    int kseq_status = 0;
    while ((kseq_status = KSEQ_READ(seq)) >= 0) {
        if (params->uppercaseDNA)
            PgHelpers::upperSequence(seq->seq.s, seq->seq.l);
        largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
        string header = readHeader(seq);
        processHeader(i + targetsShift, header);
        seqsCounts[i]++;
        char* seqPtr = seq->seq.s;
        size_t seqLeft = seq->seq.l;
#ifdef DEVELOPER_BUILD
        size_t bSize = params->splitContigsIntoBlocks() ? params->SEQ_BLOCK_SIZE : seqLeft;
#else
        size_t bSize = seqLeft;
#endif
        while (processedTargetsCount < i - allowedTargetsOutrun)
            nanosleep(SLEEP_TIME, nullptr);
        while (seqLeft > 0) {
            size_t extRefPos = targetRefExtensions[i].size();
            bSize = seqLeft < bSize ? seqLeft : bSize;
            matcher->matchTexts(resMatches, seqPtr, bSize, false, false, params->k, matchingLocksPos[i]);
            resCount += resMatches.size();
            size_t currentUnmatched = processMatches(resMatches, seqPtr, bSize, i, matchingLocksPos[i]);
            if (currentUnmatched == PROCESSING_MATCHES_SKIPPED_DUE_TO_CONTIG_DISSIMILARITY) {
                resCount -= resMatches.size();
                resMatches.clear();
                while (processedTargetsCount < i - params->allowedTargetsOutrunForDissimilarContigs)
                    nanosleep(SLEEP_TIME, nullptr);
                continue;
            }
            if (params->isContigProperForRefExtension(bSize, currentUnmatched, unmatchedFractionFactor)) {
                largestRefContigSize = bSize > largestRefContigSize ? bSize : largestRefContigSize;
                targetRefExtensions[i].append(seqPtr, bSize);
            }
            if (!params->isRCinReferenceDisabled() && params->contigsIndivduallyReversed &&
                    params->isContigProperForRefRCExtension(bSize, currentUnmatched, params->unmatchedFractionRCFactor)) {
                targetRefExtensions[i].resize(targetRefExtensions[i].size() + bSize);
                char* rcPtr = targetRefExtensions[i].data() + targetRefExtensions[i].size() - bSize;
                PgHelpers::upperReverseComplement(rcPtr - bSize, bSize, rcPtr);
            }
            seqPtr += bSize;
            seqLeft -= bSize;
        }
        processAfterSequence(i);
    }
    processAfterTarget(i);
    validate_kseq_status(fileNames[i + targetsShift], seq->name.s, kseq_status);
    dnaLineLength[i] = KSEQ_DNA_LINE_LENGTH(seq);
    kseq_destroy(seq);
#pragma omp critical(matchingStatsBlock)
    {
        this->resCount += resCount;
        this->largestRefContigSize = this->largestRefContigSize > largestRefContigSize ?
                                     this->largestRefContigSize: largestRefContigSize;
        this->largestContigSize = this->largestContigSize > largestContigSize ?
                                  this->largestContigSize: largestContigSize;
        largestFileLength = fp.size > largestFileLength ? fp.size : largestFileLength;
        totalFilesLength += fp.size;
    }
    if (!singleFastaFileMode)
        mgmpInClose(fp);
    fileNames[i + targetsShift].clear();
    if (calledByMaster)
#pragma omp atomic
        masterTargetsStats++;
    else
#pragma omp atomic
        taskTargetsStats++;
    finalizeParallelProcessingOfTarget(calledByMaster);
    applyTemplateToHeaders(i + targetsShift);
}

omp_lock_t finalizingLock;

inline int MultipleGenomeMatchingProcessor::finalizeParallelProcessingOfTarget(bool calledByMaster) {
    int counter = 0;
    int64_t tmp = processedTargetsCount;
    if (tmp < targetsCount && fileNames[tmp + targetsShift].empty())
        if (omp_test_lock(&finalizingLock)) {
            int64_t e = processedTargetsCount;
            while (e < targetsCount && fileNames[e + targetsShift].empty()) {
                const size_t startPos = matcher->getLoadedRefLength();
                matcher->loadRef(targetRefExtensions[e].data(), targetRefExtensions[e].size(),
                                 !params->isRCinReferenceDisabled() && !params->contigsIndivduallyReversed,
                                 params->refRegionSeparators, MGMP_Params::REF_REGION_SEPARATOR);
#ifdef DEVELOPER_BUILD
                currentRefExtGoal -= targetRefExtensions[e].size();
                currentRefExtGoal += refExtLengthPerFile;
                if (currentRefExtGoal > 0)
                    params->relaxUnmatchedFractionFactor();
                else
                    params->tightenUnmatchedFractionFactor();
#endif
                targetRefExtensions[e].clear();
                targetRefExtensions[e].shrink_to_fit();
                processTargetMeta(seqsCounts[e], dnaLineLength[e]);
                finalizeParallelProcessingOfTarget(e, startPos);
                matcher->releaseWorkerMatchingLockPos(matchingLocksPos[e]);
                e++;
                counter++;
                processedTargetsCount = e;
            }
            if (calledByMaster)
                masterRefExtensionsStats += counter;
            else
                taskRefExtensionsStats += counter;
            omp_unset_lock(&finalizingLock);
        }
    return counter;
}

void MultipleGenomeMatchingProcessor::tryClaimAndProcessTarget(bool calledByMaster) {
    int i = targetsCount;
#pragma omp critical(claimTargetBlock)
    {
        if (claimedTargetsCount < targetsCount) {
            int thread_no = claimedTargetsCount % readingThreadsCount;
            if (out[thread_no] != in[thread_no]) {
                i = claimedTargetsCount++;
                out[thread_no] += readingThreadsCount;
            }
        }
    }
    if (i < targetsCount)
        processTarget(i, calledByMaster);
    else if (!finalizeParallelProcessingOfTarget(calledByMaster))
        nanosleep(SLEEP_TIME, nullptr);
}

void MultipleGenomeMatchingProcessor::readFilesParallelTask(const int thread_no) {
    if (singleFastaFileMode) {
        int64_t i = 0;
        while (claimedTargetsCount < targetsCount) {
            if (thread_no == readingThreadsCount - 1) {
                int64_t t = i % readingThreadsCount;
                if (in[t] < targetsCount && in[t] < out[t] + readingThreadsCount * readingBufferSize) {
                    fps[in[t]] = mgmpInSplit_next(firstFp, params->MIN_BASIC_BLOCK_SIZE, '>');
                    if (fps[in[t]].size == 0)
                        targetsCount = in[t];
                    else {
                        in[t] += readingThreadsCount;
                        i++;
                    }
                    continue;
                }
            }
            tryClaimAndProcessTarget(false);
        }
    } else {
        while (claimedTargetsCount < targetsCount) {
            if (in[thread_no] < targetsCount &&
                in[thread_no] < out[thread_no] + readingThreadsCount * readingBufferSize) {
                fps[in[thread_no]] = (in[thread_no] || !firstFpIsFirstTarget) ?
                                     mgmpInOpen(fileNames[in[thread_no] + targetsShift].c_str()) : mgmpInReset(firstFp);
                in[thread_no] += readingThreadsCount;
            } else
                tryClaimAndProcessTarget(false);
        }
    }
}

void MultipleGenomeMatchingProcessor::processTargetsParallel() {
    omp_init_lock(&finalizingLock);
    initParallelProcessing();
#pragma omp parallel
    {
#pragma omp single
        {
            if (params->matcherWorkingThreads > PgHelpers::numberOfThreads)
                params->matcherWorkingThreads = PgHelpers::numberOfThreads;
            readingThreadsCount = params->matcherWorkingThreads - 1;
            if (readingThreadsCount > targetsCount / readingBufferSize)
                readingThreadsCount = targetsCount ? (((int64_t) targetsCount - 1) / readingBufferSize) + 1 : 0;
            allowedTargetsOutrun = params->matcherWorkingThreads * (int) params->allowedTargetsOutrunFactor;
            out.resize(readingThreadsCount);
            in.resize(readingThreadsCount);
            claimedTargetsCount = 0;
            for (int i = 0; i < readingThreadsCount; i++) {
                out[i] = i;
                in[i] = i;
#pragma omp task
                {
                    readFilesParallelTask(i);
                }
                nanosleep(SLEEP_TIME, nullptr);
            }
            while (processedTargetsCount != targetsCount) {
                tryClaimAndProcessTarget(true);
            }
        }
    }
    if (singleFastaFileMode) {
        mgmpInClose(firstFp);
        finalizeParallelProcessingInSingleFastaFileMode();
    }
    omp_destroy_lock(&finalizingLock);
}

void MultipleGenomeMatchingProcessor::processTargetsBruteParallel() {
    initParallelProcessing();
#pragma omp parallel for ordered schedule(static, 1) num_threads(params->matcherWorkingThreads)
    for(int i = 0; i < targetsCount; i++) {
        fps[i] = mgmpInOpen(fileNames[i + targetsShift].c_str());
        processTarget(i, false);
    }
    finalizeParallelProcessingOfTarget(false);
}


void MultipleGenomeMatchingProcessor::performMatching() {
    if (params->sequentialMatching)
        processTargetsWithParallelIO();
    else if (params->bruteParallel)
        processTargetsBruteParallel();
    else if (targetsCount)
        processTargetsParallel();
    else {
        fprintf(stderr, "Error selecting processing mode (no targets for parallel matching?)!\n");
        exit(EXIT_FAILURE);
    }

    *PgHelpers::devout << endl << "processed targets count: " << targetsCount << endl;
    if (PgHelpers::numberOfThreads > 1 && !params->sequentialMatching && !params->bruteParallel) {
        *PgHelpers::devout << "targets processed by - master: " << masterTargetsStats << " / (" << readingThreadsCount << ")tasks: " <<
                           taskTargetsStats << " (check: " << (masterTargetsStats + taskTargetsStats) << ")" << endl;
        *PgHelpers::devout << "ref extensions loaded by - master: " << masterRefExtensionsStats << " / (" << readingThreadsCount << ")tasks: " <<
                           taskRefExtensionsStats << " (check: " << (masterRefExtensionsStats + taskRefExtensionsStats) << ")" << endl;
    }

    refFinalTotalLength = matcher->getRefLength();
    size_t loadedRefTotalLength = matcher->getLoadedRefLength();
    delete(matcher);
    *PgHelpers::devout << "sparseEM reference length: " << refFinalTotalLength << endl;
    *PgHelpers::devout << "loaded reference total length: " << loadedRefTotalLength << endl;
    *PgHelpers::devout << "largest reference contig length: " << largestRefContigSize << endl;
    *PgHelpers::devout << "largest contig length: " << largestContigSize << endl;
    *PgHelpers::devout << "targets dna total length: " << totalDestLenAll << endl;
    *PgHelpers::devout << "exact matches total: " << resCount << endl;
    *PgHelpers::devout << "swsMEM unmatched chars: " << unmatchedCharsAll << " (removed: " <<
                       PgHelpers::getPartShareString(totalMatchedAll, totalDestLenAll) << "; " << totalDestOverlapAll << " chars in overlapped target symbol)" << endl;
    printAdditionalMatchingStats();
    *PgHelpers::appout << "matching finished - " << PgHelpers::time_millis() << " [ms]" << endl << endl;
    *PgHelpers::logout << PgHelpers::time_millis() << "       \t";
}

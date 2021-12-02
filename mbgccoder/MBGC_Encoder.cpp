#include "MBGC_Encoder.h"

#include "../utils/kseq.h"
#include "../utils/libdeflate_wrapper.h"
#include <fstream>

KSEQ_INIT(gzFile, gzread)

#include "../libs/asmlib.h"

#include "../coders/PropsLibrary.h"
#include <omp.h>

MBGC_Encoder::MBGC_Encoder(MBGC_Params *mbgcParams): params(mbgcParams) {

}

string readHeader(const kseq_t *seq) {
    string headersRes(seq->name.s, seq->name.l);
    if (seq->comment.l) {
        headersRes.push_back(' ');
        headersRes.append(seq->comment.s, seq->comment.l);
    }
    return headersRes;
}

void MBGC_Encoder::updateHeadersTemplate(uint32_t fileIndex, const string &header) {
    string& currentHeadersTemplate = fileHeadersTemplates[fileIndex];
    if (currentHeadersTemplate.empty()) {
        currentHeadersTemplate = header;
        return;
    }
    int hPos = 0;
    int templatePos = currentHeadersTemplate[0] == MBGC_Params::MATCH_MARK ? 1 : 0;
    string newTemplate(currentHeadersTemplate, 0, templatePos);
    while (templatePos + MBGC_Params::MINIMAL_PATTERN_LENGTH <= currentHeadersTemplate.size()) {
        int patternEnd = currentHeadersTemplate.find(MBGC_Params::MATCH_MARK, templatePos);
        if (patternEnd == std::string::npos)
            patternEnd = currentHeadersTemplate.size();
        if (patternEnd - templatePos < MBGC_Params::MINIMAL_PATTERN_LENGTH)
            templatePos = patternEnd + 1;
        else {
            string part = currentHeadersTemplate.substr(templatePos + MBGC_Params::MAX_PATTERN_SHIFT,
                                                        MBGC_Params::MINIMAL_PATTERN_LENGTH - MBGC_Params::MAX_PATTERN_SHIFT);
            int tmpH = header.find(part, hPos);
            if (tmpH != std::string::npos) {
                int tmpT = templatePos + MBGC_Params::MAX_PATTERN_SHIFT;
                while (tmpH != hPos && tmpT != templatePos
                       && currentHeadersTemplate[--tmpT] == header[--tmpH]);
                templatePos = tmpT + (currentHeadersTemplate[tmpT] == header[tmpH] ? 0 : 1);
                if (tmpH != hPos && newTemplate.empty())
                    newTemplate.push_back(MBGC_Params::MATCH_MARK);
                while (tmpH != header.size() && tmpT != patternEnd
                       && currentHeadersTemplate[++tmpT] == header[++tmpH]);
                if (tmpT - templatePos >= MBGC_Params::MINIMAL_PATTERN_LENGTH) {
                    newTemplate.append(currentHeadersTemplate.substr(templatePos, tmpT - templatePos));
                    if (tmpH != header.size() || tmpT != currentHeadersTemplate.size())
                        newTemplate.push_back(MBGC_Params::MATCH_MARK);
                }
                templatePos = tmpT + 1;
                hPos = tmpH + 1;
            } else {
                templatePos += MBGC_Params::MAX_PATTERN_SHIFT;
            }
        }
    }
    if (newTemplate.back() != MBGC_Params::MATCH_MARK && hPos < header.size())
        newTemplate.push_back(MBGC_Params::MATCH_MARK);
    currentHeadersTemplate = newTemplate;
}

void MBGC_Encoder::processHeader(uint32_t fileIndex, const string& header) {
    updateHeadersTemplate(fileIndex, header);
    fileHeaders[fileIndex].append(header);
    fileHeaders[fileIndex].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
}

void MBGC_Encoder::processFileName(string &fileName) {
    if(!fileName.empty() && *fileName.rbegin() == '\r') fileName.resize(fileName.length() - 1);
    if (fileName.size() > 3 && fileName.substr(fileName.size() - 3, 3) == ".gz")
        namesStr.append(fileName, 0, fileName.length() - 3);
    else
        namesStr.append(fileName);
    namesStr.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
}

gzFile refFp;

void MBGC_Encoder::loadRef(string& refName) {
    string refStr;
    int seqCounter = 0;

    refFp = gzopen(refName.c_str());
    if (!params->sequentialMatching) {
        largestFileLength = refFp.size;
        if (params->singleFastaFileMode)
            gzsplit_iter(refFp, params->MIN_REF_INIT_SIZE, '>');
        totalFilesLength += refFp.size;
    }
    kseq_t* seq = kseq_init(refFp);
    targetLiterals.resize(1);
    fileHeaders.push_back("");
    fileHeadersTemplates.push_back("");

    while (kseq_read(seq) >= 0) {
        largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
        refStr.append(seq->seq.s, seq->seq.l);
        targetLiterals[0].append(seq->seq.s, seq->seq.l);
        targetLiterals[0].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
        string header = readHeader(seq);
        if (params->sequentialMatching)
            break;
        processHeader(0, header);
        seqCounter++;
    }
    rcStart = refStr.size();
    *PgHelpers::devout << "loaded reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    if (!params->dontUseRCinReference()) {
        refStr.resize(rcStart * 2);
        char *refPtr = (char *) refStr.data();
        PgHelpers::reverseComplement(refPtr, rcStart, refPtr + rcStart);
        *PgHelpers::devout << "reversed reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    }
    size_t basicRefLength = params->sequentialMatching && !params->singleFastaFileMode ? 2 * refFp.size : refStr.size();
    basicRefLength = basicRefLength > params->MIN_REF_INIT_SIZE * 2 ? basicRefLength : params->MIN_REF_INIT_SIZE * 2;
    int64_t elementsCount = params->singleFastaFileMode ?
            (((int64_t) refFp.fileSize) - 1) / ((int64_t) basicRefLength / 2) + 1 : filesCount;
    targetsCount = elementsCount - (params->sequentialMatching ? 0 : 1);
    if (params->referenceFactor < 1) {
        int tmp = 15 - (__builtin_clz(elementsCount) / 3) +
                (params->referenceFactor == MBGC_Params::BOOSTED_ADJUSTED_REFERENCE_FACTOR_FLAG ? 1 : 0);
        tmp = tmp < 5 ? 5 : (tmp > 12 ? 12 : tmp);
        params->referenceFactor = (size_t) 1 << tmp;
    }
    size_t refLengthLimit = params->referenceFactor * basicRefLength;
    if (refLengthLimit <= UINT32_MAX)
        params->disable40bitReference = true;
    else if (params->disable40bitReference)
        refLengthLimit = UINT32_MAX;
    else if (refLengthLimit > MBGC_Params::REFERENCE_LENGTH_LIMIT)
        refLengthLimit = MBGC_Params::REFERENCE_LENGTH_LIMIT;

    matcher = new SlidingWindowSparseEMMatcher(refLengthLimit, params->k, params->k1, params->k2,
                                               params->skipMargin);
    if (params->sequentialMatching)
        matcher->disableSlidingWindow();
    else
        matcher->setSlidingWindowSize(params->referenceSlidingWindowFactor);
    matcher->loadRef(refStr.data(), refStr.size());
    kseq_destroy(seq);
    if (!params->sequentialMatching) {
        PgHelpers::writeValue<uint32_t>(seqsCountDest, seqCounter);
        if (!params->singleFastaFileMode)
            gzclose(refFp);
    }
    *PgHelpers::devout << "initial sparseEM reference length: " << matcher->getRefLength() << endl;
#ifdef DEVELOPER_BUILD
    currentRefExtGoal = rcStart * MBGC_Params::INITIAL_REF_EXT_GOAL_FACTOR;
    refExtLengthPerFile = (refLengthLimit - refStr.size() - currentRefExtGoal) /
            (filesCount > MBGC_Params::FINAL_REF_EXT_END_MARGIN_FILES ?
                (filesCount - MBGC_Params::FINAL_REF_EXT_END_MARGIN_FILES) : 1);
#endif
    *PgHelpers::appout << "processed reference dataset - " << PgHelpers::time_millis() << " [ms]" << endl;
}

size_t unmatchedCharsAll = 0;
size_t totalMatchedAll = 0;
size_t totalDestOverlapAll = 0;
size_t totalDestLenAll = 0;

string getTotalMatchStat(size_t totalMatchLength, size_t destLen) {
    return PgHelpers::toString(totalMatchLength)
           + " (" + PgHelpers::toString(((totalMatchLength * 100.0) / destLen) + 0.0005, 3) + "% == x" +
            PgHelpers::toString((destLen * 1.0 / totalMatchLength + 0.05), 1) + ")";
}

void MBGC_Encoder::appendRef(string &refExtRes, const char *extPtr, size_t length) {
    refExtRes.append(extPtr, length);
    if (!params->dontUseRCinReference()) {
        size_t tmpSize = refExtRes.size();
        refExtRes.resize(tmpSize + length);
        PgHelpers::reverseComplement(extPtr, length, (char*) refExtRes.data() + tmpSize);
    }
}

void MBGC_Encoder::processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes) {
#ifdef DEVELOPER_BUILD
    if (params->isLiteralProperForRefExtension(destLen, length)) {
        int before = params->REF_LITERAL_BEFORE_AFTER_EXT < pos ? params->REF_LITERAL_BEFORE_AFTER_EXT : pos;
        int after = (pos + length + params->REF_LITERAL_BEFORE_AFTER_EXT < destLen) ?
                    params->REF_LITERAL_BEFORE_AFTER_EXT : destLen - pos - length;
        appendRef(refExtRes, destPtr + pos - before, length + before + after);
    }
#endif
}

size_t MBGC_Encoder::processMatches(vector<PgTools::TextMatch>& textMatches, char *destPtr, size_t destLen, int i) {
    uint32_t pos = 0;
    uint32_t unmatchedChars = 0;
    uint32_t totalDestOverlap = 0;
    uint32_t totalMatched = 0;
    for (TextMatch &match: textMatches) {
        if (match.posDestText < pos) {
            uint32_t overflow = pos - match.posDestText;
            if (overflow >= match.length) {
                totalDestOverlap += match.length;
                match.length = 0;
                continue;
            }
            totalDestOverlap += overflow;
            match.length -= overflow;
            match.posDestText += overflow;
            match.posSrcText += overflow;
        }
        if (match.length < params->k) {
            totalDestOverlap += match.length;
            continue;
        }
        totalMatched += match.length;
        uint64_t length = match.posDestText - pos;
        targetLiterals[i].append(destPtr + pos, length);
        processLiteral(destPtr, pos, length, destLen, targetRefExtensions[i]);
        unmatchedChars += length;
        targetLiterals[i].push_back(MBGC_Params::MATCH_MARK);
        PgHelpers::writeValue<uint32_t>(targetMapOffDests[i], match.posSrcText);
        if (!params->disable40bitReference)
            PgHelpers::writeValue<uint8_t>(targetMapOff5thByteDests[i], match.posSrcText >> 32);
        PgHelpers::writeUIntWordFrugal(targetMapLenDests[i], match.length);
        pos = match.endPosDestText();
    }
    uint64_t length = destLen - pos;
    targetLiterals[i].append(destPtr + pos, length);
    processLiteral(destPtr, pos, length, destLen, targetRefExtensions[i]);
    unmatchedChars += length;
    textMatches.clear();

#pragma omp atomic update
    unmatchedCharsAll += unmatchedChars;
#pragma omp atomic update
    totalMatchedAll += totalMatched;
#pragma omp atomic update
    totalDestOverlapAll += totalDestOverlap;
#pragma omp atomic update
    totalDestLenAll += destLen;
    return unmatchedChars;
}

void MBGC_Encoder::loadFileNames() {
    if (params->singleFastaFileMode) {
        fileNames.push_back(params->inputFastaFileName);
    } else {
        ifstream listSrc(params->seqListFileName, ios_base::in | ios_base::binary);
        if (listSrc.fail()) {
            fprintf(stderr, "cannot open sequences list file %s\n", params->seqListFileName.c_str());
            exit(EXIT_FAILURE);
        }
        string line;
        while (getline(listSrc, line))
            fileNames.push_back(line);
    }
    filesCount = fileNames.size();
    if (!filesCount) {
        fprintf(stderr, "ERROR: filelist is empty.\n");
        exit(EXIT_FAILURE);
    }
    if (params->interleaveFileOrder)
        interleaveOrderOfFiles();
    for(string& fileName: fileNames)
        processFileName(fileName);
}

void MBGC_Encoder::interleaveOrderOfFiles() {
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
    fileNames = move(interleaveNames);
}

void MBGC_Encoder::encodeTargetsWithParallelIO() {
    int l;
    fileHeaders.resize(filesCount);
    fileHeadersTemplates.resize(filesCount);
    gzFile fp;
    kseq_t *seq;
    vector<PgTools::TextMatch> resMatches;
    targetMapOffDests.resize(1);
    targetMapOff5thByteDests.resize(1);
    targetMapLenDests.resize(1);
    targetRefExtensions.resize(1);
#ifdef DEVELOPER_BUILD
    targetRefExtensions[0].reserve(rcStart / 16);
#endif
    int threadsLimit = PgHelpers::numberOfThreads < MBGC_Params::PARALLELIO_MODE_THREADS_LIMIT ?
            PgHelpers::numberOfThreads : MBGC_Params::PARALLELIO_MODE_THREADS_LIMIT;
#pragma omp parallel for ordered schedule(static, 1) private(fp, seq) num_threads(threadsLimit)
    for(int i = 0; i < filesCount; i++) {
        int seqCounter = 0;
        fp = i ? gzopen(fileNames[i].c_str()) : gzreset(refFp);
        seq = kseq_init(fp);
#pragma omp ordered
        {
            unmatchedFractionFactors.push_back(params->currentUnmatchedFractionFactor);
            totalFilesLength += fp.size;
            largestFileLength = fp.size > largestFileLength ? fp.size : largestFileLength;
            while ((l = kseq_read(seq)) >= 0) {
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
                    if (params->isContigProperForRefExtension(bSize, currentUnmatched,
                                                              params->currentUnmatchedFractionFactor)) {
                        largestRefContigSize = bSize > largestRefContigSize ? bSize : largestRefContigSize;
                        matcher->loadRef(seqPtr, bSize);
#ifdef DEVELOPER_BUILD
                        currentRefExtGoal -= bSize;
#endif
                        if (!params->dontUseRCinReference()) {
                            PgHelpers::reverseComplementInPlace(seqPtr, bSize);
                            matcher->loadRef(seqPtr, bSize);
#ifdef DEVELOPER_BUILD
                            currentRefExtGoal -= bSize;
#endif
                        }
                    }
#ifdef DEVELOPER_BUILD
                    else if (params->useLiteralsinReference()) {
                        largestRefContigSize = targetRefExtensions[0].size() > largestRefContigSize ?
                                               targetRefExtensions[0].size() : largestRefContigSize;
                        matcher->loadRef(targetRefExtensions[0].data(), targetRefExtensions[0].size());
                        currentRefExtGoal -= bSize;
                        targetRefExtensions[0].clear();
                    }
#endif
                    seqPtr += bSize;
                    seqLeft -= bSize;
                }
                targetLiterals[0].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
            }
            PgHelpers::writeValue<uint32_t>(seqsCountDest, seqCounter);
#ifdef DEVELOPER_BUILD
            currentRefExtGoal += refExtLengthPerFile;
            if (currentRefExtGoal > 0)
                params->relaxUnmatchedFractionFactor();
            else
                params->tightenUnmatchedFractionFactor();
#endif
        }
        size_t tmp = 0;
        locksPosStream.append((char*) &tmp, sizeof(tmp));
        kseq_destroy(seq);
        gzclose(fp);
    }
    mapOffStream = targetMapOffDests[0].str();
    targetMapOffDests.pop_back();
    mapOff5thByteStream = targetMapOff5thByteDests[0].str();
    targetMapOff5thByteDests.pop_back();
    mapLenStream = targetMapLenDests[0].str();
    targetMapLenDests.pop_back();
}

vector<gzFile> fps;

void MBGC_Encoder::initParallelEncoding() {
    fileHeaders.resize(params->singleFastaFileMode? targetsCount + 1 : filesCount);
    fileHeadersTemplates.resize(params->singleFastaFileMode? targetsCount + 1 : filesCount);
    if (params->singleFastaFileMode)
        fileNames.resize(targetsCount + 1, fileNames[0] );
    unmatchedFractionFactors.resize(targetsCount);
    targetRefExtensions.resize(targetsCount);
    targetLiterals.resize(targetsCount);
    targetMapOffDests.resize(targetsCount);
    targetMapOff5thByteDests.resize(targetsCount);
    targetMapLenDests.resize(targetsCount);
    targetSeqsCounts.resize(targetsCount, 0);
    matchingLocksPos.resize(targetsCount, SIZE_MAX);
    fps.resize(targetsCount);
    processedTargetsCount = 0;
}

void MBGC_Encoder::finalizeParallelEncodingInSingleFastaFileMode() {
    fileHeaders.resize(targetsCount + 1);
    fileHeadersTemplates.resize(targetsCount + 1);
    fileNames.resize(targetsCount + 1);
    unmatchedFractionFactors.resize(targetsCount);
    targetRefExtensions.resize(targetsCount);
    targetLiterals.resize(targetsCount);
    targetMapOffDests.resize(targetsCount);
    targetMapOff5thByteDests.resize(targetsCount);
    targetMapLenDests.resize(targetsCount);
    targetSeqsCounts.resize(targetsCount);
    matchingLocksPos.resize(targetsCount);
    fps.resize(targetsCount);
}

inline void MBGC_Encoder::encodeTargetSequence(int i) {
    uint32_t resCount = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize= 0;
    gzFile& fp = fps[i];
    int l;
    vector<TextMatch> resMatches;
    targetRefExtensions[i].reserve(rcStart / 16);
    targetLiterals[i].reserve(rcStart / 32);
    kseq_t* seq = kseq_init(fp);
    uint8_t unmatchedFractionFactor = params->currentUnmatchedFractionFactor;
    unmatchedFractionFactors[i] = unmatchedFractionFactor;
#pragma omp critical(acquireWorkerMachingLockPos)
    {
        int64_t j = i;
        while (j >= 0 && matchingLocksPos[j] == SIZE_MAX)
            matchingLocksPos[j--] = matcher->acquireWorkerMatchingLockPos();
    }
    while ((l = kseq_read(seq)) >= 0) {
        largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
        string header = readHeader(seq);
        processHeader(i + 1, header);
        targetSeqsCounts[i]++;
        char* seqPtr = seq->seq.s;
        size_t seqLeft = seq->seq.l;
#ifdef DEVELOPER_BUILD
        size_t bSize = params->splitContigsIntoBlocks() ? params->SEQ_BLOCK_SIZE : seqLeft;
#else
        size_t bSize = seqLeft;
#endif
        while (seqLeft > 0) {
            bSize = seqLeft < bSize ? seqLeft : bSize;
            matcher->matchTexts(resMatches, seqPtr, bSize, false, false, params->k, matchingLocksPos[i]);
            resCount += resMatches.size();
            size_t currentUnmatched = processMatches(resMatches, seqPtr, bSize, i);
            if (params->isContigProperForRefExtension(bSize, currentUnmatched, unmatchedFractionFactor)) {
                largestRefContigSize = bSize > largestRefContigSize ? bSize : largestRefContigSize;
                targetRefExtensions[i].append(seqPtr, bSize);
                if (!params->dontUseRCinReference()) {
                    PgHelpers::reverseComplementInPlace(seqPtr, bSize);
                    targetRefExtensions[i].append(seqPtr, bSize);
                }
            };
            seqPtr += bSize;
            seqLeft -= bSize;
        }
        targetLiterals[i].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
    }
    kseq_destroy(seq);
#pragma omp critical
    {
        this->resCount += resCount;
        this->largestRefContigSize = this->largestRefContigSize > largestRefContigSize ?
                                     this->largestRefContigSize: largestRefContigSize;
        this->largestContigSize = this->largestContigSize > largestContigSize ?
                                  this->largestContigSize: largestContigSize;
        largestFileLength = fp.size > largestFileLength ? fp.size : largestFileLength;
        totalFilesLength += fp.size;
    }
    if (!params->singleFastaFileMode)
        gzclose(fp);
    fileNames[i + 1].clear();
}

inline int MBGC_Encoder::finalizeParallelEncodingOfTarget() {
    int counter = 0;
    int64_t tmp = processedTargetsCount;
    if (tmp < targetsCount && fileNames[tmp + 1].empty())
#pragma omp critical
    {
        int64_t e = processedTargetsCount;
        processedTargetsCount = targetsCount + 1;
        while (e < targetsCount && fileNames[e + 1].empty()) {
            matcher->loadRef(targetRefExtensions[e].data(),
                             targetRefExtensions[e].size());
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
            PgHelpers::writeValue<uint32_t>(seqsCountDest, targetSeqsCounts[e]);
            if (e) {
                targetLiterals[0].append(targetLiterals[e]);
                targetLiterals[e].clear();
                targetLiterals[e].shrink_to_fit();
            }
            mapOffStream.append(targetMapOffDests[e].str());
            targetMapOffDests[e].str(""); targetMapOffDests[e].clear();
            mapOff5thByteStream.append(targetMapOff5thByteDests[e].str());
            targetMapOff5thByteDests[e].str(""); targetMapOff5thByteDests[e].clear();
            mapLenStream.append(targetMapLenDests[e].str());
            targetMapLenDests[e].str(""); targetMapLenDests[e].clear();
            locksPosStream.append((char*) &matchingLocksPos[e], sizeof(matchingLocksPos[e]));
            matcher->releaseWorkerMatchingLockPos(matchingLocksPos[e]);
            e++; counter++;
        }
        processedTargetsCount = e;
    }
    return counter;
}

void MBGC_Encoder::tryClaimAndProcessTarget(bool calledByMaster) {
    int i;
#pragma omp critical
    {
        if (claimedTargetsCount < targetsCount) {
            i = claimedTargetsCount;
            int thread_no = i % readingThreadsCount;
            if (out[thread_no] != in[thread_no]) {
                claimedTargetsCount++;
                out[thread_no] += readingThreadsCount;
            } else
                i = targetsCount;
        } else
            i = targetsCount;
    }
    if (i < targetsCount) {
        encodeTargetSequence(i);
        if (calledByMaster)
#pragma omp atomic
            masterTargetsStats++;
        else
#pragma omp atomic
            taskTargetsStats++;
    } else
        nanosleep((const struct timespec[]) {{0, 100L}}, NULL);
    int counter = finalizeParallelEncodingOfTarget();
    if (calledByMaster)
#pragma omp atomic
        masterRefExtensionsStats += counter;
    else
#pragma omp atomic
        taskRefExtensionsStats += counter;
}

void MBGC_Encoder::readFilesParallelTask(const int thread_no) {
    if (params->singleFastaFileMode) {
        int64_t i = 0;
        while (claimedTargetsCount < targetsCount) {
            if (thread_no == readingThreadsCount - 1) {
                int64_t t = i % readingThreadsCount;
                if (in[t] < targetsCount && in[t] < out[t] + readingThreadsCount * READING_BUFFER_SIZE) {
                    fps[in[t]] = gzsplit_next(refFp, params->MIN_REF_INIT_SIZE, '>');
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
                in[thread_no] < out[thread_no] + readingThreadsCount * READING_BUFFER_SIZE) {
                fps[in[thread_no]] = gzopen(fileNames[in[thread_no] + 1].c_str());
                in[thread_no] += readingThreadsCount;
            } else
                tryClaimAndProcessTarget(false);
        }
    }
}

void MBGC_Encoder::encodeTargetsParallel() {
    initParallelEncoding();
#pragma omp parallel
    {
#pragma omp single
        {
            readingThreadsCount = omp_get_num_threads() - 1;
            if (readingThreadsCount > targetsCount / READING_BUFFER_SIZE)
                readingThreadsCount = targetsCount ? (((int64_t) targetsCount - 1) / READING_BUFFER_SIZE) + 1 : 0;
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
                nanosleep((const struct timespec[]) {{0, 100L}}, NULL);
            }
            while (processedTargetsCount != targetsCount) {
                tryClaimAndProcessTarget(true);
            }
        }
    }
    if (params->singleFastaFileMode) {
        gzclose(refFp);
        finalizeParallelEncodingInSingleFastaFileMode();
    }
}

void MBGC_Encoder::encodeTargetsBruteParallel() {
    initParallelEncoding();
#pragma omp parallel for ordered schedule(static, 1)
    for(int i = 0; i < targetsCount; i++) {
        fps[i] = gzopen(fileNames[i + 1].c_str());
        encodeTargetSequence(i);
        finalizeParallelEncodingOfTarget();
    }
    finalizeParallelEncodingOfTarget();
}

void MBGC_Encoder::applyTemplatesToHeaders() {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    int64_t count = params->singleFastaFileMode && !params->sequentialMatching ? targetsCount + 1 : filesCount;
    for(uint32_t f = 0; f < count; f++) {
        vector<size_t> matchesPos;
        string& hTemplate = fileHeadersTemplates[f];
        matchesPos.clear();
        size_t matchPos = 0;
        while ((matchPos = hTemplate.find(MBGC_Params::MATCH_MARK, matchPos)) != std::string::npos)
            matchesPos.push_back(matchPos++);
        size_t headersPos = 0;
        matchesPos.push_back(hTemplate.size());
        while (headersPos < fileHeaders[f].size()) {
            int tPos = 0;
            size_t headerEnd = fileHeaders[f].find(MBGC_Params::SEQ_SEPARATOR_MARK, headersPos);
            for (int i = 0; i < matchesPos.size() - 1; i++) {
                headersPos += matchesPos[i] - tPos;
                tPos = matchesPos[i] + 1;
                int hPos = matchesPos[i + 1] - tPos == 0 ? headerEnd :
                           fileHeaders[f].find(hTemplate.substr(tPos, matchesPos[i + 1] - tPos), headersPos);
                headersStr.append(fileHeaders[f], headersPos, hPos - headersPos);
                headersPos = hPos;
                headersStr.push_back(MBGC_Params::MATCH_MARK);
            }
            headersPos = headerEnd + 1;
        }
        headersTemplates.append(hTemplate);
        headersStr.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
        headersTemplates.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
    }
    fileHeaders.clear();
    fileHeadersTemplates.clear();
    *PgHelpers::devout << "applied header templates in " << PgHelpers::time_millis(start_t) << " msec." << endl;
}

size_t MBGC_Encoder::prepareAndCompressStreams() {
    string seqsCount = seqsCountDest.str();
    seqsCountDest.clear();
    ostream* out;
    string tempArchiveFileName = params->archiveFileName + MBGC_Params::TEMPORARY_FILE_SUFFIX;
    out = compressToStdout() ? &cout : new fstream(tempArchiveFileName, ios::out | ios::binary | ios::trunc);
    writeParamsAndStats(*out);
    applyTemplatesToHeaders();
    vector<CompressionJob> cJobs;
    auto namesCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 6);
    cJobs.emplace_back("file names stream... ", namesStr, namesCoderProps.get());
    auto seqCountersCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("sequences counters stream... ", seqsCount, seqCountersCoderProps.get());
    auto headersTemplateCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    cJobs.emplace_back("sequences headers templates stream... ", headersTemplates, headersTemplateCoderProps.get());
    auto headers1stCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    auto headers2ndCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_NORMAL, 3);
    auto headersCoderProps = params->headerMaxCompression ?
            getCompoundCoderProps(headers1stCoderProps.get(), headers2ndCoderProps.get()) :
            getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 16);
    cJobs.emplace_back("sequences headers stream... ", headersStr, headersCoderProps.get());
    auto refExtFactorsProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("unmatched fraction factors stream... ", unmatchedFractionFactors.data(),
                                   unmatchedFractionFactors.size(), refExtFactorsProps.get());
    auto seqCoderProps = params->fastDecoder ?
            getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_FAST, LZMA_DATAPERIODCODE_8_t)
            : getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                   params->coderLevel == CODER_LEVEL_FAST ? 4 : 7);
    int seqBlocksCount = (params->coderLevel >= CODER_LEVEL_MAX ? 1 : 2) * (params->fastDecoder ? 2 : 1);
    ParallelBlocksCoderProps blockSeqCoderProps(seqBlocksCount, seqCoderProps.get());
    cJobs.emplace_back("reference and literals stream... ", targetLiterals[0], &blockSeqCoderProps);
    auto refLocksCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("reference loading locks position stream... ", locksPosStream, refLocksCoderProps.get());
    auto nMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, params->coderLevel, LZMA_DATAPERIODCODE_32_t);
    int mapOffBlocksCount = (params->coderLevel >= CODER_LEVEL_MAX ?
            1 : (params->coderLevel == CODER_LEVEL_NORMAL ? 2 : 4));
    ParallelBlocksCoderProps blockNMapOffCoderProps(mapOffBlocksCount, nMapOffCoderProps.get());
    cJobs.emplace_back("matches offsets stream... ", mapOffStream, &blockNMapOffCoderProps);
    auto nMapOff5thByteCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    if (refFinalTotalLength > UINT32_MAX)
        cJobs.emplace_back("matches offsets 5th byte stream... ", mapOff5thByteStream, nMapOff5thByteCoderProps.get());
    auto nMapLenCoderProps = params->fastDecoder ?
                                getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_FAST, LZMA_DATAPERIODCODE_16_t)
                                : getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                                       params->coderLevel == CODER_LEVEL_FAST? 3 : 6);
    int mapLenBlocksCount = (params->coderLevel >= CODER_LEVEL_MAX ?
            1 : (params->coderLevel == CODER_LEVEL_NORMAL ? 3 : 4));
    ParallelBlocksCoderProps blockNMapLenCoderProps(mapLenBlocksCount, nMapLenCoderProps.get());
    cJobs.emplace_back("matches lengths stream... ", mapLenStream, &blockNMapLenCoderProps);
    CompressionJob::writeCompressedCollectiveParallel(*out, cJobs);
    size_t outSize = out->tellp();
    if (!compressToStdout()) {
        delete (out);
        if (std::ifstream(params->archiveFileName))
            remove(params->archiveFileName.c_str());
        if (rename(tempArchiveFileName.c_str(), params->archiveFileName.c_str()) != 0) {
            fprintf(stderr, "Error preparing output file: %s\n", params->archiveFileName.c_str());
            exit(EXIT_FAILURE);
        };
    }
    *PgHelpers::devout << endl;
    return outSize;
}

void MBGC_Encoder::writeParamsAndStats(ostream &out) const {
    out.write(MBGC_Params::MBGC_HEADER, strlen(MBGC_Params::MBGC_HEADER));
    out.put(MBGC_Params::MBGC_VERSION_MODE);
    out.put(MBGC_Params::MBGC_VERSION_MAJOR);
    out.put(MBGC_Params::MBGC_VERSION_MINOR);
    out.put(MBGC_Params::MBGC_VERSION_REVISION);
    out.put(params->coderLevel);
    out.put(PgHelpers::numberOfThreads);
    out.put(PgHelpers::numberOfThreads == 1 ? true : params->sequentialMatching);
    out.put(params->k);
    PgHelpers::writeValue(out, params->k1, false);
    PgHelpers::writeValue(out, params->k2, false);
    PgHelpers::writeValue(out, params->skipMargin, false);
    PgHelpers::writeValue(out, params->referenceFactor, false);
    out.put(params->refExtensionStrategy);
#ifdef DEVELOPER_BUILD
    if (params->useLiteralsinReference())
        PgHelpers::writeValue(out, MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_EXT, false);
    if (params->useLiteralsinReference())
        PgHelpers::writeValue(out, MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT, false);
    if (params->usesCombinedRefExtensionStrategy())
        PgHelpers::writeValue(out, params->useLiteralsinReference() ?
        params->dynamicUnmatchedFractionFactorLimit : MBGC_Params::MINIMAL_UNMATCHED_LENGTH_FACTOR, false);
    if (params->splitContigsIntoBlocks())
        PgHelpers::writeValue(out, MBGC_Params::SEQ_BLOCK_SIZE, false);
#endif
    PgHelpers::writeValue(out, filesCount, false);
    PgHelpers::writeValue(out, totalFilesLength, false);
    PgHelpers::writeValue(out, rcStart, false);
    PgHelpers::writeValue(out, largestFileLength, false);
    PgHelpers::writeValue(out, largestContigSize, false);
    PgHelpers::writeValue(out, refFinalTotalLength, false);
}

void MBGC_Encoder::encode() {
    if (!compressToStdout() && std::ifstream(params->archiveFileName))
        fprintf(stderr, "Warning: file %s already exists\n", params->archiveFileName.data());
    loadFileNames();
    if (PgHelpers::numberOfThreads == 1)
        params->sequentialMatching = true;
    loadRef(fileNames[0]);
    if (params->sequentialMatching)
        encodeTargetsWithParallelIO();
    else if (params->bruteParallel)
        encodeTargetsBruteParallel();
    else if (targetsCount)
        encodeTargetsParallel();

    *PgHelpers::devout << endl << "encoded targets count: " << targetsCount << endl;
    if (PgHelpers::numberOfThreads > 1 && !params->sequentialMatching && !params->bruteParallel) {
        *PgHelpers::devout << "targets encoded by - master: " << masterTargetsStats << " - tasks: " <<
                           taskTargetsStats << " (check: " << (masterTargetsStats + taskTargetsStats) << ")" << endl;
        *PgHelpers::devout << "ref extensions loaded by - master: " << masterRefExtensionsStats << " - tasks: " <<
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
    *PgHelpers::devout << "final unmatched chars: " << unmatchedCharsAll << " (removed: " <<
                       getTotalMatchStat(totalMatchedAll, totalDestLenAll) << "; " << totalDestOverlapAll << " chars in overlapped target symbol)" << endl;
    *PgHelpers::appout << "matching finished - " << PgHelpers::time_millis() << " [ms]" << endl << endl;
    size_t size = prepareAndCompressStreams();
    *PgHelpers::appout << "compressed " << totalFilesLength << " bytes to " << getTotalMatchStat(size, totalFilesLength) << endl;
}

bool MBGC_Encoder::compressToStdout() {
    return params->archiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS;
}

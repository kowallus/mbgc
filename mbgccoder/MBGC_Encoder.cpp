#include "MBGC_Encoder.h"

#include "../utils/kseq.h"
#include "../utils/libdeflate_wrapper.h"

KSEQ_INIT(gzFile, gzread)

#include "../libs/asmlib.h"

#include "../coders/CodersLib.h"
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
        if (patternEnd - templatePos < MBGC_Params::MINIMAL_PATTERN_LENGTH) {
            templatePos = patternEnd + 1;
            patternEnd = currentHeadersTemplate.find(MBGC_Params::MATCH_MARK, templatePos);
            if (patternEnd == std::string::npos)
                patternEnd = currentHeadersTemplate.size();
        } else {
            string part = currentHeadersTemplate.substr(templatePos + MBGC_Params::MAX_PATTERN_SHIFT,
                                                        MBGC_Params::MINIMAL_PATTERN_LENGTH - MBGC_Params::MAX_PATTERN_SHIFT);
            int tmpH = header.find(part, hPos);
            if (tmpH != std::string::npos) {
                int tmpT = templatePos + MBGC_Params::MAX_PATTERN_SHIFT;
                while (tmpH != hPos && tmpT != templatePos
                       && currentHeadersTemplate[--tmpT] == header[--tmpH]);
                templatePos = tmpT + (currentHeadersTemplate[tmpT] == header[tmpH]?0:1);
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

void MBGC_Encoder::loadRef(string& refName) {
    int seqCounter = 0;
    gzFile fp;
    kseq_t *seq;

    fp = gzopen(refName.c_str(), "r");
    seq = kseq_init(fp);
    totalFilesLength += fp.size;
    largestFileLength = fp.size;
    fileHeaders.push_back("");
    fileHeadersTemplates.push_back("");

    while (kseq_read(seq) >= 0) {
        largestContigSize = seq->seq.l > largestContigSize?seq->seq.l:largestContigSize;
        refStr.append(seq->seq.s, seq->seq.l);
        literalStr.append(seq->seq.s, seq->seq.l);
        literalStr.push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
        string header = readHeader(seq);
        processHeader(0, header);
        seqCounter++;
    }
    PgHelpers::writeValue<uint32_t>(seqsCountDest, seqCounter);
    kseq_destroy(seq);
    gzclose(fp);

    rcStart = refStr.size();
    if (rcStart == 0) {
        fprintf(stderr, "ERROR: reference sequence is empty\n");
        exit(EXIT_FAILURE);
    }
    *PgHelpers::logout << "loaded reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    if (!params->dontUseRCinReference()) {
        refStr.resize(rcStart * 2);
        char *refPtr = (char *) refStr.data();
        PgHelpers::reverseComplement(refPtr, rcStart, refPtr + rcStart);
        *PgHelpers::logout << "reversed reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    }
    size_t refLengthLimit = params->referenceFactor * refStr.size();
    if (refLengthLimit > UINT32_MAX)
        refLengthLimit = UINT32_MAX;
    matcher = new SparseEMMatcher(refLengthLimit, params->k, false, params->k1, params->k2);
    matcher->loadSrc(refStr.data(), refStr.size());
    *PgHelpers::logout << "initial sparseEM reference length: " << matcher->getSrcLength() << endl;
    currentRefExtGoal = rcStart * MBGC_Params::INITIAL_REF_EXT_GOAL_FACTOR;
    refExtLengthPerFile = (refLengthLimit - refStr.size() - currentRefExtGoal) /
            filesCount > MBGC_Params::FINAL_REF_EXT_END_MARGIN_FILES?
                (filesCount - MBGC_Params::FINAL_REF_EXT_END_MARGIN_FILES):1;
    cout << "processed reference dataset - " << PgHelpers::time_millis() << " [ms]" << endl;
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
        int after = pos + length + params->REF_LITERAL_BEFORE_AFTER_EXT < destLen ?
                    params->REF_LITERAL_BEFORE_AFTER_EXT : destLen - pos - length;
        appendRef(refExtRes, destPtr + pos - before, length + before + after);
    }
#endif
}

size_t MBGC_Encoder::processMatches(vector<PgTools::TextMatch>& textMatches, char *destPtr, size_t destLen,
                                    string& literalRes, ostringstream& mapOffDest, ostringstream& mapLenDest,
                                    string& refExtRes) {
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
        literalRes.append(destPtr + pos, length);
        processLiteral(destPtr, pos, length, destLen, refExtRes);
        unmatchedChars += length;
        literalRes.push_back(MBGC_Params::MATCH_MARK);
        PgHelpers::writeValue<uint32_t>(mapOffDest, match.posSrcText);
        PgHelpers::writeUIntWordFrugal(mapLenDest, match.length);
        pos = match.endPosDestText();
    }
    uint64_t length = destLen - pos;
    literalRes.append(destPtr + pos, length);
    processLiteral(destPtr, pos, length, destLen, refExtRes);
    unmatchedChars += length;
    literalRes.push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
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

void MBGC_Encoder::loadFileNames(ifstream &listSrc) {
    string line;
    while (getline(listSrc, line))
        fileNames.push_back(line);
    filesCount = fileNames.size();
    if (!filesCount) {
        fprintf(stderr, "ERROR: filelist is empty.\n");
        exit(EXIT_FAILURE);
    }
    targetsCount = filesCount - 1;
    if (params->interleaveFileOrder)
        interleaveOrderOfFiles();
    for(string& fileName: fileNames)
        processFileName(fileName);
    if (params->referenceFactor < 1 ) {
        int tmp = 14 - (__builtin_clz(filesCount) / 3) + (MBGC_Params::BOOSTED_ADJUSTED_REFERENCE_FACTOR_FLAG ? 1 : 0);
        tmp = tmp < 4 ? 4 : (tmp > 12 ? 12 : tmp);
        params->referenceFactor = (size_t) 1 << tmp;
    }
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
    ostringstream mapOffDest;
    ostringstream mapLenDest;
    string targetRefExtension;
#ifdef DEVELOPER_BUILD
    targetRefExtension.reserve(rcStart / 16);
#endif
#pragma omp parallel for ordered schedule(static, 1) private(fp, seq)
    for(int i = 0; i < targetsCount; i++) {
        int seqCounter = 0;
        fp = gzopen(fileNames[i + 1].c_str(), "r");
        seq = kseq_init(fp);
#pragma omp ordered
        {
            unmatchedFractionFactors.push_back(params->currentUnmatchedFractionFactor);
            totalFilesLength += fp.size;
            largestFileLength = fp.size > largestFileLength?fp.size:largestFileLength;
            while ((l = kseq_read(seq)) >= 0) {
                largestContigSize = seq->seq.l > largestContigSize?seq->seq.l:largestContigSize;
                string header = readHeader(seq);
                processHeader(i + 1, header);
                seqCounter++;
                char* seqPtr = seq->seq.s;
                size_t seqLeft = seq->seq.l;
                size_t bSize = params->splitContigsIntoBlocks() ? params->SEQ_BLOCK_SIZE : seqLeft ;
                while (seqLeft > 0) {
                    bSize = seqLeft < bSize ? seqLeft : bSize;
                    matcher->matchTexts(resMatches, seqPtr, bSize, false, false, params->k);
                    resCount += resMatches.size();
                    size_t currentUnmatched = processMatches(resMatches, seqPtr, bSize,
                                                             literalStr, mapOffDest, mapLenDest, targetRefExtension);
                    if (params->isContigProperForRefExtension(bSize, currentUnmatched,
                                                              params->currentUnmatchedFractionFactor)) {
                        largestRefContigSize = bSize > largestRefContigSize ? bSize : largestRefContigSize;
                        matcher->loadSrc(seqPtr, bSize);
                        currentRefExtGoal -= bSize;
                        if (!params->dontUseRCinReference()) {
                            PgHelpers::reverseComplementInPlace(seqPtr, bSize);
                            matcher->loadSrc(seqPtr, bSize);
                            currentRefExtGoal -= bSize;
                        }
                    } else if (params->useLiteralsinReference()) {
                        largestRefContigSize = targetRefExtension.size() > largestRefContigSize ?
                                               targetRefExtension.size() : largestRefContigSize;
                        matcher->loadSrc(targetRefExtension.data(), targetRefExtension.size());
                        currentRefExtGoal -= bSize;
                        targetRefExtension.clear();
                    }
                    seqPtr += bSize;
                    seqLeft -= bSize;
                }
                if (params->splitContigsIntoBlocks() && bSize == params->SEQ_BLOCK_SIZE)
                    literalStr.push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
            }
            PgHelpers::writeValue<uint32_t>(seqsCountDest, seqCounter);
            currentRefExtGoal += refExtLengthPerFile;
            if (currentRefExtGoal > 0)
                params->relaxUnmatchedFractionFactor();
            else
                params->tightenUnmatchedFractionFactor();
        }
        kseq_destroy(seq);
        gzclose(fp);
    }
    mapOff = mapOffDest.str();
    mapLen = mapLenDest.str();
}

vector<gzFile> fps;

void MBGC_Encoder::initParallelEncoding() {
    fileHeaders.resize(filesCount);
    fileHeadersTemplates.resize(filesCount);
    unmatchedFractionFactors.resize(targetsCount);
    targetRefExtensions.resize(targetsCount);
    targetLiterals.resize(targetsCount);
    targetMapOffDests.resize(targetsCount);
    targetMapLenDests.resize(targetsCount);
    targetSeqsCounts.resize(targetsCount, 0);
    fps.resize(targetsCount);
    processedTargetsCount = 0;
}

void MBGC_Encoder::finalizeRefExtensionsOfTarget(int i) {
    fileNames[i + 1].clear();
    if (processedTargetsCount < targetsCount && fileNames[processedTargetsCount + 1].empty())
#pragma omp critical
    {
        int e = processedTargetsCount;
        processedTargetsCount = targetsCount;
        while (e < targetsCount && fileNames[e + 1].empty()) {
            matcher->loadSrc(targetRefExtensions[e].data(),
                             targetRefExtensions[e].size());
            currentRefExtGoal -= targetRefExtensions[e].size();
            currentRefExtGoal += refExtLengthPerFile;
            if (currentRefExtGoal > 0)
                params->relaxUnmatchedFractionFactor();
            else
                params->tightenUnmatchedFractionFactor();
            targetRefExtensions[e].clear();
            targetRefExtensions[e++].shrink_to_fit();
        }
        processedTargetsCount = e;
    }
}

void MBGC_Encoder::finalizeParallelEncoding() {
    for(int i = 0; i < targetsCount; i++) {
        PgHelpers::writeValue<uint32_t>(seqsCountDest, targetSeqsCounts[i]);
        literalStr.append(targetLiterals[i]);
        mapOff.append(targetMapOffDests[i].str());
        mapLen.append(targetMapLenDests[i].str());
    }
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
    while ((l = kseq_read(seq)) >= 0) {
        largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
        string header = readHeader(seq);
        processHeader(i + 1, header);
        targetSeqsCounts[i]++;
        char* seqPtr = seq->seq.s;
        size_t seqLeft = seq->seq.l;
        size_t bSize = params->splitContigsIntoBlocks() ? params->SEQ_BLOCK_SIZE : seqLeft ;
        while (seqLeft > 0) {
            bSize = seqLeft < bSize ? seqLeft : bSize;
            matcher->matchTexts(resMatches, seqPtr, bSize, false, false, params->k);
            resCount += resMatches.size();
            size_t currentUnmatched = processMatches(resMatches, seqPtr, bSize,
                                                     targetLiterals[i], targetMapOffDests[i], targetMapLenDests[i],
                                                     targetRefExtensions[i]);
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
        if (params->splitContigsIntoBlocks() && bSize == params->SEQ_BLOCK_SIZE)
            literalStr.push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
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
    gzclose(fp);
    fileNames[i + 1].clear();
}

inline int MBGC_Encoder::finalizeParallelEncodingOfTarget() {
    int counter = 0;
    if (processedTargetsCount < targetsCount && fileNames[processedTargetsCount + 1].empty())
#pragma omp critical
    {
        int e = processedTargetsCount;
        processedTargetsCount = targetsCount + 1;
        while (e < targetsCount && fileNames[e + 1].empty()) {
            matcher->loadSrc(targetRefExtensions[e].data(),
                             targetRefExtensions[e].size());
            currentRefExtGoal -= targetRefExtensions[e].size();
            currentRefExtGoal += refExtLengthPerFile;
            if (currentRefExtGoal > 0)
                params->relaxUnmatchedFractionFactor();
            else
                params->tightenUnmatchedFractionFactor();
            targetRefExtensions[e].clear();
            targetRefExtensions[e].shrink_to_fit();
            PgHelpers::writeValue<uint32_t>(seqsCountDest, targetSeqsCounts[e]);
            literalStr.append(targetLiterals[e]);
            mapOff.append(targetMapOffDests[e].str());
            mapLen.append(targetMapLenDests[e].str());
            targetLiterals[e].clear(); targetLiterals[e].shrink_to_fit();
            targetMapOffDests[e].str(""); targetMapOffDests[e].clear();
            targetMapLenDests[e].str(""); targetMapLenDests[e].clear();
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
    while(claimedTargetsCount < targetsCount) {
        if (in[thread_no] < targetsCount && in[thread_no] < out[thread_no] + readingThreadsCount * READING_BUFFER_SIZE) {
            fps[in[thread_no]] = gzopen(fileNames[in[thread_no] + 1].c_str(), "r");
            in[thread_no] += readingThreadsCount;
        } else
            tryClaimAndProcessTarget(false);
    }
}

void MBGC_Encoder::encodeTargetsParallel() {
    int l;
    initParallelEncoding();
#pragma omp parallel
    {
#pragma omp single
        {
            readingThreadsCount = omp_get_num_threads() - 1;
            if (readingThreadsCount > targetsCount / READING_BUFFER_SIZE)
                readingThreadsCount = targetsCount ? ((targetsCount - 1) / READING_BUFFER_SIZE) + 1 : 0;
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
}

void MBGC_Encoder::encodeTargetsBruteParallel() {
    int l;
    initParallelEncoding();
#pragma omp parallel for schedule(static, 1)
    for(int i = 0; i < filesCount - 1; i++) {
        fps[i] = gzopen(fileNames[i + 1].c_str(), "r");
        encodeTargetSequence(i);
        finalizeParallelEncodingOfTarget();
    }
    finalizeParallelEncodingOfTarget();
}

void MBGC_Encoder::applyTemplatesToHeaders() {
    chrono::steady_clock::time_point start_t = chrono::steady_clock::now();
    for(uint32_t f = 0; f < filesCount; f++) {
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
    *PgHelpers::logout << "applied header templates in " << PgHelpers::time_millis(start_t) << " msec." << endl;
}

size_t MBGC_Encoder::prepareAndCompressStreams() {
    string seqsCount = seqsCountDest.str();
    seqsCountDest.clear();
    string tempArchiveFileName = params->archiveFileName + MBGC_Params::TEMPORARY_FILE_SUFFIX;
    fstream fout(tempArchiveFileName, ios::out | ios::binary | ios::trunc);
    writeParamsAndStats(fout);
    applyTemplatesToHeaders();
    vector<CompressionJob> cJobs;
    auto namesCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 6);
    cJobs.push_back(CompressionJob("file names stream... ", namesStr, namesCoderProps.get()));
    auto seqCountersCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.push_back(CompressionJob("sequences counters stream... ", seqsCount, seqCountersCoderProps.get()));
    auto headersTemplateCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    cJobs.push_back(CompressionJob("sequences headers templates stream... ", headersTemplates,
                                   headersTemplateCoderProps.get()));
    auto headersCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 16);
    cJobs.push_back(CompressionJob("sequences headers stream... ", headersStr, headersCoderProps.get()));
    auto refExtFactorsProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.push_back(CompressionJob("unmatched fraction factors stream... ", unmatchedFractionFactors.data(),
                                   unmatchedFractionFactors.size(),
                                   refExtFactorsProps.get()));
    auto seqCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                              params->coderLevel == CODER_LEVEL_FAST? 4 : 7);
    int seqBlocksCount = params->coderLevel >= CODER_LEVEL_MAX? 1 : 2;
    ParallelBlocksCoderProps blockSeqCoderProps(seqBlocksCount, seqCoderProps.get());
    cJobs.push_back(CompressionJob("reference and literals stream... ", literalStr, &blockSeqCoderProps));
    auto nMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, params->coderLevel, LZMA_DATAPERIODCODE_32_t);
    int mapOffBlocksCount = (params->coderLevel >= CODER_LEVEL_MAX? 1 : (params->coderLevel == CODER_LEVEL_NORMAL? 2 : 4));
    ParallelBlocksCoderProps blockNMapOffCoderProps(mapOffBlocksCount, nMapOffCoderProps.get());
    cJobs.push_back(CompressionJob("matches offsets stream... ", mapOff, &blockNMapOffCoderProps));
    auto nMapLenCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                                  params->coderLevel == CODER_LEVEL_FAST? 3 : 6);
    int mapLenBlocksCount = (params->coderLevel >= CODER_LEVEL_MAX? 1 : (params->coderLevel == CODER_LEVEL_NORMAL? 3 : 4));
    ParallelBlocksCoderProps blockNMapLenCoderProps(mapLenBlocksCount, nMapLenCoderProps.get());
    cJobs.push_back(CompressionJob("matches lengths stream... ", mapLen, &blockNMapLenCoderProps));
    CompressionJob::writeCompressedCollectiveParallel(fout, cJobs);
    size_t outSize = fout.tellp();
    fout.close();
    if (std::ifstream(params->archiveFileName))
        remove(params->archiveFileName.c_str());
    if (rename(tempArchiveFileName.c_str(), params->archiveFileName.c_str()) != 0) {
        fprintf(stderr, "Error preparing output file: %s\n", params->archiveFileName.c_str());
        exit(EXIT_FAILURE);
    };
    *PgHelpers::logout << endl;
    return outSize;
}

void MBGC_Encoder::writeParamsAndStats(fstream &out) const {
    out.write(params->MBGC_HEADER, strlen(params->MBGC_HEADER));
    out.put(params->MBGC_VERSION_MODE);
    out.put(params->MBGC_VERSION_MAJOR);
    out.put(params->MBGC_VERSION_MINOR);
    out.put(params->MBGC_VERSION_REVISION);
    out.put(params->coderLevel);
    out.put(PgHelpers::numberOfThreads);
    out.put(PgHelpers::numberOfThreads == 1?true:params->sequentialMatching);
    out.put(params->k);
    PgHelpers::writeValue(out, params->k1, false);
    PgHelpers::writeValue(out, params->k2, false);
    PgHelpers::writeValue(out, params->referenceFactor, false);
    out.put(params->refExtensionStrategy);
    if (params->useLiteralsinReference())
        PgHelpers::writeValue(out, params->REF_LITERAL_MINIMAL_LENGTH_EXT, false);
    if (params->useLiteralsinReference())
        PgHelpers::writeValue(out, params->REF_LITERAL_BEFORE_AFTER_EXT, false);
    if (params->usesCombinedRefExtensionStrategy())
        PgHelpers::writeValue(out, params->useLiteralsinReference() ?params->dynamicUnmatchedFractionFactorLimit :
        params->MINIMAL_UNMATCHED_LENGTH_FACTOR, false);
    if (params->splitContigsIntoBlocks())
        PgHelpers::writeValue(out, params->SEQ_BLOCK_SIZE, false);
    PgHelpers::writeValue(out, filesCount, false);
    PgHelpers::writeValue(out, totalFilesLength, false);
    PgHelpers::writeValue(out, rcStart, false);
    PgHelpers::writeValue(out, largestFileLength, false);
    PgHelpers::writeValue(out, largestContigSize, false);
    PgHelpers::writeValue(out, refTotalLength, false);
}

void MBGC_Encoder::encode() {
    if (std::ifstream(params->archiveFileName))
        fprintf(stderr, "Warning: file %s already exists\n", params->archiveFileName.data());

    ifstream listSrc(params->seqListFileName, ios_base::in | ios_base::binary);
    if (listSrc.fail()) {
        fprintf(stderr, "cannot open sequences list file %s\n", params->seqListFileName.c_str());
        exit(EXIT_FAILURE);
    }
    loadFileNames(listSrc);
    loadRef(fileNames[0]);

    if (PgHelpers::numberOfThreads == 1 || params->sequentialMatching)
        encodeTargetsWithParallelIO();
    else if (params->bruteParallel)
        encodeTargetsBruteParallel();
    else
        encodeTargetsParallel();

    *PgHelpers::logout << endl << "encoded targets count: " << (filesCount - 1) << endl;
    if (PgHelpers::numberOfThreads > 1 && !params->sequentialMatching && !params->bruteParallel) {
        *PgHelpers::logout << "targets encoded by - master: " << masterTargetsStats << " - tasks: " <<
                           taskTargetsStats << " (check: " << (masterTargetsStats + taskTargetsStats) << ")" << endl;
        *PgHelpers::logout << "ref extensions loaded by - master: " << masterRefExtensionsStats << " - tasks: " <<
                           taskRefExtensionsStats << " (check: " << (masterRefExtensionsStats + taskRefExtensionsStats) << ")" << endl;
    }

    refTotalLength = matcher->getSrcLength();
    delete(matcher);
    *PgHelpers::logout << "sparseEM reference length: " << refTotalLength << endl;
    *PgHelpers::logout << "largest reference contig length: " << largestRefContigSize << endl;
    *PgHelpers::logout << "largest contig length: " << largestContigSize << endl;
    *PgHelpers::logout << "targets dna total length: " << totalDestLenAll << endl;
    *PgHelpers::logout << "exact matches total: " << resCount << endl;
    *PgHelpers::logout << "final unmatched chars: " << unmatchedCharsAll << " (removed: " <<
                       getTotalMatchStat(totalMatchedAll, totalDestLenAll) << "; " << totalDestOverlapAll << " chars in overlapped target symbol)" << endl;
    cout << "matching finished - " << PgHelpers::time_millis() << " [ms]" << endl << endl;
    size_t size = prepareAndCompressStreams();
    cout << "compressed " << totalFilesLength << " bytes to " << getTotalMatchStat(size, totalFilesLength) << endl;
}

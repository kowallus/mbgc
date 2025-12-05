#include "MBGC_Encoder.h"

#include "../utils/kseq.h"
#include "../utils/input_with_libdeflate_wrapper.h"
#include <fstream>
#include "../matching/SimpleSequenceMatcher.h"
#include <numeric>

KSEQ_INIT(mbgcInFile, mbgcInRead)

#define KSEQ_READ(seq) (params->enableDNALineLengthEncoding ? \
    ( params->allowLossyCompression ? kseq_read_lossy(seq) : kseq_read_lossless_fasta(seq) ): kseq_read(seq))

#define KSEQ_DNA_LINE_LENGTH(seq) (params->enableDNALineLengthEncoding ? \
    ( params->allowLossyCompression ? seq->maxLastDnaLineLen \
        : (seq->dnaLineLen == DNA_NOT_WELLFORMED ? 0 : seq->dnaLineLen )) : 0)


void validate_kseq_status(string& fileName, kseq_t* seq, int kseq_status)
{
    if (kseq_status <= -2) {
        fprintf(stderr, "Error parsing file %s", fileName.c_str());
        switch (kseq_status)
        {
        case -2:
        case -3: fprintf(stderr, " - expected FASTA format.\n");
            break;
        case -4: fprintf(stderr, "\nDetected inconsistent line length in sequences while processing:\n>%s\n",
            seq->name.s);
            break;
        default: fprintf(stderr, " - unknown error (code %d).\n", kseq_status);
        }
        fprintf(stderr, "Consider using lossy compression mode (-%c option).\n", MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT);
        exit(EXIT_FAILURE);
    }
}

#if !defined(__arm__) && !defined(__aarch64__) && !defined(__ARM_ARCH)
#include "../libs/asmlib.h"
#endif
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
    switch (fileName[0]) {
    case ' ': case '\t': case '\r': case '\n':
        fprintf(stderr, "ERROR: Sequences list file %s is invalid!\n", params->seqListFileName.c_str());
        fprintf(stderr, "File name cannot start from '%c' whitespace: %s\n", fileName[0],
            fileName.c_str());
        exit(EXIT_FAILURE);
    default: ;
    }
    for (char c : fileName) {
        switch (c) {
        case '<': case '>': case ':': case '"': case '|': case '?': case '*':
            fprintf(stderr, "ERROR: Sequences list file %s is invalid!\n", params->seqListFileName.c_str());
            fprintf(stderr, "It contains invalid file name: %s\n", fileName.c_str());
            exit(EXIT_FAILURE);
        default: ;
        }
    }
    int pos, len;
    PgHelpers::normalizePath(fileName, pos, len, params->ignoreFastaFilesPath);
    if (fileNamesSet.emplace(fileName.substr(pos, len)).second == true) {
        namesStr.append(fileName, pos, len);
        namesStr.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
    } else {
        fprintf(stderr, "WARNING: skipped %s file (already in the archive).\n", fileName.substr(pos, len).c_str());
        fileName.clear();
    }
}

mbgcInFile firstFp;
bool firstFpIsFirstTarget = false;

void MBGC_Encoder::loadG0Ref(string& refName) {
    string refStr;
    int seqCounter = 0;

    firstFp = mbgcInOpen(refName.c_str());
    if (!params->sequentialMatching) {
        largestFileLength = firstFp.size;
        if (singleFastaFileMode) {
            mbgcInSplit_init(firstFp);
            mbgcInSplit_next(firstFp, params->MIN_REF_INIT_SIZE, '>');
        }
        totalFilesLength += firstFp.size;
    }
    kseq_t* seq = kseq_init(firstFp);
    targetLiterals.resize(1);
    fileHeaders.push_back("");
    fileHeadersTemplates.push_back("");
    int kseq_status = 0;
    while ((kseq_status = KSEQ_READ(seq)) >= 0) {
        if (params->uppercaseDNA)
            PgHelpers::upperSequence(seq->seq.s, seq->seq.l);
        params->probeProteinsProfile(seq->seq.s, seq->seq.l);
        largestContigSize = seq->seq.l > largestContigSize ? seq->seq.l : largestContigSize;
        refStr.append(seq->seq.s, seq->seq.l);
        targetLiterals[0].append(seq->seq.s, seq->seq.l);
        targetLiterals[0].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
        if (params->sequentialMatching)
        {
            int probed_len = seq->seq.l;
            while (probed_len < MBGC_Params::MIN_PROBE_LEN && ((kseq_status = KSEQ_READ(seq)) >= 0))
            {
                params->probeProteinsProfile(seq->seq.s, seq->seq.l);
                probed_len += seq->seq.l;
            }
            validate_kseq_status(refName, seq, kseq_status);
            break;
        }
        string header = readHeader(seq);
        processHeader(0, header);
        seqCounter++;
    }
    validate_kseq_status(refName, seq, kseq_status);
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
    if (!params->g0IsTarget && filesCount == 1 && targetsCount < MBGC_Params::SINGLEFILE_PARALLEL_MIN_TARGETS) {
        fprintf(stderr, "Switching to sequential matching mode (input file too small).\n");
        largestFileLength = 0;
        totalFilesLength = 0;
        largestContigSize = 0;
        targetLiterals.clear();
        fileHeaders.clear();
        fileHeadersTemplates.clear();
        kseq_destroy(seq);
        mbgcInClose(firstFp);
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
        PgHelpers::writeValue<uint32_t>(seqsCountsDest, seqCounter);
        PgHelpers::writeValue<uint64_t>(dnaLineLengthsDest, dnaLineLength);
        applyTemplateToHeaders(0);
        if (!singleFastaFileMode)
            mbgcInClose(firstFp);
    } else
        firstFpIsFirstTarget = true;

    if (basicRefLength > params->AVG_REF_INIT_SIZE * 2)
        readingBufferSize = 1 + readingBufferSize * params->AVG_REF_INIT_SIZE * 2 / basicRefLength;

    *PgHelpers::appout << "processed reference dataset - " << PgHelpers::time_millis() << " [ms]" << endl;
}


void MBGC_Encoder::initMatcher(const char* refStrPtr, const size_t refStrSize, size_t basicRefLength) {
    size_t refLengthLimit = params->referenceFactor * basicRefLength;
    if (!params->appendCommand && !__builtin_ctz((uint32_t) params->k1))
        params->separateRCBuffer = false;
    if (!params->appendCommand && (!params->isRCinReferenceDisabled() || (!params->separateRCBuffer)))
        refLengthLimit *= 2;
    if (refLengthLimit <= UINT32_MAX)
        params->enable40bitReference = false;
    else if (!params->enable40bitReference)
        refLengthLimit = UINT32_MAX;
    else if (refLengthLimit > MBGC_Params::REFERENCE_LENGTH_LIMIT)
        refLengthLimit = MBGC_Params::REFERENCE_LENGTH_LIMIT;
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
    matcher->loadRef(refStrPtr , refStrSize, !params->appendCommand && !params->isRCinReferenceDisabled(),
                     !params->appendCommand && params->lazyDecompressionSupport, MBGC_Params::REF_REGION_SEPARATOR);
#ifdef DEVELOPER_BUILD
    currentRefExtGoal = refG0InitPos * MBGC_Params::INITIAL_REF_EXT_GOAL_FACTOR;
    refExtLengthPerFile = (refLengthLimit - refStrSize - currentRefExtGoal) /
                          (filesCount > MBGC_Params::FINAL_REF_EXT_END_MARGIN_FILES ?
                           (filesCount - MBGC_Params::FINAL_REF_EXT_END_MARGIN_FILES) : 1);
#endif
    *PgHelpers::devout << "initial sparseEM reference length: " << matcher->getRefLength() << endl;
}

size_t unmatchedCharsAll = 0;
size_t extensionsMatchedCharsAll = 0;
size_t extensionsMismatchesAll = 0;
size_t totalMatchedAll = 0;
size_t totalDestOverlapAll = 0;
size_t totalDestLenAll = 0;
size_t removedGapBreakingMatchesAll = 0;

string getTotalMatchStat(size_t totalMatchLength, size_t destLen) {
    return PgHelpers::toString(totalMatchLength)
           + " (" + PgHelpers::toString(((totalMatchLength * 100.0) / destLen) + 0.0005, 3) + "% == x" +
            PgHelpers::toString((destLen * 1.0 / totalMatchLength + 0.05), 1) + ")";
}

void MBGC_Encoder::processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes) {
    if (params->isLiteralProperForRefExtension(length)) {
        int before = params->refLiteralBeforeAfterExt < pos ? params->refLiteralBeforeAfterExt : pos;
        int after = (pos + length + params->refLiteralBeforeAfterExt < destLen) ?
                    params->refLiteralBeforeAfterExt : destLen - pos - length;
        refExtRes.append(destPtr + pos - before, length + before + after);
    }
}

size_t MBGC_Encoder::getMatchLoadedPos(size_t pos) {
    while (pos + (matcher->getRefLength() - matcher->REF_SHIFT) < refExtLoadedPosArr.back())
        pos += (matcher->getRefLength() - matcher->REF_SHIFT);
    return pos;
}

size_t PROCESSING_MATCHES_SKIPPED_DUE_TO_CONTIG_DISSIMILARITY = SIZE_MAX;

size_t MBGC_Encoder::processMatches(vector<PgTools::TextMatch>& textMatches, char *destStart, size_t destLen, int targetIdx,
                                    size_t matchingLockPos) {
    uint32_t pos = 0;
    int64_t unmatchedChars = 0;
    uint32_t extensionsMatchedChars = 0;
    uint32_t extensionsMismatches = 0;
    uint32_t totalDestOverlap = 0;
    uint32_t totalMatched = 0;
    uint32_t removedGapBreakingMatches = 0;
    int64_t jj = 0, j;
    for (j = 0; j < textMatches.size(); j++, jj++) {
        TextMatch &match = textMatches[j];
#ifdef DEVELOPER_BUILD
        if (match.posDestText < pos) {
            fprintf(stderr, "WARNING: Unexpected matches state: overlapping matches\n");
            uint32_t overflow = pos - match.posDestText;
            if (overflow >= match.length) {
                totalDestOverlap += match.length;
                jj--;
                continue;
            }
            totalDestOverlap += overflow;
            match.shiftStartPos(overflow);
        }
        if (match.length < params->k) {
            fprintf(stderr, "WARNING: Unexpected matches state: too short match %ld < %d \n", match.length, params->k);
            totalDestOverlap += match.length;
            jj--;
            continue;
        }
#endif
        if (params->enableExtensionsWithMismatches) {
            if (j + 1 < textMatches.size()) {
                if (jj > 0 && textMatches[j + 1].pairedWith(textMatches[jj - 1]) &&
                              !match.pairedWith(textMatches[jj - 1]) &&
                              match.length < params->gapBreakingMatchMinLength) {
                    int64_t leftExtension = 0;
                    if (match.endPosDestText() == textMatches[j + 1].posDestText) {
                        const char *srcPtr = matcher->getRef() + textMatches[j + 1].posSrcText;
                        const char *ptr = destStart + textMatches[j + 1].posDestText;
                        while (*--ptr == *--srcPtr)
                            leftExtension++;
                        textMatches[j + 1].shiftStartPos(-leftExtension);
                    }
                    removedGapBreakingMatches++;
                    jj--;
                    continue;
                }
            }
        }
        totalMatched += match.length;
        uint64_t length = match.posDestText - pos;
        unmatchedChars += length;
        pos = match.endPosDestText();
        textMatches[jj] = textMatches[j];
    }
    textMatches.resize(jj);
    uint64_t literalsLeft = destLen - pos;
    unmatchedChars += literalsLeft;
    int unmatchedFractionFactor = unmatchedFractionFactors[2 * (appendedTargetsCount + targetIdx)];
    if (processedTargetsCount < targetIdx - params->allowedTargetsOutrunForDissimilarContigs &&
            params->isContigDissimilar(destLen, unmatchedChars, unmatchedFractionFactor))
        return PROCESSING_MATCHES_SKIPPED_DUE_TO_CONTIG_DISSIMILARITY;
    bool pairedGap[MBGC_Params::MAX_GAP_DEPTH] = {};
    const int NO_GAP = -1;
    int64_t gapStartIdx = NO_GAP;
    int64_t gapEndIdx = NO_GAP;
    int gapCurIdx = 0;
    bool loadLiteralsToRef = !params->isContigProperForRefExtension(destLen, unmatchedChars, unmatchedFractionFactor);
    pos = 0;
    bool isGap = false;
    for (int64_t j = 0; j < textMatches.size(); j++) {
        TextMatch &match = textMatches[j];
        literalsLeft = match.posDestText - pos;
        if (params->enableExtensionsWithMismatches) {
            if (!isGap && literalsLeft) {
                literalsLeft -= extendMatchLeft(destStart, literalsLeft, match, targetIdx, matchingLockPos, extensionsMatchedChars,
                                                extensionsMismatches);
            }
            if (j == gapEndIdx) {
                gapStartIdx = NO_GAP;
                gapEndIdx = NO_GAP;
            }
        }
        targetLiterals[targetIdx].append(destStart + pos, literalsLeft);
        if (loadLiteralsToRef) processLiteral(destStart, pos, literalsLeft, destLen, targetRefExtensions[targetIdx]);
        bool skipOffset = pairedGap[gapCurIdx];
        if (!skipOffset) {
            PgHelpers::writeValue<uint32_t>(targetMapOffDests[targetIdx], match.posSrcText);
            if (params->enable40bitReference)
                targetMapOff5thByte[targetIdx].push_back(match.posSrcText >> 32);
        }
        if (params->frugal64bitLenEncoding)
            PgHelpers::writeUInt64Frugal(targetMapLenDests[targetIdx], match.length);
        else
            PgHelpers::writeValue<uint32_t>(targetMapLenDests[targetIdx], match.length);
        targetLiterals[targetIdx].push_back(MBGC_Params::MATCH_MARK);
        pos = match.endPosDestText();

        literalsLeft = (j + 1 < textMatches.size() ? textMatches[j + 1].posDestText : destLen) - pos;
        int gapIdx = (gapCurIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;
        int gCnt = std::min((int) (textMatches.size() - j - 1), (int) params->gapDepthOffsetEncoding);
        int reduce = 0, g;
        for (g = 1; g <= gCnt; g++, gapIdx = (gapIdx + 1) % MBGC_Params::MAX_GAP_DEPTH) {
            if (pairedGap[gapIdx] || (!params->lazyDecompressionSupport && g == 1 && gapStartIdx == NO_GAP && literalsLeft == 0)) {
                reduce++;
                continue;
            }
            if (match.pairedWith(textMatches[j + g], matchingLockPos)) {
                if (params->lazyDecompressionSupport) {
                    if(!match.nextSrcRegionLoadingPos) {
                        auto itr = std::upper_bound(refExtLoadedPosArr.begin(), refExtLoadedPosArr.end(),
                                                    getMatchLoadedPos(match.posSrcText));
                        match.nextSrcRegionLoadingPos = itr == refExtLoadedPosArr.end() ? UINT64_MAX : *itr;
                    }
                    if (getMatchLoadedPos(textMatches[j + g].posSrcText) >= match.nextSrcRegionLoadingPos)
                        continue;
                    textMatches[j + g].nextSrcRegionLoadingPos = match.nextSrcRegionLoadingPos;
                }
                pairedGap[gapIdx] = true;
                if (params->enableExtensionsWithMismatches && gapEndIdx <= j + g && g <= params->gapDepthMismatchesEncoding) {
                    gapStartIdx = j;
                    gapEndIdx = j + g;
                }
                break;
            }
        }
        if (gCnt)
            targetGapDeltas[targetIdx].push_back(g <= gCnt ? g - reduce : 0);
        pairedGap[gapCurIdx] = false;
        gapCurIdx = (gapCurIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;

        bool gapStart = gapStartIdx == j;
        bool gapEnd = gapEndIdx == j + 1;
        bool gapMiddle = gapStartIdx < j && j + 1 < gapEndIdx;
        isGap = gapStart || gapMiddle || gapEnd;
        if (params->enableExtensionsWithMismatches) {
            TextMatch& corrMatch = isGap ? textMatches[gapStartIdx] : match;
            uint64_t extensionLength = extendMatchRight(destStart + pos, corrMatch, match,
                                                        literalsLeft, targetIdx, isGap, gapStart, gapMiddle, gapEnd,
                                                        extensionsMatchedChars, extensionsMismatches);
            literalsLeft -= extensionLength;
            pos += extensionLength;
        }
    }
    literalsLeft = destLen - pos;
    targetLiterals[targetIdx].append(destStart + pos, literalsLeft);
    if (loadLiteralsToRef) processLiteral(destStart, pos, literalsLeft, destLen, targetRefExtensions[targetIdx]);
    textMatches.clear();

#pragma omp atomic update
    unmatchedCharsAll += unmatchedChars;
#pragma omp atomic update
    extensionsMatchedCharsAll += extensionsMatchedChars;
#pragma omp atomic update
    extensionsMismatchesAll += extensionsMismatches;
#pragma omp atomic update
    totalMatchedAll += totalMatched;
#pragma omp atomic update
    totalDestOverlapAll += totalDestOverlap;
#pragma omp atomic update
    totalDestLenAll += destLen;
#pragma omp atomic update
    removedGapBreakingMatchesAll += removedGapBreakingMatches;
    return unmatchedChars;
}

uint64_t MBGC_Encoder::extendMatchRight(const char *gapStartPtr, const TextMatch &coreMatch, const TextMatch &match, uint64_t &length,
                                        int targetIdx, bool isGap, bool gapStart, bool gapMiddle, bool gapEnd,
                                        uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches) {
    if (length == 0) {
        if (gapMiddle)
            targetGapMismatchesFlags[targetIdx].push_back(1);
        return 0;
    }
    const char *srcStartPtr = matcher->getRef();
    const char *gapPtr = gapStartPtr;
    const char *srcPtr = matcher->getRef() + (isGap
            ? coreMatch.posSrcText + match.endPosDestText() - coreMatch.posDestText
            : match.endPosSrcText());
    size_t tmpPos = matcher->getLoadingPosition();
    const char *srcLoadingPos = matcher->getRef() + tmpPos;
    const char *srcGuard = srcPtr + length;
    const char *validSrcGuard = srcPtr + length;
    if (srcPtr == srcLoadingPos)
        validSrcGuard = srcPtr;
    if (!isGap) {
        const char *srcEnd = matcher->getRef() + matcher->getMaxRefLength();
        if (validSrcGuard > srcEnd)
            validSrcGuard = srcEnd;
        if (srcPtr <= srcLoadingPos && srcLoadingPos < validSrcGuard)
            validSrcGuard = srcLoadingPos;
    }
    if (gapStart || !isGap) {
        if (params->lazyDecompressionSupport && *srcPtr == MBGC_Params::REF_REGION_SEPARATOR)
            validSrcGuard = srcPtr;
        extensionsMismatches++;
        targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion && srcPtr < validSrcGuard ?
                                            mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
        gapPtr++;
    } else
        srcPtr--;
    int mismatches_score = params->mmsMismatchesInitialScore;
    while (++srcPtr < validSrcGuard && (!params->lazyDecompressionSupport || *srcPtr != MBGC_Params::REF_REGION_SEPARATOR) &&
        (isGap || mismatches_score < params->mmsMismatchesScoreThreshold)) {
        bool mismatch = *gapPtr != *srcPtr;
        targetGapMismatchesFlags[targetIdx].push_back(mismatch ? 1 : 0);
        if (mismatch) {
            extensionsMismatches++;
            mismatches_score += params->mmsMismatchPenalty;
            targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion ?
                    mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
        } else {
            extensionsMatchedChars++;
            mismatches_score -= params->mmsMatchBonus;
            if (mismatches_score < 0) mismatches_score = 0;
        }
        gapPtr++;
    }
    while (srcPtr++ < srcGuard && (isGap || mismatches_score < params->mmsMismatchesScoreThreshold)) {
        targetGapMismatchesFlags[targetIdx].push_back(1);
        targetLiterals[targetIdx].push_back(*gapPtr++);
        extensionsMismatches++;
        mismatches_score += params->mmsMismatchPenalty;
    }
    if ((isGap && !gapEnd) || (!isGap && mismatches_score < params->mmsMismatchesScoreThreshold))
        targetGapMismatchesFlags[targetIdx].push_back(1);
    return gapPtr - gapStartPtr;
}

uint64_t MBGC_Encoder::extendMatchLeft(const char *destStart, uint64_t length, const TextMatch &match,
                                   int targetIdx, size_t matchingLockPos,
                                   uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches) {
    const char *srcMatch = matcher->getRef() + match.posSrcText;
    const char *srcGuard = matcher->getRef() + 1;
    if (srcGuard < srcMatch - MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH)
        srcGuard = srcMatch - MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH;
    const char *srcLock = matcher->getRef() + matchingLockPos;
    if (srcGuard < srcLock && srcLock <= srcMatch)
        srcGuard = srcLock;
    bool isGuardKnownInDecompression = true;
    if (srcGuard < srcMatch - length) {
        srcGuard = srcMatch - length;
        isGuardKnownInDecompression = false;
    }
    if (srcGuard == srcMatch)
        return 0;
    const char *srcPtr = srcMatch - 1;
    const char* gapPtr = destStart + match.posDestText - 1;
    bool validSrcRegion = !params->lazyDecompressionSupport || *srcPtr != MBGC_Params::REF_REGION_SEPARATOR;
    extensionsMismatches++;
    targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion && validSrcRegion ?
            mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
    int mismatches_score = params->mmsMismatchesInitialScore;
    while (validSrcRegion && srcPtr > srcGuard && mismatches_score < params->mmsMismatchesScoreThreshold) {
        bool mismatch = *--gapPtr != *--srcPtr;
        if (params->lazyDecompressionSupport && *srcPtr == MBGC_Params::REF_REGION_SEPARATOR) {
            validSrcRegion = false;
            srcPtr++;
            gapPtr++;
            break;
        }
        targetGapMismatchesFlags[targetIdx].push_back(mismatch ? 1 : 0);
        if (mismatch) {
            mismatches_score += params->mmsMismatchPenalty;
            extensionsMismatches++;
            targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion ?
                    mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
        } else {
            mismatches_score -= params->mmsMatchBonus;
            if (mismatches_score < 0) mismatches_score = 0;
            extensionsMatchedChars++;
        }
    }
    while (!validSrcRegion && srcPtr > srcGuard && mismatches_score < params->mmsMismatchesScoreThreshold) {
        srcPtr--;
        targetGapMismatchesFlags[targetIdx].push_back(1);
        targetLiterals[targetIdx].push_back(*--gapPtr);
        extensionsMismatches++;
        mismatches_score += params->mmsMismatchPenalty;
    }
    if ((srcPtr != srcGuard || !isGuardKnownInDecompression) && mismatches_score < params->mmsMismatchesScoreThreshold)
        targetGapMismatchesFlags[targetIdx].push_back(1);
    return srcMatch - srcPtr;
}

void MBGC_Encoder::loadFileNames() {
    if (!params->repackCommand) {
        if (singleFastaFileMode) {
            fileNames.push_back(params->inputFileName);
            if (params->appendCommand) {
                filesCount = 1;
                return;
            }
        } else {
            ifstream listSrc(params->seqListFileName, ios_base::in | ios_base::binary);
            if (listSrc.fail()) {
                fprintf(stderr, "cannot open sequences list file %s\n", params->seqListFileName.c_str());
                exit(EXIT_FAILURE);
            }
            string line;
            while (getline(listSrc, line))
                fileNames.push_back(line);
            if (fileNames.empty()) {
                fprintf(stderr, "ERROR: Sequences list file %s is empty!\n", params->seqListFileName.c_str());
                exit(EXIT_FAILURE);
            }
            if (fileNames[0][0] == '>') {
                fprintf(stderr, "ERROR: Sequences list file %s is invalid!\n", params->seqListFileName.c_str());
                fprintf(stderr, "File name cannot start from '>' character: %s\n", fileNames[0].c_str());
                fprintf(stderr, "TIP: To compress a single FASTA file use -%c option.\n",
                    MBGC_Params::INPUT_FILE_OPT);
                exit(EXIT_FAILURE);
            }
        }
    }
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
    if (appendedFilesCount + filesCount == 1) {
        fprintf(stderr, "Switching to single file compression mode.\n");
        singleFastaFileMode = true;
        params->lazyDecompressionSupport = false;
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
    fileNames = std::move(interleaveNames);
}

void MBGC_Encoder::encodeTargetsWithParallelIO() {
    int l;
    fileHeaders.resize(filesCount);
    fileHeadersTemplates.resize(filesCount);
    mbgcInFile fp;
    kseq_t *seq;
    vector<PgTools::TextMatch> resMatches;
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        vector<ostringstream> &dests = p.second;
        dests.resize(1);
    }
    targetRefExtensions.resize(1);
    for(vector<string>* s: usedByteStreams)
        s->resize(1);
#ifdef DEVELOPER_BUILD
    targetRefExtensions[0].reserve(refG0InitPos / 16);
#endif
    int threadsLimit = PgHelpers::numberOfThreads < MBGC_Params::PARALLELIO_MODE_THREADS_LIMIT ?
            PgHelpers::numberOfThreads : MBGC_Params::PARALLELIO_MODE_THREADS_LIMIT;
#pragma omp parallel for ordered schedule(static, 1) private(fp, seq) num_threads(threadsLimit)
    for(int i = 0; i < filesCount; i++) {
        int seqCounter = 0;
        fp = (i || !firstFpIsFirstTarget) ? mbgcInOpen(fileNames[i].c_str()) : mbgcInReset(firstFp);
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
                                     params->lazyDecompressionSupport, MBGC_Params::REF_REGION_SEPARATOR);
#ifdef DEVELOPER_BUILD
                    currentRefExtGoal -= extSize;
#endif
                    seqPtr += bSize;
                    seqLeft -= bSize;
                }
                targetLiterals[0].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
            }
            if (params->enableExtensionsWithMismatches)
                targetGapMismatchesFlags[0].push_back(MBGC_Params::FILE_SEPARATOR_MARK);
            validate_kseq_status(fileNames[i], seq, kseq_status);

            PgHelpers::writeValue<uint32_t>(seqsCountsDest, seqCounter);
            PgHelpers::writeValue<uint64_t>(dnaLineLengthsDest, KSEQ_DNA_LINE_LENGTH(seq));
#ifdef DEVELOPER_BUILD
            currentRefExtGoal += refExtLengthPerFile;
            if (currentRefExtGoal > 0)
                params->relaxUnmatchedFractionFactor();
            else
                params->tightenUnmatchedFractionFactor();
#endif
            if (params->lazyDecompressionSupport) {
                matcher->loadSeparator(MBGC_Params::REF_REGION_SEPARATOR);
                size_t refExtSize = matcher->getLoadedRefLength() - startPos;
                PgHelpers::writeUInt64Frugal(refExtSizeDest, refExtSize);
                refExtLoadedPosArr.emplace_back(refExtLoadedPosArr.back() + refExtSize);
            }
            size_t tmp = matcher->acquireWorkerMatchingLockPos();
            locksPosStream.append((char*) &tmp, sizeof(tmp));
            matcher->releaseWorkerMatchingLockPos(tmp);
        }
        kseq_destroy(seq);
        mbgcInClose(fp);
        applyTemplateToHeaders(i);
    }
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        string &stream = p.first;
        vector<ostringstream> &dests = p.second;
        stream = dests[0].str();
        dests.pop_back();
    }
}

vector<mbgcInFile> fps;

void MBGC_Encoder::initParallelEncoding() {
    fileHeaders.resize(singleFastaFileMode ? targetsCount + targetsShift : filesCount);
    fileHeadersTemplates.resize(singleFastaFileMode ? targetsCount + targetsShift : filesCount);
    if (singleFastaFileMode)
        fileNames.resize(targetsCount + targetsShift, fileNames[0] );
    unmatchedFractionFactors.resize(2 * (appendedTargetsCount + targetsCount));
    targetRefExtensions.resize(targetsCount);
    for(vector<string>* s: usedByteStreams)
        s->resize(targetsCount);
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        vector<ostringstream> &dests = p.second;
        dests.resize(targetsCount);
    }
    seqsCounts.resize(targetsCount, 0);
    dnaLineLength.resize(targetsCount, 0);
    matchingLocksPos.resize(targetsCount, SIZE_MAX);
    fps.resize(targetsCount);
    processedTargetsCount = 0;
    if (params->lazyDecompressionSupport)
        refExtLoadedPosArr.reserve(appendedTargetsCount + targetsCount + targetsShift);
}

void MBGC_Encoder::finalizeParallelEncodingInSingleFastaFileMode() {
    fileHeaders.resize(targetsCount + targetsShift);
    fileHeadersTemplates.resize(targetsCount + targetsShift);
    fileNames.resize(targetsCount + targetsShift);
    unmatchedFractionFactors.resize(2 * (appendedTargetsCount + targetsCount));
    targetRefExtensions.resize(targetsCount);
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        vector<ostringstream> &dests = p.second;
        dests.resize(targetsCount);
    }
    for(vector<string>* s: usedByteStreams)
        s->resize(targetsCount);
    seqsCounts.resize(targetsCount);
    dnaLineLength.resize(targetsCount);
    matchingLocksPos.resize(targetsCount);
    fps.resize(targetsCount);
}

inline void MBGC_Encoder::encodeTargetSequence(int i, bool calledByMaster) {
    uint32_t resCount = 0;
    uint32_t largestRefContigSize = 0;
    uint32_t largestContigSize= 0;
    mbgcInFile& fp = fps[i];
    int l;
    vector<TextMatch> resMatches;
    targetRefExtensions[i].reserve(refG0InitPos / 16);
    targetLiterals[i].reserve(refG0InitPos / 32);
    kseq_t* seq = kseq_init(fp);
    int unmatchedFractionFactor = params->currentUnmatchedFractionFactor;
    unmatchedFractionFactors[2 * (appendedTargetsCount + i)] = unmatchedFractionFactor < 256 ? unmatchedFractionFactor : 0;
    unmatchedFractionFactors[2 * (appendedTargetsCount + i) + 1] = params->unmatchedFractionRCFactor;
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
        targetLiterals[i].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
    }
    if (params->enableExtensionsWithMismatches)
        targetGapMismatchesFlags[i].push_back(MBGC_Params::FILE_SEPARATOR_MARK);
    validate_kseq_status(fileNames[i + targetsShift], seq, kseq_status);
    dnaLineLength[i] = KSEQ_DNA_LINE_LENGTH(seq);
    kseq_destroy(seq);
#pragma omp critical(encodingStatsBlock)
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
        mbgcInClose(fp);
    fileNames[i + targetsShift].clear();
    if (calledByMaster)
#pragma omp atomic
        masterTargetsStats++;
    else
#pragma omp atomic
        taskTargetsStats++;
    finalizeParallelEncodingOfTarget(calledByMaster);
    applyTemplateToHeaders(i + targetsShift);
}

omp_lock_t finalizingLock;

inline int MBGC_Encoder::finalizeParallelEncodingOfTarget(bool calledByMaster) {
    int counter = 0;
    int64_t tmp = processedTargetsCount;
    if (tmp < targetsCount && fileNames[tmp + targetsShift].empty())
        if (omp_test_lock(&finalizingLock)) {
            int64_t e = processedTargetsCount;
            while (e < targetsCount && fileNames[e + targetsShift].empty()) {
                const size_t startPos = matcher->getLoadedRefLength();
                matcher->loadRef(targetRefExtensions[e].data(), targetRefExtensions[e].size(),
                                 !params->isRCinReferenceDisabled() && !params->contigsIndivduallyReversed,
                                 params->lazyDecompressionSupport, MBGC_Params::REF_REGION_SEPARATOR);
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
                PgHelpers::writeValue<uint32_t>(seqsCountsDest, seqsCounts[e]);
                PgHelpers::writeValue<uint64_t>(dnaLineLengthsDest, dnaLineLength[e]);
                for (vector<string> *s: usedByteStreams) {
                    if (e) {
                        (*s)[0].append((*s)[e]);
                        (*s)[e].clear();
                        (*s)[e].shrink_to_fit();
                    }
                }
                for (pair<string &, vector<ostringstream> &> p: usedStreamsWithDests) {
                    string &stream = p.first;
                    vector<ostringstream> &dests = p.second;
                    stream.append(dests[e].str());
                    dests[e].str("");
                    dests[e].clear();
                }
                if (params->lazyDecompressionSupport) {
                    matcher->loadSeparator(MBGC_Params::REF_REGION_SEPARATOR);
                    size_t refExtSize = matcher->getLoadedRefLength() - startPos;
                    PgHelpers::writeUInt64Frugal(refExtSizeDest, refExtSize);
                    refExtLoadedPosArr.emplace_back(matcher->getLoadedRefLength());
                }
                locksPosStream.append((char *) &matchingLocksPos[e], sizeof(matchingLocksPos[e]));
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

void MBGC_Encoder::tryClaimAndProcessTarget(bool calledByMaster) {
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
        encodeTargetSequence(i, calledByMaster);
    else if (!finalizeParallelEncodingOfTarget(calledByMaster))
        nanosleep(SLEEP_TIME, nullptr);
}

void MBGC_Encoder::readFilesParallelTask(const int thread_no) {
    if (singleFastaFileMode) {
        int64_t i = 0;
        while (claimedTargetsCount < targetsCount) {
            if (thread_no == readingThreadsCount - 1) {
                int64_t t = i % readingThreadsCount;
                if (in[t] < targetsCount && in[t] < out[t] + readingThreadsCount * readingBufferSize) {
                    fps[in[t]] = mbgcInSplit_next(firstFp, params->MIN_BASIC_BLOCK_SIZE, '>');
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
                                     mbgcInOpen(fileNames[in[thread_no] + targetsShift].c_str()) : mbgcInReset(firstFp);
                in[thread_no] += readingThreadsCount;
            } else
                tryClaimAndProcessTarget(false);
        }
    }
}

void MBGC_Encoder::encodeTargetsParallel() {
    omp_init_lock(&finalizingLock);
    initParallelEncoding();
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
        mbgcInClose(firstFp);
        finalizeParallelEncodingInSingleFastaFileMode();
    }
    omp_destroy_lock(&finalizingLock);
}

void MBGC_Encoder::encodeTargetsBruteParallel() {
    initParallelEncoding();
#pragma omp parallel for ordered schedule(static, 1) num_threads(params->matcherWorkingThreads)
    for(int i = 0; i < targetsCount; i++) {
        fps[i] = mbgcInOpen(fileNames[i + targetsShift].c_str());
        encodeTargetSequence(i, false);
    }
    finalizeParallelEncodingOfTarget(false);
}


void MBGC_Encoder::performMatching() {
    if (params->currentUnmatchedFractionFactor < params->unmatchedFractionRCFactor) {
        fprintf(stderr, "WARNING: %c cannot be smaller than %c. Decreased %c down to %d.\n",
                MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT, MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT,
                MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT, params->currentUnmatchedFractionFactor);
        params->unmatchedFractionRCFactor = params->currentUnmatchedFractionFactor;
    }
    if (params->matcherWorkingThreads > PgHelpers::numberOfThreads) {
        fprintf(stderr, "WARNING: %c cannot be larger than number of threads limit. Decreased %c down to %d.\n",
                MBGC_Params::MATCHER_NO_OF_THREADS_OPT,
                MBGC_Params::MATCHER_NO_OF_THREADS_OPT, PgHelpers::numberOfThreads);
        params->matcherWorkingThreads = PgHelpers::numberOfThreads;
    }
    enrollByteStream(targetLiterals);
    enrollStream(mapOffStream, targetMapOffDests);
    if (params->enable40bitReference)
        enrollByteStream(targetMapOff5thByte);
    enrollStream(mapLenStream, targetMapLenDests);
    enrollByteStream(targetGapDeltas);
    if (params->enableExtensionsWithMismatches) {
        enrollByteStream(targetGapMismatchesFlags);
        params->initMismatchesMatchingScoreParams();
    }
    if (params->lazyDecompressionSupport) {
        refExtLoadedPosArr.emplace_back(matcher->getLoadingPosition());
    }
    if (params->sequentialMatching)
        encodeTargetsWithParallelIO();
    else if (params->bruteParallel)
        encodeTargetsBruteParallel();
    else if (targetsCount)
        encodeTargetsParallel();
    else {
        fprintf(stderr, "Error selecting encoding mode (no targets for parallel matching?)!\n");
        exit(EXIT_FAILURE);
    }

    *PgHelpers::devout << endl << "encoded targets count: " << targetsCount << endl;
    if (PgHelpers::numberOfThreads > 1 && !params->sequentialMatching && !params->bruteParallel) {
        *PgHelpers::devout << "targets encoded by - master: " << masterTargetsStats << " / (" << readingThreadsCount << ")tasks: " <<
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
                       getTotalMatchStat(totalMatchedAll, totalDestLenAll) << "; " << totalDestOverlapAll << " chars in overlapped target symbol)" << endl;
    if (params->enableExtensionsWithMismatches) {
        *PgHelpers::devout << "removed matches breaking gaps total: " << removedGapBreakingMatchesAll << endl;
        *PgHelpers::devout << "extensions matched chars: " << extensionsMatchedCharsAll << endl;
        *PgHelpers::devout << "extensions mismatched chars: ~" <<
                           getTotalMatchStat(extensionsMismatchesAll,
                                             extensionsMatchedCharsAll + extensionsMismatchesAll) << endl;
        *PgHelpers::devout << "final unmatched chars: " << (unmatchedCharsAll - extensionsMatchedCharsAll)
                           << " (removed: " <<
                           getTotalMatchStat(totalMatchedAll + extensionsMatchedCharsAll, totalDestLenAll) << ")"
                           << endl;
    }
    *PgHelpers::appout << "matching finished - " << PgHelpers::time_millis() << " [ms]" << endl << endl;
    *PgHelpers::logout << PgHelpers::time_millis() << "       \t";
}

void MBGC_Encoder::applyTemplateToHeaders(uint32_t fileIdx) {
    vector<size_t> matchesPos;
    string& header = fileHeaders[fileIdx];
    string& hTemplate = fileHeadersTemplates[fileIdx];
    matchesPos.clear();
    size_t matchPos = 0;
    while ((matchPos = hTemplate.find(MBGC_Params::MATCH_MARK, matchPos)) != std::string::npos)
        matchesPos.push_back(matchPos++);
    size_t headersPos = 0;
    size_t newPos = 0;
    matchesPos.push_back(hTemplate.size());
    while (headersPos < header.size()) {
        int tPos = 0;
        size_t headerEnd = header.find(MBGC_Params::SEQ_SEPARATOR_MARK, headersPos);
        for (int i = 0; i < matchesPos.size() - 1; i++) {
            headersPos += matchesPos[i] - tPos;
            tPos = matchesPos[i] + 1;
            int hPos = matchesPos[i + 1] - tPos == 0 ? headerEnd :
                       header.find(hTemplate.substr(tPos, matchesPos[i + 1] - tPos), headersPos);
            std::copy(header.begin() + headersPos, header.begin() + hPos, header.begin() + newPos);
            newPos += hPos - headersPos;
            headersPos = hPos;
            header[newPos++] = MBGC_Params::MATCH_MARK;
        }
        headersPos = headerEnd + 1;
    }
    header.resize(newPos);
}

void MBGC_Encoder::prepareHeadersStreams() {
    int64_t count = singleFastaFileMode && !params->sequentialMatching ? targetsCount + targetsShift : filesCount;
    auto stringSize = [](int a, const std::string& b) { return a + b.size(); };
    size_t totalHeadersLength = std::accumulate(std::next(fileHeaders.begin()), fileHeaders.end(),
                                                fileHeaders[0].size(), stringSize);
    size_t totalHeaderTemplatesLength = std::accumulate(std::next(fileHeadersTemplates.begin()),
                                                        fileHeadersTemplates.end(), fileHeadersTemplates[0].size(), stringSize);
    headersStr.reserve(totalHeadersLength + count);
    headersTemplates.reserve(totalHeaderTemplatesLength + count);
    for(uint32_t f = 0; f < count; f++) {
        headersStr.append(fileHeaders[f]);
        headersTemplates.append(fileHeadersTemplates[f]);
        if (params->storeFileSeparatorMarksInHeadersStream)
            headersStr.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
        headersTemplates.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
    }
    fileHeaders.clear();
    fileHeadersTemplates.clear();
}

void MBGC_Encoder::prepareAndCompressStreams() {
    string seqsCount = seqsCountsDest.str();
    seqsCountsDest.clear();
    string dnaLineLengthStr = dnaLineLengthsDest.str();
    dnaLineLengthsDest.clear();
    ostream* out;
    string tempArchiveFileName = params->outArchiveFileName + MBGC_Params::TEMPORARY_FILE_SUFFIX;
    if (params->isStdoutMode(false)) {
        out = &cout;
    } else {
        PgHelpers::createFolders(tempArchiveFileName);
        out = new fstream(tempArchiveFileName, ios::out | ios::binary | ios::trunc);
    }
    if (!*out)
    {
        fprintf(stderr, "ERROR: Cannot create archive: %s!\n", params->outArchiveFileName.c_str());
        exit(EXIT_FAILURE);
    }
    params->write(*out);
    writeStats(*out);
    prepareHeadersStreams();
    string rcMapOff, rcMapLen;
    if (params->rcRedundancyRemoval)
        SimpleSequenceMatcher::rcMatchSequence(targetLiterals[0], rcMapOff, rcMapLen, params->rcMatchMinLength);

    params->enableOmpThreads(params->backendThreads);
    vector<CompressionJob> cJobs;
    bool fastDecoder = params->coderMode == MBGC_Params::SPEED_CODER_MODE;
    bool repoMode = params->coderMode == MBGC_Params::REPO_CODER_MODE;
    bool bestRatio = params->coderMode >= MBGC_Params::MAX_CODER_MODE;
    auto namesCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 6);
    cJobs.emplace_back("file names stream... ", namesStr, namesCoderProps.get());
    auto seqCountersCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("sequences counters stream... ", seqsCount, seqCountersCoderProps.get());
    auto headersTemplateCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    cJobs.emplace_back("sequences headers templates stream... ", headersTemplates, headersTemplateCoderProps.get());
    auto headers1stCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    auto headers2ndCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_NORMAL, 3);
    auto headersCoderProps = params->ultraStreamsCompression ?
            getCompoundCoderProps(headers1stCoderProps.get(), headers2ndCoderProps.get()) :
            getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,  fastDecoder ? 8 : 16);
    int headersBlocksCount = (bestRatio ? 1 : 2) ;
    ParallelBlocksCoderProps blockHeadersCoderProps(headersBlocksCount, headersCoderProps.get());
    cJobs.emplace_back("sequences headers stream... ", headersStr, &blockHeadersCoderProps);
    auto dnaLineLengthsProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 8);
    cJobs.emplace_back("dna line lengths stream... ", dnaLineLengthStr, dnaLineLengthsProps.get());
    auto refExtFactorsProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("unmatched fraction factors stream... ", unmatchedFractionFactors.data(),
                                   unmatchedFractionFactors.size(), refExtFactorsProps.get());
    auto seqCoderProps = params->ultraStreamsCompression ?
            getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_MAX, LZMA_DATAPERIODCODE_8_t)
            :(fastDecoder || params->k == MBGC_Params::PROTEINS_PROFILE_KMER_LENGTH ? getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_FAST, LZMA_DATAPERIODCODE_8_t)
                : getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, fastDecoder ? 4 :
                                       (params->enableExtensionsWithMismatches ?
                                           (params->mismatchesWithExclusion ? 5 : 7) :
                                    (params->sequentialMatching ? 6 : 7))));
    int seqBlocksCount = bestRatio ? 1 : (fastDecoder ? 4 : (repoMode ? 3 : 2));
    ParallelBlocksCoderProps blockSeqCoderProps(seqBlocksCount, seqCoderProps.get());
    cJobs.emplace_back("literals stream... ", targetLiterals[0], &blockSeqCoderProps);
    auto rcMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_MAX, LZMA_DATAPERIODCODE_32_t);
    auto rcMapLenCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_MAX, LZMA_DATAPERIODCODE_8_t);
    if (params->rcRedundancyRemoval) {
        cJobs.emplace_back("rc matches offsets stream... ", rcMapOff, rcMapOffCoderProps.get());
        cJobs.emplace_back("rc matches lengths stream... ", rcMapLen, rcMapLenCoderProps.get());
    }
    auto refLocksCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 8);
    cJobs.emplace_back("reference loading locks position stream... ", locksPosStream, refLocksCoderProps.get());
    auto gapDeltaCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                                   bestRatio ? 12 : (fastDecoder ? 2 : 8));
    int gapDeltaBlocksCount = repoMode ? 4 : (bestRatio ? 1 : 2);
    ParallelBlocksCoderProps blockGapDeltaCoderProps(gapDeltaBlocksCount, gapDeltaCoderProps.get());
    cJobs.emplace_back("gap between paired matches (fwd-delta) stream... ", targetGapDeltas[0], &blockGapDeltaCoderProps);
    auto gapMismatchesFlagsCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                                             fastDecoder ? 8 : 14);
    int gapMismatchesFlagsBlocksCount = fastDecoder ? 3 : 2;
    ParallelBlocksCoderProps blockGapMismatchesFlagsCoderProps(gapMismatchesFlagsBlocksCount, gapMismatchesFlagsCoderProps.get());
    uint8_t defaultCoderLevel = bestRatio ? CODER_LEVEL_MAX : (fastDecoder ? CODER_LEVEL_FAST : CODER_LEVEL_NORMAL);
    if (params->enableExtensionsWithMismatches)
        cJobs.emplace_back("mismatches in gaps between matches (flags) stream... ", targetGapMismatchesFlags[0], &blockGapMismatchesFlagsCoderProps);
    auto nMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, defaultCoderLevel, LZMA_DATAPERIODCODE_32_t);
    int mapOffBlocksCount = (bestRatio ? 1 : (fastDecoder ? 4 : 3));
    ParallelBlocksCoderProps blockNMapOffCoderProps(mapOffBlocksCount, nMapOffCoderProps.get());
    cJobs.emplace_back("matches offsets stream... ", mapOffStream, &blockNMapOffCoderProps);
    auto mapOff5thByteCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    if (refFinalTotalLength > UINT32_MAX)
        cJobs.emplace_back("matches offsets 5th byte stream... ", targetMapOff5thByte[0], mapOff5thByteCoderProps.get());
    auto nMapLenCoderProps = getDefaultCoderProps(LZMA_CODER, defaultCoderLevel,
        params->frugal64bitLenEncoding ? LZMA_DATAPERIODCODE_16_t : LZMA_DATAPERIODCODE_32_t);
    int mapLenBlocksCount = (bestRatio ? 2 : (fastDecoder ? 10 : (repoMode ? 3 : 5)));
    ParallelBlocksCoderProps blockNMapLenCoderProps(mapLenBlocksCount, nMapLenCoderProps.get());
    cJobs.emplace_back("matches lengths stream... ", mapLenStream, &blockNMapLenCoderProps);
    string refExtSizeStream = refExtSizeDest.str();
    auto refExtSizeCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 6);
    if (params->lazyDecompressionSupport)
        cJobs.emplace_back("reference extensions size stream... ", refExtSizeStream, refExtSizeCoderProps.get());
    CompressionJob::writeCompressedCollectiveParallel(*out, cJobs);
    size_t outSize = out->tellp();
    if (!params->isStdoutMode(false)) {
        delete (out);
        if (std::ifstream(params->outArchiveFileName)) {
            if (!params->appendCommand && !params->forceOverwrite) {
                fprintf(stderr, "ERROR: file %s already exists (use -%c to force overwrite)\n",
                    params->outArchiveFileName.data(), MBGC_Params::FORCE_OVERWRITE_OPT);
                exit(EXIT_FAILURE);
            }
            remove(params->outArchiveFileName.c_str());
        }
        if (rename(tempArchiveFileName.c_str(), params->outArchiveFileName.c_str()) != 0) {
            fprintf(stderr, "Error preparing output file: %s\n", params->outArchiveFileName.c_str());
            exit(EXIT_FAILURE);
        };
    }
    *PgHelpers::devout << endl;
    *PgHelpers::appout << "compressed " << totalFilesLength << " bytes to " << getTotalMatchStat(outSize, totalFilesLength) << endl;
    *PgHelpers::logout << PgHelpers::time_millis() << "       \t";
    *PgHelpers::logout << PgHelpers::toString((totalFilesLength * 1.0 / outSize + 0.05), 1) << "\t";
    *PgHelpers::logout << totalFilesLength << "\t" << outSize << "\t";
}

void MBGC_Encoder::writeStats(ostream &out) const {
    PgHelpers::writeValue(out, appendedFilesCount + filesCount);
    PgHelpers::writeValue(out, totalFilesLength);
    PgHelpers::writeValue(out, refG0InitPos);
    PgHelpers::writeValue(out, largestFileLength);
    PgHelpers::writeValue(out, largestContigSize);
    PgHelpers::writeValue(out, refFinalTotalLength);
    PgHelpers::writeValue(out, refBuffersCount);
    PgHelpers::writeValue(out, (uint64_t) targetLiterals[0].size());
}

void MBGC_Encoder::encode() {
    if (!params->noOfThreadsLimited) {
        params->coderThreads = omp_get_num_procs();
        params->backendThreads = omp_get_num_procs();
    }
    params->enableOmpThreads(params->coderThreads);
    if (!params->isStdoutMode(false) && !params->forceOverwrite
        && std::ifstream(params->outArchiveFileName)
        && (!params->appendCommand || params->inArchiveFileName != params->outArchiveFileName)) {
        fprintf(stderr, "ERROR: file %s already exists (use -%c to force overwrite)\n",
                    params->outArchiveFileName.data(), MBGC_Params::FORCE_OVERWRITE_OPT);
        exit(EXIT_FAILURE);
    }
    if (PgHelpers::numberOfThreads == 1)
        params->setSequentialMatchingMode();
    if (!params->appendCommand) {
        singleFastaFileMode = !params->inputFileName.empty();
        loadFileNames();
        loadG0Ref(fileNames[0]);
    }
    performMatching();

    prepareAndCompressStreams();
}

template<bool lazyMode>
void MBGC_Encoder::appendTemplate(MBGC_Decoder<lazyMode>* decoder, MBGC_Params& newParams) {
    decoder->decode();
    namesStr = std::move(decoder->namesStr);
    size_t namesPos = 0;
    fileNamesSet.reserve(decoder->filesCount * 2);
    for (int i = 0; i < decoder->filesCount; i++) {
        size_t tmpNamesPos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
        string name = namesStr.substr(namesPos, tmpNamesPos - namesPos);
        fileNamesSet.emplace(name);
        namesPos = tmpNamesPos + 1;
    }
    seqsCountsDest << decoder->seqsCount;
    headersTemplates = std::move(decoder->headersTemplates);
    headersStr = std::move(decoder->headersStr);
    dnaLineLengthsDest << decoder->dnaLineLengthsStr;
    unmatchedFractionFactors.assign(decoder->unmatchedFractionFactorsArray,
                                    decoder->unmatchedFractionFactorsArray + decoder->unmatchedFractionFactorsStr.size());
    params->currentUnmatchedFractionFactor = *((uint8_t*) unmatchedFractionFactors.data() + unmatchedFractionFactors.size() - 2);
    params->unmatchedFractionRCFactor = *((uint8_t*) unmatchedFractionFactors.data() + unmatchedFractionFactors.size() - 1);
    targetLiterals.push_back(std::move(decoder->literalStr));
    locksPosStream = std::move(decoder->matchingLocksPosStream);
    targetGapDeltas.push_back(std::move(decoder->gapDeltaStr));
    if (params->enableExtensionsWithMismatches)
        targetGapMismatchesFlags.push_back(std::move(decoder->gapMismatchesFlagsStr));
    targetMapOffDests.resize(1);
    targetMapOffDests[0] << decoder->mapOffStream;
    if (decoder->refTotalLength > UINT32_MAX)
        targetMapOff5thByte.push_back(std::move(decoder->mapOff5thByteStream));
    targetMapLenDests.resize(1);
    targetMapLenDests[0] << decoder->mapLenStream;
    if (params->lazyDecompressionSupport) {
        refExtSizeDest << decoder->refExtSizeStream;
        refExtLoadedPosArr.assign(decoder->refExtLoadedPosArr.begin(),
                                  decoder->refExtLoadedPosArr.end() - 1);
    }

    singleFastaFileMode = !params->inputFileName.empty();
    appendedFilesCount = singleFastaFileMode ? 0 : decoder->filesCount;
    appendedTargetsCount = decoder->targetsCount;
    totalFilesLength = decoder->totalFilesLength;
    refG0InitPos = decoder->refG0InitPos;
    largestFileLength = decoder->largestFileLength;
    largestContigSize = decoder->largestContigLength;
    refBuffersCount = decoder->refBuffersCount;

    loadFileNames();
    string refStr = std::move(decoder->refStr);
    size_t initialRefLength = refStr.size();
    firstFp = mbgcInOpen(fileNames[0].c_str());
    firstFpIsFirstTarget = true;
    if (decoder->reachedRefLengthCount == 0) {
        int64_t targetsProportion = singleFastaFileMode ?
                                (((int64_t) firstFp.fileSize) - 1) / ((int64_t) initialRefLength) + 1 :
                                filesCount / appendedFilesCount + 1;
        int tmp = 15 - (__builtin_clz(targetsProportion) / 3);
        tmp = tmp < 5 ? 5 : (tmp > 12 ? 12 : tmp);
        params->referenceFactor = (size_t) 1 << tmp;
    } else
        params->referenceFactor = 1;
    if (newParams.coderThreads > 0)
        params->setThreadsLimit(newParams.coderThreads);
    if (newParams.matcherWorkingThreads > 0)
        params->setMatcherThreads(newParams.matcherWorkingThreads);
    if (params->g0IsTarget || params->matcherWorkingThreads == 1)
        params->setSequentialMatchingMode();
    if (!params->sequentialMatching && singleFastaFileMode)
        mbgcInSplit_init(firstFp);
    size_t basicRefLength = refG0InitPos;
    basicRefLength = basicRefLength > params->MIN_BASIC_BLOCK_SIZE ? basicRefLength : params->MIN_BASIC_BLOCK_SIZE;
    targetsCount = singleFastaFileMode && !params->sequentialMatching ?
        (((int64_t) firstFp.fileSize) + basicRefLength - 1) / basicRefLength : filesCount;
    targetsShift = 0;
    initMatcher(refStr.data() + 1, refStr.size() - 1, initialRefLength);
    if (decoder->reachedRefLengthCount > 0) {
        matcher->setPosition(decoder->refPos, decoder->reachedRefLengthCount);
    }
    if (basicRefLength > params->AVG_REF_INIT_SIZE * 2)
        readingBufferSize = 1 + readingBufferSize * params->AVG_REF_INIT_SIZE * 2 / basicRefLength;
    *PgHelpers::appout << "processed reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    delete(decoder);
    encode();
}

void MBGC_Encoder::append(MBGC_Params& newParams) {
    MBGC_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(params);
    params->engageAppendCommand();
    if (params->isLazyDecompressionEnabled())
        appendTemplate<true>((MBGC_Decoder<true>*) decoder, newParams);
    else
        appendTemplate<false>((MBGC_Decoder<false>*) decoder, newParams);
}

template<bool lazyMode>
void MBGC_Encoder::repackTemplate(MBGC_Decoder<lazyMode>* decoder) {
    decoder->params->repackCommand = true;
    iterableDecoder = decoder;
    decoder->iterateInit();
    decoder->getSelectedFilesNames(fileNames);
    encode();
}

void MBGC_Encoder::repack(MBGC_Params& inParams) {
    MBGC_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&inParams);
    if (params->isLazyDecompressionEnabled())
        repackTemplate<true>((MBGC_Decoder<true>*) decoder);
    else
        repackTemplate<false>((MBGC_Decoder<false>*) decoder);
}

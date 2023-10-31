#include "MBGC_Decoder.h"

#include "../coders/CodersLib.h"
#include "../libs/asmlib.h"
#include "../matching/SimpleSequenceMatcher.h"

#include <fstream>
#include <omp.h>
#include <numeric>

#ifdef __MINGW32__
#include <fcntl.h>
#endif

#include "../utils/input_with_libdeflate_wrapper.h"

template<bool lazyMode> MBGC_Decoder<lazyMode>::MBGC_Decoder(MBGC_Params *mbgcParams, istream* inStream): MBGC_Decoder_API(), params(mbgcParams), inStream(inStream) {
    if (lazyMode && !mbgcParams->isLazyDecompressionEnabled()) {
        fprintf(stderr, "ERROR: Using unsupported lazy decompression mode\n");
        exit(EXIT_FAILURE);
    }
}

template<bool lazyMode>
MBGC_Decoder_API *MBGC_Decoder<lazyMode>::getInstance(MBGC_Params *params, string optionalWarning) {
    istream* in;
    if (params->isStdinMode(true)) {
#ifdef __MINGW32__
        if (_setmode(_fileno(stdin), _O_BINARY) == -1) {
            fprintf(stderr, "ERROR: switching cin to binary mode (errCode: %d)\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
#endif
        in = &cin;
    } else {
        in = new fstream(params->inArchiveFileName, ios_base::in | ios_base::binary);
        if (!*in) {
            fprintf(stderr, "Cannot open archive %s%s\n", params->inArchiveFileName.c_str(), optionalWarning.c_str());
            exit(EXIT_FAILURE);
        }
    }
    params->read(*in);
    if (params->isLazyDecompressionEnabled())
        return new MBGC_Decoder<true>(params, in);
    else
        return new MBGC_Decoder<false>(params, in);
}


template<bool lazyMode>
MBGC_Decoder<lazyMode>::~MBGC_Decoder() {
    if (!params->isStdinMode(true))
        delete(inStream);
}

inline void A_append(string &dest, const char *src, size_t n) {
    size_t destPos = dest.size();
    dest.resize(destPos + n);
    A_memcpy((char *) dest.data() + destPos, src, n);
}

inline void A_append(string &dest, string &src, size_t pos, size_t n) {
    size_t destPos = dest.size();
    dest.resize(destPos + n);
    A_memcpy((char *) dest.data() + destPos, src.data() + pos, n);
    //dest.append(src, pos, n);
}

inline void A_append(string &dest, string &src) {
    A_append(dest, src, 0, src.size());
    //dest.append(src);
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::writeDNA(const char *sequence, const int64_t length,
        const int64_t dnaLineLength, int threadId) {
#ifdef DEVELOPER_BUILD
    if (params->uppercaseDNA)
        PgHelpers::upperSequence((char*) sequence, length);
#endif
    uint32_t pos = 0;
    while (pos < length - dnaLineLength) {
        A_append(threadOutBuffer[threadId], sequence + pos, dnaLineLength);
        threadOutBuffer[threadId].push_back('\n');
        pos += dnaLineLength;
    }
    if (length - pos > 0) {
        A_append(threadOutBuffer[threadId], sequence + pos, length - pos);
        threadOutBuffer[threadId].push_back('\n');
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decodeHeader(string& headerTemplate, int tId) {
    int threadId = lazyMode && !sequentialDecoding ? omp_get_thread_num() : 0;
    size_t tmp, tPos = 0;
    while ((tmp = headerTemplate.find(MBGC_Params::MATCH_MARK, tPos)) != std::string::npos) {
        A_append(threadOutBuffer[threadId], headerTemplate, tPos, tmp - tPos);
        tPos = tmp + 1;
        tmp = headersStr.find(MBGC_Params::MATCH_MARK, tIdHeadersPos[tId]);
        A_append(threadOutBuffer[threadId], headersStr, tIdHeadersPos[tId], tmp - tIdHeadersPos[tId]);
        tIdHeadersPos[tId] = tmp + 1;
    }
    A_append(threadOutBuffer[threadId], headerTemplate, tPos, headerTemplate.size() - tPos);
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::moveToFile(const string& filename, string& src, const int thread_no, bool append) {
    int threadId = lazyMode && !sequentialDecoding ? omp_get_thread_num() : 0;
    string filepath = params->outputPath + (params->ignoreFastaFilesPath ? filename.substr(filename.find_last_of("/\\") + 1) : filename);
    std::replace(filepath.begin(), filepath.end(), '\\', '/');
#ifdef DEVELOPER_BUILD
    if (params->validationMode) {
        if (filesCount == 1 && targetsCount > 1 && &src != &threadGzOut[threadId]) {
            threadGzOut[threadId].append(src);
            src.clear();
            return;
        }
        bool ok = true;
        if (!params->skipActualValidation) {
            string tmpfile = filepath;
            {
                if (!fstream(filepath)) {
                    tmpfile += ".gz";
                    if (!fstream(tmpfile)) {
                        if (params->invalidFilesCount < MBGC_Params::VALIDATION_LOG_LIMIT)
                            fprintf(stderr, "Cannot find %s for validation.\n", filepath.c_str());
                        ok = false;
                    }
                }
            }
            if (ok) {
                mbgcInFile file = mbgcInOpen(tmpfile.c_str());
                if (src.size() != file.size) {
                    if (params->invalidFilesCount < MBGC_Params::VALIDATION_LOG_LIMIT)
                        *PgHelpers::devout << "Validation ERROR: ~" << (params->invalidFilesCount + params->validFilesCount) <<
                                           ". " << filepath << " size differ (" << src.size() << " instead of "
                                           << file.size << ")" << endl;
                    ok = false;
                }
                if (ok && memcmp(src.data(), file.out, file.size) != 0) {
                    if (params->invalidFilesCount < MBGC_Params::VALIDATION_LOG_LIMIT)
                        *PgHelpers::devout << "Validation ERROR: ~" << (params->invalidFilesCount + params->validFilesCount) <<
                                           ". " << filepath << " contents differ." << endl;
                    ok = false;
                }
                if (!ok && params->dumpedFilesCount++ < MBGC_Params::VALIDATION_DUMP_LIMIT) {
                    string dumpfilepath = string(MBGC_Params::VALIDATION_DUMP_PATH) + params->inArchiveFileName + "/" +
                            filepath;
                    PgHelpers::createFolders(dumpfilepath);
                    fstream fstr(dumpfilepath, ios::out | ios::binary | (append ? ios::app : ios::trunc));
                    if (fstr) {
                        PgHelpers::writeArray(fstr, (void *) src.data(), src.size());
                        *PgHelpers::devout << "Dumped invalid file to " << dumpfilepath << endl;
                    } else
                        fprintf(stderr, "Error: cannot create a validation dump file %s\n", dumpfilepath.c_str());
                    fstr.close();
                    string tmp = string(src.begin(),
                                        std::mismatch(src.begin(), src.end(), file.out, file.out + file.size).first);
                    if (tmp.find('>') != string::npos) {
                        int seqIdx = std::count(tmp.begin(), tmp.end(), '>') - 1;
                        tmp = tmp.substr(tmp.rfind('>'));
                        if (tmp.find('\n') == string::npos)
                            *PgHelpers::devout << "Error in header:\t\tseqIdx = " << seqIdx << endl;
                        else {
                            tmp = tmp.substr(tmp.find('\n'));
                            int seqPos = tmp.size() - std::count(tmp.begin(), tmp.end(), '\n');
                            *PgHelpers::devout << "Error location:\t\tseqIdx = " << seqIdx << "\tseqPos = " << seqPos << endl;
                        }
                    } else
                        *PgHelpers::devout << "Error in FASTA format - no sequences found in equal part." << endl;
                }
                mbgcInClose(file);
            }
        }
        src.clear();
        if (ok)
#pragma omp atomic update
            params->validFilesCount++;
        else
#pragma omp atomic update
            params->invalidFilesCount++;
        return;
    }
#endif
    string& outStr = params->decompressionToGzCoderLevel || params->repackCommand ? threadGzOut[threadId] : src;
    if (filesCount == 1 && (params->decompressionToGzCoderLevel || params->repackCommand) && &outStr != &src) {
        outStr.append(src);
        src.clear();
        return;
    } else if (params->repackCommand) {
        int pos, len;
        PgHelpers::normalizePath(filepath, pos, len, params->ignoreFastaFilesPath);
        iterCurFilename = filepath.substr(pos, len);
        iterCurFileContent = std::move(src);
    } else {
        if (params->decompressionToGzCoderLevel)
            outStr = gzcompress(src, params->decompressionToGzCoderLevel);
        if (params->outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS)
            PgHelpers::writeArray(cout, (void *) outStr.data(), outStr.size());
        else {
            PgHelpers::createFolders(filepath);
            if (params->decompressionToGzCoderLevel)
                filepath += ".gz";
            fstream fstr(filepath, ios::out | ios::binary | (append ? ios::app : ios::trunc));
            if (!fstr) {
                fprintf(stderr, "Error: cannot create a file %s\n", filepath.c_str());
#ifdef DEVELOPER_BUILD
#pragma omp atomic update
                params->invalidFilesCount++;
#endif
                return;
            }
            PgHelpers::writeArray(fstr, (void *) outStr.data(), outStr.size());
            fstr.close();
        }
    }
    extractedFilesCount[thread_no]++;
    src.clear();
#ifdef DEVELOPER_BUILD
    if (filesCount > 1 || &outStr == &src)
#pragma omp atomic update
        params->validFilesCount++;
#endif
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::moveToFileIfSelected(const string& filename, string& src, const int thread_no, int targetIdx, bool append) {
#ifdef DEVELOPER_BUILD
    if (!params->validationMode)
#endif
    if (!isTargetSelected(filename, targetIdx)) {
        src.clear();
        return;
    }
    moveToFile(filename, src, thread_no, append);
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::initReference(const string &name) {
    size_t tmp;
    refG0InitPos = 0;
    refStr[refPos] = 0;
    refPos = REF_SHIFT;
    tIdLiteralPos.push_back(0);
    if (literalStr.empty())
        return;
    tmp = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, tIdLiteralPos[MASTER_THREAD_ID]);
    A_memcpy((char *) refStr.data() + refPos, literalStr.data() + tIdLiteralPos[MASTER_THREAD_ID],
             tmp - tIdLiteralPos[MASTER_THREAD_ID]);
    refPos += tmp - tIdLiteralPos[MASTER_THREAD_ID];
    refG0InitPos += tmp - tIdLiteralPos[MASTER_THREAD_ID];
    tIdLiteralPos[MASTER_THREAD_ID] = tmp + 1;
    if (!params->isRCinReferenceDisabled()) {
        char *refPtr = (char *) refStr.data() + 1;
        char *rcPtr = refPtr + refG0InitPos;
        PgHelpers::upperReverseComplement(refPtr, refG0InitPos, rcPtr);
        refPos += refG0InitPos;
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decodeReference(const string &name) {
    size_t tmp;
    refG0InitPos = 0;
    refStr[refPos] = 0;
    refPos = REF_SHIFT;
    if (isTargetSelectedPtr != nullptr) isTargetSelectedPtr++;
    size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdHTemplatesPos[MASTER_THREAD_ID]);
    string headerTemplate(headersTemplates, tIdHTemplatesPos[MASTER_THREAD_ID], hTemplateEnd - tIdHTemplatesPos[MASTER_THREAD_ID]);
    uint32_t seqCount = *seqsCountArr++;
    int64_t dnaLineLength = params->enableCustomDNAformatting ? params->dnaLineLength :
            (dnaLineLengthsArr == nullptr ? 0 : *dnaLineLengthsArr++);
    dnaLineLength = dnaLineLength ? dnaLineLength : INT64_MAX;
#ifdef DEVELOPER_BUILD
    if (params->concatHeadersAndSequencesMode) {
        threadOutBuffer[MASTER_THREAD_ID].push_back('>');
        for (uint32_t i = 0; i < seqCount; i++) {
            if (i) threadOutBuffer[MASTER_THREAD_ID].push_back('&');
            decodeHeader(headerTemplate, MASTER_THREAD_ID);
        }
        threadOutBuffer[MASTER_THREAD_ID].push_back('\n');
    }
    for (uint32_t i = 0; i < seqCount; i++) {
        if (!params->concatHeadersAndSequencesMode) {
            threadOutBuffer[MASTER_THREAD_ID].push_back('>');
            decodeHeader(headerTemplate, MASTER_THREAD_ID);
            threadOutBuffer[MASTER_THREAD_ID].push_back('\n');
        }
#else
    for (uint32_t i = 0; i < seqCount; i++) {
        threadOutBuffer[MASTER_THREAD_ID].push_back('>');
        decodeHeader(headerTemplate, MASTER_THREAD_ID);
        threadOutBuffer[MASTER_THREAD_ID].push_back('\n');
#endif
        tmp = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, tIdLiteralPos[MASTER_THREAD_ID]);
        writeDNA(literalStr.data() + tIdLiteralPos[MASTER_THREAD_ID], tmp - tIdLiteralPos[MASTER_THREAD_ID],
                 dnaLineLength, MASTER_THREAD_ID);
        A_memcpy((char *) refStr.data() + refPos, literalStr.data() + tIdLiteralPos[MASTER_THREAD_ID],
                 tmp - tIdLiteralPos[MASTER_THREAD_ID]);
        refPos += tmp - tIdLiteralPos[MASTER_THREAD_ID];
        refG0InitPos += tmp - tIdLiteralPos[MASTER_THREAD_ID];
        tIdLiteralPos[MASTER_THREAD_ID] = tmp + 1;
    }
    if (params->storeFileSeparatorMarksInHeadersStream)
        tIdHeadersPos[MASTER_THREAD_ID]++;
    tIdHTemplatesPos[MASTER_THREAD_ID] = hTemplateEnd + 1;
    moveToFileIfSelected(name, threadOutBuffer[MASTER_THREAD_ID], 0);
    if (!params->isRCinReferenceDisabled()) {
        char *refPtr = (char *) refStr.data() + 1;
        char *rcPtr = refPtr + refG0InitPos;
        PgHelpers::upperReverseComplement(refPtr, refG0InitPos, rcPtr);
        refPos += refG0InitPos;
    }
}

template<bool lazyMode> uint32_t MBGC_Decoder<lazyMode>::decodeSequenceAndReturnUnmatchedChars(string &dest, size_t refLockPos, const int tId) {
    int threadId = lazyMode && !sequentialDecoding ? omp_get_thread_num() : 0;
    size_t seqEnd = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, tIdLiteralPos[tId]);
    char* extStrEndPtr = params->enableExtensionsWithMismatches ? (char*) threadMatchExtension[threadId].data() + threadMatchExtension[threadId].size() : nullptr;
#ifdef DEVELOPER_BUILD
    threadLitRefExtPos[threadId].clear();
    threadLitRefExtLen[threadId].clear();
#endif
    dest.clear();
    uint32_t unmatchedChars = 0;
    uint32_t minMatchLength = 0;

    int64_t pairedGapOffsetDelta[MBGC_Params::MAX_GAP_DEPTH];
    memcpy(pairedGapOffsetDelta, pairedGapOffsetDeltaInit, sizeof(pairedGapOffsetDeltaInit));
    const int NO_GAP = -1;
    int64_t gapStartIdx = NO_GAP;
    int64_t gapEndIdx = NO_GAP;
    int gapCurIdx = 0;
    uint64_t matchSrcPos = 0;
    uint64_t prevMatchDestPos = 0;
    int64_t offsetDelta = -1;
    size_t extLeftLen = 0;
    size_t extRightLen = 0;
    bool isGap = false;
    int64_t j = 0;
    uint64_t markPos = literalStr.find(MBGC_Params::MATCH_MARK, tIdLiteralPos[tId]);
    while (markPos != std::string::npos && markPos < seqEnd) {
        size_t literalsLeft = markPos - tIdLiteralPos[tId];

        matchSrcPos = 0;
        bool skipOffset = pairedGapOffsetDelta[gapCurIdx] != INT64_MAX;
        if (skipOffset) {
            matchSrcPos = pairedGapOffsetDelta[gapCurIdx] + dest.size() + literalsLeft;
            pairedGapOffsetDelta[gapCurIdx] = INT64_MAX;
        } else {
            PgHelpers::readValue<uint32_t>(tIdMapOffSrc[tId], (uint32_t&) matchSrcPos);
            if (refTotalLength > UINT32_MAX) {
                uint8_t tmp = *tIdMapOff5thBytePtr[tId]++;
                matchSrcPos += (((uint64_t) tmp) << 32);
            }
        }

        extLeftLen = 0;
        if (params->enableExtensionsWithMismatches) {
            if (!isGap && literalsLeft) {
                extLeftLen = extendMatchLeft(extStrEndPtr, matchSrcPos, skipOffset, refLockPos, markPos, tId);
            }
            if (gapEndIdx == j) {
                gapStartIdx = NO_GAP;
                gapEndIdx = NO_GAP;
            }
        }
        size_t literalLen = markPos - tIdLiteralPos[tId] + extLeftLen + extRightLen;
#ifdef DEVELOPER_BUILD
        processLiteral(dest.size() - extRightLen, literalLen);
#endif
        A_append(dest, literalStr, tIdLiteralPos[tId], markPos - tIdLiteralPos[tId]);
        A_append(dest, extStrEndPtr - extLeftLen, extLeftLen);
        unmatchedChars += literalLen;
        tIdLiteralPos[tId] = markPos + 1;

        uint32_t matchLength = *tIdMapLenPtr[tId]++;
        matchLength += minMatchLength;

        prevMatchDestPos = dest.size();
        A_append(dest, refStr, matchSrcPos, matchLength);

        markPos = literalStr.find(MBGC_Params::MATCH_MARK, tIdLiteralPos[tId]);
        uint8_t gapDelta = 0;
        if (params->gapDepthOffsetEncoding && markPos < seqEnd && (gapDelta = *tIdGapDeltaPtr[tId]++)) {
            int gapIdx = gapCurIdx;
            int g = gapDelta;
            if (!params->lazyDecompressionSupport && gapStartIdx == NO_GAP && markPos - tIdLiteralPos[tId] == 0) {
                gapIdx = (gapIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;
                g++;
            }
            while (gapDelta) {
                gapIdx = (gapIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;
                if (pairedGapOffsetDelta[gapIdx] == INT64_MAX)
                    gapDelta--;
                else
                    g++;
            }
            pairedGapOffsetDelta[gapIdx] = (int64_t) matchSrcPos - (int64_t) prevMatchDestPos;
            if (params->enableExtensionsWithMismatches && gapEndIdx <= j + g && g <= params->gapDepthMismatchesEncoding) {
                gapStartIdx = j;
                gapEndIdx = j + g;
            }
        }
        gapCurIdx = (gapCurIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;

        bool gapStart = gapStartIdx == j;
        bool gapEnd = gapEndIdx == j + 1;
        bool gapMiddle = gapStartIdx < j && j + 1 < gapEndIdx;
        isGap = gapStart || gapMiddle || gapEnd;
        extRightLen = 0;
        if (params->enableExtensionsWithMismatches) {
            if (!isGap || gapStart) {
                offsetDelta = (int64_t) matchSrcPos + matchLength - dest.size();
            }
            extRightLen = extendMatchRight(dest, offsetDelta, isGap, gapStart, gapMiddle, gapEnd,
                                           markPos < seqEnd ? markPos : seqEnd, tId);
        }
        j++;
    }
    size_t literalLen = seqEnd - tIdLiteralPos[tId] + extRightLen;
#ifdef DEVELOPER_BUILD
    processLiteral(dest.size() - extRightLen, literalLen);
#endif
    A_append(dest, literalStr, tIdLiteralPos[tId], seqEnd - tIdLiteralPos[tId]);
    unmatchedChars += literalLen;
    tIdLiteralPos[tId] = seqEnd + 1;
    return unmatchedChars;
}

template<bool lazyMode> size_t MBGC_Decoder<lazyMode>::extendMatchRight(string& dest, int64_t& offsetDelta,
                                      bool isGap, bool gapStart, bool gapMiddle, bool gapEnd,
                                      size_t guardLitPos, const int tId) {
    if (tIdLiteralPos[tId] == guardLitPos && !gapMiddle)
        return 0;
    size_t destStart = dest.size();
    const char *srcPtr = refStr.data() + dest.size() + offsetDelta;
    if (gapStart || !isGap)
        dest.push_back(mismatchesCoder->code2mismatch(*srcPtr, literalStr[tIdLiteralPos[tId]++]));
    else
        srcPtr--;
    int mismatches_score = params->mmsMismatchesInitialScore;
    while ((!gapEnd || tIdLiteralPos[tId] != guardLitPos)
        && (isGap || mismatches_score < params->mmsMismatchesScoreThreshold)) {
        bool mismatch = *tIdGapMismatchesFlagsPtr[tId]++;
        if (mismatch && tIdLiteralPos[tId] == guardLitPos)
            break;
        if (mismatch) {
            mismatches_score += params->mmsMismatchPenalty;
        } else {
            mismatches_score -= params->mmsMatchBonus;
            if (mismatches_score < 0) mismatches_score = 0;
        }
        dest.push_back(mismatch ? mismatchesCoder->code2mismatch(*++srcPtr, literalStr[tIdLiteralPos[tId]++]) : *++srcPtr);
    };
    return dest.size() - destStart;
}

template<bool lazyMode> size_t MBGC_Decoder<lazyMode>::extendMatchLeft(char *extStrEndPtr, uint64_t& matchSrcPos, bool skipOffset,
                                     size_t refLockPos, size_t markPos, const int tId) {
    char *ptr = extStrEndPtr;
    const char *srcMatch = refStr.data() + matchSrcPos;
    const char *srcGuard = srcMatch - MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH;
    if (!skipOffset) {
        if (srcGuard < refStr.data() + REF_SHIFT)
            srcGuard = refStr.data() + REF_SHIFT;
        const char *srcLock = refStr.data() + refLockPos;
        if (srcGuard < srcLock && srcLock <= srcMatch)
            srcGuard = srcLock;
    }
    if (srcGuard == srcMatch)
        return 0;

    if (skipOffset) {
        const char *srcPtr = srcMatch - 1;
        size_t length = 0;
        size_t mismatches = 0;
        bool isExtLenKnown = true;
        int mismatches_score = params->mmsMismatchesInitialScore;
        while (--srcPtr > srcGuard && mismatches_score < params->mmsMismatchesScoreThreshold) {
            bool mismatch = *(tIdGapMismatchesFlagsPtr[tId] + length++);
            if (mismatch && tIdLiteralPos[tId] + ++mismatches == markPos) {
                isExtLenKnown = false;
                break;
            }
            if (mismatch) {
                mismatches_score += params->mmsMismatchPenalty;
            } else {
                mismatches_score -= params->mmsMatchBonus;
                if (mismatches_score < 0) mismatches_score = 0;
            }
        }
        if (srcPtr == srcGuard && isExtLenKnown) {
            mismatches++;
            length++;
        }
        size_t matchingChars = length - mismatches;
        matchSrcPos += matchingChars;
        srcGuard += matchingChars;
        srcMatch += matchingChars;
    }

    const char *srcPtr = srcMatch - 1;
    *(--ptr) = mismatchesCoder->code2mismatch(*srcPtr, literalStr[tIdLiteralPos[tId]++]);
    int mismatches_score = params->mmsMismatchesInitialScore;
    while (--srcPtr >= srcGuard && mismatches_score < params->mmsMismatchesScoreThreshold) {
        bool mismatch = *tIdGapMismatchesFlagsPtr[tId]++;
        if (mismatch && tIdLiteralPos[tId] == markPos)
            break;
        if (mismatch) {
            mismatches_score += params->mmsMismatchPenalty;
        } else {
            mismatches_score -= params->mmsMatchBonus;
            if (mismatches_score < 0) mismatches_score = 0;
        }
        *(--ptr) = mismatch ? mismatchesCoder->code2mismatch(*srcPtr, literalStr[tIdLiteralPos[tId]++]) : *srcPtr;
    }
    return extStrEndPtr - ptr;
}

#ifdef DEVELOPER_BUILD
template<bool lazyMode> void MBGC_Decoder<lazyMode>::processLiteral(size_t pos, size_t literalLen) {
    int threadId = lazyMode && !sequentialDecoding ? omp_get_thread_num() : 0;
    if (params->isLiteralProperForRefExtension(literalLen)) {
        size_t before = params->refLiteralBeforeAfterExt < pos ? params->refLiteralBeforeAfterExt : pos;
        threadLitRefExtPos[threadId].push_back(pos - before);
        threadLitRefExtLen[threadId].push_back(before + literalLen + params->refLiteralBeforeAfterExt);
    }
}
#endif

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decodeTarget(size_t targetIdx, int tId) {
    int threadId = lazyMode && !sequentialDecoding ? omp_get_thread_num() : 0;
    threadOutBuffer[threadId].reserve(largestFileLength);
    threadSeqStr[threadId].reserve(largestContigLength);
    threadExtStr[threadId].reserve(largestContigLength);
#ifdef DEVELOPER_BUILD
    threadLitRefExtPos[threadId].resize(matchesPerTargetEstimate);
    threadLitRefExtLen[threadId].resize(matchesPerTargetEstimate);
#endif
    if (params->enableExtensionsWithMismatches)
        threadMatchExtension[threadId].resize(MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH);
    size_t locRefPos = 0;
    bool lazyActive = lazyMode && targetIdx < lazyDecompressionTargetsGuard;
    size_t& tRefPos = lazyActive ? locRefPos : refPos;
    if (lazyActive)
        tRefPos = refExtPosArr[targetIdx];
    uint8_t unmatchedFractionFactor = unmatchedFractionFactorsArray[params->mbgcVersionMajor > 1 ? 2 * targetIdx : targetIdx];
    uint8_t unmatchedFractionRCFactor = params->mbgcVersionMajor > 1 ? unmatchedFractionFactorsArray[2 * targetIdx + 1] : unmatchedFractionFactor;
    threadSeqStr[threadId].clear();
    size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdHTemplatesPos[tId]);
    string headerTemplate(headersTemplates, tIdHTemplatesPos[tId], hTemplateEnd - tIdHTemplatesPos[tId]);
    uint32_t seqCount = seqsCountArr[targetIdx];
    int64_t dnaLineLength = params->enableCustomDNAformatting ? params->dnaLineLength :
            (dnaLineLengthsArr == nullptr ? 0 : dnaLineLengthsArr[targetIdx]);
    dnaLineLength = dnaLineLength ? dnaLineLength : INT64_MAX;
    size_t startPos = tRefPos;
    size_t refLockPos = matchingLocksPosArr != nullptr ? matchingLocksPosArr[targetIdx] : this->refTotalLength;

#ifdef DEVELOPER_BUILD
    if (params->concatHeadersAndSequencesMode) {
        threadOutBuffer[threadId].push_back('>');
        for (uint32_t i = 0; i < seqCount; i++) {
            if (i) threadOutBuffer[threadId].push_back('&');
            decodeHeader(headerTemplate, tId);
        }
        threadOutBuffer[threadId].push_back('\n');
    }
    for (uint32_t i = 0; i < seqCount; i++) {
        if (!params->concatHeadersAndSequencesMode) {
            threadOutBuffer[threadId].push_back('>');
            decodeHeader(headerTemplate, tId);
            threadOutBuffer[threadId].push_back('\n');
        }
#else
    for (uint32_t i = 0; i < seqCount; i++) {
        threadOutBuffer[threadId].push_back('>');
        decodeHeader(headerTemplate, tId);
        threadOutBuffer[threadId].push_back('\n');
#endif
        uint32_t unmatchedChars;
        unmatchedChars = decodeSequenceAndReturnUnmatchedChars(threadSeqStr[threadId], refLockPos, tId);
        writeDNA(threadSeqStr[threadId].data(), threadSeqStr[threadId].size(), dnaLineLength, threadId);
        bool loadContigToRef = params->isContigProperForRefExtension(threadSeqStr[threadId].size(), unmatchedChars, unmatchedFractionFactor);
#ifdef DEVELOPER_BUILD
        if (!loadContigToRef)
            prepareLiteralRefExtension();
#endif
        char* extPtr = (char*) (loadContigToRef ? threadSeqStr[threadId].data() : threadExtStr[threadId].data());
        size_t extSize = loadContigToRef ? threadSeqStr[threadId].size() : threadExtStr[threadId].size() ;
        loadRef(extPtr, extSize, refLockPos, tRefPos, false);
        if (!params->isRCinReferenceDisabled() && params->contigsIndivduallyReversed &&
            params->isContigProperForRefRCExtension(threadSeqStr[threadId].size(), unmatchedChars, unmatchedFractionRCFactor))
            loadRef(extPtr, extSize, refLockPos, tRefPos, true);
    }
    if (!params->isRCinReferenceDisabled() && !params->contigsIndivduallyReversed) {
        if (tRefPos >= startPos) {
            loadRef(refStr.data() + startPos, tRefPos - startPos, refLockPos, tRefPos, true);
        } else {
            size_t tmpStartPos = startPos;
            loadRef(refStr.data() + 1, tRefPos - 1, refLockPos, tRefPos,  true);
            loadRef(refStr.data() + tmpStartPos, refTotalLength - tmpStartPos, refLockPos, tRefPos, true);
        }
    }
    if (params->lazyDecompressionSupport) {
        char regionSep = MBGC_Params::REF_REGION_SEPARATOR;
        loadRef(&regionSep, 1, refLockPos, tRefPos, false);
    }
    if (params->storeFileSeparatorMarksInHeadersStream)
        tIdHeadersPos[tId]++;
    tIdHTemplatesPos[tId] = hTemplateEnd + 1;
    if (params->enableExtensionsWithMismatches)
        tIdGapMismatchesFlagsPtr[tId]++;
    tIdNamesPos[tId] = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[tId]) + 1;
    if (lazyActive) {
        finalizeLazyTargetDecoding(targetIdx, true);
    } else {
        decodedTargetsCount++;
    }
    if (tId == MASTER_THREAD_ID)
        masterTargetsStats++;
}

#ifdef DEVELOPER_BUILD
template<bool lazyMode> void MBGC_Decoder<lazyMode>::prepareLiteralRefExtension() {
    int threadId = lazyMode && !sequentialDecoding ? omp_get_thread_num() : 0;
    size_t extCount = threadLitRefExtPos[threadId].size();
    const char* src = (const char*) threadSeqStr[threadId].data();
    char* guard = (char*) threadSeqStr[threadId].data() + threadSeqStr[threadId].size();
    size_t i = extCount;
    while (i-- > 0 && src + threadLitRefExtPos[threadId][i] + threadLitRefExtLen[threadId][i] > guard)
        threadLitRefExtLen[threadId][i] = guard - (src + threadLitRefExtPos[threadId][i]);
    threadExtStr[threadId].clear();
    for(size_t i = 0; i < extCount; i++)
        A_append(threadExtStr[threadId], src + threadLitRefExtPos[threadId][i], threadLitRefExtLen[threadId][i]);
}

#endif
template<bool lazyMode> void MBGC_Decoder<lazyMode>::loadRef(const char *seqText, size_t seqLength, size_t refLockPos,
                                                                size_t& refPos1, bool loadRCRef)  {
    if (seqLength == 0)
        return;
    if (refPos1 == refTotalLength && refLockPos != refTotalLength) {
        reachedRefLengthCount++;
        refPos1 = REF_SHIFT;
    }
    size_t tmpLength = seqLength;
    size_t tmpMax = refLockPos < refPos1 ? refTotalLength : refLockPos;
    if (refPos1 + tmpLength > tmpMax) {
        tmpLength = tmpMax - refPos1;
    }
    if (loadRCRef) {
        PgHelpers::upperReverseComplement(seqText + seqLength - tmpLength, tmpLength,
                                          (char *) refStr.data() + refPos1);
    } else
        A_memcpy((char*) refStr.data() + refPos1, seqText, tmpLength);
    refPos1 += tmpLength;
    seqText += loadRCRef ? 0 : tmpLength;
    seqLength = refPos1 == refLockPos ? 0 : seqLength - tmpLength;
    if (params->lazyDecompressionSupport && refPos1 == refLockPos)
        refStr[refLockPos - 1] = MBGC_Params::REF_REGION_SEPARATOR;
    this->loadRef(seqText, seqLength, refLockPos, refPos1, loadRCRef);
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::extractFilesSequentially() {
    masterThreadTargetIdx = 0;
    int tId = MASTER_THREAD_ID;
    if (namesEndGuard >= tIdNamesPos[tId] && namesEndGuard != std::string::npos) {
        while (masterThreadTargetIdx < targetsCount && namePos <= namesEndGuard) {
            if (lazyMode && (!params->masterFilterPattern.empty() || isFilterListActive)
                && masterThreadTargetIdx < lazyDecompressionTargetsGuard) {
                int extractTargetIdx = prepareStreamsForNextFileToDecompress(tId, masterThreadTargetIdx);
                if (extractTargetIdx == NO_TARGET_TO_SCHEDULE)
                    break;
                if (extractTargetIdx < lazyDecompressionTargetsGuard) {
                    tId = extractTargetIdx;
                    masterThreadTargetIdx = extractTargetIdx;
                    if (!isTargetScheduled[masterThreadTargetIdx] && lowestContributingNotDecodedTarget < masterThreadTargetIdx) {
                        scheduleDependencyChain(masterThreadTargetIdx);
                        for (int i = 0; i < masterThreadTargetIdx; i++)
                            if (isTargetScheduled[i] && !isTargetProcessed[i]) {
#ifdef DEVELOPER_BUILD
                                if (params->validationMode) {
                                    namePos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[i]);
                                    currentName = namesStr.substr(tIdNamesPos[i], namePos - tIdNamesPos[i]);
                                }
#endif
                                decodeTarget(i, i);
#ifdef DEVELOPER_BUILD
                                if (params->validationMode)
                                    moveToFile(currentName, threadOutBuffer[MASTER_THREAD_ID], 0, false);
#endif
                                threadOutBuffer[MASTER_THREAD_ID].clear();
                            }
                    }
                } else {
                    for (int i = 0; i < lazyDecompressionTargetsGuard; i++)
                        if (!isTargetProcessed[i] &&
                            (i == lazyDecompressionTargetsGuard - 1 || refExtLoadedPosArr[i + 1] - refExtLoadedPosArr[i] > 1)) {
#ifdef DEVELOPER_BUILD
                            if (params->validationMode) {
                                namePos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[i]);
                                currentName = namesStr.substr(tIdNamesPos[i], namePos - tIdNamesPos[i]);
                            }
#endif
                            decodeTarget(i, i);
#ifdef DEVELOPER_BUILD
                            if (params->validationMode)
                                moveToFile(currentName, threadOutBuffer[MASTER_THREAD_ID], 0, false);
#endif
                            threadOutBuffer[MASTER_THREAD_ID].clear();
                        }
                    copyStreamsPositions(lazyDecompressionTargetsGuard - 1, MASTER_THREAD_ID);
                    refPos = refExtPosArr[lazyDecompressionTargetsGuard];
                    masterThreadTargetIdx = lazyDecompressionTargetsGuard;
                    tId = MASTER_THREAD_ID;
                }
            }
            extractNextSequentially(tId);
        }
    }
}

template<bool lazyMode>
void MBGC_Decoder<lazyMode>::extractNextSequentially(int tId) {
    if (filesCount == 1) tIdNamesPos[tId] = 0;
    namePos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[tId]);
    currentName = namesStr.substr(tIdNamesPos[tId], namePos - tIdNamesPos[tId]);
    if (lazyMode && !isTargetSelected(currentName, masterThreadTargetIdx) && masterThreadTargetIdx < lazyDecompressionTargetsGuard &&
        refExtLoadedPosArr[masterThreadTargetIdx + 1] - refExtLoadedPosArr[masterThreadTargetIdx] <= 1) {
        refStr[refExtPosArr[masterThreadTargetIdx]] = MBGC_Params::REF_REGION_SEPARATOR;
        fastForwardTargetStreams(masterThreadTargetIdx, tId);
        finalizeLazyTargetDecoding(masterThreadTargetIdx++, false);
    } else {
        decodeTarget(masterThreadTargetIdx, tId);
        moveToFileIfSelected(currentName, threadOutBuffer[MASTER_THREAD_ID], 0, masterThreadTargetIdx,
                             filesCount == 1 && (!params->g0IsTarget || masterThreadTargetIdx > 0));
        masterThreadTargetIdx++;
    }
}

template<bool lazyMode> bool MBGC_Decoder<lazyMode>::tryClaimingTarget(int tIdx) {
    bool targetClaimed = false;
#pragma omp critical(claimTargetBlock)
    {
        if (!isTargetClaimed[tIdx]) {
            isTargetClaimed[tIdx] = true;
            targetClaimed = true;
        }
    }

    return targetClaimed;
}

template<bool lazyMode> int MBGC_Decoder<lazyMode>::tryClaimingTarget() {
    int tIdx = NO_TARGET_TO_SCHEDULE;
#pragma omp critical(claimTargetBlock)
    {
        while (lowestNotClaimedParallelIdx < targetsForParallelTaskThreadDecoding.size()
               && isTargetClaimed[targetsForParallelTaskThreadDecoding[lowestNotClaimedParallelIdx]])
            lowestNotClaimedParallelIdx++;
        if (lowestNotClaimedParallelIdx < targetsForParallelTaskThreadDecoding.size()) {
            tIdx = targetsForParallelTaskThreadDecoding[lowestNotClaimedParallelIdx];
            isTargetClaimed[tIdx] = true;
            lowestNotClaimedParallelIdx++;
        }
    }

    return tIdx;
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::workerParallelTask(const int thread_no) {
    int threadId = lazyMode ? omp_get_thread_num() : 0;
    int tIdx = NO_TARGET_TO_SCHEDULE;
    while(isDecoding || in[thread_no] != out[thread_no] ){
        if (in[thread_no] != out[thread_no]) {
            string &content = contentsBuf[thread_no][out[thread_no] % writingBufferSize];
            string &name = namesBuf[thread_no][out[thread_no] % writingBufferSize];
            moveToFile(name, content, thread_no);

            out[thread_no] += 1;
        } else {
            if (lazyMode) {
                if (tIdx == NO_TARGET_TO_SCHEDULE && lowestNotClaimedParallelIdx < targetsForParallelTaskThreadDecoding.size()) {
                    tIdx = tryClaimingTarget();
                }
                if (tIdx != NO_TARGET_TO_SCHEDULE
                    && (tIdx == lowestTargetForProcessing || highestTargetMapping[tIdx] < lowestTargetForProcessing)) {
                    if (filesCount == 1) tIdNamesPos[tIdx] = 0;
                    size_t tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[tIdx]);
                    string name = namesStr.substr(tIdNamesPos[tIdx], tmp - tIdNamesPos[tIdx]);
                    decodeTarget(tIdx, tIdx);
                    moveToFileIfSelected(name, threadOutBuffer[threadId], thread_no, tIdx, filesCount == 1);
                    tIdx = NO_TARGET_TO_SCHEDULE;
                } else
                    nanosleep(WRITE_THREAD_SLEEP_TIME, nullptr);
            } else
                nanosleep(WRITE_THREAD_SLEEP_TIME, nullptr);
        }
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::extractFilesParallel() {
#pragma omp parallel default(none)
    {
#pragma omp single
        {
            size_t tmp;
            isDecoding = true;
            masterThreadTargetIdx = 0;

            if (omp_get_num_threads() - 1 < workerThreadsCount)
                workerThreadsCount = omp_get_num_threads() - 1;
            extractedFilesCount.resize(workerThreadsCount, 0);
            in.resize(workerThreadsCount, 0);
            out.resize(workerThreadsCount, 0);
            namesBuf.resize(workerThreadsCount, vector<string>(writingBufferSize));
            contentsBuf.resize(workerThreadsCount, vector<string>(writingBufferSize));
            for (int i = 0; i < workerThreadsCount; i++) {
#pragma omp task
                {
                    if (lazyMode && i == 0)
                        scheduleParallelLazyDecompression();
                    workerParallelTask(i);
                }
            }
            int io_thread_no = workerThreadsCount > 1 ? 1 : 0;
            int threadId = lazyMode ? omp_get_thread_num() : 0;
            if (namesEndGuard >= tIdNamesPos[MASTER_THREAD_ID] && namesEndGuard != std::string::npos) {
                while (masterThreadTargetIdx < targetsCount) {
                    int ii = workerThreadsCount;
                    while (in[io_thread_no] == out[io_thread_no] + writingBufferSize) {
                        io_thread_no = (io_thread_no + 1) % workerThreadsCount;
                        if (--ii)
                            continue;
                        ii = workerThreadsCount;
                        nanosleep(SLEEP_TIME, nullptr);
                    }
                    bool targetClaimed = false;
                    if (lazyMode && masterThreadTargetIdx < scheduledTargetsGuard) {
                        while (masterThreadTargetIdx < scheduledTargetsGuard) {
                            if (isTargetScheduled[masterThreadTargetIdx] && !isTargetClaimed[masterThreadTargetIdx]) {
                                if (tryClaimingTarget(masterThreadTargetIdx)) {
                                    targetClaimed = true;
                                    break;
                                }
                            }
                            masterThreadTargetIdx++;
                        }
                        if (masterThreadTargetIdx >= targetsCount)
                            break;
                        if (masterThreadTargetIdx < scheduledTargetsGuard) {
                            copyStreamsPositions(masterThreadTargetIdx, MASTER_THREAD_ID);
                            while (masterThreadTargetIdx != lowestTargetForProcessing && highestTargetMapping[masterThreadTargetIdx] >= lowestTargetForProcessing)
                                nanosleep(SLEEP_TIME, nullptr);
                        }
                    }
                    if (!lazyMode || targetClaimed || tryClaimingTarget(masterThreadTargetIdx)) {
                        if (lazyMode && scheduledTargetsGuard && masterThreadTargetIdx == scheduledTargetsGuard) {
                            while (lowestTargetForProcessing < scheduledTargetsGuard)
                                nanosleep(SLEEP_TIME, nullptr);
                            bool masterDecodedFurthestTarget = tIdLiteralPos[MASTER_THREAD_ID] > tIdLiteralPos[scheduledTargetsGuard - 1];
                            if (!masterDecodedFurthestTarget)
                                copyStreamsPositions(scheduledTargetsGuard - 1, MASTER_THREAD_ID);
                            refPos = refExtPosArr[scheduledTargetsGuard];
                        }
                        tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[MASTER_THREAD_ID]);
                        const string &filename = namesStr.substr(tIdNamesPos[MASTER_THREAD_ID], tmp - tIdNamesPos[MASTER_THREAD_ID]);
                        bool extractFile = isTargetSelected(filename, masterThreadTargetIdx);
                        if (!extractFile && lazyMode && masterThreadTargetIdx < scheduledTargetsGuard &&
                            refExtLoadedPosArr[masterThreadTargetIdx + 1] - refExtLoadedPosArr[masterThreadTargetIdx] <= 1) {
                            refStr[refExtPosArr[masterThreadTargetIdx]] = MBGC_Params::REF_REGION_SEPARATOR;
                            fastForwardTargetStreams(masterThreadTargetIdx, MASTER_THREAD_ID);
                            finalizeLazyTargetDecoding(masterThreadTargetIdx++, false);
                        } else {
#ifdef DEVELOPER_BUILD
                            extractFile |= params->validationMode;
#endif
                            if (extractFile)
                                namesBuf[io_thread_no][in[io_thread_no] % writingBufferSize] = filename;
                            threadOutBuffer[threadId] = std::move(
                                    contentsBuf[io_thread_no][in[io_thread_no] % writingBufferSize]);
                            threadOutBuffer[threadId].reserve(largestFileLength);

                            decodeTarget(masterThreadTargetIdx++, MASTER_THREAD_ID);
                            if (extractFile) {
                                contentsBuf[io_thread_no][in[io_thread_no] % writingBufferSize] = std::move(
                                        threadOutBuffer[threadId]);
                                in[io_thread_no] += 1;
                                io_thread_no = (io_thread_no + 1) % workerThreadsCount;
                            }
                        }
                    }
                    tIdNamesPos[MASTER_THREAD_ID] = tmp + 1;
                    if (tmp > namesEndGuard)
                        break;
                }
            }
            while (lowestTargetForProcessing < scheduledTargetsGuard)
                nanosleep(SLEEP_TIME, nullptr);
            isDecoding = false;
        }
    }
    if (lazyMode && params->lazyDecompressionSupport && targetsCount == scheduledTargetsGuard) {
        refPos = refExtPosArr[targetsCount];
        reachedRefLengthCount = (refExtLoadedPosArr[targetsCount] - refPos) / (refTotalLength - 1);
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::finalizeLazyTargetDecoding(int targetIdx, bool targetActuallyDecoded) {
#pragma omp critical(finalizeTargetDecodingBlock)
    {
        isTargetProcessed[targetIdx] = true;
        if (targetActuallyDecoded)
            decodedTargetsCount++;
        while (lowestContributingNotDecodedTarget < lazyDecompressionTargetsGuard
               && (isTargetProcessed[lowestContributingNotDecodedTarget] ||
                   refExtLoadedPosArr[lowestContributingNotDecodedTarget + 1] - refExtLoadedPosArr[lowestContributingNotDecodedTarget] <= 1))
            lowestContributingNotDecodedTarget++;
        refPos = refExtPosArr[lowestContributingNotDecodedTarget];
        while (lowestTargetForProcessing < targetsCount &&
                (isTargetProcessed[lowestTargetForProcessing] ||
                (lowestTargetForProcessing < scheduledTargetsGuard && !isTargetScheduled[lowestTargetForProcessing])))
            lowestTargetForProcessing++;
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decodeRefExtSizeStream(string &refExtSizeStream) {
    lazyDecompressionTargetsGuard = 0;
    refExtLoadedPosArr.push_back(refPos);
    refExtPosArr.push_back(refPos);
    WrapperStrStream refExtSizeSrc(refExtSizeStream);
    refExtLoadedPosArr.reserve(refExtSizeStream.size() / 3);
    refExtPosArr.reserve(refExtSizeStream.size() / 3);
    uint64_t refExtSize = 0;
    size_t diff = 0;
    while (!refExtSizeSrc.isEnd()) {
        PgHelpers::readUInt64Frugal<uint64_t>(refExtSizeSrc, refExtSize);
        uint64_t refExtPos = refExtLoadedPosArr.back() + refExtSize;
        if (refExtPos <= refTotalLength)
            lazyDecompressionTargetsGuard++;
        refExtLoadedPosArr.push_back(refExtPos);
        refExtPos -= diff;
        if (refExtPos > refTotalLength) {
            diff += refTotalLength - REF_SHIFT;
            refExtPos -= refTotalLength - REF_SHIFT;
        }
        refExtPosArr.push_back(refExtPos);
        if (refExtSize > 1)
            refContributingFilesCount++;
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decodeMapLenStream(string &mapLenStream, vector<uint32_t> &mapLenArr) {
    if (!params->frugal64bitLenEncoding && params->isVersionAtLeast(2, 0)) {
        tIdMapLenPtr.push_back((uint32_t*) mapLenStream.data());
        return;
    }
    mapLenArr.reserve(mapLenStream.size() / 2);
    WrapperStrStream mapLenSrc(mapLenStream);
    uint32_t matchLength = 0;
    if (params->frugal64bitLenEncoding) {
        while (!mapLenSrc.isEnd()) {
            PgHelpers::readUInt64Frugal<uint32_t>(mapLenSrc, matchLength);
            mapLenArr.push_back(matchLength);
        }
    } else {
        while (!mapLenSrc.isEnd()) {
            PgHelpers::readUIntWordFrugal(mapLenSrc, matchLength);
            mapLenArr.push_back(matchLength);
        }
    }
    tIdMapLenPtr.push_back(mapLenArr.data());
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::processFilterPatterns() {
    namesEndGuard = namesStr.rfind(params->masterFilterPattern);
    if (namesEndGuard == std::string::npos || params->inputFileName.empty())
        return;
    isFilterListActive = true;
    isFileSelectedFlag.resize(filesCount, 0);
    isTargetSelectedPtr = isFileSelectedFlag.data();
    ifstream listSrc(params->inputFileName, ios_base::in | ios_base::binary);
    if (listSrc.fail()) {
        fprintf(stderr, "cannot open file patterns list file %s\n", params->inputFileName.c_str());
        exit(EXIT_FAILURE);
    }
    vector<string> patterns;
    string line;
    while (getline(listSrc, line)) {
        if(!line.empty() && *line.rbegin() == '\r') line.resize(line.length() - 1);
        patterns.emplace_back(std::move(line));
    }
    if (patterns.empty()) {
        fprintf(stderr, "ERROR: file patterns list file is empty.\n");
        exit(EXIT_FAILURE);
    }

    int fileIdx = -1;
    int tmp, namesPos = -1;
    size_t tmpGuard = 0;
    do {
        size_t masterMatchPos = namesStr.find(params->masterFilterPattern, namesPos + 1);
        do {
            tmp = namesPos + 1;
            fileIdx++;
            namesPos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tmp);
        } while (namesPos < masterMatchPos);
        string name = namesStr.substr(tmp, namesPos - tmp);
        for(string& pattern : patterns)
            if (name.find(pattern) != std::string::npos) {
                isFileSelectedFlag[fileIdx] = 1;
                tmpGuard = namesPos;
                break;
            }
    } while (namesPos < namesEndGuard);
    namesEndGuard = tmpGuard;
}

template<bool lazyMode> bool MBGC_Decoder<lazyMode>::isTargetSelected(const string &filename, int targetIdx) {
    if (params->appendCommand)
        return false;
    if (!isFilterListActive)
        return filename.find(params->masterFilterPattern) != string::npos;
    if (targetIdx == FIRST_FILE_TARGET_IDX)
        return isFileSelectedFlag[0];
    return isTargetSelectedPtr[targetIdx];
}

template<bool lazyMode> bool MBGC_Decoder<lazyMode>::isFileSelected(const string &filename, int fileIdx) {
    return isTargetSelected(filename, fileIdx - (params->g0IsTarget ? 0 : 1));
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::getSelectedFilesNames(vector<string> &filenames) {
    size_t tmpNamesPos;
    size_t namesPos = 0;
    for (int i = 0; i < filesCount; i++) {
        tmpNamesPos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
        string name = namesStr.substr(namesPos, tmpNamesPos - namesPos);
        if (isFileSelected(name, i)) {
            name = params->ignoreFastaFilesPath ? name.substr(name.find_last_of("/\\") + 1) : name;
            filenames.push_back(name);
        }
        namesPos = tmpNamesPos + 1;
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decodeInit() {
    readStats(*inStream);
#ifdef __MINGW32__
    if (params->isStdoutMode(true)) {
        if (_setmode(_fileno(stdout), _O_BINARY) == -1)
            fprintf(stderr, "WARNING: switching cout to binary mode failed (errCode: %d)\n", strerror(errno));
    }
#endif
    if (params->appendCommand && filesCount == 1 && params->inputFileName.empty()) {
        fprintf(stderr, "ERROR: Appending single FASTA archive requires specifying input file option (%c).\n", MBGC_Params::INPUT_FILE_OPT);
        exit(EXIT_FAILURE);
    }
    if (params->appendCommand && filesCount > 1 && !params->inputFileName.empty()) {
        fprintf(stderr, "ERROR: Appending multiple FASTA archive cannot be used with input file option (%c).\n", MBGC_Params::INPUT_FILE_OPT);
        exit(EXIT_FAILURE);
    }
    threadSeqStr.resize(1);
    threadExtStr.resize(1);
    threadGzOut.resize(1);
    refStr.resize(refTotalLength + params->maxConsecutiveMismatches);
    threadOutBuffer.resize(1);
    destStrings.push_back(&namesStr);
    destStrings.push_back(&seqsCount);
    destStrings.push_back(&headersTemplates);
    destStrings.push_back(&headersStr);
    if (params->mbgcVersionMajor > 1)
        destStrings.push_back(&dnaLineLengthsStr);
    destStrings.push_back(&unmatchedFractionFactorsStr);
    destStrings.push_back(&literalStr);
    string rcMapOff, rcMapLen;
    if (params->rcRedundancyRemoval) {
        destStrings.push_back(&rcMapOff);
        destStrings.push_back(&rcMapLen);
    }
    if (params->isVersionAtLeast(1, 1))
        destStrings.push_back(&matchingLocksPosStream);
    if (params->isVersionAtLeast(1, 3))
        destStrings.push_back(&gapDeltaStr);
    if (params->enableExtensionsWithMismatches) {
        destStrings.push_back(&gapMismatchesFlagsStr);
        threadMatchExtension.emplace_back("");
        threadMatchExtension[0].resize(MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH);
    }
    destStrings.push_back(&mapOffStream);
    if (refTotalLength > UINT32_MAX)
        destStrings.push_back(&mapOff5thByteStream);
    destStrings.push_back(&mapLenStream);
    if (lazyMode)
        destStrings.push_back(&refExtSizeStream);
    params->enableOmpThreads(params->noOfThreadsLimited ? params->coderThreads : omp_get_num_procs());
    readCompressedCollectiveParallel(*inStream, destStrings);
    unmatchedFractionFactorsArray = (uint8_t *) unmatchedFractionFactorsStr.data();
    targetsCount = unmatchedFractionFactorsStr.size() / (params->mbgcVersionMajor > 1 ? 2 : 1);
#ifdef DEVELOPER_BUILD
    matchesPerTargetEstimate = (mapOffStream.size() / sizeof(uint32_t)) / targetsCount;
    threadLitRefExtPos.resize(1);
    threadLitRefExtLen.resize(1);
#endif
    seqsCountArr = (uint32_t *) seqsCount.data();
    tIdMapOffSrc.emplace_back(mapOffStream);
    tIdMapOff5thBytePtr.emplace_back(refTotalLength > UINT32_MAX ? (uint8_t *) mapOff5thByteStream.data() : nullptr);
    tIdGapDeltaPtr.push_back(params->isVersionAtLeast(1, 3) ?
                             (uint8_t *) gapDeltaStr.data() : nullptr);
    gapDeltaGuard = (uint8_t *) gapDeltaStr.data() + gapDeltaStr.size();
    tIdGapMismatchesFlagsPtr.push_back(params->enableExtensionsWithMismatches ? (uint8_t *) gapMismatchesFlagsStr.data() : nullptr);
    gapMismatchesFlagsGuard = (uint8_t *) gapMismatchesFlagsStr.data() + gapMismatchesFlagsStr.size();

    decodeMapLenStream(mapLenStream, mapLenArr);

    dnaLineLengthsArr = dnaLineLengthsStr.size() ? (size_t *) dnaLineLengthsStr.data() : nullptr;
    matchingLocksPosArr = matchingLocksPosStream.size() ? (size_t *) matchingLocksPosStream.data() : nullptr;

    if (params->rcRedundancyRemoval)
        PgTools::SimpleSequenceMatcher::restoreRCMatchedSequence(literalStr, rcMapOff, rcMapLen, literalsLength);
    else
        literalsLength = literalStr.size();
    *PgHelpers::appout << "loaded streams - " << PgHelpers::time_millis() << " [ms]" << endl;

    processFilterPatterns();

    size_t tmp;
    tIdNamesPos.push_back(0);
    tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[MASTER_THREAD_ID]);
    g0name = namesStr.substr(tIdNamesPos[MASTER_THREAD_ID], tmp - tIdNamesPos[MASTER_THREAD_ID]);
    extractedFilesCount.push_back(0);
    tIdHeadersPos.push_back(0);
    tIdHTemplatesPos.push_back(0);
    tIdLiteralPos.push_back(0);
    if (params->g0IsTarget)
        initReference(g0name);
    else {
        tIdNamesPos[MASTER_THREAD_ID] = tmp + 1;
        decodeReference(g0name);
    }
    if (filesCount == 1 || params->outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS || PgHelpers::numberOfThreads == 1)
        sequentialDecoding = true;

    if (lazyMode) {
        if (!sequentialDecoding || !params->masterFilterPattern.empty() || params->appendCommand) {
            decodeRefExtSizeStream(refExtSizeStream);
            initLazyDecompression();
        }
    }

    for (int i = 0; i < MBGC_Params::MAX_GAP_DEPTH; i++)
        pairedGapOffsetDeltaInit[i] = INT64_MAX;
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::decode() {
    if (params->infoCommand) {
        list();
        return;
    }
    decodeInit();
    if (sequentialDecoding)
        extractFilesSequentially();
    else
        extractFilesParallel();

#ifdef DEVELOPER_BUILD
    if (filesCount == 1 && (params->decompressionToGzCoderLevel ||
    (params->validationMode && params->invalidFilesCount + params->validFilesCount == 0)))
#else
        if (filesCount == 1 && params->decompressionToGzCoderLevel)
#endif
        moveToFileIfSelected(g0name, threadGzOut[MASTER_THREAD_ID], MASTER_THREAD_ID);

#ifdef DEVELOPER_BUILD
    if (params->validationMode && !params->skipActualValidation) {
        *PgHelpers::devout << "Validation" << (params->invalidFilesCount == 0 ? "" : " ERROR") <<
                           ": correctly decoded " << params->validFilesCount << " out of " << filesCount << " files."
                           << endl;
        if (params->invalidFilesCount)
            fprintf(stderr, "Validation ERROR: errors in contents of %d decoded files.\n", params->invalidFilesCount);
    }
#endif

    if (lazyMode) {
        *PgHelpers::devout << "targets decoded by - master: " << (PgHelpers::numberOfThreads > 1 ? masterTargetsStats : decodedTargetsCount);
        if (PgHelpers::numberOfThreads > 1)
            *PgHelpers::devout << " - tasks: " << decodedTargetsCount - masterTargetsStats << " (check: " << (decodedTargetsCount)
                << "; seq: " << this->targetsDecodedBeforeScheduling << "/" << this->lowestContributingNotDecodedTargetBeforeScheduling
                << "; par: " << targetsForParallelTaskThreadDecoding.size() << ")";
        *PgHelpers::devout << endl;
    }
    refStr.resize(refTotalLength);
    *PgHelpers::devout << "reference size: " << refStr.size() << endl;
    *PgHelpers::devout << "reference loading length: " << getLoadedRefLength() << endl;
    uint32_t totalExtractedFilesCount = accumulate(extractedFilesCount.begin(), extractedFilesCount.end(), 0);
    if (filesCount == 1 && totalExtractedFilesCount)
        totalExtractedFilesCount = 1;
#ifdef DEVELOPER_BUILD
    if (!params->validationMode)
#endif
    if (!params->appendCommand)
        *PgHelpers::appout << "extracted " << totalExtractedFilesCount << (totalExtractedFilesCount == 1 ? " file" : " files") << endl;
    *PgHelpers::logout << PgHelpers::time_millis() << "       \t";
    *PgHelpers::logout << "     \t              \t           \t";
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::list() {
    readStats(*inStream);
    vector<string*> destStrings;
    destStrings.push_back(&namesStr);
    string seqsCount;
    if (params->listHeadersMode) {
        destStrings.push_back(&seqsCount);
        destStrings.push_back(&headersTemplates);
        destStrings.push_back(&headersStr);
    }
    params->enableOmpThreads(params->noOfThreadsLimited ? params->coderThreads : omp_get_num_procs());
    readCompressedCollectiveParallel(*inStream, destStrings);
    seqsCountArr = (uint32_t*) seqsCount.data();
    *PgHelpers::devout << "loaded archive - " << PgHelpers::time_millis() << " [ms]" << endl;
    processFilterPatterns();
    size_t tmpNamesPos, tmpHTemplatesPos, tmpHeadersPos;
    size_t namesPos = 0;
    size_t headersPos = 0;
    string headerTemplate;
    size_t hTemplatesPos = 0;
    int hMatchMarksCount = 0;
    for (int i = 0; i < filesCount; i++) {
        tmpNamesPos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
        string name = namesStr.substr(namesPos, tmpNamesPos - namesPos);
        if (params->listHeadersMode) {
            hMatchMarksCount = 0;
            size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, hTemplatesPos);
            headerTemplate = string(headersTemplates, hTemplatesPos, hTemplateEnd - hTemplatesPos);
            while ((hTemplatesPos = headersTemplates.find(MBGC_Params::MATCH_MARK, hTemplatesPos)) != std::string::npos
                   && hTemplatesPos++ < hTemplateEnd)
                hMatchMarksCount++;
            hTemplatesPos = hTemplateEnd + 1;
        }
        if (isFileSelected(name, i)) {
            name = params->ignoreFastaFilesPath ? name.substr(name.find_last_of("/\\") + 1) : name;
            if (params->listHeadersMode) {
                for (uint32_t s = 0; s < seqsCountArr[i]; s++) {
                    *PgHelpers::appout << ">";
                    size_t tmp, tPos = 0;
                    while ((tmp = headerTemplate.find(MBGC_Params::MATCH_MARK, tPos)) != std::string::npos) {
                        *PgHelpers::appout << headerTemplate.substr(tPos, tmp - tPos);
                        tPos = tmp + 1;
                        tmp = headersStr.find(MBGC_Params::MATCH_MARK, headersPos);
                        *PgHelpers::appout << headersStr.substr(headersPos, tmp - headersPos);
                        headersPos = tmp + 1;
                    }
                    *PgHelpers::appout << headerTemplate.substr(tPos, headerTemplate.size() - tPos);
                    *PgHelpers::appout << ">" << name << endl;
                }
            } else
                *PgHelpers::appout << name << endl;
        } else if (params->listHeadersMode) {
            for (int s = 0; s < seqsCountArr[i]; s++)
                for (int m = 0; m < hMatchMarksCount; m++)
                    headersPos = headersStr.find(MBGC_Params::MATCH_MARK, headersPos) + 1;
        }
        if (params->listHeadersMode && params->storeFileSeparatorMarksInHeadersStream)
            headersPos++;
        namesPos = tmpNamesPos + 1;
    }
    *PgHelpers::logout << PgHelpers::time_millis() << "       \t";
    *PgHelpers::logout << "     \t              \t           \t";
}


template<bool lazyMode> void MBGC_Decoder<lazyMode>::readStats(istream &in) {
    PgHelpers::readValue<uint32_t>(in, filesCount);
    PgHelpers::readValue<uint64_t>(in, totalFilesLength);
    if (params->isVersionAtLeast(2, 0)) {
        PgHelpers::readValue<uint64_t>(in, refG0InitPos);
        PgHelpers::readValue<uint64_t>(in, largestFileLength);
    }else {
        int tmpInt;
        PgHelpers::readValue<int>(in, tmpInt);
        refG0InitPos = tmpInt;
        uint32_t tmpUint32;
        PgHelpers::readValue<uint32_t>(in, tmpUint32);
        largestFileLength = tmpUint32;
    }
    PgHelpers::readValue<uint32_t>(in, largestContigLength);
    if (params->isVersionAtLeast(1, 2))
        PgHelpers::readValue<uint64_t>(in, refTotalLength);
    else {
        uint32_t tmp;
        PgHelpers::readValue<uint32_t>(in, tmp);
        refTotalLength = tmp;
    }
    if (params->isVersionAtLeast(2, 0)) {
        PgHelpers::readValue<uint8_t>(in, refBuffersCount);
        if (refBuffersCount > 1) {
            fprintf(stderr, "ERROR: Unsupported using multiple reference buffers (%d).\n", refBuffersCount);
            exit(EXIT_FAILURE);
        }
        PgHelpers::readValue<uint64_t>(in, literalsLength);
    }
    if (largestFileLength > params->AVG_REF_INIT_SIZE)
        writingBufferSize = 1 + writingBufferSize * params->AVG_REF_INIT_SIZE / largestFileLength;
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::copyStreamsPositions(const int srcId, const int destId) {
    tIdHeadersPos[destId] = tIdHeadersPos[srcId];
    tIdHTemplatesPos[destId] = tIdHTemplatesPos[srcId];
    tIdNamesPos[destId] = tIdNamesPos[srcId];
    tIdLiteralPos[destId] = tIdLiteralPos[srcId];
    tIdMapOffSrc[destId] = tIdMapOffSrc[srcId];
    if (tIdMapOff5thBytePtr[0] != nullptr)
        tIdMapOff5thBytePtr[destId] = tIdMapOff5thBytePtr[srcId];
    tIdMapLenPtr[destId] = tIdMapLenPtr[srcId];
    if (params->gapDepthOffsetEncoding)
        tIdGapDeltaPtr[destId] = tIdGapDeltaPtr[srcId];
    if (params->enableExtensionsWithMismatches)
        tIdGapMismatchesFlagsPtr[destId] = tIdGapMismatchesFlagsPtr[srcId];
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::initLazyDecompression() {
    isTargetProcessed.resize(targetsCount, false);
    isTargetScheduled.resize(targetsCount, false);
    isTargetClaimed.resize(targetsCount, false);
    highestTargetMapping.resize(targetsCount, -1);
    tIdHeadersPos.resize(targetsCount);
    tIdHTemplatesPos.resize(targetsCount);
    tIdNamesPos.resize(targetsCount);
    tIdLiteralPos.resize(targetsCount);
    tIdMapOffSrc.resize(targetsCount);
    if (tIdMapOff5thBytePtr[0] != nullptr)
        tIdMapOff5thBytePtr.resize(targetsCount);
    tIdMapLenPtr.resize(targetsCount);
    if (params->gapDepthOffsetEncoding)
        tIdGapDeltaPtr.resize(targetsCount);
    int threadsCount = PgHelpers::numberOfThreads;
    if (params->enableExtensionsWithMismatches) {
        tIdGapMismatchesFlagsPtr.resize(targetsCount);
        threadMatchExtension.resize(threadsCount);
    }
    threadOutBuffer.resize(threadsCount);
    threadExtStr.resize(threadsCount);
    threadSeqStr.resize(threadsCount);
    if (params->decompressionToGzCoderLevel)
        threadGzOut.resize(threadsCount);
#ifdef DEVELOPER_BUILD
    threadLitRefExtPos.resize(threadsCount);
    threadLitRefExtLen.resize(threadsCount);
#endif
    copyStreamsPositions(MASTER_THREAD_ID, 1);
}

template<bool lazyMode> int MBGC_Decoder<lazyMode>::findStreamsPositions(int prevTargetSrcId, int startTargetIdx, int endTargetIdx, size_t namesEndGuard) {
    copyStreamsPositions(prevTargetSrcId, startTargetIdx);
    int i = 0;
    for(i = startTargetIdx; i <= endTargetIdx && namesEndGuard >= tIdNamesPos[i]; i++) {
        fastForwardTargetStreams(i - 1, i);
        if (i + 1 < targetsCount)
            copyStreamsPositions(i, i + 1);
    }
    return i;
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::fastForwardTargetStreams(int targetIdx, int tId) {
    int hMatchMarksCount = 0;
    size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdHTemplatesPos[tId]);
    while ((tIdHTemplatesPos[tId] = headersTemplates.find(MBGC_Params::MATCH_MARK, tIdHTemplatesPos[tId])) != std::string::npos
           && tIdHTemplatesPos[tId]++ < hTemplateEnd)
        hMatchMarksCount++;
    tIdHTemplatesPos[tId] = hTemplateEnd + 1;
    size_t tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[tId]);
    tIdNamesPos[tId] = tmp + 1;
    if (params->storeFileSeparatorMarksInHeadersStream) {
        tIdHeadersPos[tId] = headersStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdHeadersPos[tId]) + 1;
    } else {
        for (int s = 0; s < seqsCountArr[targetIdx]; s++)
            for (int m = 0; m < hMatchMarksCount; m++)
                tIdHeadersPos[tId] = headersStr.find(MBGC_Params::MATCH_MARK, tIdHeadersPos[tId]) + 1;
    }
    uint8_t separatorMark = MBGC_Params::FILE_SEPARATOR_MARK;
    for (int s = 0; s < seqsCountArr[targetIdx]; s++) {
        size_t seqSeparatorPos = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, tIdLiteralPos[tId]);
        size_t markPos = tIdLiteralPos[tId];
        int matchesCount = 0;
        size_t skippedOffsetsCount = 0;
        while ((markPos = literalStr.find(MBGC_Params::MATCH_MARK, markPos)) != std::string::npos &&
               markPos++ < seqSeparatorPos) {
            matchesCount++;
            if (params->gapDepthOffsetEncoding)
                if (*tIdGapDeltaPtr[tId]++) skippedOffsetsCount++;
        }
        if (matchesCount && params->gapDepthOffsetEncoding && *--tIdGapDeltaPtr[tId])
            skippedOffsetsCount--;
        tIdMapOffSrc[tId].skipOffset((matchesCount - skippedOffsetsCount) * sizeof(uint32_t));
        if (tIdMapOff5thBytePtr[0] != nullptr)
            tIdMapOff5thBytePtr[tId] += matchesCount - skippedOffsetsCount;
        tIdMapLenPtr[tId] += matchesCount;
        tIdLiteralPos[tId] = seqSeparatorPos + 1;
    }
    if (params->enableExtensionsWithMismatches)
        tIdGapMismatchesFlagsPtr[tId] = find(tIdGapMismatchesFlagsPtr[tId], gapMismatchesFlagsGuard, separatorMark) + 1;
}

template<bool lazyMode> int MBGC_Decoder<lazyMode>::prepareStreamsForNextFileToDecompress(const int streamsSrcId, int streamsTargetIdx) {
    if (streamsSrcId >= targetsCount)
        return NO_TARGET_TO_SCHEDULE;
    int nextDecFileIdx = streamsTargetIdx + (streamsTargetIdx < streamsSrcId ? 1 : 0);
    if (params->appendCommand) {
        while (refExtLoadedPosArr[nextDecFileIdx + 1] - refExtLoadedPosArr[nextDecFileIdx] <= 1) {
            isTargetScheduled[nextDecFileIdx] = true;
            if (++nextDecFileIdx >= targetsCount)
                return NO_TARGET_TO_SCHEDULE;
        }
    } else {
        size_t tmp = streamsTargetIdx < streamsSrcId ?
                     namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tIdNamesPos[streamsSrcId]) + 1
                                                     : tIdNamesPos[streamsSrcId];
        if (isFilterListActive) {
            while (!isTargetSelectedPtr[nextDecFileIdx]) {
                if (++nextDecFileIdx >= targetsCount)
                    return NO_TARGET_TO_SCHEDULE;
            }
        } else {
            size_t nameMatch = namesStr.find(params->masterFilterPattern, tmp);
            if (nameMatch == std::string::npos)
                return NO_TARGET_TO_SCHEDULE;
            while ((tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, tmp) + 1) < nameMatch)
                nextDecFileIdx++;
        }
    }
    if (streamsTargetIdx > streamsSrcId)
        copyStreamsPositions(streamsSrcId, streamsTargetIdx);
    int endFindTargetIdx =
            nextDecFileIdx < lazyDecompressionTargetsGuard ? nextDecFileIdx : lazyDecompressionTargetsGuard - 1;
    if (streamsTargetIdx + 1 < lazyDecompressionTargetsGuard)
        findStreamsPositions(streamsSrcId, streamsTargetIdx + 1, endFindTargetIdx);
    return nextDecFileIdx;
}

vector<int> tIdQueue;

template<bool lazyMode> void MBGC_Decoder<lazyMode>::scheduleDependencyChain(int64_t tId) {
    if (tId < lowestContributingNotDecodedTarget) {
        isTargetScheduled[tId] = true;
        return;
    }
    bool pairedGapFlag[MBGC_Params::MAX_GAP_DEPTH] {};
    int gapCurIdx = 0;
    uint64_t matchSrcPos = 0;
    tIdQueue.reserve(targetsCount);
    tIdQueue.clear();
    tIdQueue.emplace_back(tId);
    for (int i = 0; i < tIdQueue.size(); i++) {
        tId = tIdQueue[i];
        isTargetScheduled[tId] = true;
        scheduledTargetsCount++;
        if (tId == 0 && (params->enableExtensionsWithMismatches || params->gapDepthOffsetEncoding))
            continue;
        uint8_t *gapDeltaPtr = tIdGapDeltaPtr[tId];
        WrapperStrStream mapOffSrc = tIdMapOffSrc[tId];
        uint8_t *mapOff5thBytePtr = tIdMapOff5thBytePtr[0] != nullptr ? tIdMapOff5thBytePtr[tId] : nullptr;
        uint32_t *mapLenPtr = tIdMapLenPtr[tId];
        size_t seqEndPos = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, tIdLiteralPos[tId]);
        uint64_t markPos = literalStr.find(MBGC_Params::MATCH_MARK, tIdLiteralPos[tId]);
        int seqLeft = seqsCountArr[tId];
        while (markPos != std::string::npos && seqLeft) {
            if (pairedGapFlag[gapCurIdx]) {
                pairedGapFlag[gapCurIdx] = false;
            } else {
                matchSrcPos = 0;
                PgHelpers::readValue<uint32_t>(mapOffSrc, (uint32_t &) matchSrcPos);
                if (refTotalLength > UINT32_MAX) {
                    uint8_t tmp = *mapOff5thBytePtr++;
                    matchSrcPos += (((uint64_t) tmp) << 32);
                }
                while (matchSrcPos + (refTotalLength - REF_SHIFT) < refExtLoadedPosArr[tId])
                    matchSrcPos += refTotalLength - REF_SHIFT;
                auto startItr = refExtLoadedPosArr.begin() + lowestContributingNotDecodedTarget;
                if (matchSrcPos >= *startItr) {
                    auto endItr = refExtLoadedPosArr.begin() + tId + 1;
                    auto itr = std::upper_bound(startItr, endItr, matchSrcPos);
                    int mapIdx = std::distance(refExtLoadedPosArr.begin(), itr) - 1;
                    if (mapIdx < tId) {
                        if (highestTargetMapping[tId] < mapIdx)
                            highestTargetMapping[tId] = mapIdx;
                        if (!isTargetProcessed[mapIdx] && !isTargetScheduled[mapIdx]) {
                            isTargetScheduled[mapIdx] = true;
                            tIdQueue.emplace_back(mapIdx);
                        }
                    }
                }
                if (!params->enableExtensionsWithMismatches && !params->gapDepthOffsetEncoding) {
                    matchSrcPos += refTotalLength - REF_SHIFT;
                    auto startItr = refExtLoadedPosArr.begin() + tId + 1;
                    auto endItr = refExtLoadedPosArr.end();
                    auto itr = std::upper_bound(startItr, endItr, matchSrcPos);
                    int depIdx = std::distance(refExtLoadedPosArr.begin(), itr) - 1;
                    size_t matchEndPos = matchSrcPos + *(mapLenPtr++);
                    while (depIdx < targetsCount && refExtLoadedPosArr[depIdx] < matchEndPos) {
                        if (tId < depIdx && highestTargetMapping[depIdx] < tId)
                            highestTargetMapping[depIdx] = tId;
                        depIdx++;
                    }
                }
            }
            markPos = literalStr.find(MBGC_Params::MATCH_MARK, markPos + 1);
            uint8_t gapDelta = 0;
            if (params->gapDepthOffsetEncoding && markPos < seqEndPos && (gapDelta = *gapDeltaPtr++)) {
                int gapIdx = gapCurIdx;
                while (gapDelta) {
                    gapIdx = (gapIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;
                    if (!pairedGapFlag[gapIdx])
                        gapDelta--;
                }
                pairedGapFlag[gapIdx] = true;
            }
            gapCurIdx = (gapCurIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;
            while (seqLeft && markPos >= seqEndPos) {
                seqEndPos = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, seqEndPos + 1);
                seqLeft--;
            }
        }
    }
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::scheduleParallelLazyDecompression() {
    int streamsFileIdx = 0;
    int tId = 0;
    int firstOverwritingRefTarget = lazyDecompressionTargetsGuard;
    lazyDecompressionTargetsGuard = targetsCount;
    if (params->appendCommand) {
        int64_t totalLoadedRefLength = refExtLoadedPosArr[targetsCount];
        int64_t recoverFinalRefStateStartPos = totalLoadedRefLength - refTotalLength;
        if (recoverFinalRefStateStartPos > 0) {
            auto startItr = refExtLoadedPosArr.begin();
            auto endItr = refExtLoadedPosArr.end();
            auto itr = std::upper_bound(startItr, endItr, recoverFinalRefStateStartPos);
            tId = std::distance(refExtLoadedPosArr.begin(), itr) - 1;
            findStreamsPositions(1, 1, tId);
        }
    }
    while (tId != NO_TARGET_TO_SCHEDULE && tId < lazyDecompressionTargetsGuard) {
        scheduleDependencyChain(tId);
        streamsFileIdx = tId;
        tId = prepareStreamsForNextFileToDecompress(tId + 1, streamsFileIdx);
    }
    int _scheduledTargetsGuard = streamsFileIdx + 1;
    if (tId >= lazyDecompressionTargetsGuard) {
        for(int i = 0; i < lazyDecompressionTargetsGuard; i++)
            if (!isTargetScheduled[i] &&
                (i == lazyDecompressionTargetsGuard - 1 || refExtLoadedPosArr[i + 1] - refExtLoadedPosArr[i] > 1))
                scheduleDependencyChain(i);
        _scheduledTargetsGuard = lazyDecompressionTargetsGuard;
    }
    lazyDecompressionTargetsGuard = _scheduledTargetsGuard;
    vector<int> scheduledTargets;
    scheduledTargets.reserve(scheduledTargetsCount);
    vector<int> revIdxOfST(_scheduledTargetsGuard, -1);
    for(int i = 0; i < _scheduledTargetsGuard; i++)
        if (isTargetScheduled[i]) {
            revIdxOfST[i] = scheduledTargets.size();
            scheduledTargets.emplace_back(i);
        }
    if (params->enableExtensionsWithMismatches || params->gapDepthOffsetEncoding) {
        for(int i = 1; i < scheduledTargets.size(); i++) {
            int tId = scheduledTargets[i];
            if (tId >= firstOverwritingRefTarget)
                highestTargetMapping[tId] = scheduledTargets[i - 1];
        }
    }
    vector<int> orderedTargets = scheduledTargets;
    std::sort(orderedTargets.begin(), orderedTargets.end(),
              [this](const int& a, const int& b) -> bool { return highestTargetMapping[a] < highestTargetMapping[b]
              || (highestTargetMapping[a] == highestTargetMapping[b] && b < a); } );
    isTargetForParallelTaskThreadDecoding.resize(_scheduledTargetsGuard, false);
    targetsForParallelTaskThreadDecoding.reserve(scheduledTargetsCount / 3);
    for(int i = 0; i < scheduledTargetsCount; i++) {
        int tIdx = orderedTargets[i];
        int htmIdx = highestTargetMapping[tIdx];
        if (tIdx > 0 && (htmIdx < 0 || revIdxOfST[tIdx] - revIdxOfST[htmIdx] > 1)) {
            isTargetForParallelTaskThreadDecoding[tIdx] = true;
            targetsForParallelTaskThreadDecoding.emplace_back(tIdx);
        };
    }
    this->targetsDecodedBeforeScheduling = this->masterTargetsStats;
    this->lowestContributingNotDecodedTargetBeforeScheduling = this->lowestContributingNotDecodedTarget;
    this->scheduledTargetsGuard = _scheduledTargetsGuard;
    this->lowestNotClaimedParallelIdx = 0;
    *PgHelpers::devout << "finished scheduling " << scheduledTargets.size() << " (up to " << _scheduledTargetsGuard << ") targets - " << PgHelpers::time_millis() << " [ms]" << endl;
}

template<bool lazyMode> void MBGC_Decoder<lazyMode>::printSchedulingStats() {
    const int HIST_MAX = 32;
    vector<int> hist(HIST_MAX, 0);
    for(int i = 0; i < targetsCount; i++) {
        if (isTargetScheduled[i]) {
            *PgHelpers::devout << i << "\t" << highestTargetMapping[i] << "\t";
            int j = i;
            int c = 0;
            while (j-- && j > highestTargetMapping[i]) {
                if (isTargetScheduled[j])
                    c++;
            }
            *PgHelpers::devout << c << endl;
            if (c > HIST_MAX)
                c = HIST_MAX - 1;
            hist[c]++;
        }
    }
    int sum = 0;
    for(int i = 0; i < HIST_MAX; i++) {
        *PgHelpers::devout << i << "\t" << hist[i] << endl;
        sum += hist[i];
    }
    *PgHelpers::devout << sum << endl;
}

template<bool lazyMode>
void MBGC_Decoder<lazyMode>::iterateInit() {
    sequentialDecoding = true;
    decodeInit();
    masterThreadTargetIdx = 0;
}

template<bool lazyMode>
bool MBGC_Decoder<lazyMode>::iterateNext(const char* name, uint8_t* &res, size_t &size) {
    bool stopFlag = false;
    bool status;
    while (!stopFlag) {
#pragma omp critical(iteration)
        {
            while (iterCurFilename.empty()) {
                if (masterThreadTargetIdx >= targetsCount) {
                    status = false;
                    stopFlag = true;
                    break;
                }
                extractNextSequentially(MASTER_THREAD_ID);
                if (filesCount == 1 && masterThreadTargetIdx == targetsCount)
                    moveToFileIfSelected(g0name, threadGzOut[MASTER_THREAD_ID], MASTER_THREAD_ID);
            }
            if (iterCurFilename == name) {
                size = iterCurFileContent.size();
                res = (uint8_t *) malloc(size);
                memcpy(res, iterCurFileContent.data(), size);
                iterCurFileContent.clear();
                iterCurFilename.clear();
                status = true;
                stopFlag = true;
            }
        }
        if (!stopFlag)
            nanosleep(SLEEP_TIME, nullptr);
    }
    return status;
}

template class MBGC_Decoder<false>;
template class MBGC_Decoder<true>;

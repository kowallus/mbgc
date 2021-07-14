#include "MBGC_Decoder.h"

#include "../coders/CodersLib.h"
#include "../libs/asmlib.h"

#include <omp.h>
#include <numeric>

#ifdef DEVELOPER_BUILD
#include "../utils/libdeflate_wrapper.h"
#endif

MBGC_Decoder::MBGC_Decoder(MBGC_Params *mbgcParams): params(mbgcParams) {

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

void MBGC_Decoder::writeDNA(const char *sequence, int64_t length) {
    uint32_t pos = 0;
#ifdef DEVELOPER_BUILD
    while (pos < length - 80 && !params->disableDNAformatting) {
#else
    while (pos < length - 80) {
#endif
        A_append(outBuffer, sequence + pos, 80);
        outBuffer.push_back('\n');
        pos += 80;
    }
    if (length - pos > 0) {
        A_append(outBuffer, sequence + pos, length - pos);
        outBuffer.push_back('\n');
    }
}

void MBGC_Decoder::decodeHeader(string& headerTemplate) {
    size_t tmp, tPos = 0;
    while ((tmp = headerTemplate.find(MBGC_Params::MATCH_MARK, tPos)) != std::string::npos) {
        A_append(outBuffer, headerTemplate, tPos, tmp - tPos);
        tPos = tmp + 1;
        tmp = headersStr.find(MBGC_Params::MATCH_MARK, headersPos);
        A_append(outBuffer, headersStr, headersPos, tmp - headersPos);
        headersPos = tmp + 1;
    }
    A_append(outBuffer, headerTemplate, tPos, headerTemplate.size() - tPos);
}

bool MBGC_Decoder::moveToFile(const string& filename, string& src, const int thread_no) {
    if (filename.find(params->filterPattern) == string::npos) {
        src.clear();
        return true;
    }
    string filepath = params->outputPath + filename;
#ifdef DEVELOPER_BUILD
    if (params->validationMode) {
        string tmpfile = filepath;
        {
            if (!fstream(filepath)) {
                tmpfile += ".gz";
                if (!fstream(tmpfile)) {
                    fprintf(stderr, "Cannot find %s for validation.\n", filepath.c_str());
                    exit(EXIT_FAILURE);
                }
            }
        }
        gzFile file = gzopen(tmpfile.c_str(), "r");
        bool ok = true;
        if (src.size() != file.size) {
            if (params->invalidFilesCount < 100)
                *PgHelpers::logout << "Validation ERROR: ~" << (params->invalidFilesCount + params->validFilesCount) <<
                ". " << filepath << " size differ (" << src.size() << " instead of "
                    << file.size << ")" << endl;
            ok = false;
        }
        if (ok && memcmp(src.data(), file.out, file.size) != 0) {
            if (params->invalidFilesCount < 100)
                *PgHelpers::logout << "Validation ERROR: ~"  << (params->invalidFilesCount + params->validFilesCount) <<
                ". " << filepath << " contents differ." << endl;
            ok = false;
        }
        src.clear();
        gzclose(file);
        return ok;
    }
#endif
    PgHelpers::createFolders(filepath);
    fstream fstr(filepath, ios::out | ios::binary | ios::trunc);
    if (!fstr)
        return false;
    PgHelpers::writeArray(fstr, (void *) src.data(), src.size());
    fstr.close();
    extractedFilesCount[thread_no]++;
    src.clear();
    return true;
}

void MBGC_Decoder::decodeReference(const string &name) {
    size_t tmp;
    rcStart = 0;
    refStr[refPos++] = 0;
    size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, hTemplatesPos);
    string headerTemplate(headersTemplates, hTemplatesPos, hTemplateEnd - hTemplatesPos);
    uint32_t seqCount;
    PgHelpers::readValue<uint32_t>(seqsCountSrc, seqCount, false);
    for (uint32_t i = 0; i < seqCount; i++) {
        outBuffer.push_back('>');
        decodeHeader(headerTemplate);
        outBuffer.push_back('\n');
        tmp = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, literalPos);
        writeDNA(literalStr.data() + literalPos, tmp - literalPos);
        A_memcpy((char *) refStr.data() + refPos, literalStr.data() + literalPos, tmp - literalPos);
        refPos += tmp - literalPos;
        rcStart += tmp - literalPos;
        literalPos = tmp + 1;
    }
    headersPos++;
    hTemplatesPos = hTemplateEnd + 1;
    bool ok = moveToFile(name, outBuffer, 0);
#ifdef DEVELOPER_BUILD
    if (ok)
        params->validFilesCount++;
    else
        params->invalidFilesCount++;
#endif
    if (!params->dontUseRCinReference()) {
        char *refPtr = (char *) refStr.data() + 1;
        PgHelpers::reverseComplement(refPtr, rcStart, refPtr + rcStart);
        refPos += rcStart;
    }
}

uint32_t MBGC_Decoder::decodeSequenceAndReturnUnmatchedChars(string &dest) {
    size_t seqEnd = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, literalPos);

    dest.clear();
    uint32_t unmatchedChars = 0;
    uint32_t minMatchLength = 0;

    uint64_t markPos = 0;
    while ((markPos = literalStr.find(MBGC_Params::MATCH_MARK, literalPos)) != std::string::npos) {
        if (markPos > seqEnd)
            break;
        A_append(dest, literalStr, literalPos, markPos - literalPos);
        unmatchedChars += markPos - literalPos;
        literalPos = markPos + 1;
        uint32_t matchSrcPos = 0;
        PgHelpers::readValue<uint32_t>(mapOffSrc, matchSrcPos, false);
        uint32_t matchLength = 0;
        PgHelpers::readUIntWordFrugal<uint32_t>(mapLenSrc, matchLength);
        matchLength += minMatchLength;
        A_append(dest, refStr, matchSrcPos, matchLength);
    }
    A_append(dest, literalStr, literalPos, seqEnd - literalPos);
    unmatchedChars += seqEnd - literalPos;
    literalPos = seqEnd + 1;
    return unmatchedChars;
}

string seqStr;

void MBGC_Decoder::decodeFile(uint8_t unmatchedFractionFactor) {
    seqStr.clear();
    size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, hTemplatesPos);
    string headerTemplate(headersTemplates, hTemplatesPos, hTemplateEnd - hTemplatesPos);
    uint32_t seqCount;
    PgHelpers::readValue<uint32_t>(seqsCountSrc, seqCount, false);
    size_t refLockPos = this->refTotalLength;
    if (matchingLocksPosSrc.rdbuf()->in_avail())
        PgHelpers::readValue<size_t>(matchingLocksPosSrc, refLockPos, false);
    for (uint32_t i = 0; i < seqCount; i++) {
        outBuffer.push_back('>');
        decodeHeader(headerTemplate);
        outBuffer.push_back('\n');
        uint32_t unmatchedChars = decodeSequenceAndReturnUnmatchedChars(seqStr);
        writeDNA(seqStr.data(), seqStr.size());
        if (params->isContigProperForRefExtension(seqStr.size(), unmatchedChars, unmatchedFractionFactor)) {
            loadRef(seqStr.data(), seqStr.size(), refLockPos);
            if (!params->dontUseRCinReference()) {
                PgHelpers::reverseComplementInPlace((char *) seqStr.data(), seqStr.size());
                loadRef(seqStr.data(), seqStr.size(), refLockPos);
            }
        }
    }
    headersPos++;
    hTemplatesPos = hTemplateEnd + 1;
}

void MBGC_Decoder::loadRef(const char *seqText, size_t seqLength, size_t refLockPos)  {
    if (seqLength == 0)
        return;
    if (refPos == refTotalLength && refLockPos != refTotalLength) {
        refPos = 1;
    }
    size_t tmpEnd = refLockPos;
    size_t tmpLength = seqLength;
    size_t tmpMax = tmpEnd < refPos ? refTotalLength : tmpEnd;
    if (refPos + tmpLength > tmpMax) {
        tmpLength = tmpMax - refPos;
    }
    A_memcpy((char*) refStr.data() + refPos, seqText, tmpLength);
    refPos += tmpLength;
    seqText += tmpLength;
    seqLength = refPos == tmpEnd ? 0 : seqLength - tmpLength;
    this->loadRef(seqText, seqLength, refLockPos);
}

void MBGC_Decoder::extractFiles() {
    size_t tmp;
    size_t endGuard = namesStr.rfind(params->filterPattern);
    if (endGuard >= namesPos && endGuard != std::string::npos) {
        while ((tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos)) != std::string::npos) {
            decodeFile(unmatchedFractionFactors[fileIdx++]);
            string name = namesStr.substr(namesPos, tmp - namesPos);
            namesPos = tmp + 1;
            bool ok = moveToFile(name, outBuffer, 0);
#ifdef DEVELOPER_BUILD
            if (ok)
                params->validFilesCount++;
            else
                params->invalidFilesCount++;
#endif
            if (tmp > endGuard)
                break;
        }
    }
}

void MBGC_Decoder::writeFilesParallelTask(const int thread_no) {
#ifdef DEVELOPER_BUILD
    uint32_t validFilesCount = 0;
    uint32_t invalidFilesCount = 0;
#endif
    while(isDecoding || in[thread_no] != out[thread_no] ){
        if (in[thread_no] != out[thread_no]) {
            string &content = contentsBuf[thread_no][out[thread_no]];
            string &name = namesBuf[thread_no][out[thread_no]];
            bool ok = moveToFile(name, content, thread_no);
#ifdef DEVELOPER_BUILD
            if (ok)
                validFilesCount++;
            else
                invalidFilesCount++;
#endif
            out[thread_no] = (out[thread_no] + 1) % WRITING_BUFFER_SIZE;
        } else {
            nanosleep((const struct timespec[]){{0, 1000L}}, NULL);
        }
    }
#ifdef DEVELOPER_BUILD
#pragma omp critical
    {
        params->validFilesCount += validFilesCount;
        params->invalidFilesCount += invalidFilesCount;
    }
#endif
}

void MBGC_Decoder::extractFilesParallel() {
    size_t tmp;
    isDecoding = true;
#pragma omp parallel shared(tmp)
    {
#pragma omp single
        {
            int writingThreadsCount = omp_get_num_threads() - 1;
            extractedFilesCount.resize(writingThreadsCount, 0);
            in.resize(writingThreadsCount, 0);
            out.resize(writingThreadsCount, 0);
            namesBuf.resize(writingThreadsCount, vector<string>(WRITING_BUFFER_SIZE));
            contentsBuf.resize(writingThreadsCount, vector<string>(WRITING_BUFFER_SIZE));
            for (int i = 0; i < writingThreadsCount; i++) {
#pragma omp task
                {
                    writeFilesParallelTask(i);
                }
            }
            int thread_no = 0;
            size_t endGuard = namesStr.rfind(params->filterPattern);
            if (endGuard >= namesPos && endGuard != std::string::npos) {
                while ((tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos)) != std::string::npos) {
                    while (((in[thread_no] + 1) % WRITING_BUFFER_SIZE) == out[thread_no])
                        nanosleep((const struct timespec[]) {{0, 100L}}, NULL);
                    namesBuf[thread_no][in[thread_no]] = namesStr.substr(namesPos, tmp - namesPos);
                    namesPos = tmp + 1;
                    outBuffer = std::move(contentsBuf[thread_no][in[thread_no]]);
                    outBuffer.reserve(largestFileLength);
                    decodeFile(unmatchedFractionFactors[fileIdx++]);
                    contentsBuf[thread_no][in[thread_no]] = std::move(outBuffer);
                    in[thread_no] = (in[thread_no] + 1) % WRITING_BUFFER_SIZE;
                    thread_no = (thread_no + 1) % writingThreadsCount;
                    if (tmp > endGuard)
                        break;
                }
            }
            isDecoding = false;
        }
    }
}

void MBGC_Decoder::decode() {
    fstream fin(params->archiveFileName);
    if (!fin) {
        fprintf(stderr, "Cannot open archive %s\n", params->archiveFileName.c_str());
        exit(EXIT_FAILURE);
    }
    readParamsAndStats(fin);
    seqStr.reserve(largestContigLength);
    refStr.resize(refTotalLength);
    outBuffer.reserve(largestFileLength);
    vector<string*> destStrings;
    destStrings.push_back(&namesStr);
    string seqsCount;
    destStrings.push_back(&seqsCount);
    destStrings.push_back(&headersTemplates);
    destStrings.push_back(&headersStr);
    string tmpRefExtStr;
    destStrings.push_back(&tmpRefExtStr);
    destStrings.push_back(&literalStr);
    string mapOffStream, mapLenStream, matchingLocksPosStream;
    if (mbgcVersionMajor > 1 || (mbgcVersionMajor == 1 && mbgcVersionMinor >= 1))
        destStrings.push_back(&matchingLocksPosStream);
    destStrings.push_back(&mapOffStream);
    destStrings.push_back(&mapLenStream);
    readCompressedCollectiveParallel(fin, destStrings);
    unmatchedFractionFactors.resize(tmpRefExtStr.size());
    memcpy((void*) unmatchedFractionFactors.data(), (void*) tmpRefExtStr.data(), tmpRefExtStr.size());
    seqsCountSrc.str(seqsCount);
    mapOffSrc.str(mapOffStream);
    mapLenSrc.str(mapLenStream);
    matchingLocksPosSrc.str(matchingLocksPosStream);
    cout << "loaded archive - " << PgHelpers::time_millis() << " [ms]" << endl;
    size_t tmp;

    tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
    string name = namesStr.substr(namesPos, tmp - namesPos);
    namesPos = tmp + 1;
    extractedFilesCount.push_back(0);
    decodeReference(name);

    if (PgHelpers::numberOfThreads == 1)
        extractFiles();
    else
        extractFilesParallel();

#ifdef DEVELOPER_BUILD
    if (params->validationMode) {
        *PgHelpers::logout << "Validation" << (params->validFilesCount == filesCount ? "" : " ERROR") <<
                           ": correctly decoded " << params->validFilesCount << " out of " << filesCount << " files."
                           << endl;
        if (params->invalidFilesCount)
            fprintf(stderr, "Validation ERROR: errors in contents of %d decoded files.\n", params->invalidFilesCount);
    }
#endif
    uint32_t totalExtractedFilesCount = accumulate(extractedFilesCount.begin(), extractedFilesCount.end(), 0);
    cout << "extracted " << totalExtractedFilesCount << (totalExtractedFilesCount == 1 ? " file" : " files") << endl;
}

void MBGC_Decoder::readParamsAndStats(fstream &fin) {
    for (int i = 0; i < strlen(params->MBGC_HEADER); i++) {
        if (params->MBGC_HEADER[i] != fin.get()) {
            fprintf(stderr, "Error processing header.\n");
            exit(EXIT_FAILURE);
        }
    }
    char ch = fin.get();
    if (ch != params->MBGC_VERSION_MODE) {
        fprintf(stderr, "Error processing header.\n");
        exit(EXIT_FAILURE);
    }
    mbgcVersionMajor = fin.get();
    mbgcVersionMinor = fin.get();
    mbgcVersionRevision = fin.get();
    ch = fin.get(); // numberOfThreads
    params->coderLevel = fin.get();
    params->sequentialMatching = (bool) fin.get();
    params->k = fin.get();
    PgHelpers::readValue<int>(fin, params->k1, false);
    PgHelpers::readValue<int>(fin, params->k2, false);
    if (mbgcVersionMajor > 1 || (mbgcVersionMajor == 1 && mbgcVersionMinor >= 1))
        PgHelpers::readValue<uint8_t>(fin, params->skipMargin, false);
    PgHelpers::readValue<int>(fin, params->referenceFactor, false);
    params->refExtensionStrategy = fin.get();
    if (params->useLiteralsinReference()) {
        fprintf(stderr, "Literals ref extension strategy decompression not implemented.\n");
        exit(EXIT_FAILURE);
    }
    if (params->splitContigsIntoBlocks()) {
        fprintf(stderr, "Split contigs into blocks strategy decompression not implemented.\n");
        exit(EXIT_FAILURE);
    }
    int tmp;
    if (params->usesCombinedRefExtensionStrategy()) {
#ifndef DEVELOPER_BUILD
        fprintf(stderr, "Combined ref extension strategy available only in developer build.\n");
        exit(EXIT_FAILURE);
#endif
        PgHelpers::readValue<int>(fin, tmp, false);
        if (tmp != params->MINIMAL_UNMATCHED_LENGTH_FACTOR) {
            fprintf(stderr, "Invalid minimal unmatched length factor: %d (expected %d).\n", tmp,
                    params->MINIMAL_UNMATCHED_LENGTH_FACTOR);
            exit(EXIT_FAILURE);
        }
    }
    PgHelpers::readValue<uint32_t>(fin, filesCount, false);
    PgHelpers::readValue<size_t>(fin, totalFilesLength, false);
    PgHelpers::readValue<int>(fin, rcStart, false);
    PgHelpers::readValue<uint32_t>(fin, largestFileLength, false);
    PgHelpers::readValue<uint32_t>(fin, largestContigLength, false);
    PgHelpers::readValue<uint32_t>(fin, refTotalLength, false);
}

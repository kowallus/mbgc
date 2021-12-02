#include "MBGC_Decoder.h"

#include "../coders/CodersLib.h"
#include "../libs/asmlib.h"

#include <fstream>
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
    const int dnaLineLength = params->dnaLineLength;
    uint32_t pos = 0;
    while (pos < length - dnaLineLength && params->enableDNAformatting) {
        A_append(outBuffer, sequence + pos, dnaLineLength);
        outBuffer.push_back('\n');
        pos += dnaLineLength;
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

bool MBGC_Decoder::moveToFile(const string& filename, string& src, const int thread_no, bool append) {
    if (filename.find(params->filterPattern) == string::npos) {
        src.clear();
        return true;
    }
    string filepath = params->outputPath + filename;
#ifdef DEVELOPER_BUILD
    if (params->validationMode) {
        if (filesCount == 1 && unmatchedFractionFactors.size() > 1) {
            fprintf(stderr, "Validation in single fasta mode archive not supported.\n");
            exit(EXIT_FAILURE);
        }
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
        gzFile file = gzopen(tmpfile.c_str());
        bool ok = true;
        if (src.size() != file.size) {
            if (params->invalidFilesCount < 100)
                *PgHelpers::devout << "Validation ERROR: ~" << (params->invalidFilesCount + params->validFilesCount) <<
                                   ". " << filepath << " size differ (" << src.size() << " instead of "
                    << file.size << ")" << endl;
            ok = false;
        }
        if (ok && memcmp(src.data(), file.out, file.size) != 0) {
            if (params->invalidFilesCount < 100)
                *PgHelpers::devout << "Validation ERROR: ~" << (params->invalidFilesCount + params->validFilesCount) <<
                                   ". " << filepath << " contents differ." << endl;
            ok = false;
        }
        src.clear();
        gzclose(file);
        return ok;
    }
#endif
    if (params->outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS)
        PgHelpers::writeArray(cout, (void *) src.data(), src.size());
    else {
        PgHelpers::createFolders(filepath);
        fstream fstr(filepath, ios::out | ios::binary | (append ? ios::app : ios::trunc));
        if (!fstr) {
            fprintf(stderr, "Error: cannot create a file %s\n", filepath.c_str());
            return false;
        }
        PgHelpers::writeArray(fstr, (void *) src.data(), src.size());
        fstr.close();
    }
    extractedFilesCount[thread_no]++;
    src.clear();
    return true;
}

void MBGC_Decoder::initReference(const string &name) {
    size_t tmp;
    rcStart = 0;
    refStr[refPos++] = 0;
    if (literalStr.empty()) {
        literalPos = 0;
        return;
    }
    tmp = literalStr.find(MBGC_Params::SEQ_SEPARATOR_MARK, literalPos);
    A_memcpy((char *) refStr.data() + refPos, literalStr.data() + literalPos, tmp - literalPos);
    refPos += tmp - literalPos;
    rcStart += tmp - literalPos;
    literalPos = tmp + 1;
    if (!params->dontUseRCinReference()) {
        char *refPtr = (char *) refStr.data() + 1;
        PgHelpers::reverseComplement(refPtr, rcStart, refPtr + rcStart);
        refPos += rcStart;
    }
    if (filesCount == 1)
        bool ok = moveToFile(name, outBuffer, 0);
}

void MBGC_Decoder::decodeReference(const string &name) {
    size_t tmp;
    rcStart = 0;
    refStr[refPos++] = 0;
    size_t hTemplateEnd = headersTemplates.find(MBGC_Params::FILE_SEPARATOR_MARK, hTemplatesPos);
    string headerTemplate(headersTemplates, hTemplatesPos, hTemplateEnd - hTemplatesPos);
    uint32_t seqCount;
    PgHelpers::readValue<uint32_t>(seqsCountSrc, seqCount, false);
#ifdef DEVELOPER_BUILD
    if (params->concatHeadersAndSequencesMode) {
        outBuffer.push_back('>');
        for (uint32_t i = 0; i < seqCount; i++) {
            if (i) outBuffer.push_back('&');
            decodeHeader(headerTemplate);
        }
        outBuffer.push_back('\n');
    }
    for (uint32_t i = 0; i < seqCount; i++) {
        if (!params->concatHeadersAndSequencesMode) {
            outBuffer.push_back('>');
            decodeHeader(headerTemplate);
            outBuffer.push_back('\n');
        }
#else
    for (uint32_t i = 0; i < seqCount; i++) {
        outBuffer.push_back('>');
        decodeHeader(headerTemplate);
        outBuffer.push_back('\n');
#endif
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
        uint32_t matchSrc32bitPos = 0;
        PgHelpers::readValue<uint32_t>(mapOffSrc, matchSrc32bitPos, false);
        uint32_t matchLength = 0;
        PgHelpers::readUIntWordFrugal<uint32_t>(mapLenSrc, matchLength);
        matchLength += minMatchLength;
        if (refTotalLength <= UINT32_MAX)
            A_append(dest, refStr, matchSrc32bitPos, matchLength);
        else {
            uint8_t tmp;
            PgHelpers::readValue<uint8_t>(mapOff5thByteSrc, tmp, false);
            uint64_t matchSrcPos = (((uint64_t) tmp) << 32) + matchSrc32bitPos;
            A_append(dest, refStr, matchSrcPos, matchLength);
        }
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
#ifdef DEVELOPER_BUILD
    if (params->concatHeadersAndSequencesMode) {
        outBuffer.push_back('>');
        for (uint32_t i = 0; i < seqCount; i++) {
            if (i) outBuffer.push_back('&');
            decodeHeader(headerTemplate);
        }
        outBuffer.push_back('\n');
    }
    for (uint32_t i = 0; i < seqCount; i++) {
        if (!params->concatHeadersAndSequencesMode) {
            outBuffer.push_back('>');
            decodeHeader(headerTemplate);
            outBuffer.push_back('\n');
        }
#else
    for (uint32_t i = 0; i < seqCount; i++) {
        outBuffer.push_back('>');
        decodeHeader(headerTemplate);
        outBuffer.push_back('\n');
#endif
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
        while (fileIdx < unmatchedFractionFactors.size()) {
            decodeFile(unmatchedFractionFactors[fileIdx++]);
            if (filesCount == 1) namesPos = 0;
            tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
            string name = namesStr.substr(namesPos, tmp - namesPos);
            namesPos = tmp + 1;
            bool ok = moveToFile(name, outBuffer, 0, filesCount == 1);
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
    istream* in;
    if (decompressFromStdin())
        in = &cin;
    else {
        in = new fstream(params->archiveFileName);
        if (!*in) {
            fprintf(stderr, "Cannot open archive %s\n", params->archiveFileName.c_str());
            exit(EXIT_FAILURE);
        }
    }
    readParamsAndStats(*in);
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
    string mapOffStream, mapOff5thByteStream, mapLenStream, matchingLocksPosStream;
    if (mbgcVersionMajor > 1 || (mbgcVersionMajor == 1 && mbgcVersionMinor >= 1))
        destStrings.push_back(&matchingLocksPosStream);
    destStrings.push_back(&mapOffStream);
    if (refTotalLength > UINT32_MAX)
        destStrings.push_back(&mapOff5thByteStream);
    destStrings.push_back(&mapLenStream);
    readCompressedCollectiveParallel(*in, destStrings);
    unmatchedFractionFactors.resize(tmpRefExtStr.size());
    memcpy((void*) unmatchedFractionFactors.data(), (void*) tmpRefExtStr.data(), tmpRefExtStr.size());
    seqsCountSrc.str(seqsCount);
    mapOffSrc.str(mapOffStream);
    mapOff5thByteSrc.str(mapOff5thByteStream);
    mapLenSrc.str(mapLenStream);
    matchingLocksPosSrc.str(matchingLocksPosStream);
    if (!decompressFromStdin())
        delete(in);
    *PgHelpers::appout << "loaded archive - " << PgHelpers::time_millis() << " [ms]" << endl;
    size_t tmp;

    tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
    string name = namesStr.substr(namesPos, tmp - namesPos);
    extractedFilesCount.push_back(0);
    if (params->sequentialMatching && (mbgcVersionMajor > 1 || (mbgcVersionMajor == 1 && mbgcVersionMinor >= 2)))
        initReference(name);
    else {
        namesPos = tmp + 1;
        decodeReference(name);
    }
    if (filesCount == 1 || params->outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS)
        PgHelpers::numberOfThreads = 1;
    if (PgHelpers::numberOfThreads == 1)
        extractFiles();
    else
        extractFilesParallel();

#ifdef DEVELOPER_BUILD
    if (params->validationMode) {
        *PgHelpers::devout << "Validation" << (params->validFilesCount == filesCount ? "" : " ERROR") <<
                           ": correctly decoded " << params->validFilesCount << " out of " << filesCount << " files."
                           << endl;
        if (params->invalidFilesCount)
            fprintf(stderr, "Validation ERROR: errors in contents of %d decoded files.\n", params->invalidFilesCount);
    }
#endif
    uint32_t totalExtractedFilesCount = accumulate(extractedFilesCount.begin(), extractedFilesCount.end(), 0);
    if (filesCount == 1 && totalExtractedFilesCount)
        totalExtractedFilesCount = 1;
    *PgHelpers::appout << "extracted " << totalExtractedFilesCount << (totalExtractedFilesCount == 1 ? " file" : " files") << endl;
}

void MBGC_Decoder::readParamsAndStats(istream &in) {
    for (int i = 0; i < strlen(MBGC_Params::MBGC_HEADER); i++) {
        if (MBGC_Params::MBGC_HEADER[i] != in.get()) {
            fprintf(stderr, "Error processing header.\n");
            exit(EXIT_FAILURE);
        }
    }
    char ch = in.get();
    if (ch != MBGC_Params::MBGC_VERSION_MODE) {
        fprintf(stderr, "Error processing header.\n");
        exit(EXIT_FAILURE);
    }
    mbgcVersionMajor = in.get();
    mbgcVersionMinor = in.get();
    mbgcVersionRevision = in.get();
    ch = in.get(); // numberOfThreads
    params->coderLevel = in.get();
    params->sequentialMatching = (bool) in.get();
    params->k = in.get();
    PgHelpers::readValue<int>(in, params->k1, false);
    PgHelpers::readValue<int>(in, params->k2, false);
    if (mbgcVersionMajor > 1 || (mbgcVersionMajor == 1 && mbgcVersionMinor >= 1))
        PgHelpers::readValue<uint8_t>(in, params->skipMargin, false);
    PgHelpers::readValue<int>(in, params->referenceFactor, false);
    params->refExtensionStrategy = in.get();
#ifdef DEVELOPER_BUILD
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
        PgHelpers::readValue<int>(in, tmp, false);
        if (tmp != MBGC_Params::MINIMAL_UNMATCHED_LENGTH_FACTOR) {
            fprintf(stderr, "Invalid minimal unmatched length factor: %d (expected %d).\n", tmp,
                    MBGC_Params::MINIMAL_UNMATCHED_LENGTH_FACTOR);
            exit(EXIT_FAILURE);
        }
    }
#endif
    PgHelpers::readValue<uint32_t>(in, filesCount, false);
    PgHelpers::readValue<size_t>(in, totalFilesLength, false);
    PgHelpers::readValue<int>(in, rcStart, false);
    PgHelpers::readValue<uint32_t>(in, largestFileLength, false);
    PgHelpers::readValue<uint32_t>(in, largestContigLength, false);
    if (mbgcVersionMajor > 1 || (mbgcVersionMajor == 1 && mbgcVersionMinor >= 2))
        PgHelpers::readValue<uint64_t>(in, refTotalLength, false);
    else {
        uint32_t tmp;
        PgHelpers::readValue<uint32_t>(in, tmp, false);
        refTotalLength = tmp;
    }

}

bool MBGC_Decoder::decompressFromStdin() {
    return params->archiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS;
}
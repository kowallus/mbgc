#include "MBGC_Decoder.h"

#include "../coders/CodersLib.h"
#include "../libs/asmlib.h"

#include <omp.h>

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

bool MBGC_Decoder::moveToFile(const string& filepath, string& src) {
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
    src.clear();
    return true;
}

void MBGC_Decoder::decodeReference(const string &name) {
    size_t tmp;
    rcStart = 0;
    refStr.push_back(0);
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
        A_append(refStr, literalStr, literalPos, tmp - literalPos);
        rcStart += tmp - literalPos;
        literalPos = tmp + 1;
    }
    headersPos++;
    hTemplatesPos = hTemplateEnd + 1;
    bool ok = moveToFile(params->outputPath + name, outBuffer);
#ifdef DEVELOPER_BUILD
    if (ok)
        params->validFilesCount++;
    else
        params->invalidFilesCount++;
#endif
    if (!params->dontUseRCinReference()) {
        refStr.resize(rcStart * 2 + 1);
        char *refPtr = (char *) refStr.data() + 1;
        PgHelpers::reverseComplement(refPtr, rcStart, refPtr + rcStart);
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
    for (uint32_t i = 0; i < seqCount; i++) {
        outBuffer.push_back('>');
        decodeHeader(headerTemplate);
        outBuffer.push_back('\n');
        uint32_t unmatchedChars = decodeSequenceAndReturnUnmatchedChars(seqStr);
        writeDNA(seqStr.data(), seqStr.size());
        if (params->isContigProperForRefExtension(seqStr.size(), unmatchedChars, unmatchedFractionFactor)) {
            size_t refExtLength = refStr.size() + seqStr.size() <= refTotalLength? seqStr.size()
                     : refTotalLength - refStr.size();
            A_append(refStr, seqStr.data(), refExtLength);
            if (!params->dontUseRCinReference()) {
                PgHelpers::reverseComplementInPlace((char *) seqStr.data(), seqStr.size());
                refExtLength = refStr.size() + seqStr.size() <= refTotalLength? seqStr.size()
                        : refTotalLength - refStr.size();
                A_append(refStr, seqStr.data(), refExtLength);
            }
        }
    }
    headersPos++;
    hTemplatesPos = hTemplateEnd + 1;
}

void MBGC_Decoder::extractFiles() {
    size_t tmp;
    while ((tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos)) != std::string::npos) {
        decodeFile(unmatchedFractionFactors[fileIdx++]);
        string name = namesStr.substr(namesPos, tmp - namesPos);
        namesPos = tmp + 1;
        bool ok = moveToFile(params->outputPath + name, outBuffer);
#ifdef DEVELOPER_BUILD
        if (ok)
            params->validFilesCount++;
        else
            params->invalidFilesCount++;
#endif
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
            bool ok = moveToFile(params->outputPath + name, content);
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
            while ((tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos)) != std::string::npos) {
                while (((in[thread_no] + 1) % WRITING_BUFFER_SIZE) == out[thread_no])
                    nanosleep((const struct timespec[]){{0, 100L}}, NULL);
                namesBuf[thread_no][in[thread_no]] = namesStr.substr(namesPos, tmp - namesPos);
                namesPos = tmp + 1;
                outBuffer = std::move(contentsBuf[thread_no][in[thread_no]]);
                outBuffer.reserve(largestFileLength);
                decodeFile(unmatchedFractionFactors[fileIdx++]);
                contentsBuf[thread_no][in[thread_no]] = std::move(outBuffer);
                in[thread_no] = (in[thread_no] + 1) % WRITING_BUFFER_SIZE;
                thread_no = (thread_no + 1) % writingThreadsCount;
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
    refStr.reserve(refTotalLength);
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
    string mapOff, mapLen;
    destStrings.push_back(&mapOff);
    destStrings.push_back(&mapLen);
    readCompressedCollectiveParallel(fin, destStrings);
    unmatchedFractionFactors.resize(tmpRefExtStr.size());
    memcpy((void*) unmatchedFractionFactors.data(), (void*) tmpRefExtStr.data(), tmpRefExtStr.size());
    seqsCountSrc.str(seqsCount);
    mapOffSrc.str(mapOff);
    mapLenSrc.str(mapLen);
    cout << "loaded archive - " << PgHelpers::time_millis() << " [ms]" << endl;
    size_t tmp;

    tmp = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
    string name = namesStr.substr(namesPos, tmp - namesPos);
    namesPos = tmp + 1;
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
    params->k = fin.get();;
    PgHelpers::readValue<int>(fin, params->k1, false);
    PgHelpers::readValue<int>(fin, params->k2, false);
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

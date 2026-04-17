#ifndef MBGC_MBGC_ENCODER_H
#define MBGC_MBGC_ENCODER_H

#include <iostream>

#include "../utils/helper.h"

#include "../matching/MultipleGenomeMatchingProcessor.h"
#include "MBGC_Params.h"
#include "../coders/ContextAwareMismatchesCoder.h"
#include "MBGC_Decoder.h"


class MBGC_Encoder : public MultipleGenomeMatchingProcessor {
private:
    MBGC_Params *params;

    ContextAwareMismatchesCoder* mismatchesCoder = &ContextAwareMismatchesCoder::defaultInstance;

    ostringstream seqsCountsDest, dnaLineLengthsDest;

    uint32_t appendedFilesCount = 0;

    string headersStr;
    string headersTemplates;

    string namesStr;

    string locksPosStream;
    string mapOffStream;
    string mapLenStream;
    string singularGapFlagsStream;
    string gapFlagsStream;
    string gapMismatchesFlagsStream;

    vector<string> fileHeadersTemplates;
    vector<string> fileHeaders;

    vector<string> targetLiterals;
    vector<ostringstream> targetMapOffDests;
    vector<string> targetMapOff5thByte;
    vector<ostringstream> targetMapLenDests;
    vector<string> targetGapDeltas;
    vector<string> targetGapMismatchesFlags;

    ostringstream refExtSizeDest;
    vector<size_t> refExtLoadedPosArr;

    vector<vector<string>*> usedByteStreams;
    void enrollByteStream(vector<string>& stream) {
        usedByteStreams.push_back(&stream);
    }
    vector<pair<string&, vector<ostringstream>&>> usedStreamsWithDests;
    void enrollStream(string& stream, vector<ostringstream>& dests) {
        usedStreamsWithDests.emplace_back(pair<string& , vector<ostringstream>&>(stream, dests));
    }

    void print_invalid_kseq_status_message(int kseq_status) override;

    void setProteinsProfile() override;
    void initStreamsForG0Ref() override;
    void processG0RefContig(const char *seq, size_t len) override;

    void initProcessTargetsWithParallelIO() override;
    void processAfterTargetWithParallelIO(size_t matcherLoaderStartPos) override;
    void finalizeProcessTargetsWithParallelIO() override;

    void processAfterTarget(uint32_t targetIdx) override;
    void processAfterSequence(uint32_t targetIdx) override;

    void initProcessTarget(uint32_t targetIdx) override;
    void initParallelProcessing() override;
    void finalizeParallelProcessingInSingleFastaFileMode() override;
    void finalizeParallelProcessingOfTarget(uint32_t targetIdx, size_t matcherLoaderStartPos) override;

    void processFileName(string &fileName) override;
    void processTargetMeta(uint32_t seqCount, uint64_t dnaLineLength) override;
    void updateHeadersTemplate(uint32_t fileIndex, const string& header);
    void processHeader(uint32_t fileIndex, const string &header) override;

    size_t getMatchLoadedPos(size_t pos);

    void processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes);
    size_t processMatches(vector<PgTools::TextMatch>& textMatches, char *destStart, size_t destLen, int targetIdx,
        size_t matchingLockPos = SIZE_MAX) override;
    uint64_t extendMatchLeft(const char *destStart, uint64_t length, const TextMatch &match,
                         int targetIdx, size_t matchingLockPos,
                         uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches);
    uint64_t extendMatchRight(const char *gapStartPtr, const TextMatch &coreMatch, const TextMatch &match,
                              uint64_t &length, int targetIdx, bool isGap, bool gapStart, bool gapMiddle, bool gapEnd,
                              uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches);

    void applyTemplateToHeaders(uint32_t fileIdx) override;
    void prepareHeadersStreams();
    void prepareAndCompressStreams();
    void printAdditionalMatchingStats() override;
    void writeStats(ostream &out) const;

    void loadFileNames() override;

    template<bool lazyMode>
    void appendTemplate(MBGC_Decoder<lazyMode>* decoder, MBGC_Params& newParams);

    template<bool lazyMode>
    void repackTemplate(MBGC_Decoder<lazyMode>* decoder);

public:

    explicit MBGC_Encoder(MBGC_Params *mbgcParams);

    void encode();

    void append(MBGC_Params& newParams);

    void repack(MBGC_Params& inParams);

};


#endif //MBGC_MBGC_ENCODER_H


#ifndef PGTOOLS_SIMPLESEQUENCEMATCHER_H
#define PGTOOLS_SIMPLESEQUENCEMATCHER_H

#include "TextMatchers.h"

namespace PgTools {

    struct PgMatch;

    class SimpleSequenceMatcher {
    private:
        uint32_t targetMatchLength;

        TextMatcher* matcher = 0;
        const string& srcSeq;

        uint64_t destSeqLength;
        vector<TextMatch> textMatches;
        bool revComplMatching;

        void exactMatchSequence(string& destSeq, bool destSeqIsRef, uint32_t minMatchLength);

        void correctDestPositionDueToRevComplMatching();
        void resolveMappingCollisionsInTheSameText();

        string getTotalMatchStat(size_t totalMatchLength);

    public:
        SimpleSequenceMatcher(const string& srcSeq, uint32_t targetMatchLength,
                              uint32_t minMatchLength = UINT32_MAX);

        virtual ~SimpleSequenceMatcher();

        void markAndRemoveExactMatches(bool destSeqIsRef,
                                       string &destSeq, string &resMapOff, string& resMapLen,
                                       bool revComplMatching, uint32_t minMatchLength = UINT32_MAX);

        static void rcMatchSequence(string &sequence, string& rcMapOff, string& rcMapLen,
                                    size_t targetMatchLength, uint32_t minMatchLength = UINT32_MAX);

        static void restoreRCMatchedSequence(string &srcSequence, string& rcMapOff, string& rcMapLen, size_t orgSrcLen);

        static char MATCH_MARK;
    };
}

#endif //PGTOOLS_SIMPLESEQUENCEMATCHER_H

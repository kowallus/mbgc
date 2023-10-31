#include "SimpleSequenceMatcher.h"

#include "copmem/CopMEMMatcher.h"
#include "SlidingWindowSparseEMMatcher.h"
#include "../mbgccoder/MBGC_Params.h"


namespace PgTools {

    using namespace PgHelpers;

    SimpleSequenceMatcher::SimpleSequenceMatcher(const string &srcSeq, uint32_t targetMatchLength,
                                                 uint32_t minMatchLength)
            : srcSeq(srcSeq), targetMatchLength(targetMatchLength) {
        *devout << "Source sequence length: " << srcSeq.length() << endl;
        if (srcSeq.size() >= targetMatchLength) {
            matcher = new CopMEMMatcher(srcSeq.data(), srcSeq.length(), targetMatchLength, minMatchLength);
        }
    }

    SimpleSequenceMatcher::~SimpleSequenceMatcher() {
        if (matcher)
            delete (matcher);
    }

    void SimpleSequenceMatcher::exactMatchSequence(string &destSeq, bool destSeqIsRef, uint32_t minMatchLength) {
        chrono::steady_clock::time_point start_t = chrono::steady_clock::now();

        if (!destSeqIsRef)
            *devout << "Destination sequence length: " << destSeqLength << endl;

        if (revComplMatching) {
            if (destSeqIsRef) {
                string querySeq = reverseComplement(destSeq);
                matcher->matchTexts(textMatches, querySeq, destSeqIsRef, revComplMatching, minMatchLength);
            } else {
                reverseComplementInPlace(destSeq);
                matcher->matchTexts(textMatches, destSeq, destSeqIsRef, revComplMatching, minMatchLength);
                upperReverseComplementInPlace(destSeq);
            }
        } else
            matcher->matchTexts(textMatches, destSeq, destSeqIsRef, revComplMatching, minMatchLength);

        *devout << "... found " << textMatches.size() << " exact matches in " << time_millis(start_t) << " msec. " << endl;

        /*        std::sort(textMatches.begin(), textMatches.end(), [](const TextMatch &match1, const TextMatch &match2) -> bool
            { return match1.length > match2.length; });
        cout << "Largest matches:" << endl;
        for (uint32_t i = 0; i < textMatches.size() && i < 10; i++)
            textMatches[i].report(cout);*/

        if (revComplMatching)
            correctDestPositionDueToRevComplMatching();
    }

    using namespace PgTools;

    void SimpleSequenceMatcher::correctDestPositionDueToRevComplMatching() {
        for (TextMatch &match: textMatches)
            match.posDestText = destSeqLength - (match.posDestText + match.length);
    }

    string SimpleSequenceMatcher::getTotalMatchStat(size_t totalMatchLength) {
        return toString(totalMatchLength) + " (" + toString((totalMatchLength * 100.0) / destSeqLength, 1) + "%)";
    }

    char SimpleSequenceMatcher::MATCH_MARK = MBGC_Params::RC_MATCH_MARK;

    void SimpleSequenceMatcher::markAndRemoveExactMatches(
            bool destSeqIsRef, string &destSeq, string &resMapOff, string &resMapLen,
            bool revComplMatching, uint32_t minMatchLength) {
        if (!matcher) {
            resMapOff.clear();
            resMapLen.clear();
            return;
        }

        this->revComplMatching = revComplMatching;
        this->destSeqLength = destSeq.length();

        if (minMatchLength == UINT32_MAX)
            minMatchLength = targetMatchLength;
        exactMatchSequence(destSeq, destSeqIsRef, minMatchLength);

        chrono::steady_clock::time_point post_start_t = chrono::steady_clock::now();
        if (destSeqIsRef)
            resolveMappingCollisionsInTheSameText();

        ostringstream mapOffDest;
        ostringstream mapLenDest;

        PgHelpers::writeUIntByteFrugal(mapLenDest, minMatchLength);

        sort(textMatches.begin(), textMatches.end());
        textMatches.erase(unique(textMatches.begin(), textMatches.end()), textMatches.end());
        *devout << "Unique exact matches: " << textMatches.size() << endl;

        char *destPtr = (char *) destSeq.data();
        size_t pos = 0;
        size_t nPos = 0;
        size_t totalDestOverlap = 0;
        size_t totalMatched = 0;
        bool isSeqLengthStd = srcSeq.length() <= UINT32_MAX;
        for (TextMatch &match: textMatches) {
            if (match.posDestText < pos) {
                size_t overflow = pos - match.posDestText;
                if (overflow >= match.length) {
                    totalDestOverlap += match.length;
                    match.length = 0;
                    continue;
                }
                totalDestOverlap += overflow;
                match.length -= overflow;
                match.posDestText += overflow;
                if (!revComplMatching)
                    match.posSrcText += overflow;
            }
            if (match.length < minMatchLength) {
                totalDestOverlap += match.length;
                continue;
            }
            totalMatched += match.length;
            uint64_t length = match.posDestText - pos;
            memmove(destPtr + nPos, destPtr + pos, length);
            nPos += length;
            destSeq[nPos++] = MATCH_MARK;
            if (isSeqLengthStd)
                PgHelpers::writeValue<uint32_t>(mapOffDest, match.posSrcText);
            else
                PgHelpers::writeValue<uint64_t>(mapOffDest, match.posSrcText);
            PgHelpers::writeUIntByteFrugal(mapLenDest, match.length - minMatchLength);
            pos = match.endPosDestText();
        }
        uint64_t length = destSeq.length() - pos;
        memmove(destPtr + nPos, destPtr + pos, length);
        nPos += length;
        destSeq.resize(nPos);

        textMatches.clear();
        resMapOff = mapOffDest.str();
        mapOffDest.clear();
        resMapLen = mapLenDest.str();
        mapLenDest.clear();

        *devout << "Preparing output time: " << time_millis(post_start_t) << " msec." << endl;
        *devout << "Final size of sequence: " << nPos << " (removed: " <<
             getTotalMatchStat(totalMatched) << "; " << totalDestOverlap << " chars in overlapped dest symbol)" << endl;
    }

    void SimpleSequenceMatcher::resolveMappingCollisionsInTheSameText() {
        for (TextMatch &match: textMatches) {
            if (match.posSrcText > match.posDestText) {
                uint64_t tmp = match.posSrcText;
                match.posSrcText = match.posDestText;
                match.posDestText = tmp;
            }
            if (revComplMatching && match.endPosSrcText() > match.posDestText) {
                uint64_t margin = (match.endPosSrcText() - match.posDestText + 1) / 2;
                match.length -= margin;
                match.posDestText += margin;
            }
        }
    }

    void SimpleSequenceMatcher::rcMatchSequence(string &sequence, string& rcMapOff, string& rcMapLen,
                                                size_t targetMatchLength, uint32_t minMatchLength) {
        chrono::steady_clock::time_point ref_start_t = chrono::steady_clock::now();
        
        PgTools::SimpleSequenceMatcher* matcher = new PgTools::SimpleSequenceMatcher(sequence, targetMatchLength, minMatchLength);
        *devout << "Feeding sequence finished in " << time_millis(ref_start_t) << " msec. " << endl;
        chrono::steady_clock::time_point start_t = chrono::steady_clock::now();

        matcher->markAndRemoveExactMatches(true, sequence, rcMapOff, rcMapLen, true, minMatchLength);
        *devout << "rc-matching sequence finished in " << time_millis(start_t) << " msec. " << endl;
        delete(matcher);
    }

    void SimpleSequenceMatcher::restoreRCMatchedSequence(string &sequence, string& rcMapOff, string& rcMapLen,
                                                    size_t orgSrcLen) {
        istringstream mapOffSrc(rcMapOff);
        istringstream mapLenSrc(rcMapLen);

        bool isSeqLengthStd = orgSrcLen <= UINT32_MAX;
        string destSeq = std::move(sequence);
        sequence.resize(0);
        string tmp;
        uint64_t posDest = 0;
        uint32_t minMatchLength = 0;

        PgHelpers::readUIntByteFrugal(mapLenSrc, minMatchLength);
        uint64_t markPos = 0;
        while ((markPos = destSeq.find(MATCH_MARK, posDest)) != std::string::npos) {
            sequence.append(destSeq, posDest, markPos - posDest);
            posDest = markPos + 1;
            uint64_t matchSrcPos = 0;
            if (isSeqLengthStd) {
                uint32_t tmp;
                PgHelpers::readValue<uint32_t>(mapOffSrc, tmp);
                matchSrcPos = tmp;
            } else
                PgHelpers::readValue<uint64_t>(mapOffSrc, matchSrcPos);
            uint64_t matchLength = 0;
            PgHelpers::readUIntByteFrugal(mapLenSrc, matchLength);
            matchLength += minMatchLength;
            sequence.append(reverseComplement(sequence.substr(matchSrcPos, matchLength)));

        }
        sequence.append(destSeq, posDest, destSeq.length() - posDest);

        *PgHelpers::devout << "Restored sequence of length: " << sequence.length() << endl;
    }
}
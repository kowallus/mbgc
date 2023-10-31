#ifndef PGTOOLS_TEXTMATCHERS_H
#define PGTOOLS_TEXTMATCHERS_H

#include <vector>
#include "../utils/helper.h"

namespace PgTools {

    struct TextMatch {
        uint64_t posSrcText;
        uint64_t length;
        uint64_t posDestText;
        uint64_t nextSrcRegionLoadingPos;

        TextMatch(uint64_t posSrcText, uint64_t length, uint64_t posDestText) :
                posSrcText(posSrcText), length(length), posDestText(posDestText), nextSrcRegionLoadingPos(0)  {}

        TextMatch() : TextMatch(0, 0, 0) {};

        bool operator==(const TextMatch &rhs) const {
            return posSrcText == rhs.posSrcText &&
                   length == rhs.length &&
                   posDestText == rhs.posDestText;
        }

        bool operator!=(const TextMatch &rhs) const {
            return !(rhs == *this);
        }

        bool operator<(const TextMatch &rhs) const {
            if (posDestText < rhs.posDestText)
                return true;
            if (rhs.posDestText < posDestText)
                return false;
            if (posSrcText < rhs.posSrcText)
                return true;
            if (rhs.posSrcText < posSrcText)
                return false;
            return length < rhs.length;
        }

        bool pairedWith(const TextMatch &rhs) const {
            return posSrcText + rhs.posDestText == rhs.posSrcText + posDestText;
        }

        bool pairedWith(const TextMatch &rhs, size_t pairBreakSrcPos) const {
            return posSrcText + rhs.posDestText == rhs.posSrcText + posDestText &&
                    ((posSrcText > pairBreakSrcPos && rhs.posSrcText > pairBreakSrcPos) ||
                            (posSrcText < pairBreakSrcPos && rhs.posSrcText < pairBreakSrcPos));
        }

        const static int8_t NOT_PAIRED_WITH_INDEL = INT8_MAX;

        int8_t pairedWithIndel(const TextMatch &rhs, int8_t maxIndel) const {
            int64_t delta = (int64_t) posSrcText - (int64_t) posDestText;
            int64_t rhsDelta = (int64_t) rhs.posSrcText - (int64_t) rhs.posDestText;
            int64_t indel = rhsDelta - delta;
            if (indel >= -maxIndel && indel <= maxIndel)
                return indel;
            return NOT_PAIRED_WITH_INDEL;
        }

        uint64_t endPosSrcText() const {
            return posSrcText + length;
        }

        uint64_t endPosDestText() const {
            return posDestText + length;
        }

        void shiftStartPos(int64_t shift) {
            this->posSrcText += shift;
            this->posDestText += shift;
            this->length -= shift;
        }

        void report(ostream &out) {
            out << length << ":<" << posSrcText << ", " << endPosSrcText() << ") in " <<
                posDestText << endl;
        }
    };

    class TextMatcher {

    public:
        virtual void matchTexts(vector<TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                                uint32_t minMatchLength) = 0;

        virtual ~TextMatcher() {};

    };

}

#endif //PGTOOLS_TEXTMATCHERS_H

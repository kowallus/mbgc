
#include <deque>
#include "TextMatchers.h"

using namespace PgHelpers;

namespace PgTools {

    void backMatchExpand(const string &text, uint64_t &textPos, const string &pattern, uint64_t &patternPos, uint64_t &length) {
        const char *textGuardPtr = text.data();
        const char *textPtr = text.data() + textPos;
        const char *patternGuardPtr = pattern.data();
        const char *patternPtr = pattern.data() + patternPos;
        uint64_t i = 0;
        while (textPtr-- != textGuardPtr && patternPtr-- != patternGuardPtr && *textPtr == *patternPtr)
            i++;
        textPos -= i;
        patternPos -= i;
        length += i;
    }

    void forwardMatchExpand(const string &text, uint64_t textPos, const string &pattern, uint64_t patternPos, uint64_t &length) {
        const char *textGuardPtr = text.data() + text.length();
        const char *textPtr = text.data() + textPos + length;
        const char *patternGuardPtr = pattern.data() + pattern.length();
        const char *patternPtr = pattern.data() + patternPos + length;
        while (*textPtr == *patternPtr && textPtr++ != textGuardPtr && patternPtr++ != patternGuardPtr)
            length++;
    }

}

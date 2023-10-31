#include "SymbolsPackingFacility.h"

#include <cassert>

namespace PgIndex {

    static const vector<char> binaryCodes { 0, 1 };
    static const vector<char> trenaryCodes { 0, 1, 2 };
    static const vector<char> quaternaryCodes { 0, 1, 2, 3 };

    static const vector<char> acgtSymbols { 'A', 'C', 'G', 'T' };
    static const vector<char> acgtnSymbols { 'A', 'C', 'G', 'N', 'T' };

    SymbolsPackingFacility SymbolsPackingFacility::ACGTPacker(acgtSymbols, true);
    SymbolsPackingFacility SymbolsPackingFacility::ACGTNPacker(acgtnSymbols, true);
    
    SymbolsPackingFacility::SymbolsPackingFacility(const vector<char> symbolsList, bool isGloballyManaged):
         symbolsCount(symbolsList.size()),
         globallyManaged(isGloballyManaged),
         symbolsPerElement(SymbolsPackingFacility::maxSymbolsPerElement(symbolsCount)) {
        size_t combinationCount = powuint(symbolsCount, symbolsPerElement);
        maxValue = combinationCount - 1;
        if (maxValue > (int) (uint8_t) - 1)
            cerr << "ERROR in symbols packaging: max value for type: " << (int) (uint8_t) - 1 << " while max " << " \n";

        reverse = new char*[combinationCount];
        reverseFlat = new char[combinationCount * symbolsPerElement];

        clear = new uint8_t*[combinationCount];
        clearFlat = new uint8_t[combinationCount * symbolsPerElement];

        std::copy(symbolsList.begin(), symbolsList.end(), std::begin(this->symbolsList));
        memset(symbolOrder, -1, UINT8_MAX);
        for (int i = 0; i < symbolsCount; i++)
            symbolOrder[(unsigned char) symbolsList[(unsigned char) i]] = i;

        buildReversePackAndClearIndexes();
    }

    SymbolsPackingFacility::~SymbolsPackingFacility() {
        delete[]reverse;
        delete[]reverseFlat;
        delete[]clear;
        delete[]clearFlat;
        delete[]packLUT0;
        delete[]packLUT1;
    }

    void SymbolsPackingFacility::buildReversePackAndClearIndexes() {
        char* rPtr = reverseFlat;
        uint8_t* cPtr = clearFlat;
        for (size_t i = 0; i <= maxValue; i++) {
            reverse[i] = rPtr;
            rPtr += symbolsPerElement;
            clear[i] = cPtr;
            cPtr += symbolsPerElement;
        }

        uint8_t* currentClear = new uint8_t[symbolsPerElement]();
        uint8_t* sequence = new uint8_t[symbolsPerElement]();
        for (size_t i = 0; i <= maxValue; i++) {
            for (uint8_t j = 0; j < symbolsPerElement; j++) {
                reverse[i][j] = symbolsList[sequence[j]];
                clear[i][j] = currentClear[j];
            }

            uint8_t j = symbolsPerElement - 1;
            while ((++sequence[j] == symbolsCount) && ((int) j > 0)) {
                sequence[j] = 0;
                currentClear[j] = i + 1;
                j--;
            }
        }
        delete[]currentClear;
        delete[]sequence;

        packLUT0 = new uint8_t[PACK_LUT_SIZE]();
        packLUT1 = new uint8_t[PACK_LUT_SIZE]();
        symbolsPerLUT1 = symbolsPerElement - SYMBOLS_PER_LUT_0;
        for (size_t i = 0; i <= maxValue; i++) {
            uint16_t temp = 0;
            if (reverse[i][SYMBOLS_PER_LUT_0] == symbolsList[0]
                && (symbolsPerLUT1 == 1 || reverse[i][SYMBOLS_PER_LUT_0 + 1] == symbolsList[0])) {
                memcpy(&temp, reverse[i], SYMBOLS_PER_LUT_0);
                packLUT0[temp & PACK_MASK] = (uint8_t) i;
            }
            if (reverse[i][0] == symbolsList[0] && reverse[i][1] == symbolsList[0]) {
                memcpy(&temp, reverse[i] + SYMBOLS_PER_LUT_0, symbolsPerLUT1);
                packLUT1[temp & PACK_MASK] = (uint8_t) i;
            }
        }
    }

    uint8_t SymbolsPackingFacility::clearSuffix(const uint8_t value, uint8_t prefixLength) {
        return clear[value][prefixLength];
    }

    uint8_t SymbolsPackingFacility::getMaxValue() {
        return maxValue;
    }

    bool SymbolsPackingFacility::isCompatible(uint8_t symbolsPerElement, uint8_t symbolsCount) {
        return powuint(symbolsCount, symbolsPerElement) - 1 <= (uint8_t) - 1;
    }

    uint8_t SymbolsPackingFacility::maxSymbolsPerElement(uint8_t symbolsCount) {
        for (int i = 0; i < UINT8_MAX; i++)
            if (!SymbolsPackingFacility::isCompatible(i + 1, symbolsCount))
                return i;
        return UINT8_MAX;
    }

    uint8_t SymbolsPackingFacility::packPrefixSymbols(const char* symbols, const size_t length) {
        uint8_t value = 0;
        for (uint8_t j = 0; j < length; j++) {
            validateSymbol((uint8_t) symbols[j]);
            value = value * symbolsCount + symbolOrder[(uint8_t) symbols[j]];
        }

        return value;
    }

    size_t SymbolsPackingFacility::packSequence(const char* source, const size_t length, uint8_t* dest) {
        size_t i = 0;

        const char* guard = source + length - symbolsPerElement;

        while (source <= guard) {
            dest[i++] = packSymbols(source);
            source += symbolsPerElement;
        }

        size_t left = guard + symbolsPerElement - source;
        if (left > 0)
            dest[i++] = packSuffixSymbols(source, left);

        return i;
    }

    string SymbolsPackingFacility::packSequence(const char *source, const size_t length) {
        size_t packedLength = (length + symbolsPerElement - 1) / symbolsPerElement * sizeof(uint8_t);
        string tmp;
        tmp.resize(packedLength);
        this->packSequence(source, length, (uint8_t*) tmp.data());
        return tmp;
    }

    uint8_t SymbolsPackingFacility::packSuffixSymbols(const char* symbols, const size_t length) {
        uint8_t value = 0;
        for (uint8_t j = 0; j < symbolsPerElement; j++) {
            value *= symbolsCount;
            if (j < length) {
                validateSymbol((uint8_t) symbols[j]);
                value += symbolOrder[(uint8_t) symbols[j]];
            }
        }
        return value;
    }

    uint8_t SymbolsPackingFacility::packSymbols(const char* symbols) {
        uint16_t temp0 = 0, temp1 = 0;
        memcpy(&temp0, symbols, SYMBOLS_PER_LUT_0);
        memcpy(&temp1, symbols + SYMBOLS_PER_LUT_0, symbolsPerLUT1);
        return packLUT0[temp0 & PACK_MASK] + packLUT1[temp1 & PACK_MASK];
    }

    const string SymbolsPackingFacility::reverseSequence(const uint8_t* sequence, const size_t pos, const size_t length) {
        string res;
        size_t i = divideBySmallInteger(pos, symbolsPerElement);
        size_t reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);
        res.append(reverse[sequence[i++]], reminder, this->symbolsPerElement - reminder);
        while (res.size() < length)
            res.append(reverse[sequence[i++]], symbolsPerElement);
        res.resize(length);
        return res;
    }

    void SymbolsPackingFacility::reverseSequence(const uint8_t* sequence, const size_t pos, const size_t length, string &res) {
        res.clear();
        size_t i = divideBySmallInteger(pos, symbolsPerElement);
        size_t reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);
        uint8_t value = sequence[i++];
        for (uint8_t j = reminder; j < symbolsPerElement; j++) {
            res.push_back(reverse[value][j]);
            if (res.size() < length)
                return;
        }
        res.append(reverse[sequence[i++]], reminder, this->symbolsPerElement - reminder);
        while (res.size() < length - symbolsPerElement)
            res.append(reverse[sequence[i++]], symbolsPerElement);
        value = sequence[i];
        for (uint8_t j = 0; res.size() < length; j++)
            res.push_back(reverse[value][j]);
    }

    void SymbolsPackingFacility::reverseSequence(const uint8_t* sequence, const size_t pos, const size_t length, char* destPtr) {
        size_t i = divideBySmallInteger(pos, symbolsPerElement);
        size_t reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);

        char* ptr = destPtr;
        const char* endPtr = destPtr + length;
        uint8_t value = sequence[i++];
        for (uint8_t j = reminder; j < symbolsPerElement; j++) {
            *ptr++ = reverse[value][j];
            if (ptr == endPtr)
                return;
        }
        while(true) {
            value = sequence[i++];
            for (uint8_t j = 0; j < symbolsPerElement; j++) {
                *ptr++ = reverse[value][j];
                if (ptr == endPtr)
                    return;
            }
        }
   }

    char SymbolsPackingFacility::reverseValue(uint8_t value, uint8_t position) {
        return reverse[value][position];
    }

    char SymbolsPackingFacility::reverseValue(uint8_t* sequence, size_t pos) {
        size_t i = divideBySmallInteger(pos, symbolsPerElement);
        size_t reminder = moduloBySmallInteger<size_t>(pos, this->symbolsPerElement, i);

        uint8_t value = sequence[i];
        return reverse[value][reminder];
    }

    string SymbolsPackingFacility::reverseValue(uint8_t value) {
        string res;
        res.resize(symbolsPerElement);
        for (uint8_t j = 0; j < symbolsPerElement; j++)
            res[j] = reverse[value][j];
        return res;
    }

    int SymbolsPackingFacility::compareSequences(uint8_t* lSeq, uint8_t* rSeq, const size_t length) {
        size_t i = length;
        while (i >= symbolsPerElement) {
            int cmp = (int) *lSeq++ - *rSeq++;
            if (cmp)
                return cmp;
            i -= symbolsPerElement;
        } 
        
        int j = 0;
        while (i--) {
            int cmp = (int) reverse[*lSeq][j] - reverse[*rSeq][j];
            if (cmp)
                return cmp;
            j++;
        }
        
        return 0;
    }

    int SymbolsPackingFacility::compareSequences(uint8_t* lSeq, uint8_t* rSeq, size_t pos, size_t length) {
        size_t i = divideBySmallInteger(pos, symbolsPerElement);
        size_t reminder = moduloBySmallInteger(pos, this->symbolsPerElement, i);

        lSeq += i;
        rSeq += i;
        for (uint8_t j = reminder; j < symbolsPerElement; j++) {
            int cmp = (int) reverse[*lSeq][j] - reverse[*rSeq][j];
            if (cmp)
                return cmp;
            if (--length == 0)
                return 0;
        }
        
        return compareSequences(lSeq + 1, rSeq + 1, length);
    }

    int SymbolsPackingFacility::compareSuffixWithPrefix(uint8_t* sufSeq, uint8_t* preSeq, size_t sufPos, size_t length) {
        size_t i = divideBySmallInteger(sufPos, symbolsPerElement);
        size_t reminder = moduloBySmallInteger(sufPos, this->symbolsPerElement, i);
        
        sufSeq += i;
        if (reminder == 0)
            return compareSequences(sufSeq, preSeq, length);
        
        uint8_t sufIdx = reminder;
        uint8_t preIdx = 0;
        while (true) {
            int cmp = (int) reverse[*sufSeq][sufIdx] - reverse[*preSeq][preIdx];
            if (cmp)
                return cmp;
            if (--length == 0)
                return 0;
            if (++preIdx == symbolsPerElement) {
                preIdx = 0; preSeq++;
            } 
            if (++sufIdx == symbolsPerElement) {
                sufIdx = 0; sufSeq++;
            }
        }
    }

    int SymbolsPackingFacility::compareSequenceWithUnpacked(uint8_t *seq,
                                                                          const char *pattern,
                                                                          size_t length) {
        size_t i = length;
        while (i >= symbolsPerElement) {
            uint8_t pSeq = packSymbols(pattern);
            int cmp = (int) *seq++ - pSeq;
            if (cmp)
                return cmp;
            pattern += symbolsPerElement;
            i -= symbolsPerElement;
        }

        int j = 0;
        while (i--) {
            int cmp = (int) reverse[*seq][j] - *pattern++;
            if (cmp)
                return cmp;
            j++;
        }

        return 0;
    }

    uint8_t SymbolsPackingFacility::countSequenceMismatchesVsUnpacked(uint8_t *seq,
                                                                                    const char *pattern,
                                                                                    size_t length,
                                                                                    uint8_t maxMismatches) {
        uint8_t res = 0;
        size_t i = length;
        while (i >= symbolsPerElement) {
            uint8_t pSeq = packSymbols(pattern);
            if (*seq != pSeq) {
                for(uint8_t j = 0; j < symbolsPerElement; j++) {
                    if (reverse[*seq][j] != *pattern++)
                        if (res++ >= maxMismatches)
                            return UINT8_MAX;
                }
            } else
                pattern += symbolsPerElement;
            seq++;
            i -= symbolsPerElement;
        }

        int j = 0;
        while (i--) {
            if (reverse[*seq][j] != *pattern++) {
                if (res++ >= maxMismatches)
                    return UINT8_MAX;
            }
            j++;
        }

        return res;
    }

    void SymbolsPackingFacility::validateSymbol(uint8_t symbol) {
        if (symbolOrder[symbol] == -1) {
            fprintf(stdout, "Unexpected symbol '%c'. Packer supports: %s.\n", symbol, symbolsList);
            exit(EXIT_FAILURE);
        }
    }
}
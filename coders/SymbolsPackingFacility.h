#ifndef SYMBOLSPACKINGFACILITY_H
#define	SYMBOLSPACKINGFACILITY_H

#include "../utils/helper.h"
#include <vector>

using namespace PgHelpers;

namespace PgIndex {

    class SymbolsPackingFacility {
    private:

        const static uint16_t PACK_LUT_SIZE = 1 << 11;
        const static uint16_t PACK_MASK = PACK_LUT_SIZE - 1;

        size_t maxValue;
        const uint8_t symbolsCount;
        const uint8_t symbolsPerElement;
        
        char symbolsList[UINT8_MAX] = {};
        int symbolOrder[UINT8_MAX] = {};

        uint8_t* packLUT0;
        uint8_t* packLUT1;
        const uint8_t SYMBOLS_PER_LUT_0 = 2;
        uint8_t symbolsPerLUT1;

        // for a given value reverse[value][pos] returns character at the given position
        char** reverse;
        char* reverseFlat;
        
        // for a given value clear[value][prefixLenght] returns value of first prefixLength symbols followed by 0 symbols.
        uint8_t** clear;
        uint8_t* clearFlat;

        const bool globallyManaged = false;

        void buildReversePackAndClearIndexes();

        inline void validateSymbol(uint8_t symbol);

    public:
        SymbolsPackingFacility(const vector<char> symbolsList, bool isGloballyManaged = false);
        
        ~SymbolsPackingFacility();

        bool isGloballyManaged() { return globallyManaged; };

        uint8_t getMaxValue();
        
        // value should be not greater then maxValue
        char reverseValue(uint8_t value, uint8_t position);
        char reverseValue(uint8_t* sequence, size_t position);

        // value should be not greater then maxValue
        string reverseValue(uint8_t value);
        
        const string reverseSequence(const uint8_t* sequence, const size_t pos, const size_t length);
        void reverseSequence(const uint8_t* sequence, const size_t pos, const size_t length, string& res);
        void reverseSequence(const uint8_t* sequence, const size_t pos, const size_t length, char* destPtr);
      
        // sequence should consist of at least symbolsPerElement symbols; result should be not greater then maxValue
        inline uint8_t packSymbols(const char* symbols);
        
        // result should be not greater then maxValue
        uint8_t packSuffixSymbols(const char *symbols, const size_t length);
        
        // result should be not greater then maxValue
        uint8_t packPrefixSymbols(const char* symbols, const size_t length);
        
        size_t packSequence(const char* source, const size_t length, uint8_t* dest);
        string packSequence(const char* source, const size_t length);

        uint8_t clearSuffix(const uint8_t value, uint8_t prefixLength);
        
        static bool isCompatible(uint8_t symbolsPerElement, uint8_t symbolsCount);
        
        static uint8_t maxSymbolsPerElement(uint8_t symbolsCount);
      
        int compareSequences(uint8_t* lSeq, uint8_t* rSeq, const size_t length);
        int compareSequences(uint8_t* lSeq, uint8_t* rSeq, size_t pos, size_t length);
        int compareSuffixWithPrefix(uint8_t* sufSeq, uint8_t* preSeq, size_t sufPos, size_t length);

        int compareSequenceWithUnpacked(uint8_t* seq, const char *pattern, size_t length);

        uint8_t countSequenceMismatchesVsUnpacked(uint8_t* seq, const char *pattern,  size_t length,
                                          uint8_t maxMismatches);

        static SymbolsPackingFacility ACGTPacker, ACGTNPacker;
    };
    
}

#endif	/* SYMBOLSPACKINGFACILITY_H */


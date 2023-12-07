#ifndef HELPER_H_INCLUDED
#define HELPER_H_INCLUDED

#include <ctime>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>
#include <set>
#include <algorithm>
#include <climits>
#include <cmath>
#include <sstream>

using namespace std;

class NullBuffer : public std::streambuf
{
public:
    int overflow(int c) { return c; }
};

class WrapperStreamBuf : public std::streambuf
{
public:
    WrapperStreamBuf() { setg(nullptr, nullptr, nullptr); }
    WrapperStreamBuf(char* s, std::size_t n) { setg(s, s, s + n); }
    void wrapData(char* s, std::size_t n) { setg(s, s, s + n); }
    size_t getPosition() { return gptr() - eback(); }
    void setPosition(size_t pos) { setg(eback(), eback() + pos, egptr()); }
    bool isEnd() { return gptr() == egptr(); };
};

class WrapperStrStream : public std::istringstream
{
private:
    WrapperStreamBuf buf;
public:
    WrapperStrStream() : std::istringstream("") { }
    WrapperStrStream(string& s) : std::istringstream("") { str(s); }
    WrapperStrStream(const WrapperStrStream& wss) : std::istringstream() { set_rdbuf(&(buf = wss.buf)); }
    WrapperStrStream(WrapperStrStream&& wss) : std::istringstream() { set_rdbuf(&(buf = wss.buf)); }
    WrapperStrStream& operator=(WrapperStrStream &wss) { set_rdbuf(&(buf = wss.buf)); return *this; }
    void str(string& str) { wrapData((char*) str.data(), str.size()); }
    size_t getPosition() { return buf.getPosition(); }
    void setPosition(size_t pos) { buf.setPosition(pos); }
    bool isEnd() { return buf.isEnd(); }
    void skipOffset(uint64_t offset) { buf.setPosition(buf.getPosition() + offset); }

    void wrapData(char* s, size_t n) {
        buf.wrapData(s, n);
        set_rdbuf(&buf);
    }
};


extern std::ostream null_stream;

namespace PgHelpers {

    extern int numberOfThreads;

    // bioinformatical routines

    char upperReverseComplement(char symbol);
    char reverseComplement(char symbol);
    string upperReverseComplement(string kmer);
    string reverseComplement(string kmer);
    void upperReverseComplement(const char* start, const std::size_t N, char* target);
    void reverseComplement(const char* start, const std::size_t N, char* target);
    void upperReverseComplementInPlace(char* start, const std::size_t N);
    void reverseComplementInPlace(char* start, const std::size_t N);
    void upperReverseComplementInPlace(string &kmer);
    void reverseComplementInPlace(string &kmer);
    void upperSequence(char* start, const std::size_t N);
    double qualityScore2approxCorrectProb(string quality);
    double qualityScore2correctProb(string quality);

    void revertToMBGC121();

    template<typename uint_read_len>
    void convertMisRevOffsets2Offsets(uint_read_len *mismatchOffsets, uint8_t mismatchesCount, uint_read_len readLength) {
        for(uint8_t i = 0; i < mismatchesCount / 2; i++) {
            uint_read_len tmp = mismatchOffsets[mismatchesCount - i - 1];
            mismatchOffsets[mismatchesCount - i - 1] = mismatchOffsets[i];
            mismatchOffsets[i] = tmp;
        }
        uint_read_len pos = readLength;
        for(uint8_t i = mismatchesCount; i-- > 0;) {
            pos -= mismatchOffsets[i] + 1;
            mismatchOffsets[i] = pos;
        }
    }

    // time routines

    clock_t clock_checkpoint();
    unsigned long long int clock_millis();
    unsigned long long int clock_millis(clock_t checkpoint);

    chrono::steady_clock::time_point time_checkpoint();
    unsigned long long int time_millis();
    unsigned long long int time_millis(chrono::steady_clock::time_point checkpoint);


    // string conversion routines

    string toString(unsigned long long value);
    string toMB(unsigned long long value, unsigned char decimalPlaces);
    string toString(long double value, unsigned char decimalPlaces);


    // mathematical routines

    unsigned long long int powuint(unsigned long long int base, int exp);

    template<typename uint>
    inline uint divideBySmallInteger(const uint dividend, const unsigned char divisor) {
        switch (divisor) {
            case 1: return dividend / 1;
            case 2: return dividend / 2;
            case 3: return dividend / 3;
            case 4: return dividend / 4;
            case 5: return dividend / 5;
            case 6: return dividend / 6;
            case 7: return dividend / 7;
            case 8: return dividend / 8;
        };
        cout << "Unsupported denominator " << divisor << "!\n";
        return 0;
    }

    template<typename uint>
    inline uint multiplyBySmallInteger(const uint value, const unsigned char smallInteger) {
        switch (smallInteger) {
            case 0: return 0;
            case 1: return value * 1;
            case 2: return value * 2;
            case 3: return value * 3;
            case 4: return value * 4;
            case 5: return value * 5;
            case 6: return value * 6;
            case 7: return value * 7;
            case 8: return value * 8;
        };
        cout << "Unsupported multiplied " << smallInteger << "!\n";
        return 0;
    }

    template<typename uint>
    inline uint moduloBySmallInteger(const uint dividend, const unsigned char divisor, const uint resultOfDivision) {
        return dividend - resultOfDivision * divisor;
    }
    
    template<typename uint>
    inline uint moduloBySmallInteger(const uint dividend, const unsigned char divisor) {
        return moduloBySmallInteger(dividend, divisor, divideBySmallInteger(dividend, divisor));
    }

    template<typename t_val>
    string transpose(string matrix, uint64_t rows, uint64_t cols) {
        string tRes;
        tRes.resize(matrix.size());
        t_val* nMatrix = (t_val*) matrix.data();
        t_val* tMatrix = (t_val*) tRes.data();
        for(uint64_t i = 0; i < rows; i++) {
            for(uint8_t m = 0; m < cols; m++)
                tMatrix[m * rows + i] = nMatrix[i * cols + m];
        }
        return tRes;
    }

    // string comparison routines

    int readsSufPreCmp(const char* suffixPart, const char* prefixRead);
    
    int strcmplcp(const char* lStrPtr, const char* rStrPtr, int length);

    // input output routines

    extern std::ostream *appout;
    extern std::ostream *devout;
    extern std::ostream *logout;

    void openLogFile(char* filename);
    void closeLogFile();

    void createFolders(string pathToFile);
    void normalizePath(string &fileName, int &finalPos, int &finalLen, bool ignorePath);

    void* readArray(std::istream&, size_t arraySizeInBytes);
    void readArray(std::istream&, void* destArray, size_t arraySizeInBytes);
    void writeArray(std::ostream&, void* srcArray, size_t arraySize, bool verbose = false);

    void* readWholeArray(std::istream&, size_t& arraySizeInBytes);
    void* readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes);
    void writeArrayToFile(string destFile, void* srcArray, size_t arraySize);
    void writeStringToFile(string destFile, const string &src);

    extern bool plainTextWriteMode;
    const static string TEXT_MODE_ID = "TXT";
    const static string BINARY_MODE_ID = "BIN";

    void writeReadMode(std::ostream &dest, bool plainTextWriteMode);
    bool confirmTextReadMode(std::istream &src);

    template<typename t_val>
    inline void writeValue(std::ostream &dest, const t_val value) {
        dest.write((char *) &value, sizeof(t_val));
    }

    template<typename t_val>
    inline void readValue(std::istream &src, t_val& value) {
        src.read((char *) &value, sizeof(t_val));
    }

    template<typename t_val>
    void readValue(unsigned char* &src, t_val& value) {
        memcpy((void*) &value, src, sizeof(value));
        src += sizeof(value);
    }

    void writeUIntByteFrugal(std::ostream &dest, uint64_t value);
    void writeUIntWordFrugal(std::ostream &dest, uint64_t value);
    void writeUInt64Frugal(std::ostream &dest, uint64_t value);

    template<typename t_val>
    void readUIntByteFrugal(std::istream &src, t_val& value) {
        value = 0;
        uint8_t yByte = 0;
        t_val base = 1;
        do {
            src.read((char *) &yByte, sizeof(uint8_t));
            value += base * (yByte % 128);
            base *= 128;
        } while (yByte >= 128);
    }

    template<typename t_val>
    void readUIntWordFrugal(std::istream &src, t_val& value) {
        value = 0;
        uint16_t yWord = 0;
        t_val base = 1;
        do {
            src.read((char *) &yWord, sizeof(uint16_t));
            value += base * (yWord % 32768);
            base *= 32768;
        } while (yWord >= 32768);
    }

    template<typename t_val>
    void readUInt64Frugal(std::istream &src, t_val& value) {
        uint16_t yWord = 0;
        src.read((char *) &yWord, sizeof(uint16_t));
        if (yWord < UINT16_MAX)
            value = yWord;
        else {
            uint32_t y32 = 0;
            src.read((char *) &y32, sizeof(uint32_t));
            if (y32 < UINT32_MAX)
                value = y32;
            else {
                uint64_t y64 = 0;
                src.read((char *) &y64, sizeof(uint64_t));
                value = y64;
            }
        }
    }

    extern bool bytePerReadLengthMode;

    void readReadLengthValue(std::istream &src, uint16_t& value);

    void writeReadLengthValue(std::ostream &dest, const uint16_t value);

    class BufferedFileIStream : public istream {

    private:

        struct membuf : std::streambuf
        {
            membuf(char* begin, char* end) {
                this->setg(begin, begin, end);
            }
        };

        membuf* sbuf = 0;
        char* readsArray = 0;

        BufferedFileIStream(membuf* sbuf, char* readsArray) : istream(sbuf) {
            this->sbuf = sbuf;
            this->readsArray = readsArray;
        }

    public:

        static BufferedFileIStream* getIStream(string filename) {
            size_t readsArraySize;
            char* readsArray = (char*) PgHelpers::readWholeArrayFromFile(filename, readsArraySize);

            membuf* sbuf = new membuf(readsArray, readsArray + readsArraySize);

            BufferedFileIStream* source = new BufferedFileIStream(sbuf, readsArray);
            if (!*source)
                std::cout << "Whoops";

            return source;
        }

        virtual ~BufferedFileIStream() {
            delete(sbuf);
            delete[] readsArray;
        }
    };

    // memory management routines

    unsigned long long getTotalSystemMemory();

    template<typename t_arr>
    size_t safeNewArrayAlloc(t_arr *&arr, size_t allocSize, bool zerofill, uint8_t totalLimitPercent = UINT8_MAX) {
        size_t totalRAM = getTotalSystemMemory();
        if (totalLimitPercent != UINT8_MAX && allocSize > totalRAM * totalLimitPercent / 100)
            allocSize = totalRAM * totalLimitPercent / 100;
        bool success = false;
        do {
            try {
                arr = zerofill ? new t_arr[allocSize]() : new t_arr[allocSize];
                success = true;
            } catch (const std::bad_alloc& e) {
                allocSize /= 2;
            }
        } while (!success);
        return allocSize;
    }

}

#if defined(__arm__) || defined(__aarch64__) || defined(__ARM_ARCH)

void A_memcpy(void *dest, const void *src, size_t n);

int A_memcmp(const void *s1, const void *s2, size_t n);

#endif

#endif // HELPER_H_INCLUDED

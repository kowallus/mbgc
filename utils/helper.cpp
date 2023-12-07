#include "helper.h"

#include "byteswap.h"
#include <sys/stat.h>
#include <regex>

#ifdef __MINGW32__
#include <sysinfoapi.h>
#else
#include <unistd.h>
#endif

int PgHelpers::numberOfThreads = 8;

std::ostream *PgHelpers::appout = &std::cout;
std::ostream *PgHelpers::devout = &std::cout;

NullBuffer null_buffer;
std::ostream null_stream(&null_buffer);

std::ostream *PgHelpers::logout = &null_stream;
std::fstream *logfileout = NULL;

void PgHelpers::openLogFile(char* filename) {
    logfileout = new fstream(filename, ios_base::in | ios_base::app);
    if (!*logfileout) {
        fprintf(stderr, "Cannot open log file: %s\n", filename);
        exit(EXIT_FAILURE);
    }
    logout = logfileout;
}

void PgHelpers::closeLogFile() {
    if (logfileout) {
        delete(logfileout);
        logfileout = NULL;
        logout = &null_stream;
    }
}

// TIME

clock_t checkpoint;

clock_t PgHelpers::clock_checkpoint() {
    checkpoint = clock();
    return checkpoint;
//    cout << "Clock reset!\n";
}

unsigned long long int PgHelpers::clock_millis(clock_t checkpoint) {
    return (clock() - checkpoint) * (unsigned long long int) 1000 / CLOCKS_PER_SEC;
}

unsigned long long int PgHelpers::clock_millis() {
    return clock_millis(checkpoint);
}


chrono::steady_clock::time_point chronocheckpoint;

chrono::steady_clock::time_point PgHelpers::time_checkpoint() {
    chronocheckpoint = chrono::steady_clock::now();
    return chronocheckpoint;
}

unsigned long long int PgHelpers::time_millis(chrono::steady_clock::time_point checkpoint) {
    chrono::nanoseconds time_span = chrono::duration_cast<chrono::nanoseconds>(chrono::steady_clock::now() - checkpoint);
    return (double)time_span.count() / 1000000.0;
}

unsigned long long int PgHelpers::time_millis() {
    return time_millis(chronocheckpoint);
}

void PgHelpers::createFolders(string pathToFile) {
    size_t pos = 0;
    while ((pos = pathToFile.find('/', pos)) != std::string::npos) {
#ifdef __MINGW32__
        mkdir(pathToFile.substr(0, pos++).c_str());
#else
        mkdir(pathToFile.substr(0, pos++).c_str(), 0777);
#endif
    }
}

void PgHelpers::normalizePath(string &fileName, int &finalPos, int &finalLen, bool ignorePath) {
    if(!fileName.empty() && *fileName.rbegin() == '\r') fileName.resize(fileName.length() - 1);
    fileName = regex_replace(fileName, regex("\\\\"), "/");
    while (fileName.find("//") != string::npos)
        fileName = regex_replace(fileName, regex("//"), "/");
    if(!fileName.empty() && fileName.substr(0, 2) == "./")
        fileName = fileName.substr(2);
    bool removePath = ignorePath || fileName.find(':') != string::npos ||
                      fileName.find("/./") != string::npos || fileName.find("/../") != string::npos ||
                      fileName.substr(0, 3) == "../" || fileName[0] == '/';
    finalPos = removePath ? fileName.find_last_of("/") + 1 : 0;
    finalLen = fileName.length() - finalPos;
    if (fileName.size() > 3 && fileName.substr(fileName.size() - 3, 3) == ".gz")
        finalLen -= 3;
}

const size_t chunkSize = 10000000;

void* PgHelpers::readArray(std::istream& in, size_t arraySizeInBytes) {

    if (in) {
        char* destArray = new char[arraySizeInBytes];
        readArray(in, destArray, arraySizeInBytes);
        return (void*) destArray;
    } else
        throw(errno);
}

void PgHelpers::readArray(std::istream& in, void* destArray, size_t arraySizeInBytes) {

    if (in) {
        size_t length = arraySizeInBytes;
        char* ptr = (char*) destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);
    } else
        throw(errno);
}

void PgHelpers::writeArray(std::ostream& out, void* srcArray, size_t arraySize, bool verbose) {
    size_t bytesLeft = arraySize;
    if (out) {
        while (bytesLeft > chunkSize) {
            out.write((char*) srcArray, chunkSize);
            srcArray = (void*) (((char*) srcArray) + chunkSize);
            bytesLeft -= chunkSize;
        }
        out.write((char*) srcArray, bytesLeft);

    } else
        throw(errno);
    if (verbose)
        cout << "Written " << arraySize << " bytes\n";
}

void* PgHelpers::readWholeArray(std::istream& in, size_t& arraySizeInBytes) {

    if (in) {
        in.seekg(0, in.end);
        size_t length = in.tellg();
        char* destArray = new char[length];
        in.seekg(0, in.beg);
        char* ptr = destArray;
        size_t bytesLeft = length;
        while (bytesLeft > chunkSize) {
            in.read(ptr, chunkSize);
            ptr = ptr + chunkSize;
            bytesLeft -= chunkSize;
        }
        in.read(ptr, bytesLeft);

        arraySizeInBytes = length;
        return (void*) destArray;
    } else
        throw(errno);
}

void* PgHelpers::readWholeArrayFromFile(string srcFile, size_t& arraySizeInBytes) {
    time_checkpoint();

    std::ifstream in(srcFile.c_str(), std::ifstream::binary);

    void* destArray = PgHelpers::readWholeArray(in, arraySizeInBytes);

    cout << "Read " << arraySizeInBytes << " bytes from " << srcFile << " in " << time_millis() << " msec \n";

    return destArray;
}

void PgHelpers::writeArrayToFile(string destFile, void* srcArray, size_t arraySize) {
    time_checkpoint();

    std::ofstream out(destFile.c_str(), std::ios::out | std::ios::binary);

    PgHelpers::writeArray(out, srcArray, arraySize);

    cout << "Write " << arraySize << " bytes to " << destFile << " in " << time_millis() << " msec \n";
}

void PgHelpers::writeStringToFile(string destFile, const string &src) {
    writeArrayToFile(destFile, (void*) src.data(), src.length());
}

void PgHelpers::writeReadMode(std::ostream &dest, bool plainTextWriteMode) {
    dest << (plainTextWriteMode?TEXT_MODE_ID:BINARY_MODE_ID) << "\n";
}

bool PgHelpers::confirmTextReadMode(std::istream &src) {
    string readMode;
    src >> readMode;
    if (readMode != TEXT_MODE_ID && readMode != BINARY_MODE_ID) {
        fprintf(stderr, "Expected READ MODE id (not: %s)\n", readMode.c_str());
        exit(EXIT_FAILURE);
    }
    char check = src.get();
    return readMode == TEXT_MODE_ID;
}

void PgHelpers::writeUIntByteFrugal(std::ostream &dest, uint64_t value) {
    while (value >= 128) {
        uint8_t yByte = 128 + (value % 128);
        dest.write((char *) &yByte, sizeof(uint8_t));
        value = value / 128;
    }
    uint8_t yByte = value;
    dest.write((char *) &yByte, sizeof(uint8_t));
}

void PgHelpers::writeUIntWordFrugal(std::ostream &dest, uint64_t value) {
    while (value >= 32768) {
        uint16_t yWord = 32768 + (value % 32768);
        dest.write((char *) &yWord, sizeof(uint16_t));
        value = value / 32768;
    }
    uint16_t yWord = value;
    dest.write((char *) &yWord, sizeof(uint16_t));
}

void PgHelpers::writeUInt64Frugal(std::ostream &dest, uint64_t value) {
    uint16_t yWord = value < UINT16_MAX ? value : UINT16_MAX;
    dest.write((char *) &yWord, sizeof(uint16_t));
    if (value >= UINT16_MAX) {
        uint32_t y32 = value < UINT32_MAX ? value : UINT32_MAX;
        dest.write((char *) &y32, sizeof(uint32_t));
        if (value >= UINT32_MAX)
            dest.write((char *) &value, sizeof(value));
    }
}

bool PgHelpers::bytePerReadLengthMode = false;

void PgHelpers::readReadLengthValue(std::istream &src, uint16_t &value) {
    if (bytePerReadLengthMode)
        readValue<uint8_t>(src, (uint8_t&) value);
    else
        readValue<uint16_t>(src, value);
}

void PgHelpers::writeReadLengthValue(std::ostream &dest, const uint16_t value) {
    if (bytePerReadLengthMode)
        writeValue<uint8_t>(dest, (uint8_t) value);
    else
        writeValue<uint16_t>(dest, value);
}

string PgHelpers::toString(unsigned long long value) {
        std::ostringstream oss;
        oss << value;
        return oss.str();
};

string PgHelpers::toMB(unsigned long long value, unsigned char decimalPlaces) {
    std::ostringstream oss;
    int power = 1000000;
    oss << value / power;
    if (decimalPlaces > 0) {
        oss << ".";
        for(int i = 0; i < decimalPlaces; i++) 
            oss << ((value / (power /= 10)) % 10);
        
    }
    return oss.str();
}

string PgHelpers::toString(long double value, unsigned char decimalPlaces) {
    std::ostringstream oss;
    oss << (long long int) value;
    if (decimalPlaces > 0) {
        oss << ".";
        int power = 1000000;        
        unsigned long long decimals = (value - (long long int) value) * power;
        for(int i = 0; i < decimalPlaces; i++) 
            oss << ((decimals / (power /= 10)) % 10);
    }
    return oss.str();
}

unsigned long long int PgHelpers::powuint(unsigned long long int base, int exp)
{
    if (exp == 0) return 1;
    if (exp == 1) return base;

    unsigned long long int tmp = PgHelpers::powuint(base, exp/2);
    if (exp%2 == 0) return tmp * tmp;
        else return base * tmp * tmp;
}

struct LUT
{
    char complementsLut[256];
    char upperComplementsLut[256];

    char* cLutPtr = complementsLut - CHAR_MIN;
    char* uLutPtr = upperComplementsLut - CHAR_MIN;

    LUT() {
        for(int i = CHAR_MIN; i < CHAR_MAX; i++)
            uLutPtr[i] = i;
        uLutPtr['A'] = 'T'; uLutPtr['a'] = 'T';
        uLutPtr['C'] = 'G'; uLutPtr['c'] = 'G';
        uLutPtr['G'] = 'C'; uLutPtr['g'] = 'C';
        uLutPtr['T'] = 'A'; uLutPtr['t'] = 'A';
        uLutPtr['N'] = 'N'; uLutPtr['n'] = 'N';
        uLutPtr['U'] = 'A'; uLutPtr['u'] = 'A';
        uLutPtr['Y'] = 'R'; uLutPtr['y'] = 'R';
        uLutPtr['R'] = 'Y'; uLutPtr['r'] = 'Y';
        uLutPtr['K'] = 'M'; uLutPtr['k'] = 'M';
        uLutPtr['M'] = 'K'; uLutPtr['m'] = 'K';
        uLutPtr['B'] = 'V'; uLutPtr['b'] = 'V';
        uLutPtr['D'] = 'H'; uLutPtr['d'] = 'H';
        uLutPtr['H'] = 'D'; uLutPtr['h'] = 'D';
        uLutPtr['V'] = 'B'; uLutPtr['v'] = 'B';
        uLutPtr['W'] = 'S'; uLutPtr['w'] = 'S';
        uLutPtr['S'] = 'W'; uLutPtr['s'] = 'W';
        for(int i = CHAR_MIN; i < CHAR_MAX; i++)
            cLutPtr[i] = uLutPtr[i];
        cLutPtr['U'] = 'U';
        cLutPtr['u'] = 'u';
        cLutPtr['a'] = 't';
        cLutPtr['c'] = 'g';
        cLutPtr['g'] = 'c';
        cLutPtr['t'] = 'a';
        cLutPtr['n'] = 'n';
        cLutPtr['y'] = 'r';
        cLutPtr['r'] = 'y';
        cLutPtr['k'] = 'm';
        cLutPtr['m'] = 'k';
        cLutPtr['b'] = 'v';
        cLutPtr['d'] = 'h';
        cLutPtr['h'] = 'd';
        cLutPtr['v'] = 'b';
        cLutPtr['w'] = 's';
        cLutPtr['s'] = 'w';
    }

    void removeUpperWScomplements() {
        upperComplementsLut['W'] = 'W'; upperComplementsLut['w'] = 'w';
        upperComplementsLut['S'] = 'S'; upperComplementsLut['s'] = 's';
    }
} instance;

char* upperComplementsLUT = instance.uLutPtr;
char* complementsLUT = instance.cLutPtr;

void PgHelpers::revertToMBGC121() {
    instance.removeUpperWScomplements();
}

char PgHelpers::upperReverseComplement(char symbol) {
    return upperComplementsLUT[symbol];
}

char PgHelpers::reverseComplement(char symbol) {
    return complementsLUT[symbol];
}

void PgHelpers::upperReverseComplementInPlace(char* start, const std::size_t N) {
    char* left = start - 1;
    char* right = start + N;
    while (--right > ++left) {
        char tmp = upperComplementsLUT[*left];
        *left = upperComplementsLUT[*right];
        *right = tmp;
    }
    if (left == right)
        *left = upperComplementsLUT[*left];
}

void PgHelpers::reverseComplementInPlace(char* start, const std::size_t N) {
    char* left = start - 1;
    char* right = start + N;
    while (--right > ++left) {
        char tmp = complementsLUT[*left];
        *left = complementsLUT[*right];
        *right = tmp;
    }
    if (left == right)
        *left = complementsLUT[*left];
}

void PgHelpers::upperReverseComplement(const char* start, const std::size_t N, char* target) {
    char* right = target + N;
    while (right-- > target) {
        *right = upperComplementsLUT[*(start++)];
    }
}

void PgHelpers::reverseComplement(const char* start, const std::size_t N, char* target) {
    char* right = target + N;
    while (right-- > target) {
        *right = complementsLUT[*(start++)];
    }
}

string PgHelpers::upperReverseComplement(string kmer) {
    size_t kmer_length = kmer.size();
    string revcomp;
    revcomp.resize(kmer_length);
    size_t j = kmer_length;
    for(size_t i = 0; i < kmer_length; i++)
        revcomp[--j] = upperComplementsLUT[kmer[i]];
    return revcomp;
}

string PgHelpers::reverseComplement(string kmer) {
    size_t kmer_length = kmer.size();
    string revcomp;
    revcomp.resize(kmer_length);
    size_t j = kmer_length;
    for(size_t i = 0; i < kmer_length; i++)
        revcomp[--j] = complementsLUT[kmer[i]];
    return revcomp;
}

void PgHelpers::upperReverseComplementInPlace(string &kmer) {
    upperReverseComplementInPlace((char *) kmer.data(), kmer.length());
}

void PgHelpers::reverseComplementInPlace(string &kmer) {
    reverseComplementInPlace((char *) kmer.data(), kmer.length());
}

void PgHelpers::upperSequence(char* start, const std::size_t N) {
    char* left = start - 1;
    char* guard = start + N;
    while (++left < guard) {
        *left = toupper(*left);
    }
}

double PgHelpers::qualityScore2approxCorrectProb(string quality) {
    double val = 1;
    for (char q : quality) {
        switch (q) {
            case 33:case 34:case 35:case 36: return 0;
            case 37: val *= 0.6018928294465028; break;
            case 38: val *= 0.683772233983162; break;
            case 39: val *= 0.748811356849042; break;
            case 40: val *= 0.800473768503112; break;
            case 41: val *= 0.8415106807538887; break;
            case 42: val *= 0.8741074588205833; break;
            case 43: val *= 0.9; break;
            case 44: val *= 0.9205671765275718; break;
            case 45: val *= 0.9369042655519807; break;
            case 46: val *= 0.9498812766372727; break;
            case 47: val *= 0.9601892829446502; break;
            case 48: val *= 0.9683772233983162; break;
            case 49: val *= 0.9748811356849042; break;
            case 50: val *= 0.9800473768503112; break;
            case 51: val *= 0.9841510680753889; break;
            case 52: val *= 0.9874107458820583; break;
            case 53: val *= 0.99; break;
            case 54: val *= 0.9920567176527572; break;
            case 55: val *= 0.993690426555198; break;
            case 56: val *= 0.9949881276637272; break;
            case 57: val *= 0.996018928294465; break;
            case 58: val *= 0.9968377223398316; break;
            case 59: val *= 0.9974881135684904; break;
            case 60: val *= 0.9980047376850312; break;
            case 61: val *= 0.9984151068075389; break;
            case 62: val *= 0.9987410745882058; break;
            case 63: val *= 0.999; break;
            case 64: val *= 0.9992056717652757; break;
            case 65: val *= 0.9993690426555198; break;
            case 66: val *= 0.9994988127663728; break;
            case 67: val *= 0.9996018928294464; break;
            case 68: val *= 0.9996837722339832; break;
            default: ;
        }
    }
    return pow(val, 1.0/quality.length());
}

double PgHelpers::qualityScore2correctProb(string quality) {
    double val = 1;
    for (char q : quality) {
        switch (q) {
            case 33: return 0;
            case 34: val *= 0.2056717652757185; break;
            case 35: val *= 0.36904265551980675; break;
            case 36: val *= 0.49881276637272776; break;
            case 37: val *= 0.6018928294465028; break;
            case 38: val *= 0.683772233983162; break;
            case 39: val *= 0.748811356849042; break;
            case 40: val *= 0.800473768503112; break;
            case 41: val *= 0.8415106807538887; break;
            case 42: val *= 0.8741074588205833; break;
            case 43: val *= 0.9; break;
            case 44: val *= 0.9205671765275718; break;
            case 45: val *= 0.9369042655519807; break;
            case 46: val *= 0.9498812766372727; break;
            case 47: val *= 0.9601892829446502; break;
            case 48: val *= 0.9683772233983162; break;
            case 49: val *= 0.9748811356849042; break;
            case 50: val *= 0.9800473768503112; break;
            case 51: val *= 0.9841510680753889; break;
            case 52: val *= 0.9874107458820583; break;
            case 53: val *= 0.99; break;
            case 54: val *= 0.9920567176527572; break;
            case 55: val *= 0.993690426555198; break;
            case 56: val *= 0.9949881276637272; break;
            case 57: val *= 0.996018928294465; break;
            case 58: val *= 0.9968377223398316; break;
            case 59: val *= 0.9974881135684904; break;
            case 60: val *= 0.9980047376850312; break;
            case 61: val *= 0.9984151068075389; break;
            case 62: val *= 0.9987410745882058; break;
            case 63: val *= 0.999; break;
            case 64: val *= 0.9992056717652757; break;
            case 65: val *= 0.9993690426555198; break;
            case 66: val *= 0.9994988127663728; break;
            case 67: val *= 0.9996018928294464; break;
            case 68: val *= 0.9996837722339832; break;
            case 69: val *= 0.999748811356849; break;
            case 70: val *= 0.9998004737685031; break;
            case 71: val *= 0.9998415106807539; break;
            case 72: val *= 0.9998741074588205; break;
            case 73: val *= 0.9999; break;
            default: val *= 1;
        }
    }
    return pow(val, 1.0/quality.length());
}

int PgHelpers::readsSufPreCmp(const char* suffixPart, const char* prefixRead) {
    while (*suffixPart) {
        if (*suffixPart > *prefixRead)
            return 1;
        if (*suffixPart++ < *prefixRead++)
            return -1;
    }
    return 0;
}

int PgHelpers::strcmplcp(const char* lStrPtr, const char* rStrPtr, int length) {
    
    int i = 0;
    while (length - i >= 4) {
        int cmp = bswap_32(*(uint32_t*) lStrPtr) - bswap_32(*(uint32_t*) rStrPtr);
        if (cmp != 0)
            break;
        lStrPtr += 4;
        rStrPtr += 4;
        i += 4;
    }

    while (i < length) {
        i++;
        int cmp = *(unsigned char*)(lStrPtr++) - *(unsigned char*)(rStrPtr++);
        if (cmp > 0)
            return 1;
        if (cmp < 0)
            return -1;
    }

    return 0;

}

unsigned long long PgHelpers::getTotalSystemMemory()
{
#ifdef __MINGW32__
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return status.ullTotalPhys;
#else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
#endif
}

#if defined(__arm__) || defined(__aarch64__) || defined(__ARM_ARCH)

void A_memcpy(void *dest, const void *src, size_t n) {
    memcpy(dest, src, n);
}

int A_memcmp(const void *s1, const void *s2, size_t n) {
    return memcmp(s1, s2, n);
}

#endif

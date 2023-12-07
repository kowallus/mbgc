/* ================================================================= *
*  SlidingWindowSparseEMMatcher implementation is based on copMEM:  *
*                                                                   *
*  copMEM is a program for efficient computation of MEMs            *
*  (Maximal Exact Matches) in a pair of genomes.                    *
*  Its main algorithmic idea requires that two internal parameters  *
*  (k1 and k2) are coprime, hence the name.                         *
*                                                                   *
*                                                                   *
*  Copyright (c) 2018, Szymon Grabowski and Wojciech Bieniecki      *
*  All rights reserved                                              *
*                                                                   *
*  This program is free software: you can redistribute it and/or    *
*  modify it under the terms of the GNU General Public License as   *
*  published by the Free Software Foundation, either version 3 of   *
*  the License, or (at your option) any later version.              *
*                                                                   *
*  This program is distributed in the hope that it will be useful,  *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
*  GNU General Public License for more details.                     *
*                                                                   *
*  You should have received a copy of the GNU General Public        *
*  License along with this program.                                 *
*                                                                   *
*  This file is subject to the terms and conditions defined in the  *
*  file 'license', which is part of this source code package.       *
* ================================================================= */

#include "SlidingWindowSparseEMMatcher.h"
#if !defined(__arm__) && !defined(__aarch64__) && !defined(__ARM_ARCH)
#include "../libs/asmlib.h"
#endif
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <cassert>
#include <cmath>

#define _prefetch(x,y) __builtin_prefetch(x,1,(4-y))

#include "../utils/Hashes.h"
#include <omp.h>

//////////////////// GLOBAL CONSTS //////////////////////////
const uint64_t NOT_MATCHED_POSITION = UINT64_MAX;

//////////////////// GLOBAL VARS ////////////////////////////

#define INIT_HASH_FUNC(n) hashFuncMatrix[n][1] = maRushPrime1HashSimplified<n>; hashFuncMatrix[n][2] = xxhash32<n>; hashFuncMatrix[n][3] = maRushPrime1HashSparsified<n>;hashFuncMatrix[n][4] = metroHash64<n>; hashFuncMatrix[n][5] = cityHash64<n>;
string hashNames[] = { "(?)", "maRushPrime1HashSimplified", "xxhash32", "maRushPrime1HashSparsified", "metroHash64", "cityHash64"};
void SlidingWindowSparseEMMatcher::initHashFuncMatrix() {
    INIT_HASH_FUNC(12);
    INIT_HASH_FUNC(16);
    INIT_HASH_FUNC(20);
    INIT_HASH_FUNC(24);
    INIT_HASH_FUNC(28);
    INIT_HASH_FUNC(32);
    INIT_HASH_FUNC(36);
    INIT_HASH_FUNC(40);
    INIT_HASH_FUNC(44);
    INIT_HASH_FUNC(56);
}


//////////////////// GLOBALS ////////////////////////////

std::ostream *v1logger;

void SlidingWindowSparseEMMatcher::initParams(uint32 minMatchLength, int _k1, int _k2) {
    if (L > 110) K = 56;
    else if (L > 62) K = 44;
    else if (L > 53) K = 40;
    else if (L > 46) K = 36;
    else if (L > 42) K = 32;
    else if (L > 32) K = 28;
    else K = (L / 4 - 1) * 4;
    if (minMatchLength < MIN_MATCH_LENGTH) {
        cerr << "Error: Minimal matching length too short!" << endl;
        exit(EXIT_FAILURE);
    }
    int KmmL = (minMatchLength / 4 - 1) * 4;
    if (KmmL < K) K = KmmL;
    calcCoprimes();
    k1 = _k1 != -1?_k1:k1;
    k2 = _k2 != -1?_k2:k2;

    hashFunc32 = hashFuncMatrix[K][H];
    v1logger = PgHelpers::devout;

    LK2 = (L - K) / 2;
    LK2_MINUS_4 = LK2 - 4;
    K_PLUS_LK24 = K + LK2_MINUS_4;

    uint8_t i = HASH_SIZE_MIN_ORDER;
    do {
        hash_size = ((uint32_t) 1) << (i++);
    } while (i <= HASH_SIZE_MAX_ORDER && hash_size < maxRefLength / k1);
    hash_size_minus_one = hash_size - 1;
}

void SlidingWindowSparseEMMatcher::displayParams() {
	*PgHelpers::devout << "SparseEM PARAMETERS: ";
	*PgHelpers::devout << "k-mer exact match length = " << L << "; ";
	*PgHelpers::devout << "hash pattern length = " << K << "; ";
	*PgHelpers::devout << "HASH_SIZE = " << hash_size << "; ";
	*PgHelpers::devout << "k1 = " << k1 << "; ";
	*PgHelpers::devout << "k2 = " << k2 << "; ";
    *PgHelpers::devout << "skipMargin = " << skipMargin << std::endl;
	*PgHelpers::devout << "Hash function: " << hashNames[H] << std::endl;
    *PgHelpers::devout << "Hash bin size: 1" << std::endl;
}

void SlidingWindowSparseEMMatcher::calcCoprimes()
{
	/* setting k1 and k2 */
	int tempVar = L - K + 1;
	if (tempVar <= 0) {
		std::cerr << "\nL and K mismatch.\n";
		exit(EXIT_FAILURE);
	}
	if (tempVar >= 20) {
        k1 = (int) (pow(tempVar, 0.5)) + 1;
        k2 = k1 - 1;
        if (k1 * k2 > tempVar) {
            --k2;
            --k1;
        }
    } else if (tempVar >= 15) {
	    k1 = 5; k2 = 3;
	} else if (tempVar >= 12) {
        k1 = 4; k2 = 3;
    } else if (tempVar >= 10) {
        k1 = 5; k2 = 2;
    } else if (tempVar >= 6) {
        k1 = 3; k2 = 2;
    } else {
        k1 = tempVar; k2 = 1;
    }
}

template<typename uint_ht_t>
void SlidingWindowSparseEMMatcher::processIgnoreCollisionsRef(uint_ht_t* ht, size_t& refSamplingPos, int64_t& refPos1) {
    const unsigned int MULTI2 = 128;
    const unsigned int k1MULTI2 = k1 * MULTI2;

    #pragma omp parallel for
    for (int64_t i1 = refSamplingPos; i1 < (refPos1 - K) - k1MULTI2; i1 += k1MULTI2) {
        uint32_t hashPositions[MULTI2];
        size_t i2 = i1;
        for (unsigned int temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            hashPositions[temp] = hashFunc(start1 + htRePos(i2));
            _prefetch((char*)(ht + hashPositions[temp]), 1);
        }

        i2 = i1;
        for (size_t temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            ht[hashPositions[temp]] = this->htEncodePos(i2);
        }
    }
    //////////////////// processing the end part of R
    int64_t i1;
    for (i1 = k1 + ((refPos1 - K - 1) / k1MULTI2) * k1MULTI2; i1 < refPos1 - K + 1; i1 += k1) {
        ht[hashFunc(start1 + htRePos(i1))] = this->htEncodePos(i1);
    }
    refSamplingPos = i1;
}

template<typename uint_ht_t>
void SlidingWindowSparseEMMatcher::processExactMatchQueryIgnoreCollisionsTightTemplate(
        uint_ht_t* ht, vector<TextMatch> &resMatches, const char *start2, size_t N2,
        bool destIsRef, bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos) {
    processExactMatchQueryIgnoreCollisionsTight<uint_ht_t, true>(ht, resMatches,
                                            start2, N2, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
}

template<typename uint_ht_t, bool checkOverlaps>
void SlidingWindowSparseEMMatcher::processExactMatchQueryIgnoreCollisionsTight(uint_ht_t* ht,
                vector<TextMatch> &resMatches, const char* start2, size_t N2, bool destIsRef, bool revComplMatching,
                uint32_t minMatchLength, size_t matchingLockPos) {
    if (destIsRef || revComplMatching) {
        fprintf(stderr, "Source as destination and reverse-complement matching modes unsupported\n\n");
        exit(EXIT_FAILURE);
    }

    uint32_t posArray;

    std::uint32_t l1 = 0, l2 = 0, r1 = 0, r2 = 0;

    size_t charExtensions = 0ULL;
    vector<TextMatch> replacedMatches;
    int64_t i2 = 0;
    const char* end2 = start2 + N2;
    const char* curr2 = start2 + i2;
    const char* guard2 = start2;
    for (; i2 + K < N2 + 1; i2 += k2, curr2 += k2) {
        posArray = hashFunc(curr2);

        if (curr2 - LK2 >= start2) memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
        if (curr2 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end2) memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

        uint32_t j = posArray;
        ++charExtensions;
        if (ht[j] == 0)
            continue;

        const char* curr1 = start1 + this->htDecodePos(ht[j]);
        const char* swStart = start1 + pos1;
        const char* swStop = start1 + (matchingLockPos != SW_END_ERASED_FLAG ? matchingLockPos : pos1);
        bool matchEndsBeforeSwStart = curr1 + K < swStart;
        bool matchStartsBeforeSwStop = curr1 < swStop;
        if (swStart <= swStop) {
            if (!matchEndsBeforeSwStart && matchStartsBeforeSwStop)
                continue;
        } else if (!matchEndsBeforeSwStart || matchStartsBeforeSwStop)
            continue;
        const char* tmpEnd1 = matchEndsBeforeSwStart ? swStart : start1 + this->getRefLength();
        const char* tmpStart1 = matchStartsBeforeSwStop ? start1 : swStop;
        memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
        memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

        if (r1 == r2 || l1 == l2) {
            const char* p1 = curr1 + K;
            const char* p2 = curr2 + K;
            const int MEMCMP_BLOCK = 64;
            size_t blocks1 = (tmpEnd1 - p1) / MEMCMP_BLOCK;
            size_t blocks2 = (end2 - p2) / MEMCMP_BLOCK;
            size_t blocks = blocks1 < blocks2 ? blocks1 : blocks2;
            blocks = blocks >> 1 << 1;
            if (blocks) blocks--;
            for (int i = 0; i < blocks; i += 2, p1 += 2 * MEMCMP_BLOCK, p2 += 2 * MEMCMP_BLOCK) {
                if (A_memcmp(p1, p2, MEMCMP_BLOCK) != 0)
                    break;
                if (A_memcmp(p1 + MEMCMP_BLOCK, p2 + MEMCMP_BLOCK, MEMCMP_BLOCK) != 0) {
                    p1 += MEMCMP_BLOCK, p2 += MEMCMP_BLOCK;
                    break;
                }
            }
            p1--; p2--;
            while (++p1 != tmpEnd1 && ++p2 != end2 && *p1 == *p2);
            const char *right1 = p1;
            const char *right2 = p2;
            p1 = curr1;
            p2 = curr2;

            const int PAIRED_MATCH_LENGTH_LOSS_LIMIT = 0; // 0 - DISABLED
            int64_t resSizeWithoutOverlapped = resMatches.size();
            if (checkOverlaps) {
                while (resSizeWithoutOverlapped--) {
                    const TextMatch& lastMatch = resMatches[resSizeWithoutOverlapped];
                    if (lastMatch.endPosDestText() < p2 - start2) {
                        const char* tmpguard2 = start2 + lastMatch.endPosDestText();
                        while (p1 != tmpStart1 && p2 > tmpguard2 - 1 && *p1 == *p2) {
                            p1--; p2--;
                        }
                        if (p2 > tmpguard2 - 1)
                            break;
                        p1++; p2++;
                    }
                    int64_t lastDelta = (p2 - start2) - lastMatch.posDestText;
                    if (p1 - tmpStart1 < lastDelta || lastMatch.length > OVERLAP_MATCH_MAX_LENGTH ||
                        PgHelpers::strcmplcp(p2 - lastDelta, p1 - lastDelta, lastDelta) != 0) {
                        p1--; p2--;
                        break;
                    }
                    p2 -= lastDelta; p1 -= lastDelta;
                    if (PAIRED_MATCH_LENGTH_LOSS_LIMIT && !replacedMatches.empty() &&
                        (right2 - p2 - 1) - replacedMatches.back().length >= PAIRED_MATCH_LENGTH_LOSS_LIMIT)
                        replacedMatches.clear();
                }
                if (PAIRED_MATCH_LENGTH_LOSS_LIMIT && resSizeWithoutOverlapped + 2 < resMatches.size())
                    replacedMatches.clear();
                if (resSizeWithoutOverlapped < 0) {
                    while (p1 != tmpStart1 && p2 > guard2 - 1 && *p1 == *p2) {
                        p1--; p2--;
                    }
                } else {
                    if (PAIRED_MATCH_LENGTH_LOSS_LIMIT && resSizeWithoutOverlapped + 1 < resMatches.size() && resMatches[resSizeWithoutOverlapped + 1].pairedWith(resMatches[resSizeWithoutOverlapped]))
                        if ((right2 - p2 - 1) - resMatches[resSizeWithoutOverlapped + 1].length < PAIRED_MATCH_LENGTH_LOSS_LIMIT)
                            replacedMatches.push_back(resMatches[resSizeWithoutOverlapped + 1]);
                    int64_t overlap = ((int64_t) resMatches[resSizeWithoutOverlapped].endPosDestText()) - (p2 + 1 - start2);
                    if (overlap > 0) {
                        p1 += overlap;
                        p2 += overlap;
                    }
                }
                ++resSizeWithoutOverlapped;
            } else {
                while (p1 != tmpStart1 && p2 > guard2 - 1 && *p1 == *p2) {
                    p1--; p2--;
                }
//                while (p2 < guard2 - 1) { p1++; p2++; } // to definitely remove overlapping matches
            }
            if (right1 - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                resMatches.resize(resSizeWithoutOverlapped);
                resMatches.push_back(TextMatch(p1 + 1 - start1, right1 - p1 - 1, p2 + 1 - start2));
                if (!replacedMatches.empty() && resMatches.size() > 1 &&
                    resMatches.back().posDestText > replacedMatches.back().endPosDestText()) {
                    int64_t m = resMatches.size() - 2;
                    if (!resMatches.back().pairedWith(resMatches[m]))
                        resMatches[m] = replacedMatches.back();
                    replacedMatches.clear();
                }
                int skip = ((p2 - start2 + right1 - p1) - (curr2 - start2)) / k2 * k2;
                skip -= skip > skipMargin ? skipMargin : skip;
                if (skip) {
                    curr2 += skip - k2;
                    i2 += skip - k2;
                }
                if (!checkOverlaps) guard2 = p2 + (right1 - p1); // to avoid large overlaps
            }

        }
    }
    //////////////////// processing the end part of Q  //////////////////////
//    *v1logger << "Character extensions = " << charExtensions <<  "\n";
}

using namespace std;

SlidingWindowSparseEMMatcher::SlidingWindowSparseEMMatcher(const size_t refLengthLimit,
                                                           const uint32_t targetMatchLength,
                                                           int _k1, int _k2, int skipMargin, uint32_t minMatchLength,
                                                           bool skipHtInit)
        : maxRefLength(refLengthLimit), L(targetMatchLength), skipMargin(skipMargin) {
    size_t allocSize;
    allocSize = PgHelpers::safeNewArrayAlloc(start1, this->maxRefLength, false,
                                                    TOTAL_RAM_REF_LIMIT_PERCENT);
    if (this->maxRefLength != allocSize)
        this->maxRefLength = allocSize;
    start1[0] = 0;

    pos1 = REF_SHIFT;
    swEnd = this->maxRefLength;
    setSlidingWindowSize(SW_WIDTH_FACTOR);

    initHashFuncMatrix();
    if (minMatchLength > targetMatchLength)
        minMatchLength = targetMatchLength;
    initParams(minMatchLength, _k1, _k2);
    if (skipHtInit)
        return;
    if (refLengthLimit >= (1ULL << 32)) {
        ht64bitFlag = true;  // huge Reference
        allocSize = PgHelpers::safeNewArrayAlloc<uint64_t>(ht64bit, hash_size, true);
        *v1logger << "WARNING - BIG reference file (>4GB), 64-bit arrays\n";
    } else {
        allocSize = PgHelpers::safeNewArrayAlloc<uint32_t>(ht32bit, hash_size, true);
    }
    if (hash_size != allocSize) {
        hash_size = allocSize;
        hash_size_minus_one = hash_size - 1;
    }
    displayParams();
}

size_t SlidingWindowSparseEMMatcher::acquireWorkerMatchingLockPos() {
    if (swSize == 0 || !circularBuffer)
        return swEnd;
    size_t workerSwEndPos;
    #pragma omp critical(workerLocksManagement)
    {
        workerSwEndPos = pos1 + swSize;
        if (reachedRefLengthCount || workerSwEndPos > this->maxRefLength) {
            if (workerSwEndPos > maxRefLength)
                workerSwEndPos -= maxRefLength - REF_SHIFT;
        } else
            workerSwEndPos = maxRefLength;
        if (workersSwEndPositions.empty())
            swEnd = workerSwEndPos;
        workersSwEndPositions.push_back(workerSwEndPos);
    }
    return workerSwEndPos;
}

void SlidingWindowSparseEMMatcher::releaseWorkerMatchingLockPos(size_t lockValue) {
    if (swSize == 0 || !circularBuffer)
        return;
    #pragma omp critical(workerLocksManagement)
    {
        int i = 0;
        while (i < workersSwEndPositions.size() && workersSwEndPositions[i] != lockValue)
            i++;
        if (i == workersSwEndPositions.size()) {
            fprintf(stderr, "ERROR: Invalid worker lock value (%ld)\n\n", lockValue);
            exit(EXIT_FAILURE);
        }
        if (i == 0) {
            do {
                workersSwEndPositions.pop_front();
            } while (!workersSwEndPositions.empty() && workersSwEndPositions.front() == SW_END_ERASED_FLAG);
            swEnd = workersSwEndPositions.empty() ? swEnd : workersSwEndPositions.front();
        } else
            workersSwEndPositions[i] = SW_END_ERASED_FLAG;
    }
}

size_t SlidingWindowSparseEMMatcher::loadRef(const char *refText, size_t refLength, size_t& refSamplingPos,
                                           int64_t& refPos1, bool rcRef, bool addRegionSeparators,
                                           char regionSeparator) {
    if (refLength == 0)
        return 0;
    if (refPos1 == this->maxRefLength && swEnd != this->maxRefLength) {
        reachedRefLengthCount++;
        refPos1 = REF_SHIFT;
        samplingPos = REF_SHIFT;
    }
    size_t tmpEnd = swEnd;
    size_t tmpLength = refLength;
    size_t tmpMax = tmpEnd < refPos1 ? maxRefLength : tmpEnd;
    if (refPos1 + tmpLength > tmpMax) {
        tmpLength = tmpMax - refPos1;
    }
    if (rcRef) {
        PgHelpers::upperReverseComplement(refText + refLength - tmpLength, tmpLength, start1 + refPos1);
    } else
        A_memcpy(start1 + refPos1, refText, tmpLength);
    if (addRegionSeparators && refPos1 + tmpLength == swEnd)
        start1[swEnd - 1] = regionSeparator;
    refPos1 += tmpLength;

    if (ht64bitFlag) {
        processIgnoreCollisionsRef<std::uint64_t>(ht64bit, refSamplingPos, refPos1);
    } else {
        processIgnoreCollisionsRef<std::uint32_t>(ht32bit, refSamplingPos, refPos1);
    }

    refText += rcRef ? 0 : tmpLength;
    refLength = refPos1 == tmpEnd ? 0 : refLength - tmpLength;

    return tmpLength + this->loadRef(refText, refLength, refSamplingPos, refPos1, rcRef, addRegionSeparators,
                                     regionSeparator);
}

void SlidingWindowSparseEMMatcher::loadSeparator(char regionSeparator) {
    if (pos1 == this->maxRefLength && swEnd != this->maxRefLength) {
        reachedRefLengthCount++;
        pos1 = REF_SHIFT;
        samplingPos = REF_SHIFT;
    }
    if (pos1 == this->maxRefLength)
        return;
    if (pos1 == swEnd) {
        start1[pos1 - 1] = regionSeparator;
    } else
        start1[pos1++] = regionSeparator;
}

void SlidingWindowSparseEMMatcher::loadRef(const char *refText, size_t refLength, bool loadRCRef,
                                           bool addRegionSeparators, char regionSeparator) {
    this->loadRef(refText, refLength, samplingPos, pos1, false, addRegionSeparators, regionSeparator);
    if (loadRCRef)
        this->loadRef(refText, refLength, samplingPos, pos1, true, addRegionSeparators, regionSeparator);
}

SlidingWindowSparseEMMatcher::~SlidingWindowSparseEMMatcher() {
    if (ht64bitFlag)
        delete[] ht64bit;
    else
        delete[] ht32bit;

    delete[] start1;
}


void
SlidingWindowSparseEMMatcher::matchTexts(vector <TextMatch> &resMatches, const string &destText, bool destIsRef, bool revComplMatching,
                                         uint32_t minMatchLength, size_t matchingLockPos) {
    matchTexts(resMatches, destText.data(), destText.size(), destIsRef, revComplMatching, minMatchLength,
               matchingLockPos);
}

void
SlidingWindowSparseEMMatcher::matchTexts(vector<TextMatch> &resMatches, const char* destText, size_t destLen, bool destIsRef,
                                         bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos) {
    if (minMatchLength < K) {
        fprintf(stderr, "Minimal matching length cannot be smaller than K (%d < %d)\n\n", minMatchLength, K);
        exit(EXIT_FAILURE);
    }
    resMatches.clear();
    if (ht64bitFlag) {
        processExactMatchQueryIgnoreCollisionsTightTemplate<std::uint64_t>(ht64bit, resMatches,
                            destText, destLen, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
    } else {
        processExactMatchQueryIgnoreCollisionsTightTemplate<std::uint32_t>(ht32bit, resMatches,
                            destText, destLen, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
    }
}

SlidingWindowExpSparseEMMatcher::SlidingWindowExpSparseEMMatcher(const size_t refLengthLimit,
                                                                 const uint32_t targetMatchLength, int _k1, int _k2,
                                                                 int skipMargin, uint32_t minMatchLength) : SlidingWindowSparseEMMatcher(
        refLengthLimit, targetMatchLength, _k1, _k2, skipMargin, minMatchLength, true),
        k1ord(__builtin_ctz((uint32_t) _k1)) {
    if (_k1 % 2) {
        fprintf(stderr, "Error initializing ExpSparseMEM: incorrect k1 (%d)\n\n", _k1);
        exit(EXIT_FAILURE);
    }
    this->samplingPos = _k1;
    size_t allocSize;
    if (this->htEncodePos(refLengthLimit) >= (1ULL << 32)) {
        ht64bitFlag = true;  // large Reference
        allocSize = PgHelpers::safeNewArrayAlloc<uint64_t>(ht64bit, hash_size, true);
        *v1logger << "WARNING - LARGE reference file (SIZE / k1ord > 4GB), 64-bit arrays\n";
    } else {
        allocSize = PgHelpers::safeNewArrayAlloc<uint32_t>(ht32bit, hash_size, true);
    }
    if (hash_size != allocSize) {
        hash_size = allocSize;
        hash_size_minus_one = hash_size - 1;
    }
    SlidingWindowSparseEMMatcher::displayParams();
    *PgHelpers::devout << "Exponential mode: ";
    *PgHelpers::devout << "k1 order = " << k1ord << std::endl;
}

SlidingWindowExpRandSparseEMMatcher::SlidingWindowExpRandSparseEMMatcher(const size_t refLengthLimit,
                                                                 const uint32_t targetMatchLength, int _k1, int _k2,
                                                                 int skipMargin, uint32_t minMatchLength) : SlidingWindowExpSparseEMMatcher(
        refLengthLimit, targetMatchLength, _k1, _k2, skipMargin, minMatchLength) {
    if (_k1 < __builtin_popcount(POPCOUNT_MASK)) {
        fprintf(stderr, "Error initializing RandExpSparseMEM: Too small k1 (%d < %d)\n\n", _k1, __builtin_popcount(POPCOUNT_MASK));
        exit(EXIT_FAILURE);
    }
    *PgHelpers::appout << "Randomized mode!" << std::endl;
}
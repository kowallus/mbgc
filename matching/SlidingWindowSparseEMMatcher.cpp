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

#include "../libs/asmlib.h"

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
    if (minMatchLength < 24) {
        cout << "Error: Minimal matching length too short!" << endl;
        exit(EXIT_FAILURE);
    }
    int KmmL = (minMatchLength / 4 - 1) * 4;
    if (KmmL < K) K = KmmL;
    calcCoprimes();
    k1 = _k1 != -1?_k1:k1;
    k2 = _k2 != -1?_k2:k2;

    hashFunc32 = hashFuncMatrix[K][H];
    v1logger = PgHelpers::logout;

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
	std::cout << "SparseEM PARAMETERS: ";
	std::cout << "k-mer exact match length = " << L << "; ";
	std::cout << "hash pattern length = " << K << "; ";
	std::cout << "HASH_SIZE = " << hash_size << "; ";
	std::cout << "k1 = " << k1 << "; ";
	std::cout << "k2 = " << k2 << "; ";
    std::cout << "skipMargin = " << skipMargin << std::endl;
	*PgHelpers::logout << "Hash function: " << hashNames[H] << std::endl;
    *PgHelpers::logout << "Hash bin size: 1" << std::endl;
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

template<typename MyUINT1, typename MyUINT2>
void SlidingWindowSparseEMMatcher::processIgnoreCollisionsRef(HashBuffer<MyUINT1, MyUINT2>& buffer) {
    const unsigned int MULTI2 = 128;
    const unsigned int k1MULTI2 = k1 * MULTI2;

    MyUINT1* sampledPositions = buffer.first;

    #pragma omp parallel for
    for (int64_t i1 = samplingPos; i1 < (pos1 - K) - k1MULTI2; i1 += k1MULTI2) {
        uint32_t hashPositions[MULTI2];
        const char* tempPointer = start1 + i1;
        for (unsigned int temp = 0; temp < MULTI2; ++temp) {
            hashPositions[temp] = hashFunc(tempPointer);
            tempPointer += k1;
            _prefetch((char*)(sampledPositions + hashPositions[temp]), 1);
        }

        int64_t i2 = i1;
        for (size_t temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            sampledPositions[hashPositions[temp]] = i2;
        }
    }

    //////////////////// processing the end part of R
    int64_t i1;
    for (i1 = REF_SHIFT + ((pos1 - K - 1) / k1MULTI2) * k1MULTI2; i1 < pos1 - K + 1; i1 += k1) {
        sampledPositions[hashFunc(start1 + i1)] = i1;
    }
    samplingPos = i1;
}

template <class MyUINT1, class MyUINT2>
void SlidingWindowSparseEMMatcher::deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf) {
    delete[] buf.first;
    delete[] buf.second;
}

template<typename MyUINT1, typename MyUINT2>
void SlidingWindowSparseEMMatcher::processExactMatchQueryIgnoreCollisionsTightTemplate(
        HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches, const char *start2, size_t N2,
        bool destIsRef, bool revComplMatching, uint32_t minMatchLength, size_t matchingLockPos) {
    processExactMatchQueryIgnoreCollisionsTight<MyUINT1, MyUINT2, true>(buffer, resMatches,
            start2, N2, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
}


template<typename MyUINT1, typename MyUINT2, bool checkOverlaps>
void SlidingWindowSparseEMMatcher::processExactMatchQueryIgnoreCollisionsTight(HashBuffer<MyUINT1, MyUINT2> buffer,
           vector<TextMatch> &resMatches, const char* start2, size_t N2, bool destIsRef, bool revComplMatching,
           uint32_t minMatchLength, size_t matchingLockPos) {
    if (destIsRef || revComplMatching) {
        fprintf(stderr, "Source as destination and reverse-complement matching modes unsupported\n\n");
        exit(EXIT_FAILURE);
    }

    MyUINT1* sampledPositions = buffer.first;
    MyUINT2 posArray;

    std::uint32_t l1 = 0, l2 = 0, r1 = 0, r2 = 0;

    size_t charExtensions = 0ULL;

    int64_t i2 = 0;
    const char* end2 = start2 + N2;
    const char* curr2 = start2 + i2;
    const char* guard2 = start2;
    for (; i2 + K < N2 + 1; i2 += k2, curr2 += k2) {
        posArray = hashFunc(curr2);

        if (curr2 - LK2 >= start2) memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
        if (curr2 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end2) memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

        MyUINT1 j = posArray;
        ++charExtensions;
        if (sampledPositions[j] == 0)
            continue;

        const char* curr1 = start1 + sampledPositions[j];
        const char* swStart = start1 + pos1;
        const char* swStop = start1 + (matchingLockPos != SW_END_ERASED_FLAG ? matchingLockPos : pos1);
        bool matchBeforeSwStart = curr1 < swStart;
        bool matchBeforeSwStop = curr1 < swStop;
        if (swStart <= swStop) {
            if (!matchBeforeSwStart && matchBeforeSwStop)
                continue;
        } else if (!matchBeforeSwStart || matchBeforeSwStop)
            continue;
        const char* tmpEnd1 = matchBeforeSwStart ? swStart : start1 + this->getRefLength();
        const char* tmpStart1 = matchBeforeSwStop ? start1 : swStop;
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
            const char *right = p1;
            p1 = curr1;
            p2 = curr2;

            if (checkOverlaps) {
                while (!resMatches.empty()) {
                    const TextMatch& lastMatch = resMatches.back();
                    if (lastMatch.endPosDestText() < p2 - start2) {
                        const char* tmpguard2 = start2 + lastMatch.endPosDestText();
                        while (p1 != tmpStart1 && p2 > tmpguard2 - 1 && *--p1 == *--p2);
                        if (p2 > tmpguard2 - 1)
                            break;
                        p1++; p2++;
                    }
                    int64_t lastDelta = (p2 - start2) - lastMatch.posDestText;
                    if (p1 - start1 < lastDelta || lastMatch.length > OVERLAP_MATCH_MAX_LENGTH ||
                        PgHelpers::strcmplcp(p2 - lastDelta, p1 - lastDelta, lastDelta) != 0) {
                        p1--; p2--;
                        break;
                    }
                    p2 -= lastDelta; p1 -= lastDelta;
                    resMatches.pop_back();
                }
                if (resMatches.empty())
                    while (p1 != tmpStart1 && p2 > guard2 - 1 && *--p1 == *--p2);
                else {
                    int64_t overlap = ((int64_t) resMatches.back().endPosDestText()) - (p2 + 1 - start2);
                    if (overlap > 0) {
                        p1 += overlap;
                        p2 += overlap;
                    }
                }
            } else {
                while (p1 != tmpStart1 && p2 > guard2 - 1 && *--p1 == *--p2);
//                while (p2 < guard2 - 1) { p1++; p2++; } // to definitely remove overlapping matches
            }
            if (right - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                resMatches.push_back(TextMatch(p1 + 1 - start1, right - p1 - 1, p2 + 1 - start2));
                int skip = ((p2 - start2 + right - p1) - (curr2 - start2)) / k2 * k2;
                skip -= skip > skipMargin ? skipMargin : skip;
                if (skip) {
                    curr2 += skip - k2;
                    i2 += skip - k2;
                }
                if (!checkOverlaps) guard2 = p2 + (right - p1); // to avoid large overlaps
            }

        }
    }
    //////////////////// processing the end part of Q  //////////////////////

//    *v1logger << "Character extensions = " << charExtensions <<  "\n";
}

using namespace std;

SlidingWindowSparseEMMatcher::SlidingWindowSparseEMMatcher(const size_t refLengthLimit, const uint32_t targetMatchLength,
                                                           int _k1, int _k2, int skipMargin, uint32_t minMatchLength)
    : maxRefLength(refLengthLimit), L(targetMatchLength), skipMargin(skipMargin) {
    start1 = new char[this->maxRefLength];
    start1[0] = 0;
    pos1 = REF_SHIFT;
    swEnd = this->maxRefLength;
    setSlidingWindowSize(SW_WIDTH_FACTOR);

    initHashFuncMatrix();
    if (minMatchLength > targetMatchLength)
        minMatchLength = targetMatchLength;
    initParams(minMatchLength, _k1, _k2);
    displayParams();
    if (refLengthLimit / k1 >= (1ULL << 32)) {
        bigRef = 2;  // huge Reference
        buffer2.first = new uint64_t[hash_size]();
        *v1logger  << "WARNING - LARGE reference file (SIZE / k1 > 4GB), 64-bit arrays\n";
    } else if (refLengthLimit >= (1ULL << 32)) {
        bigRef = 1;  // large Reference
        buffer1.first = new uint64_t[hash_size]();
        *v1logger << "WARNING - BIG reference file (>4GB), 64-bit arrays\n";
    } else {
        bigRef = 0;  // small Reference
        buffer0.first = new uint32_t[hash_size]();
    }
}

size_t SlidingWindowSparseEMMatcher::acquireWorkerMatchingLockPos() {
    if (swSize == 0)
        return this->maxRefLength;
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
    if (swSize == 0)
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

void SlidingWindowSparseEMMatcher::loadRef(const char *refText, size_t refLength) {
    if (refLength == 0)
        return;
    if (pos1 == this->maxRefLength && swEnd != this->maxRefLength) {
        reachedRefLengthCount++;
        pos1 = REF_SHIFT;
    }
    size_t tmpEnd = swEnd;
    size_t tmpLength = refLength;
    size_t tmpMax = tmpEnd < pos1 ? maxRefLength : tmpEnd;
    if (pos1 + tmpLength > tmpMax) {
        tmpLength = tmpMax - pos1;
    }
    A_memcpy(start1 + pos1, refText, tmpLength);
    pos1 += tmpLength;
    if (bigRef == 2) {
        processIgnoreCollisionsRef<std::uint64_t, std::uint64_t>(buffer2);
    } else if (bigRef == 1) {
        processIgnoreCollisionsRef<std::uint64_t, std::uint32_t>(buffer1);
    } else {
        processIgnoreCollisionsRef<uint32_t, uint32_t>(buffer0);
    }

    refText += tmpLength;
    refLength = pos1 == tmpEnd ? 0 : refLength - tmpLength;

    this->loadRef(refText, refLength);
}

SlidingWindowSparseEMMatcher::~SlidingWindowSparseEMMatcher() {
    if (bigRef == 2)
        deleteHashBuffer(buffer2);
    else if (bigRef == 1)
        deleteHashBuffer(buffer1);
    else if (bigRef == 0)
        deleteHashBuffer(buffer0);

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
    if (bigRef == 2) {
        processExactMatchQueryIgnoreCollisionsTightTemplate<std::uint64_t, std::uint64_t>(buffer2, resMatches,
                  destText, destLen, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
    } else if (bigRef == 1) {
        processExactMatchQueryIgnoreCollisionsTightTemplate<std::uint64_t, std::uint32_t>(buffer1, resMatches,
                  destText, destLen, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
    }
    else {
        processExactMatchQueryIgnoreCollisionsTightTemplate<std::uint32_t, std::uint32_t>(buffer0, resMatches,
                  destText, destLen, destIsRef, revComplMatching, minMatchLength, matchingLockPos);
    }
}

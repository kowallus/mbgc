/* ================================================================= *
*  SparseEMMatcher is  based on copMEM:                               *
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

#include "SparseEMMatcher.h"

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
void SparseEMMatcher::initHashFuncMatrix() {
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

void SparseEMMatcher::initParams(uint32 minMatchLength, int _k1, int _k2) {
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
    } while (i <= HASH_SIZE_MAX_ORDER && hash_size < maxSrcLength / k1);
    hash_size_minus_one = hash_size - 1;
}

void SparseEMMatcher::displayParams() {
	std::cout << "SparseEM PARAMETERS: ";
	std::cout << "k-mer exact match length = " << L << "; ";
	std::cout << "hash pattern length = " << K << "; ";
	std::cout << "HASH_SIZE = " << hash_size << "; ";
	std::cout << "k1 = " << k1 << "; ";
	std::cout << "k2 = " << k2 << std::endl;
	*PgHelpers::logout << "Hash function: " << hashNames[H] << std::endl;
    *PgHelpers::logout << "Hash bin size: "
        << (collisionSupport?HASH_COLLISIONS_PER_POSITION_LIMIT:1) << std::endl;
    if (collisionSupport) {
        *PgHelpers::logout << "Average hash collisions per position limit (in approx mode): "
                           << AVERAGE_HASH_COLLISIONS_PER_POSITION_LIMIT << std::endl;
        *PgHelpers::logout << "Unlimited number of hash collisions per position (in approx mode): "
                           << UNLIMITED_NUMBER_OF_HASH_COLLISIONS_PER_POSITION << std::endl;
    }
}

void SparseEMMatcher::calcCoprimes()
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
void SparseEMMatcher::genCumm(size_t N, const char* gen, MyUINT2* cumm, vector<MyUINT1> &skippedList) {
    const size_t MULTI1 = 128;
    const size_t k1MULTI1 = k1 * MULTI1;

    uint32_t hashPositions[MULTI1];
    size_t i;

    for (i = 0; i + K + k1MULTI1 < N ; i += k1MULTI1) {
        const char* tempPointer = gen + i;
        for (size_t temp = 0; temp < MULTI1; ++temp) {
            hashPositions[temp] = hashFunc(tempPointer) + 2;
            tempPointer += k1;
            _prefetch((char*)(cumm + hashPositions[temp]), 1);
        }

        for (size_t temp = 0; temp < MULTI1; ++temp) {
            if (cumm[hashPositions[temp]] <= HASH_COLLISIONS_PER_POSITION_LIMIT)
                ++cumm[hashPositions[temp]];
            else
                skippedList.push_back(i + k1 * temp);
        }
    }

    //////////////////// processing the end part of R  //////////////////////
    for (; i < N - K + 1; i += k1) {
        uint32_t h = hashFunc(gen + i) + 2;
        if (cumm[h] <= HASH_COLLISIONS_PER_POSITION_LIMIT)
            ++cumm[h];
        else
            skippedList.push_back(i);
    }
    //////////////////// processing the end part of R //////////////////////
    std::partial_sum(cumm, cumm + hash_size + 2, cumm);
    skippedList.push_back(N);
}

template<typename MyUINT1, typename MyUINT2>
void SparseEMMatcher::processRef(HashBuffer<MyUINT1, MyUINT2>& buffer) {
    if (!collisionSupport)
        return processIgnoreCollisionsRef<MyUINT1, MyUINT2>(buffer);

    if (PgHelpers::numberOfThreads > 1)
        return processRefMultithreaded<MyUINT1, MyUINT2>(buffer);

    const unsigned int MULTI2 = 128;
    const unsigned int k1MULTI2 = k1 * MULTI2;

    MyUINT2* cumm = new MyUINT2[hash_size + 2]();
    vector<MyUINT1> skippedList;
    const char* start1 = src.data() + SRC_SHIFT;
    const size_t N = src.size() - SRC_SHIFT;
    genCumm(N, start1, cumm, skippedList);
    const size_t hashCount = cumm[hash_size + 1];
    MyUINT1* sampledPositions = new MyUINT1[hashCount + 2];
    *v1logger << "Hash count = " << hashCount << std::endl;

    uint32_t hashPositions[MULTI2];
    MyUINT1 i1;

    size_t s = 0;
    for (i1 = 0; i1 + K + k1MULTI2 < N ; i1 += k1MULTI2) {
        const char* tempPointer = start1 + i1;
        for (unsigned int temp = 0; temp < MULTI2; ++temp) {
            hashPositions[temp] = hashFunc(tempPointer) + 1;
            tempPointer += k1;
            _prefetch((char*)(cumm + hashPositions[temp]), 1);
        }

        for (unsigned int temp = 0; temp < MULTI2; ++temp) {
            _prefetch((char*)(sampledPositions + *(cumm + hashPositions[temp])), 1);
        }

        MyUINT1 i2 = i1;
        for (size_t temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            if (skippedList[s] == i2) {
                s++;
                continue;
            }
            sampledPositions[cumm[hashPositions[temp]]] = i2;
            ++cumm[hashPositions[temp]];
        }
    }

    //////////////////// processing the end part of R
    for (; i1 < N - K + 1; i1 += k1) {
        if (skippedList[s] == i1) {
            s++;
            continue;
        }
        uint32_t h = hashFunc(start1 + i1) + 1;
        sampledPositions[cumm[h]] = i1;
        ++cumm[h];
    }

    buffer = { sampledPositions, cumm };
}

template<typename MyUINT2>
void SparseEMMatcher::genCummMultithreaded(size_t N, const char* gen, uint8_t* counts, MyUINT2* cumm) {
	const size_t MULTI1 = 128;
	const size_t k1MULTI1 = k1 * MULTI1;

	#pragma omp parallel for
    for (size_t i = 0; i < N - K - k1MULTI1; i += k1MULTI1) {
        uint32_t hashPositions[MULTI1];
		const char* tempPointer = gen + i;
		for (size_t temp = 0; temp < MULTI1; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer);
			tempPointer += k1;
			_prefetch((char*)(counts + hashPositions[temp]), 1);
		}

        for (size_t temp = 0; temp < MULTI1; ++temp) {
            {
                if (counts[hashPositions[temp]] <= HASH_COLLISIONS_PER_POSITION_LIMIT)
                    ++counts[hashPositions[temp]];
            }
        }
	}

	//////////////////// processing the end part of R  //////////////////////
	for (size_t i = ((N - K - 1) / k1MULTI1) * k1MULTI1; i < N - K + 1; i += k1) {
		uint32_t h = hashFunc(gen + i);
		if (counts[h] <= HASH_COLLISIONS_PER_POSITION_LIMIT)
		    ++counts[h];
	}
	//////////////////// processing the end part of R //////////////////////
    MyUINT2 sum = 0;
    for(size_t i = 0; i < hash_size + 1; i++) {
        cumm[i] = sum;
        sum += counts[i];
    }
    cumm[hash_size + 1] = sum;
}

template<typename MyUINT1, typename MyUINT2>
void SparseEMMatcher::processRefMultithreaded(HashBuffer<MyUINT1, MyUINT2>& buffer) {

	const unsigned int MULTI2 = 128;
	const unsigned int k1MULTI2 = k1 * MULTI2;

	MyUINT2* cumm = new MyUINT2[hash_size + 2]();
    uint8_t* counts = new uint8_t[hash_size + 2]();
    const char* start1 = src.data() + SRC_SHIFT;
    const size_t N = src.size() - SRC_SHIFT;
    genCummMultithreaded(N, start1, counts, cumm);
    const size_t hashCount = cumm[hash_size + 1];
    MyUINT1* sampledPositions = new MyUINT1[hashCount + 2];
    *v1logger << "Hash count = " << hashCount << std::endl;

    #pragma omp parallel for
	for (MyUINT1 i1 = 0; i1 < N - K - k1MULTI2 ; i1 += k1MULTI2) {
        uint32_t hashPositions[MULTI2];
		const char* tempPointer = start1 + i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
		    hashPositions[temp] = hashFunc(tempPointer);
			tempPointer += k1;
            _prefetch((char*)(counts + hashPositions[temp]), 1);
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			_prefetch((char*)(sampledPositions + *(cumm + hashPositions[temp])), 1);
		}

		MyUINT1 i2 = i1;
		for (size_t temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            if (!counts[hashPositions[temp]])
                continue;
			sampledPositions[cumm[hashPositions[temp]] + (--counts[hashPositions[temp]])] = i2;
		}
	}

	//////////////////// processing the end part of R
	for (MyUINT1 i1 = ((N - K - 1) / k1MULTI2) * k1MULTI2; i1 < N - K + 1; i1 += k1) {
        uint32_t h = hashFunc(start1 + i1);
	    if (!counts[h])
            continue;
		sampledPositions[cumm[h] + (--counts[h])] = i1;
	}
	delete[] counts;

    #pragma omp parallel for
	for (size_t i = 0; i < hashCount + 2; i++) {
	    if (sampledPositions[i] > N - K)
	        sampledPositions[i] = 0;
	}

	buffer = { sampledPositions, cumm };
}

template<typename MyUINT1, typename MyUINT2>
void SparseEMMatcher::processIgnoreCollisionsRef(HashBuffer<MyUINT1, MyUINT2>& buffer) {
    const unsigned int MULTI2 = 128;
    const unsigned int k1MULTI2 = k1 * MULTI2;

    const char* start1 = src.data();
    const size_t N = src.size();

    MyUINT1* sampledPositions = (processedSrc - SRC_SHIFT)?buffer.first:new MyUINT1[hash_size]();

    #pragma omp parallel for
    for (MyUINT1 i1 = processedSrc; i1 < N - K - k1MULTI2; i1 += k1MULTI2) {
        uint32_t hashPositions[MULTI2];
        const char* tempPointer = start1 + i1;
        for (unsigned int temp = 0; temp < MULTI2; ++temp) {
            hashPositions[temp] = hashFunc(tempPointer);
            tempPointer += k1;
            _prefetch((char*)(sampledPositions + hashPositions[temp]), 1);
        }

        MyUINT1 i2 = i1;
        for (size_t temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            sampledPositions[hashPositions[temp]] = i2;
        }
    }

    //////////////////// processing the end part of R
    MyUINT1 i1;
    for (i1 = SRC_SHIFT + ((N - K - 1) / k1MULTI2) * k1MULTI2; i1 < N - K + 1; i1 += k1) {
        sampledPositions[hashFunc(start1 + i1)] = i1;
    }
    processedSrc = i1;
    buffer.first = sampledPositions;
}

template <class MyUINT1, class MyUINT2>
void SparseEMMatcher::deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf) {
    delete[] buf.first;
    delete[] buf.second;
}

template<typename MyUINT1, typename MyUINT2>
void SparseEMMatcher::processExactMatchQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches,
                                                const char* start2, size_t N2, bool destIsSrc, bool revComplMatching,
                                                uint32_t minMatchLength){
    if (!collisionSupport) {
        processExactMatchQueryIgnoreCollisionsTight<MyUINT1, MyUINT2>(buffer, resMatches, start2, N2, destIsSrc,
                revComplMatching, minMatchLength);
        return;
    }

	const unsigned int MULTI = 256;
	const unsigned int k2MULTI = k2 * MULTI;

	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	std::uint32_t l1 = 0, l2 = 0, r1 = 0, r2 = 0;

	size_t charExtensions = 0ULL;

    const int skip = K / k1 - 1;
    const int skipK2 = skip * k2;
    const bool MULTI_MODE = true;
    *v1logger << "Minimal matching length = " << minMatchLength << "; ";
    *v1logger << "Skip factor = " << skip << "; ";
    *v1logger << "Multi-mode = " << (MULTI_MODE?"true":"false") << std::endl;

    const char* start1 = src.data() + SRC_SHIFT;
    const size_t N = src.size() - SRC_SHIFT;

    size_t i1 = 0;
    const char* end1 = start1 + N;
    const char* end2 = start2 + N2;
    if (MULTI_MODE) {
        for (i1 = 0; i1 + K + k2MULTI < N2 + 1; i1 += k2MULTI) {
            const char *curr2 = start2 + i1;
            size_t tempCount = 0;
            for (size_t i2 = 0; i2 < MULTI; ++i2) {
                hArray[tempCount++] = hashFunc(curr2);
                curr2 += k2;
            }
            for (size_t i2 = 0; i2 < tempCount; ++i2) {
                memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT2) * 2);
            }

            curr2 = start2 + i1;
            for (size_t i2 = 0; i2 < tempCount; ++i2) {
                if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
                    curr2 += k2;
                    continue;
                }

                if (curr2 - LK2 >= start2) memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
                if (curr2 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end2) memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

                for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
                    ++charExtensions;
                    const char *curr1 = start1 + sampledPositions[j];

                    uint64_t tmpMatchSrcPos = sampledPositions[j];
                    uint64_t tmpMatchDestPos = curr2 - start2;
                    if (destIsSrc && (revComplMatching ? N2 - tmpMatchSrcPos < tmpMatchDestPos
                                                       : curr2 - start2 >= tmpMatchSrcPos))
                        continue;
                    if (resMatches.size() > 0 &&
                        tmpMatchDestPos - tmpMatchSrcPos == resMatches.back().posDestText - resMatches.back().posSrcText
                        && tmpMatchDestPos + K < resMatches.back().posDestText + resMatches.back().length) {
                        curr2 += skipK2;
                        i2 += skip;
                        break;
                    }

                    if (curr1 - LK2 >= start1) memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
                    if (curr1 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end1) memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

                    if (r1 == r2 || l1 == l2) {
                        const char *p1 = curr1 + K - 1;
                        const char *p2 = curr2 + K - 1;
                        while (++p1 != end1 && ++p2 != end2 && *p1 == *p2);
                        const char *right = p1;
                        p1 = curr1;
                        p2 = curr2;
                        while (p1 != start1 && p2 != start2 && *--p1 == *--p2);

                        if (right - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                            resMatches.push_back(TextMatch(p1 + 1 - start1, right - p1 - 1, (p2 + 1 - start2)));

                            curr2 += skipK2;
                            i2 += skip;
                            break;
                        }
                    }
                }
                curr2 += k2;
            }
        }
    }
    //////////////////// processing the end part of Q  //////////////////////
    const char* curr2 = start2 + i1;
    for (; i1 + K < N2 + 1; i1 += k2) {
        memcpy(posArray, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);

        if (posArray[0] == posArray[1]) {
            curr2 += k2;
            continue;
        }

        if (curr2 - LK2 >= start2) memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
        if (curr2 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end2) memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

        for (MyUINT1 j = posArray[0]; j < posArray[1]; ++j) {
            ++charExtensions;
            const char* curr1 = start1 + sampledPositions[j];

            uint64_t tmpMatchSrcPos = sampledPositions[j];
            uint64_t tmpMatchDestPos = curr2 - start2;
            if (destIsSrc && (revComplMatching ? N2 - tmpMatchSrcPos < tmpMatchDestPos
                                               : curr2 - start2 >= tmpMatchSrcPos))
                continue;
            if (resMatches.size() > 0 &&
                tmpMatchDestPos - tmpMatchSrcPos == resMatches.back().posDestText - resMatches.back().posSrcText
                && tmpMatchDestPos + K < resMatches.back().posDestText + resMatches.back().length) {
                curr2 += skipK2;
                i1 += skipK2;
                break;
            }

            memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
            memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

            if (r1 == r2 || l1 == l2) {
                const char* p1 = curr1 + K - 1;
                const char* p2 = curr2 + K - 1;
                while (++p1 != end1 && ++p2 != end2 && *p1 == *p2);
                const char *right = p1;
                p1 = curr1;
                p2 = curr2;
                while (p1 != start1 && p2 != start2 && *--p1 == *--p2);

                if (right - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                    resMatches.push_back(TextMatch(p1 + 1 - start1, right - p1 - 1, (p2 + 1 - start2)));
                    curr2 += skipK2;
                    i1 += skipK2;
                    break;
                }
            }
        }
        curr2 += k2;
    }
    //////////////////// processing the end part of Q  //////////////////////

	*v1logger << "Character extensions = " << charExtensions <<  "\n";
}

template<typename MyUINT1, typename MyUINT2>
void SparseEMMatcher::processExactMatchQueryIgnoreCollisionsTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches,
        const char* start2, size_t N2, bool destIsSrc, bool revComplMatching, uint32_t minMatchLength){
    if (destIsSrc || revComplMatching) {
        fprintf(stderr, "Source as destination and reverse-complement matching modes unsupported\n\n");
        exit(EXIT_FAILURE);
    }

    MyUINT1* sampledPositions = buffer.first;
    MyUINT2 posArray;

    std::uint32_t l1 = 0, l2 = 0, r1 = 0, r2 = 0;

    size_t charExtensions = 0ULL;

    const char* start1 = src.data();
    const size_t N = src.size();

    size_t i2 = 0;
    const char* end1 = start1 + N;
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

        memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
        memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

        if (r1 == r2 || l1 == l2) {
            const char* p1 = curr1 + K;
            const char* p2 = curr2 + K;
            const int MEMCMP_BLOCK = 64;
            size_t blocks1 = (end1 - p1) / MEMCMP_BLOCK;
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
            while (++p1 != end1 && ++p2 != end2 && *p1 == *p2);
            const char *right = p1;
            p1 = curr1;
            p2 = curr2;
            while (p1 != start1 && p2 > guard2 - 1 && *--p1 == *--p2);
//            while (p2 < guard2 - 1) { p1++; p2++; } // to definitely remove overlapping matches

            if (right - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                resMatches.push_back(TextMatch(p1 + 1 - start1, right - p1 - 1, (p2 + 1 - start2)));
                const int skip = ((p2 - start2 + right - p1) - (curr2 - start2)) / k2 * k2;
                if (skip) {
                    curr2 += skip - k2;
                    i2 += skip - k2;
                }
                guard2 = p2 + (right - p1); // to avoid large overlaps
                continue;
            }
        }
    }
    //////////////////// processing the end part of Q  //////////////////////

//    *v1logger << "Character extensions = " << charExtensions <<  "\n";
}

using namespace std;

SparseEMMatcher::SparseEMMatcher(const size_t srcLengthLimit, const uint32_t targetMatchLength,
                                 bool collisionSupport, int _k1, int _k2, uint32_t minMatchLength)
    : maxSrcLength(srcLengthLimit), L(targetMatchLength), collisionSupport(collisionSupport) {
    src.reserve(this->maxSrcLength);
    src.push_back(0);
    initHashFuncMatrix();
    if (minMatchLength > targetMatchLength)
        minMatchLength = targetMatchLength;
    initParams(minMatchLength, _k1, _k2);
    displayParams();
    if (srcLengthLimit / k1 >= (1ULL << 32)) {
        bigRef = 2;  // huge Reference
        *v1logger  << "WARNING - LARGE reference file (SIZE / k1 > 4GB), 64-bit arrays\n";
    } else if (srcLengthLimit >= (1ULL << 32)) {
        bigRef = 1;  // large Reference
        *v1logger << "WARNING - BIG reference file (>4GB), 64-bit arrays\n";
    } else {
        bigRef = 0;  // small Reference
    }
}

void SparseEMMatcher::loadSrc(const char *srcText, size_t srcLength) {
    if (src.size() + srcLength > maxSrcLength)
        srcLength = maxSrcLength - src.size();

    src.append(srcText, srcLength);

    if (bigRef == 2) {
        processRef<std::uint64_t, std::uint64_t>(buffer2);
    } else if (bigRef == 1) {
        processRef<std::uint64_t, std::uint32_t>(buffer1);
    } else {
        processRef<uint32_t, uint32_t>(buffer0);
    }
}

SparseEMMatcher::~SparseEMMatcher() {
    if (bigRef == 2)
        deleteHashBuffer(buffer2);
    else if (bigRef == 1)
        deleteHashBuffer(buffer1);
    else if (bigRef == 0)
        deleteHashBuffer(buffer0);
}


void
SparseEMMatcher::matchTexts(vector <TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                          uint32_t minMatchLength) {
    matchTexts(resMatches, destText.data(), destText.size(), destIsSrc, revComplMatching, minMatchLength);
}

void
SparseEMMatcher::matchTexts(vector<TextMatch> &resMatches, const char* destText, size_t destLen, bool destIsSrc,
                bool revComplMatching, uint32_t minMatchLength) {
    if (minMatchLength < K) {
        fprintf(stderr, "Minimal matching length cannot be smaller than K (%d < %d)\n\n", minMatchLength, K);
        exit(EXIT_FAILURE);
    }
    resMatches.clear();
    if (bigRef == 2) {
        processExactMatchQueryTight<std::uint64_t, std::uint64_t>(buffer2, resMatches, destText, destLen, destIsSrc,
                                                                  revComplMatching, minMatchLength);
    } else if (bigRef == 1) {
        processExactMatchQueryTight<std::uint64_t, std::uint32_t>(buffer1, resMatches, destText, destLen, destIsSrc,
                                                                  revComplMatching, minMatchLength);
    }
    else {
        processExactMatchQueryTight<std::uint32_t, std::uint32_t>(buffer0, resMatches, destText, destLen, destIsSrc,
                                                                  revComplMatching, minMatchLength);
    }
}


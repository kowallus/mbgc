cmake_minimum_required(VERSION 3.5)
project(mbgc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2")

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(LIBDEFLATE_SOURCES
        coders/libdeflate/common_defs.h
        coders/libdeflate/libdeflate.h
        coders/libdeflate/lib/arm/cpu_features.c
        coders/libdeflate/lib/arm/cpu_features.h
        coders/libdeflate/lib/cpu_features_common.h
        coders/libdeflate/lib/deflate_constants.h
        coders/libdeflate/lib/lib_common.h
        coders/libdeflate/lib/utils.c
        coders/libdeflate/lib/x86/cpu_features.c
        coders/libdeflate/lib/x86/cpu_features.h
        coders/libdeflate/lib/arm/matchfinder_impl.h
        coders/libdeflate/lib/bt_matchfinder.h
        coders/libdeflate/lib/deflate_compress.c
        coders/libdeflate/lib/deflate_compress.h
        coders/libdeflate/lib/hc_matchfinder.h
        coders/libdeflate/lib/ht_matchfinder.h
        coders/libdeflate/lib/matchfinder_common.h
        coders/libdeflate/lib/x86/matchfinder_impl.h
        coders/libdeflate/lib/decompress_template.h
        coders/libdeflate/lib/deflate_decompress.c
        coders/libdeflate/lib/x86/decompress_impl.h
        coders/libdeflate/lib/arm/crc32_impl.h
        coders/libdeflate/lib/arm/crc32_pmull_helpers.h
        coders/libdeflate/lib/arm/crc32_pmull_wide.h
        coders/libdeflate/lib/crc32.c
        coders/libdeflate/lib/crc32_multipliers.h
        coders/libdeflate/lib/crc32_tables.h
        coders/libdeflate/lib/gzip_constants.h
        coders/libdeflate/lib/x86/crc32_impl.h
        coders/libdeflate/lib/x86/crc32_pclmul_template.h
        coders/libdeflate/lib/gzip_compress.c
        coders/libdeflate/lib/gzip_decompress.c)

set(HELPER_FILES
        utils/kseq.h utils/byteswap.h utils/helper.cpp utils/helper.h
        utils/input_with_libdeflate_wrapper.h utils/igzip_wrapper.h utils/input_with_libdeflate_wrapper.cpp
        coders/SymbolsPackingFacility.cpp coders/SymbolsPackingFacility.h
        coders/ContextAwareMismatchesCoder.cpp coders/ContextAwareMismatchesCoder.h)

set(LZMA_FILES
        coders/CodersLib.h coders/CodersLib.cpp
        coders/lzma/LzmaEnc.h coders/lzma/LzmaEnc.c
        coders/lzma/LzmaDec.h coders/lzma/LzmaDec.c
        coders/lzma/Alloc.h coders/lzma/Alloc.c
        coders/lzma/LzFindMt.h coders/lzma/LzFindMt.c
        coders/lzma/LzFind.h coders/lzma/LzFind.c
        coders/lzma/Ppmd7.h coders/lzma/Ppmd7.c
        coders/lzma/Ppmd7Enc.c coders/lzma/Ppmd7Dec.c
        coders/lzma/Threads.h coders/lzma/Threads.c
        coders/VarLenDNACoder.cpp coders/VarLenDNACoder.h
        coders/LzmaCoder.cpp coders/LzmaCoder.h
        coders/PpmdCoder.cpp coders/PpmdCoder.h
        coders/PropsLibrary.cpp coders/PropsLibrary.h)

set(HASHES_FILES
        utils/Hashes.h utils/hashes/xxhash.c utils/hashes/xxhash.h
        utils/hashes/metrohash64.cpp utils/hashes/metrohash64.h utils/hashes/platform.h
        utils/hashes/city.cpp utils/hashes/city.h utils/hashes/citycrc.h)

set(MATCHER_FILES
        ${HELPER_FILES}
        ${HASHES_FILES}
        matching/TextMatchers.cpp matching/TextMatchers.h
        matching/SlidingWindowSparseEMMatcher.cpp matching/SlidingWindowSparseEMMatcher.h
        matching/SimpleSequenceMatcher.cpp matching/SimpleSequenceMatcher.h
        matching/copmem/CopMEMMatcher.cpp matching/copmem/CopMEMMatcher.h)

set(MBGC_FILES
        ${MATCHER_FILES}
        ${LZMA_FILES}
        ${LIBDEFLATE_SOURCES}
        mbgccoder/MBGC_API.h
        mbgccoder/MBGC_Params.h
        mbgccoder/MBGC_Encoder.cpp mbgccoder/MBGC_Encoder.h
        mbgccoder/MBGC_Decoder.cpp mbgccoder/MBGC_Decoder.h)

add_executable(mbgc main.cpp ${MBGC_FILES})
add_executable(mbgc-dev main.cpp ${MBGC_FILES})
target_compile_definitions(mbgc-dev PUBLIC "-DDEVELOPER_BUILD")

if (CYGWIN OR MINGW)
    target_link_libraries(mbgc -static-libgcc)
    target_link_libraries(mbgc -static-libstdc++)
    target_link_libraries(mbgc-dev -static-libgcc)
    target_link_libraries(mbgc-dev -static-libstdc++)
    target_link_libraries(mbgc ${CMAKE_SOURCE_DIR}/libs/libacof64.lib)
    target_link_libraries(mbgc-dev ${CMAKE_SOURCE_DIR}/libs/libacof64.lib)
else()
    target_link_libraries(mbgc ${CMAKE_SOURCE_DIR}/libs/libaelf64.a)
    target_link_libraries(mbgc-dev ${CMAKE_SOURCE_DIR}/libs/libaelf64.a)
endif()

target_link_libraries(mbgc pthread)
target_link_libraries(mbgc-dev pthread)

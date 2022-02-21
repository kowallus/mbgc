cmake_minimum_required(VERSION 3.4)
project(mbgc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2")

include(CheckIncludeFileCXX)
CHECK_INCLUDE_FILE_CXX("libdeflate.h" HAVE_LIBDEFLATE)
IF(NOT HAVE_LIBDEFLATE)
    message( "Libdeflate is not found" )
ENDIF()

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(HELPER_FILES
        utils/kseq.h utils/byteswap.h utils/helper.cpp utils/helper.h
        utils/libdeflate_wrapper.h utils/igzip_wrapper.h utils/libdeflate_wrapper.cpp)

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
        coders/VarLenDNACoder.cpp coders/VarLenDNACoder.h coders/LzmaCoder.cpp coders/LzmaCoder.h coders/PpmdCoder.cpp coders/PpmdCoder.h coders/PropsLibrary.cpp coders/PropsLibrary.h)

set(HASHES_FILES
        utils/Hashes.h utils/hashes/xxhash.c utils/hashes/xxhash.h
        utils/hashes/metrohash64.cpp utils/hashes/metrohash64.h utils/hashes/platform.h
        utils/hashes/city.cpp utils/hashes/city.h utils/hashes/citycrc.h)

set(MATCHER_FILES
        ${HELPER_FILES}
        ${HASHES_FILES}
        matching/TextMatchers.cpp matching/TextMatchers.h
        matching/SlidingWindowSparseEMMatcher.cpp matching/SlidingWindowSparseEMMatcher.h)

set(MBGC_FILES
        ${MATCHER_FILES}
        ${LZMA_FILES}
        mbgccoder/MBGC_Params.h
        mbgccoder/MBGC_Encoder.cpp mbgccoder/MBGC_Encoder.h mbgccoder/MBGC_Decoder.cpp mbgccoder/MBGC_Decoder.h)

add_executable(mbgc main.cpp ${MBGC_FILES})
add_executable(mbgc-dev main.cpp ${MBGC_FILES})
target_compile_definitions(mbgc-dev PUBLIC "-DDEVELOPER_BUILD")

if (CYGWIN OR MINGW)
    target_link_libraries(mbgc ${CMAKE_SOURCE_DIR}/libs/libacof64.lib)
    target_link_libraries(mbgc-dev ${CMAKE_SOURCE_DIR}/libs/libacof64.lib)
else()
    target_link_libraries(mbgc ${CMAKE_SOURCE_DIR}/libs/libaelf64.a)
    target_link_libraries(mbgc-dev ${CMAKE_SOURCE_DIR}/libs/libaelf64.a)
endif()
if (HAVE_LIBDEFLATE)
    target_link_libraries(mbgc libdeflate.a)
    target_link_libraries(mbgc-dev libdeflate.a)
else()
    target_compile_definitions(mbgc PUBLIC "-DNO_GZ_SUPPORT")
    target_compile_definitions(mbgc-dev PUBLIC "-DNO_GZ_SUPPORT")
endif()
target_link_libraries(mbgc pthread)
target_link_libraries(mbgc-dev pthread)

#target_link_libraries( mbgc libisal.so )

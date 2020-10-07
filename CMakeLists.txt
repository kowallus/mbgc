cmake_minimum_required(VERSION 3.6)
project(fnacoder)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG -mavx2")

if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    message("Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(HELPER_FILES
        utils/kseq.h utils/byteswap.h utils/helper.cpp utils/helper.h
        utils/libdeflate_wrapper.h utils/igzip_wrapper.h)

set(HASHES_FILES
        utils/Hashes.h utils/hashes/xxhash.c utils/hashes/xxhash.h
        utils/hashes/metrohash64.cpp utils/hashes/metrohash64.h utils/hashes/platform.h
        utils/hashes/city.cpp utils/hashes/city.h utils/hashes/citycrc.h)

set(MATCHER_FILES
        ${HELPER_FILES}
        ${HASHES_FILES}
        matching/TextMatchers.cpp matching/TextMatchers.h
        matching/CopMEMMatcher.cpp matching/CopMEMMatcher.h )

add_executable(fnacoder main.cpp ${MATCHER_FILES})

target_link_libraries( fnacoder libdeflate.a )
#target_link_libraries( fnacoder libisal.so )

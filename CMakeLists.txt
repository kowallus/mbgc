cmake_minimum_required(VERSION 3.5)
project(mbgc)

set(CMAKE_CXX_STANDARD 17)

if(APPLE)
    message("Processing Apple build")
    if(DEFINED ENV{CONDA_PREFIX})
        set(LIBOMP_PATHS "$ENV{CONDA_PREFIX}/lib/libomp.dylib")
        if(NOT EXISTS ${LIBOMP_PATHS})
            message("libomp not found in ${LIBOMP_PATHS}")
            file(GLOB LIBOMP_PATHS $ENV{CONDA_PREFIX}/*/libomp.dylib)
        endif()
        set(CLANG_PATH $ENV{CONDA_PREFIX}/bin/clang)
        if(NOT EXISTS ${CLANG_PATH})
            message("clang not found in ${CLANG_PATH}")
            unset(CLANG_PATH)
        endif()
        set(CLANGXX_PATH $ENV{CONDA_PREFIX}/bin/clang++)
        if(NOT EXISTS ${CLANGXX_PATH})
            message("clang++ not found in ${CLANGXX_PATH}")
            unset(CLANGXX_PATH)
        endif()
    else()
        message("Conda not found")
        message("Proceeding with backup configuration")
        # Backup configuration may require setting proper paths in
        #     LIBOMP_PATHS AND CUSTOM_COMPILER_PATH variables.
        set(LIBOMP_PATHS "/usr/local/opt/libomp/lib/libomp.dylib")
        if(NOT EXISTS ${LIBOMP_PATHS})
            message("libomp not found in ${LIBOMP_PATHS}")
            unset(LIBOMP_PATHS)
            file(GLOB LIBOMP_PATHS /usr/local/opt/*/libomp.dylib)
            if(NOT LIBOMP_PATHS)
                file(GLOB LIBOMP_PATHS /opt/*/libomp.dylib)
            endif()
            if(NOT LIBOMP_PATHS)
                file(GLOB LIBOMP_PATHS ~/*/libomp.dylib)
            endif()
        endif()
        if((NOT CLANG_PATH) OR (NOT CLANGXX_PATH))
            set(LLVM_PATH "/usr/local/opt/llvm/bin")
            set(CUSTOM_COMPILER_PATH "/usr/local/opt/llvm@13/bin")
            find_program(CLANG_PATH clang ${CUSTOM_COMPILER_PATH} ${LLVM_PATH} NO_DEFAULT_PATH)
            find_program(CLANGXX_PATH clang++ ${CUSTOM_COMPILER_PATH} ${LLVM_PATH} NO_DEFAULT_PATH)
        endif()
    endif()
    if(LIBOMP_PATHS)
        list(GET LIBOMP_PATHS 0 LIBOMP_PATH)
        get_filename_component(LIBOMP_DIRECTORY "${LIBOMP_PATH}" DIRECTORY)
        message("Found LIBOMP: ${LIBOMP_DIRECTORY}")
    endif()
    if((NOT CLANG_PATH) OR (NOT CLANGXX_PATH))
        message("Could not find clang or clang++ in specified paths: ${LLVM_PATH}; ${LLVM16_PATH}; ${LLVM14_PATH}; ${BIN_PATH};")
    else()
        set(CMAKE_C_COMPILER ${CLANG_PATH})
        set(CMAKE_CXX_COMPILER ${CLANGXX_PATH})
    endif()
    message("compilers: ${CMAKE_C_COMPILER};${CMAKE_CXX_COMPILER}")
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.10" CACHE STRING "Minimum OS X deployment version")
    get_filename_component(CLANG_DIR "${CLANG_PATH}" DIRECTORY)
    get_filename_component(COMPILER_DIR "${CLANG_DIR}" DIRECTORY)
    message("${CLANG_PATH};${CLANG_DIR};${COMPILER_DIR}")
    set(OPENMP_LIBRARIES "${COMPILER_DIR}/lib")
    set(OPENMP_INCLUDES "${COMPILER_DIR}/include")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(OpenMP_C "${CMAKE_C_COMPILER}")
        set(OpenMP_C_FLAGS "-fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
    endif()
    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
        set(OpenMP_CXX_FLAGS "-fopenmp -Wno-unused-command-line-argument")
        set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
        set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
        set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
    endif()
    find_package(OpenMP REQUIRED)
    if (OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
    link_directories("${LIBOMP_DIRECTORY}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-format -Wl,-rpath ${LIBOMP_DIRECTORY} -L${LIBOMP_DIRECTORY}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|AMD64)")
    if(CYGWIN OR MINGW)
        set(ASM_LIB ${CMAKE_SOURCE_DIR}/libs/libacof64.lib)
    elseif(UNIX AND NOT APPLE)
        set(ASM_LIB ${CMAKE_SOURCE_DIR}/libs/libaelf64.a)
    elseif(APPLE)
        set(ASM_LIB ${CMAKE_SOURCE_DIR}/libs/libamac64.a)
    endif()
    link_libraries(${ASM_LIB})
endif()

link_libraries(pthread)

if(CYGWIN OR MINGW)
    message("Processing Windows build")
    link_libraries(-static-libgcc)
    link_libraries(-static-libstdc++)
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
add_executable(mbgc-dev EXCLUDE_FROM_ALL main.cpp ${MBGC_FILES})

target_compile_definitions(mbgc-dev PUBLIC -DDEVELOPER_BUILD)

if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm|aarch64)")
    add_executable(mbgc-noavx main.cpp ${MBGC_FILES})
    target_compile_options(mbgc PUBLIC -mavx2)
    target_compile_options(mbgc-dev PUBLIC -mavx2)
endif()


cmake_minimum_required(VERSION 3.7)
project(mtVariant VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_COMPILER "/share/app/gcc-7.4.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/share/app/gcc-7.4.0/bin/g++")

if (NOT CMAKE_BUILD_TYPE)
    set(DEFAULT_CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE ${DEFAULT_CMAKE_BUILD_TYPE} CACHE STRING
            "Choose the type of build (default: ${DEFAULT_CMAKE_BUILD_TYPE})" FORCE)
endif ()


find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(CURL REQUIRED)
find_library(htslib libhts.a "/share/FGI2017B/users/huangfei/software/bin/lib")
find_library(htslib hts "/share/FGI2017B/users/huangfei/software/bin/lib")

add_executable(mtVariant
        mtVariant.h mtVariant.cc
        Bam.h Bam.cc
        CoverageCounter.h CoverageCounter.cc
        EventScanner.h EventScanner.cc
        IndelEvent.h IndelEvent.cc
        ReferenceSequence.h ReferenceSequence.cc
        SnpEvent.h SnpEvent.cc
        MapStat.h MapStat.cc
        Primer.h Primer.cc
        VcfWriter.h VcfWriter.cc)

target_include_directories(mtVariant PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${LIBLZMA_INCLUDE_DIRS}
        ${CURL_INCLUDE_DIRS}
        "/share/FGI2017B/users/huangfei/software/bin/include"
        )

# Set static linking of gcc standard libraries to simplify binary distribution
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(STATIC_FLAGS -static-libgcc -static-libstdc++)
endif ()

target_link_libraries(mtVariant PUBLIC
        ${STATIC_FLAGS}
        ${htslib} -lpthread -lm  -lssl -lcrypto
        ${LIBLZMA_LIBRARIES}
        ${CURL_LIBRARIES}
        ZLIB::ZLIB
        BZip2::BZip2)


install(TARGETS mtVariant RUNTIME DESTINATION bin)

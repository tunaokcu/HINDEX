
cmake_minimum_required (VERSION 3.16)

include (CheckCXXCompilerFlag)

project (fusor CXX C)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build using static libraries" OFF)

find_library(JEMALLOC_LIB jemalloc HINTS ${CMAKE_CURRENT_SOURCE_DIR}/../jemalloc_install/lib)
if (JEMALLOC_LIB)
    set(CMAKE_CXX_STANDARD_LIBRARIES "-static-libgcc -static-libstdc++ ${JEMALLOC_LIB}" ${CMAKE_CXX_STANDARD_LIBRARIES})
else()
    message(FATAL_ERROR "jemalloc not found! You MUST run ./install_jemalloc.sh in the project root before running CMake.")
endif()

check_cxx_compiler_flag ("-fopenmp" HAS_OPENMP)
if (HAS_OPENMP)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif ()

check_cxx_compiler_flag ("-fopenmp-simd" HAS_OPENMP_SIMD)
if (HAS_OPENMP_SIMD)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp-simd -DSIMDE_ENABLE_OPENMP")
endif ()

# Production flags: full optimization, native arch, assertions disabled
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -DNDEBUG")
add_compile_options("-static-libstdc++")
add_compile_options("-static-libgcc")

# add all subdirectories

add_subdirectory(seqan)
add_subdirectory(hashing)
add_subdirectory(taxonomy)
add_subdirectory(hixf)
add_subdirectory(main)
add_subdirectory(ganon-classify)


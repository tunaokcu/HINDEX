cmake_minimum_required(VERSION 3.16)

message(STATUS "${CMAKE_GENERATOR}")

#if (WIN32)
#  if (MSVC)
#    set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
#                        -DMSVC_COMPILER")
#  else()
#    set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite
#                          -Dgtest_disable_pthreads=ON -Dgoogle_disable_pthreads=ON")
#  endif()
#else()
#  set(sdsl_cache_args  "-DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite")
#endif()

# Download seqan3 for SDSL dependency, but use local headers
configure_file(CMakeLists.txt.in ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download)
if(result)
  message(FATAL_ERROR "CMake step for seqan failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . --config Release
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/seqan3-download )

#add_library(sdsl STATIC IMPORTED GLOBAL)
#if (WIN32)
#  set_target_properties(sdsl PROPERTIES 
#    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/sdsl-lite-prefix/src/sdsl-lite-build/lib/Release/sdsl.lib")
#else()
#  set_target_properties(sdsl PROPERTIES 
#    IMPORTED_LOCATION "${CMAKE_CURRENT_BINARY_DIR}/seqan2-download/sdsl-lite-prefix/src/sdsl-lite-build/lib/libsdsl.a")
#endif()
#target_include_directories( sdsl SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/sdsl-lite/include )


#set(ZLIB_ROOT "${CMAKE_CURRENT_BINARY_DIR}/zlib")
#add_library(zlib STATIC IMPORTED GLOBAL)
#set_target_properties(zlib PROPERTIES IMPORTED_LOCATION "${ZLIB_ROOT}/lib/libz.a")
#set(ZLIB_LIBRARIES "${ZLIB_ROOT}/lib/libz.a")
#set(ZLIB_FOUND True)
#set(ZLIB_INCLUDE_DIRS "${ZLIB_ROOT}/include")

#set(BZIP2_ROOT "${CMAKE_CURRENT_BINARY_DIR}/bzip2")
#add_library(bzip2 STATIC IMPORTED GLOBAL)
#set_target_properties(bzip2 PROPERTIES IMPORTED_LOCATION "${BZIP2_ROOT}/lib/libbz2_static.a")
#set(BZIP2_LIBRARIES bzip2)
#set(BZIP2_INCLUDE_DIR "${BZIP2_ROOT}/include")

#find_library(zlib)

find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)


# Use local modified seqan3 copy instead of downloaded version
set( SEQAN3_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/external_libs )

# Set C++20 requirement (was done by find_package before)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Skip find_package - we're using local headers directly
# find_package (SeqAn3 3.3.0 REQUIRED HINTS ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/build_system)

add_library( seqan3 INTERFACE )
# Use our local copy for includes
target_include_directories( seqan3 SYSTEM INTERFACE ${SEQAN3_INCLUDE_PATH} )
# Add SDSL include path (needed by seqan3)
target_include_directories( seqan3 SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/submodules/sdsl-lite/include )
# Add Cereal include path (needed by seqan3 and chopper)
target_include_directories( seqan3 SYSTEM INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/submodules/cereal/include )
# Create alias to match the namespaced target that other CMakeLists expect
add_library( seqan3::seqan3 ALIAS seqan3 )
# Link ZLIB and BZip2 (needed by seqan3)
target_link_libraries( seqan3 INTERFACE ZLIB::ZLIB BZip2::BZip2 )
# Define that ZLIB and BZIP2 are available (needed by seqan3 headers)
target_compile_definitions( seqan3 INTERFACE SEQAN3_HAS_ZLIB=1 SEQAN3_HAS_BZIP2=1 )


# No SEQAN3_DEFINITIONS needed
# add_compile_options( ${SEQAN3_DEFINITIONS} )
# Don't use SEQAN3_CXX_FLAGS since we don't have them
# set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN3_CXX_FLAGS}" )


#list(APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_BINARY_DIR}/seqan3/build_system")
#set(seqan3_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/build_system)
#find_package (seqan3 3.0 REQUIRED)
#set(CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)
#set(SEQAN_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/seqan3/seqan3-src/include)
#include_directories (${SEQAN_INCLUDE_PATH})
#set(SeqAn_DIR ${CMAKE_CURRENT_BINARY_DIR}/seqan2/seqan2-src/util/cmake)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

#find_package( SeqAn CONFIG REQUIRED PATHS ${SeqAn_DIR} )

#add_library( seqan INTERFACE )
#target_include_directories( seqan SYSTEM INTERFACE ${SEQAN_INCLUDE_DIRS} )
#target_link_libraries( seqan INTERFACE ${SEQAN_LIBRARIES} )

# Add include directories.

#message(FATAL_ERROR ${SEQAN_INCLUDE_DIRS})
#include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
#add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")


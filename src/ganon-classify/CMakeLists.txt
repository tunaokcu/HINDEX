file(GLOB GanonClassify_SRC "*.cpp")
file(GLOB GanonClassify_HDRS "include/ganon-classify/*.hpp")
add_executable (ganon-classify ${GanonClassify_HDRS} ${GanonClassify_SRC})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include"
                    "${CMAKE_CURRENT_SOURCE_DIR}/include/ganon-classify"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../ganon-utils"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../hashing"
                    "${CMAKE_CURRENT_SOURCE_DIR}/../taxonomy"
                    "${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/lib/simde/simde")

target_include_directories(ganon-classify INTERFACE 
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${CMAKE_CURRENT_BINARY_DIR}/../hashing
    ${CMAKE_CURRENT_BINARY_DIR}/../taxonomy
    ${CMAKE_CURRENT_BINARY_DIR}/../hixf
    ${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/lib/simde/simde
    )

if (WIN32)
    target_link_libraries (ganon-classify seqan3::seqan3)
else()
    target_link_libraries (ganon-classify seqan3::seqan3 hashing taxonomy hixf)
endif()

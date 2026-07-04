file(GLOB Main_SRC "*.cpp" "*.c")
list(REMOVE_ITEM Main_SRC "${CMAKE_CURRENT_SOURCE_DIR}/verify_fusor.cpp")
list(REMOVE_ITEM Main_SRC "${CMAKE_CURRENT_SOURCE_DIR}/build_hibff.cpp")
file(GLOB MAIN_HDRS "*.hpp" "*.h")
add_executable (${PROJECT_NAME} ${MAIN_HDRS} ${Main_SRC})
add_executable (verify_fusor verify_fusor.cpp)
target_link_libraries(verify_fusor hashing taxonomy hixf seqan3::seqan3)

add_executable(build_hibff build_hibff.cpp)
#find_package (Threads)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../hashing"
					"${CMAKE_CURRENT_SOURCE_DIR}/../taxonomy"
					"${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/lib/simde/simde")

target_include_directories(${PROJECT_NAME} INTERFACE 
	${CMAKE_CURRENT_BINARY_DIR}/../hashing
	${CMAKE_CURRENT_BINARY_DIR}/../taxonomy
	${CMAKE_CURRENT_BINARY_DIR}/../hixf
	${CMAKE_CURRENT_BINARY_DIR}/../seqan/chopper/lib/simde/simde
	)

# using seqan3
if (WIN32)
	#add_definitions(-DMSVC_COMPILER)
	target_link_libraries (${PROJECT_NAME} seqan3::seqan3)
else()
	target_link_libraries (${PROJECT_NAME} seqan3::seqan3 hashing taxonomy hixf)
endif()

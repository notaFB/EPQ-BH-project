# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name
project(bssn_simulation)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler flags for handling large objects
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fPIC -Wall -Wextra -Wpedantic -Wno-unknown-pragmas -march=native")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /O2 /EHsc")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
endif()

# Set the output directories for build artifacts
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Automatically gather all source files in src/ directory
file(GLOB_RECURSE SRC_FILES ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Create the executable
add_executable(bssn_simulation ${SRC_FILES})

# Link libraries if any (example: -lm)
# target_link_libraries(bssn_simulation <library_name>)

# Unit tests
file(GLOB_RECURSE TEST_FILES ${CMAKE_SOURCE_DIR}/test/*.cpp)

foreach(TEST_FILE ${TEST_FILES})
    # Extract test name from file name
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_FILE} ${SRC_FILES})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()

# Set build type to Release by default
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

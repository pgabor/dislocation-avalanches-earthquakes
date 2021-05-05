cmake_minimum_required(VERSION 3.5.1)
project(sdddstEQ)

# Check for build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Setting flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -g")
if (CMAKE_BUILD_TYPE EQUAL Release)
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Handling dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
message(STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH})

find_package("UMFPACK" MODULE REQUIRED)

include_directories("include")
include_directories(${UMFPACK_INCLUDES})

# Start of compilation
add_subdirectory(src)

cmake_minimum_required(VERSION 3.10)
project(Gatete)

set(CMAKE_BUILD_TYPE Release)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(GSL REQUIRED)
find_package(Threads REQUIRED)

# Only set up to -O2 in release mode as -O3 may break
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG"
    CACHE STRING "Flags used by the compiler during release builds." FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG"
    CACHE STRING "Flags used by the compiler during release builds." FORCE)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -ggdb -Wall -Wextra -pedantic -fno-omit-frame-pointer")

include(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
endif()

add_definitions(-DUSE_THREADS)
add_definitions(-DUSE_THREAD_POOL)

add_definitions(-DLEPTON_STEADY_STATE=0)
add_definitions(-DPROTON_STEADY_STATE=0)

add_subdirectory(src)
add_subdirectory(tables)
add_subdirectory(examples)

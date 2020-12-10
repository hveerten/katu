cmake_minimum_required(VERSION 3.2)
project(Katu)

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

if(${CMAKE_VERSION} VERSION_LESS "3.9.0")
    set(IPO_SUPPORTED 0)
else()
    cmake_policy(SET CMP0069 NEW)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED)
endif()

macro(prepare_executable target source)
    add_executable(${target} ${source})
    target_link_libraries(${target} m gsl gslcblas pthread model_static)
    target_include_directories(${target} PUBLIC ${MODEL_STATIC_INCLUDE_DIR})
    if(IPO_SUPPORTED)
        set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
    endif()
endmacro(prepare_executable)

add_definitions(-DUSE_THREADS)
add_definitions(-DUSE_THREAD_POOL)

add_definitions(-DLEPTON_STEADY_STATE=0)
add_definitions(-DPROTON_STEADY_STATE=0)

add_subdirectory(src)
add_subdirectory(tables)
add_subdirectory(examples)

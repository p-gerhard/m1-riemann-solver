cmake_minimum_required(VERSION 3.5)
project(m1-riemann-solver VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()

# Build the library
add_library(m1_riemann STATIC src/m1_riemann.c)

target_include_directories(m1_riemann
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_compile_options(m1_riemann PRIVATE $<$<C_COMPILER_ID:GNU>:-Wall -O3>)

find_package(GSL REQUIRED)
target_link_libraries(m1_riemann GSL::gsl)

add_subdirectory(test)
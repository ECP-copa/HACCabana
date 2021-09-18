cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX)

# Use C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Cabana)

option(ENABLE_SERIAL OFF)
option(ENABLE_GPU ON)

if (ENABLE_SERIAL)
    add_definitions(-DUSE_SERIAL_CPU=ON)
endif(ENABLE_SERIAL)

if (ENABLE_GPU)
    add_definitions(-DUSE_GPU=ON)
    find_package(CUDA)
    if (CUDA_FOUND)
        message("Found CUDA")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --expt-extended-lambda -g")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xcudafe --diag_suppress=esa_on_defaulted_function_ignored") # supress express cuda warnings
    else()
        message( FATAL_ERROR "Requested CUDA but cant find it!")
    endif(CUDA_FOUND)
endif(ENABLE_GPU)

INCLUDE_DIRECTORIES(./src)
file(GLOB SOURCES "src/*.cxx")

add_executable(driver_short-range ${SOURCES})
target_link_libraries(driver_short-range PUBLIC Cabana::cabanacore dl)

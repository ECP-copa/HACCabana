cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX)

# Use C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Cabana)

file(GLOB SOURCES "src/*.cxx")

add_executable(driver_short-range ${SOURCES})
target_link_libraries(driver_short-range PUBLIC Cabana::cabanacore dl)

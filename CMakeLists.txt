cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX C)

include(GNUInstallDirs)

find_package(Cabana)

add_subdirectory(src)

cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX C)

include(GNUInstallDirs)

find_package(Cabana REQUIRED VERSION 0.6.0)

add_subdirectory(src)

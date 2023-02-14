cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX C)

include(GNUInstallDirs)

find_package(Cabana)

option(HACCabana_REQUIRE_ARBORX "Build HACCabana with Cabana-ArborX neighbor list support" ${Cabana_ENABLE_ARBORX})
set(HACCabana_ENABLE_ARBORX ${HACCabana_REQUIRE_ARBORX})

add_subdirectory(src)

cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX C)

include(GNUInstallDirs)

find_package(Cabana)

option(HACCabana_REQUIRE_ARBORX "Build HACCabana with Cabana-ArborX neighbor list support" ${Cabana_ENABLE_ARBORX})
if(${HACCabana_REQUIRE_ARBORX} AND NOT ${Cabana_ENABLE_ARBORX})
  message(FATAL_ERROR "Cabana must be compiled with ArborX support to use HACCabana_REQUIRE_ARBORX")
endif()
set(HACCabana_ENABLE_ARBORX ${HACCabana_REQUIRE_ARBORX})

add_subdirectory(src)

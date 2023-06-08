cmake_minimum_required(VERSION 3.9)

project(HACCabana LANGUAGES CXX C)

include(GNUInstallDirs)

find_package(Cabana REQUIRED)
if(Kokkos_VERSION_MAJOR LESS 4 AND Kokkos_VERSION_MINOR LESS 7)
  message( FATAL_ERROR "HACCabana requires Kokkos >= 3.7" )
endif()

add_subdirectory(src)

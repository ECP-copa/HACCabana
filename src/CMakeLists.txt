file(GLOB HACCABANA_HEADERS *.h)
file(GLOB HACCABANA_SOURCES *.cxx)

install(FILES ${HACCABANA_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(HACCabana ${HACCABANA_SOURCES})

target_include_directories(HACCabana PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

install(TARGETS HACCabana
  EXPORT HACCabana_Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

target_link_libraries(HACCabana Cabana::Core)

install(TARGETS HACCabana DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Require minimum of C++14
set_target_properties(HACCabana PROPERTIES INTERFACE_COMPILE_FEATURES cxx_std_14)

add_executable(driver_short-range driver_gpu.cxx)
target_link_libraries(driver_short-range PUBLIC HACCabana)
install(TARGETS driver_short-range DESTINATION ${CMAKE_INSTALL_BINDIR})

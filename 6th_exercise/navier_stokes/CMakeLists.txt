cmake_minimum_required(VERSION 2.8)

PROJECT(navier_stokes)

set(THIS_DIR ${CMAKE_CURRENT_LIST_DIR})
message(${THIS_DIR})
set(SRC_DIR ${THIS_DIR}/src)
set(INCLUDE_DIR ${THIS_DIR}/include)
set(LIB_DIR ${THIS_DIR}/lib)
set(VTK_DIR D:/source/Libs/VTK-8.1.2/build)

set(SOURCE
    ${SRC_DIR}/main.cpp
    ${SRC_DIR}/navier_stokes.c
)
message(SOURCE: ${SOURCE})

set(HEADERS
    ${INCLUDE_DIR}/navier_stokes.h
    ${INCLUDE_DIR}/fftw3.h
    ${INCLUDE_DIR}/viridis.h
)
message(HEADERS: ${HEADERS})

find_package(VTK 8.1.2 REQUIRED)
include(${VTK_USE_FILE})
include_directories(INCLUDE_DIR)

add_library(LIBFFTW SHARED IMPORTED)
set_target_properties(LIBFFTW
    PROPERTIES
    IMPORTED_LOCATION ${LIBFFTW_DIR}/libfftw3-3.lib
)

add_executable(main ${SOURCE} ${HEADERS})
target_link_libraries(main ${LIBFFTW} ${VTK_LIBRARIES})

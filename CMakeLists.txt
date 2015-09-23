cmake_minimum_required(VERSION 2.8)
project(ingest_cm1 Fortran)

include (GNUInstallDirs)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  if (APPLE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -pipe -O2")
  else()
    set (CMAKE_Fortran_FLAGS_RELEASE "-std=f2008 -march=native -pipe -O2")
  endif()
  set (CMAKE_Fortran_FLAGS_DEBUG   "-std=f2008 -O0 -ggdb -g3 -fbacktrace -Wall -Wextra
  -pedantic -fsanitize=undefined -fsanitize=address -fsanitize=leak -fcheck=all")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-stand f08 -xhost -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-stand f08 -O0 -g")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # default
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

#set (CMAKE_Fortran_MODDIR_FLAG ${CMAKE_CURRENT_BINARY_DIR}/include)

if (WITH_HDF5)
  find_package(HDF5 COMPONENTS Fortran Fortran_HL REQUIRED)
  include_directories(${HDF5_INCLUDE_DIR})
  set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
  add_definitions(-DHDF5)
endif (WITH_HDF5)

if (WITH_HDF5_CMAKE)
  find_package(HDF5 NAMES hdf5 REQUIRED)
  include_directories(${HDF5_INCLUDE_DIR})
  set (LINK_LIBS ${LINK_LIBS} ${HDF5_LIBRARIES})
endif (WITH_HDF5_CMAKE)



add_subdirectory(src)
add_subdirectory(bin)

#http://www.cmake.org/Wiki/CMakeForFortranExample
#http://stackoverflow.com/questions/2186110/cmake-tutorial

message ("Installation prefix:       ${CMAKE_INSTALL_PREFIX}")
message ("Library installation path: ${CMAKE_INSTALL_PREFIX}/lib")
message ("Module file install path:  ${CMAKE_INSTALL_PREFIX}/include")

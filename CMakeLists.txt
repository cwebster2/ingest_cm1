cmake_minimum_required(VERSION 2.8)
project(ingest_cm1 Fortran)

include (GNUInstallDirs)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
#get_filename_component (default_prefix "." ABSOLUTE)
#set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#      "Choose the installation directory; by default it installs in the ingest-cm1 directory."
#      FORCE)


#set (INGESTCM1_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  if (APPLE)
    set (CMAKE_Fortran_FLAGS_RELEASE "-pipe -O2")
  else()
    set (CMAKE_Fortran_FLAGS_RELEASE "-march=native -pipe -O2")
  endif()
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -ggdb")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-xhost -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
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

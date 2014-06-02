#MAKEFILE

# SET your compiler and any compile / link options here:
#
# GNU Fortran build options.  Requires gcc >= 4.8
FC=gfortran
CC=gcc
FFLAGS= -std=f2008 -Wall -march=native -O2 -fimplicit-none -I../include

# Intel Fortran build options.
#FC=ifort
#CC=icc
#FFLAGS= -stand f08 -xHost -assume byterecl
##FFLAGS= -assume byterecl -g -debug full
##-diag-enable sc3

HDFINC=-I/usr/include
HDFLIB=-lhdf5_fortran -lhdf5hl_fortran -lhdf5 

#HDFINC=-I/twist/s2/ccw5079/hdf5/include
#HDFLIB=-L/twist/s2/ccw5079/hdf5/lib -lhdf5_fortran -lhdf5hl_fortran -lhdf5 -Wl,-rpath,/twist/s2/ccw5079/hdf5/lib

# You should not need to edit anything below this line.
# -----------------------------------------------------
#  ----------------------------------------------------
#  ----------------------------------------------------

INGESTOBJS =  cm1_base.o cm1_grads.o cm1_grads_single.o cm1_grads_mpi.o cm1_hdf5.o ingest_cm1.o
INGESTLIB = ../lib/ingest_cm1.a

all: test 

test: test_cm1_hdf 

test_cm1: $(INGESTLIB) test_cm1.o
	$(FC) $(FFLAGS) -o test_cm1 test_cm1.o $(INGESTLIB) $(HDFLIB)

test_cm1_mpi: $(INGESTLIB) test_cm1_mpi.o
	$(FC) test_cm1_mpi.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_mpi

test_cm1_hdf: $(INGESTLIB) test_cm1_hdf.o
	$(FC) test_cm1_hdf.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_hdf

test_cm1_gsingle: $(INGESTLIB) test_cm1_gsingle.o
	$(FC) test_cm1_gsingle.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_gsingle

#doc:
#	doxygen ../Doxygen.conf

%.o: %.f
	$(FC) $(FFLAGS) $(HDFINC) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(HDFINC) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<
clean: 
	rm -f $(TESTOBJS) test_cm1.o test_cm1_mpi.o test_cm1_hdf.o test_cm1 test_cm1_mpi test_cm1_hdf *.mod

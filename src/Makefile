# SET your compiler and any compile / link options here:

# GNU compiler --------------------------------------------
# GNU Fortran build options.  Requires gcc >= 4.8
FC=gfortran
CC=gcc
FFLAGS= -std=f2008 -Wall -march=native -O2 -fimplicit-none -fcheck-array-temporaries

# Intel compiler -----------------------------------------
# Intel Fortran build options.
#FC=ifort
#CC=icc
#FFLAGS= -stand f08 -xHost -assume byterecl

# HDF5 Support -------------------------------------------
# Uncomment the following lines to enable HDF5 support
# Make sure that the path in HDFINC points to your HDF5 install
HDFFLAGS=-DHDF5
HDFINC=-I/usr/include
HDFLIB=-lhdf5_fortran -lhdf5hl_fortran -lhdf5 


# You should not need to edit anything below this line.
# -----------------------------------------------------
#  ----------------------------------------------------
#  ----------------------------------------------------

INGESTOBJS =  cm1_base.o cm1_grads.o cm1_grads_single.o cm1_grads_mpi.o cm1_hdf5.o ingest_cm1.o
INGESTLIB = ingest_cm1.a

all: $(INGESTLIB) test install

test: test_cm1 test_cm1_mpi test_cm1_gsingle

test_cm1: $(INGESTLIB) test_cm1.o
	$(FC) $(FFLAGS) -o test_cm1 test_cm1.o $(INGESTLIB) $(HDFLIB)

test_cm1_mpi: $(INGESTLIB) test_cm1_mpi.o
	$(FC) test_cm1_mpi.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_mpi

test_cm1_gsingle: $(INGESTLIB) test_cm1_gsingle.o
	$(FC) test_cm1_gsingle.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_gsingle

ingest_cm1_lib: $(INGESTLIB)

$(INGESTLIB): $(INGESTOBJS)
	ar rcs $(INGESTLIB) $(INGESTOBJS)

install: $(INGESTLIB) test_cm1
	mkdir -p ../lib ../include
	cp $(INGESTLIB) ../lib
	cp *.mod ../include

tags: ingest_cm1.f90
#	ctags -R --fortran-kinds=+i
	etags -R *.f90


doc:
	doxygen ../Doxygen.conf

%.o: %.F90
	$(FC) -cpp $(FFLAGS) ${HDFFLAGS} $(HDFINC) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(HDFINC) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<
clean: 
	rm -f $(INGESTOBJS) $(TESTOBJS) $(INGESTLIB) test_cm1.o test_cm1_mpi.o test_cm1_hdf.o test_cm1 test_cm1_mpi test_cm1_hdf *.mod

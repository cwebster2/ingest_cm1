#MAKEFILE

# SET your compiler and any compile / link options here:
#
FC=gfortran
CC=gcc
FFLAGS= -std=f2008 -Wall -march=native -O2 -fimplicit-none
#FFLAGS= -assume byterecl -g -debug full
#-diag-enable sc3

HDFINC=-I/usr/include
HDFLIB=-lhdf5_fortran -lhdf5hl_fortran -lhdf5 

# You should not need to edit anything below this line.
# -----------------------------------------------------
#  ----------------------------------------------------
#  ----------------------------------------------------

INGESTOBJS =  cm1_base.o cm1_grads.o cm1_grads_mpi.o cm1_hdf5.o ingest_cm1.o
INGESTLIB = ingest_cm1.a

all: $(INGESTLIB) test install

test: test_cm1 test_cm1_mpi test_cm1_hdf

test_cm1: $(INGESTLIB) test_cm1.o
	$(FC) $(FFLAGS) -o test_cm1 test_cm1.o $(INGESTLIB) $(HDFLIB)

test_cm1_mpi: $(INGESTLIB) test_cm1_mpi.o
	$(FC) test_cm1_mpi.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_mpi

test_cm1_hdf: $(INGESTLIB) test_cm1_hdf.o
	$(FC) test_cm1_hdf.o $(FFLAGS) $(INGESTLIB) $(HDFLIB) -o test_cm1_hdf

ingest_cm1_lib: $(INGESTLIB)

$(INGESTLIB): $(INGESTOBJS)
	ar rcs $(INGESTLIB) $(INGESTOBJS)

install: $(INGESTLIB) test_cm1
	mkdir -p ../lib ../include
	cp $(INGESTLIB) ../lib
	cp *.mod ../include

tags: ingest_cm1.f90
	ctags -R --fortran-kinds=+i

#doc:
#	doxygen ../Doxygen.conf

%.o: %.f
	$(FC) $(FFLAGS) $(HDFINC) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(HDFINC) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<
clean: 
	rm -f $(INGESTOBJS) $(TESTOBJS) $(INGESTLIB) test_cm1.o test_cm1_mpi.o test_cm1_hdf.o test_cm1 test_cm1_mpi test_cm1_hdf *.mod

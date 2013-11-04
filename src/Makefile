#MAKEFILE

# SET your compiler and any compile / link options here:
#
FC=ifort
CC=icc
FFLAGS= -warn -stand 08 -traceback -assume byterecl -xHost 
#FFLAGS= -assume byterecl -g -debug full
#-diag-enable sc3

# You should not need to edit anything below this line.
# -----------------------------------------------------
#  ----------------------------------------------------
#  ----------------------------------------------------

INGESTOBJS = ingest_cm1.o 
INGESTLIB = ingest_cm1.a

all: $(INGESTLIB) test install

test: test_cm1 test_cm1_mpi

test_cm1: $(INGESTLIB) test_cm1.o
	$(FC) $(FFLAGS) -o test_cm1 test_cm1.o $(INGESTLIB)

test_cm1_mpi: $(INGESTLIB) test_cm1_mpi.o
	$(FC) test_cm1_mpi.o $(FFLAGS) $(INGESTLIB) -o test_cm1_mpi

ingest_cm1_lib: $(INGESTLIB)

$(INGESTLIB): $(INGESTOBJS)
	ar rcs $(INGESTLIB) $(INGESTOBJS)

install: $(INGESTLIB) test_cm1
	mkdir -p ../lib ../include
	cp $(INGESTLIB) ../lib
	cp ingest_cm1.mod ../include

tags: ingest_cm1.f90
	ctags -R

#doc:
#	doxygen ../Doxygen.conf

%.o: %.f
	$(FC) $(FFLAGS) $(OUTPUTINC) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) $(OUTPUTINC) -c $<

%.o: %.c
	$(CC) $(CFLAGS) -c $<
clean: 
	rm -f $(INGESTOBJS) $(TESTOBJS) $(INGESTLIB) test_cm1 test_cm1_mpi *.mod

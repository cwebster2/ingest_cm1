# Ingest_cm1 #

## README ##

Ingest_cm1 is a Fortran module to abstract access to model output 
files from George Bryan's CM1 cloud model.  

The current version of this module supports access to GrADS style files
with one output file per timestep or one output file per timestep per
MPI node; HDF5 output. 

These represent files created by the following options in the
CM1 namelist.input file (GRADS, GRADSMPI):

    output_format = 1
    output_filetpe = 2 or 3

And for HDF data:

    output_format = 3, 4, 5
    output_filetype = 3

NOTE: Support for the native tiled output is not yet implemented.  You must untile the output with the included python script.

Future output format support is expected for NetCDF4 and native MPI tiled HDF output.

## GETTING THE SOFTWARE ##

    $ git clone https://github.com/cwebster2/ingest_cm1.git

## INSTALLING ##

A Fortran 2003 compiler is required to build this software.  If using GNU `gfortran`, 
use of version 4.8 or later is required due to the use of allocatable arrays of 
polymorphic types.  Intel `ifort` should work but is untested, as are other fortran compilers.

    $ cd ingest_cm1/src; make

## USING ##

The cm1_dataset type is the public documented API, but feel free to use the 
particular backend classes directly if so included. 

### open_dataset ###

```fortran
    integer function open_dataset(self, dsetpath, dsetbasename, dsettype, grids, nodex, nodey, hdfmetadatatime)
      class(cm1_dataset) :: self
      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in) :: dsetbasename
      integer, intent(in)          :: dsettype
      character, dimension(:)      :: grids
      integer, optional :: nodex, nodey, hdfmetadatatime
```

### close_dataset ###
### read_3d and read_2d ###
### get_nx, get_ny, get_nz ###
### get_x, get_y, get_z ###

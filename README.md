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
integer function open_dataset(self, dsetpath, dsetbasename, dsettype, grids, nodex, nodey)
   class(cm1_dataset) :: self
   character(len=*), intent(in) :: dsetpath
   character(len=*), intent(in) :: dsetbasename
   integer, intent(in)          :: dsettype
   character, dimension(:)      :: grids
   integer, optional :: nodex, nodey
```
This opens the dataset located at `dsetpath` with basename `dsetbasename`.  
`dsettype` is one of `GRADS`, `GRADSMPI` or `HDF`.
`grids` is an array of grids to open, e.g. ['s', 'u', 'v', 'w'] for a full dataset or ['s'] if only the scalar grid is desired.
The variables `nodex` and `nodey` are only used for the `GRADSMPI` dsettype.  These are the same values used in the namelist.input for the MPI run.

### close_dataset ###

```fortran
integer function close_dataset(self)
   class(cm1_dataset) :: self
```
This closes the dataset.

### read_3d and read_2d ###

```fortran
integer function read_3d(self, time, grid, varname, Field3D)
   implicit none
   class(cm1_dataset) :: self
   integer            :: time, gridno
   character          :: grid
   character(len=*)   :: varname
   real, dimension(:,:,:) :: Field3D
```
These read a 2/3D variable `varname` from grid `grid` at time `time`.

### get_nx, get_ny, get_nz ###

```fortran
integer function get_nx(self, grid)
   implicit none
   class(cm1_dataset) :: self
   character          :: grid
```
These get dimensions of the specified grid `grid`.

### get_x, get_y, get_z ###

```fortran
integer function get_x(self, grid, x)
   implicit none
   class(cm1_dataset) :: self
   character          :: grid
   integer            :: gridno
   real, dimension(:) :: x
```
These get the mesh of grid `grid` along the specified dimension

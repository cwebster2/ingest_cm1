# Ingest_cm1 #

[![Build Status](https://travis-ci.org/cwebster2/ingest_cm1.svg?branch=master)](https://travis-ci.org/cwebster2/ingest_cm1)



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


* [Getting](#getting-the-software)
* [Installing](#installing)
* [Using](#using)
  - [The cm1_dataset front-end](#The-cm1_dataset-front-end)
    + [Open_dataset](#open_dataset)
    + [Close_dataset](#close_dataset)
    + [Reading data](#reading-data)
    + [Accessing grid information](#accessing_grid_information)
  - The backend classes
* Contributing
* Author
* License


## GETTING THE SOFTWARE ##

Download a current snapshot of ingest_cm1 with:

    git clone https://github.com/cwebster2/ingest_cm1.git


## INSTALLING ##

A Fortran 2003 compiler is required to build this software.  If using GNU `gfortran`, 
use of version 4.8 or later is required due to the use of allocatable arrays of 
polymorphic types.  Intel `ifort` should work but is untested, as are other fortran compilers.

IngestCM1 uses cmake to generate makefiles and you have a few options when building.  The default
options are to build a static library without HDF5 support and install into ingest_cm1's directory.
To use this default, from the ingest_cm1 directory, do:

    cmake .

To enable HDF5 support, use:

    cmake . -DWITH_HDF5=1

If your HDF5 library is not found, you can specify a search path with the `HDF5_ROOT` environment variable.

    HDF5_ROOT="/path/to/hdf5" cmake . -DWITH_HDF5=1

If your HDF5 library is built with cmake support, you can try this instead:

    cmake . -DWITH_HDF5_CMAKE

To build a shared library instead of a static library, add the flag `-DBUILD_SHARED_LIBS=1` flag to cmake and
to change the installation location, e.g. /usr/local, use `-DCMAKE_INSTALL_PREFIX:PATH=/usr/local`.  An example
of all of the above is:

    HDF5_ROOT="/path/to/hdf5" cmake . -DWITH_HDF5=1 -DBUILD_SHARED_LIBS=1 -DCMAKE_INSTALL_PREFIX:PATH=/usr

which will generate makefiles to build a shared library, install into /usr/lib and /usr/include and
include HDF5 support.

To build and install ingest_cm1 do:

    make
    make install

## USING ##

### The cm1_dataset front end ###

The cm1_dataset front-end is able to load a dataset with multiple grids spread across multiple files.
This paradigm assumes the output style of CM1 where there are 4 grids each for scalar, u, v, and w
variables. Within each grid the output may be a single file per timestep, a single file for all timesteps
or a file for each MPI rank at each timestep.  This is defined by the `dsettype` variable below.

#### open_dataset ####

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
`dsettype` is one of `GRADS`, `GRADSMPI`, `GRADSSINGLE`, or `HDF`.
`grids` is an array of grids to open, e.g. `['s', 'u', 'v', 'w']` for a full dataset or `['s']` if only the scalar grid is desired.
The variables `nodex` and `nodey` are only used for the `GRADSMPI` dsettype.  These are the same values used in the namelist.input for the MPI run.

#### close_dataset ####

```fortran
integer function close_dataset(self)
   class(cm1_dataset) :: self
```
This closes the dataset.

#### Reading data ####

You can read a 3d grid using the `read_3d` function:

```fortran
integer function read_3d(self, time, grid, varname, Field3D)
   implicit none
   class(cm1_dataset) :: self
   integer            :: time, gridno
   character          :: grid
   character(len=*)   :: varname
   real, dimension(:,:,:) :: Field3D
```
These read a 2/3D variable `varname` from grid `grid` at time `time`.  To read 2D variables,
use the similar function `read_2d`.

#### Accessing grid information ####
##### get_nx, get_ny, get_nz #####

```fortran
integer function get_nx(self, grid)
   implicit none
   class(cm1_dataset) :: self
   character          :: grid
```
These get dimensions of the specified grid `grid`.

##### get_x, get_y, get_z #####

```fortran
integer function get_x(self, grid, x)
   implicit none
   class(cm1_dataset) :: self
   character          :: grid
   integer            :: gridno
   real, dimension(:) :: x
```
These get the mesh of grid `grid` along the specified dimension

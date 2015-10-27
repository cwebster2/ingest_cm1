!==========================================================================
!
!   Ingest_CM1
!
!   Casey Webster, Dept of Meteorology, Penn State
!   <ccw5079 at psu dot edu>
!
!   DESCRIPTION:
!    A Fortran 2003 module to access model output by the CM1 cloud model
!    by George Bryan.
!
!    Formats supported :
!       * grads style flat files
!         * multi file per time step (MPI)
!         * single file per time step
!
!                        Requirement: T dimension values must conincide
!                                     with filenames.
!
!    Future formats:
!      * grads style flat files
!        * all timesteps in one file
!      * netcdf
!      * native cm1hdf5 cdir output
!
!    TODO: for HDF support:
!     derived types (if no qvpert, use qv-q00 from basestate)
!
!   Copyright (C) 2014 Casey Webster - All Rights Reserved
!   You may use, distribute and modiy this code under the
!   terms of the BSD (3-clause) license.
!
!   You should have received a copy of the BSD (3-clause)
!   license with this file.  If not, please visit :
!   http://opensource.org/licenses/BSD-3-Clause or 
!   https://github.com/cwebster2/ingest_cm1/blob/master/LICENSE
!
!==========================================================================!

module ingest_cm1

   ! This module defines the base API for ingest_cm1
   use ingest_cm1_base

   ! These modules implement the backend access for output types
   use ingest_cm1_grads
   use ingest_cm1_grads_mpi
   use ingest_cm1_grads_single
#ifdef HDF5
   use ingest_cm1_hdf5
#endif

   implicit none

   ! Object to abstract access to multiple datasets (same type, different grids)
   ! CM1 datasets are output on 4 grids (Arakawa C grid) that vary for scalar, u, v, and w points
   ! This type will allow access to all from one interface
   type cm1_dataset
      private
      class(cm1_base), allocatable, dimension(:) :: cm1
      character, allocatable, dimension(:)       :: grids
      integer                       :: ngrids, dsettype
      logical :: isopen = .false.

      contains
         private
         procedure, pass(self)         :: get_gridno
         procedure, pass(self)         :: check_valid

         procedure, public, pass(self) :: read_3d
         procedure, public, pass(self) :: read_3d_slice
         procedure, public, pass(self) :: read_2d

         procedure, public, pass(self) :: open_dataset
         procedure, public, pass(self) :: close_dataset

         procedure, public, pass(self) :: get_nx
         procedure, public, pass(self) :: get_ny
         procedure, public, pass(self) :: get_nz
         procedure, public, pass(self) :: get_nt
         procedure, public, pass(self) :: get_x
         procedure, public, pass(self) :: get_y
         procedure, public, pass(self) :: get_z
         procedure, public, pass(self) :: get_t
         procedure, public, pass(self) :: get_var_exists

   end type cm1_dataset

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Opens a CM1 output dataset for reading.
  !! @param dsetpath The path to the CM1 dataset
  !! @param dsetbasename The basename of the CM1 dataset (e.g. cm1out)
  !! @param dsettype The dataset type
  !! @see enum
  !! @param grids An array of dataset grids to load (e.g. ['u', 'v', 'w', 's'])
  !! @param nodex Number of MPI nodes in the x direction for the dataset
  !! @param nodey Number of MPI nodes in the y direction for the dataset
  !! @return 1 on sucess, 0 on failure to open the dataset
  
   integer function open_dataset(self, dsetpath, dsetbasename, dsettype, grids, nodex, nodey)
      class(cm1_dataset) :: self
      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in) :: dsetbasename
      integer, intent(in)          :: dsettype
      character, dimension(:)      :: grids
      integer, optional :: nodex, nodey

      integer :: status, i

      if (self%isopen) then
         call cm1log(LOG_ERROR,'open_dataset','Error, dataset already open')
         open_dataset = 0
         return
      end if

      self%dsettype = dsettype
      self%ngrids = size(grids)
      self%grids = grids
      self%isopen = .false.
      open_dataset = 0

      select case (self%dsettype)
         case (GRADS)
            allocate(cm1_grads :: self%cm1(self%ngrids))
         case (GRADSMPI)
            allocate(cm1_grads_mpi :: self%cm1(self%ngrids))
         case (GRADSSINGLE)
            allocate(cm1_grads_single :: self%cm1(self%ngrids))
         case (HDF)
#ifdef HDF5
            allocate(cm1_hdf5 :: self%cm1(self%ngrids))
#else
            call cm1log(LOG_ERROR, 'open_dataset',"Ingest_CM1 not built with HDF support")
#endif
         case default
            call cm1log(LOG_ERROR, 'open_dataset',"Unrecognized dataset type")
            open_dataset = 0
            return
      end select

      do i=1,self%ngrids
         call cm1log(LOG_INFO, 'open_dataset',"Opening grid "//grids(i))
         select case (self%dsettype)
#ifdef HDF5
            ! HDF needs some special handling to intialize the library
            case (HDF)
               call initiate_hdf()
               status = self%cm1(i)%open_cm1(dsetpath, dsetbasename, dsettype, grids(i), nodex, nodey, .false.)
#endif
            case default
               status = self%cm1(i)%open_cm1(dsetpath, dsetbasename, dsettype, grids(i), nodex, nodey)
         end select
         if (status.ne.1) then
            call cm1log(LOG_ERROR,'open_dataset',"Error opening dataet")
            stop
         end if
      end do

      self%isopen = .true.
      open_dataset = 1

   end function open_dataset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_dataset(self)
      class(cm1_dataset) :: self
      integer :: i, status

      if (.not. self%isopen) then
         call cm1log(LOG_ERROR,'close_dataset','No dataset open')
         close_dataset = 0
         return
      end if

      do i=1,self%ngrids
         status = self%cm1(i)%close_cm1()
      end do

#ifdef HDF5
      if (self%dsettype .eq. HDF) then
         call deinitiate_hdf()
      endif
#endif

      deallocate(self%grids)
      deallocate(self%cm1)

      self%isopen = .false.
      close_dataset = 1

   end function close_dataset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_gridno(self, grid)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer :: i

      get_gridno = -1

      do i=1,self%ngrids
         if (grid .eq. self%grids(i)) then
            get_gridno = i
            return
         end if
      end do
   end function get_gridno

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function check_valid(self, funcname, grid)
      implicit none
      class(cm1_dataset) :: self
      character(len=*)   :: funcname
      character          :: grid

      check_valid = -1

      if (.not. self%isopen) then
         call cm1log(LOG_ERROR, trim(funcname), 'No dataset open')
         return
      end if

      check_valid = self%get_gridno(grid)

      if (check_valid .eq. -1) then
         call cm1log(LOG_ERROR, trim(funcname), 'Invalid grid specified')
         return
      end if

   end function check_valid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read_3d(self, time, grid, varname, Field3D)
      implicit none
      class(cm1_dataset) :: self
      integer            :: time, gridno
      character          :: grid
      character(len=*)   :: varname
      real, dimension(:,:,:) :: Field3D

      read_3d = 0

      gridno = self%check_valid('read_3d', grid)
      if (gridno .eq. -1) return

      read_3d = self%cm1(gridno)%read3D(varname, time, Field3D)
      
   end function read_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read_3d_slice(self, time, grid, varname, Field3D, ib, ie, jb, je, kb, ke)
      implicit none
      class(cm1_dataset) :: self
      integer            :: time, gridno, ib, ie, jb, je, kb, ke
      character          :: grid
      character(len=*)   :: varname
      real, dimension(:,:,:) :: Field3D

      read_3d_slice = 0

      gridno = self%check_valid('read_3d', grid)
      if (gridno .eq. -1) return

      read_3d_slice = self%cm1(gridno)%read3DSlice(varname, time, Field3D, ib, ie, jb, je, kb, ke)
      
   end function read_3d_slice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read_2d(self, time, grid, varname, Field2D)
      implicit none
      class(cm1_dataset) :: self
      integer            :: time, gridno
      character          :: grid
      character(len=*)   :: varname
      real, dimension(:,:) :: Field2D

      read_2d = 0

      gridno = self%check_valid('read_2d', grid)
      if (gridno .eq. -1) return

      read_2d = self%cm1(gridno)%read2D(varname, time, Field2D)

   end function read_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nx(self, grid)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno

      get_nx = 0

      gridno = self%check_valid('get_nx', grid)
      if (gridno .eq. -1) return

      get_nx = self%cm1(gridno)%get_nx()
      
   end function get_nx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_ny(self, grid)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno

      get_ny = 0

      gridno = self%check_valid('get_ny', grid)
      if (gridno .eq. -1) return

      get_ny = self%cm1(gridno)%get_ny()
      
   end function get_ny

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nz(self, grid)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno

      get_nz = 0

      gridno = self%check_valid('get_nz', grid)
      if (gridno .eq. -1) return

      get_nz = self%cm1(gridno)%get_nz()
      
   end function get_nz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nt(self, grid)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno

      get_nt = 0

      gridno = self%check_valid('get_nz', grid)
      if (gridno .eq. -1) return

      get_nt = self%cm1(gridno)%get_nt()
      
   end function get_nt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function get_x(self, grid, cm1err) result(x)
      implicit none
      class(cm1_dataset)       :: self
      character                :: grid
      integer                  :: gridno
      logical, optional        :: cm1err
      real, allocatable, dimension(:) :: x

      if (present(cm1err)) cm1err = .true.
      
      gridno = self%check_valid('get_x', grid)
      if (gridno .eq. -1) return
      
      x = self%cm1(gridno)%get_x(cm1err)
      if (present(cm1err)) cm1err = .false.
   end function get_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function get_y(self, grid, cm1err) result(y)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno
      logical, optional  :: cm1err
      real, allocatable, dimension(:) :: y

      if (present(cm1err)) cm1err = .true.
      
      gridno = self%check_valid('get_y', grid)
      if (gridno .eq. -1) return
      
      y = self%cm1(gridno)%get_y(cm1err)
      if (present(cm1err)) cm1err = .false.
   end function get_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function get_z(self, grid, cm1err) result(z)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno
      logical, optional  :: cm1err
      real, allocatable, dimension(:) :: z

      if (present(cm1err)) cm1err = .true.
      
      gridno = self%check_valid('get_z', grid)
      if (gridno .eq. -1) return
      
      z = self%cm1(gridno)%get_z(cm1err)
      if (present(cm1err)) cm1err = .false.
   end function get_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function get_t(self, grid, cm1err) result(t)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno
      logical, optional  :: cm1err
      real, allocatable, dimension(:) :: t

      if (present(cm1err)) cm1err = .true.
      
      gridno = self%check_valid('get_z', grid)
      if (gridno .eq. -1) return
      
      t= self%cm1(gridno)%get_t(cm1err)
      if (present(cm1err)) cm1err = .false.
   end function get_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function get_var_exists(self, grid, varname)
      implicit none
      class(cm1_dataset) :: self
      character          :: grid
      integer            :: gridno, varid
      character(len=*)   :: varname

      get_var_exists = .false.

      gridno = self%check_valid('get_var_exists', grid)
      if (gridno .eq. -1) return

      varid = self%cm1(gridno)%get_var_byname(varname)
      if (varid .eq. 0) return

      get_var_exists = .true.

   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine cm1log(loglevel, funcname, message)
      use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
      integer, intent(in) :: loglevel
      character(len=*), intent(in) :: funcname, message

      1001 format('[ingestcm1:',A,']: Error: ',A)
      1002 format('[ingestcm1:',A,']: Warning: ',A)
      1004 format('[ingestcm1:',A,']: ',A)

      select case (loglevel)

         case (LOG_ERROR)
            write(error_unit,1001) funcname, message

         case (LOG_WARN)
            write(error_unit,1002) funcname, message

         case (LOG_INFO)
            write(error_unit,1004) funcname, message

         case (LOG_MSG)
            write(output_unit,1004) funcname, message

      end select

   end subroutine

end module ingest_cm1

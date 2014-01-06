!==========================================================================!
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
!==========================================================================!

module ingest_cm1

   ! This module defines the base API for ingest_cm1
   use ingest_cm1_base

   ! These modules implement the backend access for output types
   use ingest_cm1_grads
   use ingest_cm1_grads_mpi
   use ingest_cm1_hdf5

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
      logical, allocatable, dimension(:) :: ismult

      contains
         private
         procedure, public, pass(self) :: open_dataset
         procedure, public, pass(self) :: close_dataset

   end type cm1_dataset

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_dataset(self, dsetpath, dsetbasename, dsettype, grids, nodex, nodey, hdfmetadatatime)
      class(cm1_dataset) :: self
      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in) :: dsetbasename
      integer, intent(in)          :: dsettype
      character, dimension(:)      :: grids
      integer, optional :: nodex, nodey, hdfmetadatatime

      integer :: status, i

      if (self%isopen) then
         print *,'Error, dataset already open'
         open_dataset = 0
         return
      end if

      self%dsettype = dsettype
      self%ngrids = size(grids)
      self%grids = grids
      allocate(self%ismult(self%ngrids))
      self%ismult = .false.
      self%isopen = .false.
      open_dataset = 0

      select case (self%dsettype)
         case (GRADS)
            allocate(cm1_grads :: self%cm1(self%ngrids))
         case (GRADSMPI)
            allocate(cm1_grads_mpi :: self%cm1(self%ngrids))
         case (HDF)
            print *,"HDFALLOC", self%ngrids
            allocate(cm1_hdf5 :: self%cm1(self%ngrids))
         case default
            print *,"Unrecognized dataset type"
            open_dataset = 0
            return
      end select

      do i=1,self%ngrids
         print *,"Opening grid "//grids(i), i
         status = self%cm1(i)%open_cm1(dsetpath, dsetbasename, dsettype, grids(i), nodex, nodey, hdfmetadatatime, .true.)
         if (status.ne.1) then
            print *,"Error opening dataet"
            !stop
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
         print *,'Error, no dataset open'
         close_dataset = 0
         return
      end if

      do i=1,self%ngrids
         status = self%cm1(i)%close_cm1()
      end do

      deallocate(self%grids)
      deallocate(self%cm1)

      self%isopen = .false.
      close_dataset = 1

   end function close_dataset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ingest_cm1

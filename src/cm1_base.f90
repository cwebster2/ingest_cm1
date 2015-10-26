!==========================================================================!
!
!   Ingest_CM1_base
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

module ingest_cm1_base

   implicit none

   ! module variables
   private
   public :: variable, cm1_base, min_loglevel

! Parameters to specify dataset type
   integer, public, parameter :: GRADSSINGLE = 2
   integer, public, parameter :: GRADS = 3
   integer, public, parameter :: GRADSMPI = 4
   integer, public, parameter :: HDF = 5

   type variable
      character(len=32) :: varname
      integer           :: levs
      character(len=64) :: units
      integer           :: dims
   end type variable


   ! Abstract polymorphic type to establish interface for all dataset types
   type, abstract :: cm1_base
      !private
         character(len=128) :: path
         character(len=64)  :: basename
         integer            :: dtype
         character          :: grid
         integer :: nx, ny, nz, nt, nv, dt, time_open
         real, dimension(:), allocatable :: x, y, z
         integer, dimension(:), allocatable :: times
         type (variable), dimension(:), allocatable  :: vars
         logical :: is_dataset_open = .false., is_time_open = .false.
      contains
         ! Private procedures
         procedure, pass(self) :: check_dataset_open
         procedure, pass(self) :: check_time_open
         procedure, pass(self) :: cm1log
         procedure(open_dataset_time_i), deferred, public ,pass(self) :: open_dataset_time
         procedure(close_dataset_time_i), deferred, public ,pass(self) :: close_dataset_time
         procedure(read3D_at_time_i), deferred, public ,pass(self) :: read3D_at_time
         procedure(read2D_at_time_i), deferred, public ,pass(self) :: read2D_at_time
         
         ! Public getter functions for generic dataset attributes
         procedure, public ,pass(self) :: get_nx
         procedure, public ,pass(self) :: get_ny
         procedure, public ,pass(self) :: get_nz
         procedure, public ,pass(self) :: get_nt
         procedure, public ,pass(self) :: get_nv
         procedure, public ,pass(self) :: get_x
         procedure, public ,pass(self) :: get_y
         procedure, public ,pass(self) :: get_z
         procedure, public ,pass(self) :: get_t
         procedure, public ,pass(self) :: get_var_byname
         procedure, public ,pass(self) :: read3D
         procedure, public ,pass(self) :: read3DSlice
         procedure, public ,pass(self) :: read2D
         procedure(open_cm1_i), deferred, public ,pass(self) :: open_cm1
         procedure(close_cm1_i),deferred, public ,pass(self) :: close_cm1
         
   end type cm1_base

   abstract interface
      integer function open_cm1_i(self, dsetpath, dsetbasename, dsettype, grid, nodex, nodey, hdfmanage)
         import cm1_base
         class(cm1_base) :: self
         character(len=*), intent(in) :: dsetpath
         character(len=*), intent(in) :: dsetbasename
         integer, intent(in)          :: dsettype
         character, optional :: grid
         integer, optional :: nodex, nodey
         logical, optional :: hdfmanage
      end function open_cm1_i

      integer function close_cm1_i(self)
         import cm1_base
         class(cm1_base) :: self
      end function close_cm1_i

      integer function open_dataset_time_i(self,time)
         import cm1_base
         class(cm1_base)          :: self
         integer, intent(in) :: time
      end function open_dataset_time_i

      integer function close_dataset_time_i(self)
         import cm1_base
         class(cm1_base) :: self
      end function close_dataset_time_i

      integer function read2D_at_time_i(self, varname, Field2D)
         import cm1_base
         class(cm1_base) :: self
         character(len=*), intent(in) :: varname
         real, dimension(self%nx,self%ny) :: Field2D
      end function read2D_at_time_i

      integer function read3D_at_time_i(self, varname, Field3D)
         import cm1_base
         class(cm1_base) :: self
         character(len=*), intent(in) :: varname
         real, dimension(self%nx,self%ny,self%nz) :: Field3D
      end function read3D_at_time_i
   end interface

   !for logging
   integer, public, parameter :: LOG_ERROR = 1005
   integer, public, parameter :: LOG_WARN  = 1004
   integer, public, parameter :: LOG_INFO  = 1003
   integer, public, parameter :: LOG_MSG   = 1002
   integer, public, parameter :: LOG_DEBUG = 1001

   integer :: min_loglevel = LOG_INFO

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_var_byname(self,varname)
      implicit none
      class(cm1_base)                   :: self
      character(len=*), intent(in) :: varname
      integer :: i

      get_var_byname = 0
      if (.not. self%check_dataset_open('get_var_byname')) return

      do i = 1, self%nv
         if (varname.eq.self%vars(i)%varname) then
            get_var_byname = i
            return
         endif
      end do

      call self%cm1log(LOG_WARN, 'get_var_byname', 'Requested variable does not exist: '//varname)
      
   end function get_var_byname

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2D(self, varname, time, Field2D)
      implicit none
      class(cm1_base), intent(in)   :: self
      character(len=*), intent(in) :: varname
      integer, intent(in) :: time
      real, dimension(self%nx,self%ny) :: Field2D
      integer :: status

      if (.not. self%check_dataset_open('read2D')) then
         read2D = 0
         return
      end if

      status = self%open_dataset_time(time)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read2D', 'Data open failure, aborting: ')
         read2D = 0
         return
      endif

      status = self%read2D_at_time(varname, Field2D)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read2D', 'Data read failure, aborting: ')
         read2D = 0
         status = self%close_dataset_time()
         return
      endif

      status = self%close_dataset_time()
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read2D', 'Data close failure, aborting: ')
         read2D = 0
         return
      endif

      read2D = 1

   end function read2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3D(self, varname, time, Field3D)
      implicit none
      class(cm1_base), intent(in)   :: self
      character(len=*), intent(in) :: varname
      integer, intent(in) :: time
      real, dimension(self%nx,self%ny,self%nz) :: Field3D
      integer :: status

      if (.not. self%check_dataset_open('read3D')) then
         read3D = 0
         return
      end if

      status = self%open_dataset_time(time)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3D', 'Data open failure, aborting: ')
         read3D = 0
         return
      endif

      status = self%read3D_at_time(varname, Field3D)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3D', 'Data read failure, aborting: ')
         read3D = 0
         status = self%close_dataset_time()
         return
      endif

      status = self%close_dataset_time()
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3D', 'Data close failure, aborting: ')
         read3D = 0
         return
      endif

      read3D = 1
   end function read3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DSlice(self, varname, time, Field3D, ib, ie, jb, je, kb, ke)
      implicit none
      class(cm1_base), intent(in)   :: self
      character(len=*), intent(in) :: varname
      integer, intent(in) :: time, ib, ie, jb, je, kb, ke
      real, dimension(ib:ie, jb:je, kb:ke) :: Field3D
      integer :: status
      real, dimension(self%nx,self%ny,self%nz) :: Field3DFull
      
      if (.not. self%check_dataset_open('read3DSlice')) then
         read3DSlice = 0
         return
      end if

      if ((ib < 0) .or. (jb < 0) .or. (kb < 0) .or.  &
          (ie > self%nx) .or. (je > self%ny) .or. (ke > self%nz)) then
         call self%cm1log(LOG_ERROR, 'read3DSice', 'Slice out of bounds, aborting')
         read3DSlice = 0
         return
      endif
 
      status = self%read3D(varname, time, Field3DFull)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3DSice', 'Data read failure, aborting.')
         read3DSlice = 0
         return
      endif

      Field3D = Field3DFull(ib:ie, jb:je, kb:ke)

      read3DSlice = 1
   end function read3DSlice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nx(self)
      implicit none
      class(cm1_base), intent(in) :: self

      get_nx = 0
      if (.not. self%check_dataset_open('get_nx')) return

      get_nx = self%nx
   end function get_nx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_ny(self)
      implicit none
      class(cm1_base), intent(in) :: self

      get_ny = 0
      if (.not. self%check_dataset_open('get_ny')) return

      get_ny = self%ny
   end function get_ny

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nz(self)
      implicit none
      class(cm1_base), intent(in) :: self

      get_nz = 0
      if (.not. self%check_dataset_open('get_nz')) return

      get_nz = self%nz
   end function get_nz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nt(self)
      implicit none
      class(cm1_base), intent(in) :: self

      get_nt = 0
      if (.not. self%check_dataset_open('get_nt')) return

      get_nt = self%nt
   end function get_nt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_nv(self)
      implicit none
      class(cm1_base), intent(in) :: self

      get_nv = 0
      if (.not. self%check_dataset_open('get_nv')) return

      get_nv = self%nv
   end function get_nv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_dt(self)
      implicit none
      class(cm1_base), intent(in) :: self

      get_dt = 0
      if (.not. self%check_dataset_open('get_dt')) return

      get_dt = self%dt
   end function get_dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_x(self, cm1x)
      implicit none
      class(cm1_base), intent(in) :: self
      real, dimension(self%nx) :: cm1x

      get_x = 0
      if (.not. self%check_dataset_open('get_x')) return

      cm1x(:) = self%x(:)
      get_x = 1
   end function get_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_y(self, cm1y)
      implicit none
      class(cm1_base), intent(in) :: self
      real, dimension(self%ny) :: cm1y

      get_y = 0
      if (.not. self%check_dataset_open('get_y')) return

      cm1y(:) = self%y(:)
      get_y = 1
   end function get_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_z(self, cm1z)
      implicit none
      class(cm1_base), intent(in) :: self
      real, dimension(self%nz) :: cm1z

      get_z = 0
      if (.not. self%check_dataset_open('get_z')) return

      cm1z(:) = self%z(:)
      get_z = 1
   end function get_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_t(self, cm1t)
      implicit none
      class(cm1_base), intent(in) :: self
      real, dimension(self%nt) :: cm1t

      get_t = 0 
      if (.not. self%check_dataset_open('get_x')) return

      cm1t(:) = self%times(:)
      get_t = 1
   end function get_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function get_time_open(self)
     implicit none
     class(cm1_base), intent(in) :: self

     get_time_open = -1
     if (.not. self%check_dataset_open('get_time_open')) return
     if (.not. self%check_time_open('get_time_open')) return

     get_time_open = self%time_open
   end function get_time_open
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function check_dataset_open(self, funcname)
      use, intrinsic :: iso_fortran_env, only: error_unit
      class(cm1_base), intent(in) :: self
      character(len=*), intent(in) :: funcname

      ! Is the dataset open (has the ctl file been scanned)
      if (.not. self%is_dataset_open) then
         call self%cm1log(LOG_ERROR, funcname, 'No dataset open, aborting')
         check_dataset_open = .false.
         return
      end if
      check_dataset_open = .true.

   end function check_dataset_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function check_time_open(self, funcname)
      use, intrinsic :: iso_fortran_env, only: error_unit
      class(cm1_base),intent(in) :: self
      character(len=*), intent(in) :: funcname

      ! Is the dataset open for reading (is the dat file open)
      if (.not. self%is_time_open) then
         call cm1log(self, LOG_ERROR, funcname, 'No multiread started, aborting')
         check_time_open = .false.
         return
      end if
      check_time_open = .true.

   end function check_time_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine cm1log(self, loglevel, funcname, message)
     use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
     implicit none
     class(cm1_base), intent(in) :: self
     integer, intent(in) :: loglevel
     character(len=*), intent(in) :: funcname, message
     
     character(len=6) :: proto
     
     character(len=*),parameter :: fmt_error = "('[ingestcm1',3(':',A),']: Error: ',A)"
     character(len=*),parameter :: fmt_warning = "('[ingestcm1',3(':',A),']: Warning: ',A)"
     character(len=*),parameter :: fmt_message = "('[ingestcm1',3(':',A),']: ',A)"
     
     if (loglevel .lt. min_loglevel) return
     
     select case (self%dtype)
        
        case (GRADS)
           proto = 'GRADS'
        
        case (GRADSMPI)
           proto = 'GRMPI'
        
        case (HDF)
           proto = "HDF5"
        
        case (GRADSSINGLE)
           proto = "GRADS"
        
        case default
           proto = 'X'
        
     end select
     
     select case (loglevel)

        case (LOG_ERROR)
           write(error_unit,fmt_error) trim(proto), self%grid, funcname, message
        
        case (LOG_WARN)
           write(error_unit,fmt_warning) trim(proto), self%grid, funcname, message
        
        case (LOG_INFO)
           write(output_unit,fmt_message) trim(proto), self%grid, funcname, message
        
        case (LOG_MSG)
           write(output_unit,fmt_message) trim(proto), self%grid, funcname, message
        
        case (LOG_DEBUG)
           write(output_unit,fmt_message) trim(proto), self%grid, funcname, message
        
     end select
     
   end subroutine cm1log
   
 end module ingest_cm1_base

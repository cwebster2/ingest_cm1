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
         integer :: nx, ny, nz, nt, nv, dt
         real, dimension(:), allocatable :: x, y, z
         integer, dimension(:), allocatable :: times
         type (variable), dimension(:), allocatable  :: vars
         logical :: isopen = .false., ismult = .false.
      contains
         ! Private procedures
         procedure, pass(self) :: check_open
         procedure, pass(self) :: check_mult
         procedure, pass(self) :: cm1log

         ! Public getter functions for generic dataset attributes
         procedure, public ,pass(self) :: cm1_nx
         procedure, public ,pass(self) :: cm1_ny
         procedure, public ,pass(self) :: cm1_nz
         procedure, public ,pass(self) :: cm1_nt
         procedure, public ,pass(self) :: cm1_nv
         procedure, public ,pass(self) :: cm1_x
         procedure, public ,pass(self) :: cm1_y
         procedure, public ,pass(self) :: cm1_z
         procedure, public ,pass(self) :: cm1_t
         procedure, public ,pass(self) :: getVarByName
         procedure, public ,pass(self) :: read3D
         procedure, public ,pass(self) :: read2D

         ! deferred procedures
         procedure(open_cm1_i), deferred, public ,pass(self) :: open_cm1
         procedure(close_cm1_i),deferred, public ,pass(self) :: close_cm1
         procedure(readMultStart_i), deferred, public ,pass(self) :: readMultStart
         procedure(readMultStop_i), deferred, public ,pass(self) :: readMultStop
         procedure(read3DMult_i), deferred, public ,pass(self) :: read3DMult
         procedure(read2DMult_i), deferred, public ,pass(self) :: read2DMult
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

      integer function readMultStart_i(self,time)
         import cm1_base
         class(cm1_base)          :: self
         integer, intent(in) :: time
      end function readMultStart_i

      integer function readMultStop_i(self)
         import cm1_base
         class(cm1_base) :: self
      end function readMultStop_i

      integer function read2DMult_i(self, varname, Field2D)
         import cm1_base
         class(cm1_base) :: self
         character(len=*), intent(in) :: varname
         real, dimension(self%nx,self%ny) :: Field2D
      end function read2DMult_i

      integer function read3DMult_i(self, varname, Field3D)
         import cm1_base
         class(cm1_base) :: self
         character(len=*), intent(in) :: varname
         real, dimension(self%nx,self%ny,self%nz) :: Field3D
      end function read3DMult_i
   end interface

   !for logging
   integer, public, parameter :: LOG_ERROR = 1001
   integer, public, parameter :: LOG_WARN  = 1002
   integer, public, parameter :: LOG_INFO  = 1003
   integer, public, parameter :: LOG_MSG   = 1004

   integer :: min_loglevel = LOG_ERROR

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function getVarByName(self,varname)
      implicit none
      class(cm1_base)                   :: self
      character(len=*), intent(in) :: varname
      integer :: i

      if (.not. self%check_open('getVarByName')) then
         getVarByName = 0
         return
      end if

      do i = 1, self%nv
         if (varname.eq.self%vars(i)%varname) then
            getVarByName = i
            return
         endif
      end do

      getVarByName = 0

   end function getVarByName

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2D(self, varname, time, Field2D)
      implicit none
      class(cm1_base), intent(in)   :: self
      character(len=*), intent(in) :: varname
      integer, intent(in) :: time
      real, dimension(self%nx,self%ny) :: Field2D
      integer :: status

      if (.not. self%check_open('read2D')) then
         read2D = 0
         return
      end if

      status = self%readMultStart(time)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read2D', 'Data open failure, aborting: ')
         read2D = 0
         return
      endif

      status = self%read2DMult(varname, Field2D)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read2D', 'Data read failure, aborting: ')
         read2D = 0
         status = self%readMultStop()
         return
      endif

      status = self%readMultStop()
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

      if (.not. self%check_open('read3D')) then
         read3D = 0
         return
      end if

      status = self%readMultStart(time)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3D', 'Data open failure, aborting: ')
         read3D = 0
         return
      endif

      status = self%read3DMult(varname, Field3D)
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3D', 'Data read failure, aborting: ')
         read3D = 0
         status = self%readMultStop()
         return
      endif

      status = self%readMultStop()
      if (status.eq.0) then
         call self%cm1log(LOG_ERROR, 'read3D', 'Data close failure, aborting: ')
         read3D = 0
         return
      endif

      read3D = 1
   end function read3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nx(self)
      implicit none
      class(cm1_base),intent(in) :: self

      if (.not. self%check_open('cm1_nx')) then
         cm1_nx = 0
         return
      end if

      cm1_nx = self%nx
   end function cm1_nx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_ny(self)
      implicit none
      class(cm1_base),intent(in) :: self

      if (.not. self%check_open('cm1_ny')) then
         cm1_ny = 0
         return
      end if

      cm1_ny = self%ny
   end function cm1_ny

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nz(self)
      implicit none
      class(cm1_base),intent(in) :: self

      if (.not. self%check_open('cm1_nz')) then
         cm1_nz = 0
         return
      end if

      cm1_nz = self%nz
   end function cm1_nz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nt(self)
      implicit none
      class(cm1_base),intent(in) :: self

      if (.not. self%check_open('cm1_nt')) then
         cm1_nt = 0
         return
      end if

      cm1_nt = self%nt
   end function cm1_nt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nv(self)
      implicit none
      class(cm1_base),intent(in) :: self

      if (.not. self%check_open('cm1_nv')) then
         cm1_nv = 0
         return
      end if

      cm1_nv = self%nv
   end function cm1_nv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_dt(self)
      implicit none
      class(cm1_base),intent(in) :: self

      if (.not. self%check_open('cm1_dt')) then
         cm1_dt = 0
         return
      end if

      cm1_dt = self%dt
   end function cm1_dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_x(self, cm1x)
      implicit none
      class(cm1_base),intent(in) :: self
      real, dimension(self%nx) :: cm1x

      if (.not. self%check_open('cm1_x')) then
         cm1_x = 0
         return
      end if

      cm1x(:) = self%x(:)
      cm1_x = 1
   end function cm1_x

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_y(self, cm1y)
      implicit none
      class(cm1_base),intent(in) :: self
      real, dimension(self%ny) :: cm1y

      if (.not. self%check_open('cm1_x')) then
         cm1_y = 0
         return
      end if

      cm1y(:) = self%y(:)
      cm1_y = 1
   end function cm1_y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_z(self, cm1z)
      implicit none
      class(cm1_base),intent(in) :: self
      real, dimension(self%nz) :: cm1z

      if (.not. self%check_open('cm1_x')) then
         cm1_z = 0
         return
      end if

      cm1z(:) = self%z(:)
      cm1_z = 1
   end function cm1_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_t(self, cm1t)
      implicit none
      class(cm1_base),intent(in) :: self
      real, dimension(self%nt) :: cm1t

      if (.not. self%check_open('cm1_x')) then
         cm1_t = 0
         return
      end if

      cm1t(:) = self%times(:)
      cm1_t = 1
   end function cm1_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function check_open(self, funcname)
      use, intrinsic :: iso_fortran_env, only: error_unit
      class(cm1_base),intent(in) :: self
      character(len=*), intent(in) :: funcname

      ! Is the dataset open (has the ctl file been scanned)
      if (.not. self%isopen) then
         call self%cm1log(LOG_ERROR, funcname, 'No dataset open, aborting')
         check_open = .false.
         return
      end if
      check_open = .true.

   end function check_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function check_mult(self, funcname)
      use, intrinsic :: iso_fortran_env, only: error_unit
      class(cm1_base),intent(in) :: self
      character(len=*), intent(in) :: funcname

      ! Is the dataset open for reading (is the dat file open)
      if (.not. self%ismult) then
         call cm1log(self, LOG_ERROR, funcname, 'No multiread started, aborting')
         check_mult = .false.
         return
      end if
      check_mult = .true.

   end function check_mult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine cm1log(self, loglevel, funcname, message)
      use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
      class(cm1_base), intent(in) :: self
      integer, intent(in) :: loglevel
      character(len=*), intent(in) :: funcname, message

      character(len=6) :: proto

      1001 format('[ingestcm1',3(':',A),']: Error: ',A)
      1002 format('[ingestcm1',3(':',A),']: Warning: ',A)
      1004 format('[ingestcm1',3(':',A),']: ',A)

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
            write(error_unit,1001) trim(proto), self%grid, funcname, message

         case (LOG_WARN)
            write(error_unit,1002) trim(proto), self%grid, funcname, message

         case (LOG_INFO)
            write(error_unit,1004) trim(proto), self%grid, funcname, message

         case (LOG_MSG)
            write(output_unit,1004) trim(proto), self%grid, funcname, message

      end select

   end subroutine

end module ingest_cm1_base

!==========================================================================!
!
!   Ingest_CM1_grads_single
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

module ingest_cm1_grads_single

   use ingest_cm1_grads

   implicit none

   type, extends(cm1_grads) :: cm1_grads_single
      integer t ! keep track of what time level is "selected"
    contains
      !private

      procedure, pass(self) :: read3DXYSlice => read3DXYSlice_grads_single
      procedure, public ,pass(self) :: open_cm1 => open_cm1_grads_single
      procedure, public ,pass(self) :: close_cm1 => close_cm1_grads_single
      procedure, public ,pass(self) :: open_dataset_time => open_dataset_time_grads_single
      procedure, public ,pass(self) :: close_dataset_time => close_dataset_time_grads_single
      
      ! This ensures that the filehandles and arrays are closed if the variable goes out of scope
      !TODO: finalization needs gfortran 4.9 or ifort
      !!!final :: final_cm1

   end type cm1_grads_single

contains

   integer function open_cm1_grads_single(self, dsetpath, dsetbasename, dsettype, &
                                   grid, nodex, nodey, hdfmanage) result(open_cm1)

      implicit none
      class(cm1_grads_single) :: self

      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in)  :: dsetbasename
      integer, intent(in)            :: dsettype
      character, optional :: grid
      integer, optional :: nodex, nodey
      logical, optional :: hdfmanage
      character(len=256) :: datfile


      if (self%is_dataset_open) then
         call self%cm1log(LOG_WARN, 'open_cm1', 'Already open, aborting')
         open_cm1 = 0
         return
      end if

      self%path = dsetpath
      self%basename = dsetbasename
      self%dtype = dsettype

      if (present(grid)) then
        self%grid = grid
      else
        self%grid = 's'
      endif
      call self%cm1log(LOG_MSG, 'open_cm1', 'Grid ('//self%grid//') selected.')

      call self%cm1log(LOG_MSG, 'open_cm1', 'Reading GRADS control file.')
      open_cm1 = self%read_ctl()
      self%nunits = 1
      allocate(self%dat_units(self%nunits))
      
      datfile = trim(self%path)//'/'//trim(self%basename)//'_'//self%grid//'.dat'
      call self%cm1log(LOG_MSG, 'open_cm1', 'Opening: '//trim(datfile))

      ! TODO: check if time is in times for vailidy
      !       Make sure file successfully opens!

      ! open dat file
      open(newunit=self%dat_units(1),file=datfile,form='unformatted',access='direct',recl=self%reclen,status='old')
      self%is_dataset_open = .true.

   end function open_cm1_grads_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   subroutine final_cm1(self)
!      implicit none
!      type(cm1) :: self
!      integer :: stat
!
!      stat = close_cm1(self)
!
!   end subroutine final_cm1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_cm1_grads_single(self) result(close_cm1)
      implicit none
      class(cm1_grads_single) :: self

      close_cm1 = 0
      ! Check if the dataset is open
      if (.not. self%check_dataset_open('close_cm1')) return

      ! If data files are open, close them
      if (self%is_time_open) then
         close_cm1 = self%close_dataset_time()
      end if
      close(self%dat_units(1))
      
      deallocate(self%x)
      deallocate(self%y)
      deallocate(self%z)
      deallocate(self%times)
      deallocate(self%vars)
      deallocate(self%dat_units)
      self%nx = 0
      self%ny = 0
      self%nz = 0
      self%nt = 0
      self%nv = 0

      self%is_dataset_open = .false.
      close_cm1 = 1
      call self%cm1log(LOG_MSG, 'close_cm1', 'Dataset closed:')

   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DXYSlice_grads_single(self, varid, level, slice) result(read3DXYSlice)
      implicit none
      class(cm1_grads_single)             :: self
      integer, intent(in)    :: varid, level
      real, dimension(self%nx,self%ny) :: slice
      integer :: idx

      if (.not. self%check_dataset_open('read3DXYSlice')) then
         call self%cm1log(LOG_ERROR, 'read3DXYSlice', 'No datasef open, aborting')
         read3DXYSlice = 0
         return
      end if

      ! NOTE this offset is for one timestep per file!
      ! Calculate record index for 2D or 3D field
      if (self%vars(varid)%levs == 0) then ! 2D
          idx = ((self%t-1)*(self%n2d+(self%nv-self%n2d)*self%nz)) + varid-1
       else ! 3D
          idx = ((self%t-1)*(self%n2d+(self%nv-self%n2d)*self%nz)) + (self%n2d+(varid-1-self%n2d)*(self%nz)+level)
       endif

       ! Read Slice
       read(self%dat_units(1), rec=idx) slice
       read3DXYSlice = 1

   end function read3DXYSlice_grads_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_dataset_time_grads_single(self,time) result(open_dataset_time)
      implicit none
      class(cm1_grads_single) :: self
      integer, intent(in) :: time
      character(len=6) :: dtime

      open_dataset_time = 0
      if (.not. self%check_dataset_open('open_dataset_time')) return

      if (self%is_time_open) then
         call self%cm1log(LOG_WARN, 'open_dataset_time', 'Internal Error: timelevel inconsistency')
         stop
      end if

      self%t = time
      505 format(I6.6)
      ! filename?
      write(dtime,505) self%t
      call self%cm1log(LOG_MSG, 'read3DMultStart', 'Multiread started for time: '//trim(dtime))
      open_dataset_time = 1
      self%time_open = time
      self%is_time_open = .true.

    end function open_dataset_time_grads_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_dataset_time_grads_single(self) result(close_dataset_time)
      implicit none
      class(cm1_grads_single) :: self

      close_dataset_time = 0
      if ((.not. self%check_dataset_open('close_dataset_time')) .or. &
           (.not. self%check_time_open('close_dataset_time'))) return

      self%t = 0
      call self%cm1log(LOG_MSG, 'read3DMultStop', 'Multiread stopped.')
      close_dataset_time = 1
      self%time_open = -1
      self%is_time_open = .false.

    end function close_dataset_time_grads_single


end module ingest_cm1_grads_single

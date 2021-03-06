!==========================================================================!
!
!   Ingest_CM1_grads
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

module ingest_cm1_grads

   use ingest_cm1_base

   implicit none

   type, extends(cm1_base) :: cm1_grads
      !private
      integer            :: t0

      ! GRADS stuff
      integer            :: n2d, reclen
      integer            :: ctl_unit
      integer,dimension(:),allocatable :: dat_units
      integer            :: nunits

      contains
         !private

         procedure, pass(self) :: read_ctl
         procedure, pass(self) :: read3DXYSlice => read3DXYSlice_grads
         procedure, public ,pass(self) :: open_cm1 => open_cm1_grads
         procedure, public ,pass(self) :: close_cm1 => close_cm1_grads
         procedure, public ,pass(self) :: open_dataset_time => open_dataset_time_grads
         procedure, public ,pass(self) :: close_dataset_time => close_dataset_time_grads
         procedure, public ,pass(self) :: read3D_at_time => read3D_at_time_grads
         procedure, public ,pass(self) :: read2D_at_time => read2D_at_time_grads

         ! This ensures that the filehandles and arrays are closed if the variable goes out of scope
         !TODO: finalization needs gfortran 4.9 or ifort
         !!!final :: final_cm1

   end type cm1_grads

contains

   integer function open_cm1_grads(self, dsetpath, dsetbasename, dsettype, &
                                   grid, nodex, nodey, hdfmanage) result(open_cm1)

      implicit none
      class(cm1_grads) :: self

      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in)  :: dsetbasename
      integer, intent(in)            :: dsettype
      character, optional :: grid
      integer, optional :: nodex, nodey
      logical, optional :: hdfmanage

      if (self%is_dataset_open) then
         call self%cm1log(LOG_ERROR, 'open_cm1', 'Already open, aborting')
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
      self%is_dataset_open = .true.

   end function open_cm1_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read_ctl(self)
      implicit none
      class(cm1_grads) :: self
      character(len=256) :: dset
      character(len=128) :: tmp, output
      integer            :: i,j,k
      real, allocatable, dimension(:,:) :: recltest

      dset = trim(self%path)//'/'//trim(self%basename)//'_'//self%grid//'.ctl'
      call self%cm1log(LOG_MSG, 'read_ctl', 'Opening: '//trim(dset))
      open(newunit=self%ctl_unit, file=dset, status='old')
      if (self%dtype .ne. GRADSSINGLE) then
         read(self%ctl_unit,*) ! This line wont exist in single
      endif
      read(self%ctl_unit,*)
      read(self%ctl_unit,*)
      read(self%ctl_unit,*)

      500 format('... Found ',A1,' dimension with ',i5,' points')
      read(self%ctl_unit,*) tmp, self%nx  ! for stretched grids only.  detect and calc linear.
      allocate(self%x(self%nx))
      do i = 1, self%nx
         read(self%ctl_unit,*) self%x(i)
      end do
      write (output,500) 'X',self%nx
      call self%cm1log(LOG_INFO, 'read_ctl', trim(output))

      read(self%ctl_unit,*) tmp, self%ny  ! for stretched grids only.  detect and calc linear.
      allocate(self%y(self%ny))
      do j = 1, self%ny
         read(self%ctl_unit,*) self%y(j)
      end do
      write (output,500) 'Y',self%ny
      call self%cm1log(LOG_INFO, 'read_ctl', trim(output))

      read(self%ctl_unit,*) tmp, self%nz  ! for stretched grids only.  detect and calc linear.
      allocate(self%z(self%nz))
      do k = 1, self%nz
         read(self%ctl_unit,*) self%z(k)
      end do
      write (output,500) 'Z',self%nz
      call self%cm1log(LOG_INFO, 'read_ctl', trim(output))

      ! calculate timelevels.  note - supports a specific timeformat
      ! and only an integer dt (minimum timestep 1 s).
      ! tdef     NT linear hh:mmZddMMMYYYY SSSSSYR
      ! where NT is number of timelevels
      !       YYYY is the starting time in seconds
      !       SSSSS is the timestep in seconds
      !
      ! Modify cm1 to write your ctl file this way, or hand edit them
      read(self%ctl_unit,501) self%nt, self%t0, self%dt
      !!! format('tdef ',I10,' linear ',11x,I4,1x,I5,'YR')
      501 format( 5x    ,I10, 8x,       11x,I4,1x,I5, 2x)
      allocate(self%times(self%nt))
      do i = 1, self%nt
         self%times(i) = self%t0 + ((i-1)*self%dt)
      end do
      write (output,500) 'T',self%nt
      call self%cm1log(LOG_INFO, 'read_ctl', trim(output))

      ! variables
      read(self%ctl_unit,*) tmp, self%nv
      allocate(self%vars(self%nv))
      self%n2d = 0
      do i = 1,self%nv
         read(self%ctl_unit,502) self%vars(i)%varname, self%vars(i)%levs, self%vars(i)%units
         502 format(A10,1x,I5,6x,A)
         if (self%vars(i)%levs.eq.0) then
            self%n2d = self%n2d+1
         end if
      end do
      503 format('... Found ',i3,' variables')
      write (output,503) self%nv
      call self%cm1log(LOG_INFO, 'read_ctl', trim(output))


      ! This determines the recl to pass to open()
      allocate(recltest(self%nx,self%ny))
      inquire(iolength=self%reclen) recltest
      deallocate(recltest)
      504 format('... Using record length = ', I10)
      write (output,504) self%reclen
      call self%cm1log(LOG_DEBUG, 'read_ctl', trim(output))

      call self%cm1log(LOG_MSG, 'read_ctl', 'Closing Grads control file: '//trim(dset))
      close(self%ctl_unit)
      read_ctl = 1

   end function read_ctl

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

   integer function close_cm1_grads(self) result(close_cm1)
      implicit none
      class(cm1_grads) :: self

      ! Check if the dataset is open
      close_cm1 = 0
      if (.not. self%check_dataset_open('close_cm1')) return

      ! If data files are open, close them
      if (self%is_time_open) then
         close_cm1 = self%close_dataset_time()
      end if

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

   integer function read3DXYSlice_grads(self, varid, level, slice) result(read3DXYSlice)
      implicit none
      class(cm1_grads)             :: self
      integer, intent(in)    :: varid, level
      real, dimension(self%nx,self%ny) :: slice
      integer :: idx

      if (.not. self%check_dataset_open('read3DXYSlice')) then
         call self%cm1log(LOG_ERROR, 'read3DXYSlice', 'No dataset open, aborting')
         read3DXYSlice = 0
         return
      end if

      ! NOTE this offset is for one timestep per file!
      ! Calculate record index for 2D or 3D field
      if (self%vars(varid)%levs == 0) then ! 2D
          idx = varid-1
       else ! 3D
          idx = (self%n2d+(varid-1-self%n2d)*(self%nz)+level)
       endif

       ! Read Slice
       read(self%dat_units(1), rec=idx) slice
       read3DXYSlice = 1

   end function read3DXYSlice_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_dataset_time_grads(self,time) result(open_dataset_time)
     implicit none
     class(cm1_grads)          :: self
     integer, intent(in) :: time
     character(len=6) :: dtime
     character(len=256) :: datfile
     
     open_dataset_time = 0
     if (.not. self%check_dataset_open('open_dataset_time')) return
     
     if (self%is_time_open) then
        call self%cm1log(LOG_ERROR, 'open_dataset_time', 'Internal Error: different timelevel is already open')
        stop
     end if
     
     ! filename?
     write(dtime,'(I6.6)') time
     datfile = trim(self%path)//'/'//trim(self%basename)//'_'//trim(dtime)//'_'//self%grid//'.dat'
     call self%cm1log(LOG_MSG, 'open_dataset_time', 'Opening: '//trim(datfile))
     
     ! TODO: check if time is in times for vailidy
     !       Make sure file successfully opens!
     
     ! open dat file
     open(newunit=self%dat_units(1),file=datfile,form='unformatted',access='direct',recl=self%reclen,status='old')
     
     call self%cm1log(LOG_MSG, 'open_dataset_time', 'Timelevel open for reading: '//trim(dtime))
     open_dataset_time = 1
     self%time_open = time
     self%is_time_open = .true.
     
   end function open_dataset_time_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_dataset_time_grads(self) result(close_dataset_time)
     implicit none
     class(cm1_grads) :: self
     
     close_dataset_time = 0
     if ((.not. self%check_dataset_open('close_dataset_time')) .or. &
          (.not. self%check_time_open('close_dataset_time'))) return
     
     close(self%dat_units(1))
     call self%cm1log(LOG_MSG, 'close_dataset_time', 'Timelevel closed for reading')
     close_dataset_time = 1
     self%time_open = -1
     self%is_time_open = .false.
     
   end function close_dataset_time_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2D_at_time_grads(self, varname, Field2D) result(read2D_at_time)
     implicit none
     class(cm1_grads) :: self
     character(len=*), intent(in) :: varname
     real, dimension(self%nx,self%ny) :: Field2D
     integer :: varid,status
     
     read2D_at_time = 0
     if ((.not. self%check_dataset_open('read2D_at_time')) .or. &
          (.not. self%check_time_open('read2D_at_time'))) return
     
     ! Does the variable exist in this dataset?
     varid = self%get_var_byname(varname)
     if (varid == 0) then
        call self%cm1log(LOG_MSG, 'read2D_at_time', 'Variable not found: '//trim(varname))
        return
     end if
     
     call self%cm1log(LOG_INFO, 'read2D_at_time', 'Reading: '//trim(varname))
     ! Read the variable from the dataset
     status = self%read3DXYSlice(varid, 0, Field2D(:,:))
     
     read2D_at_time = status
     
   end function read2D_at_time_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3D_at_time_grads(self, varname, Field3D) result(read3D_at_time)
     implicit none
     class(cm1_grads) :: self
     character(len=*), intent(in) :: varname
     real, dimension(self%nx,self%ny,self%nz) :: Field3D
     integer :: k,varid,status
     
     ! Is the dataset open (has the ctl file been scanned)
     ! Is the dataset open for reading (is the dat file open)
     read3D_at_time = 0
     if ((.not. self%check_dataset_open('read3D_at_time')) .or. &
          (.not. self%check_time_open('read3D_at_time'))) return
     
     ! Does the variable exist in this dataset?
     varid = self%get_var_byname(varname)
     if (varid.eq.0) then
        call self%cm1log(LOG_WARN, 'read3D_at_time', 'Variable not found: '//trim(varname))
        return
     end if
     
     call self%cm1log(LOG_INFO, 'read3D_at_time', 'Reading: '//trim(varname))
     ! Read the variable from the dataset
     do k = 1,self%nz
        status = self%read3DXYSlice(varid, k, Field3D(:,:,k))
        if (status == 0) then
           call self%cm1log(LOG_ERROR, 'read3D_at_time', 'Error reading variable, aborting')
           return
        end if
     end do
     
     read3D_at_time = status
     
   end function read3D_at_time_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ingest_cm1_grads

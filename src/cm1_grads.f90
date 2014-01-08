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
         procedure, public ,pass(self) :: readMultStart => readMultStart_grads
         procedure, public ,pass(self) :: readMultStop => readMultStop_grads
         procedure, public ,pass(self) :: read3DMult => read3DMult_grads
         procedure, public ,pass(self) :: read2DMult => read2DMult_grads

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

      if (self%isopen) then
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
      call self%cm1log(LOG_INFO, 'open_cm1', 'Grid ('//self%grid//') selected.')

      call self%cm1log(LOG_INFO, 'open_cm1', 'Reading GRADS control file.')
      open_cm1 = self%read_ctl()
      self%nunits = 1
      allocate(self%dat_units(self%nunits))
      self%isopen = .true.

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
      read(self%ctl_unit,*)
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
      call self%cm1log(LOG_MSG, 'read_ctl', trim(output))

      read(self%ctl_unit,*) tmp, self%ny  ! for stretched grids only.  detect and calc linear.
      allocate(self%y(self%ny))
      do j = 1, self%ny
         read(self%ctl_unit,*) self%y(j)
      end do
      write (output,500) 'Y',self%ny
      call self%cm1log(LOG_MSG, 'read_ctl', trim(output))

      read(self%ctl_unit,*) tmp, self%nz  ! for stretched grids only.  detect and calc linear.
      allocate(self%z(self%nz))
      do k = 1, self%nz
         read(self%ctl_unit,*) self%z(k)
      end do
      write (output,500) 'Z',self%nz
      call self%cm1log(LOG_MSG, 'read_ctl', trim(output))

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
      call self%cm1log(LOG_MSG, 'read_ctl', trim(output))

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
      call self%cm1log(LOG_MSG, 'read_ctl', trim(output))


      ! This determines the recl to pass to open()
      allocate(recltest(self%nx,self%ny))
      inquire(iolength=self%reclen) recltest
      deallocate(recltest)
      504 format('... Using record length = ', I10)
      write (output,504) self%reclen
      call self%cm1log(LOG_INFO, 'read_ctl', trim(output))

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
      if (.not. self%check_open('close_cm1')) then
         close_cm1 = 0
         return
      end if

      ! If data files are open, close them
      if (self%ismult) then
         close_cm1 = self%readMultStop()
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

      self%isopen = .false.
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

      if (.not. self%check_open('read3DXYSlice')) then
         call self%cm1log(LOG_ERROR, 'read3DXYSlice', 'No datasef open, aborting')
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

   integer function readMultStart_grads(self,time) result(readMultStart)
      implicit none
      class(cm1_grads)          :: self
      integer, intent(in) :: time
      character(len=6) :: dtime
      character(len=256) :: datfile

      if (.not. self%check_open('read3DMultStart')) then
         readMultStart = 0
         return
      end if

      if (self%ismult) then
         call self%cm1log(LOG_WARN, 'read3DMultStart', 'Multiread already started, aborting')
         readMultStart = 0
         return
      end if

      505 format(I6.6)
      ! filename?
      write(dtime,505) time
      datfile = trim(self%path)//'/'//trim(self%basename)//'_'//trim(dtime)//'_'//self%grid//'.dat'
      call self%cm1log(LOG_MSG, 'read3DMultStart', 'Opening: '//trim(datfile))

      ! TODO: check if time is in times for vailidy
      !       Make sure file successfully opens!

      ! open dat file
      open(newunit=self%dat_units(1),file=datfile,form='unformatted',access='direct',recl=self%reclen,status='old')

      call self%cm1log(LOG_INFO, 'read3DMultStart', 'Multiread started for time: '//trim(dtime))
      readMultStart = 1
      self%ismult = .true.

   end function readMultStart_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStop_grads(self) result(readMultStop)
      implicit none
      class(cm1_grads) :: self

      if ((.not. self%check_open('readMultStop')) .or. (.not. self%check_mult('readMultStop'))) then
         readMultStop = 0
         return
      end if

      close(self%dat_units(1))
      call self%cm1log(LOG_INFO, 'read3DMultStop', 'Multiread stopped.')
      readMultStop = 1
      self%ismult = .false.

   end function readMultStop_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2DMult_grads(self, varname, Field2D) result(read2DMult)
      implicit none
      class(cm1_grads) :: self
      character(len=*), intent(in) :: varname
      real, dimension(self%nx,self%ny) :: Field2D
      integer :: varid,status

      if ((.not. self%check_open('read2DMult')) .or. (.not. self%check_mult('read2DMult'))) then
         read2DMult = 0
         return
      end if

      ! Does the variable exist in this dataset?
      varid = self%getVarByName(varname)
      if (varid.eq.0) then
         call self%cm1log(LOG_WARN, 'read2DMult', 'Variable not found: '//trim(varname))
         read2DMult = 0
         return
      end if

      call self%cm1log(LOG_INFO, 'read2DMult', 'Reading: '//trim(varname))
      ! Read the variable from the dataset
      status = self%read3DXYSlice(varid, 0, Field2D(:,:))

      read2DMult = 1

   end function read2DMult_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DMult_grads(self, varname, Field3D) result(read3DMult)
      implicit none
      class(cm1_grads) :: self
      character(len=*), intent(in) :: varname
      real, dimension(self%nx,self%ny,self%nz) :: Field3D
      integer :: k,varid,status

      ! Is the dataset open (has the ctl file been scanned)
      ! Is the dataset open for reading (is the dat file open)
      if ((.not. self%check_open('read3DMult')) .or. (.not. self%check_mult('read3DMult'))) then
         read3DMult = 0
         return
      end if

      ! Does the variable exist in this dataset?
      varid = self%getVarByName(varname)
      if (varid.eq.0) then
         call self%cm1log(LOG_WARN, 'read3DMult', 'Variable not found: '//trim(varname))
         read3DMult = 0
         return
      end if

      call self%cm1log(LOG_INFO, 'read3DMult', 'Reading: '//trim(varname))
      ! Read the variable from the dataset
      do k = 1,self%nz
         status = self%read3DXYSlice(varid, k, Field3D(:,:,k))
      end do

      read3DMult = 1

   end function read3DMult_grads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ingest_cm1_grads

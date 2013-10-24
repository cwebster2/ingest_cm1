!==========================================================================!
!                                                                          !
!   Ingest_CM1                                                             !
!                                                                          !
!                                                                          !
!   A Fortran 2003 module to access model output by the CM1 cloud model    !
!   by George Bryan.                                                       !
!                                                                          !
!   Formats supported : grads style flat files                             !
!                       - multi file per time step (MPI)                   !
!                       - single file per time step                        !
!                                                                          !
!                       Requirement: T dimension values must conincide     !
!                                    with filenames.                       !
!                                                                          !
!   Future formats    : grads style flat files                             !
!                       - all timesteps in one file                        !
!                                                                          !
!                       netcdf / hdf5 files                                !
!                                                                          !
!==========================================================================!



module ingest_cm1

   implicit none

   ! module variables
   private
   public :: variable, cm1

! Parameters to specify dataset type
   integer, public, parameter :: GRADS = 3
   integer, public, parameter :: GRADSMPI = 4
   integer, public, parameter :: HDF = 5

   type variable
      character(len=32) :: varname
      integer           :: levs
      character(len=64) :: units
      integer           :: dims
   end type variable

!TODO: Split these into multipe types? Or keep single type?
   type cm1
      private
      character(len=128) :: path
      character(len=64)  :: basename
      integer            :: dtype
      integer            :: nx, ny, nz, nt, nv, dt
      real, dimension(:), allocatable    :: x, y, z
      integer, dimension(:), allocatable :: times
      type (variable), dimension(:), allocatable :: vars
      integer            :: t0
      logical            :: isopen = .false., ismult = .false.

      ! GRADS stuff
      integer            :: n2d, reclen
      character          :: grid
      integer            :: ctl_unit
      integer,dimension(:),allocatable :: dat_units
      integer            :: nunits

      ! GRADSMPI stuff
      integer            :: ni, nj
      integer            :: nodex, nodey
      integer            :: mpireclen


      contains
         private

         procedure, pass(self) :: cm1_set_nodes
         procedure, pass(self) :: read_ctl
         procedure, pass(self) :: read3DXYSlice
         procedure, pass(self) :: check_open
         procedure, pass(self) :: check_mult
         procedure, pass(self) :: cm1log
!TODO: TYPE BOUND PROCEDURE HERE
         procedure, public ,pass(self) :: open_cm1
         procedure, public ,pass(self) :: close_cm1
         procedure, public ,pass(self) :: getVarByName
         procedure, public ,pass(self) :: readMultStart
         procedure, public ,pass(self) :: readMultStop
         procedure, public ,pass(self) :: read3DMult
         procedure, public ,pass(self) :: read2DMult
         procedure, public ,pass(self) :: read3D
         procedure, public ,pass(self) :: cm1_nx
         procedure, public ,pass(self) :: cm1_ny
         procedure, public ,pass(self) :: cm1_nz
         procedure, public ,pass(self) :: cm1_nt
         procedure, public ,pass(self) :: cm1_nv
         procedure, public ,pass(self) :: cm1_x
         procedure, public ,pass(self) :: cm1_y
         procedure, public ,pass(self) :: cm1_z
         procedure, public ,pass(self) :: cm1_t

!         final :: close_cm1


   end type cm1

   !for logging
   integer, parameter :: LOG_ERROR = 1001
   integer, parameter :: LOG_WARN  = 1002
   integer, parameter :: LOG_INFO  = 1003
   integer, parameter :: LOG_MSG   = 1004

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_cm1(self, dsetpath, dsetbasename, dsettype, grid, nodex, nodey)
      implicit none
      class(cm1) :: self
      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in)  :: dsetbasename
      integer, intent(in)            :: dsettype
      character, optional :: grid
      integer, optional :: nodex, nodey

      if (self%isopen) then
         call cm1log(self, LOG_WARN, 'open_cm1', 'Already open, aborting')
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
      call cm1log(self, LOG_INFO, 'open_cm1', 'Grid ('//self%grid//') selected.')

      select case (self%dtype)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        case (GRADS)
          call cm1log(self, LOG_INFO, 'open_cm1', 'Reading GRADS control file.')
          open_cm1 = self%read_ctl()
          self%nunits = 1
          allocate(self%dat_units(self%nunits))
          self%isopen = 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        case (GRADSMPI)
          if ((.not.present(nodex)) .or. (.not.present(nodey))) then
            call cm1log(self, LOG_ERROR, 'open_cm1', 'GRADSMPI open requires nodex and nodey parameters')
            open_cm1 = 0
            return
          endif

          call cm1log(self, LOG_INFO, 'open_cm1', 'Reading GRADS control file.')
          open_cm1 = self%read_ctl()
          self%nunits = self%cm1_set_nodes(nodex,nodey)
          allocate(self%dat_units(self%nunits))
          self%isopen = 1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        case (HDF)
          call cm1log(self, LOG_ERROR, 'open_cm1', 'HDF ingest not implemented.')
          stop
        case default
          call cm1log(self, LOG_ERROR, 'open_cm1', 'unknown dset type.')
          stop
      end select

   end function open_cm1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read_ctl(self)
      implicit none
      class(cm1) :: self
      character(len=256) :: dset
      character(len=128) :: tmp
      integer            :: i,j,k
      real, allocatable, dimension(:,:) :: recltest

      dset = trim(self%path)//'/'//trim(self%basename)//'_'//self%grid//'.ctl'
      print *, ('[ingest_cm1::open_cm1]: Opening ',trim(dset))
      open(newunit=self%ctl_unit, file=dset, status='old')
      read(self%ctl_unit,*)
      read(self%ctl_unit,*)
      read(self%ctl_unit,*)
      read(self%ctl_unit,*)

      500 format(' [ingest_cm1::open_cm1]: ... Found ',A1,' dimension with ',i5,' points')
      read(self%ctl_unit,*) tmp, self%nx  ! for stretched grids only.  detect and calc linear.
      allocate(self%x(self%nx))
      do i = 1, self%nx
         read(self%ctl_unit,*) self%x(i)
      end do
      print 500,'X',self%nx

      read(self%ctl_unit,*) tmp, self%ny  ! for stretched grids only.  detect and calc linear.
      allocate(self%y(self%ny))
      do j = 1, self%ny
         read(self%ctl_unit,*) self%y(j)
      end do
      print 500,'Y',self%ny

      read(self%ctl_unit,*) tmp, self%nz  ! for stretched grids only.  detect and calc linear.
      allocate(self%z(self%nz))
      do k = 1, self%nz
         read(self%ctl_unit,*) self%z(k)
      end do
      print 500,'Z',self%nz

      ! calculate timelevels.  note - supports a specific timeformat
      ! and only an integer dt (minimum timestep 1 s).
      ! tdef     NT linear hh:mmZddMMMYYYY SSSSSYR
      ! where NT is number of timelevels
      !       YYYY is the starting time in seconds
      !       SSSSS is the timestep in seconds
      !
      ! Modify cm1 to write your ctl file this way, or hand edit them
      read(self%ctl_unit,501) self%nt, self%t0, self%dt
      501 format('tdef ',I10,' linear ',11x,I4,1x,I5,'YR')
      allocate(self%times(self%nt))
      do i = 1, self%nt
         self%times(i) = self%t0 + ((i-1)*self%dt)
      end do
      print 500,'T',self%nt

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
      503 format(' [ingest_cm1::open_cm1]: ... Found ',i3,' variables')
      print 503,self%nv

      print *, ('[ingest_cm1::open_cm1]: Closing Grads control file',trim(dset))

      ! This determines the recl to pass to open()
      allocate(recltest(self%nx,self%ny))
      inquire(iolength=self%reclen) recltest
      deallocate(recltest)
      print *,'[ingest_cm1::open_cm1]: Using record length = ', self%reclen

      close(self%ctl_unit)
      read_ctl = 1

   end function read_ctl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_cm1(self)
      implicit none
      class(cm1) :: self

      ! Check if the dataset is open
      if (.not. check_open(self,'close_cm1')) then
         close_cm1 = 0
         return
      end if

      ! If data files are open, close them
      select case(self%dtype)
         case(GRADS,GRADSMPI)
            if (self%ismult) then
               close_cm1 = self%readMultStop()
            end if
      end select

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

      self%isopen = 0
      close_cm1 = 1
      print *,'[ingest_cm1::close_cm1]: Dataset closed'

   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function getVarByName(self,varname)
      implicit none
      class(cm1)                   :: self
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

   integer function read3DXYSlice(self, varid, level, slice)
      implicit none
      class(cm1)             :: self
      integer, intent(in)    :: varid, level
      real, dimension(self%nx,self%ny) :: slice
      integer :: idx
      real, allocatable, dimension(:,:) :: nodeslice
      integer :: nf, si, sj, i, j

      if (.not. self%check_open('read3DXYSlice')) then
         print *,'[ingest_cm1::read3DXYSlice]: No dataset open, aborting'
         read3DXYSlice = 0
         return
      end if

      select case (self%dtype)
         case (GRADS, GRADSMPI)
            ! NOTE this offset is for one timestep per file!
            ! Calculate record index for 2D or 3D field
            if (self%vars(varid)%levs == 0) then ! 2D
               idx = varid-1
            else ! 3D
               idx = (self%n2d+(varid-1-self%n2d)*(self%nz)+level)
            endif
      end select

      select case (self%dtype)
         case (GRADS)
            ! Read Slice
            read(self%dat_units(1), rec=idx) slice
            read3DXYSlice = 1

         case (GRADSMPI)
            allocate(nodeslice(self%ni,self%nj))
            do nf=0, self%nodex*self%nodey-1

              nodeslice = 0.0
              read(self%dat_units(nf+1), rec=idx) nodeslice
              sj = nf / self%nodex + 1
              si = nf - (sj-1)*self%nodex + 1

              select case(self%grid)
                 case ('s', 'w')
                   do j=1,self%nj
                   do i=1,self%ni
                     slice((si-1)*self%ni+i,(sj-1)*self%nj+j) = nodeslice(i,j)
                   end do
                   end do
                 case ('u')
                   do j=1,self%nj
                   do i=1,self%ni
                     slice((si-1)*(self%ni-1)+i,(sj-1)*self%nj+j) = nodeslice(i,j)
                   end do
                   end do
                 case ('v')
                   do j=1,self%nj
                   do i=1,self%ni
                     slice((si-1)*self%ni+i,(sj-1)*(self%nj-1)+j) = nodeslice(i,j)
                   end do
                   end do
                 case default
                   print *,'[ingest_cm1::read3DXYSlice]: Unsupported MPI grid, something bad happened.'
              end select

            end do
            deallocate(nodeslice)
            read3DXYSlice = 1

         case (HDF)
            print*,'Not implemented'

      end select



   end function read3DXYSlice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStart(self,time)
      implicit none
      class(cm1)          :: self
      integer, intent(in) :: time
      character(len=6) :: dtime
      character(len=256) :: datfile
      integer :: fu
      character(len=6) :: node

      1001 format('[ingest_cm1::readMultStart]: ',A,A,A,A,A)

      if (.not. self%check_open('read3DMultStart')) then
         readMultStart = 0
         return
      end if

      if (self%ismult) then
         call cm1log(self, LOG_WARN, 'read3DMultStart', 'Multiread already started, aborting')
         readMultStart = 0
         return
      end if

      505 format(I6.6)
      select case (self%dtype)
         case (GRADS)
            ! filename?
            write(dtime,505) time
            datfile = trim(self%path)//'/'//trim(self%basename)//'_'//trim(dtime)//'_'//self%grid//'.dat'
            print *, ('[ingest_cm1::read3DMultStart]: Opening ',trim(datfile))

            ! TODO: check if time is in times for vailidy
            !       Make sure file successfully opens!

            ! open dat file
            open(newunit=self%dat_units(1),file=datfile,form='unformatted',access='direct',recl=self%reclen,status='old')

            print *,'[ingest_cm1::read3DMultStart]: Multiread started for time ',time
            readMultStart = 1
            self%ismult = .true.

         case (GRADSMPI)
            if (self%nodex == 0 .or. self%nodey == 0) then
               print *,'[ingest_cm1::read3DMultStart]: Need to set nodex and nodey'
               readMultStart = 0
               self%ismult = .false.
               return
            endif
            write(dtime,505) time
            do fu=0, self%nodex*self%nodey-1
               write(node,505) fu
               datfile = trim(self%path)//'/'//trim(self%basename)//'_'//trim(node)//'_'//trim(dtime)//'_'//self%grid//'.dat'
               call cm1log(self, LOG_MSG, 'read3DMultStart', 'Opening: '//trim(datfile))

               ! open dat file
               open(newunit=self%dat_units(fu+1),file=datfile,form='unformatted',access='direct',recl=self%mpireclen,status='old')
            end do
            print *,'[ingest_cm1::read3DMultStart]: Multiread started for time ',time
            readMultStart = 1
            self%ismult = .true.

         case (HDF)
            print*,'Not implemented'

      end select

   end function readMultStart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStop(self)
      implicit none
      class(cm1) :: self
      integer :: fu

      if ((.not. self%check_open('readMultStop')) .or. (.not. self%check_mult('readMultStop'))) then
         readMultStop = 0
         return
      end if

      select case (self%dtype)
         case (GRADS)
            close(self%dat_units(1))
            print *,'[ingest_cm1::read3DMultStop]: Multiread stopped'
            readMultStop = 1
            self%ismult = .false.

         case (GRADSMPI)
            do fu=1, self%nodex*self%nodey
               close(self%dat_units(fu))
            end do
            print *,'[ingest_cm1::read3DMultStop]: Multiread stopped'
            readMultStop = 1
            self%ismult = .false.

         case (HDF)
            print*,'Not implemented'

      end select

   end function readMultStop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2DMult(self, varname, Field2D)
      implicit none
      class(cm1), intent(in) :: self
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
         print *,'[ingest_cm1::read2DMult]: variable not found: ',trim(varname)
         read2DMult = 0
         return
      end if

      print *,'[ingest_cm1::read2DMult]: Reading: ',trim(varname)
      ! Read the variable from the dataset
      status = self%read3DXYSlice(varid, 0, Field2D(:,:))

      read2DMult = 1

   end function read2DMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   integer function read3DMult(varname, Field4D)
!      implicit none
!      character*(*), intent(in) :: varname
!      real, dimension(nx,ny,nz,1) :: Field4D
!      real, dimension(nx,ny,nz) :: Field3D
!      integer :: k,varid,status

!      status = read3dMult(varname, Field3D)
!      Field4D(:,:,:,1) = Field3d(:,:,:)

!      read3DMult = status

!   end function read3DMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer function read3DMult(self, varname, Field3D)
      implicit none
      class(cm1), intent(in) :: self
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
         print *,'[ingest_cm1::read3DMult]: variable not found: ',trim(varname)
         read3DMult = 0
         return
      end if

      print *,'[ingest_cm1::read3DMult]: Reading: ',trim(varname)
      ! Read the variable from the dataset
      do k = 1,self%nz
         status = self%read3DXYSlice(varid, k, Field3D(:,:,k))
      end do

      read3DMult = 1

   end function read3DMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3D(self, varname, time, Field3D)
      implicit none
      class(cm1), intent(in)   :: self
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
         print *,'[ingest_cm1:read3D]: data open failure, aborting'
         read3D = 0
         return
      endif

      status = self%read3DMult(varname, Field3D)
      if (status.eq.0) then
         print *,'[ingest_cm1:read3D]: data read failure, aborting'
         read3D = 0
         status = self%readMultStop()
         return
      endif

      status = self%readMultStop()
      if (status.eq.0) then
         print *,'[ingest_cm1:read3D]: data close failure, aborting'
         read3D = 0
         return
      endif

      read3D = 1
   end function read3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nx(self)
      implicit none
      class(cm1),intent(in) :: self

      if (.not. self%check_open('cm1_nx')) then
         cm1_nx = 0
         return
      end if

      cm1_nx = self%nx
   end function cm1_nx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_set_nodes(self, mpix, mpiy)
      implicit none
      class(cm1) :: self
      integer, intent(in) :: mpix, mpiy
      real, allocatable, dimension(:,:) :: recltest

!      if (.not. self%check_open('cm1_set_nodes')) then
!         cm1_set_nodes = 0
!         return
!      end if

      if (self%dtype /= GRADSMPI) then
         print *,'[ingest_cm1::cm1_set_nodes]: Dataset not of MPI type, aborting'
         cm1_set_nodes = 0
         return
      end if

      self%nodex = mpix
      self%nodey = mpiy
      print *,'[ingest_cm1::cm1_set_nodes]: nodex=',self%nodex,' nodey=',self%nodey

      select case(self%grid)
        case ('s','w')
           self%ni = self%nx/self%nodex
           self%nj = self%ny/self%nodey
        case ('u')
           self%ni = 1+(self%nx-1)/self%nodex
           self%nj = self%ny/self%nodey
        case ('v')
           self%ni = self%nx/self%nodex
           self%nj = 1+(self%ny-1)/self%nodey
        case default
           print *,'Unsupported grid'
           stop
      end select
      print *,'[ingest_cm1::cm1_set_nodes]: ni=',self%ni,' nj=',self%nj

      allocate(recltest(self%ni,self%nj))
      inquire(iolength=self%mpireclen) recltest
      deallocate(recltest)
      print *,'[ingest_cm1::cm1_set_nodes]: mpi recl = ',self%mpireclen

      cm1_set_nodes = self%nodex*self%nodey
   end function cm1_set_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_ny(self)
      implicit none
      class(cm1),intent(in) :: self

      if (.not. self%check_open('cm1_ny')) then
         cm1_ny = 0
         return
      end if

      cm1_ny = self%ny
   end function cm1_ny

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nz(self)
      implicit none
      class(cm1),intent(in) :: self

      if (.not. self%check_open('cm1_nz')) then
         cm1_nz = 0
         return
      end if

      cm1_nz = self%nz
   end function cm1_nz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nt(self)
      implicit none
      class(cm1),intent(in) :: self

      if (.not. self%check_open('cm1_nt')) then
         cm1_nt = 0
         return
      end if

      cm1_nt = self%nt
   end function cm1_nt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_nv(self)
      implicit none
      class(cm1),intent(in) :: self

      if (.not. self%check_open('cm1_nv')) then
         cm1_nv = 0
         return
      end if

      cm1_nv = self%nv
   end function cm1_nv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_dt(self)
      implicit none
      class(cm1),intent(in) :: self

      if (.not. self%check_open('cm1_dt')) then
         cm1_dt = 0
         return
      end if

      cm1_dt = self%dt
   end function cm1_dt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_x(self, cm1x)
      implicit none
      class(cm1),intent(in) :: self
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
      class(cm1),intent(in) :: self
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
      class(cm1),intent(in) :: self
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
      class(cm1),intent(in) :: self
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
      class(cm1),intent(in) :: self
      character(len=*), intent(in) :: funcname

      ! Is the dataset open (has the ctl file been scanned)
      if (.not. self%isopen) then
         call cm1log(self, LOG_ERROR, funcname, 'No dataset open, aborting')
         check_open = .false.
         return
      end if
      check_open = .true.

   end function check_open

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   logical function check_mult(self, funcname)
      use, intrinsic :: iso_fortran_env, only: error_unit
      class(cm1),intent(in) :: self
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
      class(cm1), intent(in) :: self
      integer, intent(in) :: loglevel
      character(len=*), intent(in) :: funcname, message

      character(len=6) :: proto

      1001 format('[ingestcm1',3(':',A),']: Error: ',A)
      1002 format('[ingestcm1',3(':',A),']: Warning: ',A)
      1003 format('[ingestcm1',3(':',A),']: Info: ',A)
      1004 format('[ingestcm1',3(':',A),']: ',A)

      select case (self%dtype)

         case (GRADS)
            proto = 'GRADS'

         case (GRADSMPI)
            proto = 'GRMPI'

         case default
            proto = 'X'

      end select

      select case (loglevel)

         case (LOG_ERROR)
            write(error_unit,1001) trim(proto), self%grid, funcname, message

         case (LOG_WARN)
            write(error_unit,1002) trim(proto), self%grid, funcname, message

         case (LOG_INFO)
            write(error_unit,1003) trim(proto), self%grid, funcname, message

         case (LOG_MSG)
            write(output_unit,1004) trim(proto), self%grid, funcname, message

      end select

   end subroutine

end module ingest_cm1

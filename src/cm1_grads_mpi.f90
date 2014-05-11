!==========================================================================!
!
!   Ingest_CM1_grads_mpi
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

module ingest_cm1_grads_mpi
 
    use ingest_cm1_grads
 
    implicit none
 
    type, extends(cm1_grads) :: cm1_grads_mpi
       !private
          ! GRADSMPI stuff
          integer            :: ni, nj
          integer            :: nodex, nodey
          integer            :: mpireclen
       contains
          !private
 
          !procedure, pass(self) :: read_ctl
          procedure, pass(self) :: cm1_set_nodes => cm1_set_nodes_grads_mpi
          procedure, pass(self) :: read3DXYSlice => read3DXYSlice_grads_mpi
          procedure, public ,pass(self) :: open_cm1 => open_cm1_grads_mpi
          procedure, public ,pass(self) :: readMultStart => readMultStart_grads_mpi
          procedure, public ,pass(self) :: readMultStop => readMultStop_grads_mpi
 
          ! This ensures that the filehandles and arrays are closed if the variable goes out of     scope
          !TODO: finalization needs gfortran 4.9 or ifort
          !!!final :: final_cm1
 
    end type cm1_grads_mpi
 
contains

   integer function open_cm1_grads_mpi(self, dsetpath, dsetbasename, dsettype, &
                                       grid, nodex, nodey, hdfmanage) result(open_cm1)
      implicit none
      class(cm1_grads_mpi) :: self
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

      if ((.not.present(nodex)) .or. (.not.present(nodey))) then
        call self%cm1log(LOG_ERROR, 'open_cm1', 'GRADSMPI open requires nodex and nodey parameters')
        open_cm1 = 0
        return
      endif

      call self%cm1log(LOG_INFO, 'open_cm1', 'Reading GRADS control file.')
      open_cm1 = self%read_ctl()
      self%nunits = self%cm1_set_nodes(nodex,nodey)
      allocate(self%dat_units(self%nunits))
      self%isopen = .true.

   end function open_cm1_grads_mpi


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

   integer function read3DXYSlice_grads_mpi(self, varid, level, slice) result(read3DXYSlice)
      implicit none
      class(cm1_grads_mpi)             :: self
      integer, intent(in)    :: varid, level
      real, dimension(self%nx,self%ny) :: slice
      integer :: idx
      real, allocatable, dimension(:,:) :: nodeslice
      integer :: nf, si, sj, i, j

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
             call self%cm1log(LOG_ERROR, 'read3DXYSlice', 'Unsupported MPI grid, something bad happened.')
        end select

      end do
      deallocate(nodeslice)
      read3DXYSlice = 1

   end function read3DXYSlice_grads_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStart_grads_mpi(self,time) result(readMultStart)
      implicit none
      class(cm1_grads_mpi)          :: self
      integer, intent(in) :: time
      character(len=6) :: dtime
      character(len=256) :: datfile
      integer :: fu
      character(len=6) :: node

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
      if (self%nodex == 0 .or. self%nodey == 0) then
         call self%cm1log(LOG_WARN, 'read3DMultStart', 'Need to set nodex and nodey')
         readMultStart = 0
         self%ismult = .false.
         return
      endif
      write(dtime,505) time
      do fu=0, self%nodex*self%nodey-1
         write(node,505) fu
         datfile = trim(self%path)//'/'//trim(self%basename)//'_'//trim(node)//'_'//trim(dtime)//'_'//self%grid//'.dat'
         call self%cm1log(LOG_MSG, 'read3DMultStart', 'Opening: '//trim(datfile))

         ! open dat file
         open(newunit=self%dat_units(fu+1),file=datfile,form='unformatted',access='direct',recl=self%mpireclen,status='old')
      end do
      call self%cm1log(LOG_INFO, 'read3DMultStart', 'Multiread started for time: '//trim(dtime))
      readMultStart = 1
      self%ismult = .true.

   end function readMultStart_grads_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStop_grads_mpi(self) result(readMultStop)
      implicit none
      class(cm1_grads_mpi) :: self
      integer :: fu

      if ((.not. self%check_open('readMultStop')) .or. (.not. self%check_mult('readMultStop'))) then
         readMultStop = 0
         return
      end if

      do fu=1, self%nodex*self%nodey
         close(self%dat_units(fu))
      end do
      call self%cm1log(LOG_INFO, 'read3DMultStop', 'Multiread stopped.')
      readMultStop = 1
      self%ismult = .false.

   end function readMultStop_grads_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function cm1_set_nodes_grads_mpi(self, mpix, mpiy) result(cm1_set_nodes)
      implicit none
      class(cm1_grads_mpi) :: self
      integer, intent(in) :: mpix, mpiy
      real, allocatable, dimension(:,:) :: recltest
      character(len=128) :: output

!      Don't check open because this is now called automatically before the open flag is set.
!      if (.not. self%check_open('cm1_set_nodes')) then
!         cm1_set_nodes = 0
!         return
!      end if

      101 format('nodex= ',I45,', nodey= ',I5)
      102 format('ni   = ',I45,', nj   = ',I5)
      103 format('Using mpi record length: ',I8)

      if (self%dtype /= GRADSMPI) then
         call self%cm1log(LOG_WARN, 'cm1_set_nodes', 'Dataset not of MPI type, aborting.')
         cm1_set_nodes = 0
         return
      end if

      self%nodex = mpix
      self%nodey = mpiy
      write(output,101) self%nodex, self%nodey
      call self%cm1log(LOG_INFO, 'cm1_set_nodes', trim(output))

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
           call self%cm1log(LOG_ERROR, 'cm1_set_nodes', 'Unsupported grid')
           stop
      end select
      write(output,102) self%ni, self%nj
      call self%cm1log(LOG_INFO, 'cm1_set_nodes', trim(output))

      allocate(recltest(self%ni,self%nj))
      inquire(iolength=self%mpireclen) recltest
      deallocate(recltest)
      write(output,103) recltest
      call self%cm1log(LOG_INFO, 'cm1_set_nodes', trim(output))

      cm1_set_nodes = self%nodex*self%nodey
   end function cm1_set_nodes_grads_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ingest_cm1_grads_mpi

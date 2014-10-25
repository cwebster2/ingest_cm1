!==========================================================================!
!
!   Ingest_CM1_hdf
!
!   Implementation of the ingest_cm1 routines to support reading from hd5
!   datafiles.  Note that this code does not yet read the native cm1hdf5 
!   output and requires un-tiling using a provided python program.
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


#ifdef HDF5
module ingest_cm1_hdf5

   use ingest_cm1_base
   use hdf5

   implicit none

   type, extends(cm1_base) :: cm1_hdf5
      !private
         integer :: h5err
         integer(HID_T) :: file_id
         logical :: manage_hdf
      contains
         !private
         !procedure, pass(self) :: initiate_hdf
         !procedure, pass(self) :: deinitiate_hdf

         procedure, pass(self) :: get_times
         procedure, pass(self) :: scan_hdf
         procedure, public ,pass(self) :: getVarByName => getVarByName_hdf5
         procedure, public ,pass(self) :: open_cm1 => open_cm1_hdf5
         procedure, public ,pass(self) :: close_cm1 => close_cm1_hdf5
         procedure, public ,pass(self) :: readMultStart => readMultStart_hdf5
         procedure, public ,pass(self) :: readMultStop => readMultStop_hdf5
         procedure, public ,pass(self) :: read3DMult => read3DMult_hdf5
         procedure, public ,pass(self) :: read2DMult => read2DMult_hdf5

         procedure, private, pass(self) :: read3DDerived

         !TODO: finalization needs gfortran 4.9 or ifort
         !!!final :: final_cm1
   end type cm1_hdf5

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine initiate_hdf()
      implicit none
      integer :: h5err

!      call cm1log(LOG_DEBUG, 'initiate_hdf', 'Initializing HDF5 library')
      call h5open_f(h5err)
   end subroutine initiate_hdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine deinitiate_hdf()
      implicit none
      integer :: h5err

!      call cm1log(LOG_DEBUG, 'deinitiate_hdf', 'De-initializing HDF5 library')
      call h5close_f(h5err)
   end subroutine deinitiate_hdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_cm1_hdf5(self, dsetpath, dsetbasename, dsettype, & 
                                  grid, nodex, nodey, hdfmanage) &
           result(open_cm1)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in)  :: dsetbasename
      integer, intent(in)            :: dsettype
      character, optional :: grid
      integer, optional :: nodex, nodey
      logical, optional :: hdfmanage


      self%dtype = dsettype
      if (present(grid)) then
        self%grid = grid
      else
        self%grid = 's'
      endif
      call self%cm1log(LOG_MSG, 'open_cm1', 'Opening HDF dataset')
      call self%cm1log(LOG_MSG, 'open_cm1', 'Grid ('//self%grid//') selected.')

      if (self%isopen) then
         call self%cm1log(LOG_WARN, 'open_cm1', 'Already open, aborting')
         open_cm1 = 0
         return
      end if

      if (.not. present(hdfmanage)) then
         self%manage_hdf = .true.
      else
         self%manage_hdf = hdfmanage
      end if
      if (self%manage_hdf) call initiate_hdf

      self%path = dsetpath
      self%basename = dsetbasename


      call self%cm1log(LOG_MSG, 'open_cm1', 'Scanning HDF files.')
      open_cm1 = self%scan_hdf()
      if (open_cm1 .eq. 1) then
         self%isopen = .true.
      else
         self%isopen = .false.
      end if

   end function open_cm1_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TODO: protect this with ifdef for non-gnu compiler?
!subroutine execute_command_line(cmd)
!   character(len=*) :: cmd
!   call system(cmd)
!end subroutine execute_command_line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! this is UGLY
!TODO: is there a better way?

   subroutine get_times(self)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), parameter :: tmpfile = '/tmp/ingest_cm1_flist'
      integer :: u, i, reason
      real :: r

!TODO: re-write this as a call to find, not ls
      call execute_command_line('find '//trim(self%path)//'/'//trim(self%basename)//&
                                "*.h5 -printf '%f\n' | sed -n 's/.*"//trim(self%basename)//&
                                "\.//;s/\.h5//p' > "//trim(tmpfile))
      !print *,'find '//trim(self%path)//'/'//trim(self%basename)//"*.h5 -printf '%f\n' | sed -n 's/.*"//trim(self%basename)//"\.//;s/\.h5//p' > "//trim(tmpfile)

      open(newunit=u, file=tmpfile, action="read")
      i = 0
      do
        read(u,fmt='(A)', iostat=reason) r
        if (reason /= 0) exit
        i = i + 1
      end do

      self%nt = i
      allocate(self%times(self%nt))
      !print *, self%nt

      rewind(u)

      do i=1,self%nt
         read(u,'(I5.5)') self%times(i)
         !print *, self%times(i)
      end do

      close(u)

      self%dt = self%times(2) - self%times(1)

      call execute_command_line('rm '//trim(tmpfile))

   end subroutine get_times

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function get_h5_scalar_int(file_id, dset, err)
    implicit none
    integer(HID_T), intent(inout) :: file_id
    integer, intent(inout) :: err
    character(len=*) :: dset
    integer(HID_T) :: dset_id
    integer(HSIZE_T), dimension(1) :: dims
    integer, dimension(1) :: data_out

    call h5dopen_f(file_id, dset, dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, data_out, dims, err)
    call h5dclose_f(dset_id, err)

    get_h5_scalar_int = data_out(1)

  end function get_h5_scalar_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_h5_1d_float(file_id, dset, err, data_out, data_len)
    implicit none
    integer(HID_T), intent(inout) :: file_id
    integer, intent(inout) :: err, data_len
    character(len=*) :: dset
    integer(HID_T) :: dset_id
    integer(HSIZE_T), dimension(1) :: dims
    real, dimension(data_len), intent(out) :: data_out

    dims(1) = data_len
    call h5dopen_f(file_id, dset, dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, data_out, dims, err)
    call h5dclose_f(dset_id, err)

  end subroutine get_h5_1d_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_h5_2d_float(file_id, dset, err, var, nx, ny)
    implicit none
    integer(HID_T), intent(inout) :: file_id
    integer, intent(inout) :: err, nx, ny
    character(len=*) :: dset
    integer(HID_T) :: dset_id
    integer(HSIZE_T), dimension(2) :: dims
    real, dimension(nx,ny), intent(out) :: var

    dims(1) = nx
    dims(2) = ny
    call h5dopen_f(file_id, dset, dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, var, dims, err)
    call h5dclose_f(dset_id, err)

  end subroutine get_h5_2d_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine get_h5_3d_float(file_id, dset, err, var, nx, ny, nz)
    implicit none
    integer(HID_T), intent(inout) :: file_id
    integer, intent(inout) :: err, nx, ny, nz
    character(len=*) :: dset
    integer(HID_T) :: dset_id
    integer(HSIZE_T), dimension(3) :: dims
    real, dimension(nx,ny,nz), intent(out) :: var

    dims(1) = nx
    dims(2) = ny
    dims(3) = nz
    call h5dopen_f(file_id, dset, dset_id, err)
    call h5dread_f(dset_id, H5T_NATIVE_REAL, var, dims, err)
    call h5dclose_f(dset_id, err)

  end subroutine get_h5_3d_float

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function scan_hdf(self)
    implicit none
    class(cm1_hdf5)     :: self
    integer             :: hdfmetatime
    character(len=256)  :: filename, output
    character(len=5)    :: dtime
    integer(HID_T) :: gid_3d, gid_2d, gid_base
    integer :: i, n3d, n2d, nbase, unused, unused2

    ! GET TIMES
    call self%get_times()
    hdfmetatime = self%times(1)

    500 format(I5.5)
    write(dtime, 500) hdfmetatime

    filename = trim(self%path)//'/'//trim(self%basename)//'.'//dtime//'.h5'

    call self%cm1log(LOG_MSG, 'scan_hdf', 'Scanning metadata in '//trim(filename))
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, self%file_id, self%h5err)

    if (self%h5err .eq. -1) then
       call self%cm1log(LOG_ERROR, 'scan_hdf', 'Unable to open: '//trim(filename))
       scan_hdf = 0
       return
    end if

    self%nx = get_h5_scalar_int(self%file_id, '/grid/nx', self%h5err)
    self%ny = get_h5_scalar_int(self%file_id, '/grid/ny', self%h5err)
    self%nz = get_h5_scalar_int(self%file_id, '/grid/nz', self%h5err)

    select case(self%grid)
       case ('u')
         self%nx = self%nx + 1
       case ('v')
         self%ny = self%ny + 1
       case ('w')
         self%nz = self%nz + 1
    end select

    allocate(self%x(self%nx))
    allocate(self%y(self%ny))
    allocate(self%z(self%nz))

    select case(self%grid)
       case ('s')
         call get_h5_1d_float(self%file_id, '/mesh/xh', self%h5err, self%x, self%nx)
         call get_h5_1d_float(self%file_id, '/mesh/yh', self%h5err, self%y, self%ny)
         call get_h5_1d_float(self%file_id, '/mesh/zh', self%h5err, self%z, self%nz)
       case ('u')
         call get_h5_1d_float(self%file_id, '/mesh/xf', self%h5err, self%x, self%nx)
         call get_h5_1d_float(self%file_id, '/mesh/yh', self%h5err, self%y, self%ny)
         call get_h5_1d_float(self%file_id, '/mesh/zh', self%h5err, self%z, self%nz)
       case ('v')
         call get_h5_1d_float(self%file_id, '/mesh/xh', self%h5err, self%x, self%nx)
         call get_h5_1d_float(self%file_id, '/mesh/yf', self%h5err, self%y, self%ny)
         call get_h5_1d_float(self%file_id, '/mesh/zh', self%h5err, self%z, self%nz)
       case ('w')
         call get_h5_1d_float(self%file_id, '/mesh/xh', self%h5err, self%x, self%nx)
         call get_h5_1d_float(self%file_id, '/mesh/yh', self%h5err, self%y, self%ny)
         call get_h5_1d_float(self%file_id, '/mesh/zf', self%h5err, self%z, self%nz)
    end select

    ! GET VARS
!TODO: read units / description attributes
    select case(self%grid)
       case ('s')
          ! scan 3d, 2d, and basestate variables
          call h5gopen_f(self%file_id, '/3d_s', gid_3d, self%h5err)
          call h5gopen_f(self%file_id, '/2d', gid_2d, self%h5err)
          call h5gopen_f(self%file_id, '/basestate', gid_base, self%h5err)

          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          call h5gget_info_f(gid_2d, unused, n2d, unused2, self%h5err)
          call h5gget_info_f(gid_base, unused, nbase, unused2, self%h5err)

          self%nv = n3d + n2d + nbase
          allocate(self%vars(self%nv))

          do i = 1,n3d
             call h5gget_obj_info_idx_f(gid_3d, '/3d_s', i-1, self%vars(i)%varname, unused, self%h5err)
             self%vars(i)%levs = self%nz
             self%vars(i)%dims = 3
          end do
          do i = 1,n2d
             call h5gget_obj_info_idx_f(gid_2d, '/2d', i-1, self%vars(i+n3d)%varname, unused, self%h5err)
             self%vars(i+n3d)%levs = 1
             self%vars(i+n3d)%dims = 2
          end do
          do i = 1,nbase
             call h5gget_obj_info_idx_f(gid_base, '/basestate', i-1, self%vars(i+n3d+n2d)%varname, unused, self%h5err)
             self%vars(i+n2d+n3d)%levs = self%nz
             self%vars(i+n2d+n3d)%dims = 1
          end do

          call h5gclose_f(gid_3d, self%h5err)
          call h5gclose_f(gid_2d, self%h5err)
          call h5gclose_f(gid_base, self%h5err)

       case ('u')
          call h5gopen_f(self%file_id, '/3d_u', gid_3d, self%h5err)
          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          self%nv = n3d
          allocate(self%vars(self%nv))
          do i = 1,n3d
             call h5gget_obj_info_idx_f(gid_3d, '/3d_u', i-1, self%vars(i)%varname, unused, self%h5err)
             self%vars(i)%levs = self%nz
             self%vars(i)%dims = 3
          end do
          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          call h5gclose_f(gid_3d, self%h5err)

       case ('v')
          call h5gopen_f(self%file_id, '/3d_v', gid_3d, self%h5err)
          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          self%nv = n3d
          allocate(self%vars(self%nv))
          do i = 1,n3d
             call h5gget_obj_info_idx_f(gid_3d, '/3d_v', i-1, self%vars(i)%varname, unused, self%h5err)
             self%vars(i)%levs = self%nz
             self%vars(i)%dims = 3
          end do
          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          call h5gclose_f(gid_3d, self%h5err)

       case ('w')
          call h5gopen_f(self%file_id, '/3d_w', gid_3d, self%h5err)
          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          self%nv = n3d
          allocate(self%vars(self%nv))
          do i = 1,n3d
             call h5gget_obj_info_idx_f(gid_3d, '/3d_w', i-1, self%vars(i)%varname, unused, self%h5err)
             self%vars(i)%levs = self%nz
             self%vars(i)%dims = 3
          end do
          call h5gget_info_f(gid_3d, unused, n3d, unused2, self%h5err)
          call h5gclose_f(gid_3d, self%h5err)
!TODO: need default case

    end select

!    do i=1,self%nv
!       print *,i,self%vars(i)%varname, self%vars(i)%levs, self%vars(i)%dims
!    end do


    501 format('... Found ',A1,' dimension with ',i5,' points')
    503 format('... Found ',i3,' variables')

    write (output,501) 'X',self%nx
    call self%cm1log(LOG_INFO, 'scan_hdf', trim(output))
    write (output,501) 'Y',self%ny
    call self%cm1log(LOG_INFO, 'scan_hdf', trim(output))
    write (output,501) 'Z',self%nz
    call self%cm1log(LOG_INFO, 'scan_hdf', trim(output))
    write (output,501) 'T',self%nt
    call self%cm1log(LOG_INFO, 'scan_hdf', trim(output))
    write (output,503) self%nv
    call self%cm1log(LOG_INFO, 'scan_hdf', trim(output))

    scan_hdf = 1
    call h5fclose_f(self%file_id, self%h5err)

  end function scan_hdf

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

   integer function close_cm1_hdf5(self) result(close_cm1)
      implicit none
      class(cm1_hdf5) :: self

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
      self%nx = 0
      self%ny = 0
      self%nz = 0
      self%nt = 0
      self%nv = 0

      self%isopen = .false.
      if (self%manage_hdf) call deinitiate_hdf
      close_cm1 = 1
      call self%cm1log(LOG_MSG, 'close_cm1', 'Dataset closed')

   end function close_cm1_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStart_hdf5(self,time) result(readMultStart)
      implicit none
      class(cm1_hdf5)     :: self
      integer, intent(in) :: time
      character(len=5) :: dtime
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

      500 format(I5.5)
      write(dtime, 500) time

      datfile = trim(self%path)//'/'//trim(self%basename)//'.'//dtime//'.h5'

      call self%cm1log(LOG_MSG, 'read3DMultStart', 'Opening for Multiread: '//trim(datfile))
      call h5fopen_f(datfile, H5F_ACC_RDONLY_F, self%file_id, self%h5err)

      self%ismult = .true.
      readMultStart = 1

   end function readMultStart_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStop_hdf5(self) result(readMultStop)
      implicit none
      class(cm1_hdf5) :: self

      if ((.not. self%check_open('readMultStop')) .or. (.not. self%check_mult('readMultStop'))) then
         readMultStop = 0
         return
      end if

      call h5fclose_f(self%file_id, self%h5err)
      call self%cm1log(LOG_MSG, 'read3DMultStop', 'Closed file for Multiread')

      self%ismult = .false.
      readMultStop = 1

   end function readMultStop_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer function getVarByName_hdf5(self,varname) result(getVarByName)
       implicit none
       class(cm1_hdf5)                   :: self
       character(len=*), intent(in) :: varname
       integer :: i
 
       if (.not. self%check_open('getVarByName')) then
          getVarByName = 0
          return
       end if
 
       do i = 1, self%nv
          if (varname.eq.self%vars(i)%varname) then
             getVarByName = self%vars(i)%dims
             return
          endif
       end do
 
       getVarByName = 0
 
    end function getVarByName_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2DMult_hdf5(self, varname, Field2D) result(read2DMult)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), intent(in) :: varname
      real, dimension(self%nx,self%ny) :: Field2D
      integer :: varid

      if ((.not. self%check_open('read2DMult')) .or. (.not. self%check_mult('read2DMult'))) then
         read2DMult = 0
         return
      end if

      ! Does the variable exist in this dataset?
      varid = self%getVarByName(varname)
      if (varid.ne.2) then
         call self%cm1log(LOG_WARN, 'read2DMult', 'Variable not found: '//trim(varname))
         read2DMult = 0
         return
      end if

      call self%cm1log(LOG_INFO, 'read2DMult', 'Reading: '//trim(varname))
      call get_h5_2d_float(self%file_id, '/2d/'//trim(varname), self%h5err, Field2D, self%nx, self%ny)

      read2DMult = 1

   end function read2DMult_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DMult_hdf5(self, varname, Field3D) result(read3DMult)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), intent(in) :: varname
      real, dimension(self%nx,self%ny,self%nz) :: Field3D
      integer :: varid

      ! Is the dataset open (has the ctl file been scanned)
      ! Is the dataset open for reading (is the dat file open)
      if ((.not. self%check_open('read3DMult')) .or. (.not. self%check_mult('read3DMult'))) then
         read3DMult = 0
         return
      end if

      ! Does the variable exist in this dataset?
      varid = self%getVarByName(varname)
      if (varid.ne.3) then
         call self%cm1log(LOG_WARN, 'read3DMult', 'Variable not found: '//trim(varname))
         ! lets see if this variable can be reconstructed before giving up.
         read3DMult = self%read3DDerived(varname, Field3D)
         return
      end if

      call self%cm1log(LOG_INFO, 'read3DMult', 'Reading: '//trim(varname))
      call get_h5_3d_float(self%file_id, '/3d_'//self%grid//'/'//trim(varname), self%h5err, Field3D, self%nx, self%ny, self%nz)

      read3DMult = 1

   end function read3DMult_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DDerived(self, varname, Field3D) result(read3D)
     implicit none
     class(cm1_hdf5) :: self
     character(len=*), intent(in) :: varname
     real, dimension(self%nx, self%ny, self%nz), intent(inout) :: Field3D
     
     integer :: varid0, varid
     character(len=128) :: varname0, varname_pert
     real, dimension(self%nz) :: Field0
     integer :: k

     ! Is the dataset open (has the ctl file been scanned)
     ! Is the dataset open for reading (is the dat file open)
     if ((.not. self%check_open('read3DDerived')) .or. (.not. self%check_mult('read3DDerived'))) then
        read3D = 0
        return
     end if

     ! This function reads a basestate variable and a perturbation variable to
     ! return the requested variable.  e.g. p = p_0 + p', read p_0 and p' to 
     ! return p

     varname0 = trim(varname)//'0'
     varname_pert = trim(varname)//'pert'

     ! fix some quirky names
     if (varname0 .eq. 'p0') varname0 = 'pres0'
     if (varname0 .eq. 'rho0') varname0 = 'rh0'

     ! Verify variables exist. Returns 0 for non-existance or the number of dimensions if found
     varid0 = self%getVarByName(varname0)
     varid = self%getVarByName(varname_pert)
     if (varid0.ne.1) then ! we are expecting a 1-D variable
        call self%cm1log(LOG_WARN, 'read3DDerived', 'Variable not found: '//trim(varname0))
        read3D = 0
        return
     end if
     if (varid.ne.3) then ! we are expecting a 3-D variable
        call self%cm1log(LOG_WARN, 'read3DDerived', 'Variable not found: '//trim(varname_pert))
        read3D = 0
        return
     end if

     ! Get perturbation field into 3d variable
     call self%cm1log(LOG_INFO, 'read3DDerived', 'Reading: '//trim(varname_pert))
     call get_h5_3d_float(self%file_id, '/3d_'//self%grid//'/'//trim(varname_pert), self%h5err, Field3D, self%nx, self%ny, self%nz)

     ! Get basestate into 1D variable
     call self%cm1log(LOG_INFO, 'read3DDerived', 'Reading: '//trim(varname0))
     call get_h5_1d_float(self%file_id, '/basestate/'//trim(varname0), self%h5err, Field0, self%nz)

     ! Add them
     !$omp parallel do private(k) shared(Field3D, Field0)
     do k = 1, self%nz
        Field3D(:,:,k) = Field3D(:,:,k) + Field0(k)
     end do

     read3D = 1

   end function read3DDerived

end module ingest_cm1_hdf5
#endif

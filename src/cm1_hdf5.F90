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
         procedure, pass(self) :: get_times
         procedure, pass(self) :: scan_hdf
         procedure, public ,pass(self) :: get_var_byname => get_var_byname_hdf5
         procedure, public ,pass(self) :: open_cm1 => open_cm1_hdf5
         procedure, public ,pass(self) :: close_cm1 => close_cm1_hdf5
         procedure, public ,pass(self) :: open_dataset_time => open_dataset_time_hdf5
         procedure, public ,pass(self) :: close_dataset_time => close_dataset_time_hdf5
         procedure, public ,pass(self) :: read3D_at_time => read3D_at_time_hdf5
         procedure, public ,pass(self) :: read2D_at_time => read2D_at_time_hdf5

         procedure, private, pass(self) :: read3D_derived
!         procedure, pass(self) :: final_cm1
         
         !TODO: abstract this into a makefile.in macro or something (autotools?)
         !TODO: require a minimum version of ifort
#if (defined(__GFORTRAN__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR > 8))) || defined(__INTEL_COMPILER)
         final :: final_cm1
#endif
   end type cm1_hdf5

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This function is split from the main HDF code to handle the case that this
  ! library is used within a program that also makes use of HDF.  The HDF library
  ! is not friendly to multiple open/close operations, so we only want to call
  ! them if the program calling the library allows it.
  
   subroutine initiate_hdf()
      implicit none
      integer :: h5err

      call h5open_f(h5err)
   end subroutine initiate_hdf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This function is split from the main HDF code to handle the case that this
   ! library is used within a program that also makes use of HDF.  The HDF library
   ! is not friendly to multiple open/close operations, so we only want to call
   ! them if the program calling the library allows it.
 
   subroutine deinitiate_hdf()
      implicit none
      integer :: h5err

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

      if (self%is_dataset_open) then
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
         self%is_dataset_open = .true.
      else
         self%is_dataset_open = .false.
      end if

   end function open_cm1_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TODO: protect this with ifdef for non-gnu compiler?
!subroutine execute_command_line(cmd)
!   character(len=*) :: cmd
!   call system(cmd)
!end subroutine execute_command_line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!TODO: randomize the filename with C's "char* mkdtmp(char*)"  

   subroutine get_times(self)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), parameter :: tmpfile = '/tmp/ingest_cm1_flist'
      integer :: u, i, reason
      real :: r

      call execute_command_line('find '//trim(self%path)//'/'//trim(self%basename)//&
                                "*.h5 -printf '%f\n' | sed -n 's/.*"//trim(self%basename)//&
                                "\.//;s/\.h5//p' > "//trim(tmpfile))

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
! https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58175

  subroutine final_cm1(cm1hdf5)
    implicit none
    type(cm1_hdf5) :: cm1hdf5
    integer :: stat
    
    stat = cm1hdf5%close_cm1()

  end subroutine final_cm1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_cm1_hdf5(self) result(close_cm1)
      implicit none
      class(cm1_hdf5) :: self

      close_cm1 = 0
      ! Check if the dataset is open
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
      self%nx = 0
      self%ny = 0
      self%nz = 0
      self%nt = 0
      self%nv = 0

      self%is_dataset_open = .false.
      if (self%manage_hdf) call deinitiate_hdf
      close_cm1 = 1
      call self%cm1log(LOG_MSG, 'close_cm1', 'Dataset closed')

   end function close_cm1_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_dataset_time_hdf5(self,time) result(open_dataset_time)
      implicit none
      class(cm1_hdf5)     :: self
      integer, intent(in) :: time
      character(len=5) :: dtime
      character(len=256) :: datfile

      open_dataset_time = 0
      if (.not. self%check_dataset_open('open_dataset_time')) return

      if (self%is_time_open) then
         call self%cm1log(LOG_WARN, 'open_dataset_time', 'Internal Error: timelevel inconsistency')
         stop
      end if

      500 format(I5.5)
      write(dtime, 500) time

      datfile = trim(self%path)//'/'//trim(self%basename)//'.'//dtime//'.h5'

      call self%cm1log(LOG_MSG, 'open_dataset_time', 'Opening for Multiread: '//trim(datfile))
      call h5fopen_f(datfile, H5F_ACC_RDONLY_F, self%file_id, self%h5err)

      self%time_open = time
      self%is_time_open = .true.
      open_dataset_time = 1

    end function open_dataset_time_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function close_dataset_time_hdf5(self) result(close_dataset_time)
      implicit none
      class(cm1_hdf5) :: self

      close_dataset_time = 0
      if ((.not. self%check_dataset_open('close_dataset_time')) .or. &
           (.not. self%check_time_open('close_dataset_time'))) return

      call h5fclose_f(self%file_id, self%h5err)
      call self%cm1log(LOG_MSG, 'close_dataset_time', 'Timelevel closed for reading')

      self%time_open = -1
      self%is_time_open = .false.
      close_dataset_time = 1

    end function close_dataset_time_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!TODO: rename this to get_var_dims ?
    integer function get_var_byname_hdf5(self,varname) result(get_var_byname)
       implicit none
       class(cm1_hdf5)                   :: self
       character(len=*), intent(in) :: varname
       integer :: i

       get_var_byname = 0
       if (.not. self%check_dataset_open('get_var_byname')) return
 
       do i = 1, self%nv
          if (varname.eq.self%vars(i)%varname) then
             get_var_byname = self%vars(i)%dims
             return
          endif
       end do
 
       call self%cm1log(LOG_WARN, 'get_var_byname', 'Requested variable does not exist: '//varname)
 
     end function get_var_byname_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2D_at_time_hdf5(self, varname, Field2D) result(read2D_at_time)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), intent(in) :: varname
      real, dimension(self%nx,self%ny) :: Field2D
      integer :: varid

      read2D_at_time = 0
      if ((.not. self%check_dataset_open('read2D_at_time')) .or. &
           (.not. self%check_time_open('read2D_at_time'))) return

      ! Does the variable exist in this dataset?
      varid = self%get_var_byname(varname)
      if (varid.ne.2) then
         call self%cm1log(LOG_WARN, 'read2D_at_time', 'Variable not found: '//trim(varname))
         return
      end if

      call self%cm1log(LOG_INFO, 'read2D_at_time', 'Reading: '//trim(varname))
      call get_h5_2d_float(self%file_id, '/2d/'//trim(varname), self%h5err, Field2D, self%nx, self%ny)

      read2D_at_time = 1

    end function read2D_at_time_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3D_at_time_hdf5(self, varname, Field3D) result(read3D_at_time)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), intent(in) :: varname
      real, dimension(self%nx,self%ny,self%nz) :: Field3D
      integer :: varid

      read3D_at_time = 0
      ! Is the dataset open (has the ctl file been scanned)
      ! Is the dataset open for reading (is the dat file open)
      if ((.not. self%check_dataset_open('read3D_at_time')) .or. &
           (.not. self%check_time_open('read3D_at_time'))) return

      ! Does the variable exist in this dataset?
      varid = self%get_var_byname(varname)
      if (varid.ne.3) then
         call self%cm1log(LOG_WARN, 'read3D_at_time', 'Variable not found: '//trim(varname))
         ! lets see if this variable can be reconstructed before giving up.
         read3D_at_time = self%read3D_derived(varname, Field3D)
         return
      end if

      call self%cm1log(LOG_INFO, 'read3D_at_time', 'Reading: '//trim(varname))
      call get_h5_3d_float(self%file_id, '/3d_'//self%grid//'/'//trim(varname), self%h5err, Field3D, self%nx, self%ny, self%nz)

      read3D_at_time = 1

    end function read3D_at_time_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3D_derived(self, varname, Field3D) result(read3D)
     implicit none
     class(cm1_hdf5) :: self
     character(len=*), intent(in) :: varname
     real, dimension(self%nx, self%ny, self%nz), intent(inout) :: Field3D
     
     integer :: varid0, varid
     character(len=128) :: varname0, varname_pert
     real, dimension(self%nz) :: Field0
     integer :: k

     read3D = 0
     ! Is the dataset open (has the ctl file been scanned)
     ! Is the dataset open for reading (is the dat file open)
     if ((.not. self%check_dataset_open('read3DDerived')) .or. &
          (.not. self%check_time_open('read3DDerived'))) return

     ! This function reads a basestate variable and a perturbation variable to
     ! return the requested variable.  e.g. p = p_0 + p', read p_0 and p' to 
     ! return p

     varname0 = trim(varname)//'0'
     varname_pert = trim(varname)//'pert'

     ! fix some quirky names
     if (varname0 .eq. 'p0') varname0 = 'pres0'
     if (varname0 .eq. 'rho0') varname0 = 'rh0'

     ! Verify variables exist. Returns 0 for non-existance or the number of dimensions if found
     varid0 = self%get_var_byname(varname0)
     varid = self%get_var_byname(varname_pert)
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

   end function read3D_derived

end module ingest_cm1_hdf5
#endif

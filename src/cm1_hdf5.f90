!==========================================================================!
!
!   Ingest_CM1_hdf
!
!   Casey Webster, Dept of Meteorology, Penn State
!
!==========================================================================!



module ingest_cm1_hdf5

   use ingest_cm1_base
   use hdf5

   implicit none

   type, extends(cm1_base) :: cm1_hdf5
      !private
         integer :: h5err
         integer(HID_T) :: file_id
      contains
         !private

         procedure, pass(self) :: scan_hdf
         procedure, public ,pass(self) :: open_cm1 => open_cm1_hdf5
         procedure, public ,pass(self) :: close_cm1 => close_cm1_hdf5
         procedure, public ,pass(self) :: readMultStart => readMultStart_hdf5
         procedure, public ,pass(self) :: readMultStop => readMultStop_hdf5
         procedure, public ,pass(self) :: read3DMult => read3DMult_hdf5
         procedure, public ,pass(self) :: read2DMult => read2DMult_hdf5

         !TODO: finalization needs gfortran 4.9 or ifort
         !!!final :: final_cm1
   end type cm1_hdf5

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function open_cm1_hdf5(self, dsetpath, dsetbasename, dsettype, grid, nodex, nodey, hdfmetadatatime) result(open_cm1)
      implicit none
      class(cm1_hdf5) :: self
      character(len=*), intent(in) :: dsetpath
      character(len=*), intent(in)  :: dsetbasename
      integer, intent(in)            :: dsettype
      character, optional :: grid
      integer, optional :: nodex, nodey, hdfmetadatatime

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

      if (.not.present(hdfmetadatatime)) then
        call self%cm1log(LOG_ERROR, 'open_cm1', 'HDF open requires a valid time to get metadata for the dataset')
        open_cm1 = 0
        return
      endif

      call self%cm1log(LOG_INFO, 'open_cm1', 'Scanning HDF files.')
      open_cm1 = self%scan_hdf(hdfmetadatatime)
      call h5open_f(self%h5err)
      self%isopen = .true.

   end function open_cm1_hdf5


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

  integer function scan_hdf(self, hdfmetatime)
    implicit none
    class(cm1_hdf5)     :: self
    integer, intent(in) :: hdfmetatime
    character(len=256)  :: filename, output
    character(len=5)    :: dtime
    integer(HSIZE_T), dimension(1) :: dims

    !TODO: open the hdf file and read metadata
    ! vars
    ! nt and times
         !integer :: nt, nv
         !type (variable), dimension(:), allocatable  :: vars
         !integer, dimension(:), allocatable :: times

    500 format(I5.5) 
    write(dtime, 500) hdfmetatime

    filename = trim(self%path)//'/'//trim(self%basename)//'.'//dtime//'.h5'

    call self%cm1log(LOG_INFO, 'scan_hdf', 'Scanning metadata in '//trim(filename))
    call h5open_f(self%h5err)
    call h5fopen_f(filename, H5F_ACC_RDONLY_F, self%file_id, self%h5err)

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

!/grid                    Group
!/grid/x0                 Dataset {1}
!/grid/x1                 Dataset {1}
!/grid/y0                 Dataset {1}
!/grid/y1                 Dataset {1}
!/mesh                    Group
!/mesh/dx                 Dataset {1}
!/mesh/dy                 Dataset {1}
!/time                    Dataset {1}

    ! GET VARS

    501 format('... Found ',A1,' dimension with ',i5,' points')
    503 format('... Found ',i3,' variables')

    write (output,501) 'X',self%nx
    call self%cm1log(LOG_MSG, 'scan_hdf', trim(output))
    write (output,501) 'Y',self%ny
    call self%cm1log(LOG_MSG, 'scan_hdf', trim(output))
    write (output,501) 'Z',self%nz
    call self%cm1log(LOG_MSG, 'scan_hdf', trim(output))
    !write (output,501) 'T',self%nt
    !call self%cm1log(LOG_MSG, 'read_ctl', trim(output))
    !write (output,503) self%nv
    !call self%cm1log(LOG_MSG, 'read_ctl', trim(output))

    ! GET TIMES

    scan_hdf = 1
    call h5fclose_f(self%file_id, self%h5err)
    call h5close_f(self%h5err)

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
      call h5close_f(self%h5err)
      close_cm1 = 1
      call self%cm1log(LOG_MSG, 'close_cm1', 'Dataset closed:')

   end function close_cm1_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   integer function getVarByName(self,varname)
!      implicit none
!      class(cm1)                   :: self
!      character(len=*), intent(in) :: varname
!      integer :: i
!
!      if (.not. self%check_open('getVarByName')) then
!         getVarByName = 0
!         return
!      end if
!
!      do i = 1, self%nv
!         if (varname.eq.self%vars(i)%varname) then
!            getVarByName = i
!            return
!         endif
!      end do
!
!      getVarByName = 0
!
!   end function getVarByName

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStart_hdf5(self,time) result(readMultStart)
      implicit none
      class(cm1_hdf5)     :: self
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
      select case (self%dtype)
         case (HDF)
            call self%cm1log(LOG_ERROR, 'read3DMultStart', 'Not implemented.')
      end select

   end function readMultStart_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function readMultStop_hdf5(self) result(readMultStop)
      implicit none
      class(cm1_hdf5) :: self

      if ((.not. self%check_open('readMultStop')) .or. (.not. self%check_mult('readMultStop'))) then
         readMultStop = 0
         return
      end if

      select case (self%dtype)
         case (HDF)
            call self%cm1log(LOG_ERROR, 'read3DMultStop', 'Not implemented.')
      end select

   end function readMultStop_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read2DMult_hdf5(self, varname, Field2D) result(read2DMult)
      implicit none
      class(cm1_hdf5), intent(in) :: self
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
      !status = self%read3DXYSlice(varid, 0, Field2D(:,:))

      read2DMult = 1

   end function read2DMult_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DMult_hdf5(self, varname, Field3D) result(read3DMult)
      implicit none
      class(cm1_hdf5), intent(in) :: self
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
      !do k = 1,self%nz
      !   status = self%read3DXYSlice(varid, k, Field3D(:,:,k))
      !end do

      read3DMult = 1

   end function read3DMult_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module ingest_cm1_hdf5

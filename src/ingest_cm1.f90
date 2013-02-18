! module for ingesting cm1 model output.   
!
! Formats supported : grads style flat files
!                       - single file per time step
!                         - filenames by time index
!                         - filenames by time in seconds
! Future formats    : grads style flat files
!                       - everything in one file
!                     netcdf / hdf5 files


module ingest_cm1

   implicit none

   ! module variables
   character(len=128) :: path
   character(len=64)  :: basename
   integer            :: dtype
   integer            :: nx, ny, nz, nt, nv, dt
   real, dimension(:), allocatable :: x, y, z
   integer, dimension(:), allocatable :: t
   integer            :: n2d, reclen
   integer            :: i,j,k
   integer            :: t0, isopen = 0, ismult = 0
   type variable
      character(len=32) :: varname
      integer           :: levs
      character(len=54) :: units
   end type variable
   type (variable), dimension(:), allocatable :: vars


contains

   integer function open_cm1(dsetpath, dsetbasename, dsettype)
      implicit none
      character(len=256) :: dset
      character(len=128), intent(in) :: dsetpath
      character(len=64), intent(in)  :: dsetbasename
      integer, intent(in)            :: dsettype
      character(len=128) :: tmp

      if (isopen.eq.1) then
         print *,'[ingest_cm1::open_cm1]: Already open, aborting'
         open_cm1 = 0
         return
      end if

      path = dsetpath
      basename = dsetbasename
      dtype = dsettype

      dset = trim(path)//'/'//trim(basename)//'_s.ctl'
      print *, ('[ingest_cm1::open_cm1]: Opening ',trim(dset))

      open(unit=41, file=dset, status='old')
      read(41,*)
      read(41,*)
      read(41,*)
      read(41,*)

      500 format(' [ingest_cm1::open_cm1]: ... Found ',A1,' dimension with ',i5,' points')
      read(41,*) tmp, nx  ! for stretched grids only.  detect and calc linear.
      allocate(x(nx))
      do i = 1, nx
         read(41,*) x(i)
      end do
      print 500,'X',nx

      read(41,*) tmp, ny  ! for stretched grids only.  detect and calc linear.
      allocate(y(ny))
      do j = 1, ny
         read(41,*) y(j)
      end do
      print 500,'Y',ny

      reclen = nx*ny*4

      read(41,*) tmp, nz  ! for stretched grids only.  detect and calc linear.
      allocate(z(nz))
      do k = 1, nz
         read(41,*) z(k)
      end do
      print 500,'Z',nz

      ! calculate timelevels.  note - supports a specific timeformat
      ! and only an integer dt (minimum timestep 1 s). 
      ! tdef     NT linear hh:mmZddMMMYYYY SSSSSYR
      ! where NT is number of timelevels
      !       YYYY is the starting time in seconds
      !       SSSSS is the timestep in seconds
      !
      ! Modify cm1 to write your ctl file this way, or hand edit them 
      read(41,501) nt, t0, dt
      501 format('tdef ',I10,' linear ',11x,I4,1x,I5,'YR')
      allocate(t(nt))
      do i = 1, nt
         t(i) = t0 + ((i-1)*dt)
      end do
      print 500,'T',nt

      ! variables
      read(41,*) tmp, nv
      allocate(vars(nv))
      n2d = 0
      do i = 1,nv
         read(41,502) vars(i)%varname, vars(i)%levs, vars(i)%units
         502 format(A10,1x,I5,6x,A)
         if (vars(i)%levs.eq.0) then
            n2d = n2d+1
         end if
      end do
      503 format(' [ingest_cm1::open_cm1]: ... Found ',i3,' variables')
      print 503,nv

      print *, ('[ingest_cm1::open_cm1]: Closing ',trim(dset))
      close(41)
      isopen = 1
      open_cm1 = 1
   end function open_cm1

   integer function close_cm1()
      implicit none

      if (isopen.eq.0) then
         print *,'[ingest_cm1::close_cm1]: No dataset open, aborting'
         close_cm1 = 0
         return
      end if
      nx = 0
      ny = 0
      nz = 0
      nt = 0
      nv = 0
      deallocate(x)
      deallocate(y)
      deallocate(z)
      deallocate(t)
      deallocate(vars)

      isopen = 0
      close_cm1 = 1
      print *,'[ingest_cm1::close_cm1]: Dataset closed'

   end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   integer function getVarByName(varname)
      implicit none
      character(len=32), intent(in) :: varname
      integer :: i

      if (isopen.eq.0) then
         print *,'[ingest_cm1::getVarByName]: No dataset open, aborting'
         getVarByName = 0
         return
      end if

      do i = 1, nv 
         if (varname.eq.vars(i)%varname) then
            getVarByName = i
            return
         endif
      end do

      getVarByName = 0

   end function getVarByName

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DXYSlice(varid, level, slice)
      implicit none
      integer, intent(in)                 :: varid, level
      real, dimension(nx,ny), intent(out) :: slice
      integer :: idx

      if (isopen.eq.0) then
         print *,'[ingest_cm1::read3DXYSlice]: No dataset open, aborting'
         read3DXYSlice = 0
         return
      end if

      ! Read Slice
      ! NOTE this offset is for one timestep per file!
      idx = (n2d+(varid-1-n2d)*(nz)+1)

      read(42, rec=i) slice

      read3DXYSlice = 1

   end function read3DXYSlice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DMultStart(time)
      implicit none
      integer, intent(in) :: time
      character(len=10) :: dtime
      character(len=256) :: datfile

      if (isopen.eq.0) then
         print *,'[ingest_cm1::read3DMultStart]: No dataset open, aborting'
         read3DMultStart = 0
         return
      end if

      if (ismult.eq.1) then
         print *,'[ingest_cm1::read3DMultStart]: Multiread already started, aborting'
         read3DMultStart = 0
         return
      end if

      ! filename?
      write(dtime,505) time
      505 format(I6.6)
      datfile = trim(path)//'/'//trim(basename)//'_'//trim(dtime)//'_s.dat'
      print *, ('[ingest_cm1::read3DMultStart]: Opening ',trim(datfile))
      
      ! open dat file
      open(42,file=datfile,form='unformatted',access='direct',recl=reclen,status='old')

      print *,'[ingest_cm1::read3DMultStart]: Multiread started for time ',time
      read3DMultStart = 1
      ismult = 1

   end function read3DMultStart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DMultStop()
      implicit none

      if (isopen.eq.0) then
         print *,'[ingest_cm1::read3DMultStop]: No dataset open, aborting'
         read3DMultStop = 0
         return
      end if

      if (ismult.eq.0) then
         print *,'[ingest_cm1::read3DMultStop]: No multiread started, aborting'
         read3DMultStop = 0
         return
      end if

      close(42)
      print *,'[ingest_cm1::read3DMultStop]: Multiread stopped'
      read3DMultStop = 1
      ismult = 0

   end function read3DMultStop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3DMult(varname, Field3D)
      implicit none
      character(len=32), intent(in) :: varname
      real, dimension(nx,ny,nz), intent(out) :: Field3D
      integer :: k,varid,status

      ! Is the dataset open (has the ctl file been scanned)
      if (isopen.eq.0) then
         print *,'[ingest_cm1::read3DMult]: No dataset open, aborting'
         read3DMult = 0
         return
      end if

      ! Is the dataset open for reading (is the dat file open)
      if (isopen.eq.0) then
         print *,'[ingest_cm1::read3DMult]: No multread started, aborting'
         read3DMult = 0
         return
      end if

      ! Does the variable exist in this dataset?
      varid = getVarByName(varname)
      if (varid.eq.0) then
         print *,'[ingest_cm1::read3DMult]: variable not found: ',trim(varname)
         read3DMult = 0
         return
      end if

      print *,'[ingest_cm1::read3DMult]: Reading: ',trim(varname)
      ! Read the variable from the dataset
      do k = 1,nz
         status = read3DXYSlice(varid, k, Field3D(:,:,k))
      end do

      read3DMult = 1

   end function read3DMult

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function read3D(varname, time, Field3D)
      implicit none
      character(len=32), intent(in) :: varname
      integer, intent(in) :: time
      real, dimension(nx,ny,nz), intent(out) :: Field3D
      integer :: status

      if (isopen.eq.0) then
         print *,'[ingest_cm1::read3D]: No dataset open, aborting'
         read3D = 0
         return
      end if

      status = read3DMultStart(time)
      if (status.eq.0) then
         print *,'[ingest_cm1:read3D]: data open failure, aborting'
         read3D = 0
         return
      endif

      status = read3DMult(varname, Field3D)
      if (status.eq.0) then
         print *,'[ingest_cm1:read3D]: data read failure, aborting'
         read3D = 0
         status = read3DMultStop()
         return
      endif

      status = read3DMultStop()
      if (status.eq.0) then
         print *,'[ingest_cm1:read3D]: data close failure, aborting'
         read3D = 0
         return
      endif

      read3D = 1
   end function read3D

end module ingest_cm1

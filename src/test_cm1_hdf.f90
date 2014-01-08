! Test program for the ingest cm1 module

program test_cm1

   use ingest_cm1
!   use testing

   implicit none
   integer :: status

   character(len=128) :: dsetpath
   character(len=64)  :: dsetbasename
   integer            :: dsettype
   character(len=32)  :: varname
   real, dimension(:,:,:), allocatable :: dbz, dum3
   real, dimension(:), allocatable :: x_s, y_s, z_s
   real, dimension(:), allocatable :: x_u, y_v, z_w
   integer :: nx_s, nz_s, ny_s
   integer :: nx_u, nz_w, ny_v
   integer :: hdf_meta_time

   type(cm1_dataset)         :: cm1

   dsetpath = '/home/casey/Research/circuitAnalysis/data/reference/'
   dsetbasename = 'curved90-qv14'
   dsettype = HDF
   hdf_meta_time = 0

   !status = test("Open dataset",1, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype, grids=['s','u','v','w'],hdfmetadatatime=hdf_meta_time)
   call check(status, 1)

   !status = test("Open dataset already open",0, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype, grids=['s','u','v','w'], hdfmetadatatime=hdf_meta_time)
   call check(status, 0)

! TEST some stuff

  print *,'Getting dimensions of data'
  nx_s = cm1%get_nx('s')
  call check(nx_s, 448)
  allocate(x_s(nx_s))
  status = cm1%get_x('s',x_s)

  ny_s = cm1%get_ny('s')
  call check(ny_s, 448)
  allocate(y_s(ny_s))
  status = cm1%get_y('s',y_s)

  nz_s = cm1%get_nz('s')
  call check(nz_s, 68)
  allocate(z_s(nz_s))
  status = cm1%get_z('s',z_s)

  nx_u = cm1%get_nx('u')
  call check(nx_u, 449)
  allocate(x_u(nx_u))
  status = cm1%get_x('u',x_u)

  ny_v = cm1%get_ny('v')
  call check(ny_v, 449)
  allocate(y_v(ny_v))
  status = cm1%get_y('v',y_v)

  nz_w = cm1%get_nz('w')
  call check(nz_w, 69)
  allocate(z_w(nz_w))
  status = cm1%get_z('w',z_w)

  varname = 'dbz'
  print *, 'Allocating dbz dims (s)=',nx_s,ny_s,nz_s
  allocate (dbz(nx_s,ny_s,nz_s))
  print *, 'Fetching dbz'
  status = cm1%read_3d(2400, 's', varname, dbz(:,:,:))
  print *,'status = ',status

  deallocate(dbz)
  print *,'Reallocating and changing time'
  allocate (dbz(nx_s,ny_s,nz_s))
  allocate (dum3(nx_u,ny_s,nz_s))
  print *, 'Fetching dbz'
  status = cm1%read_3d(2100, 's', 'dbz', dbz(:,:,:))
  print *,'status = ',status
  print *, 'Fetching u'
  status = cm1%read_3d(1800, 'u', 'u', dum3(:,:,:))
  print *,'status = ',status
  print *, 'Fetching u'
  status = cm1%read_3d(2100, 'u', 'u', dum3(:,:,:))
  print *,'status = ',status

! Close dataset

   status = cm1%close_dataset()

   if (status.eq.0) then
      print *,'Error'
   else
      print *,'Success'
   endif

   status = cm1%close_dataset()

   if (status.eq.0) then
      print *,'Error'
   else
      print *,'Success'
   endif

   contains

   subroutine check(val, expected)
      implicit none
      integer :: val, expected

      print *,'--------------------------------------------------------'
      if (val == expected) then
         print *,'SUCCESS: ', val
      else
         print *,'FAILURE: ', val
      end if
      print *,'--------------------------------------------------------'
   end subroutine check

end program test_cm1

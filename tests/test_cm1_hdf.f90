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
   integer :: k

   type(cm1_dataset)         :: cm1

   dsetpath = 'data/'
   dsetbasename = 'curved90-qv14'
   dsettype = HDF

   !status = test("Open dataset",1, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype, grids=['s','u','v','w'])
   call check(status, 1)

   !status = test("Open dataset already open",0, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype, grids=['s','u','v','w'])
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
  status = cm1%read_3d(3600, 's', varname, dbz(:,:,:))
  call check(status, 1)

  deallocate(dbz)
  allocate(dbz(nx_s,ny_s,nz_s))
  allocate(dum3(nx_s,ny_s,nz_s))
  print *,'Fettching p'
  status = cm1%read_3d(3600, 's', 'p', dbz(:,:,:))
  call check(status, 1)
  print *,'Fetching ppert'
  status = cm1%read_3d(3600, 's', 'ppert', dum3(:,:,:))
  call check(status, 1)

  ! calculate p0
!TODO need a file with p in it to test
  print *,'Vertical profile of p0 (p, ppert, p0)'
  do k = 1, nz_s
     print *, k, dbz(100,100,k), dum3(100,100,k), dbz(100,100,k)-dum3(100,100,k)
  end do

  ! test slicing
  deallocate(dum3)
  allocate(dum3(1:50,1:50,1:4))
  print *, 'Fetching slice of p'
  status = cm1%read_3d_slice(3600,'s','p', dum3, 1,50,1,50,1,4)
  call check(status,1)
  print *, 'Checking slice of p'
  call checkl(all(dum3-dbz(1:50,1:50,1:4)<0.0001), .true.)
  do k = 1,4
    print *, dum3(1,1,k), dbz(1,1,k), dum3(1,1,k)-dbz(1,1,k)
    print *, dum3(1,50,k), dbz(1,50,k), dum3(1,50,k)-dbz(1,50,k)
    print *, dum3(50,1,k), dbz(50,1,k), dum3(50,1,k)-dbz(50,1,k)
    print *, dum3(50,50,k), dbz(50,50,k), dum3(50,50,k)-dbz(50,50,k)
  end do

! Close dataset

  status = cm1%close_dataset()
  call check(status, 1)

  status = cm1%close_dataset()
  call check(status, 0)

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

   subroutine checkl(val, expected)
      implicit none
      logical :: val, expected

      print *,'--------------------------------------------------------'
      if (val .eqv. expected) then
         print *,'SUCCESS: ', val
      else
         print *,'FAILURE: ', val
      end if
      print *,'--------------------------------------------------------'
   end subroutine checkl

end program test_cm1

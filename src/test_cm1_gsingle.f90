! Test program for the ingest cm1 module

program test_cm1

   use ingest_cm1
!   use testing

   implicit none
   integer :: status

   character(len=128) :: dsetpath
   character(len=64)  :: dsetbasename
   integer            :: dsettype
   integer            :: varid
   character(len=32)  :: varname
   real, dimension(:,:,:), allocatable :: dbz, dum3
   real, dimension(:), allocatable :: myx, myy, myz
   integer :: mynx, mynz, myny

   type(cm1_dataset)         :: cm1

   dsetpath = '/home/casey/Research/test'
   dsetbasename = 'curved90-qv14-2'
   dsettype = GRADSSINGLE

   !status = test("Open dataset",1, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype,grids=['s'])
   print *,'status = ',status
   call check(status, 1)

   !status = test("Open dataset already open",0, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype,grids=['s'])
   print *,'status = ',status
   call check(status, 0)

! TEST some stuff

  print *,'Getting dimensions of data'
  mynx = cm1%get_nx('s')
  call check(mynx, 448)
  allocate(myx(mynx))
  status = cm1%get_x('s',myx)

  myny = cm1%get_ny('s')
  call check(myny, 448)
  allocate(myy(myny))
  status = cm1%get_y('s',myy)

  mynz = cm1%get_nz('s')
  call check(mynz, 68)
  allocate(myz(mynz))
  status = cm1%get_z('s',myz)

  varname = 'dbz'
  print *, 'Allocating dbz dims=',mynx,myny,mynz
  allocate (dbz(mynx,myny,mynz))
  print *, 'Fetching dbz'
  status = cm1%read_3d(4500, 's', varname, dbz(:,:,:))
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

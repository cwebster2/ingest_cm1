! Test program for the ingest cm1 module

program test_cm1

   use ingest_cm1

   implicit none
   integer :: status

   character(len=128) :: dsetpath
   character(len=64)  :: dsetbasename
   integer            :: dsettype
   integer            :: varid
   character(len=32)  :: varname
   real, dimension(:,:,:), allocatable :: dbz, dum3
   real, dimension(:,:,:,:), allocatable :: dbz2
   real, dimension(:), allocatable :: myx, myy, myz, myt
   integer :: mydt, mynx, mynz, myny

   dsetpath = '/home/casey/Research/test'
   dsetbasename = 'curved90-qv14-2'
   dsettype = GRADS

   status = open_cm1(dsetpath,dsetbasename,dsettype)

   if (status.eq.0) then
      print *,'Error'
   else
      print *,'Success'
   endif

   status = open_cm1(dsetpath,dsetbasename,dsettype)

   if (status.eq.0) then
      print *,'Error'
   else
      print *,'Success'
   endif

! TEST some stuff

  varid = getVarByName('uinterp') 
  print *,'Variable uinterp = ',varid

  varname = 'qg'
  varid = getVarByName(varname) 
  print *,'Variable ',varname,' = ',varid

  varname = 'xyzvort'
  varid = getVarByName(varname) 
  print *,'Variable ',varname,' = ',varid

  print *,'Getting dimensions of data'
  mynx = cm1_nx()
  allocate(myx(mynx))
  status = cm1_x(myx)

  myny = cm1_ny()
  allocate(myy(myny))
  status = cm1_y(myy)

  mynz = cm1_nz() + 1
  allocate(myz(mynz))
  status = cm1_z(myz(2:mynz))
  myz(1) = 0.0

  !print *,myx
  !print *,' '
  !print *,myy
  !print *,' '
  !print *,myz
  !print *,' '

  varname = 'dbz'
  print *, 'Allocating dbz dims=',mynx,myny,mynz
  allocate (dbz(mynx,myny,mynz))
  print *, 'Fetching dbz'
  status = read3D(varname, 4500, dbz(:,:,2:mynz))
  print *,'status = ',status

!  print *, 'Allocating dbz dims=',mynx,myny,mynz,1
!  allocate (dbz2(nx,ny,nz,10))
!  print *, 'Fetching dbz2 w/ 4D array'
!  status = read3D(varname, 4500, dbz2(:,:,2:mynz,1))
!  print *,'status = ',status
  
  !deallocate(dbz,dbz2)
  deallocate(dbz)
  print *,'Reallocating and doing multiread'
  allocate (dbz(mynx,myny,mynz))
  allocate (dbz2(mynx,myny,mynz,10))
  allocate (dum3(mynx,myny,mynz))
  status = readMultStart(4500)
  print *,'status = ',status
  status = read3DMult('dbz', dbz(:,:,2:mynz))
  print *,'status = ',status
  status = read3DMult('dbz', dum3)
  dbz2(:,:,2:mynz,1) = dum3(:,:,:)
  print *,'status = ',status
  print *,'Stopping multiread'
  status = readMultStop()
  print *,'status = ',status

! Close dataset

   status = close_cm1()

   if (status.eq.0) then
      print *,'Error'
   else
      print *,'Success'
   endif

   status = close_cm1()

   if (status.eq.0) then
      print *,'Error'
   else
      print *,'Success'
   endif

end program test_cm1

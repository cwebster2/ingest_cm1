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
   real, dimension(:,:,:,:), allocatable :: dbz2
   real, dimension(:), allocatable :: myx, myy, myz, myt
   integer :: mydt, mynx, mynz, myny
   integer :: hdf_meta_time

   type(cm1_dataset)         :: cm1
   type(cm1_hdf5)         :: cm1s

   dsetpath = '/Users/casey/Research/test/'
   dsetbasename = 'curved90-qv14'
   dsettype = HDF
   hdf_meta_time = 0

   !status = test("Open dataset",1, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype, grids=['s','u','v','w'],hdfmetadatatime=hdf_meta_time)
   print *,'status = ',status
   call check(status, 1)

   !status = test("Open dataset already open",0, cm1s%open_cm1, cm1s, dsetpath, dsetbasename, dsettype)
   status = cm1%open_dataset(dsetpath,dsetbasename,dsettype, grids=['s','u','v','w'], hdfmetadatatime=hdf_meta_time)
   print *,'status = ',status
   call check(status, 0)

! TEST some stuff

  varname = 'p'
  !status = test('Get defined variable:'//trim(varname), 32, cm1s%getVarByName, varname)
  varid = cm1s%getVarByName('p')
   call check(varid, 3)

  varname = 'qg'
  !status =  test('Get defined variable:'//trim(varname), 26, cm1s%getVarByName, varname)
  varid = cm1s%getVarByName(varname)
   call check(varid, 3)

  varname = 'xyzvort'
  !status =  test('Get UNdefined variable:'//trim(varname), 0, cm1s%getVarByName, varname)
  varid = cm1s%getVarByName(varname)
   call check(varid, 0)


  print *,'Getting dimensions of data'
  !mynx = test('NX', 300, cm1s%cm1_nx)
  mynx = cm1s%cm1_nx()
  allocate(myx(mynx))
  status = cm1s%cm1_x(myx)

  !myny = test('NY', 300, cm1s%cm1_ny)
  myny = cm1s%cm1_ny()
  allocate(myy(myny))
  status = cm1s%cm1_y(myy)

  !mynz =  test('NZ', 72, cm1s%cm1_nz) + 1
  mynz = cm1s%cm1_nz() + 1
  allocate(myz(mynz))
  status = cm1s%cm1_z(myz(2:mynz))
  myz(1) = 0.0

  varname = 'dbz'
  print *, 'Allocating dbz dims=',mynx,myny,mynz
  allocate (dbz(mynx,myny,mynz))
  print *, 'Fetching dbz'
  status = cm1s%read3D(varname, 2400, dbz(:,:,2:mynz))
  print *,'status = ',status

  deallocate(dbz)
  print *,'Reallocating and doing multiread'
  allocate (dbz(mynx,myny,mynz))
  allocate (dbz2(mynx,myny,mynz,10))
  allocate (dum3(mynx,myny,mynz))
  status = cm1s%readMultStart(2400)
  print *,'status = ',status
  status = cm1s%read3DMult('dbz', dbz(:,:,2:mynz))
  print *,'status = ',status
  status = cm1s%read3DMult('dbz', dum3)
  dbz2(:,:,2:mynz,1) = dum3(:,:,:)
  print *,'status = ',status
  print *,'Stopping multiread'
  status = cm1s%readMultStop()
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

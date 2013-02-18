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
   real, dimension(:,:,:), allocatable :: dbz

   dsetpath = '/home/casey/Research/test'
   dsetbasename = 'curved90-qv14-2'
   dsettype = 3

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

  varname = 'uinterp'
  varid = getVarByName(varname) 
  print *,'Variable ',varname,' = ',varid

  varname = 'qg'
  varid = getVarByName(varname) 
  print *,'Variable ',varname,' = ',varid

  varname = 'xyzvort'
  varid = getVarByName(varname) 
  print *,'Variable ',varname,' = ',varid

  varname = 'dbz'
  allocate (dbz(nx,ny,nz+1))
  status = read3D(varname, 4500, dbz(:,:,2:nz+1))

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

!==========================================================================!
!
!   Ingest_CM1
!
!   Casey Webster, Dept of Meteorology, Penn State
!
!   DESCRIPTION:
!    A Fortran 2003 module to access model output by the CM1 cloud model
!    by George Bryan.
!
!    Formats supported :
!       * grads style flat files
!         * multi file per time step (MPI)
!         * single file per time step
!
!                        Requirement: T dimension values must conincide
!                                     with filenames.
!
!    Future formats:
!      * grads style flat files
!        * all timesteps in one file
!      * netcdf / hdf5 files
!
!    TODO: for HDF support:
!     scan_hdf
!     derived types (if no qvpert, use qv-q00 from basestate)
!==========================================================================!

module ingest_cm1

   use ingest_cm1_base
   use ingest_cm1_grads
   use ingest_cm1_grads_mpi
   use ingest_cm1_hdf5

   implicit none

   ! Object to abstract access to multiple datasets (same type, different grids)
   !   type cm1_dataset

   !   end type cm1_dataset

end module ingest_cm1

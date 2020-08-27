! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_vel_self_attraction_loading
!
!> \brief MPAS ocean module
!> \author Kristin Barton 
!> \date   August 2020
!> \details
!>  This module contains routines for the tidal potential forcing.
!>  Design document located in :
!>    MPAS-Model/docs/ocean/design_docs
!
!-----------------------------------------------------------------------

module ocn_vel_self_attraction_loading

   use mpas_kind_types
   use mpas_constants
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_timer
   use mpas_threading, only : mpas_threading_get_thread_num
   use ocn_constants
   use netcdf
#ifdef _MPI
   use mpi
#endif
   
   
   implicit none
   private
   save

   include 'shtns.f' 

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_compute_self_attraction_loading, &
             ocn_vel_self_attraction_loading_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: tidalPotentialOn
   logical :: selfAttractionLoadingOn
   type :: char_array
     character(:), allocatable :: constituent
   end type
   type(char_array), dimension(37) :: constituentList
   public :: char_array

!***********************************************************************

contains


!***********************************************************************
!
!  routine ocn_compute_self_attraction_loading
!
!> \brief   Computes self-attraction and loading 
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine computes the self-attraction and loading
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_self_attraction_loading(meshPool, forcingPool, diagnosticsPool, dminfo, ssh, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool 
      type (dm_info), intent(in) :: dminfo
      real (kind=RKIND), dimension(:), pointer :: ssh

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: forcingPool 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      
      integer :: process_Size, process_Rank
      integer, dimension(:), pointer :: indexToCellID
      type (field1DInteger), pointer :: toColValuesPtr, fromColValuesPtr
      type (field1DInteger), pointer :: toRowValuesPtr, fromRowValuesPtr
      type (field1DReal), pointer :: toSValuesPtr, fromSValuesPtr
      integer, pointer :: nMpas, nGrid 
      integer, dimension(:), allocatable :: scount, stride ! For MPI_Gather/Scatter
      real, dimension(:), allocatable :: sshGlobal
      integer :: sizeOnProcess

      err = 0
      
!      if (.NOT. selfAttractionLoadingOn) return

      ! MPI commands
      sizeOnProcess = size(ssh)
      call MPI_COMM_RANK( dminfo % comm, process_Rank, err)
      call MPI_COMM_SIZE( dminfo % comm, process_Size, err)
      
      ! Allocate arrays on root and gather sizes
      allocate( scount( process_Size ) )
!      allocate( stride( process_Size ) )
!      stride(process_Rank) = sizeOnProcess
       
      call MPI_GATHER( sizeOnProcess, 1, MPI_INTEGER, scount, 1, MPI_INTEGER, &
                       0, dminfo % comm, err)
      call MPI_BCAST( scount, process_Size, MPI_INTEGER, 0, dminfo % comm, err)

      if (process_Rank .EQ. 0) then
        allocate( sshGlobal( sum(scount) ) )
      end if

      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)

      call MPI_BARRIER(dminfo % comm, err)
      call MPI_GATHERV( ssh, sizeOnProcess, MPI_REAL, sshGlobal, scount, scount, &
                MPI_REAL, 0, dminfo % comm, err)
      
      print *, 'process ', process_Rank, ' size: ', sizeOnProcess

      call MPI_BARRIER(dminfo % comm, err) 
      if (process_Rank .EQ. 0) then

        ! Get weight data from the pool
        call mpas_pool_get_field(forcingPool, 'toColValues', toColValuesPtr)
        call mpas_pool_get_field(forcingPool, 'fromColValues', fromColValuesPtr)
        call mpas_pool_get_field(forcingPool, 'toRowValues', toRowValuesPtr)
        call mpas_pool_get_field(forcingPool, 'fromRowValues', fromRowValuesPtr)
        call mpas_pool_get_dimension(forcingPool, 'nMpas', nMpas)
        call mpas_pool_get_dimension(forcingPool, 'nGrid', nGrid)

        print *, 'scount: ', scount
        print *, '---------------------------------'

      end if

      call MPI_BARRIER(dminfo % comm, err)
      call MPI_SCATTERV( sshGlobal, scount, scount, MPI_REAL, ssh, sizeOnProcess, &
                MPI_REAL, 0, dminfo % comm, err)

      deallocate( scount )
      if (process_Rank .EQ. 0) then
        deallocate( sshGlobal )
      end if

   end subroutine ocn_compute_self_attraction_loading!}}}


!***********************************************************************
!
!  routine ocn_vel_self_attraction_loading_init
!
!> \brief   Initializes ocean tidal protential forcing module.
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine initializes the ocean self-attraction and loading module 
!>  and variables.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_self_attraction_loading_init(domain,err)!{{{

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_use_self_attraction_loading
      type (block_type), pointer :: block_ptr
      type (mpas_pool_type), pointer :: forcingPool

      integer :: toNcId, toStatus, toNsDimId, toRowId, toColId, toSId
      integer :: fromNcId, fromStatus, fromNsDimId, fromRowId, fromColId, fromSId
      integer:: nMpas, nGrid, nMpasDimId, nGridDimId
      integer:: toNsLen, fromNsLen 
      integer :: toRowDims, toRowAtts, toColdDims, toColAtts, toNsDims, toNsAtts
      integer :: fromRowDims, fromRowAtts, fromColdDims, fromColAtts, fromNsDims, fromNsAtts
      integer, dimension(:), target, allocatable :: toRowValues, toColValues
      integer, dimension(:), target, allocatable :: fromRowValues, fromColValues
      real, dimension(:), target, allocatable :: toSValues, fromSValues
      character (len = NF90_MAX_NAME) :: toNsName, fromNsName
      character (len = NF90_MAX_NAME) :: nMpasName, nGridName
      integer, pointer :: n_s

      type (field1DInteger), pointer :: toRowValuesPtr, toColValuesPtr
      type (field1Dinteger), pointer :: fromRowValuesPtr, fromColValuesPtr
      type (field1DReal), pointer :: toSValuesPtr, fromSValuesPtr

      integer :: size_Of_Cluster

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_self_attraction_loading', config_use_self_attraction_loading)

      selfAttractionLoadingOn = .false.
      if (config_use_self_attraction_loading) then
         selfAttractionLoadingOn = .true.
      else
! Remove comment when a testcase is setup with this true
!         return
      end if

      ! Open netcdf weights files
      call check( nf90_open(path = "mpas_to_grid.nc", mode = nf90_nowrite, ncid = toNcId) )
      call check( nf90_open(path = "grid_to_mpas.nc", mode = nf90_nowrite, ncid = fromNcId) )

      ! Get dimension ID
      call check( nf90_inq_dimid(toNcId, "n_s", toNsDimId) )
      call check( nf90_inq_dimid(fromNcId, "n_s", fromNsDimId) )
      call check( nf90_inq_dimid(toNcId, "n_a", nMpasDimId) )
      call check( nf90_inq_dimid(fromNcId, "n_a", nGridDimId) )

      ! Get Variable IDs
      call check( nf90_inq_varid(toNcId, "row", toRowId) )
      call check( nf90_inq_varid(toNcId, "col", toColId) )
      call check( nf90_inq_varid(toNcId, "S", toSId) )
      call check( nf90_inq_varid(fromNcId, "row", fromRowId) )
      call check( nf90_inq_varid(fromNcId, "col", fromColId) )
      call check( nf90_inq_varid(fromNcId, "S", fromSId) )

      ! Get Dimension Length
      call check( nf90_inquire_dimension(toNcId, toNsDimId, toNsName, toNsLen) )
      call check( nf90_inquire_dimension(fromNcId, fromNsDimId, fromNsName, fromNsLen) )
      call check( nf90_inquire_dimension(toNcId, nMpasDimId, nMpasName, nMpas) )
      call check( nf90_inquire_dimension(fromNcId, nGridDimId, nGridName, nGrid) )

      ! Allocate matrices to read data into
      allocate (toRowValues (toNsLen) )
      allocate ( toColValues (toNsLen) )
      allocate ( toSValues (toNsLen) )
      allocate (fromRowValues (fromNsLen) )
      allocate ( fromColValues (fromNsLen) )
      allocate ( fromSValues (fromNsLen) )

      ! Retrieve data
      call check( nf90_get_var(toNcId, toColId, toColValues(:) ) )
      call check( nf90_get_var(toNcId, toRowId, toRowValues(:) ) )
      call check( nf90_get_var(toNcId, toSId, toSValues(:) ) )
      call check( nf90_get_var(fromNcId, fromColId, fromColValues(:) ) )
      call check( nf90_get_var(fromNcId, fromRowId, fromRowValues(:) ) )
      call check( nf90_get_var(fromNcId, fromSId, fromSValues(:) ) )

      ! Add values into a field_type
      allocate( toColValuesPtr )
      allocate( toRowValuesPtr )
      allocate( toSValuesPtr )
      allocate( fromColValuesPtr )
      allocate( fromRowValuesPtr )
      allocate( fromSValuesPtr )

      allocate( toColValuesPtr % array(toNsLen) )
      allocate( toRowValuesPtr % array(toNsLen) )
      allocate( toSValuesPtr % array(toNsLen) )
      allocate( fromColValuesPtr % array(fromNsLen) )
      allocate( fromRowValuesPtr % array(fromNsLen) )
      allocate( fromSValuesPtr % array(fromNsLen) )

      toColValuesPtr % array = toColValues
      toRowValuesPtr % array = toRowValues
      toSValuesPtr % array = toSValues
      fromColValuesPtr % array = fromColValues
      fromRowValuesPtr % array = fromRowValues
      fromSValuesPtr % array = fromSValues

      ! Add fields into the forcing pool
      block_ptr => domain % blocklist
      do while ( associated( block_ptr ) )
         call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
 
        call mpas_pool_add_dimension(forcingPool, 'nMpas', nMpas)
        call mpas_pool_add_dimension(forcingPool, 'nGrid', nGrid)
        call mpas_pool_add_field(forcingPool, 'toColValues', toColValuesPtr)
        call mpas_pool_add_field(forcingPool, 'fromColValues', fromColValuesPtr)
        call mpas_pool_add_field(forcingPool, 'toRowValues', toRowValuesPtr)
        call mpas_pool_add_field(forcingPool, 'fromRowValues', fromRowValuesPtr)
        call mpas_pool_add_field(forcingPool, 'toSValues', toSValuesPtr)
        call mpas_pool_add_field(forcingPool, 'fromSValues', fromSValuesPtr)
        
        block_ptr => block_ptr % next
      end do

      ! Deallocate space
      deallocate( toRowValues )
      deallocate( toColValues )
      deallocate( toSValues )
      deallocate( fromRowValues )
      deallocate( fromColValues )
      deallocate( fromSValues )

      if (domain % dminfo % my_proc_id .EQ. 0) then
        print *, "Init"
      end if

   end subroutine ocn_vel_self_attraction_loading_init!}}}

!***********************************************************************
!
!  routine check
!
!> \brief   Check status of netcdf operations
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine checks to status of the netcdf file
!
!-----------------------------------------------------------------------

  subroutine check(status) !{{{
     integer, intent ( in) :: status

     if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        stop "Stopped"
     end if
  end subroutine!}}}

!***********************************************************************
!
!  routine interpolate
!
!> \brief   Perform interpolation
!> \author  Kristin Barton
!> \date    August 2020
!> \details
!>  This routine checks contains the sparse matrix multiplication
!>    algorithm to interpolate between MPAS and Gaussian Grid
!
!-----------------------------------------------------------------------

  subroutine interpolate( colValues, rowValues, sValues, dataIn, dataOut) !{{{

      integer, dimension(:), intent(in) :: colValues, rowValues
      real, dimension(:), intent(in) :: sValues, dataIn
      real, dimension(:), intent(out) :: dataOut
      real :: rhs = 0
      integer :: n_S, n, nRow

      n_S = size(sValues)
      n = 1

      do while (n .LE. n_S)
        nRow = rowValues( n )
        do while ( rowValues( n ) .EQ. nRow )
          rhs = rhs + dataIn(n) * sValues(n)
          n = n + 1
        end do
        dataOut(nRow) = rhs
        rhs = 0
      end do

  end subroutine!}}}

end module ocn_vel_self_attraction_loading!}}}

!***********************************************************************

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

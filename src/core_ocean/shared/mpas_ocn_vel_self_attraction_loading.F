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
!> \date   September 2019
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

   public :: ocn_vel_self_attraction_loading_tend, &
             ocn_compute_self_attraction_loading, &
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
!  routine ocn_vel_self_attraction_loading_tend
!
!> \brief   Computes tendency term for tidal potential
!> \author  Kristin Barton
!> \date    April 2020
!> \details
!>  This routine computes the tidal potential tendency for momentum
!>  based on current state.
!
!-----------------------------------------------------------------------

   subroutine ocn_vel_self_attraction_loading_tend(meshPool, forcingPool, ssh, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), intent(in) :: &
         ssh !< Input: Sea surface height

      type (mpas_pool_type), intent(in) :: meshPool          !< Input: mesh information
      type (mpas_pool_type), intent(in) :: forcingPool       !< Input: forcinginformation

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iEdge, k, cell1, cell2, nEdges
      integer, dimension(:), pointer :: nEdgesArray
      integer, dimension(:), pointer :: maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask

      real (kind=RKIND), dimension(:), pointer :: dcEdge
      real (kind=RKIND), dimension(:), pointer :: tidalPotentialEta
      real (kind=RKIND), pointer :: config_self_attraction_and_loading_beta
      real (kind=RKIND) :: invdcEdge
      real (kind=RKIND) :: potentialGrad
      logical, pointer :: config_use_tidal_potential_forcing

      err = 0

   end subroutine ocn_vel_self_attraction_loading_tend!}}}

!***********************************************************************
!
!  routine ocn_compute_self_attraction_loading
!
!> \brief   Computes equilibrium tidal potential 
!> \author  Kristin Barton
!> \date    September 2019
!> \details
!>  This routine computes the equilibrium tidal potential due to all 
!>  tidal constituents 
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_self_attraction_loading(meshPool, forcingPool, diagnosticsPool, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: diagnosticsPool 

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
      
      integer :: threadNum

      err = 0

   end subroutine ocn_compute_self_attraction_loading!}}}


!***********************************************************************
!
!  routine ocn_vel_self_attraction_loading_init
!
!> \brief   Initializes ocean tidal protential forcing module.
!> \author  Kristin Barton
!> \date    September 2019
!> \details
!>  This routine initializes the ocean tidal potential forcing module 
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
      integer:: toNsLen, fromNsLen
      integer :: toRowDims, toRowAtts, toColdDims, toColAtts, toNsDims, toNsAtts
      integer :: fromRowDims, fromRowAtts, fromColdDims, fromColAtts, fromNsDims, fromNsAtts
      integer, dimension(:), target, allocatable :: toRowValues, toColValues
      integer, dimension(:), target, allocatable :: fromRowValues, fromColValues
      real, dimension(:), target, allocatable :: toSValues, fromSValues
      character (len = NF90_MAX_NAME) :: toNsName, fromNsName
      integer, pointer :: n_s

      type (field1DInteger), pointer :: toRowValuesPtr, toColValuesPtr
      type (field1Dinteger), pointer :: fromRowValuesPtr, fromColValuesPtr
      type (field1DReal), pointer :: toSValuesPtr, fromSValuesPtr

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
 
        call mpas_pool_add_dimension(forcingPool, 'toNsLen', toNsLen)
        call mpas_pool_add_dimension(forcingPool, 'fromNsLen', fromNsLen)
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

      print *, "init"

      err = 0


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

end module ocn_vel_self_attraction_loading!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker

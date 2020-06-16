!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Computes the localised moments of subsystems given by Kim-Gordon embedding
!> \author Lukas Schreder
! **************************************************************************************************
module kg_moments

   use cp_control_types,      only: dft_control_type
   use cp_dbcsr_operations,   only: dbcsr_allocate_matrix_set
   use cp_log_handling,       only: cp_logger_type
   use cp_output_handling,    only: cp_print_key_unit_nr
   use dbcsr_api,             only: dbcsr_p_type
   use input_section_types,   only: section_vals_type
   use kinds,                 only: default_string_length, &
                                    dp
   use moments_utils,         only: get_reference_point
   use orbital_pointers,      only: current_maxl
   use qs_environment_types,  only: get_qs_env, &
                                    qs_environment_type
   use qs_rho_types,          only: qs_rho_get, &
                                    qs_rho_type
  
#include "./base/base_uses.f90"

   implicit none

   private

! *** globals ***

   character(len=*), parameter, private :: moduleN = "kg_moments"

! *** public subroutines ***

   public :: calculate_kg_moments

! *** body ***

contains

! **************************************************************************************************
!> \brief caculates and prints the expectation values of electronic moments of KG subsystems
! **************************************************************************************************
   subroutine calculate_kg_moments(qs_env, nmoments, reference, ref_point, input, logger)

      type(qs_environment_type), pointer  :: qs_env
      integer, intent(in)                 :: nmoments, reference
      real(dp), intent(in), pointer       :: ref_point(:)
      type(section_vals_type), pointer    :: input
      type(cp_logger_type), pointer       :: logger

      character(len=*), parameter  :: routineN = "calculate_kg_moments", &
           routineP = moduleN//":"//routineN

      integer                          :: handle, ispin, isubsys, nm, nspins, nsubsys, unit_nr
      real(dp)                         :: rcc(3)
      real(dp), allocatable            :: moments(:,:)
      type(dbcsr_p_type), pointer      :: moments_op(:)
      type(dft_control_type), pointer  :: dft_control
      type(qs_rho_type), pointer       :: rho

      call timeset(routineN, handle)
      
      ! sanity checks
      CPASSERT(nmoments <= current_maxl)

      ! output handling
      unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MOMENTS", extension="KGMom", &
           middle_name="KG_MOMENTS")
      
      call get_qs_env(qs_env, dft_control=dft_control, rho=rho)
      call get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      ! number of all moments (e.g. 3 for dipole)
      nm = (6 + 11 * nmoments + 6 * nmoments**2 + nmoments**3) / 6 - 1
      nspins = dft_control%nspins

      ! slice and dice the block-diagonal density matrix

      ! build the moments operator
      call dbcsr_allocate_matrix_set(moments_op, nm)
      
      ! get expectation values
      do ispin = 1, nspins
         
         moments = 0.0_dp
         do isubsys = 1, nsubsys

         end do
         
         
      end do
         
      ! print results

      ! cleanup

      call timestop(handle)
      
   end subroutine calculate_kg_moments

end module kg_moments

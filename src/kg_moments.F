!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Computes the localised moments of subsystems given by Kim-Gordon embedding
!> \author Lukas Schreder
! **************************************************************************************************
module kg_moments

   use qs_environment_types,  only: qs_environment_type
  
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
   subroutine calculate_kg_moments(qs_env)

      type(qs_environment_type), pointer  :: qs_env

      character(len=*), parameter  :: routineN = "calculate_kg_moments", &
           routineP = moduleN//":"//routineN

      integer :: handle

      call timeset(routineN, handle)

      call get_qs_env(qs_env)

      ! slice and dice the block-diagonal density matrix

      ! get expectation values

      ! print results

      call timestop(handle)
      
   end subroutine calculate_kg_moments

end module kg_moments

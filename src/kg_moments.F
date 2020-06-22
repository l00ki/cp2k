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
   use cp_output_handling,    only: cp_iter_string, &
                                    cp_p_file, &
                                    cp_print_key_finished_output, &
                                    cp_print_key_should_output, &
                                    cp_print_key_unit_nr
   use dbcsr_api,             only: dbcsr_copy, &
                                    dbcsr_dot, &
                                    dbcsr_p_type, &
                                    dbcsr_print, &
                                    dbcsr_release, &
                                    dbcsr_set
   use input_section_types,   only: section_vals_type
   use kg_environment_types,  only: kg_environment_type
   use kinds,                 only: default_string_length, &
                                    dp
   use molecule_types,        only: molecule_type
   use moments_utils,         only: get_reference_point
   use orbital_pointers,      only: current_maxl
   use pw_types,              only: pw_p_type
   use qs_environment_types,  only: get_qs_env, &
                                    qs_environment_type
   use qs_moments,            only: build_local_moment_matrix
   use qs_rho_methods,        only: qs_rho_rebuild, &
                                    qs_rho_update_rho
   use qs_rho_types,          only: qs_rho_create, &
                                    qs_rho_get, &
                                    qs_rho_set, &
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

      character(len=default_string_length)  :: iter
      integer                               :: handle, im, ispin, isubsys, nm, nspins, nsubsys, &
                                               unit_nr
      logical                               :: first_time
      real(dp)                              :: rcc(3), strace, trace
      real(dp), allocatable                 :: moments(:,:)
      type(dbcsr_p_type), pointer           :: density_matrix(:), matrix_s(:), moments_op(:), rho_ao(:)
      type(dft_control_type), pointer       :: dft_control
      type(kg_environment_type), pointer    :: kg_env
      type(molecule_type), pointer          :: molecules(:)
      type(pw_p_type), pointer              :: rho_r(:)
      type(qs_rho_type), pointer            :: rho_old, rho_struct

      call timeset(routineN, handle)
      
      ! sanity checks
      CPASSERT(nmoments <= current_maxl)

      ! output handling
      unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MOMENTS", extension=".KGMom", &
           middle_name="KG_MOMENTS")
      iter = cp_iter_string(logger%iter_info)
      
      call get_qs_env(qs_env, dft_control=dft_control, kg_env=kg_env, matrix_s=matrix_s, molecule_set=molecules, &
           rho=rho_old)
      call get_reference_point(rcc, qs_env=qs_env, reference=reference, ref_point=ref_point)
      ! number of all moments (e.g. 3 for dipole)
      nm = (6 + 11 * nmoments + 6 * nmoments**2 + nmoments**3) / 6 - 1
      nspins = dft_control%nspins
      nsubsys = kg_env%nsubsets

      ! get density matrix
      call qs_rho_get(rho_old, rho_ao=density_matrix)
      ! initialize new density
      call qs_rho_create(rho_struct)
      ! set density matrix to the blocked matrix
      call qs_rho_set(rho_struct, rho_ao=density_matrix)
      call qs_rho_rebuild(rho_struct, qs_env, rebuild_ao=.false., rebuild_grids=.true.)
      ! full density kinetic energy term
      call qs_rho_update_rho(rho_struct, qs_env)
      ! get blocked density on the grid
      call qs_rho_get(rho_struct, rho_r=rho_r)

      if (btest(cp_print_key_should_output(logger%iter_info, input, "DFT%PRINT%MOMENTS", &
           first_time=first_time), cp_p_file)) then

         allocate(moments(nsubsys+1, nm))
         allocate(moments_op(nm))
         ! build the moments operator
         do im = 1, nm
            allocate(moments_op(im)%matrix)
            call dbcsr_copy(moments_op(im)%matrix, matrix_s(1)%matrix, "moments operator")
            call dbcsr_set(moments_op(im)%matrix, 0.0_dp)
         end do
         call build_local_moment_matrix(qs_env, moments_op, nmoments, ref_point=rcc)

         do isubsys = 1, nsubsys
            ! calculates the density of the subsystem
            call qs_rho_update_rho(rho_struct, qs_env, &
                 task_list_external=kg_env%subset(isubsys)%task_list)
            call qs_rho_get(rho_struct, rho_ao=rho_ao)
            
            ! calculate eigenvalues
            do im = 1, nm
               trace = 0.0_dp
               do ispin = 1, nspins
                  call dbcsr_print(rho_ao(ispin)%matrix, unit_nr=666)
                  call dbcsr_dot(rho_ao(ispin)%matrix, moments_op(im)%matrix, strace)
                  trace = trace + strace
               end do
               moments(isubsys, im) = trace
            end do
         end do
         
         do ispin = 1, nspins
            ! print results
            if (unit_nr > 0) then
               if (first_time) then
                  write (unit_nr, "(a)") "# iter_level, ispin, subsystem, moments [a.u.]"
               end if
               do isubsys = 1, nsubsys
                  write (unit_nr, "(a10,i3,1X,a20,39es26.16)") iter(1:10), ispin, &
                       molecules(isubsys)%molecule_kind%name, moments(isubsys, :)
               end do
            end if
         end do

         ! cleanup
         do im = 1, nm
            call dbcsr_release(moments_op(im)%matrix)
            deallocate(moments_op(im)%matrix)
         end do
         deallocate(moments_op)
         deallocate(moments)
      end if

      call cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MOMENTS")
      call timestop(handle)
      
   end subroutine calculate_kg_moments

end module kg_moments

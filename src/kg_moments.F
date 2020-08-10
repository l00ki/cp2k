!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Computes the localised moments of subsystems given by Kim-Gordon embedding
!> \author Lukas Schreder
! **************************************************************************************************
module kg_moments

   use atomic_kind_types,       only: atomic_kind_type, &
                                      get_atomic_kind
   use cell_types,              only: cell_type, &
                                      pbc
   use cp_control_types,        only: dft_control_type
   use cp_dbcsr_operations,     only: copy_dbcsr_to_fm, &
                                      dbcsr_allocate_matrix_set
   use cp_fm_basic_linalg,      only: cp_fm_trace
   use cp_fm_struct,            only: cp_fm_struct_create, &
                                      cp_fm_struct_release, &
                                      cp_fm_struct_type
   use cp_fm_types,             only: cp_fm_create, &
                                      cp_fm_get_info, &
                                      cp_fm_p_type, &
                                      cp_fm_release, &
                                      cp_fm_set_all, &
                                      cp_fm_to_fm, &
                                      cp_fm_to_fm_submat, &
                                      cp_fm_type, &
                                      cp_fm_write_formatted
   use cp_log_handling,         only: cp_logger_type, &
                                      cp_to_string
   use cp_output_handling,      only: cp_iter_string, &
                                      cp_p_file, &
                                      cp_print_key_finished_output, &
                                      cp_print_key_should_output, &
                                      cp_print_key_unit_nr
   use dbcsr_api,               only: dbcsr_copy, &
                                      dbcsr_deallocate_matrix, &
                                      dbcsr_dot, &
                                      dbcsr_get_block_p, &
                                      dbcsr_get_info, &
                                      dbcsr_p_type, &
                                      dbcsr_print, &
                                      dbcsr_release, &
                                      dbcsr_set, &
                                      dbcsr_type
   use input_constants,         only: use_mom_ref_com, &
                                      use_mom_ref_zero
   use input_section_types,     only: section_vals_type
   use kg_environment,          only: kg_build_neighborlist
   use kg_environment_types,    only: kg_environment_type
   use kinds,                   only: default_string_length, &
                                      dp
   use molecule_types,          only: get_molecule, &
                                      get_molecule_set_info, &
                                      molecule_type
   use moments_utils,           only: get_reference_point
   use orbital_pointers,        only: current_maxl, &
                                      indco
   use particle_types,          onlY: particle_type
   use pw_types,                only: pw_p_type
   use qs_environment_types,    only: get_qs_env, &
                                      qs_environment_type
   use qs_kind_types,           only: qs_kind_type, &
                                      get_qs_kind
   use qs_moments,              only: build_local_moment_matrix
   use qs_neighbor_list_types,  only: neighbor_list_set_p_type
   use qs_rho_methods,          only: qs_rho_rebuild, &
                                      qs_rho_update_rho
   use qs_rho_types,            only: qs_rho_create, &
                                      qs_rho_get, &
                                      qs_rho_set, &
                                      qs_rho_type
   use wannier_moments,         only: wrap
  
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
      integer                               :: akind, handle, iatom, ibasis, im, ispin, lx, ly, lz, isubsys, nao, nions, &
                                               nm, nspins, nsubsys, unit_nr
      integer, allocatable                  :: nbasis(:)
      logical                               :: first_time, floating, ghost
      real(dp)                              :: r(3), rcc(3), strace, trace, zeff, zmom
      real(dp), allocatable                 :: moments(:, :), nuclear_moments(:)
      type(atomic_kind_type), pointer       :: atomic_kind
      type(cell_type), pointer              :: cell
      type(cp_fm_p_type), pointer           :: density_fm(:), density_slices(:, :), moments_op_fm(:)
      type(cp_fm_type), pointer             :: zero_density
      type(cp_fm_struct_type), pointer      :: density_fmstruct, mo_coeff_struct
      type(dbcsr_p_type), pointer           :: density_matrix(:), matrix_s(:), moments_op(:)
      type(dft_control_type), pointer       :: dft_control
      type(kg_environment_type), pointer    :: kg_env
      type(molecule_type), pointer          :: molecules(:)
      type(particle_type), pointer          :: particles(:)
      type(qs_kind_type), pointer           :: kinds(:)
      type(qs_rho_type), pointer            :: rho

      call timeset(routineN, handle)
      print *, "reference", reference 
      ! sanity checks
      CPASSERT(nmoments <= current_maxl)

      ! output handling
      unit_nr = cp_print_key_unit_nr(logger, input, "DFT%PRINT%MOMENTS", extension=".KGMom", &
           middle_name="KG_MOMENTS")
      iter = cp_iter_string(logger%iter_info)
      
      call get_qs_env(qs_env, cell=cell, dft_control=dft_control, kg_env=kg_env, &
           matrix_s=matrix_s, molecule_set=molecules, particle_set=particles, qs_kind_set=kinds, &
           rho=rho, natom=nions)
      call get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, ref_point=ref_point)
      !rcc = wrap(rcc, cell)
      call cp_fm_get_info(qs_env%mos(1)%mo_set%mo_coeff, matrix_struct=mo_coeff_struct, &
           nrow_global=nao)
      call cp_fm_struct_create(density_fmstruct, template_fmstruct=mo_coeff_struct, &
           ncol_global=nao)
      call cp_fm_create(zero_density, density_fmstruct)
      call cp_fm_set_all(zero_density, 0.0_dp)
      ! number of all moments (e.g. 3 for dipole)
      nm = (6 + 11 * nmoments + 6 * nmoments**2 + nmoments**3) / 6 - 1
      nspins = dft_control%nspins
      nsubsys = kg_env%nsubsets

      if (btest(cp_print_key_should_output(logger%iter_info, input, "DFT%PRINT%MOMENTS", &
           first_time=first_time), cp_p_file)) then

         ! get the block-diagonal density matrix
         call qs_rho_get(rho, rho_ao=density_matrix)
         allocate(density_fm(nspins))
         do ispin = 1, nspins
            call cp_fm_create(density_fm(ispin)%matrix, density_fmstruct)
            call copy_dbcsr_to_fm(density_matrix(ispin)%matrix, density_fm(ispin)%matrix)
         end do

         ! build the moments operator
         allocate(moments(nsubsys, nm))
         allocate(moments_op(nm))
         allocate(moments_op_fm(nm))
         allocate(nbasis(nsubsys))
         do im = 1, nm
            allocate(moments_op(im)%matrix)
            call dbcsr_copy(moments_op(im)%matrix, matrix_s(1)%matrix, "moments operator")
            call dbcsr_set(moments_op(im)%matrix, 0.0_dp)
         end do
         call build_local_moment_matrix(qs_env, moments_op, nmoments, ref_point=rcc)
         do im = 1, nm
            !allocate(moments_op_fm(im)%matrix)
            call cp_fm_create(moments_op_fm(im)%matrix, density_fmstruct)
            call copy_dbcsr_to_fm(moments_op(im)%matrix, moments_op_fm(im)%matrix)
         end do

         ! slice and dice the density matrix
         allocate(density_slices(nsubsys, nspins))
         call get_molecule_set_info(molecules, mol_to_nbasis=nbasis)
         ibasis = 1
         do isubsys = 1, nsubsys
            ! get the Kim-Gordon blocks
            do ispin = 1, nspins
               !allocate(density_slices(isubsys, ispin)%matrix)
               call cp_fm_create(density_slices(isubsys, ispin)%matrix, density_fmstruct)
               call cp_fm_to_fm(zero_density, density_slices(isubsys, ispin)%matrix)
               call cp_fm_to_fm_submat(density_fm(ispin)%matrix, &
                    density_slices(isubsys, ispin)%matrix, nbasis(isubsys), nbasis(isubsys), &
                    ibasis, ibasis, ibasis, ibasis)
            end do
            ibasis = ibasis + nbasis(isubsys)
         end do

         ! calculate eigenvalues
         do isubsys = 1, nsubsys
            do im = 1, nm
               strace = 0.0_dp
               do ispin = 1, nspins
                  call cp_fm_trace(density_slices(isubsys, ispin)%matrix, &
                       moments_op_fm(im)%matrix, trace)
                  strace = strace - trace
               end do
               moments(isubsys, im) = strace
            end do
         end do
         
         ! total nuclear contribution
         allocate(nuclear_moments(nm))
         nuclear_moments = 0.0_dp
         rcc = pbc(rcc, cell)
         do im = 1, nm
            lx = indco(1, im+1)
            ly = indco(2, im+1)
            lz = indco(3, im+1)
            do iatom = 1, size(particles)
               atomic_kind => particles(iatom)%atomic_kind
               call get_atomic_kind(atomic_kind, kind_number=akind)
               call get_qs_kind(kinds(akind), ghost=ghost, floating=floating)
               if (.not. ghost .and. .not. floating) then
                  call get_qs_kind(kinds(akind), core_charge=zeff)
                  !r = wrap(particles(iatom)%r - rcc, cell)
                  r = pbc(particles(iatom)%r - rcc, cell)
                  zmom = zeff
                  if (lx /= 0) zmom = zmom * r(1)**lx
                  if (ly /= 0) zmom = zmom * r(2)**ly
                  if (lz /= 0) zmom = zmom * r(3)**lz
                  nuclear_moments(im) = nuclear_moments(im) + zmom
               end if
            end do  ! atoms
         end do  ! moments

         do ispin = 1, nspins
            ! print results
            if (unit_nr > 0) then
               if (first_time) then
                  write (unit_nr, "(a)") "# iter_level, ispin, subsystem, moments [a.u.]"
               end if
               write (unit_nr, "(a10,i3,1X,a20,39es26.16)") iter(1:10), ispin, "nuclear", nuclear_moments(:)
               do isubsys = 1, nsubsys
                  write (unit_nr, "(a10,i3,1X,a20,39es26.16)") iter(1:10), ispin, &
                       molecules(isubsys)%molecule_kind%name, moments(isubsys, :)
               end do
            end if
         end do

         ! cleanup
         do isubsys = 1, nsubsys
            do ispin = 1, nspins
               call cp_fm_release(density_slices(isubsys, ispin)%matrix)
            end do
         end do
         deallocate(density_slices)
         
         do ispin = 1, nspins
            call cp_fm_release(density_fm(ispin)%matrix)
         end do
         deallocate(density_fm)
         do im = 1, nm
            call dbcsr_deallocate_matrix(moments_op(im)%matrix)
            call cp_fm_release(moments_op_fm(im)%matrix)
         end do
         deallocate(moments_op)
         deallocate(moments_op_fm)
         deallocate(moments)
      end if

      ! cleanup
      call cp_fm_struct_release(density_fmstruct)
      call cp_fm_release(zero_density)

      call cp_print_key_finished_output(unit_nr, logger, input, "DFT%PRINT%MOMENTS")
      call timestop(handle)
      
   end subroutine calculate_kg_moments

end module kg_moments

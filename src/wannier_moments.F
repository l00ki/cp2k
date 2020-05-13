!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the calculation of moments from Wannier functions
!> \author Lukas Schreder
! **************************************************************************************************
module wannier_moments

   use cell_types,            only: cell_type, &
                                    get_cell
   use cp_cfm_basic_linalg,   only: cp_cfm_gemm
   use cp_cfm_types,          only: cp_cfm_create, &
                                    cp_cfm_p_type, &
                                    cp_cfm_release, &
                                    cp_cfm_set_all, &
                                    cp_cfm_to_cfm, &
                                    cp_cfm_type, &
                                    cp_cfm_write_formatted, &
                                    cp_fm_to_cfm
   use cp_control_types,      only: dft_control_type
   use cp_dbcsr_operations,   only: copy_dbcsr_to_fm
   use cp_fm_struct,          only: cp_fm_struct_create, &
                                    cp_fm_struct_release, &
                                    cp_fm_struct_type
   use cp_fm_types,           only: cp_fm_create, &
                                    cp_fm_p_type, &
                                    cp_fm_release, &
                                    cp_fm_set_all, &
                                    cp_fm_type
   use cp_log_handling,       only: cp_logger_type, &
                                    cp_get_default_logger
   use cp_output_handling,    only: cp_iter_string, &
                                    cp_print_key_unit_nr, &
                                    cp_print_key_finished_output
   use cp_para_types,         only: cp_para_env_type
   use dbcsr_api,             only: dbcsr_copy, &
                                    dbcsr_deallocate_matrix, &
                                    dbcsr_get_info, &
                                    dbcsr_p_type, &
                                    dbcsr_release, &
                                    dbcsr_set
   use input_section_types,   only: section_vals_type, &
                                    section_vals_val_get
   use kinds,                 only: default_string_length, &
                                    dp
   use particle_types,        only: particle_type
   use qs_environment_types,  only: qs_environment_type, &
                                    get_qs_env
   use qs_loc_types,          only: get_qs_loc_env, &
                                    qs_loc_env_new_type
   use qs_moments,            only: build_local_moment_matrix
   use wannier_states_types,  only: wannier_centres_type
#include "./base/base_uses.f90"

   implicit none

   private

! *** Globals ***

   character(len=*), parameter, private :: moduleN = "wannier_moments"
  
! *** Public subroutines ***

   public :: calculate_wannier_moments

! *** body ***
   
contains

! **************************************************************************************************
!> \brief calculates and prints the expectation values of electric moments of Wannier functions
!> \note  only does dipoles for now, but may be easily extended  
! **************************************************************************************************
   subroutine calculate_wannier_moments(qs_env, qs_loc_env, loc_print_key, mos_are_complex)

      type(qs_environment_type), pointer   :: qs_env
      type(qs_loc_env_new_type), pointer   :: qs_loc_env
      type(section_vals_type), pointer     :: loc_print_key
      logical                              :: mos_are_complex
     
      character(len=*), parameter :: routineN = "calculate_wannier_moments", &
                                     routineP = moduleN//":"//routineN

      integer                              :: handle, idim, iion, ispin, istate, naos, norder, &
                                              ndims, nions, nspins, nstates, order, unit_nr, a, b
      integer, allocatable                 :: wcs_ions_map(:)
      real(dp), pointer                    :: rions(:, :), rcentres(:, :)
      character(len=default_string_length) :: iter
      type(cell_type), pointer             :: cell
      type(cp_cfm_type), pointer           :: op_cfm, tmp_coeff, tmp_opvec
      type(cp_cfm_p_type), pointer         :: moments(:), opvec(:)
      type(cp_fm_struct_type), pointer     :: naos_fm_struct, nstates_fm_struct
      type(cp_fm_type), pointer            :: null_fm, null_naos_fm, op_fm
      type(cp_logger_type), pointer        :: logger
      type(cp_para_env_type), pointer      :: para_env
      type(dbcsr_p_type), pointer          :: matrix_s(:), moments_op(:)
      type(dft_control_type), pointer      :: dft_control
      type(particle_type), pointer         :: ions_list(:)
      type(wannier_centres_type), pointer  :: wcs(:)

      call timeset(routineN, handle)

      logger => cp_get_default_logger()
      unit_nr = cp_print_key_unit_nr(logger, loc_print_key, "WANNIER_MOMENTS", extension=".WanMom", &
                                    middle_name="WANNIER_MOMENTS")
      iter = cp_iter_string(logger%iter_info)
      
      ! header line of output
      write (unit=unit_nr, fmt="(A)") "# iter_level, ispin, istate, mapped ion, dipole(x,y,z)[a.u.]"

      call get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, &
                      para_env=para_env, particle_set=ions_list, wanniercentres=wcs)
      nspins = dft_control%nspins
      !print *, "mos_are_complex", mos_are_complex
      nstates = size(wcs(1)%centres, 2)  ! no. of wannier states
      !print *, "nstates", nstates
      naos = qs_loc_env%moloc_coeff(1)%matrix%matrix_struct%nrow_global
      !print *, "naos", naos
      ndims = size(wcs(1)%centres, 1)
      !print *, "ndims", ndims
      nions = size(ions_list)
      !print *, "nions", nions
      call section_vals_val_get(loc_print_key, "WANNIER_MOMENTS%ORDER", i_val=order)
      !print *, "order", order
      CPASSERT(order > 0)

      ! create a temp matrix to restore complex MO coeffs
      call cp_cfm_create(tmp_coeff, qs_loc_env%moloc_coeff(1)%matrix%matrix_struct)
      if (.not. mos_are_complex) then
         ! imaginary part full of zeroes for real inputs
         call cp_fm_create(null_fm, qs_loc_env%moloc_coeff(1)%matrix%matrix_struct)
         call cp_fm_set_all(null_fm, 0.0_dp, 0.0_dp)
      end if

      call cp_fm_struct_create(naos_fm_struct, &
                               template_fmstruct=qs_loc_env%moloc_coeff(1)%matrix%matrix_struct, &
                               nrow_global=naos, ncol_global=naos)
      call cp_fm_struct_create(nstates_fm_struct, &
                               template_fmstruct=qs_loc_env%moloc_coeff(1)%matrix%matrix_struct, &
                               nrow_global=nstates, ncol_global=nstates)
      call cp_fm_create(op_fm, naos_fm_struct)
      call cp_cfm_create(op_cfm, naos_fm_struct)
      call cp_cfm_create(tmp_opvec, qs_loc_env%moloc_coeff(1)%matrix%matrix_struct)
      allocate(moments(ndims), opvec(ndims))
      do idim = 1, ndims
         call cp_cfm_create(moments(idim)%matrix, nstates_fm_struct)
         call cp_cfm_create(opvec(idim)%matrix, qs_loc_env%moloc_coeff(1)%matrix%matrix_struct)
      end do
      ! imaginary part of the operator is 0 since this is in the Schrodinger picture
      call cp_fm_create(null_naos_fm, naos_fm_struct)
      call cp_fm_set_all(null_naos_fm, 0.0_dp, 0.0_dp)

      allocate(rions(ndims, nions), moments_op(ndims))
      do iion = 1, nions
         rions(:, iion) = ions_list(iion)%r
      end do

     do ispin = 1, nspins
        ! prepare coeffs
        if (mos_are_complex) then
           tmp_coeff => qs_loc_env%complex_moloc_coeff(ispin)%matrix
        else
           call cp_fm_to_cfm(qs_loc_env%moloc_coeff(ispin)%matrix, null_fm, tmp_coeff)
        end if
        call cp_cfm_write_formatted(tmp_coeff, unit_nr, "tmp_coeff")

        ! wrap wcs into the cell and map them to "their" ions
        call wrap_wcs_to_cell(wcs(ispin)%centres, cell, rcentres)
        call map_wcs_to_ions(rcentres, rions, wcs_ions_map)

        do istate = 1, nstates
           ! build the moments operator
           do idim = 1, ndims
              allocate(moments_op(idim)%matrix)
              call dbcsr_copy(moments_op(idim)%matrix, matrix_s(1)%matrix, "Wannier moments op")
              call dbcsr_set(moments_op(idim)%matrix, 0.0_dp)
           end do
           call build_local_moment_matrix(qs_env, moments_op, order, ref_point=rcentres(:, istate))

           ! calculate expectation values
           do idim = 1, ndims
              call copy_dbcsr_to_fm(moments_op(idim)%matrix, op_fm)
              call cp_fm_to_cfm(op_fm, null_naos_fm, op_cfm)
              call cp_cfm_gemm("N", "N", naos, nstates, naos, (1.0_dp, 0.0_dp), op_cfm, &
                               tmp_coeff, (0.0_dp, 0.0_dp), tmp_opvec)
              ! only keep the istate-th column
              call cp_cfm_to_cfm(tmp_opvec, opvec(idim)%matrix, 1, istate, istate)
              call cp_cfm_set_all(tmp_opvec, (0.0_dp, 0.0_dp))
           end do
           do idim = 1, ndims
              call cp_cfm_gemm("C", "N", nstates, nstates, naos, (1.0_dp, 0.0_dp), tmp_coeff, &
                               opvec(idim)%matrix, (0.0_dp, 0.0_dp), moments(idim)%matrix)
              call dbcsr_deallocate_matrix(moments_op(idim)%matrix)
           end do

           ! output
           write (unit=unit_nr, fmt="(a,3(i7),3(es18.8))") &
                iter(1:15), ispin, istate, wcs_ions_map(istate), &
                real(moments(1)%matrix%local_data(istate,istate), dp), &
                real(moments(2)%matrix%local_data(istate,istate), dp), &
                real(moments(3)%matrix%local_data(istate,istate), dp)
         end do
      end do

      call cp_print_key_finished_output(unit_nr, logger, loc_print_key, "WANNIER_MOMENTS")

      ! cleanup
      call cp_fm_struct_release(naos_fm_struct)
      call cp_fm_struct_release(nstates_fm_struct)
      call cp_fm_release(op_fm)
      call cp_cfm_release(op_cfm)
      call cp_fm_release(null_fm)
      call cp_fm_release(null_naos_fm)
      do idim = 1, ndims
         call cp_cfm_release(moments(idim)%matrix)
         call cp_cfm_release(opvec(idim)%matrix)
      end do
      deallocate(moments, moments_op, opvec)

      call timestop(handle)

   end subroutine calculate_wannier_moments

! **************************************************************************************************
!> \brief shifts a set of Wannier centres into the simulation cell, only orthorhombic cells
! **************************************************************************************************
   subroutine wrap_wcs_to_cell(r, cell, r_wrapped)

     real(dp), pointer, intent(in)   :: r(:, :)
     type(cell_type), pointer        :: cell
     real(dp), pointer, intent(out)  :: r_wrapped(:, :)

     integer                         :: istate, idim
     real(dp)                        :: abc(3)

     ! assumes an orthorhombic cell
     call get_cell(cell, abc=abc)

     allocate(r_wrapped(size(r, 1), size(r, 2)))

     do istate = 1, size(r, 2)  ! loop over centers
        do idim = 1, size(r, 1)  ! loop over dimensions
           if(r(idim, istate) < 0.0) then
              r_wrapped(idim, istate) = r(idim, istate) + abc(idim)
           else if(r(idim, istate) > abc(idim)) then
              r_wrapped(idim, istate) = r(idim, istate) - abc(idim)
           else
              r_wrapped(idim, istate) = r(idim, istate)
           end if
        end do  ! idim
     end do  ! istate

   end subroutine wrap_wcs_to_cell

! **************************************************************************************************
!> \brief maps the Wannier centres to their nearest ions
! **************************************************************************************************
   subroutine map_wcs_to_ions(rcentres, rions, wcs_ions_map)

     real(dp), pointer, intent(in)      :: rcentres(:, :)
     real(dp), pointer, intent(in)      :: rions(:, :)
     integer, allocatable, intent(out)  :: wcs_ions_map(:)
     
     integer                            :: icentre, idim, iion, ncentres, ndims, nions
     real(dp), pointer                  :: r(:), d(:)

     ndims = size(rcentres, 1)
     ncentres = size(rcentres, 2)
     nions = size(rions, 2)

     allocate(wcs_ions_map(ncentres), r(ndims), d(nions))

     do icentre = 1, ncentres
        do iion = 1, nions
           r = rcentres(:, icentre) - rions(:, iion)
           d(iion) = norm2(r)
        end do
        wcs_ions_map(icentre) = minloc(d, 1)
     end do
     
     deallocate(r, d)
     
   end subroutine map_wcs_to_ions
   
end module wannier_moments

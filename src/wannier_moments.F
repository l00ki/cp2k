!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Routines for the calculation of moments from Wannier functions
!> \author Lukas Schreder
! **************************************************************************************************
module wannier_moments

   use atomic_kind_types,     only: atomic_kind_type, &
                                    get_atomic_kind
   use cell_types,            only: cell_type, &
                                    get_cell, &
                                    pbc
   use cp_cfm_basic_linalg,   only: cp_cfm_gemm, &
                                    cp_cfm_schur_product, &
                                    cp_cfm_transpose
   use cp_cfm_types,          only: cp_cfm_create, &
                                    cp_cfm_get_info, &
                                    cp_cfm_p_type, &
                                    cp_cfm_release, &
                                    cp_cfm_set_all, &
                                    cp_cfm_to_cfm, &
                                    cp_cfm_type, &
                                    cp_cfm_write_formatted, &
                                    cp_fm_to_cfm
   use cp_control_types,      only: dft_control_type
   use cp_dbcsr_operations,   only: copy_dbcsr_to_cfm, &
                                    copy_dbcsr_to_fm, &
                                    cp_dbcsr_sm_cfm_multiply, &
                                    dbcsr_allocate_matrix_set
   use cp_fm_struct,          only: cp_fm_struct_create, &
                                    cp_fm_struct_release, &
                                    cp_fm_struct_type
   use cp_fm_types,           only: cp_fm_create, &
                                    cp_fm_get_info, &
                                    cp_fm_p_type, &
                                    cp_fm_release, &
                                    cp_fm_set_all, &
                                    cp_fm_type
   use cp_log_handling,       only: cp_logger_type, &
                                    cp_get_default_logger
   use cp_output_handling,    only: cp_iter_string, &
                                    cp_p_file, &
                                    cp_print_key_finished_output, &
                                    cp_print_key_should_output, &
                                    cp_print_key_unit_nr
   use cp_para_types,         only: cp_para_env_type
   use dbcsr_api,             only: dbcsr_copy, &
                                    dbcsr_deallocate_matrix, &
                                    dbcsr_get_info, &
                                    dbcsr_p_type, &
                                    dbcsr_release, &
                                    dbcsr_set
   use input_constants,       only: use_mom_ref_com, &
                                    use_mom_ref_coac, &
                                    use_mom_ref_user, &
                                    use_mom_ref_zero
   use input_section_types,   only: section_get_ival, &
                                    section_vals_type, &
                                    section_vals_val_get
   use kinds,                 only: default_string_length, &
                                    dp
   use mathconstants,         only: twopi
   use molecule_types,        only: molecule_type
   use moments_utils,         only: get_reference_point
   use orbital_pointers,      only: current_maxl, &
                                    indco, &
                                    ncoset
   use particle_types,        only: particle_type
   use qs_environment_types,  only: qs_environment_type, &
                                    get_qs_env, &
                                    set_qs_env
   use qs_kind_types,         only: get_qs_kind, &
                                    qs_kind_type
   use qs_loc_molecules,      only: wfc_to_molecule
   use qs_loc_types,          only: get_qs_loc_env, &
                                    qs_loc_env_new_type
   use qs_moments,            only: build_berry_moment_matrix, &
                                    build_local_moment_matrix
   use rt_propagation_types,  only: get_rtp, &
                                    rt_prop_type
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
   subroutine calculate_wannier_moments(qs_env, qs_loc_env, print_key, molecules, mos_are_complex)

      type(qs_environment_type), pointer  :: qs_env
      type(qs_loc_env_new_type), pointer  :: qs_loc_env
      type(section_vals_type), pointer    :: print_key
      type(molecule_type), pointer        :: molecules(:)
      logical                             :: mos_are_complex

      character(len=*), parameter :: routineN = "calculate_wannier_moments", &
           routineP = moduleN//":"//routineN

      character(len=default_string_length)  :: iter
      integer                               :: handle, iion, im, ispin, istate, isubsys, nao, nions, nm, nspins, &
           nstates, nsubsys, order, unit_nr, akind, lx, ly, lz
      integer, allocatable                  :: wcs_subsys_map(:)
      integer, pointer                      :: wcs_ions_map(:)
      logical                               :: first_time, floating, ghost
      real(dp)                              :: zwfc, rcc(3), r(3), ria(3), zeff
      real(dp), allocatable                 :: moments(:, :)
      real(dp), pointer                     :: rcentres(:, :), rions(:, :), ref_point(:)
      type(atomic_kind_type), pointer       :: atomic_kind
      type(cell_type), pointer              :: cell
      type(cp_cfm_type), pointer            :: opvec, vec, momv, coeffs, cfm_moments_op
      type(cp_cfm_p_type), pointer          :: complex_moloc_coeff(:)
      type(cp_fm_struct_type), pointer      :: vec_struct, one_struct, op_struct
      type(cp_fm_p_type), pointer           :: moloc_coeff(:)
      type(cp_fm_type), pointer             :: fm_moments_op
      type(cp_logger_type), pointer         :: logger
      type(dbcsr_p_type), pointer           :: matrix_s(:), moments_op(:)
      type(dft_control_type), pointer       :: dft_control
      type(particle_type), pointer          :: particles(:)
      type(rt_prop_type), pointer           :: rtp
      type(qs_kind_type), pointer           :: qs_kind_set(:)
      type(wannier_centres_type), pointer   :: wcs(:), wos(:)
      
      call timeset(routineN, handle)

      ! output handles
      logger => cp_get_default_logger()
      unit_nr = cp_print_key_unit_nr(logger, print_key, "WANNIER_MOMENTS", extension=".WanMom", &
           middle_name="WANNIER_MOMENTS")
      iter = cp_iter_string(logger%iter_info)
      call section_vals_val_get(print_key, "WANNIER_MOMENTS%ORDER", i_val=order)
      call get_qs_env(qs_env, wanniercentres=wcs)
      if (.not. associated(qs_env%wannierorigins)) then
         call set_qs_env(qs_env, wannierorigins=wcs)
      end if
      call get_qs_env(qs_env, wannierorigins=wos)
      ! sanity checks
      CPASSERT(order > 0)
      CPASSERT(order <= current_maxl)
      CPASSERT(associated(wcs))

      if (BTEST(cp_print_key_should_output(logger%iter_info, print_key, "WANNIER_MOMENTS", &
           first_time=first_time), cp_p_file)) then
         call get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, &
              qs_kind_set=qs_kind_set, rtp=rtp, wcs_ions_map=wcs_ions_map)
         call get_qs_loc_env(qs_loc_env, complex_moloc_coeff=complex_moloc_coeff, &
              moloc_coeff=moloc_coeff, particle_set=particles)
         nions = size(particles)
         nm = ncoset(order) - 1  ! total number of moments (e.g. 3 for dipoles, 12 for qupoles)
         nspins = dft_control%nspins
         nstates = size(wcs(1)%centres, 2)
         nsubsys = size(molecules)
         zwfc = 3.0_dp - real(nspins, dp)

         allocate(rcentres(3, nstates))  ! this holds the wrapped wcs
         call wrap_to_cell(wcs(1)%centres, cell, rcentres)

         allocate(wcs_subsys_map(nstates), rions(3, nions))
         do iion = 1, nions
            rions(:, iion) = particles(iion)%r
         end do
         if (.not. associated(wcs_ions_map)) then
            print *, "new map!"
            call map_wcs_to_ions(rcentres, rions, wcs_ions_map)
            call set_qs_env(qs_env, wcs_ions_map=wcs_ions_map)
         end if
         
         ! map states to subsystems
         do istate = 1, nstates
            do isubsys = 1, nsubsys
               if (wcs_ions_map(istate) >= molecules(isubsys)%first_atom .and. &
                    wcs_ions_map(istate) <= molecules(isubsys)%last_atom) then
                  wcs_subsys_map(istate) = isubsys
               end if
            end do
         end do

         ! prepare coeff and expectation value vectors
         call cp_cfm_create(coeffs, moloc_coeff(1)%matrix%matrix_struct)
         call cp_cfm_get_info(coeffs, nrow_global=nao)
         call cp_fm_struct_create(op_struct, template_fmstruct=coeffs%matrix_struct, ncol_global=nao)
         call cp_fm_struct_create(vec_struct, template_fmstruct=coeffs%matrix_struct, ncol_global=1)
         call cp_fm_struct_create(one_struct, template_fmstruct=vec_struct, nrow_global=1)
         call cp_fm_create(fm_moments_op, op_struct)
         call cp_cfm_create(cfm_moments_op, op_struct)
         call cp_cfm_create(vec, vec_struct)
         call cp_cfm_create(opvec, vec_struct)
         call cp_cfm_create(momv, one_struct)
         
         ! get expectation values
         allocate(moments(nsubsys+1, nm))
         ! moments operator
         allocate(moments_op(nm))
         do im = 1, nm
            allocate(moments_op(im)%matrix)
            call dbcsr_copy(moments_op(im)%matrix, matrix_s(1)%matrix, "moments operator")
            call dbcsr_set(moments_op(im)%matrix, 0.0_dp)
         end do

         do ispin = 1, nspins
            ! get coeffs
            if (mos_are_complex) then
               !call cp_fm_to_cfm(moloc_coeff(2*ispin-1)%matrix, moloc_coeff(2*ispin)%matrix, coeffs)
               call cp_cfm_to_cfm(complex_moloc_coeff(ispin)%matrix, coeffs)
            else
               call cp_fm_to_cfm(moloc_coeff(ispin)%matrix, mtarget=coeffs)
            end if
               
            moments = 0.0_dp
            do istate = 1, nstates
               call get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, ifirst=1, ilast=nions, ref_point=ref_point)
               rcc = pbc(rcc, cell)
               call build_local_moment_matrix(qs_env, moments_op, order, ref_point=rcc)
               ! get the single Wannier orbital
               call cp_cfm_to_cfm(coeffs, vec, 1, istate, 1)
               !print *, "Wannier orbital coeffs", istate
               !call cp_cfm_write_formatted(vec, 6, "vec")
               r = pbc(wcs(ispin)%centres(:, istate), cell)
               do im = 1, nm
                  call copy_dbcsr_to_fm(moments_op(im)%matrix, fm_moments_op)
                  call cp_fm_to_cfm(fm_moments_op, mtarget=cfm_moments_op)
                  call cp_cfm_gemm("N", "N", nao, 1, nao, (1.0_dp, 0.0_dp), cfm_moments_op, vec, &
                       (0.0_dp, 0.0_dp), opvec)
                  call cp_cfm_gemm("C", "N", 1, 1, nao, (1.0_dp, 0.0_dp), vec, opvec, &
                       (0.0_dp, 0.0_dp), momv)
                  moments(wcs_subsys_map(istate), im) = moments(wcs_subsys_map(istate), im) &
                       - zwfc * real(sum(momv%local_data), dp)
                       !- zwfc * r(im)
               end do
            end do

            ! nuclear contribution
            do im = 1, nm
               lx = indco(1, im+1)
               ly = indco(2, im+1)
               lz = indco(3, im+1)
               do iion = 1, nions
                  atomic_kind => particles(iion)%atomic_kind
                  call get_atomic_kind(atomic_kind, kind_number=akind)
                  call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                  if (.not. ghost .and. .not. floating) then
                     call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                     ria = particles(iion)%r - rcc
                     ria = pbc(ria, cell)
                     moments(nsubsys+1, im) = moments(nsubsys+1, im) + zeff * ria(1)**lx * ria(2)**ly * ria(3)**lz
                  end if
               end do
            end do
            
            if (unit_nr > 0) then
               if (first_time) then
                  write(unit_nr, "(a)") "# iter_level, ispin, molecule, moments [a.u.]"
                  !write(unit_nr, fmt="(a10,i3,1X,a20,39es26.16)") iter(1:10), ispin, &
                  !     "nuclear", moments(nsubsys+1, :)
               end if
               do isubsys = 1, nsubsys
                  write(unit_nr, fmt="(a10,i3,1X,a20,39es26.16)") iter(1:10), ispin, &
                       molecules(isubsys)%molecule_kind%name, moments(isubsys, :)
               end do
            end if
         end do

         do im = 1, nm
            call dbcsr_release(moments_op(im)%matrix)
            deallocate(moments_op(im)%matrix)
         end do
         deallocate(moments_op)
         call cp_fm_struct_release(op_struct)
         call cp_fm_struct_release(vec_struct)
         call cp_fm_struct_release(one_struct)
         call cp_fm_release(fm_moments_op)
         call cp_cfm_release(coeffs)
         call cp_cfm_release(cfm_moments_op)
         call cp_cfm_release(vec)
         call cp_cfm_release(opvec)
         call cp_cfm_release(momv)
         deallocate(rcentres, wcs_subsys_map, rions, moments)
      end if

      call cp_print_key_finished_output(unit_nr, logger, print_key, "WANNIER_MOMENTS")
      call timestop(handle)
     
   end subroutine calculate_wannier_moments

! **************************************************************************************************
!> \brief shifts a coordinate into the simulation cell, only orthorhombic cells
! **************************************************************************************************
   subroutine wrap_to_cell(r, cell, r_wrapped)

     real(dp), pointer, intent(in)   :: r(:, :)
     type(cell_type), pointer        :: cell
     real(dp), pointer, intent(out)  :: r_wrapped(:, :)

     integer                         :: idim, istate
     real(dp)                        :: abc(3)

     ! assumes an orthorhombic cell
     call get_cell(cell, abc=abc)

     do istate = 1, size(r, 2)
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

   end subroutine wrap_to_cell

   subroutine wrap_to_cell_1d(r, cell, r_wrapped)

      real(dp), pointer, intent(in)   :: r(:)
      type(cell_type), pointer        :: cell
      real(dp), intent(out)           :: r_wrapped(:)

      integer                         :: idim
      real(dp)                        :: abc(3)

      ! assumes an orthorhombic cell
      call get_cell(cell, abc=abc)

      do idim = 1, size(r, 1)  ! loop over dimensions
         if(r(idim) < 0.0) then
            r_wrapped(idim) = r(idim) + abc(idim)
         else if(r(idim) > abc(idim)) then
            r_wrapped(idim) = r(idim) - abc(idim)
         else
            r_wrapped(idim) = r(idim)
         end if
      end do  ! idim

   end subroutine wrap_to_cell_1d

   subroutine wrap_to_cell_d(r, cell, r_wrapped)

      real(dp), intent(in)      :: r(:)
      type(cell_type), pointer  :: cell
      real(dp), intent(out)     :: r_wrapped(:)

      integer                         :: idim
      real(dp)                        :: abc(3)

      ! assumes an orthorhombic cell
      call get_cell(cell, abc=abc)

      do idim = 1, size(r, 1)  ! loop over dimensions
         if(r(idim) < 0.0) then
            r_wrapped(idim) = r(idim) + abc(idim)
         else if(r(idim) > abc(idim)) then
            r_wrapped(idim) = r(idim) - abc(idim)
         else
            r_wrapped(idim) = r(idim)
         end if
      end do  ! idim

   end subroutine wrap_to_cell_d

! **************************************************************************************************
!> \brief maps the Wannier centres to their nearest ions
! **************************************************************************************************
   subroutine map_wcs_to_ions(rcentres, rions, wcs_ions_map)

     real(dp), pointer, intent(in)      :: rcentres(:, :)
     real(dp), pointer, intent(in)      :: rions(:, :)
     integer, pointer, intent(out)      :: wcs_ions_map(:)
     
     integer                            :: icentre, iion, ncentres, ndims, nions
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

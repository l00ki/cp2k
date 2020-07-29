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
   use cp_cfm_basic_linalg,   only: cp_cfm_gemm
   use cp_cfm_types,          only: cp_cfm_create, &
                                    cp_cfm_p_type, &
                                    cp_cfm_release, &
                                    cp_cfm_to_cfm, &
                                    cp_cfm_type, &
                                    cp_cfm_write_formatted, &
                                    cp_fm_to_cfm
   use cp_control_types,      only: dft_control_type
   use cp_dbcsr_operations,   only: copy_dbcsr_to_fm, &
                                    cp_dbcsr_sm_fm_multiply
   use cp_fm_basic_linalg,    only: cp_fm_schur_product
   use cp_fm_struct,          only: cp_fm_struct_create, &
                                    cp_fm_struct_release, &
                                    cp_fm_struct_type
   use cp_fm_types,           only: cp_fm_create, &
                                    cp_fm_get_info, &
                                    cp_fm_p_type, &
                                    cp_fm_release, &
                                    cp_fm_to_fm, &
                                    cp_fm_type, &
                                    cp_fm_write_formatted
   use cp_log_handling,       only: cp_get_default_logger, &
                                    cp_logger_type, &
                                    cp_to_string
   use cp_output_handling,    only: cp_iter_string, &
                                    cp_p_file, &
                                    cp_print_key_finished_output, &
                                    cp_print_key_should_output, &
                                    cp_print_key_unit_nr
   use cp_para_types,         only: cp_para_env_type
   use dbcsr_api,             only: dbcsr_copy, &
                                    dbcsr_deallocate_matrix,&
                                    dbcsr_p_type, &
                                    dbcsr_set
   use input_constants,       only: use_mom_ref_com
   use input_section_types,   only: section_vals_type, &
                                    section_vals_val_get
   use kinds,                 only: default_string_length, &
                                    dp
   use message_passing,       only: mp_bcast, &
                                    mp_maxloc, &
                                    mp_sum
   use molecule_types,        only: get_molecule_set_info, &
                                    molecule_type
   use moments_utils,         only: get_reference_point
   use orbital_pointers,      only: indco, &
                                    ncoset
   use particle_types,        only: particle_type
   use qs_environment_types,  only: get_qs_env, &
                                    qs_environment_type, &
                                    set_qs_env
   use qs_kind_types,         only: get_qs_kind, &
                                    qs_kind_type
   use qs_loc_types,          only: get_qs_loc_env, &
                                    qs_loc_env_new_type
   use qs_moments,            only: build_local_moment_matrix
   use rt_propagation_types,  only: get_rtp, &
                                    rt_prop_type, &
                                    rtp_states_type
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
   subroutine calculate_wannier_moments(qs_env, qs_loc_env, print_key, molecules, do_rtp)

      type(qs_environment_type), pointer  :: qs_env
      type(qs_loc_env_new_type), pointer  :: qs_loc_env
      type(section_vals_type), pointer    :: print_key
      type(molecule_type), pointer        :: molecules(:)
      logical                             :: do_rtp

      character(len=*), parameter :: routineN = "calculate_wannier_moments", &
                                     routineP = moduleN//":"//routineN

      character(len=default_string_length)  :: iter, name, nameheader
      integer                               :: akind, handle, i, j, k, l, lx, ly, lz, m, nao, natoms, &
           nlocalstates(2), nmo, nmoments, norders, nspins, nstates, nsubsets, unit_nr
      integer, allocatable                  :: first_atoms(:), last_atoms(:), nbasis(:), states(:)
      integer, pointer                      :: current_states(:)
      logical                               :: do_periodic, do_orbitals, first_time, floating, ghost
      real(dp)                              :: abc(3), r(3), rcc(3), zeff, zwfc
      real(dp), allocatable                 :: moments(:, :)
      type(atomic_kind_type), pointer       :: atomic_kind
      type(cell_type), pointer              :: cell
      type(cp_cfm_p_type), pointer          :: coeffs(:), moments_op(:)
      type(cp_cfm_type), pointer            :: vec, opvec, vecopvec
      type(cp_fm_struct_type), pointer      :: op_struct, vec_struct, vecopvec_struct
      type(cp_fm_p_type), pointer           :: tmp_coeffs(:)
      type(cp_fm_type), pointer             :: opbuffer
      type(cp_logger_type), pointer         :: logger
      type(cp_para_env_type), pointer       :: para_env
      type(dbcsr_p_type), pointer           :: matrix_s(:), tmp_moments_op(:)
      type(dft_control_type), pointer       :: dft_control
      type(particle_type), pointer          :: particles(:)
      type(qs_kind_type), pointer           :: kinds(:)
      type(rt_prop_type), pointer           :: rtp
      !type(rtp_states_type), pointer        :: states(:)
      type(wannier_centres_type), pointer   :: origins(:), wcs(:), wos(:)

      call timeset(routineN, handle)

      logger => cp_get_default_logger()

      call section_vals_val_get(print_key, "WANNIER_MOMENTS%ORDER", i_val=norders)
      call section_vals_val_get(print_key, "WANNIER_MOMENTS%PERIODIC", l_val=do_periodic)
      call section_vals_val_get(print_key, "WANNIER_MOMENTS%ORBITALS", l_val=do_orbitals)

      ! sanity checks
      CPASSERT(norders > 0)

      ! get parameters
      call get_qs_env(qs_env, cell=cell, dft_control=dft_control, matrix_s=matrix_s, &
           particle_set=particles, qs_kind_set=kinds, wanniercentres=wcs, wannierorigins=wos)
      if (do_rtp) then
         call get_qs_env(qs_env, rtp=rtp)
      end if
      
      nspins = dft_control%nspins
      zwfc = 3.0_dp - real(nspins, dp)  ! occupation of one Wannier function
      nmoments = ncoset(norders)  ! includes the 0th order i.e. charge
      nstates = size(wcs(1)%centres, 2)
      if (do_orbitals) then
         nsubsets = size(wcs(1)%centres(1, :))
         nameheader = "Wannier function"
      else
         nsubsets = size(molecules)
         nameheader = "molecule"
         allocate(first_atoms(nsubsets), last_atoms(nsubsets), nbasis(nsubsets))
         call get_molecule_set_info(molecules, mol_to_first_atom=first_atoms, &
              mol_to_last_atom=last_atoms, mol_to_nbasis=nbasis)
      end if
      if (.not. associated(wos)) then
         ! save the SCF Wannier centres
         allocate(origins(nspins))
         do i = 1, nspins
            allocate(origins(i)%centres(3, nstates))
            print *, "i, wcs(i) dimensions, nstates", i, size(wcs(i)%centres, 1), size(wcs(i)%centres, 2), nstates
            origins(i)%centres = wcs(i)%centres
            allocate(origins(i)%wannierhamdiag(nstates))
            origins(i)%wannierhamdiag = wcs(i)%wannierhamdiag
         end do
         call set_qs_env(qs_env, wannierorigins=origins)
         wos => origins
      end if
      call get_cell(cell, abc=abc)

      call get_qs_loc_env(qs_loc_env, para_env=para_env, moloc_coeff=tmp_coeffs)
      call cp_fm_get_info(tmp_coeffs(1)%matrix, ncol_global=nmo, nrow_global=nao)
      allocate(coeffs(nspins))
      do i = 1, nspins
         call cp_cfm_create(coeffs(i)%matrix, tmp_coeffs(1)%matrix%matrix_struct)
         if (do_rtp) then
            call cp_fm_to_cfm(tmp_coeffs(2*i-1)%matrix, tmp_coeffs(2*i)%matrix, coeffs(i)%matrix)
         else
            call cp_fm_to_cfm(tmp_coeffs(i)%matrix, mtarget=coeffs(i)%matrix)
         end if
      end do

      allocate(moments(nmoments, nsubsets))
      moments = 0.0_dp

      if (btest(cp_print_key_should_output(logger%iter_info, print_key, "WANNIER_MOMENTS", &
           first_time=first_time), cp_p_file)) then
         call cp_fm_struct_create(vec_struct, nrow_global=nao, ncol_global=1, &
              para_env=coeffs(1)%matrix%matrix_struct%para_env, &
              context=coeffs(1)%matrix%matrix_struct%context)
         call cp_fm_struct_create(op_struct, nrow_global=nao, ncol_global=nao, &
              para_env=coeffs(1)%matrix%matrix_struct%para_env, &
              context=coeffs(1)%matrix%matrix_struct%context)
         call cp_fm_struct_create(vecopvec_struct, nrow_global=1, ncol_global=1, &
              para_env=coeffs(1)%matrix%matrix_struct%para_env, &
              context=coeffs(1)%matrix%matrix_struct%context)
         call cp_cfm_create(vec, vec_struct)
         call cp_cfm_create(opvec, vec_struct)
         call cp_cfm_create(vecopvec, vecopvec_struct)
         allocate(tmp_moments_op(nmoments-1))
         allocate(moments_op(nmoments-1))
         call cp_fm_create(opbuffer, op_struct)
         do i = 1, nmoments - 1
            allocate(tmp_moments_op(i)%matrix)
            call dbcsr_copy(tmp_moments_op(i)%matrix, matrix_s(1)%matrix, "Moments operator")
            call dbcsr_set(tmp_moments_op(i)%matrix, 0.0_dp)
            call cp_cfm_create(moments_op(i)%matrix, op_struct)
         end do
         
         do i = 1, nspins
            do j = 1, nsubsets
               if (do_orbitals) then
                  allocate(states(1))
                  states(1) = j
                  nstates = 1
               else
                  ! mpi-safe state list
                  if (associated(molecules(j)%lmi)) then
                     nlocalstates(1) = molecules(j)%lmi(i)%nstates
                  else
                     nlocalstates(1) = 0
                  end if
                  nlocalstates(2) = para_env%mepos
                  call mp_maxloc(nlocalstates, para_env%group)
                  k = nlocalstates(1)
                  if (k == 0) cycle
                  allocate(states(k))
                  nstates = k
                  if (nlocalstates(2) == para_env%mepos) then
                     states(:) = molecules(k)%lmi(i)%states(:)
                  else
                     states(:) = 0
                  end if
                  call mp_bcast(states, nlocalstates(2), para_env%group)
               end if
 
               ! get reference point
               if (do_orbitals) then
                  rcc = wrap(wos(i)%centres(:, j), cell)
               else
                  call get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
                       ifirst=first_atoms(j), ilast=last_atoms(j))
                       natoms = last_atoms(j) - first_atoms(j) + 1
               end if

               ! electronic contribution
               if (do_periodic) then
                  ! HACK
                  do k = 1, nstates
                     moments(2:4, j) = moments(2:4, j) &
                          - zwfc * (wrap(wcs(i)%centres(:, states(k)), cell) - wrap(rcc, cell))
                  end do
               else
                  do k = 1, nstates 
                     call cp_cfm_to_cfm(coeffs(i)%matrix, vec, 1, states(k), 1)
                     call build_local_moment_matrix(qs_env, tmp_moments_op, norders, wrap(rcc, cell))
                     do l = 1, nmoments - 1
                        call copy_dbcsr_to_fm(tmp_moments_op(l)%matrix, opbuffer)
                        call cp_fm_to_cfm(opbuffer, mtarget=moments_op(l)%matrix)
                        call cp_cfm_gemm("N", "N", nao, 1, nao, (1.0_dp, 0.0_dp), moments_op(l)%matrix, vec, &
                             (0.0_dp, 0.0_dp), opvec)
                        call cp_cfm_gemm("C", "N", 1, 1, nao, (1.0_dp, 0.0_dp), vec, &
                             opvec, (0.0_dp, 0.0_dp), vecopvec)
                        call cp_cfm_write_formatted(vecopvec, 6, "vecopvec"//cp_to_string(l))
                        moments(l+1, j) = moments(l+1, j) - zwfc * sum(real(vecopvec%local_data))
                     end do
                  end do
               end if

               ! local nuclear contribution
               if (.not. do_orbitals) then
                  do k = 2, nmoments
                     lx = indco(1, k)
                     ly = indco(2, k)
                     lz = indco(3, k)
                     do l = 1, natoms
                        m = first_atoms(j) + l - 1
                        atomic_kind => particles(m)%atomic_kind
                        call get_atomic_kind(atomic_kind, kind_number=akind)
                        call get_qs_kind(kinds(akind), ghost=ghost, floating=floating)
                        if (.not. ghost .and. .not. floating) then
                           call get_qs_kind(kinds(akind), core_charge=zeff)
                           r = particles(m)%r - rcc
                           r = pbc(r, cell)
                           if (lx /= 0) zeff = zeff * r(1)**lx
                           if (ly /= 0) zeff = zeff * r(2)**ly
                           if (lz /= 0) zeff = zeff * r(3)**lz
                           moments(k, j) = moments(k, j) + zeff
                        end if
                     end do  ! atoms
                  end do  ! moments
               end if

               ! cleanup
               deallocate(states)
            end do  ! subsets

            !! total nuclear contribution
            !if (do_orbitals) then
            !   do k = 2, nmoments
            !      lx = indco(1, k)
            !      ly = indco(2, k)
            !      lz = indco(3, k)
            !      do l = 1, natoms
            !         m = first_atoms(j) + l - 1
            !         atomic_kind => particles(m)%atomic_kind
            !         call get_atomic_kind(atomic_kind, kind_number=akind)
            !         call get_qs_kind(kinds(akind), ghost=ghost, floating=floating)
            !         if (.not. ghost .and. .not. floating) then
            !            call get_qs_kind(kinds(akind), core_charge=zeff)
            !            r = particles(m)%r - rcc
            !            r = pbc(r, cell)
            !            if (lx /= 0) zeff = zeff * r(1)**lx
            !            if (ly /= 0) zeff = zeff * r(2)**ly
            !            if (lz /= 0) zeff = zeff * r(3)**lz
            !            moments(k, j) = moments(k, j) + zeff
            !         end if
            !      end do  ! atoms
            !   end do  ! moments
            !end if
         end do  ! spins

         ! total nuclear contribution
         if (do_orbitals) then
         end if

         ! cleanup
         call cp_fm_struct_release(vec_struct)
         call cp_fm_struct_release(op_struct)
         call cp_fm_struct_release(vecopvec_struct)
         call cp_cfm_release(vec)
         call cp_cfm_release(opvec)
         call cp_cfm_release(vecopvec)
         call cp_fm_release(opbuffer)
         do i = 1, nmoments - 1
            call dbcsr_deallocate_matrix(tmp_moments_op(i)%matrix)
            call cp_cfm_release(moments_op(i)%matrix)
         end do
         deallocate(tmp_moments_op)
         deallocate(moments_op)
      end if

      unit_nr = cp_print_key_unit_nr(logger, print_key, "WANNIER_MOMENTS", extension=".WanMom", &
            middle_name="WANNIER_MOMENTS")
      iter = cp_iter_string(logger%iter_info)
      if (unit_nr > 0) then
         if (first_time) then
            write(unit_nr, "(a)") "# iter_level, "//trim(nameheader)//", isubset, moments [a.u.]"
         end if
         do i = 1, nsubsets
            if (do_orbitals) then
               name = cp_to_string(i)
            else
               name = molecules(i)%molecule_kind%name
            end if
            write (unit_nr, "(a10,a20,39es26.16)") iter(1:10), name, moments(2:nmoments, i)
         end do
      end if
      call cp_print_key_finished_output(unit_nr, logger, print_key, "WANNIER_MOMENTS")

      ! cleanup
      deallocate(moments)
      do i = 1, size(coeffs)
         call cp_cfm_release(coeffs(i)%matrix)
      end do
      deallocate(coeffs)
      if (.not. do_orbitals) deallocate(first_atoms, last_atoms, nbasis)
      
      call timestop(handle)

   end subroutine calculate_wannier_moments

! **************************************************************************************************
! **************************************************************************************************
   function wrap(rin, cell) result(rout)

     real(dp), intent(in)      :: rin(3)
     type(cell_type), pointer  :: cell
     real(dp)                  :: rout(3)

     integer   :: i
     real(dp)  :: abc(3)

     CPASSERT(cell%orthorhombic)
     call get_cell(cell, abc=abc)
     rout = rin
     do i = 1, 3
        do while (rout(i) < 0.0_dp)
           rout(i) = rout(i) + abc(i)
        end do
        do while (rout(i) > abc(i))
           rout(i) = rout(i) - abc(i)
        end do
     end do

   end function wrap
   
end module wannier_moments

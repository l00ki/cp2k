!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Routines for the calculation of moments from Wannier functions
!> \author Lukas Schreder
! **************************************************************************************************
module localized_moments

   use atomic_kind_types,      only: atomic_kind_type, &
                                     get_atomic_kind
   use cell_types,             only: cell_type, &
                                     pbc
   use cp_cfm_basic_linalg,    only: cp_cfm_gemm, &
                                     cp_cfm_schur_product
   use cp_cfm_types,           only: cp_cfm_create, &
                                     cp_cfm_get_element, &
                                     cp_cfm_p_type, &
                                     cp_cfm_release, &
                                     cp_cfm_set_all, &
                                     cp_cfm_to_cfm, &
                                     cp_cfm_type, &
                                     cp_fm_to_cfm
   use cp_control_types,       only: dft_control_type
   use cp_dbcsr_operations,    only: copy_dbcsr_to_fm, &
                                     cp_dbcsr_sm_fm_multiply
   use cp_fm_basic_linalg,     only: cp_fm_gemm, &
                                     cp_fm_schur_product
   use cp_fm_struct,           only: cp_fm_struct_create, &
                                     cp_fm_struct_release, &
                                     cp_fm_struct_type
   use cp_fm_types,            only: cp_fm_create, &
                                     cp_fm_get_element, &
                                     cp_fm_get_info, &
                                     cp_fm_p_type, &
                                     cp_fm_write_formatted, &
                                     cp_fm_release, &
                                     cp_fm_set_all, &
                                     cp_fm_to_fm, &
                                     cp_fm_type
   use cp_log_handling,        only: cp_get_default_logger, &
                                     cp_logger_type
   use cp_output_handling,     only: cp_iter_string, &
                                     cp_p_file, &
                                     cp_print_key_finished_output, &
                                     cp_print_key_should_output, &
                                     cp_print_key_unit_nr
   use cp_para_types,          only: cp_para_env_type
   use dbcsr_api,              only: dbcsr_copy, &
                                     dbcsr_dot, &
                                     dbcsr_get_block_p, &
                                     dbcsr_get_diag, &
                                     dbcsr_iterator_start, &
                                     dbcsr_iterator_next_block, &
                                     dbcsr_iterator_blocks_left, &
                                     dbcsr_iterator_type, &
                                     dbcsr_multiply, &
                                     dbcsr_p_type, &
                                     dbcsr_release, &
                                     dbcsr_set, &
                                     dbcsr_type
   use distribution_1d_types,  only: distribution_1d_type
   use input_constants,        only: use_mom_ref_com
   use input_section_types,    only: section_vals_type, &
                                     section_vals_val_get
   use kinds,                  only: default_string_length, &
                                     dp
   use mathconstants,          only: gaussi, &
                                     twopi
   use message_passing,        only: mp_bcast, &
                                     mp_maxloc, &
                                     mp_sum
   use molecule_kind_types,    only: molecule_kind_type
   use molecule_types,         only: get_molecule_set_info, &
                                     molecule_type
   use moments_utils,          only: get_reference_point
   use orbital_pointers,       only: indco, &
                                     ncoset
   use particle_types,         only: particle_type
   use qs_environment_types,   only: get_qs_env, &
                                     persistent_states_type, &
                                     qs_environment_type
   use qs_kind_types,          only: get_qs_kind, &
                                     qs_kind_type
   use qs_loc_types,           only: get_qs_loc_env, &
                                     qs_loc_env_new_type
   use qs_moments,             only: build_local_moment_matrix
   use qs_rho_types,           only: qs_rho_get, &
                                     qs_rho_type
   use wannier_states_types,   only: wannier_centres_type

#include "./base/base_uses.f90"

   implicit none

   private

! *** globals ***

   character(len=*), parameter, private :: moduleN = "localized_moments"

! *** public subroutines ***

   public :: calculate_kg_moments, &
             calculate_wannier_moments, &
             calculate_wc_moments

! *** body ***

contains

! **************************************************************************************************
!> \brief Calculates localised electric moments from Kim-Gordon subdensities
!> \param qs_env the QuickStep environment
!> \param norders the number of multipole orders
!> \param print_key the force_eval section
! **************************************************************************************************
   subroutine calculate_kg_moments(qs_env, norders, print_key)

     type(qs_environment_type), pointer  :: qs_env
     integer                             :: norders
     type(section_vals_type), pointer    :: print_key

     character(len=*), parameter  :: routineN = "calculate_kg_moments", &
                                     routineP = moduleN//":"//routineN

     character(len=default_string_length)  :: iter, name
     integer                               :: akind, first_atom, group, handle, iatom, jatom, blk, &
          ilocmol, imolkind, imoment, &
          ispin, isubset, last_atom, mepos, natoms, &
          nlocmols, nmolkinds, nmoments, nspins, nsubsets, output_unit
     integer, allocatable                  :: first_atoms(:), first_basis(:), last_atoms(:), &
          last_basis(:), lx, ly, lz, nbasis(:), nelectrons(:)
     logical                               :: first_time, floating, ghost, found, tr
     real(kind=dp)                         :: focc, m, mtot, r(3), zeff, zmom, d(3)
     real(kind=dp), allocatable            :: diag(:), moments(:, :), rccs(:, :)
     real(kind=dp), pointer                :: ref_point(:), blk_p(:)
     type(atomic_kind_type), pointer       :: atomic_kind
     type(cell_type), pointer              :: cell
     type(cp_logger_type), pointer         :: logger
     type(cp_para_env_type), pointer       :: para_env
     type(dbcsr_p_type), pointer           :: densmat(:), matrix_s(:), moments_op(:)
     type(dbcsr_type), pointer             :: opmat, local_densmat
     type(dbcsr_iterator_type), pointer    :: iterator
     type(dft_control_type), pointer       :: dft_control
     type(distribution_1d_type), pointer   :: local_molecules
     type(molecule_kind_type), pointer     :: molkind
     type(molecule_type), pointer          :: molecules(:)
     type(particle_type), pointer          :: particles(:)
     type(qs_kind_type), pointer           :: qs_kind_set(:)
     type(qs_rho_type), pointer            :: rho
    
     call timeset(routineN, handle)

     logger => cp_get_default_logger()

     if (btest(cp_print_key_should_output(logger%iter_info, print_key, "DFT%PRINT%MOMENTS", first_time=first_time), &
                                          cp_p_file)) then

        ! user input
        CPASSERT(norders >= 0)
        nmoments = ncoset(norders) - 1

        ! environment infos
        call get_qs_env(qs_env, cell=cell, dft_control=dft_control, &
             matrix_s=matrix_s, para_env=para_env, particle_set=particles, &
             qs_kind_set=qs_kind_set, local_molecules=local_molecules, &
             molecule_set=molecules, rho=rho)
        mepos = para_env%mepos
        group = para_env%group
        nspins = dft_control%nspins
        nsubsets = size(molecules)
        call qs_rho_get(rho, rho_ao=densmat)
        allocate(moments(nmoments, nsubsets), first_atoms(nsubsets), first_basis(nsubsets), &
             last_atoms(nsubsets), last_basis(nsubsets), &
             rccs(3, nsubsets), moments_op(nmoments), nbasis(nsubsets), ref_point(3), nelectrons(nsubsets))
        moments = 0._dp
        call get_molecule_set_info(molecules, &
             mol_to_first_atom=first_atoms, mol_to_last_atom=last_atoms, mol_to_nbasis=nbasis, &
             mol_to_nelectrons=nelectrons)
        ! get lists of the first and last basis functions of each molecule
        first_basis(1) = 1
        last_basis(1) = nbasis(1)
        do isubset = 2, nsubsets
           first_basis(isubset) = sum(nbasis(1:isubset-1)) + 1
           last_basis(isubset) = sum(nbasis(1:isubset))
        end do
        allocate(diag(sum(nbasis)))
        ! create the operator
        do imoment = 1, nmoments
           allocate(moments_op(imoment)%matrix)
           call dbcsr_copy(moments_op(imoment)%matrix, matrix_s(1)%matrix, "moments operator")
           call dbcsr_set(moments_op(imoment)%matrix, 0._dp)
        end do
        allocate(local_densmat)
        call dbcsr_copy(local_densmat, matrix_s(1)%matrix, "local densmat")
        call dbcsr_set(local_densmat, 0._dp)
        nmolkinds = size(local_molecules%n_el)

        do isubset = 1, nsubsets
           ! gives wrong results for some reason
           call get_reference_point(rccs(:, isubset), qs_env=qs_env, reference=use_mom_ref_com, &
                ref_point=ref_point, &!ifirst=1, ilast=3)!, &
                ifirst=first_atoms(isubset), ilast=last_atoms(isubset))

           call build_local_moment_matrix(qs_env, moments_op, norders, pbc(rccs(:, isubset), cell))

           ! nuclear dipole
           do imoment = 1, nmoments
              lx = indco(1, imoment + 1)
              ly = indco(2, imoment + 1)
              lz = indco(3, imoment + 1)
              do iatom = first_atoms(isubset), last_atoms(isubset)
                 atomic_kind => particles(iatom)%atomic_kind
                 call get_atomic_kind(atomic_kind, kind_number=akind)
                 call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating = floating)
                 if (.not. ghost .and. .not. floating) then
                    call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                    r = particles(iatom)%r - rccs(:, isubset)
                    r = pbc(r, cell)
                    zmom = zeff
                    if (lx /= 0) zmom = zmom * r(1)**lx
                    if (ly /= 0) zmom = zmom * r(2)**ly
                    if (lz /= 0) zmom = zmom * r(3)**lz
                    moments(imoment, isubset) = moments(imoment, isubset) + zmom
                 end if
              end do
           end do

           ! Kim-Gordon subdensity moments
           focc = (3.0_dp - real(nspins, dp)) / 2.0_dp  ! uniform occupation
           do ispin = 1, nspins
              call dbcsr_copy(local_densmat, densmat(ispin)%matrix)
              allocate(iterator)
              call dbcsr_iterator_start(iterator, local_densmat)
              do while (dbcsr_iterator_blocks_left(iterator))
                 call dbcsr_iterator_next_block(iterator, iatom, jatom, blk)
                 if (((iatom < first_atoms(isubset)) .or. (iatom > last_atoms(isubset))) .or. &
                      (jatom < first_atoms(isubset)) .or. (jatom > last_atoms(isubset))) then
                    call dbcsr_get_block_p(local_densmat, iatom, jatom, blk_p, tr, found)
                    blk_p = 0.0_dp
                 end if
              end do
              d = 0.0_dp
              do imoment = 1, nmoments
                 call dbcsr_dot(local_densmat, moments_op(imoment)%matrix, d(imoment))
                 d(imoment) = cell%hmat(imoment, imoment) / twopi * &
                      aimag(log(exp(-gaussi * twopi / cell%hmat(imoment, imoment) * focc * d(imoment))))
                 moments(imoment, isubset) = moments(imoment, isubset) + d(imoment)
              end do
           end do
        end do

        ! output
        output_unit = cp_print_key_unit_nr(logger, print_key, "DFT%PRINT%MOMENTS", &
             extension=".KGMom", middle_name="KG_MOMENTS")
        iter = cp_iter_string(logger%iter_info)
        if (output_unit > 0) then
           if (first_time) then
              write(output_unit, "(a)") "# iter_level, subset, moments [a.u.]"
           end if
           do isubset = 1, nsubsets
              name = molecules(isubset)%molecule_kind%name
              write(output_unit, "(a10,a10,39f12.6)") iter(1:10), name, moments(:, isubset)
           end do
        end if

        call cp_print_key_finished_output(output_unit, logger, print_key, "DFT%PRINT%MOMENTS")

     end if
  
     call timestop(handle)

   end subroutine calculate_kg_moments

! **************************************************************************************************
!> \brief Calculates localised eletric moments from Wannier functions
!> \param qs_env the QuickStep environment
!> \param qs_loc_env the localisation environment
!> \param mo_local the Wannier function coefficients
!> \param print_key the print section
!> \param molecules the molecule list of the system
! **************************************************************************************************
   subroutine calculate_wannier_moments(qs_env, qs_loc_env, mo_local, print_key, molecules)

     type(qs_environment_type), pointer  :: qs_env
     type(qs_loc_env_new_type), pointer  :: qs_loc_env
     type(cp_fm_p_type), pointer         :: mo_local(:)
     type(section_vals_type), pointer    :: print_key
     type(molecule_type), pointer        :: molecules(:)

     character(len=*), parameter  :: routineN = "calculate_wannier_moments", &
                                     routineP = moduleN//":"//routineN

     character(len=default_string_length)  :: iter, name
     complex(kind=dp)                      :: clocd
     integer                               :: akind, group, handle, i, iatom, &
          imoment, imo, &
          iproc, is, ispin, isubset, jatom, lx, ly, lz, &
          mepos, &
          nao, natoms, nmo, nmolkinds, nmoments, norders, ns, &
          nspins, nstates(2), nsubsets, nwc, output_unit
     integer, allocatable                   :: first_atoms(:), last_atoms(:), atom_to_mol(:), &
          mol_to_nelectrons(:), wcs_to_atoms(:, :)
     integer, pointer                       :: states(:)
     logical                                :: do_rtp, first_time, floating, ghost
     real(kind=dp)                          :: d(3), dlocal, r(3), zeff, zmom, zwfc, &
          m, mtot, rcc(3), locd, r_min
     real(kind=dp), allocatable             :: moments(:, :), rccs(:, :)
     real(kind=dp), pointer                 :: ref_point(:)
     type(atomic_kind_type), pointer        :: atomic_kind
     type(cell_type), pointer               :: cell
     type(cp_cfm_p_type), pointer           :: cmoments_op_fm(:)
     type(cp_cfm_type), pointer             :: ccoeffs, cvec, copvec, cvecopvec
     type(cp_fm_p_type), pointer            :: moments_op_fm(:)
     type(cp_fm_struct_type), pointer       :: fm_struct, mstruct, momstruct, opstruct, nmostruct
     type(cp_fm_type), pointer              :: coeffs, momvector, mvector, &
                                               omvector, opvec, vec, vecopvec
     type(cp_logger_type), pointer          :: logger
     type(cp_para_env_type), pointer        :: para_env
     type(dbcsr_p_type), pointer            :: matrix_s(:), moments_op(:)
     type(dft_control_type), pointer        :: dft_control
     type(distribution_1d_type), pointer    :: local_molecules
     type(particle_type), pointer           :: particles(:)
     type(persistent_states_type), pointer  :: persistent_states(:, :)
     type(qs_kind_type), pointer            :: qs_kind_set(:)
     type(wannier_centres_type), pointer    :: wcs(:)

     call timeset(routineN, handle)

     logger => cp_get_default_logger()

     if (btest(cp_print_key_should_output(logger%iter_info, print_key, "WANNIER_MOMENTS", &
                                          first_time=first_time), cp_p_file)) then

        ! user input
        call section_vals_val_get(print_key, "WANNIER_MOMENTS%ORDER", i_val=norders)
        CPASSERT(norders >= 0)
        nmoments = ncoset(norders) - 1

        ! environment infos
        call get_qs_env(qs_env, cell=cell, dft_control=dft_control, &
             matrix_s=matrix_s, para_env=para_env, particle_set=particles, qs_kind_set=qs_kind_set, &
             wanniercentres=wcs)
        mepos = para_env%mepos
        group = para_env%group
        nspins = dft_control%nspins
        do_rtp = (size(mo_local) == 2 * nspins)
        nsubsets = size(molecules)
        natoms = size(particles)
        do_rtp = (size(mo_local) == 2 * nspins)  ! check for time-dependent orbitals
        zwfc = 3.0_dp - real(nspins, dp)  ! occupation of Wannier functions (uniform)
        allocate(moments(nmoments, nsubsets), first_atoms(nsubsets), &
             last_atoms(nsubsets), rccs(3, nsubsets), moments_op(nmoments), moments_op_fm(nmoments), ref_point(3), &
             cmoments_op_fm(nmoments), atom_to_mol(natoms), mol_to_nelectrons(nsubsets))
        moments = 0.0_dp
        call get_molecule_set_info(molecules, &
             mol_to_first_atom=first_atoms, mol_to_last_atom=last_atoms, atom_to_mol=atom_to_mol, &
             mol_to_nelectrons=mol_to_nelectrons)
        coeffs => mo_local(1)%matrix

        call cp_fm_get_info(coeffs, nrow_global=nao, ncol_global=nmo)
        call cp_fm_struct_create(nmostruct, ncol_global=nmo, nrow_global=nmo, &
             para_env=coeffs%matrix_struct%para_env, context=coeffs%matrix_struct%context)
        call cp_fm_struct_create(opstruct, ncol_global=nao, nrow_global=nao, &
             para_env=coeffs%matrix_struct%para_env, context=coeffs%matrix_struct%context)
        if (do_rtp) then
           call cp_cfm_create(cvec, coeffs%matrix_struct)
           call cp_cfm_create(copvec, coeffs%matrix_struct)
           call cp_cfm_create(cvecopvec, nmostruct)
           call cp_cfm_create(ccoeffs, coeffs%matrix_struct)
        else
           call cp_fm_create(vec, coeffs%matrix_struct)
           call cp_fm_create(opvec, coeffs%matrix_struct)
           call cp_fm_create(vecopvec, nmostruct)
        end if
        

        do imoment = 1, nmoments
           allocate(moments_op(imoment)%matrix)
           call dbcsr_copy(moments_op(imoment)%matrix, matrix_s(1)%matrix, "moments operator")
           call dbcsr_set(moments_op(imoment)%matrix, 0.0_dp)
           call cp_fm_create(moments_op_fm(imoment)%matrix, opstruct)
           if (do_rtp) then
              call cp_cfm_create(cmoments_op_fm(imoment)%matrix, opstruct)
           end if
        end do

        call get_qs_loc_env(qs_loc_env, local_molecules=local_molecules)
        nmolkinds = size(local_molecules%n_el)

        ! get the states list
        persistent_states => qs_env%persistent_states
        allocate(persistent_states(nsubsets, nspins), wcs_to_atoms(nmo, nspins))

        !! this might be garbage
        !do isubset = 1, nsubsets
        !   do ispin = 1, nspins
        !      ! get states belonging to this molecule
        !      if (associated(molecules(isubset)%lmi)) then
        !         nstates(1) = molecules(isubset)%lmi(ispin)%nstates
        !      else
        !         nstates(1) = 0
        !      end if
        !      nstates(2) = mepos
        !      call mp_maxloc(nstates, group)
        !      if (nstates(1) == 0) cycle
        !      ns = nstates(1)
        !      iproc = nstates(2)
        !      allocate(persistent_states(isubset, ispin)%states(ns))
        !      if (iproc == mepos) then
        !         persistent_states(isubset, ispin)%states(:) = molecules(isubset)%lmi(ispin)%states(:)
        !      else
        !         persistent_states(isubset, ispin)%states(:) = 0
        !      end if
        !      call mp_bcast(persistent_states(isubset, ispin)%states, iproc, group)
        !   end do
        !end do

        ! correct the atoms -> states lists, probably inefficient
        do ispin = 1, nspins
           do imo = 1, nmo
              r_min = huge(0.0_dp)
              jatom = 0
              do iatom = 1, natoms
                 if (norm2(pbc(particles(iatom)%r, wcs(ispin)%centres(:, imo), cell)) < r_min) then
                    jatom = iatom
                    !if (output_unit > 0) then
                    !   print *, imo, jatom, r_min, norm2(pbc(particles(iatom)%r, wcs(ispin)%centres(:, imo), cell))
                    !end if
                    r_min = norm2(pbc(particles(iatom)%r, wcs(ispin)%centres(:, imo), cell))
                 end if
              end do
              wcs_to_atoms(imo, ispin) = jatom
           end do
        end do

        !if (output_unit > 0) then
        !   print *, "wcs_to_atoms", wcs_to_atoms
        !   print *, "mol_to_nelectrons", mol_to_nelectrons
        !end if
        
        ! create the subsets -> states list
        do ispin = 1, nspins
           do isubset = 1, nsubsets
              allocate(persistent_states(isubset, ispin)%states(mol_to_nelectrons(isubset) / 2))
              persistent_states(isubset, ispin)%states = 0
              nwc = 1
              do imo = 1, nmo
                 if ((first_atoms(isubset) <= wcs_to_atoms(imo, ispin)) .and. &
                     (last_atoms(isubset) >= wcs_to_atoms(imo, ispin))) then
                    persistent_states(isubset, ispin)%states(nwc) = imo
                    nwc = nwc + 1
                 end if
              end do
           end do
        end do

        call get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
             ref_point=ref_point)

        do isubset = 1, nsubsets
           ! sometimes gives wrong results?
           call get_reference_point(rccs(:, isubset), qs_env=qs_env, reference=use_mom_ref_com, &
                ref_point=ref_point, ifirst=first_atoms(isubset), ilast=last_atoms(isubset))

           ! nuclear dipole
           do imoment = 1, nmoments
              lx = indco(1, imoment + 1)
              ly = indco(2, imoment + 1)
              lz = indco(3, imoment + 1)
              do iatom = first_atoms(isubset), last_atoms(isubset)
                 atomic_kind => particles(iatom)%atomic_kind
                 call get_atomic_kind(atomic_kind, kind_number=akind)
                 call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                 if (.not. ghost .and. .not. floating) then
                    call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                    r = particles(iatom)%r - rccs(:, isubset)
                    r = pbc(r, cell)
                    zmom = zeff
                    if (lx /= 0) zmom = zmom * r(1)**lx
                    if (ly /= 0) zmom = zmom * r(2)**ly
                    if (lz /= 0) zmom = zmom * r(3)**lz
                    moments(imoment, isubset) = moments(imoment, isubset) + zmom
                 end if
              end do
           end do

           ! Wannier function dipole
           do ispin = 1, nspins
              states => persistent_states(isubset, ispin)%states
              ns = size(states)
              !if (output_unit > 0) then
              !   print *, states
              !end if
              if (do_rtp) then
                 call cp_fm_to_cfm(mo_local(2*ispin-1)%matrix, mo_local(2*ispin)%matrix, ccoeffs)
                 call cp_cfm_set_all(cvec, (0.0_dp, 0.0_dp))
                 call cp_cfm_set_all(copvec, (0.0_dp, 0.0_dp))
                 call cp_cfm_set_all(cvecopvec, (0.0_dp, 0.0_dp))
                 do is = 1, ns
                    call cp_cfm_to_cfm(ccoeffs, cvec, 1, states(is), states(is))
                 end do

                 call build_local_moment_matrix(qs_env, moments_op, norders, &
                      ref_point=pbc(rccs(:, isubset), cell))
                 do imoment = 1, nmoments
                    call copy_dbcsr_to_fm(moments_op(imoment)%matrix, moments_op_fm(imoment)%matrix)
                    call cp_fm_to_cfm(moments_op_fm(imoment)%matrix, mtarget=cmoments_op_fm(imoment)%matrix)
                 end do

                 d = 0._dp
                 do imoment = 1, nmoments
                    call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), cmoments_op_fm(imoment)%matrix, &
                         cvec, (0._dp, 0._dp), copvec)
                    call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), cvec, copvec, (0._dp, 0._dp), cvecopvec)
                    do is = 1, ns
                       call cp_cfm_get_element(cvecopvec, states(is), states(is), clocd)
                       d(imoment) = d(imoment) + real(clocd)
                    end do
                    d(imoment) = cell%hmat(imoment, imoment) / twopi * &
                         aimag(log(exp(-gaussi * twopi / cell%hmat(imoment, imoment) * zwfc * d(imoment))))
                    moments(imoment, isubset) = moments(imoment, isubset) + d(imoment)
                 end do
              else
                 coeffs => mo_local(ispin)%matrix
                 
                 call cp_fm_set_all(vec, 0.0_dp)
                 call cp_fm_set_all(opvec, 0.0_dp)
                 call cp_fm_set_all(vecopvec, 0.0_dp)
                 do is = 1, ns
                    call cp_fm_to_fm(coeffs, vec, 1, states(is), states(is))
                 end do
                 
                 call build_local_moment_matrix(qs_env, moments_op, norders, &
                      ref_point=pbc(rccs(:, isubset), cell))
                 do imoment = 1, nmoments
                    call copy_dbcsr_to_fm(moments_op(imoment)%matrix, moments_op_fm(imoment)%matrix)
                 end do

                 d = 0.0_dp
                 do imoment = 1, nmoments
                    call cp_fm_gemm("N", "N", nao, nmo, nao, 1.0_dp, moments_op_fm(imoment)%matrix, &
                         vec, 0.0_dp, opvec)
                    call cp_fm_gemm("T", "N", nmo, nmo, nao, 1.0_dp, vec, opvec, 0.0_dp, vecopvec)
                    do is = 1, ns
                       call cp_fm_get_element(vecopvec, states(is), states(is), locd)
                       d(imoment) = d(imoment) + locd
                    end do
                    !d(imoment) = sum(vecopvec%local_data)
                    !call mp_sum(d(imoment), group)
                    d(imoment) = cell%hmat(imoment, imoment) / twopi * &
                         aimag(log(exp(-gaussi * twopi / cell%hmat(imoment, imoment) * zwfc * d(imoment))))
                    moments(imoment, isubset) = moments(imoment, isubset) + d(imoment)
                 end do
              end if
              deallocate(persistent_states(isubset, ispin)%states)
           end do
        end do
        deallocate(persistent_states)


        ! output
        output_unit = cp_print_key_unit_nr(logger, print_key, "WANNIER_MOMENTS", &
             extension=".WanMom", middle_name="WANNIER_MOMENTS")
        iter = cp_iter_string(logger%iter_info)
        if (output_unit > 0) then
           if (first_time) then
              write(output_unit, "(a)") "# iter_level, subset, moments [a.u.]"
           end if
           do isubset = 1, nsubsets
              name = molecules(isubset)%molecule_kind%name
              write(output_unit, "(a10,a10,39f12.6)") iter(1:10), name, moments(:, isubset)
           end do
        end if
        call cp_print_key_finished_output(output_unit, logger, print_key, "WANNIER_MOMENTS")

        ! cleanup
        do imoment = 1, nmoments
           call dbcsr_release(moments_op(imoment)%matrix)
           call cp_fm_release(moments_op_fm(imoment)%matrix)
           deallocate(moments_op(imoment)%matrix)
        end do
        deallocate(moments_op, moments_op_fm)
        deallocate(moments, first_atoms, last_atoms, rccs, ref_point, atom_to_mol, mol_to_nelectrons)

        call cp_fm_struct_release(opstruct)
        call cp_fm_struct_release(nmostruct)
        if (do_rtp) then
           do imoment = 1, nmoments
              call cp_cfm_release(cmoments_op_fm(imoment)%matrix)
           end do
           deallocate(cmoments_op_fm)
           call cp_cfm_release(cvec)
           call cp_cfm_release(copvec)
           call cp_cfm_release(cvecopvec)
           call cp_cfm_release(ccoeffs)
        else
           call cp_fm_release(vec)
           call cp_fm_release(opvec)
           call cp_fm_release(vecopvec)
        end if
     end if
  
     call timestop(handle)

   end subroutine calculate_wannier_moments


! **************************************************************************************************
!> \brief Calculates localised electric moments from Wannier centres
!> \param qs_env the QuickStep environment
!> \param qs_loc_env the localisation environment
!> \param print_key the print section
!> \param molecules the molecule list of the system
! **************************************************************************************************
   subroutine calculate_wc_moments(qs_env, qs_loc_env, print_key, molecules)

     type(qs_environment_type), pointer  :: qs_env
     type(qs_loc_env_new_type), pointer  :: qs_loc_env
     type(section_vals_type), pointer    :: print_key
     type(molecule_type), pointer        :: molecules(:)

     character(len=*), parameter  :: routineN = "calculate_wannier_moments", &
                                     routineP = moduleN//":"//routineN

     character(len=default_string_length)  :: iter, name
     integer                               :: akind, group, handle, i, iatom, imo, &
          imoment, jatom, &
          iproc, is, ispin, isubset, lx, ly, lz, &
          mepos, nwc, &
          nao, natoms, nmo, nmolkinds, nmoments, norders, ns, &
          nspins, nstates(2), nsubsets, output_unit
     integer, allocatable                   :: mol_to_nelectrons(:), first_atoms(:), last_atoms(:), &
          wcs_to_atoms(:, :), atom_to_mol(:)
     integer, pointer                       :: states(:)
     logical                                :: do_rtp, first_time, floating, ghost
     real(kind=dp)                          :: d(3), dlocal, r(3), zeff, zmom, zwfc, &
          m, mtot, rcc(3), r_min
     real(kind=dp), allocatable             :: moments(:, :), rccs(:, :)
     real(kind=dp), pointer                 :: ref_point(:)
     type(atomic_kind_type), pointer        :: atomic_kind
     type(cell_type), pointer               :: cell
     type(cp_fm_p_type), pointer            :: moments_op_fm(:)
     type(cp_fm_struct_type), pointer       :: fm_struct, mstruct, momstruct, opstruct
     type(cp_fm_type), pointer              :: coeffs, momvector, mvector, &
                                               omvector
     type(cp_logger_type), pointer          :: logger
     type(cp_para_env_type), pointer        :: para_env
     type(dbcsr_p_type), pointer            :: matrix_s(:), moments_op(:)
     type(dft_control_type), pointer        :: dft_control
     type(distribution_1d_type), pointer    :: local_molecules
     type(particle_type), pointer           :: particles(:)
     type(persistent_states_type), pointer  :: persistent_states(:, :)
     type(qs_kind_type), pointer            :: qs_kind_set(:)
     type(wannier_centres_type), pointer    :: wcs(:)

     call timeset(routineN, handle)

     logger => cp_get_default_logger()

     if (btest(cp_print_key_should_output(logger%iter_info, print_key, "WC_MOMENTS", &
                                          first_time=first_time), cp_p_file)) then

        ! user input
        call section_vals_val_get(print_key, "WC_MOMENTS%ORDER", i_val=norders)
        CPASSERT(norders >= 0)
        nmoments = ncoset(norders) - 1

        ! environment infos
        call get_qs_env(qs_env, cell=cell, dft_control=dft_control, &
             matrix_s=matrix_s, para_env=para_env, particle_set=particles, qs_kind_set=qs_kind_set, &
             wanniercentres=wcs)
        mepos = para_env%mepos
        group = para_env%group
        nspins = dft_control%nspins
        nsubsets = size(molecules)
        natoms = size(particles)
        zwfc = 3.0_dp - real(nspins, dp)  ! occupation of Wannier functions (uniform)
        allocate(moments(nmoments, nsubsets), first_atoms(nsubsets), &
             last_atoms(nsubsets), rccs(3, nsubsets), moments_op(nmoments), moments_op_fm(nmoments), &
             ref_point(3), mol_to_nelectrons(nsubsets), atom_to_mol(natoms))
        moments = 0.0_dp
        call get_molecule_set_info(molecules, &
             mol_to_first_atom=first_atoms, mol_to_last_atom=last_atoms, &
             mol_to_nelectrons=mol_to_nelectrons, atom_to_mol=atom_to_mol)
        nmo = sum(mol_to_nelectrons) / 2
        allocate(wcs_to_atoms(nmo, nspins))

        call get_qs_loc_env(qs_loc_env, local_molecules=local_molecules)
        nmolkinds = size(local_molecules%n_el)

        ! get the states list
        persistent_states => qs_env%persistent_states
        allocate(persistent_states(nsubsets, nspins))
        !do isubset = 1, nsubsets
        !   do ispin = 1, nspins
        !      ! get states belonging to this molecule
        !      if (associated(molecules(isubset)%lmi)) then
        !         nstates(1) = molecules(isubset)%lmi(ispin)%nstates
        !      else
        !         nstates(1) = 0
        !      end if
        !      nstates(2) = mepos
        !      call mp_maxloc(nstates, group)
        !      if (nstates(1) == 0) cycle
        !      ns = nstates(1)
        !      iproc = nstates(2)
        !      allocate(persistent_states(isubset, ispin)%states(ns))
        !      if (iproc == mepos) then
        !         persistent_states(isubset, ispin)%states(:) = molecules(isubset)%lmi(ispin)%states(:)
        !      else
        !         persistent_states(isubset, ispin)%states(:) = 0
        !      end if
        !      call mp_bcast(persistent_states(isubset, ispin)%states, iproc, group)
        !   end do
        !end do

        ! correct the atoms -> states lists, probably inefficient
        do ispin = 1, nspins
           do imo = 1, nmo
              r_min = huge(0.0_dp)
              jatom = 0
              do iatom = 1, natoms
                 if (norm2(pbc(particles(iatom)%r, wcs(ispin)%centres(:, imo), cell)) < r_min) then
                    jatom = iatom
                    r_min = norm2(pbc(particles(iatom)%r, wcs(ispin)%centres(:, imo), cell))
                 end if
              end do
              wcs_to_atoms(imo, ispin) = jatom
           end do
        end do

        ! create the subsets -> states list
        do ispin = 1, nspins
           do isubset = 1, nsubsets
              allocate(persistent_states(isubset, ispin)%states(mol_to_nelectrons(isubset) / 2))
              persistent_states(isubset, ispin)%states = 0
              nwc = 1
              do imo = 1, nmo
                 if ((first_atoms(isubset) <= wcs_to_atoms(imo, ispin)) .and. &
                     (last_atoms(isubset) >= wcs_to_atoms(imo, ispin))) then
                    persistent_states(isubset, ispin)%states(nwc) = imo
                    nwc = nwc + 1
                 end if
              end do
           end do
        end do

        call get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, &
             ref_point=ref_point)

        do isubset = 1, nsubsets
           ! sometimes gives wrong results?
           call get_reference_point(rccs(:, isubset), qs_env=qs_env, reference=use_mom_ref_com, &
                ref_point=ref_point, ifirst=first_atoms(isubset), ilast=last_atoms(isubset))

           ! nuclear dipole
           do imoment = 1, nmoments
              lx = indco(1, imoment + 1)
              ly = indco(2, imoment + 1)
              lz = indco(3, imoment + 1)
              do iatom = first_atoms(isubset), last_atoms(isubset)
                 atomic_kind => particles(iatom)%atomic_kind
                 call get_atomic_kind(atomic_kind, kind_number=akind)
                 call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                 if (.not. ghost .and. .not. floating) then
                    call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                    r = particles(iatom)%r - rccs(:, isubset)
                    r = pbc(r, cell)
                    zmom = zeff
                    if (lx /= 0) zmom = zmom * r(1)**lx
                    if (ly /= 0) zmom = zmom * r(2)**ly
                    if (lz /= 0) zmom = zmom * r(3)**lz
                    moments(imoment, isubset) = moments(imoment, isubset) + zmom
                 end if
              end do
           end do

           do ispin = 1, nspins

              states => persistent_states(isubset, ispin)%states
              ns = size(states)
              
              ! Wannier centre dipole
              do is = 1, ns
                 r = wcs(ispin)%centres(:, states(is)) - rccs(:, isubset)
                 r = pbc(r, cell)
                 do imoment = 1, nmoments
                    d(imoment) = cell%hmat(imoment, imoment) / twopi * &
                         aimag(log(exp(-gaussi * twopi / cell%hmat(imoment, imoment) * zwfc * r(imoment))))
                    moments(imoment, isubset) = moments(imoment, isubset) + d(imoment)
                 end do
              end do
              deallocate(persistent_states(isubset, ispin)%states)
           end do
        end do
        deallocate(persistent_states)

        ! output
        output_unit = cp_print_key_unit_nr(logger, print_key, "WC_MOMENTS", &
             extension=".WCMom", middle_name="WC_MOMENTS")
        iter = cp_iter_string(logger%iter_info)
        if (output_unit > 0) then
           if (first_time) then
              write(output_unit, "(a)") "# iter_level, subset, moments [a.u.]"
           end if
           do isubset = 1, nsubsets
              name = molecules(isubset)%molecule_kind%name
              write(output_unit, "(a10,a10,39f12.6)") iter(1:10), name, moments(:, isubset)
           end do
        end if
        call cp_print_key_finished_output(output_unit, logger, print_key, "WC_MOMENTS")

        ! cleanup
        deallocate(moments_op, moments_op_fm)
        deallocate(moments, first_atoms, last_atoms, rccs, ref_point, mol_to_nelectrons, atom_to_mol)
     end if
  
     call timestop(handle)

   end subroutine calculate_wc_moments

end module localized_moments
 

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
   use cp_control_types,       only: dft_control_type
   use cp_dbcsr_operations,    only: cp_dbcsr_sm_fm_multiply
   use cp_fm_basic_linalg,     only: cp_fm_schur_product
   use cp_fm_struct,           only: cp_fm_struct_create, &
                                     cp_fm_struct_release, &
                                     cp_fm_struct_type
   use cp_fm_types,            only: cp_fm_create, &
                                     cp_fm_get_info, &
                                     cp_fm_p_type, &
                                     cp_fm_release, &
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
                                     dbcsr_get_diag, &
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
     integer                               :: akind, first_atom, group, handle, iatom, &
          ilocmol, imolkind, imoment, &
          ispin, isubset, last_atom, mepos, natoms, &
          nlocmols, nmolkinds, nmoments, nspins, nsubsets, output_unit
     integer, allocatable                  :: first_atoms(:), first_basis(:), last_atoms(:), &
          last_basis(:), lx, ly, lz, nbasis(:)
     logical                               :: first_time, floating, ghost
     real(kind=dp)                         :: focc, m, mtot, r(3), zeff, zmom
     real(kind=dp), allocatable            :: diag(:), moments(:, :), rccs(:, :)
     type(atomic_kind_type), pointer       :: atomic_kind
     type(cell_type), pointer              :: cell
     type(cp_logger_type), pointer         :: logger
     type(cp_para_env_type), pointer       :: para_env
     type(dbcsr_p_type), pointer           :: densmat(:), matrix_s(:), moments_op(:)
     type(dbcsr_type), pointer             :: opmat
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
        focc = 3._dp - real(nspins, dp)
        call qs_rho_get(rho, rho_ao=densmat)
        allocate(moments(nmoments, nsubsets), first_atoms(nsubsets), first_basis(nsubsets), &
             last_atoms(nsubsets), last_basis(nsubsets), &
             rccs(3, nsubsets), moments_op(nmoments))
        moments = 0._dp
        call get_molecule_set_info(molecules, &
             mol_to_first_atom=first_atoms, mol_to_last_atom=last_atoms, mol_to_nbasis=nbasis)
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
        call dbcsr_copy(opmat, matrix_s(1)%matrix, "operator * densmat")
        nmolkinds = size(local_molecules%n_el)

        do isubset = 1, nsubsets
           ! gives wrong results for some reason
           !call get_reference_point(rccs(:, isubset), qs_env=qs_env, reference=use_mom_ref_com, &
           !     ref_point=ref_point, ifirst=1, ilast=3)!, &
           !     !ifirst=first_atoms(isubset), ilast=last_atoms(isubset))

           ! custom COM
           rccs(:, isubset) = 0._dp
           mtot = 0._dp
           do iatom = first_atoms(isubset), last_atoms(isubset)
              r = particles(iatom)%r
              r = pbc(r, cell)
              atomic_kind => particles(iatom)%atomic_kind
              call get_atomic_kind(atomic_kind, mass=m)
              rccs(:, isubset) = rccs(:, isubset) + m * r
              mtot = mtot + m
           end do
           if (abs(mtot) > 0._dp) then
              rccs(:, isubset) = rccs(:, isubset) / mtot
           end if
           rccs(:, isubset) = pbc(rccs(:, isubset), cell)

           call build_local_moment_matrix(qs_env, moments_op, norders, rccs(:, isubset))

           do ispin = 1, nspins
              do imoment = 1, nmoments
                 call dbcsr_multiply("N", "N", 1._dp, densmat(ispin)%matrix, &
                      moments_op(imoment)%matrix, 0._dp, opmat)
                 call dbcsr_get_diag(opmat, diag)
                 moments(imoment, isubset) = moments(imoment, isubset) + &
                      sum(diag(first_basis(isubset):last_basis(isubset)))
              end do
           end do
        end do

        ! nuclear moments
        do imolkind = 1, nmolkinds
           nlocmols = size(local_molecules%list(imolkind)%array)
           do ilocmol = 1, nlocmols
              isubset = local_molecules%list(imolkind)%array(ilocmol)  ! global index
              molkind => molecules(isubset)%molecule_kind
              first_atom = first_atoms(isubset)
              last_atom = last_atoms(isubset)
              natoms = last_atom - first_atom + 1
              
              do imoment = 1, nmoments
                 lx = indco(1, imoment + 1)
                 ly = indco(2, imoment + 1)
                 lz = indco(3, imoment + 1)
                 do iatom = first_atom, last_atom
                    atomic_kind => particles(iatom)%atomic_kind
                    call get_atomic_kind(atomic_kind, kind_number=akind)
                    call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                    if (.not. ghost .and. .not. floating) then
                       call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                       r = particles(iatom)%r - rccs(:, isubset)
                       r = pbc(r, cell)
                       zmom = zeff
                       if (lx /= 0) zmom = zmom*r(1)**lx
                       if (ly /= 0) zmom = zmom*r(2)**ly
                       if (lz /= 0) zmom = zmom*r(3)**lz
                       moments(imoment, isubset) = moments(imoment, isubset) + zmom
                    end if
                 end do
              end do
           end do
        end do
        call mp_sum(moments, group)

        do isubset = 1, nsubsets
           do imoment = 1, nmoments
              moments(imoment, isubset) = -cell%hmat(imoment, imoment) / twopi * &
                   aimag(log(exp(-gaussi * twopi * cell%h_inv(imoment, imoment) * moments(imoment, isubset))))
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
     integer                               :: akind, group, handle, i, iatom, &
          imoment, &
          iproc, ispin, isubset, lx, ly, lz, &
          mepos, &
          nao, nmo, nmolkinds, nmoments, norders, ns, &
          nspins, nstates(2), nsubsets, output_unit
     integer, allocatable                  :: first_atoms(:), last_atoms(:), states(:)
     logical                               :: first_time, floating, ghost
     real(kind=dp)                         :: m, mtot, r(3), zeff, zmom, zwfc
     real(kind=dp), allocatable            :: moments(:, :), rccs(:, :)
     type(atomic_kind_type), pointer       :: atomic_kind
     type(cell_type), pointer              :: cell
     type(cp_fm_struct_type), pointer      :: fm_struct
     type(cp_fm_type), pointer             :: coeffs, momvector, mvector, omvector
     type(cp_logger_type), pointer         :: logger
     type(cp_para_env_type), pointer       :: para_env
     type(dbcsr_p_type), pointer           :: matrix_s(:), moments_op(:)
     type(dft_control_type), pointer       :: dft_control
     type(distribution_1d_type), pointer   :: local_molecules
     type(particle_type), pointer          :: particles(:)
     type(qs_kind_type), pointer           :: qs_kind_set(:)
     type(wannier_centres_type), pointer   :: wcs(:)

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
        nsubsets = size(molecules)
        zwfc = 3.0_dp - real(nspins, dp)  ! occupation of Wannier functions (uniform)
        allocate(moments(nmoments, nsubsets), first_atoms(nsubsets), &
             last_atoms(nsubsets), rccs(3, nsubsets), moments_op(nmoments))
        moments = 0.0_dp
        call get_molecule_set_info(molecules, &
             mol_to_first_atom=first_atoms, mol_to_last_atom=last_atoms)
        do imoment = 1, nmoments
           allocate(moments_op(imoment)%matrix)
           call dbcsr_copy(moments_op(imoment)%matrix, matrix_s(1)%matrix, "moments operator")
           call dbcsr_set(moments_op(imoment)%matrix, 0.0_dp)
        end do

        call get_qs_loc_env(qs_loc_env, local_molecules=local_molecules)
        nmolkinds = size(local_molecules%n_el)

        !! nuclear moments
        !do imolkind = 1, nmolkinds
        !   nlocmols = size(local_molecules%list(imolkind)%array)
        !   do ilocmol = 1, nlocmols
        !      isubset = local_molecules%list(imolkind)%array(ilocmol)  ! global index
        !      molkind => molecules(isubset)%molecule_kind
        !      first_atom = first_atoms(isubset)
        !      last_atom = last_atoms(isubset)
        !      natoms = last_atom - first_atom + 1
        !      
        !      do imoment = 1, nmoments
        !         lx = indco(1, imoment + 1)
        !         ly = indco(2, imoment + 1)
        !         lz = indco(3, imoment + 1)
        !         do iatom = first_atom, last_atom
        !            atomic_kind => particles(iatom)%atomic_kind
        !            call get_atomic_kind(atomic_kind, kind_number=akind)
        !            call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
        !            if (.not. ghost .and. .not. floating) then
        !               call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
        !               r = particles(iatom)%r - rccs(:, isubset)
        !               r = pbc(r, cell)
        !               zmom = zeff
        !               if (lx /= 0) zmom = zmom*r(1)**lx
        !               if (ly /= 0) zmom = zmom*r(2)**ly
        !               if (lz /= 0) zmom = zmom*r(3)**lz
        !               moments(imoment, isubset) = moments(imoment, isubset) + zmom
        !            end if
        !         end do
        !      end do
        !   end do
        !end do
        !call mp_sum(moments, group)

        ! electronic moments
        do isubset = 1, nsubsets
           ! gives wrong results for some reason
           !call get_reference_point(rccs(:, isubset), qs_env=qs_env, reference=use_mom_ref_com, &
           !     ref_point=ref_point, ifirst=1, ilast=3)!, &
           !     !ifirst=first_atoms(isubset), ilast=last_atoms(isubset))

           ! custom COM
           rccs(:, isubset) = 0._dp
           mtot = 0._dp
           do iatom = first_atoms(isubset), last_atoms(isubset)
              r = particles(iatom)%r
              r = pbc(r, cell)
              atomic_kind => particles(iatom)%atomic_kind
              call get_atomic_kind(atomic_kind, mass=m)
              rccs(:, isubset) = rccs(:, isubset) + m * r
              mtot = mtot + m
           end do
           if (abs(mtot) > 0._dp) then
              rccs(:, isubset) = rccs(:, isubset) / mtot
           end if
           rccs(:, isubset) = pbc(rccs(:, isubset), cell)

           do imoment = 1, nmoments
              call dbcsr_set(moments_op(imoment)%matrix, 0.0_dp)
           end do

           call build_local_moment_matrix(qs_env, moments_op, norders, rccs(:, isubset))

           ! nuclear moments
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
                    if (lx /= 0) zmom = zmom*r(1)**lx
                    if (ly /= 0) zmom = zmom*r(2)**ly
                    if (lz /= 0) zmom = zmom*r(3)**lz
                    moments(imoment, isubset) = moments(imoment, isubset) + zmom
                 end if
              end do
           end do

           do ispin = 1, nspins
              ! get states belonging to this molecule
              if (associated(molecules(isubset)%lmi)) then
                 nstates(1) = molecules(isubset)%lmi(ispin)%nstates
              else
                 nstates(1) = 0
              end if
              nstates(2) = mepos
              call mp_maxloc(nstates, group)
              if (nstates(1) == 0) cycle
              ns = nstates(1)
              iproc = nstates(2)
              allocate(states(ns))
              if (iproc == mepos) then
                 states(:) = molecules(isubset)%lmi(ispin)%states(:)
              else
                 states(:) = 0
              end if
              call mp_bcast(states, iproc, group)

              coeffs => mo_local(ispin)%matrix
              call cp_fm_get_info(coeffs, ncol_global=nmo, nrow_global=nao)
              call cp_fm_struct_create(fm_struct, ncol_global=ns, nrow_global=nao, &
                   para_env=coeffs%matrix_struct%para_env, context=coeffs%matrix_struct%context)
              call cp_fm_create(mvector, fm_struct, name="mvector")
              call cp_fm_create(omvector, fm_struct, name="omvector")
              call cp_fm_create(momvector, fm_struct, name="momvector")
              call cp_fm_struct_release(fm_struct)

              do i = 1, ns
                 call cp_fm_to_fm(coeffs, mvector, 1, states(i), i)
              end do

              do imoment = 1, nmoments
                 call cp_dbcsr_sm_fm_multiply(moments_op(imoment)%matrix, mvector, omvector, ns)
                 call cp_fm_schur_product(mvector, omvector, momvector)
                 moments(imoment, isubset) = moments(imoment, isubset) - zwfc * sum(momvector%local_data)
              end do
              
              call cp_fm_release(mvector)
              call cp_fm_release(omvector)
              call cp_fm_release(momvector)
              deallocate(states)
           end do
        end do

        do isubset = 1, nsubsets
           do imoment = 1, nmoments
              moments(imoment, isubset) = -cell%hmat(imoment, imoment) / twopi * &
                   aimag(log(exp(-gaussi * twopi * cell%h_inv(imoment, imoment) * moments(imoment, isubset))))
           end do
        end do

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
           deallocate(moments_op(imoment)%matrix)
        end do
        deallocate(moments_op)
        deallocate(moments, first_atoms, last_atoms, rccs)

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

     character(len=*), parameter  :: routineN = "calculate_wc_moments", &
                                     routineP = moduleN//":"//routineN

     character(len=default_string_length)  :: iter, name
     integer                               :: akind, first_atom, &
          group, handle, iatom, ilocmol, imolkind, imoment, iproc, ispin, istate, isubset, &
          last_atom, lx, ly, lz, &
          mepos, natoms, nlocmols, &
          nmolkinds, nmoments, &
          norders, ns, nspins, nstates(2), nsubsets, output_unit
     integer, allocatable                  :: first_atoms(:), last_atoms(:), states(:)
     logical                               :: first_time, floating, ghost
     real(kind=dp)                         :: m, mtot, r(3), zeff, zmom, zwfc
     real(kind=dp), allocatable            :: moments(:, :), rccs(:, :)
     type(atomic_kind_type), pointer       :: atomic_kind
     type(cell_type), pointer              :: cell
     type(cp_logger_type), pointer         :: logger
     type(cp_para_env_type), pointer       :: para_env
     type(dft_control_type), pointer       :: dft_control
     type(distribution_1d_type), pointer   :: local_molecules
     type(molecule_kind_type), pointer     :: molkind
     type(particle_type), pointer          :: particles(:)
     type(qs_kind_type), pointer           :: qs_kind_set(:)
     type(wannier_centres_type), pointer   :: wcs(:)

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
             para_env=para_env, particle_set=particles, qs_kind_set=qs_kind_set, &
             wanniercentres=wcs)
        mepos = para_env%mepos
        group = para_env%group
        nspins = dft_control%nspins
        nsubsets = size(molecules)
        zwfc = 3._dp - real(nspins, dp)
        allocate(moments(nmoments, nsubsets), first_atoms(nsubsets), &
             last_atoms(nsubsets), rccs(3, nsubsets))
        moments = 0._dp
        call get_molecule_set_info(molecules, &
             mol_to_first_atom=first_atoms, mol_to_last_atom=last_atoms)

        call get_qs_loc_env(qs_loc_env, local_molecules=local_molecules)
        nmolkinds = size(local_molecules%n_el)

        ! nuclear moments
        do imolkind = 1, nmolkinds
           nlocmols = size(local_molecules%list(imolkind)%array)
           do ilocmol = 1, nlocmols
              isubset = local_molecules%list(imolkind)%array(ilocmol)  ! global index
              molkind => molecules(isubset)%molecule_kind
              first_atom = first_atoms(isubset)
              last_atom = last_atoms(isubset)
              natoms = last_atom - first_atom + 1
              
              do imoment = 1, nmoments
                 lx = indco(1, imoment + 1)
                 ly = indco(2, imoment + 1)
                 lz = indco(3, imoment + 1)
                 do iatom = first_atom, last_atom
                    atomic_kind => particles(iatom)%atomic_kind
                    call get_atomic_kind(atomic_kind, kind_number=akind)
                    call get_qs_kind(qs_kind_set(akind), ghost=ghost, floating=floating)
                    if (.not. ghost .and. .not. floating) then
                       call get_qs_kind(qs_kind_set(akind), core_charge=zeff)
                       r = particles(iatom)%r - rccs(:, isubset)
                       r = pbc(r, cell)
                       zmom = zeff
                       if (lx /= 0) zmom = zmom*r(1)**lx
                       if (ly /= 0) zmom = zmom*r(2)**ly
                       if (lz /= 0) zmom = zmom*r(3)**lz
                       moments(imoment, isubset) = moments(imoment, isubset) + zmom
                    end if
                 end do
              end do
           end do
        end do
        call mp_sum(moments, group)
        
        ! wannier centre moments
        do isubset = 1, nsubsets
           ! gives wrong results for some reason
           !call get_reference_point(rccs(:, isubset), qs_env=qs_env, reference=use_mom_ref_com, &
           !     ref_point=ref_point, ifirst=1, ilast=3)!, &
           !     !ifirst=first_atoms(isubset), ilast=last_atoms(isubset))

           ! custom COM
           rccs(:, isubset) = 0._dp
           mtot = 0._dp
           do iatom = first_atoms(isubset), last_atoms(isubset)
              r = particles(iatom)%r
              r = pbc(r, cell)
              !print *, "r", r
              atomic_kind => particles(iatom)%atomic_kind
              call get_atomic_kind(atomic_kind, mass=m)
              rccs(:, isubset) = rccs(:, isubset) + m * r
              mtot = mtot + m
           end do
           if (abs(mtot) > 0._dp) then
              rccs(:, isubset) = rccs(:, isubset) / mtot
           end if
           
           rccs(:, isubset) = pbc(rccs(:, isubset), cell)
           !print *, isubset, "rcc", rccs(:, isubset)

           do ispin = 1, nspins
              ! get states belonging to this molecule
              if (associated(molecules(isubset)%lmi)) then
                 nstates(1) = molecules(isubset)%lmi(ispin)%nstates
              else
                 nstates(1) = 0
              end if
              nstates(2) = mepos
              call mp_maxloc(nstates, group)
              if (nstates(1) == 0) cycle
              ns = nstates(1)
              iproc = nstates(2)
              allocate(states(ns))
              if (iproc == mepos) then
                 states(:) = molecules(isubset)%lmi(ispin)%states(:)
              else
                 states(:) = 0
              end if
              call mp_bcast(states, iproc, group)

              do istate = 1, ns
                 do imoment = 1, nmoments
                    lx = indco(1, imoment + 1)
                    ly = indco(2, imoment + 1)
                    lz = indco(3, imoment + 1)
                    !print *, "states", states
                    !print *, size(wcs)
                    !print *, size(wcs(1)%centres, 1)
                    !print *, size(wcs(1)%centres, 2)
                    r = wcs(ispin)%centres(:, states(istate)) - rccs(:, isubset)
                    r = pbc(r, cell)
                    zmom = zwfc
                    if (lx /= 0) zmom = zmom * r(1)**lx
                    if (ly /= 0) zmom = zmom * r(2)**ly
                    if (lz /= 0) zmom = zmom * r(3)**lz
                    moments(imoment, isubset) = moments(imoment, isubset) - zmom
                 end do
              end do
              deallocate(states)
           end do
        end do

        do isubset = 1, nsubsets
           do imoment = 1, nmoments
              moments(imoment, isubset) = -cell%hmat(imoment, imoment) / twopi * &
                   aimag(log(exp(-gaussi * twopi * cell%h_inv(imoment, imoment) * moments(imoment, isubset))))
           end do
        end do

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
        
     end if

     call timestop(handle)

   end subroutine calculate_wc_moments

end module localized_moments

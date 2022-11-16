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
   use basis_set_types,        only: gto_basis_set_p_type, &
                                     gto_basis_set_type
   use cell_types,             only: cell_type, &
                                     pbc
   use cp_cfm_basic_linalg,    only: cp_cfm_gemm, &
                                     cp_cfm_lu_decompose, &
                                     cp_cfm_schur_product
   use cp_cfm_types,           only: cp_cfm_create, &
                                     cp_cfm_get_element, &
                                     cp_cfm_get_info, &
                                     cp_cfm_p_type, &
                                     cp_cfm_release, &
                                     cp_cfm_set_all, &
                                     cp_cfm_to_cfm, &
                                     cp_cfm_type, &
                                     cp_cfm_write_formatted, &
                                     cp_fm_to_cfm
   use cp_control_types,       only: dft_control_type
   use cp_dbcsr_cp2k_link,     only: cp_dbcsr_alloc_block_from_nbl
   use cp_dbcsr_operations,    only: copy_dbcsr_to_cfm, &
                                     copy_dbcsr_to_fm, &
                                     cp_dbcsr_sm_fm_multiply, &
                                     dbcsr_allocate_matrix_set, &
                                     dbcsr_deallocate_matrix_set
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
                                     dbcsr_create, &
                                     dbcsr_desymmetrize, &
                                     dbcsr_dot, &
                                     dbcsr_get_block_p, &
                                     dbcsr_get_diag, &
                                     dbcsr_get_info, &
                                     dbcsr_init_p, &
                                     dbcsr_iterator_blocks_left, &
                                     dbcsr_iterator_next_block, &
                                     dbcsr_iterator_start, &
                                     dbcsr_iterator_stop, &
                                     dbcsr_iterator_type, &
                                     dbcsr_multiply, &
                                     dbcsr_p_type, &
                                     dbcsr_print, &
                                     dbcsr_release, &
                                     dbcsr_set, &
                                     dbcsr_trace, &
                                     dbcsr_type, &
                                     dbcsr_type_antisymmetric, &
                                     dbcsr_type_no_symmetry
   use distribution_1d_types,  only: distribution_1d_type
   use input_constants,        only: use_mom_ref_coac, &
                                     use_mom_ref_com, &
                                     use_mom_ref_zero
   use input_section_types,    only: section_vals_type, &
                                     section_vals_val_get
   use kinds,                  only: default_string_length, &
                                     dp
   use mathconstants,          only: gaussi, &
                                     pi, &
                                     twopi
   use message_passing,        only: mp_bcast, &
                                     mp_maxloc, &
                                     mp_sum
   use molecule_kind_types,    only: get_molecule_kind, &
                                     molecule_kind_type
   use molecule_types,         only: get_molecule_set_info, &
                                     molecule_type
   use moments_utils,          only: get_reference_point
   use orbital_pointers,       only: indco, &
                                     ncoset
   use particle_types,         only: particle_type
   use physcon,                only: c_light_au
   use qs_environment_types,   only: get_qs_env, &
                                     persistent_states_type, &
                                     qs_environment_type
   use qs_kind_types,          only: get_qs_kind, &
                                     get_qs_kind_set, &
                                     qs_kind_type
   use qs_loc_types,           only: get_qs_loc_env, &
                                     qs_loc_env_new_type
   use qs_mo_types,            only: mo_set_p_type
   use qs_moments,             only: build_berry_moment_matrix, &
                                     build_local_magmom_matrix, &
                                     build_local_moment_matrix, &
                                     build_local_moments_der_matrix, &
                                     calculate_commutator_nl_terms, &
                                     op_orbbas, &
                                     op_orbbas_rtp, &
                                     print_moments, &
                                     print_moments_nl, &
                                     set_label
   use qs_neighbor_list_types, only: neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_set_p_type
   use qs_operators_ao,        only: build_lin_mom_matrix
   use qs_rho_types,           only: qs_rho_get, &
                                     qs_rho_type
   use wannier_states_types,   only: wannier_centres_type

#include "./base/base_uses.f90"

   implicit none

   private

! *** globals ***

   character(len=*), parameter, private :: &
        moduleN = "localized_moments", &
        header = "# iteration, iteration_level, " // &
        "nucl. charge, el. charge, total charge, " // &
        "eldip x, eldip y, eldip z, " // &
        "elqup xx, elqup xy, elqup xz, elqup yy, elqup yz, elqup zz, " // &
        "mdip x, mdip y, mdip z, " // &
        "mom x, mom y, mom z, " // &
        "elqupvel xx, elqupvel xy, elqupvel xz, elqupvel yx, elqupvel yy, elqupvel zz, " // &
        "mdipnl x, mdipnl y, mdipnl z, momnl x, momnl y, momnl z, " // &
        "elqupvelnl xx, elqupvelnl xy, elqupvelnl xz, elqupvelnl yx, elqupvelnl yy, elqupvelnl zz, " // &
        "ref x, ref y, refz"
   integer, private            :: iteration = 0
   integer, parameter, private :: nm = 9, nmm = 3
   type(persistent_states_type), pointer, private :: persistent_states(:, :)

! *** public subroutines ***

   public :: calculate_localized_moments

! *** body ***

contains

! **************************************************************************************************
!> \brief Calculates localised electric and magnetic moments from Wannier functions or Kim-Gordon
!>        subdensities
! **************************************************************************************************
   subroutine calculate_localized_moments(qs_env, mo_local)

      type(qs_environment_type), pointer     :: qs_env
      type(cp_fm_p_type), pointer, optional  :: mo_local(:)

      character(len=*), parameter  :: routineN = "calculate_localized_moments", &
                                      routineP = moduleN//":"//routineN

      character(len=default_string_length)   :: iter, name                    
      integer                                :: handle, iatom, iproc, jatom, imo, isubset, itimes, &
                                                natoms, &
                                                nao, nmo, nmolkinds, nmoments, &
                                                norders, ns, &
                                                nspins, nstates(2), nsubsets, nwcs, output_unit
      integer, allocatable                   :: first_atoms(:), first_basis(:), last_atoms(:), &
                                                last_basis(:), nbasis(:), &
                                                nelectrons(:), wcs_to_atoms( :)
      integer, pointer                       :: states(:)
      logical                                :: do_kg, do_perd, do_periodic
      real(dp)                               :: dist, r_min, rcc(3), rcc0(3), ra(3), rb(3), rpoint(3), r(3), &
           mtot, mass_low, center(3), mass, qupprefactor
      real(dp), allocatable                  :: moments(:), rccs(:, :)
      real(dp), pointer                      :: ref_point(:)
      type(atomic_kind_type), pointer        :: akind
      type(cell_type), pointer               :: cell
      type(cp_logger_type), pointer          :: logger
      type(cp_para_env_type), pointer        :: para_env
      type(dbcsr_p_type), pointer            :: densmat(:), matrix_s(:), dbop(:)
      type(dft_control_type), pointer        :: dft_control
      type(molecule_type), pointer           :: molecules(:)
      type(particle_type), pointer           :: particles(:)
      !type(persistent_states_type), pointer  :: persistent_states(:, :)
      type(qs_kind_type), pointer            :: qs_kinds(:)
      type(qs_rho_type), pointer             :: rho
      type(section_vals_type), pointer       :: input
      type(wannier_centres_type), pointer    :: wcs(:)

      call timeset(routineN, handle)

      logger => cp_get_default_logger()
      iteration = iteration + 1

      ! user input
      do_periodic = .false.
      norders = 2

      ! environment infos
      call get_qs_env(qs_env, &
           cell=cell, &
           dft_control=dft_control, &
           input=input, &
           matrix_s=matrix_s, &
           molecule_set=molecules, &
           para_env=para_env, &
           particle_set=particles, &
           qs_kind_set=qs_kinds, &
           rho=rho)
      ! charge(3), eldip(3), elqup(6), mdip(3), mom(3), elqupvel(6), mdipnl(3), momnl(3), elqupvelnl(6)
      nmoments = 36
      natoms = size(particles)
      nsubsets = size(molecules)
      allocate(first_atoms(nsubsets), first_basis(nsubsets), last_atoms(nsubsets), &
           last_basis(nsubsets), nbasis(nsubsets), nelectrons(nsubsets), moments(nmoments), &
           rccs(nsubsets, 3), ref_point(3))
      call get_molecule_set_info(molecules, &
           mol_to_first_atom=first_atoms, &
           mol_to_last_atom=last_atoms, &
           mol_to_nbasis=nbasis, &
           mol_to_nelectrons=nelectrons)

      ! reference points
      call get_reference_point(rcc0, qs_env=qs_env, reference=use_mom_ref_zero, ifirst=1, &
           ilast=natoms)
      call get_reference_point(rcc, qs_env=qs_env, reference=use_mom_ref_com, ifirst=1, &
           ilast=natoms)
      do isubset = 1, nsubsets
         call get_reference_point(rccs(isubset, :), &
              qs_env=qs_env, &
              reference=use_mom_ref_com, &
              ifirst=first_atoms(isubset), &
              ilast=last_atoms(isubset))
      end do

      if (iteration == 1) then
         call cp_fm_get_info(mo_local(1)%matrix, ncol_global=nmo, nrow_global=nao)
         call get_qs_env(qs_env, wanniercentres=wcs)
         allocate(persistent_states(nsubsets, 1), wcs_to_atoms(nmo))
         ! states -> subset
         do imo = 1, nmo
            r_min = huge(0._dp)
            jatom = 0
            do iatom = 1, natoms
               dist = norm2(pbc(wcs(1)%centres(:, imo) - particles(iatom)%r, cell))
               if (dist < r_min) then
                  jatom = iatom
                  r_min = dist
               end if
            end do
            wcs_to_atoms(imo) = jatom
         end do
         
         ! create the subsets -> states list
         nelectrons = 0
         do isubset = 1, nsubsets
            do imo = 1, nmo
               if ((first_atoms(isubset) <= wcs_to_atoms(imo)) .and. &
                    (last_atoms(isubset) >= wcs_to_atoms(imo))) then
                  nelectrons(isubset) = nelectrons(isubset) + 1
               end if
            end do
         end do
         ! subset -> states
         do isubset = 1, nsubsets
            allocate(persistent_states(isubset, 1)%states(nelectrons(isubset)))
            nwcs = 1
            do imo = 1, nmo
               if ((first_atoms(isubset) <= wcs_to_atoms(imo)) .and. &
                    (last_atoms(isubset) >= wcs_to_atoms(imo))) then
                  persistent_states(isubset, 1)%states(nwcs) = imo
                  nwcs = nwcs + 1
               end if
            end do
            print *, "dbg| states(isubset)", isubset, size(persistent_states(isubset, 1)%states), &
                 persistent_states(isubset, 1)%states
         end do
         deallocate(wcs_to_atoms)
      end if

      moments = 0._dp
      call calculate_loc_moments(moments, qs_env, rcc, mo_local)
      call print_localized_moments(logger, input, moments, rcc, "ref")
      moments = 0._dp
      call calculate_loc_moments(moments, qs_env, rcc0, mo_local)
      call print_localized_moments(logger, input, moments, rcc0, "ref0")

      ! subsystems
      do isubset = 1, nsubsets
         states => persistent_states(isubset, 1)%states
         name = molecules(isubset)%molecule_kind%name

         ! sub
         moments = 0._dp
         call calculate_loc_moments(moments, qs_env, rccs(isubset, :), mo_local, states)
         call print_localized_moments(logger, input, moments, rccs(isubset, :), name)
      end do
      
      deallocate(first_atoms, first_basis, last_atoms, last_basis, moments, nbasis, nelectrons, &
           rccs, ref_point)

      call timestop(handle)
     
   end subroutine calculate_localized_moments

! **************************************************************************************************
   subroutine calculate_loc_moments(moments, qs_env, rcc, mo_local, statelist)

     real(dp)                            :: moments(:)
     type(qs_environment_type), pointer  :: qs_env
     real(dp)                            :: rcc(3)
     type(cp_fm_p_type), pointer         :: mo_local(:)
     integer, pointer, optional          :: statelist(:)

     call cp_fm_get_info(mo_local(1)%matrix, nrow_global=nao, ncol_global=nmo)
     zwfc = 2._dp

     if (present(statelist)) then
        ns = size(statelist)
        do i = 1, ns
           states(i) = statelist(i)
        end do
     else
        ns = nmo
        do i = 1, ns
           states(i) = i
        end do
     end if

     ! charge

     ! mu

     ! muv

     ! m

     ! theta

     ! thetav

   end subroutine calculate_loc_moments


! **************************************************************************************************
! **************************************************************************************************
   subroutine build_wheeler_moment_matrix(qs_env, o, rcc, m)

      type(qs_environment_type), pointer  :: qs_env
      type(cp_cfm_p_type), pointer        :: o(:)
      real(dp)  :: rcc(3)
      type(cp_cfm_type), pointer  :: m

      integer  :: i, j, k, nao, nmo
      real(dp)  :: kvec(3)
      type(cp_fm_struct_type), pointer  :: o_struct
      type(cp_fm_type), pointer  :: cmat, smat
      type(dbcsr_p_type), pointer  :: matrix_s(:), qo(:)
      type(dbcsr_type), pointer  :: cosmat, sinmat
      type(cell_type), pointer  :: cell

      call get_qs_env(qs_env, cell=cell, matrix_s=matrix_s)

      call cp_cfm_get_info(m, nrow_global=nao)
      call cp_fm_struct_create(o_struct, ncol_global=nao, template_fmstruct=m%matrix_struct)
      call cp_fm_create(cmat, o_struct)
      call cp_fm_create(smat, o_struct)
      nullify(cosmat, sinmat)
      allocate(cosmat, sinmat)
      call dbcsr_copy(cosmat, matrix_s(1)%matrix)
      call dbcsr_copy(sinmat, matrix_s(1)%matrix)
      allocate(o(9))
      do i = 1, 9
         call cp_cfm_create(o(i)%matrix, o_struct)
      end do
      
      !do i = 1, 3
      !   kvec = twopi * cell%h_inv(i, :)
      !   call build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec)
      !   call copy_dbcsr_to_fm(cosmat, cmat)
      !   call copy_dbcsr_to_fm(sinmat, smat)
      !   smat%local_data = -smat%local_data
      !   call cp_fm_to_cfm(cmat, smat, o(i)%matrix)
      !end do
      !k = 4
      !do i = 1, 3
      !   do j = i, 3
      !      kvec = 0.0_dp
      !      kvec(i) = twopi * cell%h_inv(i, i) * cell%h_inv(j, j)
      !      !print *, "dbg| kvec", kvec
      !      call build_berry_moment_matrix(qs_env, cosmat, sinmat, kvec, q=.true.)
      !      call copy_dbcsr_to_fm(cosmat, cmat)
      !      call copy_dbcsr_to_fm(sinmat, smat)
      !      smat%local_data = -smat%local_data
      !      call cp_fm_to_cfm(cmat, smat, o(k)%matrix)
      !      k = k + 1
      !   end do
      !end do
      nullify(qo)
      call dbcsr_allocate_matrix_set(qo, 9)
      do i = 1, 9
         allocate(qo(i)%matrix)
         call dbcsr_copy(qo(i)%matrix, matrix_s(1)%matrix, "qo")
      end do
      call build_local_moment_matrix(qs_env, qo, 2, ref_point=rcc)
      do i = 4, 9
         call copy_dbcsr_to_fm(qo(i)%matrix, cmat)
         call cp_fm_to_cfm(cmat, mtarget=o(i)%matrix)
      end do
      !print *, "dbg| o(4)", o(4)%matrix%local_data
      !qo => o(4:9)
      !call build_wheeler_quadrupole_matrix(qs_env, qo)
      call cp_fm_release(cmat)
      call cp_fm_release(smat)
      call dbcsr_release(cosmat)
      call dbcsr_release(sinmat)
      call cp_fm_struct_release(o_struct)
      deallocate(cosmat, sinmat)
      
   end subroutine build_wheeler_moment_matrix


! **************************************************************************************************
! **************************************************************************************************
   subroutine build_modernmag_matrix(qs_env, o, rcc, m)

      type(qs_environment_type), pointer  :: qs_env
      type(cp_cfm_p_type), pointer        :: o(:)
      real(dp)                            :: rcc(3)
      type(cp_cfm_type), pointer  :: m

      integer  :: i, nao
      type(cp_fm_struct_type), pointer  :: o_struct
      type(cp_fm_type), pointer  :: or
      type(dbcsr_p_type), pointer  :: op(:), matrix_s(:)
      type(neighbor_list_set_p_type), pointer  :: sab_orb(:)
      
      call get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
      call cp_cfm_get_info(m, nrow_global=nao)
      call cp_fm_struct_create(o_struct, ncol_global=nao, template_fmstruct=m%matrix_struct)
      call cp_fm_create(or, o_struct)
      
      nullify(op)
      call dbcsr_allocate_matrix_set(op, 3)
      do i = 1, 3
         call dbcsr_init_p(op(i)%matrix)
         call dbcsr_create(op(i)%matrix, template=matrix_s(1)%matrix, &
              matrix_type=dbcsr_type_antisymmetric)
         call cp_dbcsr_alloc_block_from_nbl(op(i)%matrix, sab_orb)
         call dbcsr_set(op(i)%matrix, 0._dp)
      end do
      call build_local_magmom_matrix(qs_env, op, 2, ref_point=rcc)
      !print *, "dbg| rcc build_modernmag_matrix", rcc
      !call dbcsr_print(op(1)%matrix, unit_nr=7)
      allocate(o(3))
      do i = 1, 3
        call cp_cfm_create(o(i)%matrix, o_struct, name="mago")
        call copy_dbcsr_to_fm(op(i)%matrix, or)
        call dbcsr_release(op(i)%matrix)
        deallocate(op(i)%matrix)
        call cp_fm_to_cfm(or, mtarget=o(i)%matrix)
      end do
      deallocate(op)
      call cp_fm_release(or)
      call cp_fm_struct_release(o_struct)

   end subroutine build_modernmag_matrix


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_reference_nuc_moments(moments, qqi, qqij, particles, qs_kinds, cell, rcc)
!
!      real(dp)                      :: moments(:), qqi(3), qqij(3, 3)
!      type(particle_type), pointer  :: particles(:)
!      type(qs_kind_type), pointer   :: qs_kinds(:)
!      type(cell_type), pointer      :: cell
!      real(dp)                      :: rcc(3)
!
!      complex(dp)                      :: d, di(3), dij(3, 3), iphase(3), ijphase(3, 3)
!      integer                          :: akind, i, ix, iy, iz, j, k, l
!      real(dp)                         :: dir(3), kvec(3), r(3), z, dd
!      type(atomic_kind_type), pointer  :: atomic_kind
!
!      do k = 0, 2
!         select case (k)
!         case (0)
!            ! charge
!            do i = 1, size(particles)
!               atomic_kind => particles(i)%atomic_kind
!               call get_atomic_kind(atomic_kind, kind_number=akind)
!               call get_qs_kind(qs_kinds(akind), core_charge=z)
!               moments(1) = moments(1) + z
!            end do
!            r = twopi * matmul(cell%h_inv, rcc)
!            iphase = exp(-gaussi * r)**moments(1)
!            do i = 1, 3
!               do j = 1, 3
!                  ijphase(i, j) = exp(-gaussi * r(i) * r(j))**moments(1)
!               end do
!            end do
!         case (1)
!            ! dipole
!            di = (1._dp, 0._dp)
!            do l = 1, size(particles)
!               r = pbc(particles(l)%r, cell)
!               atomic_kind => particles(l)%atomic_kind
!               call get_atomic_kind(atomic_kind, kind_number=akind)
!               call get_qs_kind(qs_kinds(akind), core_charge=z)
!               di = di * (exp(gaussi * twopi * matmul(cell%h_inv, r))**z)**z
!            end do
!            di = di * iphase
!            qqi = aimag(log(di))  ! for phase periodicity
!            moments(4:6) = matmul(cell%hmat / twopi, aimag(log(di)))
!         case (2)
!            ! quadrupole
!            dij = (1._dp, 0._dp)
!            do l = 1, size(particles)
!               r = pbc(particles(l)%r, cell)
!               atomic_kind => particles(l)%atomic_kind
!               call get_atomic_kind(atomic_kind, kind_number=akind)
!               call get_qs_kind(qs_kinds(akind), core_charge=z)
!               do i = 1, 3
!                  do j = 1, 3
!                     kvec(1) = twopi * cell%h_inv(i, i) * cell%h_inv(j, j)
!                     dd = kvec(1) * r(i) * r(j)
!                     d = exp(gaussi * dd)**z
!                     dij(i, j) = dij(i, j) * d
!                  end do
!               end do
!            end do
!            dij = dij * ijphase
!            qqij = aimag(log(dij))
!            l = 1
!            do i = 1, 3
!               do j = i, 3
!                  moments(l + 6) = cell%hmat(i, i) * cell%hmat(j, j) / twopi * aimag(log(dij(i, j)))
!                  l = l + 1
!               end do
!            end do
!         end select
!      end do
!
!   end subroutine calculate_reference_nuc_moments


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_reference_moments(moments, qs_env, rcc, particles, qs_kinds, cell, rho, &
!        mo_local)
!
!     real(dp)                            :: moments(:)
!     type(qs_environment_type), pointer  :: qs_env
!     real(dp)                            :: rcc(3)
!     type(particle_type), pointer        :: particles(:)
!     type(qs_kind_type), pointer         :: qs_kinds(:)
!     type(cell_type), pointer            :: cell
!     type(qs_rho_type), pointer          :: rho
!     type(cp_fm_p_type), pointer         :: mo_local(:)
!
!     complex(dp)  :: cd, cda, iphase(3), ijphase(3, 3), di(3), dij(3, 3)
!     integer                                  :: i, j, k, l, nao, nmo, tmp_dim
!     logical                                  :: do_rtp, found, tr
!     real(dp)                                 :: focc, kvec(3), qupoleder(9), strace, trace, &
!          r(3), zwfc, qqi(3), ci(3), qqij(3, 3), cij(3, 3), d
!     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
!     real(dp), pointer                        :: blk_p(:)
!     type(cp_cfm_type), pointer               :: m, om, mom
!     type(cp_cfm_p_type), pointer             :: o(:), mago(:), eigrmat(:)
!     type(dbcsr_p_type), pointer              :: cosmat(:), sinmat(:)
!     type(cp_fm_struct_type), pointer         :: m_struct, o_struct, mom_struct
!     type(cp_fm_p_type), pointer              :: op_fm_set(:, :), opvec(:)
!     type(cp_fm_type), pointer                :: mr, mi, or
!     type(dbcsr_iterator_type), pointer       :: iterator
!     type(dbcsr_p_type), pointer              :: oder(:, :), matrix_s(:), &
!          rho_ao(:), rho_ao_im(:)
!     type(dbcsr_type), pointer                :: densmat, densmat_im, ocos, osin, tmp
!     type(mo_set_p_type), pointer             :: mos(:)
!     type(neighbor_list_set_p_type), pointer  :: sab_all(:), sab_orb(:)
!     type(dbcsr_p_type), pointer  :: magop(:)
!
!     print *, "dbg| ref"
!     print *, "dbg| rcc calc_ref_mom", rcc
!     
!     ! info
!     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb, mos=mos)
!     do_rtp = (size(mo_local) == 2)
!     call cp_fm_get_info(mo_local(1)%matrix, nrow_global=nao, ncol_global=nmo)
!     focc = 1._dp
!     zwfc = 2._dp
!     trace = 0._dp
!
!     ! structs
!     m_struct => mo_local(1)%matrix%matrix_struct
!     call cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
!          template_fmstruct=mo_local(1)%matrix%matrix_struct)
!     call cp_fm_struct_create(mom_struct, nrow_global=nmo, ncol_global=nmo, &
!          template_fmstruct=mo_local(1)%matrix%matrix_struct)
!
!     ! get the mo coeffs
!     call cp_fm_create(mr, m_struct, name="mr")
!     if (do_rtp) then
!        call cp_fm_create(mi, m_struct, name="mi")
!     end if
!     call cp_cfm_create(m, m_struct, name="m")
!     call cp_fm_to_fm(mo_local(1)%matrix, mr)
!     if (do_rtp) then
!        call cp_fm_to_fm(mo_local(2)%matrix, mi)
!        call cp_fm_to_cfm(mr, mi, m)
!     else
!        call cp_fm_to_cfm(mr, mtarget=m)
!     end if
!     call cp_fm_release(mr)
!     if (do_rtp) then
!        call cp_fm_release(mi)
!     end if
!     call cp_fm_create(or, o_struct, name="mi")
!     call cp_fm_struct_release(o_struct)
!
!     ! operators
!     call build_wheeler_moment_matrix(qs_env, o, rcc, m)
!     call build_modernmag_matrix(qs_env, mago, rcc, m)
!     
!     ! get the nonlocal terms
!     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
!     nlrv = 0._dp
!     nlrrv = 0._dp
!     nlrxrv = 0._dp
!     call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, ref_point=rcc)
!     
!     ! nuclear contribution
!     call calculate_reference_nuc_moments(moments, qqi, qqij, particles, qs_kinds, cell, rcc)
!     
!     ! work matrices
!     call cp_cfm_create(om, m_struct)
!     call cp_cfm_create(mom, mom_struct)
!     call cp_fm_struct_release(mom_struct)
!
!     ! charge
!     moments(2) = -zwfc * nmo
!     moments(3) = moments(1) + moments(2)
!     ! initial phase from the reference point
!     r = twopi * (-moments(2)) * matmul(cell%h_inv, rcc)
!     iphase = exp(gaussi * r)
!     do i = 1, 3
!        do j = 1, 3
!           !ijphase(i, j) = exp(twopi * (-moments(2)) * cell%h_inv(i, i) * cell%h_inv(j, j) * rcc(i) * rcc(j))
!        end do
!     end do
!
!     ! edip
!     di = (1._dp, 0._dp)
!     do i = 1, 3
!        call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
!        call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!        d = 0._dp
!        do j = 1, nmo
!           call cp_cfm_get_element(mom, j, j, cd)
!           d = d + real(cd)
!        end do
!        moments(i+3) = moments(i+3) + d
!     end do
!     
!     ! mdip
!     do i = 1, 3
!        call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), mago(i)%matrix, m, (0._dp, 0._dp), om)
!        call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!        d = 0._dp
!        do j = 1, nmo
!           call cp_cfm_get_element(mom, j, j, cd)
!           d = d + aimag(cd)
!        end do
!        moments(i+12) = moments(i+12) - zwfc * d
!     end do
!     ! nl
!     moments(25:27) = moments(13:15) + nlrxrv
!
!     ! equp
!     dij = (1._dp, 0._dp)
!     k = 4
!     do i = 1, 3
!        do j = i, 3
!           call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), o(k)%matrix, m, (0._dp, 0._dp), om)
!           call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!           d = 0._dp
!           do k = 1, nmo
!              call cp_cfm_get_element(mom, k, k, cd)
!              d = d + real(cd)
!           end do
!        end do
!     end do
!     cij = aimag(log(dij))
!     l = 1
!     do i = 1, 3
!        do j = i, 3
!           if (qqij(i, j) + cij(i, j) > pi) cij(i, j) = cij(i, j) - twopi
!           if (qqij(i, j) + cij(i, j) < -pi) cij(i, j) = cij(i, j) + twopi
!           cij(i, j) = cell%hmat(i, i) * cell%hmat(j, j) / twopi
!           moments(l + 6) = moments(l + 6) + cij(i, j)
!           l = l + 1
!        end do
!     end do
!     do i = 4, 9
!        call cp_cfm_gemm("N", "N", nao, nmo, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
!        call cp_cfm_gemm("C", "N", nmo, nmo, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
!        d = 0.
!        do j = 1, nmo
!           call cp_cfm_get_element(mom, j, j, cd)
!           d = d + zwfc * real(cd)
!        end do
!        moments(i+3) = moments(i+3) - d
!     end do
!
!     do i = 1, size(o)
!        call cp_cfm_release(o(i)%matrix)
!     end do
!     do i = 1, size(mago)
!        call cp_cfm_release(mago(i)%matrix)
!     end do
!     call cp_cfm_release(m)
!     call cp_cfm_release(om)
!     call cp_cfm_release(mom)
!     call cp_fm_release(or)
!     deallocate(mago, o)
!     moments = 0._dp
!     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, 1, &
!          size(particles), rcc)
!     call calculate_localized_loc_moments(moments, qs_env, rcc, mo_local, &
!          particles, qs_kinds, cell, 1, size(particles), ref=.true.)
!
!   end subroutine calculate_reference_moments


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
!        last_atom, rcc)
!
!      real(dp)                      :: moments(:)
!      type(particle_type), pointer  :: particles(:)
!      type(qs_kind_type), pointer   :: qs_kinds(:)
!      type(cell_type), pointer      :: cell
!      integer                       :: first_atom, last_atom
!      real(dp)                      :: rcc(3)
!
!      logical                          :: do_periodic
!      integer                          :: akind, iatom, l, lx, ly, lz
!      logical                          :: ghost, floating, do_ref
!      real(dp)                         :: d, r(3), z
!      type(atomic_kind_type), pointer  :: atomic_kind
!      
!      do iatom = first_atom, last_atom
!         r = pbc(particles(iatom)%r - rcc, cell) + rcc
!         r = r - rcc
!         atomic_kind => particles(iatom)%atomic_kind
!         call get_atomic_kind(atomic_kind, kind_number=akind)
!         call get_qs_kind(qs_kinds(akind), core_charge=z, ghost=ghost, floating=floating)
!         if (.not. ghost .and. .not. floating) then
!            r = pbc(particles(iatom)%r - rcc, cell)
!            moments(1) = moments(1) + z  ! nuclear charge
!            do l = 1, 9
!               lx = indco(1, l+1)
!               ly = indco(2, l+1)
!               lz = indco(3, l+1)
!               d = 1._dp
!               if (lx > 0) d = d * r(1)**lx
!               if (ly > 0) d = d * r(2)**ly
!               if (lz > 0) d = d * r(3)**lz
!               moments(l+3) = moments(l+3) + z * d
!            end do
!         end if
!      end do
!      moments(3) = moments(3) + moments(2) ! total charge
!
!   end subroutine calculate_localized_nuc_moments


! **************************************************************************************************
! **************************************************************************************************
   subroutine calculate_localized_loc_moments(moments, qs_env, rcc, mo_local, particles, qs_kinds, &
        cell, first_atom, last_atom, statelist, periodic, ref)

     real(dp)                            :: moments(:)
     type(qs_environment_type), pointer  :: qs_env
     real(dp)                            :: rcc(3)
     type(cp_fm_p_type), pointer         :: mo_local(:)
     type(particle_type), pointer        :: particles(:)
     type(qs_kind_type), pointer         :: qs_kinds(:)
     type(cell_type), pointer            :: cell
     integer                             :: first_atom, last_atom
     integer, pointer, optional          :: statelist(:)
     logical, optional                   :: periodic, ref

     complex(dp)                              :: cd
     integer                                  :: i, j, k, nao, nmo, ns
     integer, pointer                         :: states(:)
     logical                                  :: do_rtp, do_periodic, do_ref
     real(dp)                                 :: d, zwfc, qupoleder(9)
     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
     type(cp_cfm_p_type), pointer             :: o(:), oder(:, :), mago(:), momo(:)
     type(cp_cfm_type), pointer               :: m, om, mom
     type(dbcsr_p_type), pointer              :: matrix_s(:), op(:), opder(:, :), magop(:), momop(:)
     type(cp_fm_struct_type), pointer         :: m_struct, o_struct, mom_struct
     type(cp_fm_type), pointer                :: mr, mi, or
     type(neighbor_list_set_p_type), pointer  :: sab_all(:), sab_orb(:)

     !print *, "dbg| loc"
     !print *, "dbg| calc_loc_mom rcc", rcc
      
     ! info
     if (present(periodic)) then
        do_periodic = periodic
     else
        do_periodic = .false.
     end if
     if (present(ref)) then
        do_ref = ref
     else
        do_ref = .false.
     end if
     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
     do_rtp = (size(mo_local) == 2)
     call cp_fm_get_info(mo_local(1)%matrix, nrow_global=nao, ncol_global=nmo)
     zwfc = 2._dp
     allocate(states(nmo))
     if (present(statelist)) then
        ns = size(statelist)
        do i = 1, ns
           states(i) = statelist(i)
        end do
     else
        ns = nmo
        do i = 1, ns
           states(i) = i
        end do
     end if

     ! structs
     call cp_fm_struct_create(m_struct, nrow_global=nao, ncol_global=ns, &
          template_fmstruct=mo_local(1)%matrix%matrix_struct)
     call cp_fm_struct_create(o_struct, nrow_global=nao, ncol_global=nao, &
          template_fmstruct=mo_local(1)%matrix%matrix_struct)
     call cp_fm_struct_create(mom_struct, nrow_global=ns, ncol_global=ns, &
          template_fmstruct=mo_local(1)%matrix%matrix_struct)

     ! get only this subsystem's states
     call cp_fm_create(mr, m_struct, name="mr")
     call cp_fm_set_all(mr, 0._dp)
     if (do_rtp) then
        call cp_fm_create(mi, m_struct, name="mi")
        call cp_fm_set_all(mi, 0._dp)
     end if
     call cp_cfm_create(m, m_struct, name="m")
     do i = 1, ns
        call cp_fm_to_fm(mo_local(1)%matrix, mr, 1, states(i), i)
        if (do_rtp) then
           call cp_fm_to_fm(mo_local(2)%matrix, mi, 1, states(i), i)
           call cp_fm_to_cfm(mr, mi, m)
        else
           call cp_fm_to_cfm(mr, mtarget=m)
        end if
     end do
     call cp_fm_release(mr)
     if (do_rtp) then
        call cp_fm_release(mi)
     end if

     !! get the operators
     call cp_fm_create(or, o_struct, name="or")
     ! el moments with derivatives
     nullify(op)
     call dbcsr_allocate_matrix_set(op, nm)
     do i = 1, nm
        call dbcsr_init_p(op(i)%matrix)
        call dbcsr_create(op(i)%matrix, template=matrix_s(1)%matrix, &
             matrix_type=dbcsr_type_no_symmetry)
        call cp_dbcsr_alloc_block_from_nbl(op(i)%matrix, sab_all)
        call dbcsr_set(op(i)%matrix, 0._dp)
     end do
     nullify(opder)
     call dbcsr_allocate_matrix_set(opder, 3, 3)
     do i = 1, 3
        do j = 1, 3
           call dbcsr_init_p(opder(i, j)%matrix)
           call dbcsr_create(opder(i, j)%matrix, template=matrix_s(1)%matrix, &
                matrix_type=dbcsr_type_no_symmetry)
           call cp_dbcsr_alloc_block_from_nbl(opder(i, j)%matrix, sab_all)
           call dbcsr_set(opder(i, j)%matrix, 0._dp)
        end do
     end do
     allocate(o(9))
     if (ref) then
        call build_local_moments_der_matrix(qs_env, op, opder, 2, 1, rcc)
     else
        call build_local_moment_matrix(qs_env, op, 2, rcc)
     end if
     do i = 1, nm
        call cp_cfm_create(o(i)%matrix, o_struct)
        if (ref) then
           call copy_dbcsr_to_fm(op(i)%matrix, or)
        end if
        call dbcsr_release(op(i)%matrix)
        deallocate(op(i)%matrix)
        if (ref) then
           call cp_fm_to_cfm(or, mtarget=o(i)%matrix)
        else
           call cp_cfm_set_all(o(i)%matrix, (0._dp, 0._dp))
        end if
     end do
     deallocate(op)
     allocate(oder(3, 3))
     do i = 1, 3
        do j = 1, 3
           call cp_cfm_create(oder(i, j)%matrix, o_struct)
           if (ref) then
              call copy_dbcsr_to_fm(opder(i, j)%matrix, or)
           end if
           call dbcsr_release(opder(i, j)%matrix)
           deallocate(opder(i, j)%matrix)
           if (ref) then
              call cp_fm_to_cfm(or, mtarget=oder(i, j)%matrix)
           else
              call cp_cfm_set_all(oder(i, j)%matrix, (0._dp, 0._dp))
           end if
        end do
     end do
     deallocate(opder)

     ! momentum
     nullify(momop)
     call dbcsr_allocate_matrix_set(momop, 3)
     do i = 1, 3
        call dbcsr_init_p(momop(i)%matrix)
        call dbcsr_create(momop(i)%matrix, template=matrix_s(1)%matrix, &
             matrix_type=dbcsr_type_antisymmetric)
        call cp_dbcsr_alloc_block_from_nbl(momop(i)%matrix, sab_orb)
        call dbcsr_set(momop(i)%matrix, 0._dp)
     end do
     allocate(momo(3))
     call build_lin_mom_matrix(qs_env, momop)
     do i = 1, 3
        call cp_cfm_create(momo(i)%matrix, o_struct)
        call copy_dbcsr_to_fm(momop(i)%matrix, or)
        call dbcsr_release(momop(i)%matrix)
        deallocate(momop(i)%matrix)
        call cp_fm_to_cfm(or, mtarget=momo(i)%matrix)
     end do
     deallocate(momop)

     ! mag moments
     nullify(magop)
     call dbcsr_allocate_matrix_set(magop, 3)
     do i = 1, 3
        call dbcsr_init_p(magop(i)%matrix)
        call dbcsr_create(magop(i)%matrix, template=matrix_s(1)%matrix, &
             matrix_type=dbcsr_type_antisymmetric)
        call cp_dbcsr_alloc_block_from_nbl(magop(i)%matrix, sab_orb)
        call dbcsr_set(magop(i)%matrix, 0._dp)
     end do
     nullify(op)
     allocate(mago(3))
     call build_local_magmom_matrix(qs_env, magop, nmm, ref_point=rcc)
     do i = 1, nmm
        call cp_cfm_create(mago(i)%matrix, o_struct, name="mago")
        if (ref) then
           call copy_dbcsr_to_fm(magop(i)%matrix, or)
        end if
        call dbcsr_release(magop(i)%matrix)
        deallocate(magop(i)%matrix)
        if (ref) then
           call cp_fm_to_cfm(or, mtarget=mago(i)%matrix)
        else
           call cp_cfm_set_all(mago(i)%matrix, (0._dp, 0._dp))
        end if
     end do
     deallocate(magop)
     call cp_fm_release(or)

     ! get the nonlocal terms
     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
     nlrv = 0._dp
     nlrrv = 0._dp
     nlrxrv = 0._dp
     if (ref) then
        call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, rcc)
     end if

     ! nuclear contribution
     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
          last_atom, rcc)

     ! expectation values
     call cp_cfm_create(om, m_struct, name="om")
     call cp_cfm_create(mom, mom_struct, name="mom")
     call cp_fm_struct_release(m_struct)
     call cp_fm_struct_release(o_struct)
     call cp_fm_struct_release(mom_struct)

     ! charge
     call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, m, (0._dp, 0._dp), mom)
     d = 0.
     do i = 1, ns
        call cp_cfm_get_element(mom, i, i, cd)
        moments(2) = moments(2) + real(cd)
     end do
     moments(2) = -zwfc * ns! moments(2)
     moments(3) = moments(1) + moments(2)

     ! edip
     do i = 1, 3
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
        d = 0.
        do j = 1, ns
           call cp_cfm_get_element(mom, j, j, cd)
           d = d + real(cd)
        end do
        moments(i+3) = moments(i+3) - zwfc * d
     end do
     ! vel
     do i = 1, 3
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), momo(i)%matrix, m, (0._dp, 0._dp), om)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
        d = 0.
        do j = 1, ns
           call cp_cfm_get_element(mom, j, j, cd)
           d = d + aimag(cd)
        end do
        moments(i+15) = moments(i+15) - zwfc * d
     end do
     ! nl
     moments(28:30) = moments(16:18) + nlrv

     ! mdip
     do i = 1, 3
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), mago(i)%matrix, m, (0._dp, 0._dp), om)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
        d = 0.
        do j = 1, ns
           call cp_cfm_get_element(mom, j, j, cd)
           d = d + aimag(cd)
        end do
        moments(i+12) = moments(i+12) - zwfc * d
     end do
     ! nl
     moments(25:27) = moments(13:15) + nlrxrv

     ! equp
     do i = 4, 9
        call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), o(i)%matrix, m, (0._dp, 0._dp), om)
        call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
        d = 0.
        do j = 1, ns
           call cp_cfm_get_element(mom, j, j, cd)
           d = d + zwfc * real(cd)
        end do
        moments(i+3) = moments(i+3) - d
     end do
     ! vel
     qupoleder = 0._dp
     do i = 1, 3
        do j = 1, 3
           call cp_cfm_gemm("N", "N", nao, ns, nao, (1._dp, 0._dp), oder(i, j)%matrix, m, (0._dp, 0._dp), om)
           call cp_cfm_gemm("C", "N", ns, ns, nao, (1._dp, 0._dp), m, om, (0._dp, 0._dp), mom)
           d = 0.
           do k = 1, ns
              call cp_cfm_get_element(mom, k, k, cd)
              d = d - aimag(cd)
              if (i == j) then
                 d = d + real(cd)
              end if
           end do
           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + zwfc * d
        end do
     end do
     moments(19) = -2*qupoleder(1) - moments(2)
     moments(20) = -qupoleder(2) - qupoleder(4)
     moments(21) = -qupoleder(3) - qupoleder(7)
     moments(22) = -2*qupoleder(5) - moments(2)
     moments(23) = -qupoleder(6) - qupoleder(8)
     moments(24) = -2*qupoleder(9) - moments(2)
     ! nl
     moments(31:36) = moments(19:24) + nlrrv

     ! cleanup
     call cp_cfm_release(m)
     call cp_cfm_release(om)
     call cp_cfm_release(mom)
     do i = 1, size(o)
        call cp_cfm_release(o(i)%matrix)
     end do
     deallocate(o)
     do i = 1, size(oder, 1)
        do j = 1, size(oder, 2)
           call cp_cfm_release(oder(i, j)%matrix)
        end do
     end do
     deallocate(oder)
     do i = 1, size(momo)
        call cp_cfm_release(momo(i)%matrix)
     end do
     deallocate(momo)
     do i = 1, size(mago)
        call cp_cfm_release(mago(i)%matrix)
     end do
     deallocate(mago)
     deallocate(nlrv, nlrrv, nlrxrv, states)
     
   end subroutine calculate_localized_loc_moments


!! **************************************************************************************************
!! **************************************************************************************************
!   subroutine calculate_localized_kg_moments(moments, qs_env, rcc, rho, particles, qs_kinds, cell, &
!        first_atom, last_atom, do_periodic)
!
!     real(dp)                            :: moments(:)
!     type(qs_environment_type), pointer  :: qs_env
!     real(dp)                            :: rcc(3)
!     type(qs_rho_type), pointer          :: rho
!     type(particle_type), pointer        :: particles(:)
!     type(qs_kind_type), pointer         :: qs_kinds(:)
!     type(cell_type), pointer            :: cell
!     integer                             :: first_atom, last_atom
!     logical  :: do_periodic
!
!     integer                                  :: i, j
!     logical                                  :: found, tr
!     real(dp)                                 :: focc, qupoleder(9), strace, trace
!     real(dp), allocatable                    :: nlrv(:), nlrrv(:), nlrxrv(:)
!     real(dp), pointer                        :: blk_p(:)
!     type(dbcsr_iterator_type), pointer       :: iterator
!     type(dbcsr_p_type), pointer              :: o(:), oder(:, :), mago(:), momo(:), matrix_s(:), &
!          rho_ao(:), rho_ao_im(:)
!     type(dbcsr_type), pointer                :: densmat, densmat_im, tmp
!     type(neighbor_list_set_p_type), pointer  :: sab_all(:), sab_orb(:)
!
!     ! info
!     call get_qs_env(qs_env, matrix_s=matrix_s, sab_all=sab_all, sab_orb=sab_orb)
!     call qs_rho_get(rho, rho_ao=rho_ao, rho_ao_im=rho_ao_im)
!     focc = 1._dp
!
!     ! get real and imaginary part of the subdensity
!     allocate(densmat)
!     call dbcsr_copy(densmat, rho_ao(1)%matrix)
!     allocate(iterator)
!     call dbcsr_iterator_start(iterator, densmat)
!     do while (dbcsr_iterator_blocks_left(iterator))
!        call dbcsr_iterator_next_block(iterator, i, j, blk_p)
!        if ((i < first_atom) .or. (i > last_atom) .or. (j < first_atom) .or. (j > last_atom)) then
!           call dbcsr_get_block_p(densmat, i, j, blk_p, tr, found)
!           blk_p = 0._dp
!        end if
!     end do
!     call dbcsr_iterator_stop(iterator)
!     allocate(densmat_im)
!     call dbcsr_copy(densmat_im, rho_ao_im(1)%matrix)
!     call dbcsr_iterator_start(iterator, densmat_im)
!     do while (dbcsr_iterator_blocks_left(iterator))
!        call dbcsr_iterator_next_block(iterator, i, j, blk_p)
!        if ((i < first_atom) .or. (i > last_atom) .or. (j < first_atom) .or. (j > last_atom)) then
!           call dbcsr_get_block_p(densmat_im, i, j, blk_p, tr, found)
!           blk_p = 0._dp
!        end if
!     end do
!     call dbcsr_iterator_stop(iterator)
!     deallocate(iterator)
!
!     ! get the operators
!     ! el moments with derivatives
!     nullify(o)
!     call dbcsr_allocate_matrix_set(o, nm)
!     do i = 1, nm
!        call dbcsr_init_p(o(i)%matrix)
!        call dbcsr_create(o(i)%matrix, template=matrix_s(1)%matrix, &
!             matrix_type=dbcsr_type_no_symmetry)
!        call cp_dbcsr_alloc_block_from_nbl(o(i)%matrix, sab_all)
!        call dbcsr_set(o(i)%matrix, 0._dp)
!     end do
!     nullify(oder)
!     call dbcsr_allocate_matrix_set(oder, 3, 3)
!     do i = 1, 3
!        do j = 1, 3
!           call dbcsr_init_p(oder(i, j)%matrix)
!           call dbcsr_create(oder(i, j)%matrix, template=matrix_s(1)%matrix, &
!                matrix_type=dbcsr_type_no_symmetry)
!           call cp_dbcsr_alloc_block_from_nbl(oder(i, j)%matrix, sab_all)
!           call dbcsr_set(oder(i, j)%matrix, 0._dp)
!        end do
!     end do
!     call build_local_moments_der_matrix(qs_env, o, oder, 2, 1, rcc)
!
!     ! momentum
!     nullify(momo)
!     call dbcsr_allocate_matrix_set(momo, 3)
!     do i = 1, 3
!        call dbcsr_init_p(momo(i)%matrix)
!        call dbcsr_create(momo(i)%matrix, template=matrix_s(1)%matrix, &
!             matrix_type=dbcsr_type_antisymmetric)
!        call cp_dbcsr_alloc_block_from_nbl(momo(i)%matrix, sab_orb)
!        call dbcsr_set(momo(i)%matrix, 0._dp)
!     end do
!     call build_lin_mom_matrix(qs_env, momo)
!
!     ! mag moments
!     nullify(mago)
!     call dbcsr_allocate_matrix_set(mago, 3)
!     do i = 1, 3
!        call dbcsr_init_p(mago(i)%matrix)
!        call dbcsr_create(mago(i)%matrix, template=matrix_s(1)%matrix, &
!             matrix_type=dbcsr_type_antisymmetric)
!        call cp_dbcsr_alloc_block_from_nbl(mago(i)%matrix, sab_orb)
!        call dbcsr_set(mago(i)%matrix, 0._dp)
!     end do
!     call build_local_magmom_matrix(qs_env, mago, nmm, ref_point=rcc)
!
!     ! get the nonlocal terms
!     allocate(nlrv(3), nlrrv(6), nlrxrv(3))
!     nlrv = 0._dp
!     nlrrv = 0._dp
!     nlrxrv = 0._dp
!     call calculate_commutator_nl_terms(qs_env, nlrv, nlrxrv, nlrrv, rcc)
!     
!     ! nuclear contribution
!     call calculate_localized_nuc_moments(moments, particles, qs_kinds, cell, first_atom, &
!          last_atom, rcc)
!
!     ! temp matrix
!     allocate(tmp)
!     call dbcsr_desymmetrize(matrix_s(1)%matrix, tmp)
!     
!     ! charge
!     call dbcsr_dot(densmat, matrix_s(1)%matrix, trace)
!     moments(2) = - focc * trace
!     moments(3) = moments(1) + moments(2)
!
!     ! edip
!     do i = 1, 3
!        call dbcsr_multiply("T", "N", 1._dp, densmat, o(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+3) = moments(i+3) - focc * trace
!     end do
!     ! vel
!     do i = 1, 3
!        call dbcsr_multiply("T", "N", 1._dp, densmat_im, momo(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+15) = moments(i+15) + focc * trace
!     end do
!     ! nl
!     moments(28:30) = moments(16:18) + nlrv
!
!     ! mdip
!     do i = 1, 3
!        call dbcsr_multiply("T", "N", 1._dp, densmat_im, mago(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+12) = moments(i+12) + focc * trace
!     end do
!     ! nl
!     moments(25:27) = moments(13:15) + nlrxrv
!
!     ! equp
!     do i = 4, 9
!        call dbcsr_multiply("T", "N", 1._dp, densmat, o(i)%matrix, 0._dp, tmp)
!        call dbcsr_trace(tmp, trace)
!        moments(i+3) = moments(i+3) - focc * trace
!     end do
!     ! vel
!     qupoleder = 0._dp
!     do i = 1, 3
!        do j = 1, 3
!           call dbcsr_multiply("T", "N", 1._dp, densmat, oder(i, j)%matrix, 0._dp, tmp)
!           call dbcsr_trace(tmp, trace)
!           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + focc * trace
!           call dbcsr_multiply("T", "N", 1._dp, densmat_im, oder(i, j)%matrix, 0._dp, tmp)
!           call dbcsr_trace(tmp, trace)
!           qupoleder((i-1)*3+j) = qupoleder((i-1)*3+j) + focc * trace
!        end do
!     end do
!     moments(19) = -2*qupoleder(1) - moments(2)
!     moments(20) = -qupoleder(2) - qupoleder(4)
!     moments(21) = -qupoleder(3) - qupoleder(7)
!     moments(22) = -2*qupoleder(5) - moments(2)
!     moments(23) = -qupoleder(6) - qupoleder(8)
!     moments(24) = -2*qupoleder(9) - moments(2)
!     ! nl
!     moments(31:36) = moments(19:24) + nlrrv
!
!     do i = 1, size(o)
!        call dbcsr_release(o(i)%matrix)
!        deallocate(o(i)%matrix)
!     end do
!     deallocate(o)
!     do i = 1, size(oder, 1)
!        do j = 1, size(oder, 2)
!           call dbcsr_release(oder(i, j)%matrix)
!           deallocate(oder(i, j)%matrix)
!        end do
!     end do
!     deallocate(oder)
!     do i = 1, size(momo)
!        call dbcsr_release(momo(i)%matrix)
!        deallocate(momo(i)%matrix)
!     end do
!     deallocate(momo)
!     do i = 1, size(mago)
!        call dbcsr_release(mago(i)%matrix)
!        deallocate(mago(i)%matrix)
!     end do
!     deallocate(mago)
!     call dbcsr_release(densmat)
!     deallocate(densmat)
!     call dbcsr_release(densmat_im)
!     deallocate(densmat_im)
!     call dbcsr_release(tmp)
!     deallocate(tmp)
!     deallocate(nlrv, nlrrv, nlrxrv)
!
!   end subroutine calculate_localized_kg_moments


! **************************************************************************************************
! **************************************************************************************************
   subroutine print_localized_moments(logger, input, moments, rcc, molname, suffix)

      type(cp_logger_type), pointer     :: logger
      type(section_vals_type), pointer  :: input
      real(dp)                          :: moments(:), rcc(:)
      character(*)                      :: molname
      character(*), optional            :: suffix

      integer                               :: output_unit
      character(len=default_string_length)  :: iter

      if (present(suffix)) then
         output_unit = cp_print_key_unit_nr(logger=logger, basis_section=input, &
              print_key_path="DFT%PRINT%MOMENTS", extension=".csv", &
              middle_name="locmoments_"//trim(molname)//"_"//trim(suffix), log_filename=.false.)
      else
         output_unit = cp_print_key_unit_nr(logger=logger, basis_section=input, &
              print_key_path="DFT%PRINT%MOMENTS", extension=".csv", &
              middle_name="locmoments_"//trim(molname), log_filename=.false.)
      end if
      iter = cp_iter_string(logger%iter_info)
      
      if (output_unit > 0) then
         if (iteration == 1) then
            write (output_unit, "(A)") header
         end if
         write (output_unit, "(I6,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A," // &
              "ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3,A,ES20.12E3)") &
              iteration,   ", "//trim(iter)//", ", &
              moments( 1), ", ", moments( 2), ", ", moments( 3), ", ", &
              moments( 4), ", ", moments( 5), ", ", moments( 6), ", ", moments( 7), ", ", &
              moments( 8), ", ", moments( 9), ", ", moments(10), ", ", moments(11), ", ", &
              moments(12), ", ", moments(13), ", ", moments(14), ", ", moments(15), ", ", &
              moments(16), ", ", moments(17), ", ", moments(18), ", ", moments(19), ", ", &
              moments(20), ", ", moments(21), ", ", moments(22), ", ", moments(23), ", ", &
              moments(24), ", ", moments(25), ", ", moments(26), ", ", moments(27), ", ", &
              moments(28), ", ", moments(29), ", ", moments(30), ", ", moments(31), ", ", &
              moments(32), ", ", moments(33), ", ", moments(34), ", ", moments(35), ", ", &
              moments(36), ", ", rcc(1), ", ", rcc(2), ",", rcc(3)
      end if

      call cp_print_key_finished_output(output_unit, logger, input, "DFT%PRINT%MOMENTS")

   end subroutine print_localized_moments

end module localized_moments
 

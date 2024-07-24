module fat_methods
  use atomic_kind_types
  use cp_control_types
  use cp_files
  use cp_log_handling
  use cp_realspace_grid_cube
  use dbcsr_api
  use fat_types
  use force_env_types
  use input_section_types
  use kinds
  use particle_list_types
  use particle_types
  use physcon
  use pw_env_types
  use pw_poisson_methods
  use pw_poisson_types
  use pw_pool_types
  use pw_methods
  use pw_types
  use qs_collocate_density
  use qs_energy_types
  use qs_environment_types
  use qs_force
  use qs_force_types
  use qs_rho_types
  use qs_subsys_types
  use realspace_grid_cube
  use string_utilities
  use xc

#include "./base/base_uses.f90"

  implicit none

  private

  public :: &
       fat_finalize, &
       fat_init, &
       fat_is_done, &
       fat_iterate, &
       fat_results

  logical, private, parameter :: is_debug = .true.

  character(len=*), private, parameter :: &
       modulen = "fat_methods", &
       fat_title_fmt = '(/, T2, "FAT| ", 4(A6, ", "), 9(A26, ", "), A26)', &
       fat_data_fmt = '(/, T2, "FAT| ", 4(I6, ", "), 9(F26.15, ", "), F26.15)'

contains

! **************************************************************************************************
!> \brief Creates a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_create_frag(fat_env, pw, i_frag, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw
    integer :: i_frag
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pw_create_frag"

    integer :: use_data, in_space
    logical :: my_rspace

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    if (present(rspace)) then
       my_rspace = rspace
    else
       my_rspace = .true.
    end if

    if (my_rspace) then
       use_data = realdata3d
       in_space = realspace
    else
       use_data = complexdata1d
       in_space = reciprocalspace
    end if

    allocate(pw)
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(pw, use_data=use_data, in_space=in_space)
    call pw_zero(pw)

  end subroutine fat_pw_create_frag

! **************************************************************************************************
!> \brief Creates a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_create_ref(fat_env, pw, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pw_create_ref"

    integer :: use_data, in_space
    logical :: my_rspace

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    if (present(rspace)) then
       my_rspace = rspace
    else
       my_rspace = .true.
    end if

    if (my_rspace) then
       use_data = realdata3d
       in_space = realspace
    else
       use_data = complexdata1d
       in_space = reciprocalspace
    end if

    allocate(pw)
    call fat_env%pw_pool_ref%create_pw(pw, use_data=use_data, in_space=in_space)
    call pw_zero(pw)

  end subroutine fat_pw_create_ref
  
! **************************************************************************************************
!> \brief Releases a pw for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_release_frag(fat_env, pw, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_pw_release_frag"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(pw)
    deallocate(pw)

    if (is_debug) then
       print *, "DBG| released"
    end if

  end subroutine fat_pw_release_frag

! **************************************************************************************************
!> \brief Releases a pw for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pw_release_ref(fat_env, pw)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw

    character(len=*), parameter :: routinen = "fat_pw_release_ref"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_env%pw_pool_ref%give_back_pw(pw)
    deallocate(pw)

    if (is_debug) then
       print *, "DBG| released"
    end if

  end subroutine fat_pw_release_ref

! **************************************************************************************************
!> \brief Print a plane wave to a cube file
! **************************************************************************************************
  subroutine fat_pw_print(fat_env, v, label)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v
    character(len=*) :: label

  end subroutine fat_pw_print
  
! **************************************************************************************************
!> \brief Creates an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_create_frag(fat_env, pws, i_frag, n, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)
    integer :: i_frag, n
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pws_create_frag"

    integer :: i
    type(pw_type), pointer :: pw_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    allocate(pws(n))
    do i = 1, n
       pw_tmp = pws(i)
       call fat_pw_create_frag(fat_env, pw_tmp, i_frag, rspace)
    end do

  end subroutine fat_pws_create_frag
  
! **************************************************************************************************
!> \brief Creates an array of pws for the reference system in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_create_ref(fat_env, pws, n, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)
    integer :: n
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_pws_create_ref"

    integer :: i
    type(pw_p_type), pointer :: pw_tmp(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    allocate(pws(n))
    allocate(pw_tmp(n))
    do i = 1, n
       call fat_pw_create_ref(fat_env, pw_tmp(i)%pw, rspace)
       pws(i) = pw_tmp(i)%pw
    end do
    deallocate(pw_tmp)
    
  end subroutine fat_pws_create_ref

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_release_frag(fat_env, pws, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_pws_release_frag"

    integer :: i, n
    type(pw_type), pointer :: pw_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    n = size(pws)
    do i = 1, n
       pw_tmp => pws(i)
       call fat_pw_release_frag(fat_env, pw_tmp, i_frag)
    end do
    
    if (is_debug) then
       print *, "DBG| released"
    end if

  end subroutine fat_pws_release_frag

! **************************************************************************************************
!> \brief Releases an array of pws for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_pws_release_ref(fat_env, pws)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pws(:)

    character(len=*), parameter :: routinen = "fat_pws_release_ref"

    integer :: i, n
    type(pw_p_type), pointer :: pw_tmp(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    n = size(pws)
    allocate(pw_tmp(n))
    do i = 1, n
       pw_tmp(i)%pw => pws(i)
       call fat_pw_release_ref(fat_env, pw_tmp(i)%pw)
    end do
    deallocate(pw_tmp)
    
    if (is_debug) then
       print *, "DBG| released"
    end if

  end subroutine fat_pws_release_ref

! **************************************************************************************************
!> \brief Gets the electron density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_rho_frag(fat_env, rho, i_frag, rspace)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: rho(:)
    integer :: i_frag
    logical, optional :: rspace

    character(len=*), parameter :: routinen = "fat_rho_frag"

    character(len=default_string_length) :: label
    logical :: my_rspace
    integer :: i, in_space, use_data
    type(pw_type), pointer :: rho_r(:), rho_tmp
    type(qs_environment_type), pointer :: qs_env
    type(qs_rho_type), pointer :: qs_rho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    if (present(rspace)) then
       my_rspace = rspace
    else
       my_rspace = .true.
    end if

    if (my_rspace) then
       use_data = realdata3d
       in_space = realspace
    else
       use_data = complexdata1d
       in_space = reciprocalspace
    end if

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho=qs_rho)
    if (my_rspace) then
       call qs_rho_get(qs_rho, rho_r=rho)
    else
       call qs_rho_get(qs_rho, rho_g=rho)
    end if

    if (fat_env%extern(i_frag)) then
       call integer_to_string(i_frag, label)
       if (my_rspace) then
          call cp_cube_to_pw(rho(1), "extern_"//trim(label)//".cube", 1.0_dp, .false.)
       else
          call fat_pw_create_frag(fat_env, rho_tmp, i_frag)
          call cp_cube_to_pw(rho_tmp, "extern_"//trim(label)//".cube", 1.0_dp, .false.)
          call pw_transfer(rho_tmp, rho(1))
          call fat_pw_release_frag(fat_env, rho_tmp, i_frag)
       end if
    end if

    if (is_debug) then
       if (my_rspace) then
          print *, "DBG| rho ", i_frag, " sum ", sum(rho(1)%cr3d)
       else
          print *, "DBG| rho ", i_frag, " sum ", sum(rho(1)%cc)
       end if
    end if

  end subroutine fat_rho_frag

! **************************************************************************************************
!> \brief Gets the nuclear density for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_rhon_frag(fat_env, rhon, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: rhon
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_rhon_frag"

    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: pw_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env=qs_env, rho_core=rhon)

  end subroutine fat_rhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the nuclear density w.r.t. a nuclear position for one fragment in the
!>        freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_drhon_frag(fat_env, drhon, i_frag, i_nuc, i_dir)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: drhon
    integer :: i_frag, i_nuc, i_dir

    character(len=*), parameter :: routinen = "fat_rhon_frag"

    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call calculate_drho_core(drhon, qs_env, i_dir, i_nuc)

  end subroutine fat_drhon_frag

! **************************************************************************************************
!> \brief Gets the gradient of the electronic density w.r.t. a nuclear position for one fragment in
!>        the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_drho_frag(fat_env, drho, i_frag, i_nuc, i_dir)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: drho
    integer :: i_frag, i_nuc, i_dir

    character(len=*), parameter :: routinen = "fat_drho_frag"

    type(dbcsr_p_type), pointer :: p_mat(:), p_mat_kp(:)
    type(pw_type), pointer :: drho_g
    type(qs_environment_type), pointer :: qs_env
    type(qs_rho_type), pointer :: rho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho=rho)
    call qs_rho_get(rho, rho_ao=p_mat)
    call fat_pw_create_frag(fat_env, drho_g, i_frag, .false.)
    call calculate_drho_elec_dR( &
         p_mat(1)%matrix, drho=drho, drho_gspace=drho_g, qs_env=qs_env, beta=i_dir, lambda=i_nuc)
    call fat_pw_release_frag(fat_env, drho_g, i_frag)

  end subroutine fat_drho_frag

! **************************************************************************************************
!> \brief Gets the total electron density for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_rho_ref(fat_env, pw)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: pw(:)

    character(len=*), parameter :: routinen = "fat_rho_ref"

    integer :: i, j, n
    type(pw_type), pointer :: pw_tmp(:), rho_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    do i = 1, fat_env%n_frag
       call fat_rho_frag(fat_env, pw_tmp, i)
       do j = 1, size(pw_tmp)
          call pw_axpy(pw_tmp(j), pw(j), allow_noncompatible_grids=.true.)
       end do
    end do

  end subroutine fat_rho_ref

! **************************************************************************************************
!> \brief Calculate kinetic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_ts_frag(fat_env, v, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_ts_frag"

    real(dp) :: e, virial(3, 3)
    type(pw_env_type), pointer :: pw_env
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(pw_pool_type), pointer :: pw_pool
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, section

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")
    pw_pool => fat_env%pw_pool_sub(i_frag)%pool

    call fat_rho_frag(fat_env, rho_r, i_frag)
    nullify(v, v_tau)
    call xc_vxc_pw_create(v, v_tau, e, rho_r, rho_g, tau, section, pw_pool, .false., virial)
    fat_env%ekin_sub(i_frag) = e

  end subroutine fat_v_ts_frag

! **************************************************************************************************
!> \brief Calculate reference kinetic potential for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_ts_ref(fat_env, v)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)

    character(len=*), parameter :: routinen = "fat_v_ts_ref"

    real(dp) :: e, virial(3, 3)
    type(pw_env_type), pointer :: pw_env
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(pw_pool_type), pointer :: pw_pool
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, section

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    section => section_vals_get_subs_vals(input, "DFT%QS%OPT_EMBED%XC")
    pw_pool => fat_env%pw_pool_ref

    call fat_pws_create_ref(fat_env, rho_r, 1)
    call fat_rho_ref(fat_env, rho_r)
    nullify(v, v_tau)
    call xc_vxc_pw_create(v, v_tau, e, rho_r, rho_g, tau, section, pw_pool, .false., virial)
    fat_env%ekin_ref = e
    call fat_pws_release_ref(fat_env, rho_r)

  end subroutine fat_v_ts_ref

! **************************************************************************************************
!> \brief Calculate xc potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_xc_frag(fat_env, v, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_xc_frag"

    real(dp) :: e, virial(3, 3)
    type(pw_env_type), pointer :: pw_env
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(pw_pool_type), pointer :: pw_pool
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, section

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    section => section_vals_get_subs_vals(input, "DFT%XC")
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, auxbas_pw_pool=pw_pool)

    call fat_rho_frag(fat_env, rho_r, i_frag)
    nullify(v, v_tau)
    call xc_vxc_pw_create(v, v_tau, e, rho_r, rho_g, tau, section, pw_pool, .false., virial)
    fat_env%exc_sub(i_frag) = e

  end subroutine fat_v_xc_frag

! **************************************************************************************************
!> \brief Calculate reference xc potential for the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_xc_ref(fat_env, v)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v(:)

    character(len=*), parameter :: routinen = "fat_v_xc_ref"

    real(dp) :: e, virial(3, 3)
    type(pw_env_type), pointer :: pw_env
    type(pw_type), pointer :: rho_g(:), rho_r(:), tau(:), v_tau(:)
    type(pw_pool_type), pointer :: pw_pool
    type(qs_environment_type), pointer :: qs_env
    type(section_vals_type), pointer :: input, section

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, input=input)
    section => section_vals_get_subs_vals(input, "DFT%XC")
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, auxbas_pw_pool=pw_pool)

    call fat_pws_create_ref(fat_env, rho_r, 1)
    call fat_rho_ref(fat_env, rho_r)
    nullify(v, v_tau)
    call xc_vxc_pw_create(v, v_tau, e, rho_r, rho_g, tau, section, pw_pool, .false., virial)
    fat_env%exc_ref = e
    call fat_pws_release_ref(fat_env, rho_r)

  end subroutine fat_v_xc_ref

! **************************************************************************************************
!> \brief Calculate electronic potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_e_frag(fat_env, v, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_e_frag"

    real(dp) :: e
    type(pw_env_type), pointer :: pw_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_pool_type), pointer :: pw_pool
    type(pw_type), pointer :: rho(:), rho_g(:), v_g(:), v_sum
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_rho_frag(fat_env, rho_g, i_frag, .false.)

    allocate(v_g(1))
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_g(1), use_data=complexdata1d, in_space=reciprocalspace)
    call pw_zero(v_g(1))

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, poisson_env=poisson_env)
    call pw_poisson_solve(poisson_env, rho_g(1), e, v_g(1))

    if (is_debug) then
       print *, "DBG| ", e
    end if

    call pw_transfer(v_g(1), v)
    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(v_g(1))
    deallocate(v_g)

  end subroutine fat_v_e_frag

! **************************************************************************************************
!> \brief Calculate nuclear potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_v_n_frag(fat_env, v, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: v
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_v_n_frag"

    real(dp) :: e
    type(pw_env_type), pointer :: pw_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_pool_type), pointer :: pw_pool
    type(pw_type), pointer :: rho(:), rho_g(:), v_g(:), v_sum, rho_tmp
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    allocate(rho_g(1))
    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, rho_core=rho_tmp)
    rho_g(1) = rho_tmp

    allocate(v_g(1))
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_g(1), use_data=complexdata1d, in_space=reciprocalspace)
    call pw_zero(v_g(1))

    qs_env => fat_env%force_env_sub(i_frag)%force_env%qs_env
    call get_qs_env(qs_env, pw_env=pw_env)
    call pw_env_get(pw_env, poisson_env=poisson_env)
    call pw_poisson_solve(poisson_env, rho_g(1), e, v_g(1))

    if (is_debug) then
       print *, "DBG| ", e
    end if

    call pw_transfer(v_g(1), v)
    call fat_env%pw_pool_sub(i_frag)%pool%give_back_pw(v_g(1))
    deallocate(rho_g)
    deallocate(v_g)

  end subroutine fat_v_n_frag

! **************************************************************************************************
!> \brief Calculate non-additive kinetic embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
  subroutine fat_vemb_tsna_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_tsna_frag"

    integer :: i
    type(pw_type), pointer :: v(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call pw_zero(vemb)

    ! total
    nullify(v)
    call fat_v_ts_ref(fat_env, v)
    call pw_axpy(v(1), vemb, allow_noncompatible_grids=.true.)
    call fat_pws_release_ref(fat_env, v)

    ! fragment
    nullify(v)
    call fat_v_ts_frag(fat_env, v, i_frag)
    call pw_axpy(v(1), vemb, -1.0_dp, allow_noncompatible_grids=.true.)
    call fat_pws_release_frag(fat_env, v, i_frag)

  end subroutine fat_vemb_tsna_frag
  
! **************************************************************************************************
!> \brief Calculate non-additive xc embedding potential for one fragment in the freeze-and-thaw
!>        environment
! **************************************************************************************************
  subroutine fat_vemb_xcna_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_xcna_frag"

    integer :: i
    type(pw_type), pointer :: v(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call pw_zero(vemb)

    ! total
    nullify(v)
    call fat_v_xc_ref(fat_env, v)
    call pw_axpy(v(1), vemb, allow_noncompatible_grids=.true.)
    call fat_pws_release_ref(fat_env, v)

    ! fragment
    nullify(v)
    call fat_v_xc_frag(fat_env, v, i_frag)
    call pw_axpy(v(1), vemb, -1.0_dp, allow_noncompatible_grids=.true.)
    call fat_pws_release_frag(fat_env, v, i_frag)

  end subroutine fat_vemb_xcna_frag

! **************************************************************************************************
!> \brief Calculate electronic embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_e_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_e_frag"

    integer :: i
    real(dp) :: e
    type(pw_env_type), pointer :: pw_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_pool_type), pointer :: pw_pool
    type(pw_type), pointer :: rho(:), rho_g(:), v, v_g(:), v_sum
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call pw_zero(vemb)

    allocate(v_sum)
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_sum, use_data=realdata3d, in_space=realspace)
    call pw_zero(v_sum)

    do i = 1, fat_env%n_frag
       if (i == i_frag) cycle

       call fat_pw_create_frag(fat_env, v, i)
       call fat_v_e_frag(fat_env, v, i)
       call pw_axpy(v, v_sum, allow_noncompatible_grids=.true.)
       call fat_env%pw_pool_sub(i)%pool%give_back_pw(v)
       call fat_pw_release_frag(fat_env, v, i)
    end do

    call pw_transfer(v_sum, vemb)
    call fat_pw_release_frag(fat_env, v_sum, i_frag)

  end subroutine fat_vemb_e_frag

! **************************************************************************************************
!> \brief Calculate nuclear embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_n_frag(fat_env, vemb, i_frag)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_vemb_n_frag"

    integer :: i
    real(dp) :: e
    type(pw_env_type), pointer :: pw_env
    type(pw_poisson_type), pointer :: poisson_env
    type(pw_pool_type), pointer :: pw_pool
    type(pw_type), pointer :: rho(:), rho_g(:), v, v_g(:), v_sum
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call pw_zero(vemb)

    allocate(v_sum)
    call fat_env%pw_pool_sub(i_frag)%pool%create_pw(v_sum, use_data=realdata3d, in_space=realspace)
    call pw_zero(v_sum)

    do i = 1, fat_env%n_frag
       if (i == i_frag) cycle

       call fat_pw_create_frag(fat_env, v, i)
       call fat_v_n_frag(fat_env, v, i)
       call pw_axpy(v, v_sum, allow_noncompatible_grids=.true.)
       call fat_pw_release_frag(fat_env, v, i)
    end do

    call pw_transfer(v_sum, vemb)
    call fat_pw_release_frag(fat_env, v_sum, i_frag)

  end subroutine fat_vemb_n_frag
  
! **************************************************************************************************
!> \brief Calculate embedding potential for one fragment in the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_vemb_frag(fat_env, vemb, i_frag, to_file)
    type(fat_env_type), pointer :: fat_env
    type(pw_type), pointer :: vemb
    integer :: i_frag
    logical :: to_file

    character(len=*), parameter :: routinen = "fat_vemb_frag"

    character(len=default_string_length) :: label
    integer :: i, output_unit
    integer, pointer :: stride(:)
    type(particle_list_type), pointer :: particle_list
    type(pw_type), pointer :: rho(:), vemb_tsna, vemb_xcna, vemb_ee, vemb_ne, weights
    type(qs_environment_type), pointer :: qs_env
    type(qs_subsys_type), pointer :: subsys

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_rho_frag(fat_env, rho, i_frag)

    call fat_pw_create_frag(fat_env, vemb_tsna, i_frag)
    call fat_vemb_tsna_frag(fat_env, vemb_tsna, i_frag)
    call pw_axpy(vemb_tsna, vemb)
    call fat_pw_release_frag(fat_env, vemb_tsna, i_frag)

    call fat_pw_create_frag(fat_env, vemb_xcna, i_frag)
    call fat_vemb_xcna_frag(fat_env, vemb_xcna, i_frag)
    call pw_axpy(vemb_xcna, vemb)
    call fat_pw_release_frag(fat_env, vemb_xcna, i_frag)

    call fat_pw_create_frag(fat_env, vemb_ee, i_frag)
    call fat_vemb_e_frag(fat_env, vemb_ee, i_frag)
    call pw_axpy(vemb_ee, vemb)
    call fat_pw_release_frag(fat_env, vemb_ee, i_frag)

    call fat_pw_create_frag(fat_env, vemb_ne, i_frag)
    call fat_vemb_n_frag(fat_env, vemb_ne, i_frag)
    call pw_axpy(vemb_ne, vemb)
    call fat_pw_release_frag(fat_env, vemb_ne, i_frag)

    if (to_file) then
       call integer_to_string(i_frag, label)

       call open_file( &
            "vemb_"//trim(label)//".cube", file_status="REPLACE",&
            & file_action="WRITE", unit_number=output_unit)
       qs_env => fat_env%force_env_ref%qs_env
       call get_qs_env(qs_env, subsys=subsys)
       call qs_subsys_get(subsys, particles=particle_list)
       allocate(stride(3))
       stride = [1, 1, 1]
       call cp_pw_to_cube(vemb, output_unit, "vemb_"//trim(label), particle_list, stride=stride)
       deallocate(stride)

       call open_file( &
            "vemb_"//trim(label)//".dat", file_status="REPLACE", &
            file_action="WRITE", unit_number=output_unit)
       call fat_pw_create_frag(fat_env, weights, i_frag)
       call pw_set(weights, vemb%pw_grid%dr(1) * vemb%pw_grid%dr(2) * vemb%pw_grid%dr(3))
       call pw_to_simple_volumetric(weights, output_unit, [1, 1, 1], vemb)
       call fat_pw_release_frag(fat_env, weights, i_frag)
    end if

  end subroutine fat_vemb_frag

! **************************************************************************************************
!> \brief Calculate the kinetic energy for on fragment of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_ts_frag(fat_env, i_frag) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_ts_frag"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    e = fat_env%ekin_sub(i_frag)
    
  end function fat_e_ts_frag

! **************************************************************************************************
!> \brief Calculate the kinetic energy for the reference system of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_ts_ref(fat_env) result(e)
    type(fat_env_type), pointer :: fat_env
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_ts_ref"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    e = fat_env%ekin_ref
    
  end function fat_e_ts_ref

! **************************************************************************************************
!> \brief Calculate the xc energy for on fragment of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_xc_frag(fat_env, i_frag) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_xc_frag"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    e = fat_env%exc_sub(i_frag)
    
  end function fat_e_xc_frag

! **************************************************************************************************
!> \brief Calculate the xc energy for the reference system of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_xc_ref(fat_env) result(e)
    type(fat_env_type), pointer :: fat_env
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_xc_ref"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    e = fat_env%exc_ref
    
  end function fat_e_xc_ref
  
! **************************************************************************************************
!> \brief Calculate the ee interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_ee(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_ee"

    type(pw_type), pointer :: v, rho, rho_tmp(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_pw_create_frag(fat_env, v, i)
    call fat_v_e_frag(fat_env, v, i)
    call fat_pw_create_frag(fat_env, rho, i)
    call fat_rho_frag(fat_env, rho_tmp, j)

    call pw_transfer(rho_tmp(1), rho)
    e = pw_integral_ab(v, rho)

    call fat_pw_release_frag(fat_env, v, i)
    call fat_pw_release_frag(fat_env, rho, i)
    
    if (is_debug) then
       print *, "DBG| ", i, j, e
    end if
    
  end function fat_eint_ee

! **************************************************************************************************
!> \brief Calculate the ne interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_ne(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_ne"

    type(pw_type), pointer :: v, rho, rho_tmp(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_pw_create_frag(fat_env, v, i)
    call fat_v_n_frag(fat_env, v, i)

    call fat_pw_create_frag(fat_env, rho, i)
    call fat_rho_frag(fat_env, rho_tmp, j)

    call pw_transfer(rho_tmp(1), rho)
    e = pw_integral_ab(v, rho)

    call fat_pw_release_frag(fat_env, v, i)
    call fat_pw_release_frag(fat_env, rho, i)

    if (is_debug) then
       print *, "DBG| ", i, j, e
    end if

  end function fat_eint_ne

! **************************************************************************************************
!> \brief Calculate the nn interaction energy for two fragments of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_eint_nn(fat_env, i, j) result(e)
    type(fat_env_type), pointer :: fat_env
    integer :: i, j
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_eint_nn"

    type(pw_type), pointer :: v, rho, rho_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_pw_create_frag(fat_env, v, i)
    call fat_v_n_frag(fat_env, v, i)

    call fat_pw_create_frag(fat_env, rho, i)
    call fat_rhon_frag(fat_env, rho_tmp, j)

    call pw_transfer(rho_tmp, rho)
    e = pw_integral_ab(v, rho)

    call fat_pw_release_frag(fat_env, v, i)
    call fat_pw_release_frag(fat_env, rho, i)

    if (is_debug) then
       print *, "DBG| ", i, j, e
    end if

  end function fat_eint_nn

! **************************************************************************************************
!> \brief Calculate the total energy of a freeze-and-thaw optimization
! **************************************************************************************************
  function fat_e_tot(fat_env) result(e)
    type(fat_env_type), pointer :: fat_env
    real(dp) :: e

    character(len=*), parameter :: routinen = "fat_e_tot"

    integer :: i, j
    real(dp) :: e_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    e = 0.0_dp

    do i = 1, fat_env%n_frag
       e = e + fat_env%e_sub(i)
    end do

    do i = 1, fat_env%n_frag
       e = e - fat_env%e_emb(i)
    end do

    e = e + fat_e_ts_ref(fat_env)
    do i = 1, fat_env%n_frag
       e = e - fat_e_ts_frag(fat_env, i)
    end do

    e = e + fat_e_xc_ref(fat_env)
    do i = 1, fat_env%n_frag
       e = e - fat_e_xc_frag(fat_env, i)
    end do

    if (fat_env%i_iter > 0) then
       fat_env%eee = 0.0
       do i = 1, fat_env%n_frag
          do j = 1, fat_env%n_frag
             if (i == j) cycle
             e_tmp = fat_eint_ee(fat_env, i, j)
             e = e + .5*e_tmp
             fat_env%eee = fat_env%eee + .5*e_tmp
          end do
       end do

       fat_env%ene = 0.0
       do i = 1, fat_env%n_frag
          do j = 1, fat_env%n_frag
             if (i == j) cycle
             e_tmp = fat_eint_ne(fat_env, i, j)
             e = e + e_tmp
             fat_env%ene = fat_env%ene + e_tmp
          end do
       end do

       fat_env%enn = 0.0
       do i = 1, fat_env%n_frag
          do j = 1, fat_env%n_frag
             if (i == j) cycle
             e_tmp = fat_eint_nn(fat_env, i, j)
             e = e + .5*e_tmp
             fat_env%enn = fat_env%enn + .5*e_tmp
          end do
       end do
    end if

    fat_env%e_tot = e

  end function fat_e_tot

! ************************************************************************************************
!> \brief Print the intermediate results of a freeze-and-thaw iteration
! **************************************************************************************************
  subroutine fat_status(fat_env, i_frag)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_status"

    real(dp) :: e_tot

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    e_tot = fat_e_tot(fat_env)

    write (fat_env%output_unit, fmt=fat_data_fmt) &
         fat_env%i_iter, &
         fat_env%n_iter, &
         i_frag, &
         fat_env%n_frag, &
         fat_env%e_sub(i_frag), &
         fat_env%e_emb(i_frag), &
         fat_env%ekin_ref, &
         fat_env%ekin_sub(i_frag), &
         fat_env%exc_ref, &
         fat_env%exc_sub(i_frag), &
         fat_env%eee, &
         fat_env%ene, &
         fat_env%enn, &
         e_tot
       
  end subroutine fat_status

! **************************************************************************************************
!> \brief Check the freeze-and-thaw optimization for convergence
! **************************************************************************************************
  subroutine fat_iterate_frag(fat_env, i_frag)
    type(fat_env_type), pointer :: fat_env
    integer :: i_frag

    character(len=*), parameter :: routinen = "fat_iterate_frag"

    character :: cmdmsg(100)
    character(default_string_length) :: label
    integer:: cmdstat, exitstat, input_unit, output_unit
    type(dft_control_type), pointer :: dft_control
    type(force_env_type), pointer :: force_env_tmp
    type(pw_type), pointer :: rho(:), vemb, weights
    type(qs_energy_type), pointer :: energy
    type(qs_environment_type), pointer :: qs_env_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen, i_frag
    end if

    force_env_tmp => fat_env%force_env_sub(i_frag)%force_env
    qs_env_tmp => force_env_tmp%qs_env

    call get_qs_env(qs_env_tmp, dft_control=dft_control)
    dft_control%apply_embed_pot = .true.
    call fat_pw_create_frag(fat_env, vemb, i_frag)
    if (fat_env%i_iter > 0) then
       call fat_vemb_frag(fat_env, vemb, i_frag, .false.)
    end if
    call set_qs_env(qs_env_tmp, embed_pot=vemb)

    if (fat_env%extern(i_frag)) then

       call integer_to_string(i_frag, label)
       call open_file( &
            "vemb_"//trim(label)//".dat", file_status="REPLACE", &
            file_action="WRITE", unit_number=output_unit)
       call fat_pw_create_frag(fat_env, weights, i_frag)
       call pw_set(weights, vemb%pw_grid%dr(1) * vemb%pw_grid%dr(2) * vemb%pw_grid%dr(3))
       call pw_to_simple_volumetric(weights, output_unit, [1, 1, 1], vemb)
       call fat_pw_release_frag(fat_env, weights, i_frag)
      
       call execute_command_line("sh extern_"//trim(label)//".sh", .true., exitstat, cmdstat, cmdmsg)
       if (cmdstat /= 0) CPABORT("External command did not execute successfully")

       call get_qs_env(qs_env_tmp, energy=energy)
       call open_file( &
            "extern_"//trim(label)//".e", file_status="OLD", file_action="READ", unit_number=input_unit)
       read (input_unit, *) fat_env%e_sub(i_frag)

       call fat_rho_frag(fat_env, rho, i_frag)
       fat_env%e_emb(i_frag) = pw_integral_ab(vemb, rho(1))

    else
       call qs_calc_energy_force(qs_env_tmp, .true., .false., .false.)
       call get_qs_env(qs_env_tmp, energy=energy)
       fat_env%e_sub(i_frag) = energy%total
       fat_env%e_emb(i_frag) = energy%embed_corr
    end if


    call fat_pw_release_frag(fat_env, vemb, i_frag)

  end subroutine fat_iterate_frag

! **************************************************************************************************
!> \brief Performs one iteration of the freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_iterate(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_iterate"

    integer :: i_frag

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    do i_frag = 1, fat_env%n_frag
       call fat_iterate_frag(fat_env, i_frag)
       call fat_status(fat_env, i_frag)
    end do
    
    fat_env%i_iter = fat_env%i_iter + 1

  end subroutine fat_iterate

! **************************************************************************************************
!> \brief Check the freeze-and-thaw optimization for convergence
! **************************************************************************************************
  function fat_is_converged(fat_env) result(is_converged)
    type(fat_env_type), pointer :: fat_env
    logical :: is_converged

    character(len=*), parameter :: routinen = "fat_is_converged"

    real(dp) :: diff

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    diff = fat_env%e_tot - fat_env%e_old
    is_converged = abs(diff) .le. fat_env%eps
    fat_env%is_converged = is_converged
    
    if (is_debug) then
       print *, "DBG| ", is_converged
       print *, "DBG| ", fat_env%e_tot, fat_env%e_old
       print *, "DBG| diff", diff
    end if

    fat_env%e_old = fat_env%e_tot

  end function

! **************************************************************************************************
!> \brief Check if the freeze-and-thaw optimization has converged or run into the iteration limit
! **************************************************************************************************
  function fat_is_done(fat_env) result(is_done)
    type(fat_env_type), pointer :: fat_env
    logical :: is_done

    character(len=*), parameter :: routinen = "fat_is_done"

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    is_done = fat_env%i_iter > fat_env%n_iter .or. fat_is_converged(fat_env)
  end function

! **************************************************************************************************
!> \brief Allocate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_create(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable, intent(out) :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_create"

    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, particle_set=particles)
    allocate(forces(size(particles), 3))
    forces = 0.0

  end subroutine fat_force_create

! **************************************************************************************************
!> \brief Calculate the self force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_self(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_self"

    integer :: i, i_atom, i_kind, j, k, n
    integer, allocatable :: atom_of_kind(:), kind_of(:)
    type(atomic_kind_type), pointer :: kinds(:)
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(qs_force_type), pointer :: qs_force(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, atomic_kind_set=kinds, force=qs_force)
       call get_atomic_kind_set(kinds, natom=n, kind_of=kind_of, atom_of_kind=atom_of_kind)
       do j = 1, n
          k = k + 1
          i_kind = kind_of(j)
          i_atom = atom_of_kind(j)
          forces(k, :) = qs_force(i_kind)%total(1:3, i_atom)
          if (is_debug) then
             print *, "DBG| k, i_kind, i_atom", k, i_kind, i_atom
          end if
       end do
       deallocate(atom_of_kind, kind_of)
    end do

  end subroutine fat_force_self

! **************************************************************************************************
!> \brief Calculate the embedding force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_emb(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_emb"

    integer :: a, i, j, k, n
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_frag(fat_env, v, i, .true.)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_emb

! **************************************************************************************************
!> \brief Calculate the non-additive kinetic force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_tsna(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_tsna"

    integer :: a, i, j, k, n
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_tsna_frag(fat_env, v, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_tsna

! **************************************************************************************************
!> \brief Calculate the non-additive xc force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_xcna(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_xcna"

    integer :: a, i, j, k, n
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_xcna_frag(fat_env, v, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_xcna

! **************************************************************************************************
!> \brief Calculate the ee force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_ee(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_ee"

    integer :: a, i, j, k, n
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_e_frag(fat_env, v, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_ee

! **************************************************************************************************
!> \brief Calculate the ne force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_ne(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_ee"

    integer :: a, i, j, k, l, n
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: ve, vn, drho, drhon

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, ve, i)
       call fat_vemb_e_frag(fat_env, ve, i)

       call fat_pw_create_frag(fat_env, vn, i)
       call fat_vemb_n_frag(fat_env, vn, i)

       qs_env => fat_env%force_env_sub(i)%force_env%qs_env

       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drho_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(vn, drho)
             call fat_pw_release_frag(fat_env, drho, i)

             call fat_pw_create_frag(fat_env, drhon, i)
             call fat_drhon_frag(fat_env, drhon, i, j, a)
             forces(k, a) = pw_integral_ab(ve, drhon)
             call fat_pw_release_frag(fat_env, drhon, i)
          end do
       end do

       call fat_pw_release_frag(fat_env, ve, i)
       call fat_pw_release_frag(fat_env, vn, i)
    end do

  end subroutine fat_force_ne

! **************************************************************************************************
!> \brief Calculate the nn force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_nn(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_nn"

    integer :: a, i, j, k, n
    type(particle_type), pointer :: particles(:)
    type(qs_environment_type), pointer :: qs_env
    type(pw_type), pointer :: v, drho

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    k = 0
    do i = 1, fat_env%n_frag  ! fragments
       call fat_pw_create_frag(fat_env, v, i)
       call fat_vemb_n_frag(fat_env, v, i)
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, particle_set=particles)
       do j = 1, size(particles)  ! atoms in the fragment
          k = k + 1
          do a = 1, 3  ! xyz
             call fat_pw_create_frag(fat_env, drho, i)
             call fat_drhon_frag(fat_env, drho, i, j, a)
             forces(k, a) = pw_integral_ab(v, drho)
             call fat_pw_release_frag(fat_env, drho, i)
          end do
       end do
       call fat_pw_release_frag(fat_env, v, i)
    end do

  end subroutine fat_force_nn

! ************************************************************************************************
!> \brief Calculate the overlap force for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_overlap(fat_env, forces)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable :: forces(:, :)

    character(len=*), parameter :: routinen = "fat_force_overlap"

    integer :: i, iatom, ikind, natom
    integer, allocatable :: atom_of_kind(:), kind_of(:)
    type(atomic_kind_type), pointer :: atomic_kinds(:)
    type(qs_environment_type), pointer :: qs_env
    type(qs_force_type), pointer :: force(:)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, atomic_kind_set=atomic_kinds, force=force)
    call get_atomic_kind_set( &
         atomic_kinds, atom_of_kind=atom_of_kind, kind_of=kind_of, natom=natom)

    do iatom = 1, natom
       ikind = kind_of(iatom)
       i = atom_of_kind(iatom)
       forces(iatom, :) = force(ikind)%overlap(1:3, i)
    end do

  end subroutine fat_force_overlap

! **************************************************************************************************
!> \brief Print the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force_print(fat_env, forces, label)
    type(fat_env_type), pointer :: fat_env
    real(dp), allocatable, intent(in) :: forces(:, :)
    character(len=*) :: label

    character(len=*), parameter :: routinen = "fat_force_print"

    integer :: i, n, output_unit
    character(len=2) :: symbol
    type(particle_type), pointer :: particles(:)
    type(atomic_kind_type), pointer :: kind
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    qs_env => fat_env%force_env_ref%qs_env
    call get_qs_env(qs_env, particle_set=particles)
    n = size(particles)

    call open_file( &
         label//".xyz", file_status="REPLACE", file_action="WRITE", unit_number=output_unit)
    write (output_unit, '(I0)') size(particles)
    write (output_unit, *)
    do i = 1, n
       kind => particles(i)%atomic_kind
       call get_atomic_kind(kind, element_symbol=symbol)
       write (output_unit, '(A2,3(F26.15," "))') symbol, forces(i, 1), forces(i, 2), forces(i, 3)
    end do
    call close_file(output_unit)

  end subroutine fat_force_print

! **************************************************************************************************
!> \brief Calculate the forces for a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_force(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_force"

    integer :: i, j
    real(dp), allocatable :: &
         & forces(:, :), forces_self(:, :), forces_emb(:, :),&
         & forces_tsna(:, :), forces_xcna(:, :), forces_ee(:, :),&
         & forces_ne(:, :), forces_nn(:, :), forces_overlap(:, :)

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call fat_force_create(fat_env, forces)

    call fat_force_create(fat_env, forces_self)
    call fat_force_self(fat_env, forces_self)
    call fat_force_print(fat_env, forces_self, "force_self")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_self(i, j)
       end do
    end do
    deallocate(forces_self)

    call fat_force_create(fat_env, forces_emb)
    call fat_force_emb(fat_env, forces_emb)
    call fat_force_print(fat_env, forces_emb, "force_emb")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) - forces_emb(i, j)
       end do
    end do
    deallocate(forces_emb)

    call fat_force_create(fat_env, forces_tsna)
    call fat_force_tsna(fat_env, forces_tsna)
    call fat_force_print(fat_env, forces_tsna, "force_tsna")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_tsna(i, j)
       end do
    end do
    deallocate(forces_tsna)

    call fat_force_create(fat_env, forces_xcna)
    call fat_force_xcna(fat_env, forces_xcna)
    call fat_force_print(fat_env, forces_xcna, "force_xcna")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_xcna(i, j)
       end do
    end do
    deallocate(forces_xcna)

    call fat_force_create(fat_env, forces_ee)
    call fat_force_ee(fat_env, forces_ee)
    call fat_force_print(fat_env, forces_ee, "force_ee")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_ee(i, j)
       end do
    end do
    deallocate(forces_ee)

    call fat_force_create(fat_env, forces_ne)
    call fat_force_ne(fat_env, forces_ne)
    call fat_force_print(fat_env, forces_ne, "force_ne")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + 2*forces_ne(i, j)
       end do
    end do
    deallocate(forces_ne)

    call fat_force_create(fat_env, forces_nn)
    call fat_force_nn(fat_env, forces_nn)
    call fat_force_print(fat_env, forces_nn, "force_nn")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_nn(i, j)
       end do
    end do
    deallocate(forces_nn)

    call fat_force_create(fat_env, forces_overlap)
    call fat_force_overlap(fat_env, forces_overlap)
    call fat_force_print(fat_env, forces_overlap, "force_overlap")
    do i = 1, size(forces, 1)
       do j = 1, 3
          forces(i, j) = forces(i, j) + forces_overlap(i, j)
       end do
    end do
    deallocate(forces_overlap)

    call fat_force_print(fat_env, forces, "force")

    deallocate(forces)

  end subroutine fat_force

! **************************************************************************************************
!> \brief Print the end results of a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_results(fat_env, do_forces)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: do_forces

    character(len=*), parameter :: routinen = "fat_results"

    logical :: my_do_forces
    integer :: i, j
    real(dp) :: e_tmp

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if
    
    if (present(do_forces)) then
       my_do_forces = do_forces
    else
       my_do_forces = .false.
    end if

    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| self energy fragment   ",I2,":",F26.15)') i, fat_env%e_sub(i)
    end do

    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| emb energy fragment    ",I2,":",F26.15)') i, fat_env%e_emb(i)
    end do

    write (fat_env%output_unit, '(T2,"FAT| kin energy reference     :",F26.15)') fat_env%ekin_ref
    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| kin energy fragment    ",I2,":",F26.15)') i, fat_env%ekin_sub(i)
    end do

    write (fat_env%output_unit, '(T2,"FAT| xc energy reference      :",F26.15)') fat_env%exc_ref
    do i = 1, fat_env%n_frag
       write (fat_env%output_unit, '(T2,"FAT| xc energy fragment     ",I2,":",F26.15)') i, fat_env%exc_sub(i)
    end do

    do i = 1, fat_env%n_frag
       do j = 1, fat_env%n_frag
          if (i == j) cycle
          e_tmp = fat_eint_ee(fat_env, i, j)
          write (fat_env%output_unit, '(T2,"FAT| ee energy fragments ",I2," ",I2,":",F26.15)') i, j, e_tmp
       end do
    end do

    do i = 1, fat_env%n_frag
       do j = 1, fat_env%n_frag
          if (i == j) cycle
          e_tmp = fat_eint_ne(fat_env, i, j)
          write (fat_env%output_unit, '(T2,"FAT| ne energy fragments ",I2," ",I2,":",F26.15)') i, j, e_tmp
       end do
    end do

    fat_env%enn = 0.0
    do i = 1, fat_env%n_frag
       do j = 1, fat_env%n_frag
          if (i == j) cycle
          e_tmp = fat_eint_nn(fat_env, i, j)
          write (fat_env%output_unit, '(T2,"FAT| nn energy fragments ",I2," ",I2,":",F26.15)') i, j, e_tmp
       end do
    end do
       
    write (fat_env%output_unit, '(/,T2,"ENERGY| Total FORCE_EVAL ( ",A," ) energy [a.u.]: ",T55,F26.15)') &
       "FAT", fat_env%e_tot

    if (my_do_forces) then
       call fat_force(fat_env)
    end if

  end subroutine fat_results

! **************************************************************************************************
!> \brief Cleans up a freeze-and-thaw optimization
! **************************************************************************************************
  subroutine fat_finalize(fat_env)
    type(fat_env_type), pointer :: fat_env

    character(len=*), parameter :: routinen = "fat_finalize"

    integer :: i
    type(dft_control_type), pointer :: dft_control
    type(pw_type), pointer :: vemb
    type(qs_environment_type), pointer :: qs_env

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    do i = 1, fat_env%n_frag
       qs_env => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env, dft_control=dft_control)
       dft_control%apply_embed_pot = .false.
    end do

  end subroutine fat_finalize

! **************************************************************************************************
!> \brief Initialize the freeze-and-thaw environment
! **************************************************************************************************
  subroutine fat_init(fat_env, force_env, ref_subsys_number)
    type(fat_env_type), pointer   :: fat_env
    type(force_env_type), pointer :: force_env
    integer                       :: ref_subsys_number

    character(len=*), parameter :: routinen = "fat_init"

    integer :: handle, i, i_frag, j, n_iter, n_frag, cluster_subsys_num
    real(dp) :: eps
    type(cp_logger_type), pointer :: logger
    type(dft_control_type), pointer :: dft_control
    type(force_env_type), pointer :: force_env_tmp
    type(pw_env_type), pointer :: pw_env_tmp
    type(pw_pool_type), pointer :: pw_pool_tmp
    type(pw_type), pointer :: rho_r_tmp(:), vemb
    type(qs_environment_type), pointer :: qs_env_tmp
    type(qs_energy_type), pointer :: energy_tmp
    type(qs_rho_type), pointer :: rho_tmp
    type(section_vals_type), pointer :: input

    if (is_debug) then
       print *, "DBG| ", modulen//"::"//routinen
    end if

    call timeset(routinen, handle)

    nullify(qs_env_tmp)

    ! this is a limitation of the EMBED method
    CPASSERT(ref_subsys_number == 1)

    fat_env%force_env_ref => force_env%sub_force_env(ref_subsys_number)%force_env
    fat_env%n_frag = size(force_env%sub_force_env) - 1
    fat_env%i_iter = 0
    fat_env%output_unit = cp_logger_get_default_io_unit()

    ! force environments
    allocate(fat_env%force_env_sub(fat_env%n_frag))
    j = 1
    do i = 1, size(force_env%sub_force_env)
       if (i == ref_subsys_number) cycle
       fat_env%force_env_sub(j)%force_env => force_env%sub_force_env(i)%force_env
       j = j + 1
    end do

    ! pw pools
    force_env_tmp => fat_env%force_env_ref
    qs_env_tmp => force_env_tmp%qs_env
    call get_qs_env(qs_env_tmp, pw_env=pw_env_tmp)
    call pw_env_get(pw_env_tmp, auxbas_pw_pool=fat_env%pw_pool_ref)
    allocate(fat_env%pw_pool_sub(fat_env%n_frag))
    do i = 1, fat_env%n_frag
       force_env_tmp => fat_env%force_env_sub(i)%force_env
       qs_env_tmp => fat_env%force_env_sub(i)%force_env%qs_env
       call get_qs_env(qs_env_tmp, pw_env=pw_env_tmp)
       call pw_env_get(pw_env_tmp, auxbas_pw_pool=fat_env%pw_pool_sub(i)%pool)
    end do

    ! mark the extern fragments, if specified
    allocate(fat_env%extern(fat_env%n_frag))
    do i = 1, fat_env%n_frag
       force_env_tmp => fat_env%force_env_sub(i)%force_env
       qs_env_tmp => force_env_tmp%qs_env
       call get_qs_env(qs_env_tmp, input=input)
       call section_vals_val_get(input, "DFT%QS%FAT_EXTERN", l_val=fat_env%extern(i))
       if (is_debug) then
          print *, "DBG| fat_env%extern(i)", i, fat_env%extern(i)
       end if
    end do

    do i = 1, fat_env%n_frag
       force_env_tmp => fat_env%force_env_sub(i)%force_env
       qs_env_tmp => force_env_tmp%qs_env
       call qs_calc_energy_force(qs_env_tmp, .true., .false., .false.)
    end do

    ! energies
    allocate(fat_env%e_sub(fat_env%n_frag))
    fat_env%e_sub = 0.0
    allocate(fat_env%e_emb(fat_env%n_frag))
    fat_env%e_emb = 0.0

    fat_env%ekin_ref = 0.0
    allocate(fat_env%ekin_sub(fat_env%n_frag))
    fat_env%ekin_sub = 0.0

    fat_env%exc_ref = 0.0
    allocate(fat_env%exc_sub(fat_env%n_frag))
    fat_env%exc_sub = 0.0

    fat_env%eee = 0.0
    fat_env%ene = 0.0
    fat_env%enn = 0.0

    fat_env%e_old = 0.0
    do i = 1, fat_env%n_frag
       call fat_iterate_frag(fat_env, i)
    end do

    ! read input
    call get_qs_env(fat_env%force_env_ref%qs_env, input=input)
    call section_vals_val_get(input, "DFT%QS%OPT_EMBED%N_ITER", i_val=n_iter)
    fat_env%n_iter = n_iter
    call section_vals_val_get(input, "DFT%QS%OPT_EMBED%DENS_CONV_MAX", r_val=eps)
    fat_env%eps = eps

    write (fat_env%output_unit, fmt=fat_title_fmt) &
         "i_iter", &
         "n_iter", &
         "i_frag", &
         "n_frag", &
         "self", &
         "emb", &
         "kin ref", &
         "kin frag", &
         "xc ref", &
         "xc frag", &
         "ee", &
         "ne", &
         "nn", &
         "tot"

    ! run the isolated fragments
    call fat_iterate(fat_env)

    call timestop(handle)

  end subroutine fat_init

end module fat_methods

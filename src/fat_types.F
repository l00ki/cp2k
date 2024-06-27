module fat_types
  use force_env_types
  use kinds
  use pw_env_types
  use pw_pool_types
  use pw_types
  use qs_environment_types

#include "./base/base_uses.f90"

  implicit none

  private

  public :: &
       fat_env_type, &
       fat_env_get, &
       fat_env_set, &
       fat_env_create, &
       fat_env_release

  logical, private, parameter :: is_debug = .true.
  character(len=*), private, parameter :: modulen = "fat_types"

  type fat_env_type
     type(force_env_type), pointer   :: force_env_ref
     type(force_env_p_type), pointer :: force_env_sub(:)
     integer                         :: i_iter, n_iter, n_frag
     logical                         :: is_converged
     real(dp)                        :: eps, e_old, ekin_ref, exc_ref, e_tot, eee, ene, enn
     real(dp), allocatable           :: e_emb(:), e_sub(:), ekin_sub(:), exc_sub(:)
     type(pw_type), pointer          :: rho_r_sub(:), vemb(:)
     type(pw_pool_type), pointer     :: pw_pool_ref
     type(pw_pool_p_type), pointer   :: pw_pool_sub(:)
     integer                         :: output_unit
  end type fat_env_type
  
contains

! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_get(fat_env, is_converged)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: is_converged

    if (present(is_converged)) is_converged = fat_env%is_converged

  end subroutine fat_env_get
  
! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_set(fat_env, is_converged)
    type(fat_env_type), pointer :: fat_env
    logical, optional :: is_converged

    if (present(is_converged)) fat_env%is_converged = is_converged

  end subroutine fat_env_set
  
! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_create(fat_env)
    type(fat_env_type), pointer :: fat_env

    nullify(fat_env)
    allocate(fat_env)

  end subroutine fat_env_create
  
! **************************************************************************************************
! **************************************************************************************************
  subroutine fat_env_release(fat_env)
    type(fat_env_type), pointer :: fat_env

    deallocate(fat_env)

  end subroutine fat_env_release

end module fat_types

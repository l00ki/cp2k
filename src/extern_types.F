! **************************************************************************************************
!> \brief Types for external calculations
!> \author Lukas Schreder
! **************************************************************************************************
module extern_types
  use cell_types
  use cp_subsys_types
  use distribution_1d_types
  use input_section_types
  use kinds
  use message_passing
#include "./base/base_uses.f90"

  implicit none

  private

  logical, parameter :: debug_this_module = .true.
  character(len=*), parameter :: moduleN = "extern_types"

  public :: extern_env_type
  
  public :: extern_env_create, &
            extern_env_release, &
            get_extern_env, &
            set_extern_env

! **************************************************************************************************
!> \brief Contains the data of an external calculation.
!> \author Lukas Schreder
! **************************************************************************************************
  type extern_env_type
    real(dp) :: extern_energy = 0.0_dp
    type(cell_type), pointer :: cell_ref => NULL()
    type(section_vals_type), pointer :: input => NULL()
    type(mp_para_env_type), pointer :: para_env
    type(cp_subsys_type), pointer :: subsys => NULL()
  end type extern_env_type

contains

! **************************************************************************************************
!> \brief Gets something from the extern environment.
!> \param extern_env ...
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine get_extern_env( &
       extern_env, cell, cell_ref, extern_energy, input, local_molecules, local_particles, subsys)
    type(extern_env_type), pointer :: extern_env
    type(cell_type), pointer, optional :: cell, cell_ref
    real(dp), optional :: extern_energy
    type(section_vals_type), pointer, optional :: input
    type(distribution_1d_type), pointer, optional :: local_molecules, local_particles
    type(cp_subsys_type), pointer, optional :: subsys

    if (present(cell)) call cp_subsys_get(extern_env%subsys, cell=cell)
    if (present(cell_ref)) cell_ref => extern_env%cell_ref
    if (present(extern_energy)) extern_energy = extern_env%extern_energy
    if (present(input)) input => extern_env%input
    if (present(subsys)) subsys => extern_env%subsys
    call cp_subsys_get(extern_env%subsys, &
         local_molecules=local_molecules, &
         local_particles=local_particles)

  end subroutine get_extern_env

! **************************************************************************************************
!> \brief Sets something in the extern environment.
!> \param extern_env ...
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine set_extern_env( &
       extern_env, cell, cell_ref, extern_energy, input, local_molecules, local_particles, subsys)
    type(extern_env_type), pointer :: extern_env
    type(cell_type), pointer, optional :: cell, cell_ref
    real(dp), optional :: extern_energy
    type(section_vals_type), pointer, optional :: input
    type(distribution_1d_type), pointer, optional :: local_molecules, local_particles
    type(cp_subsys_type), pointer, optional :: subsys

    if (present(cell)) call cp_subsys_set(extern_env%subsys, cell=cell)
    if (present(cell_ref)) extern_env%cell_ref => cell_ref
    if (present(extern_energy)) extern_env%extern_energy = extern_energy
    if (present(input)) extern_env%input => input
    if (present(subsys)) extern_env%subsys => subsys
    if (present(local_molecules)) call cp_subsys_set( &
         extern_env%subsys, local_molecules=local_molecules)
    if (present(local_particles)) call cp_subsys_set( &
         extern_env%subsys, local_particles=local_particles)

  end subroutine set_extern_env

! **************************************************************************************************
!> \brief Creates the extern environment.
!> \param extern_env ...
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine extern_env_create(extern_env, para_env)
    type(extern_env_type), pointer, intent(out) :: extern_env
    type(mp_para_env_type), intent(in), target :: para_env

    CPASSERT(.not. associated(extern_env))

    allocate(extern_env)
    extern_env%para_env => para_env
    call extern_env%para_env%retain()

  end subroutine extern_env_create

! **************************************************************************************************
!> \brief Releases the extern environment.
!> \param extern_env ...
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine extern_env_release(extern_env)
    type(extern_env_type), pointer :: extern_env

    deallocate(extern_env)

  end subroutine extern_env_release

end module extern_types

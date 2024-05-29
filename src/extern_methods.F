! **************************************************************************************************
!> \brief Methods for external calculations
!> \author Lukas Schreder
! **************************************************************************************************
module extern_methods
  use atomic_kind_types
  use cell_methods
  use cell_types
  use cp_log_handling
  use cp_subsys_methods
  use cp_subsys_types
  use distribution_1d_types
  use distribution_methods
  use extern_types
  use force_env_types
  use input_section_types
  use kinds
  use message_passing
  use molecule_kind_types
  use molecule_types
  use particle_types
#include "./base/base_uses.f90"

  implicit none

  private

  logical, parameter :: debug_this_module = .true.
  character(len=*), parameter :: moduleN = "extern_methods"

  public :: extern_calc_energy_force, &
            extern_init

contains

! **************************************************************************************************
!> \brief Calculates energy and forces with an external calculation
!> \param extern_env ...
!> \param calc_force ...
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine extern_calc_energy_force(force_env, calc_force)
    type(force_env_type), pointer :: force_env
    logical :: calc_force

    character(len=*), parameter :: routineN = "extern_calc_energy_force"
    
    character :: cmdmsg(100)
    integer :: cmdstat, exitstat, handle, output_unit

    call timeset(routineN, handle)

    call execute_command_line("pymolcas -i molcas.inp", .true., exitstat, cmdstat, cmdmsg)
    if (cmdstat /= 0) CPABORT("External command did not execute")

    output_unit = cp_logger_get_default_io_unit()
    write (output_unit, "(T2,A,I3)") "EXTERN| External command returned status code", exitstat
    
    call timestop(handle)

  end subroutine extern_calc_energy_force

! **************************************************************************************************
!> \brief Initialises an external calculation.
!> \param extern_env ...
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine extern_init(extern_env, root_section, para_env, force_env_section, use_motion_section)
    type(extern_env_type), pointer, intent(inout) :: extern_env
    type(section_vals_type), pointer :: root_section
    type(mp_para_env_type), pointer :: para_env
    type(section_vals_type), pointer :: force_env_section
    logical, intent(in) :: use_motion_section

    character(len=*), parameter :: routineN = "extern_init"

    integer :: handle
    logical :: use_ref_cell
    real(dp) :: abc(3)
    type(atomic_kind_type), pointer :: atomic_kind_set(:)
    type(cell_type), pointer :: cell, cell_ref
    type(cp_subsys_type), pointer :: subsys
    type(distribution_1d_type), pointer :: local_molecules, local_particles
    type(molecule_kind_type), pointer :: molecule_kind_set(:)
    type(molecule_type), pointer :: molecule_set(:)
    type(particle_type), pointer :: particle_set(:)
    type(section_vals_type), pointer :: cell_section, subsys_section

    call timeset(routineN, handle)

    subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
    cell_section => section_vals_get_subs_vals(subsys_section, "CELL")

    call extern_env_create(extern_env, para_env)
    call set_extern_env(extern_env, input=force_env_section)
    nullify(subsys)
    call cp_subsys_create(subsys, para_env, root_section, force_env_section=force_env_section, &
         use_motion_section=use_motion_section)

    nullify(cell, cell_ref)
    call read_cell( &
         cell, cell_ref, use_ref_cell=use_ref_cell, cell_section=cell_section, para_env=para_env)
    call get_cell(cell, abc=abc)

    call write_cell(cell, subsys_section)
    call write_cell(cell_ref, subsys_section)

    molecule_kind_set => subsys%molecule_kinds%els
    call write_molecule_kind_set(molecule_kind_set, subsys_section)

    atomic_kind_set => subsys%atomic_kinds%els
    molecule_set => subsys%molecules%els
    particle_set => subsys%particles%els
    call distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, &
                                 local_particles=local_particles, &
                                 molecule_kind_set=molecule_kind_set, &
                                 molecule_set=molecule_set, &
                                 local_molecules=local_molecules, &
                                 force_env_section=force_env_section)
    
    call cp_subsys_set(subsys, cell=cell)

    call set_extern_env(extern_env, subsys=subsys)
    call set_extern_env( &
         extern_env, cell_ref=cell_ref, local_molecules=local_molecules, &
         local_particles=local_particles)

    call distribution_1d_release(local_molecules)
    call distribution_1d_release(local_particles)

    call cell_release(cell)
    call cell_release(cell_ref)

    call timestop(handle)

  end subroutine extern_init

end module extern_methods

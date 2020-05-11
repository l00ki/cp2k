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
   use cp_control_types,      only: dft_control_type
   use cp_log_handling,       only: cp_logger_type, &
                                    cp_get_default_logger
   use cp_output_handling,    only: cp_print_key_unit_nr, &
                                    cp_print_key_finished_output
   use cp_para_types,         only: cp_para_env_type
   use input_section_types,   only: section_vals_type, &
                                    section_vals_val_get
   use kinds,                 only: dp
   use particle_types,        only: particle_type
   use qs_environment_types,  only: qs_environment_type, &
                                    get_qs_env
   use qs_loc_types,          only: get_qs_loc_env, &
                                    qs_loc_env_new_type
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
! **************************************************************************************************
   subroutine calculate_wannier_moments(qs_env, qs_loc_env, loc_print_key)

     type(qs_environment_type), pointer   :: qs_env
     type(qs_loc_env_new_type), pointer   :: qs_loc_env
     type(section_vals_type), pointer     :: loc_print_key
     
     character(len=*), parameter :: routineN = "calculate_wannier_moments", &
        routineP = moduleN//":"//routineN

     integer                              :: handle, iion, ispin, istate, ndims, nions, nspins, &
                                             nstates, order, unit_nr
     real(dp), pointer                    :: rcentres(:, :), rions(:, :)
     type(dft_control_type), pointer      :: dft_control
     type(wannier_centres_type), pointer  :: wcs(:)
     type(cp_logger_type), pointer        :: logger
     type(cell_type), pointer             :: cell
     type(cp_para_env_type), pointer      :: para_env
     type(particle_type), pointer         :: particle_list(:)

     call timeset(routineN, handle)

     logger => cp_get_default_logger()
     unit_nr = cp_print_key_unit_nr(logger, loc_print_key, "WANNIER_MOMENTS", extension=".WanMom", &
                                    middle_name="WANNIER_MOMENTS")

     call get_qs_env(qs_env, cell=cell, dft_control=dft_control, wanniercentres=wcs, particle_set=particle_list)
     nspins = dft_control%nspins
     nstates = size(wcs(1)%centres, 2)
     ndims = size(wcs(1)%centres, 1)
     nions = size(particle_list)

     call get_qs_loc_env(qs_loc_env, para_env=para_env)

     call section_vals_val_get(loc_print_key, "WANNIER_MOMENTS%ORDER", i_val=order)

     allocate(rions(ndims, nions))
     do iion = 1, nions
        rions(:, iion) = particle_list(iion)%r
     end do
     
     do ispin = 1, nspins
        ! wrap wcs to cell
        call wrap_wcs_to_cell(wcs(ispin)%centres, cell, rcentres)
        ! map wcs to ions
        call map_wcs_to_ions(rcentres, rions, unit_nr)
     end do

     call cp_print_key_finished_output(unit_nr, logger, loc_print_key, "WANNIER_MOMENTS")

     deallocate(rcentres, rions)
     
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
   subroutine map_wcs_to_ions(rcentres, rions, unit_nr)

     real(dp), pointer, intent(in)      :: rcentres(:, :)
     real(dp), pointer, intent(in)      :: rions(:, :)
     integer, intent(in)                :: unit_nr
     
     integer                            :: icentre, idim, iion, ncentres, ndims, nions
     integer, allocatable               :: wcs_ions_map(:)
     real(dp), pointer                  :: r(:), d(:)

     ndims = size(rcentres, 1)
     ncentres = size(rcentres, 2)
     nions = size(rions, 2)

     allocate(wcs_ions_map(ncentres), r(ndims), d(nions))

     do icentre = 1, ncentres
        !print *, "centre", icentre, rcentres(:, icentre)
        do iion = 1, nions
           !print *, "ion", iion, rions(:, iion)
           r = rcentres(:, icentre) - rions(:, iion)
           !print *, "r", r
           d(iion) = norm2(r)
        end do
        wcs_ions_map(icentre) = minloc(d, 1)
     end do
     
     if (unit_nr > 0) then
        write(unit=unit_nr, fmt="(A)") "# wannier centre index, mapped ion index"
        do icentre = 1, ncentres
           write(unit=unit_nr, fmt="(T2,I7,A,I7)") icentre, " ", wcs_ions_map(icentre)
        end do
     end if

     deallocate(wcs_ions_map, r, d)
     
   end subroutine map_wcs_to_ions
   
end module wannier_moments

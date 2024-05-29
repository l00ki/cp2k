! **************************************************************************************************
!> \brief Subroutines to create the EXTERN section and its subsection
!> \author Lukas Schreder
! **************************************************************************************************
module input_cp2k_extern
  use input_constants
  use input_keyword_types
  use input_section_types
  use string_utilities
#include "./base/base_uses.f90"

  implicit none

  private

  character(len=*), parameter, private :: moduleN = "input_cp2k_extern"

  public :: create_extern_section

contains

! **************************************************************************************************
!> \brief Create the EXTERN section.
!> \param section the section to create
!> \author Lukas Schreder
! **************************************************************************************************
  subroutine create_extern_section(section)
    type(section_type), pointer :: section

    type(keyword_type), pointer :: keyword

    CPASSERT(.not. associated(section))

    call section_create( &
         section, __LOCATION__, name="EXTERN", &
         description="Input for external calculations", repeats=.false.)

    ! keywords
    nullify(keyword)

    call keyword_create( &
         keyword, __LOCATION__, name="COMMAND", &
         description="The command to run.", &
         usage="COMMAND sh extern.sh", &
         default_c_val="sh extern.sh")
    call section_add_keyword(section, keyword)
    call keyword_release(keyword)

    call keyword_create( &
         keyword, __LOCATION__, name="ENERGY_FILENAME", &
         description="The file containing the energy in Hartree.", &
         usage="ENERGY_FILENAME extern.energy", &
         default_c_val="extern.energy")
    call section_add_keyword(section, keyword)
    call keyword_release(keyword)

    call keyword_create( &
         keyword, __LOCATION__, name="E_DENSITY_FILENAME", &
         description="The file containing the electron density. Must be on the same grid.", &
         usage="E_DENSITY_FILENAME extern.cube", &
         default_c_val="extern.cube")
    call section_add_keyword(section, keyword)
    call keyword_release(keyword)
    ! subsections

  end subroutine create_extern_section

end module input_cp2k_extern

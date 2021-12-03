module Connection_module


  implicit none

private
  type, public :: connection_set_type
    Integer, pointer :: id_up(:)      ! list of ids of upwind cells
    Integer, pointer :: id_dn(:)      ! list of ids of downwind cells
    Real, pointer :: dist(:)    ! list of distance vectors
    Real, pointer :: area(:)      ! list of areas of faces normal to distance vectors

  end type connection_set_type

public :: get_natveg_column_id 
 
contains

! ************************************************************************** !
function get_natveg_column_id(id) result(id_out)

  implicit none
  integer, intent(in) :: id
  integer, id_out
  id_out=id ! for 2D transect only 


end function get_natveg_column_id



end module Connection_module

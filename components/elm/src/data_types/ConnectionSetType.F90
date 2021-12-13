module ConnectionSetType

use shr_kind_mod   , only : r8 => shr_kind_r8
use shr_infnan_mod  , only : isnan => shr_infnan_isnan,nan => shr_infnan_nan, assignment(=)
use decompMod      , only : bounds_type
implicit none
save
public
  
type, public :: connection_set_type
    Integer, pointer :: nconn  => null()         ! number of connections
    Integer, pointer :: grid_id_up(:)  => null()    ! list of ids of upwind cells
    Integer, pointer :: grid_id_dn(:)  => null()    ! list of ids of downwind cells
    Real(r8), pointer :: dist(:)   => null()      ! list of distance vectors
    Real(r8), pointer :: area(:)   => null()      ! list of areas of faces normal to distance vectors
contains
    procedure, public :: Init => col_connect_init
end type connection_set_type

public :: get_natveg_column_id 
 
type (connection_set_type), public, target :: conn   ! connection type

contains

! ************************************************************************** !
subroutine col_connect_init(this, bounds)

   type(bounds_type), intent(in)    :: bounds
   class(connection_set_type)       :: this 
   Integer                          :: n, iconn,begc,endc,begg, endg  
   Integer                          :: dx, dz, g
   begc = bounds%begc;  endc = bounds%endc
   begg = bounds%begg;  endg = bounds%endg
   n = endg-begg
   ! number of connections in each layer
   allocate(this%nconn)             ;  this%nconn = n                   
   allocate(this%grid_id_up(n)) ;  this%grid_id_up(:) = 0
   allocate(this%grid_id_dn(n)) ;  this%grid_id_dn(:) = 0
   allocate(this%area(n))       ;  this%area(:) = 0
   allocate(this%dist(n))       ;  this%dist(:) = 0
   dx = 1000
   dz = 1000
   do iconn = 1,n
     g = begg+iconn-1
     this%grid_id_up(iconn) = g    !  Step-2: Eventually will need to read from surface dataset
     this%grid_id_dn(iconn) = g+1  !  There is already some code that we will be able to
     this%area(iconn)       = 1000      !  use to fill this data structure
     this%dist(iconn)       = 1000      ! triangle law
   enddo
end subroutine col_connect_init

function get_natveg_column_id(id) result(id_out)

  implicit none
  integer, intent(in) :: id
  integer :: id_out
  id_out=id ! for 2D transect only 


end function get_natveg_column_id
end module ConnectionSetType

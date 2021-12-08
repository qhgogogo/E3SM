module Connection_module

use shr_kind_mod   , only : r8 => shr_kind_r8
use decompMod      , only : bounds_type
implicit none
save
public
  
type, public :: connection_set_type
    Integer, pointer :: grid_id_up(:)   => null()      ! list of ids of upwind cells
    Integer, pointer :: grid_id_dn(:)   => null()      ! list of ids of downwind cells
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
   Integer                          :: nconn, iconn,begc,endc,begg, endg  
   Integer                          :: dx, dz
   begc = bounds%begc;  endc = bounds%endc!
   begg = bounds%begg;  endg = bounds%endg
   nconn = endc-begc                                           ! number of connections in each layer
   
   allocate(this%grid_id_up(nconn)) ;  this%grid_id_up(:) = nan
   allocate(this%grid_id_dn(nconn)) ;  this%grid_id_dn(:) = nan
   allocate(this%area(nconn))       ;  this%area(:) = nan
   allocate(this%dist(nconn))       ;  this%dist(:) = nan
   dx = 1000  ! temporary
   dz = 1000  ! temporary
   do iconn = 1,nconn
     g = bounds%begg(iconn)
     this%grid_id_up(iconn) = g    !  Step-2: Eventually will need to read from surface dataset
     this%grid_id_dn(iconn) = g+1  !  There is already some code that we will be able to
     this%area(iconn)       = 1000     !  use to fill this data structure, temporary value now
     this%dist(iconn)       = 1000     ! triangle law, temporary value now
   enddo
end subroutine col_connect

  function get_natveg_column_id(id) result(id_out)

    implicit none
    integer, intent(in) :: id
    integer ::id_out
    id_out=id ! for 2D transect only 
  
  end function get_natveg_column_id

end module Connection_module



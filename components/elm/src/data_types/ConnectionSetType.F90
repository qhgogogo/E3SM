module Connection_module

use shr_kind_mod   , only : r8 => shr_kind_r8
implicit none
save
public
  
type, public :: connection_set_type
    Integer, pointer :: id_up(:)      ! list of ids of upwind cells
    Integer, pointer :: id_dn(:)      ! list of ids of downwind cells
    Real(r8), pointer :: dist(:)    ! list of distance vectors
    Real, pointer :: area(:)      ! list of areas of faces normal to distance vectors
contains
    procedure, public :: connect => col_connect
end type connection_set_type

public :: get_natveg_column_id 
 
type (connection_set_type), public, target :: conn   ! connection type

contains

! ************************************************************************** !
subroutine col_connect(this, begc, endc)

   Integer, intent(in) begc, endc
   Integer, nconn
   class(connection_set_type) :: this  
   allocate(this%grid_id_up(nconn)) ;  this%grid_id_up(:) = nan
   allocate(this%grid_id_dn(nconn)) ;  this%grid_id_dn(:) = nan
   allocate(this%area(nconn))       ;  this%area(:) = nan
   allocate(this%dist(nconn))       ;  this%dist(:) = nan
   do iconn = 1,nconn
     g = begg(iconn)
     conn%grid_id_up(iconn) = g    !... Step-2: Eventually will need to read from surface dataset
     conn%grid_id_dn(iconn) = g+1  !...         There is already some code that we will be able to
                                !            use to fill this data structure
   enddo
end subroutine col_connect

function get_natveg_column_id(id) result(id_out)

  implicit none
  integer, intent(in) :: id
  integer, id_out
  id_out=id ! for 2D transect only 


end function get_natveg_column_id
end module Connection_module



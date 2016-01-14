subroutine plotsamples(npt,x,y)
use precision
implicit none
integer :: npt
real, allocatable, dimension(:) :: xp,yp
real(double), dimension(:) :: x,y

allocate(xp(npt),yp(npt))
xp=real(x(1:npt))
yp=real(y(1:npt))

call pgsci(2)
call pgline(npt,xp,yp)

deallocate(xp,yp)

return
end

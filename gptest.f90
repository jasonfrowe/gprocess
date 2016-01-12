program gptest
use precision
implicit none
integer :: iargc,nmax,npt
real(double), allocatable, dimension(:) :: x,y,yerr
character(80) :: filename

interface
   subroutine getdata(filename,npt,nmax,x,y,yerr)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), dimension(:), intent(inout) :: x,y,yerr
      character(80), intent(in) :: filename
   end subroutine getdata
end interface

if(iargc().lt.1)then
   write(0,*) "Usage: gptest filename"
   stop
endif

!read in data (2 columns)
call getarg(1,filename)

nmax=1000 !initial guess for number of datapoints.
allocate(x(nmax),y(nmax),yerr(nmax))
call getdata(filename,npt,nmax,x,y,yerr)
write(0,*) "Number of points read: ",npt

end program gptest

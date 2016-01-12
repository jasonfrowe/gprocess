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
interface
   subroutine plotdata(npt,x,y,yerr)
      use precision
      implicit none
      integer, intent(in) :: npt
      real(double), dimension(:), intent(in) :: x,y,yerr
   end subroutine plotdata
end interface

if(iargc().lt.1)then
   write(0,*) "Usage: gptest filename"
   stop
endif

!read in data (3 columns)
call getarg(1,filename)

nmax=1000 !initial guess for number of datapoints.
allocate(x(nmax),y(nmax),yerr(nmax))
call getdata(filename,npt,nmax,x,y,yerr)
write(0,*) "Number of points read: ",npt

!open PGPLOT device
call pgopen('?')
call PGPAP (8.0 ,1.0) !use a square 8" across
!call pgsubp(1,4)
call pgpage()

!plot the data
call plotdata(npt,x,y,yerr)

!fit a straight line
call fitline(npt,x,y,yerr)

call pgclos()

end program gptest

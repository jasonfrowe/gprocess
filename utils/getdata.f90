subroutine getdata(filename,npt,nmax,x,y,yerr)
use precision
implicit none
integer :: npt,nmax,nunit,filestatus
real(double), dimension(:) :: x,y,yerr
character(80) :: filename

nunit=10
open(unit=nunit,file=filename,iostat=filestatus,status='old')
if(filestatus>0)then !trap missing file errors
   write(0,*) "Cannot open ",filename
   stop
endif


end subroutine getdata

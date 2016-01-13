program gptest
use precision
implicit none
integer :: iargc,nmax,npt,npars,info
real(double) :: chisq,bpix
real, allocatable, dimension(:) :: bb !contains dimensions of plot
real(double), allocatable, dimension(:) :: x,y,yerr,ans,eans,pars
real(double), allocatable, dimension(:,:) :: Kernel,Kfac
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
   subroutine plotdata(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(in) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(in) :: x,y,yerr
   end subroutine plotdata
end interface
interface
   subroutine fitline(npt,x,y,yerr,ans,eans,chisq)
      use precision
      implicit none
      integer, intent(in) :: npt
      real(double), intent(inout) :: chisq
      real(double), dimension(:), intent(in) :: x,y,yerr
      real(double), dimension(:), intent(inout) :: ans,eans
   end subroutine fitline
end interface
interface
   subroutine plotline(bb,ans,eans)
      use precision
      implicit none
      real, dimension(:), intent(in) :: bb
      real(double), dimension(:), intent(in) :: ans,eans
   end subroutine plotline
end interface
interface
   subroutine makekernel(Kernel,npt,x,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(in) :: npt,npars
      real(double), dimension(:), intent(in) :: x,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface
interface
   subroutine displaykernel(nx,ny,Kernel,bpix)
      use precision
      implicit none
      integer, intent(in) :: nx,ny
      real(double), intent(in) :: bpix
      real(double), dimension(:,:), intent(in) :: Kernel
   end subroutine displaykernel
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
call pgslw(3) !thicker lines

!plot the data
allocate(bb(4)) !contains plot boundaries
call plotdata(npt,x,y,yerr,bb)

!fit a straight line
allocate(ans(2),eans(2)) !contains fit for a straight line
call fitline(npt,x,y,yerr,ans,eans,chisq)
call plotline(bb,ans,eans)

!repeat but increase error via chi-squared
call pgpage()
yerr=yerr*sqrt(chisq/dble(npt))

!plot the data
call plotdata(npt,x,y,yerr,bb)

!fit a straight line
call fitline(npt,x,y,yerr,ans,eans,chisq)
call plotline(bb,ans,eans)

!lets make a Kernel for the Gaussian process
allocate(Kernel(npt,npt))
npars=2
allocate(pars(npars))
pars(1)=1.0d0
pars(2)=7.0d0
call makekernel(Kernel,npt,x,yerr,npars,pars)

call pgpage()
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt,npt,Kernel,bpix)

!Cholesky factorization
allocate(Kfac(npt,npt))
Kfac=Kernel
call dpotrf('U',npt,Kfac,npt,info)


call pgclos()

end program gptest

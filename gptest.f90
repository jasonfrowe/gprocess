program gptest
!Jason F Rowe 2016 (jasonfrowe@gmail.com)
use precision
implicit none
integer :: iargc,nmax,npt,npars,info,i,npt1,nrhs
real(double) :: chisq,bpix,ldet
real, allocatable, dimension(:) :: bb !contains dimensions of plot
real(double), allocatable, dimension(:) :: x,y,yerr,ans,eans,pars,alpha,&
   mu,x1,dbb,yerr2,std,alphanew
real(double), allocatable, dimension(:,:) :: Kernel,Kfac,newKernel,     &
   newKernelT,cov
character(80) :: filename

!These are F90 interfaces that allow one to pass assumed sized arrays
!to subroutines.
interface !reads in a three-column ascii space seperated file
   subroutine getdata(filename,npt,nmax,x,y,yerr)
      use precision
      implicit none
      integer, intent(inout) :: npt,nmax
      real(double), dimension(:), intent(inout) :: x,y,yerr
      character(80), intent(in) :: filename
   end subroutine getdata
end interface
interface !makes a plot of your data.
   subroutine plotdata(npt,x,y,yerr,bb)
      use precision
      implicit none
      integer, intent(in) :: npt
      real, dimension(:), intent(inout) :: bb
      real(double), dimension(:), intent(in) :: x,y,yerr
   end subroutine plotdata
end interface
interface !fits a straight line to data
   subroutine fitline(npt,x,y,yerr,ans,eans,chisq)
      use precision
      implicit none
      integer, intent(in) :: npt
      real(double), intent(inout) :: chisq
      real(double), dimension(:), intent(in) :: x,y,yerr
      real(double), dimension(:), intent(inout) :: ans,eans
   end subroutine fitline
end interface
interface !plots a line based on output from fitline
   subroutine plotline(bb,ans,eans)
      use precision
      implicit none
      real, dimension(:), intent(in) :: bb
      real(double), dimension(:), intent(in) :: ans,eans
   end subroutine plotline
end interface
interface !creates a co-variance matrix
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(in) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(in) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface
interface !displays a 2D array as a picture
   subroutine displaykernel(nx,ny,Kernel,bpix)
      use precision
      implicit none
      integer, intent(inout) :: nx,ny
      real(double), intent(inout) :: bpix
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine displaykernel
end interface
interface !plots samples and uncertainties
   subroutine plotsamples(npt1,x1,mu,std)
      use precision
      implicit none
      integer, intent(in) :: npt1
      real(double), dimension(:), intent(in) :: x1,mu,std
   end subroutine plotsamples
end interface

!check that we have enough information from the commandline
if(iargc().lt.1)then !if not, spit out the Usage info and stop.
   write(0,*) "Usage: gptest filename"
   stop
endif

!read in filename containing data (3 columns)
call getarg(1,filename)

nmax=1000 !initial guess for number of datapoints.
allocate(x(nmax),y(nmax),yerr(nmax)) !allocate arrays
call getdata(filename,npt,nmax,x,y,yerr) !subroutine to read in data
write(0,*) "Number of points read: ",npt !report how much data was read in

!open PGPLOT device
call pgopen('?')  !'?' lets the user choose the device.
call PGPAP (8.0 ,1.0) !use a square 8" across
call pgpage() !create a fresh page
call pgslw(3) !thicker lines

!plot the data
allocate(bb(4)) !contains plot boundaries
bb=0.0e0 !tell code to generate scale for plot
bb(1)=-25.0  !if bb(1)!=bb(2) then these values are used for plotting
bb(2)=200.0
bb(3)=-30.0  !if bb(3)!=bb(4) then these values are used for plotting
bb(4)= 14.0
call plotdata(npt,x,y,yerr,bb) !plot data

!fit a straight line
allocate(ans(2),eans(2)) !contains fit for a straight line
call fitline(npt,x,y,yerr,ans,eans,chisq) !calculates bestfit
call plotline(bb,ans,eans) !plot the best-fit line on top of our data.

!repeat but increase error via chi-squared
call pgpage() !make a fresh page for plotting
allocate(yerr2(npt)) !allocate array to hold scaled errors
yerr2(1:npt)=yerr(1:npt)*sqrt(chisq/dble(npt)) !scale errors such that k=1

!plot the data
call plotdata(npt,x,y,yerr2,bb) !plot the data with new error bars

!fit a straight line
call fitline(npt,x,y,yerr2,ans,eans,chisq) !fit data with new errorbars
call plotline(bb,ans,eans) !plot the new fit and uncertainties

!lets make a Kernel/co-variance for the Gaussian process
allocate(Kernel(npt,npt)) !allocate space
npars=4 !number of parameters used for model of the matrix
allocate(pars(npars))
pars(1)=1.0d0 !amp scale for exp
pars(2)=7.0d0 !length scale for exp
pars(3)=100.0 !second amp scale
pars(4)=500.0 !second length scale
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars) !create Kernel

call pgpage() !create fresh page for plotting
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt,npt,Kernel,bpix) !show what the Kernel looks like

!Cholesky factorization
allocate(Kfac(npt,npt)) !allocate array to contain Cholesky factorization
Kfac=Kernel  !make a copy of Kernel as dpotrf overwrites input
call dpotrf('U',npt,Kfac,npt,info) !LAPACK routine for Cholesky
if (info.ne.0) then !check for errors
   write(0,*) "Cholesky factorization failed"
   write(0,*) "dpotrf info: ",info
   stop
endif

!!example: calculate log determinant
!ldet=0.0
!do i=1,npt
!   ldet=ldet+log(Kfac(i,i))
!enddo
!ldet=2.0d0*ldet


!calculate solution
! We have to solve,
! Kernel*X=Y,
! for X.
! LAPACK routine dpotrs uses Kfac from dpotrf.  We copy y into a new
! array called alpha, which is overridden with the solution on completion.
allocate(alpha(npt))
alpha=y(1:npt) !dpotrs takes alpha as input and output.
nrhs=1 !how man columns does alpha contain - just one
call dpotrs('U',npt,nrhs,Kfac,npt,alpha,npt,info) !call LAPACK
if (info.ne.0) then !check for errors
   write(0,*) "Solver failed.."
   write(0,*) "dpotrs info: ",info
   stop
endif

call pgpage() !new page for plotting
npt1=1000 !lets generate a resampled predicted dataset based on our solution.
allocate(mu(npt1)) !allocate space for new predicted dataset
allocate(newKernel(npt1,npt),x1(npt1),dbb(4)) !need a new co-variance to match mu
dbb=dble(bb) !we are going to sample around the 'X'-axis of the plot. Need dbles
do i=1,npt1  !create new sampling intervals
   x1(i)=dbb(1)+(dbb(2)-dbb(1))/dble(npt1-1.0d0)*dble(i-1)
enddo
yerr2=0.0d0
call makekernel(newKernel,npt1,npt,x1,x,npt,yerr2,npars,pars) !make new Kernel
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt1,npt,newKernel,bpix) !display our new Kernel - not square!

!yerr2=0.0d0
!call makekernel(Kernel,npt,npt,x,x,npt,yerr2,npars,pars)
!mu=matmul(Kernel,alpha)
!do i=1,npt
!   write(0,*) y(i),mu(i)
!enddo

mu=matmul(newKernel,alpha) !to get the predicted dataset we multiply matrices



!lets calculate the uncertainty on the predicted dataset
call pgpage() !fresh page for plotting
allocate(newKernelT(npt,npt1))  !allocate space for transposed Kernel
newKernelT=transpose(newKernel) !create transpose
allocate(cov(npt1,npt1)) !allocate space for co-variance
!make Kernel based on predicted dataset
call makekernel(cov,npt1,npt1,x1,x1,npt,yerr,npars,pars)
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt1,npt1,cov,bpix) !display the Kernel
!now we call the Solver again.  We will feed in the new transposed Kernel
nrhs=npt1 !note, more than one column!
call dpotrs('U',npt,nrhs,Kfac,npt,newKernelT,npt,info) !Solve
if (info.ne.0) then !error check
   write(0,*) "Solver failed with newKernelT"
   write(0,*) "dpotrs info: ",info
   stop
endif
!calculate covariance matrix
cov=cov-matmul(newKernel,newKernelT)
allocate(std(npt1)) !allocate space for estimated uncertainties
do i=1,npt1
   std(i)=sqrt(cov(i,i)) !use diagonals as estimate of standard deviation
enddo

call pgpage()!fresh page for plotting
call plotdata(npt,x,y,yerr,bb) !plot our original dataset
call plotsamples(npt1,x1,mu,std) !plot our predicted sample set on top.
!call plotline(bb,ans,eans)
deallocate(mu,std,x1,yerr2,newKernel,cov,newKernelT)

npt1=2
allocate(mu(npt1),std(npt1),x1(npt1),newKernel(npt1,npt),yerr2(npt1))
x1(1)=0.0d0
x1(2)=x(npt)
yerr2=0.0d0
call makekernel(newKernel,npt1,npt,x1,x,npt,yerr2,npars,pars)
mu=matmul(newKernel,alpha)
allocate(newKernelT(npt,npt1))  !allocate space for transposed Kernel
newKernelT=transpose(newKernel) !create transpose
allocate(cov(npt1,npt1))
!make Kernel based on predicted dataset
call makekernel(cov,npt1,npt1,x1,x1,npt,yerr2,npars,pars)
!now we call the Solver again.  We will feed in the new transposed Kernel
nrhs=npt1 !note, more than one column!
call dpotrs('U',npt,nrhs,Kfac,npt,newKernelT,npt,info) !Solve
if (info.ne.0) then !error check
   write(0,*) "Solver failed with newKernelT"
   write(0,*) "dpotrs info: ",info
   stop
endif
!calculate covariance matrix
cov=cov-matmul(newKernel,newKernelT)
do i=1,npt1
   std(i)=sqrt(cov(i,i)) !use diagonals as estimate of standard deviation
enddo
500 format(A6,F6.2,A5,F6.2)
write(0,500) "mu1:  ",mu(1)," +/- ",std(1)
write(0,500) "mu2:  ",mu(2)," +/- ",std(2)
write(0,500) "diff: ",mu(1)-mu(2)," +/- ",                              &
   sqrt(std(1)*std(1)+std(2)+std(2))
write(0,500) "rate: ",(mu(1)-mu(2))/(x1(2)-x1(1))," +/- ",                &
   sqrt(std(1)*std(1)+std(2)+std(2))/(x1(2)-x1(1))

call pgclos() !close plotter

end program gptest

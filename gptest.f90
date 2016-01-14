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
!   subroutine makekernel(Kernel,npt,x,yerr,npars,pars)
   subroutine makekernel(Kernel,npt1,npt2,x1,x2,npt,yerr,npars,pars)
      use precision
      implicit none
      integer, intent(in) :: npt1,npt2,npars,npt
      real(double), dimension(:), intent(in) :: x1,x2,yerr,pars
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine makekernel
end interface
interface
   subroutine displaykernel(nx,ny,Kernel,bpix)
      use precision
      implicit none
      integer, intent(inout) :: nx,ny
      real(double), intent(inout) :: bpix
      real(double), dimension(:,:), intent(inout) :: Kernel
   end subroutine displaykernel
end interface
interface
   subroutine plotsamples(npt1,x1,mu,std)
      use precision
      implicit none
      integer, intent(in) :: npt1
      real(double), dimension(:), intent(in) :: x1,mu,std
   end subroutine plotsamples
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
bb=0.0e0 !tell code to generate scale for plot
bb(1)=-10.0
bb(2)=200.0
bb(3)=-30.0
bb(4)= 14.0
call plotdata(npt,x,y,yerr,bb)

!fit a straight line
allocate(ans(2),eans(2)) !contains fit for a straight line
call fitline(npt,x,y,yerr,ans,eans,chisq)
call plotline(bb,ans,eans)

!repeat but increase error via chi-squared
call pgpage()
allocate(yerr2(npt))
yerr2(1:npt)=yerr(1:npt)*sqrt(chisq/dble(npt))

!plot the data
call plotdata(npt,x,y,yerr2,bb)

!fit a straight line
call fitline(npt,x,y,yerr,ans,eans,chisq)
call plotline(bb,ans,eans)

!lets make a Kernel for the Gaussian process
allocate(Kernel(npt,npt))
npars=5
allocate(pars(npars))
pars(1)=1.0d0 !amp scale for exp
pars(2)=7.0d0 !length scale for exp
pars(3)=100.0    !m
pars(4)=500.0   !b
call makekernel(Kernel,npt,npt,x,x,npt,yerr,npars,pars)

call pgpage()
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt,npt,Kernel,bpix)

!Cholesky factorization
allocate(Kfac(npt,npt))
Kfac=Kernel
call dpotrf('U',npt,Kfac,npt,info)
if (info.ne.0) write(0,*) "dpotrf info: ",info

!calculate log determinant
ldet=0.0
do i=1,npt
   ldet=ldet+log(Kfac(i,i))
enddo
ldet=2.0d0*ldet


!calculate A*alpha=y
allocate(alpha(npt))
alpha=y(1:npt) !dpotrs takes alpha as input and output.
nrhs=1
call dpotrs('U',npt,nrhs,Kfac,npt,alpha,npt,info)
if (info.ne.0) write(0,*) "dpotrs info: ",info

npt1=100
allocate(mu(npt1))
allocate(newKernel(npt1,npt),x1(npt1),dbb(4))
dbb=dble(bb)
do i=1,npt1
   x1(i)=dbb(1)+(dbb(2)-dbb(1))/dble(npt1-1.0d0)*dble(i-1)
enddo
call makekernel(newKernel,npt1,npt,x1,x,npt,yerr,npars,pars)
call pgpage()
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt1,npt,newKernel,bpix)

mu=matmul(newKernel,alpha)
!do i=1,npt1
!   write(0,*) i,x1(i),mu(i)
!enddo

allocate(newKernelT(npt,npt1))
newKernelT=transpose(newKernel)
allocate(cov(npt1,npt1))
call makekernel(cov,npt1,npt1,x1,x1,npt,yerr,npars,pars)
call pgpage()
bpix=1.0e30 !cut off for bright pixels
call displaykernel(npt1,npt1,cov,bpix)
!cov -= np.dot(Kxs, self.solver.apply_inverse(KxsT, in_place=True))
nrhs=npt1
call dpotrs('U',npt,nrhs,Kfac,npt,newKernelT,npt,info)
if (info.ne.0) write(0,*) "dpotrs info: ",info

cov=cov-matmul(newKernel,newKernelT)
allocate(std(npt1))
do i=1,npt1
   std(i)=sqrt(cov(i,i))
!   write(0,*) "std: ",mu(i),std(i)
enddo


call pgpage()
call plotdata(npt,x,y,yerr,bb)
call plotsamples(npt1,x1,mu,std)

call pgclos()

end program gptest

subroutine makekernel(Kernel,npt,x,yerr,npars,pars)
use precision
implicit none
integer :: npt,npars,i,j
real(double) :: twol2
real(double), dimension(:) :: x,yerr,pars
real(double), dimension(:,:) :: Kernel

!initialize Kernel to zero
Kernel=0.0d0

!pre-compute
twol2=2.0d0*pars(2)*pars(2)

!add in diagonal from yerr
do i=1,npt
   Kernel(i,i)=yerr(i)*yerr(i)
enddo

do i=1,npt
   do j=1,npt
      Kernel(i,j)=Kernel(i,j)+pars(1)*exp(-(x(i)-x(j))*(x(i)-x(j))/twol2)
   enddo
enddo

return
end subroutine makekernel

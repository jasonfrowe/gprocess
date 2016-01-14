subroutine makekernel(Kernel,npt1,npt2,x1,x2,yerr,npars,pars)
use precision
implicit none
integer :: npt1,npt2,npars,i,j
real(double) :: twol2
real(double), dimension(:) :: x1,x2,yerr,pars
real(double), dimension(:,:) :: Kernel

!initialize Kernel to zero
Kernel=0.0d0

!pre-compute
twol2=2.0d0*pars(2)*pars(2)

do i=1,npt1
   do j=1,npt2
      Kernel(i,j)=Kernel(i,j)+pars(1)*exp(-(x1(i)-x2(j))*(x1(i)-x2(j))/twol2)
!      write(0,*) i,j,x1(i)-x2(j)
      if(abs(x1(i)-x2(j)).le.1.0d-15)then
!         write(0,*) "yo"
         Kernel(i,j)=Kernel(i,j)+yerr(i)*yerr(j)
      endif
   enddo
enddo

return
end subroutine makekernel

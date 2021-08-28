!======================================================================
!
! SUBROUTINE: Solve linear system with SOR method
!             SOR: Successive Over-Relaxation
!
!      [in ] nval  : size of coefficient matrix.
!            adim  : coefficient matrix.
!            bdim  : right-hand-side vector.
!            omega : acceleration factor.
!            itmax : maximum number of iteration.
!            eplson: convergence test value.
!            iopt  : only if iopt = .ture., display residual history.
!
!      [out] xdim  : solution of linear system.
!            itnum : number of iterations.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine sor(nval,adim,bdim,omega,itmax,epsln,iopt,xdim,itnum)
  implicit none
  integer,intent(in ) :: nval,itmax
  real(8),intent(in ) :: adim(nval,nval),bdim(nval),omega,epsln
  logical,intent(in ) :: iopt
  integer,intent(out) :: itnum
  real(8),intent(out) :: xdim(nval)
  integer :: i,j,k
  real(8) :: sum1,sum2,dnorm2,xold(nval),xi
  real(8),parameter :: zero=0.0d0
  !
  if( iopt ) write(*,'(/a)') '.... Residual history ...............'
  !
  dnorm2 = zero
  do k=1, itmax
     do i=1, nval
        sum1 = zero
        do j=1, i-1
           sum1 = sum1 + adim(i,j) * xdim(j)
        end do
        !
        sum2 = zero
        do j = i+1, nval
           sum2 = sum2 + adim(i,j) * xdim(j)
        end do
        !
        xi = (bdim(i) - sum2 - sum1) / adim(i,i)
        xdim(i) = xdim(i) + omega * (xi - xdim(i))
     end do
     !
     dnorm2 = zero
     do i=1, nval
        dnorm2 = dnorm2 + (xdim(i) - xold(i))**2
     end do
     !
     if( iopt ) write(*,'(i4,2x,1pe12.5)') k, sqrt(dnorm2)
     !
     if( sqrt(dnorm2) .lt. epsln ) then
        exit
     else
        xold(:) = xdim(:)
     end if
  end do
  itnum = k
  !
  if( sqrt(dnorm2) >  epsln )then
     write(*,'(a)') '... Not converge! ...'
  else if( iopt )then
     write(*,'(a/)') '....................................'
  end if
  !
  return
end subroutine sor

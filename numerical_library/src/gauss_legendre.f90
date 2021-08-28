!======================================================================
!
! SUBROUTINE: Calculate integral point and its weighting coefficient
!             for Gauss-Ledendre quadrature
!
!      [in ] nval  : number of integration points.
!            lmax  : maximum number of iteration.
!            eplson: convergence test value.
!
!      [out] gzi   : integration point.
!            weight: weight on integration point.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine gauss_legendre(nval,lmax,epsln,gzi,weight)
  implicit none
  integer,intent(in ) :: nval,lmax
  real(8),intent(in ) :: epsln
  real(8),intent(out) :: gzi(nval),weight(nval)
  integer :: i,j
  real(8) :: pi,x0,dx,pn,dp
  real(8),parameter   :: zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0
  pi=four*atan(one)
  !
  do i=1, nval
     ! ... initial value
     x0 = dble(i-1) + 0.75d0
     x0 = x0 / (dble(nval) + 0.5d0)
     x0 = cos(x0*pi)
     ! ... Newton method
     do j=1, lmax
        call legendre_func(x0,nval,pn,dp)
        dx = pn / dp
        x0 = x0 - dx
        if( abs(dx) < epsln ) exit
     end do
     gzi(i) = x0
     ! ... weighting coefficient
     call legendre_func(gzi(i),nval-1,pn,dp)
     weight(i) = dble(nval)*pn
     weight(i) = two*(one - gzi(i)**2) / weight(i)**2 
  end do
  !
  return
contains
  subroutine legendre_func(x,n,pn,dp)
    implicit none
    integer,intent(in ) :: n
    real(8),intent(in ) :: x
    real(8),intent(out) :: pn,dp
    integer :: i
    real(8) :: p0,p1
    real(8),parameter   :: zero=0.0d0,one=1.0d0
    !
    pn = zero
    dp = zero
    if( n == 0 ) then
       pn = one
       dp = zero
    else if( n == 1 ) then
       pn = x
       dp = one
    else
       p0 = one
       p1 = x
       do i=2, n
          pn = dble(2*i-1)*x*p1 - dble(i-1)*p0
          pn = pn / dble(i)
          p0 = p1
          p1 = pn
       end do
       dp = dble(n)*(p0 - x*p1)
       dp = dp / (one - x**2)
    end if
    !
    return
  end subroutine legendre_func
end subroutine gauss_legendre

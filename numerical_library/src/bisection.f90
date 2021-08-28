!======================================================================
!
! SUBROUTINE: Solve transcendental function with bisectoin method
!
!      [in ] y0    : lower limit of interval including root.
!            y1    : upper limit of interval including root.
!            val   : value of right-hand side in Eq. f(x) = v.
!            eplson: convergence test value.
!            lmax  : maximum number of iterations.
!            func  : subroutine for calcutaing function f(x).
!
!      [out] xval  : root of Eq. f(x) = v.
!            icode : only if converged, return icode = 0. 
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine bisection(y0,y1,xval,val,epsln,lmax,func,icode)
  integer,intent(in ) :: lmax
  integer,intent(out) :: icode
  real(8),intent(in ) :: y0,y1,val,epsln
  real(8),intent(out) :: xval
  integer :: i
  real(8) :: x0,x1,x2,f0,f1,f2
  real(8),parameter   :: zero=0.0d0,two=2.0d0
  interface
     subroutine func(x,f)
       real(8),intent(in ) :: x
       real(8),intent(out) :: f
     end subroutine func
  end interface
  icode = 100
  !
  x0 = y0
  x1 = y1
  do i=1, lmax
     call func(x0,f0)
     call func(x1,f1)
     f0 = f0 - val
     f1 = f1 - val
     if( f0*f1 > zero ) then
        write(*,'(/a/)') 'Bisection: sgn(f0) == sgn(f1)!'
        exit
     end if
     !
     x2 = (x0 + x1) / two
     call func(x2,f2)
     f2 = f2 - val
     if( f2 == zero ) then
        x0 = x2
        exit
     else if( f0*f2 < zero ) then
        x1 = x2
        f1 = f2
     else
        x0 = x2
        f0 = f2
     end if
     !
     if( abs(x0-x1) < epsln ) then
        icode = 0
        exit
     end if
  end do
  !
  xval = x0
  !
  return
end subroutine bisection

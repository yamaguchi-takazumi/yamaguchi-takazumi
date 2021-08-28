!======================================================================
!
! SUBROUTINE: Solve transcendental function with Newton method
!
!      [in ] y0    : initial value of root of Eq. f(x) = v.
!            val   : value of right-hand side in Eq. f(x) = v.
!            eplson: convergence test value.
!            lmax  : maximum number of iteration.
!            func  : subroutine for calcutaing function f(x).
!
!      [out] xval  : root of Eq. f(x) = v.
!            icode : only if converged, return icode = 0. 
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine newton(y0,xval,val,epsln,lmax,func,icode)
  integer,intent(in ) :: lmax
  integer,intent(out) :: icode
  real(8),intent(in ) :: y0,val,epsln
  real(8),intent(out) :: xval
  integer :: i
  real(8) :: x0,f0,dx,df
  interface
     subroutine func(x,f,df)
       real(8),intent(in ) :: x
       real(8),intent(out) :: f,df
     end subroutine func
  end interface
  icode = 100
  !
  x0 = y0
  do i=1, lmax
     call func(x0,f0,df)
     f0 = f0 - val
     dx = f0 / df
     !
     x0 = x0 - dx
     !
     if( abs(dx) < epsln ) then
        icode = 0
        exit
     end if
  end do
  !
  xval = x0
  !
  return
end subroutine newton

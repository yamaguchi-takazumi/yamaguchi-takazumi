!======================================================================
!
! SUBROUTINE: Solve linear system with LU decomposition
!
!      [in ] nval: size of square matrix.
!            adim: LU decomposed matrix.
!            bdim: right-hand-side vector.
!            ip  : array for controlling pivoting.
!
!      [out] xdim: solution of linear system.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine lu_solve(nval,adim,bdim,ip,xdim)
  implicit none 
  integer,intent(in ) :: nval,ip(nval)
  real(8),intent(in ) :: adim(nval,nval),bdim(nval)
  real(8),intent(out) :: xdim(nval)
  integer :: i,k
  real(8) :: sum,ydim(nval)
  real(8),parameter :: zero=0.0d0
  ! ... Solve Ly=Pb
  ydim(1) = bdim(ip(1))
  do i=2, nval
     sum = zero
     do k=1, i-1
        sum = sum + adim(i,k) * ydim(k)
     end do
     ydim(i) = bdim(ip(i)) - sum
  end do
  ! ... Solve Ux=y
  xdim(nval) = ydim(nval) / adim(nval,nval)
  do i=nval-1, 1, -1
     sum = zero
     do k=i+1, nval
        sum = sum + adim(i,k) * xdim(k)
     end do
     xdim(i) = (ydim(i) - sum) / adim(i,i)
  end do
  !
  return
end subroutine lu_solve

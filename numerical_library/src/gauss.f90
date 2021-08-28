!======================================================================
!
! SUBROUTINE: Solve linear system by using Gauss elimination with
!             LU decomposition
!
!      [in ] nval  : size of coefficient matrix
!            adim  : coefficient matrix
!            bdim  : right-hand-side vector
!
!      [out] xdim  : solution of linear system.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine gauss(nval,adim,bdim,xdim)
  implicit none
  integer,intent(in ) :: nval
  real(8),intent(in ) :: adim(nval,nval),bdim(nval)
  real(8),intent(out) :: xdim(nval)
  integer :: ip(nval)
  real(8) :: adim0(nval,nval)
  !
  adim0(:,:) = adim(:,:)
  call lu_decomp(nval,adim0,ip)
  call lu_solve(nval,adim0,bdim,ip,xdim)
  !
  return
end subroutine gauss

!======================================================================
!
! SUBROUTINE: Solve linear system with conjugate gradient method
!
!      [in ] nval  : size of coefficient matrix.
!            adim  : coefficient matrix.
!            bdim  : right-hand-side vector.
!            itmax : maximum number of iteration.
!            eplson: convergence test value.
!            iopt  : only if iopt = .ture., display residual history.
!
!      [out] xdim  : solution of linear system.
!            itnum : number of iterations.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine cgm(nval,adim,bdim,itmax,epsln,iopt,xdim,itnum)
  implicit none
  integer,intent(in ) :: nval,itmax
  real(8),intent(in ) :: adim(nval,nval),bdim(nval),epsln
  logical,intent(in ) :: iopt
  integer,intent(out) :: itnum
  real(8),intent(out) :: xdim(nval)
  integer :: k
  real(8) :: ydim(nval),alph,beta,rdim(nval),pdim(nval),rerror
  real(8),parameter :: zero=0.0d0
  !
  rdim(:) = bdim(:) - matmul(adim,xdim)
  pdim(:) = rdim(:)
  if( iopt ) write(*,'(/a)') '.... Residual history ...............'
  !
  do k=1, itmax
     ydim(:) = matmul(adim,pdim)
     beta = dot_product(rdim,rdim)
     alph = beta/dot_product(pdim,ydim)
     !
     xdim(:) = xdim(:) + alph*pdim(:)
     rdim(:) = rdim(:) - alph*ydim(:)
     !
     rerror = dot_product(rdim,rdim)
     if( iopt ) write(*,'(i4,2x,1pe12.5)') k, rerror
     !
     if( rerror .lt. epsln ) exit
     !
     beta = dot_product(rdim,rdim)/beta
     pdim(:) = rdim(:) + beta*pdim(:)
  end do
  itnum = k
  !
  if( rerror > epsln ) write(*,'(a)') '... Not converge! ...'
  if( iopt ) write(*,'(a/)') '....................................'
  !
  return
end subroutine cgm

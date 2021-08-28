!======================================================================
!
! SUBROUTINE: Convert sparse matrix into CRS type matrix.
!
!      [in ] m      : number of rows.
!            n      : number of columns.
!            nsize  : estimated number of nonzero elements.
!            adim   : real matrix.
!
!      [out] nonzero: actual number of nonzero elements.
!            val    : value of nonzero element.
!            col_ind: column number of nonzero element.
!            row_prt: number of nonzero elements per row.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine crs_conv(m,n,nsize,adim,nonzero,val,col_ind,row_ptr)
  implicit none
  integer,intent(in ) :: m,n,nsize
  real(8),intent(in ) :: adim(m,n)
  integer,intent(out) :: nonzero,col_ind(nsize),row_ptr(m+1)
  real(8),intent(out) :: val(nsize)
  integer :: i,j,ic
  real(8),parameter   :: zero=0.0d0
  !
  val = zero
  col_ind(:) = 0
  row_ptr(:) = 0
  !
  ic = 0
  do i=1, m
     row_ptr(i) = ic + 1
     do j=1, n
        if( abs(adim(i,j)) > zero ) then
           ic = ic + 1
           val(ic) = adim(i,j)
           col_ind(ic) = j
        end if
     end do
  end do
  row_ptr(i) = ic + 1
  nonzero = ic
  !
  return
end subroutine crs_conv

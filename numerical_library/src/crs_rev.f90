!======================================================================
!
! SUBROUTINE: Convert CRS type matrix into sparse matrix.
!
!      [in ] m      : number of rows.
!            n      : number of columns.
!            nonzero: actual number of nonzero elements.
!            val    : value of nonzero element.
!            col_ind: column number of nonzero element.
!            row_prt: number of nonzero elements per row.
!
!      [out] adim   : real matrix.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine crs_rev(m,n,nonzero,val,col_ind,row_ptr,adim)
  implicit none
  integer,intent(in ) :: m,n,nonzero,col_ind(nonzero),row_ptr(m+1)
  real(8),intent(in ) :: val(nonzero)
  real(8),intent(out) :: adim(m,n)
  integer :: i,j
  real(8),parameter   :: zero=0.0d0
  !
  adim(:,:) = zero
  !
  do i=1, m
     do j=row_ptr(i), row_ptr(i+1)-1
        adim(i,col_ind(j)) = val(j)
     end do
  end do
  !
  return
end subroutine crs_rev

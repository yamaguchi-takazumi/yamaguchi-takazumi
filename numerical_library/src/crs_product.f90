!======================================================================
!
! SUBROUTINE: Calculate product of matrix A and vector x.
!
!      [in ] m      : number of rows.
!            n      : number of columns.
!            nonzero: number of nonzero elements.
!            val    : value of nonzero element.
!            col_ind: column number of nonzero element.
!            row_prt: number of nonzero elements per row.
!            xvec   : real vector.
!
!      [out] ax     : product of matrix A and vector x.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine crs_product(m,n,nonzero,val,col_ind,row_ptr,xvec,ax)
  implicit none
  integer,intent(in ) :: m,n,nonzero,col_ind(nonzero),row_ptr(m+1)
  real(8),intent(in ) :: val(nonzero),xvec(n)
  real(8),intent(out) :: ax(m)
  integer :: i,j
  real(8),parameter   :: zero=0.0d0
  !
  ax(:) = zero
  do i=1, m
     do j=row_ptr(i),row_ptr(i+1)-1
        ax(i) = ax(i) + val(j)*xvec(col_ind(j))
     end do
  end do
  !
  return
end subroutine crs_product

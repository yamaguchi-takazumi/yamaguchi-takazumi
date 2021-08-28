!======================================================================
!
! SUBROUTINE: Calculate add C of matrixes A and B.
!
!      [in ] m      : number of rows.
!            n      : number of columns.
!            nzero_A: number of nonzero elements in matrix A.
!            nzero_B: number of nonzero elements in matrix B.
!            val_A  : value of nonzero element in matrix A.
!            val_B  : value of nonzero element in matrix B.
!            col_A  : column number of nonzero element in matrix A.
!            col_B  : column number of nonzero element in matrix B.
!            row_A  : number of nonzero elements per row in matrix A.
!            row_B  : number of nonzero elements per row in matrix B.
!            nsize  : estimated number of nonzero elements in matrix C.
!
!      [out] nzero_C: number of nonzero elements in matrix C.
!            val_C  : value of nonzero element in matrix C.
!            col_C  : column number of nonzero element in matrix C.
!            row_C  : number of nonzero elements per row in matrix C.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine crs_add(m,n,nzero_A,nzero_B,val_A,val_B,col_A,col_B, &
     row_A,row_B,nsize,nzero_C,col_C,row_C,val_C)
  implicit none
  integer,intent(in ) :: m,n,nzero_A,nzero_B,col_A(nzero_A), &
       col_B(nzero_B),row_A(m+1),row_B(m+1),nsize
  real(8),intent(in ) :: val_A(nzero_A),val_B(nzero_B)
  integer,intent(out) :: nzero_C,col_C(nsize),row_C(m+1)
  real(8),intent(out) :: val_C(nsize)
  integer :: i,j,ic
  real(8) :: val_tmp(n)
  real(8),parameter   :: zero=0.0d0
  !
  ic = 0
  do i=1, m
     val_tmp(:) = zero
     do j=row_A(i), row_A(i+1)-1
        val_tmp(col_A(j)) = val_A(j)
     end do
     do j=row_B(i), row_B(i+1)-1
        val_tmp(col_B(j)) = val_tmp(col_B(j)) + val_B(j)
     end do
     !
     row_C(i) = ic + 1
     do j=1, n
        if( abs(val_tmp(j)) > zero ) then
           ic = ic + 1
           val_C(ic) = val_tmp(j)
           col_C(ic) = j
        end if
     end do
  end do
  nzero_C  = ic
  row_C(i) = ic + 1
  !
  return
end subroutine crs_add

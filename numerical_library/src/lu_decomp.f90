!======================================================================
!
! SUBROUTINE to apply LU decomposition to square matrix
!
!      [in ] nval: size of square matrix.
!            adim: square matrix.
!
!      [out] adim: LU decomposed matrix.
!            ip  : array for controlling pivoting.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine lu_decomp(nval,adim,ip)
  integer,intent(in   ) :: nval
  integer,intent(  out) :: ip(nval)
  real(8),intent(inout) :: adim(nval,nval)
  integer :: i,j,k,maxLine,temp
  real(8) :: reserve
  !
  do k=1, nval
     ip(k) = k
  end do
  ! ... Pivoting
  do k=1, nval-1
     maxLine = k
     do i = k+1, nval
        if (abs(adim(i,k)) .gt. abs(adim(k,k))) then
           maxLine = i
        end if
     end do
     if( maxLine .ne. k )then
        temp = ip(k)
        ip(k) = ip(maxLine)
        ip(maxLine) = temp
        do j=1, nval
           reserve = adim(k,j)
           adim(k,j) = adim(maxLine,j)
           adim(maxLine,j) = reserve
        end do
     end if
     ! ... LU decomposion
     do i=k+1, nval
        adim(i,k) = adim(i,k) / adim(k,k)
     end do
     do j=k+1, nval
        do i=k+1, nval
           adim(i,j) = adim(i,j) - adim(i,k) * adim(k,j)
        end do
     end do
  end do
  !
  return
end subroutine lu_decomp

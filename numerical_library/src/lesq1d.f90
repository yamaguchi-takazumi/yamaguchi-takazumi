!======================================================================
!
! SUBROUTINE: Approximate data set as linear function by using
!             least squares method
!
!      [in ] ndata: number of data.
!            x, y : data.
!
!      [out] aval : slope of linear equation y = a*x + b.
!            bval : y-intercept of linear equation y = a*x + b.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine lesq1d(ndata,x,y,aval,bval)
  implicit none
  integer,intent(in ) :: ndata
  real(8),intent(in ) :: x(ndata),y(ndata)
  real(8),intent(out) :: aval,bval
  integer :: i
  real(8) :: sum1,sum2,sum3,sum4
  real(8),parameter   :: zero=0.0d0
  !
  sum1 = zero
  sum2 = zero
  sum3 = zero
  sum4 = zero
  do i=1, ndata
     sum1 = sum1 + x(i)
     sum2 = sum2 + y(i)
     sum3 = sum3 + x(i)*y(i)
     sum4 = sum4 + x(i)**2
  end do
  !
  aval = (dble(ndata)*sum3 - sum1*sum2) / (dble(ndata)*sum4 - sum1**2)
  bval = (sum2*sum4 - sum1*sum3) / (dble(ndata)*sum4 - sum1**2)
  !
  return
end subroutine lesq1d
